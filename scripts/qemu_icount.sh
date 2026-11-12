#!/bin/bash
# Dynamic instruction counts per operation with the QEMU TCG plugin "libinsn" (QEMU built with --enable-plugins).
#   qemu_icount.sh aarch64 <bench binary>              (user-mode: /opt/qemu-plugins/bin/qemu-aarch64)
#   qemu_icount.sh cortex-m <bench.elf> [cpu] [machine] (system-mode: /opt/qemu-plugins/bin/qemu-system-arm)
# Each operation is run once in the bench's --insn mode; the count of the setup-only run ("none") is subtracted.
set -euo pipefail
QEMU_DIR=${QEMU_DIR:-/opt/qemu-plugins/bin}
PLUGIN=${PLUGIN:-$(find /opt/qemu-src -name libinsn.so -print -quit 2>/dev/null || true)}
MODE=$1; BIN=$2; CPU=${3:-cortex-m4}; MACH=${4:-mps2-an386}
if [[ "$MODE" != aarch64 && "$MODE" != cortex-m ]]; then
    echo "Unknown QEMU mode: $MODE" >&2
    exit 2
fi
[[ -r "$PLUGIN" && -r "$BIN" ]] || { echo "Cannot read benchmark or instruction-count plugin" >&2; exit 2; }
OPS="none keypair enc dec dec_exp expkey decode_trap ntt257 intt257 basemul257 ntt3329 intt3329 basemul3329 ntt12289 intt12289 rot16 sha3_256"
declare -A CNT
out=$(mktemp)
trap 'rm -f -- "$out"' EXIT
for op in $OPS; do
    : > "$out"
    if [ "$MODE" = "aarch64" ]; then
        timeout "${QEMU_TIMEOUT:-1800}" "$QEMU_DIR/qemu-aarch64" -L /usr/aarch64-linux-gnu -plugin "$PLUGIN,inline=true" -d plugin -D "$out" "$BIN" 1 --insn "$op" > /dev/null
    else
        timeout "${QEMU_TIMEOUT:-1800}" "$QEMU_DIR/qemu-system-arm" -M "$MACH" -cpu "$CPU" -nographic -monitor none -serial none \
            -semihosting-config "enable=on,target=native,arg=bench,arg=1,arg=--insn,arg=$op" \
            -plugin "$PLUGIN,inline=true" -d plugin -D "$out" -kernel "$BIN" > /dev/null
    fi
    n=$(grep -oE '[0-9]+' "$out" | tail -1) || { echo "Missing instruction count for $op" >&2; exit 1; }
    CNT[$op]=$n
done
base=${CNT[none]}
printf "%-14s %14s\n" "operation" "instructions"
for op in $OPS; do
    [ "$op" = "none" ] && continue
    printf "%-14s %14d\n" "$op" $(( ${CNT[$op]} - base ))
done
