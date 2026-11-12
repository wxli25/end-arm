#!/bin/bash
# Dynamic instruction counts of KEMs inside pqm4 (QEMU mps2-an386 + TCG plugin libinsn).
#
#   bash pqm4/icount_pqm4.sh <pqm4-dir> <impl-path>...     e.g. crypto_kem/end512/m4 crypto_kem/ml-kem-512/m4fspeed
#
# Needs: pqm4/export_pqm4.py <pqm4> --icount (installs mupq/crypto_kem/icount_*.c), and a QEMU built with
# --enable-plugins (QEMU_DIR, default /opt/qemu-plugins/bin; PLUGIN = path of libinsn.so).
# Prints: impl, keypair, encaps, decaps (instructions), one line per implementation.
set -euo pipefail
P=${1:?usage: icount_pqm4.sh <pqm4-dir> <impl-path>...}
shift
QEMU_DIR=${QEMU_DIR:-/opt/qemu-plugins/bin}
PLUGIN=${PLUGIN:-$(find /opt/qemu-src -name libinsn.so -print -quit 2>/dev/null || true)}
[[ -r "$PLUGIN" ]] || { echo "libinsn.so not found (set PLUGIN=...)" >&2; exit 2; }
cd "$P"
out=$(mktemp)
trap 'rm -f -- "$out"' EXIT
count() {  # count <elf> -> total executed guest instructions
    : > "$out"
    timeout "${QEMU_TIMEOUT:-3600}" "$QEMU_DIR/qemu-system-arm" -M mps2-an386 -nographic -semihosting \
        -plugin "$PLUGIN,inline=true" -d plugin -D "$out" -kernel "$1" < /dev/null | grep -q "icount OK" \
        || { echo "run failed: $1" >&2; exit 1; }
    grep -oE '[0-9]+' "$out" | tail -1
}
printf "%-36s %14s %14s %14s\n" implementation keypair encaps decaps
for impl in "$@"; do
    name=${impl//\//_}
    for op in none keypair enc dec; do
        make -s PLATFORM=mps2-an386 IMPLEMENTATION_PATH="$impl" "elf/${name}_icount_${op}.elf" > /dev/null 2>&1 \
            || { echo "build failed: $impl icount_$op" >&2; exit 1; }
    done
    n0=$(count "elf/${name}_icount_none.elf")
    n1=$(count "elf/${name}_icount_keypair.elf")
    n2=$(count "elf/${name}_icount_enc.elf")
    n3=$(count "elf/${name}_icount_dec.elf")
    printf "%-36s %14d %14d %14d\n" "$impl" $((n1 - n0)) $((n2 - n1)) $((n3 - n2))
done
