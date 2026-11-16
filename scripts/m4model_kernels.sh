#!/bin/bash
# Cortex-M4 cycle model (pqm4/m4model/m4cycles.c) for every operation of this repository's own bench:
# KEM operations including the expanded-key decapsulation, and the hot kernels, for both backends.
#   bash scripts/m4model_kernels.sh [out.txt]
# Builds platform/cortex-m/build/<backend>-<N>-cortex-m4/bench.elf if missing and runs each operation
# once in --insn mode under QEMU mps2-an386; the setup-only run ("none") is subtracted.
# QEMU / QEMU_INC / PLUGIN_ARGS as in pqm4/m4model/run_m4model.sh.
# BACKENDS (default "ref armv7em m4f") and SIZES (default "512 1024") restrict what is measured.
set -euo pipefail
cd "$(dirname "$0")/.." || exit 1
OUT=${1:-build/m4model_kernels.txt}
QEMU=${QEMU:-/opt/qemu-plugins/bin/qemu-system-arm}
QEMU_INC=${QEMU_INC:-/opt/qemu-src/qemu-9.1.3/include/qemu}
PLUGIN_ARGS=${PLUGIN_ARGS:-}
OPS="none keypair enc dec dec_exp expkey decode_trap ntt257 intt257 basemul257 ntt3329 intt3329 basemul3329 ntt12289 intt12289 rot16 sha3_256"
mkdir -p build "$(dirname "$OUT")"
SO=$(pwd)/build/m4cycles.so
gcc -O2 -shared -fPIC -Wall -I"$QEMU_INC" $(pkg-config --cflags glib-2.0) -o "$SO" pqm4/m4model/m4cycles.c
tmp=$(mktemp -d)
trap 'rm -rf -- "$tmp"' EXIT
results=$tmp/results.txt
for n in ${SIZES:-512 1024}; do
    for bk in ${BACKENDS:-ref armv7em m4f}; do
        elf=platform/cortex-m/build/$bk-$n-cortex-m4/bench.elf
        [ -f "$elf" ] || make -s -C platform/cortex-m END_N=$n BACKEND=$bk > /dev/null
        for op in $OPS; do
            ( timeout "${QEMU_TIMEOUT:-1800}" "$QEMU" -M mps2-an386 -cpu cortex-m4 -nographic -monitor none -serial none \
                -semihosting-config "enable=on,target=native,arg=bench,arg=1,arg=--insn,arg=$op" \
                -plugin "$SO$PLUGIN_ARGS" -d plugin -D "$tmp/$bk-$n-$op.log" -kernel "$elf" > /dev/null ) &
        done
        wait
        {
            echo "### Cortex-M4 model END-$n backend=$bk"
            printf "%-14s %14s %14s %7s\n" operation cycles instructions CPI
            c0=$(awk '/^cycles/{print $2}' "$tmp/$bk-$n-none.log"); i0=$(awk '/^insns/{print $2}' "$tmp/$bk-$n-none.log")
            for op in $OPS; do
                [ "$op" = none ] && continue
                c=$(awk '/^cycles/{print $2}' "$tmp/$bk-$n-$op.log"); i=$(awk '/^insns/{print $2}' "$tmp/$bk-$n-$op.log")
                [ -n "$c" ] || { echo "missing result for $bk-$n $op" >&2; exit 1; }
                awk -v op="$op" -v c=$((c - c0)) -v i=$((i - i0)) 'BEGIN { printf "%-14s %14d %14d %7.3f\n", op, c, i, c / i }'
            done
        } | tee -a "$results"
    done
done
mv -- "$results" "$OUT"
