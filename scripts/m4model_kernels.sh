#!/bin/bash
# Cortex-M4 cycle model (pqm4/m4model/m4cycles.c) for every operation of this repository's own bench:
# KEM operations including the expanded-key decapsulation, and the hot kernels, for both backends.
#   bash scripts/m4model_kernels.sh [out.txt]
# Builds the current sources and runs each operation once under QEMU mps2-an386;
# the setup-only run ("none") is subtracted. The DSP prefix prototype runs separately.
# QEMU / QEMU_INC / PLUGIN_ARGS as in pqm4/m4model/run_m4model.sh.
# BACKENDS (default "ref armv7em m4f") and SIZES (default "512 1024") restrict what is measured.
set -euo pipefail
cd "$(dirname "$0")/.." || exit 1
OUT=${1:-build/m4model_kernels.txt}
QEMU=${QEMU:-/opt/qemu-plugins/bin/qemu-system-arm}
QEMU_INC=${QEMU_INC:-/opt/qemu-src/qemu-9.1.3/include/qemu}
PLUGIN_ARGS=${PLUGIN_ARGS:-}
OPS="none keypair enc dec dec_exp expkey decode_trap ntt257 intt257 basemul257 ntt3329 intt3329 basemul3329 ntt12289 intt12289 rot16 sha3_256 crt_c2 pk_decode"
mkdir -p build "$(dirname "$OUT")"
SO=$(pwd)/build/m4cycles.so
gcc -O2 -shared -fPIC -Wall -I"$QEMU_INC" $(pkg-config --cflags glib-2.0) -o "$SO" pqm4/m4model/m4cycles.c
tmp=$(mktemp -d)
trap 'rm -rf -- "$tmp"' EXIT
results=$tmp/results.txt

measure_ops() {
    local elf=$1 args=$2 op pid cycles insns base_cycles base_insns
    shift 2
    local -a jobs=()
    for op in "$@"; do
        timeout "${QEMU_TIMEOUT:-1800}" "$QEMU" -M mps2-an386 -cpu cortex-m4 -nographic -monitor none -serial none \
            -semihosting-config "enable=on,target=native,$args,arg=$op" \
            -plugin "$SO$PLUGIN_ARGS" -d plugin -D "$tmp/$op.log" -kernel "$elf" > /dev/null &
        jobs+=($!)
    done
    for pid in "${jobs[@]}"; do wait "$pid"; done
    base_cycles=$(awk '/^cycles/{print $2}' "$tmp/none.log")
    base_insns=$(awk '/^insns/{print $2}' "$tmp/none.log")
    [[ -n "$base_cycles" && -n "$base_insns" ]] || return 1
    printf "%-14s %14s %14s %7s\n" operation cycles instructions CPI
    for op in "$@"; do
        [[ "$op" != none ]] || continue
        cycles=$(awk '/^cycles/{print $2}' "$tmp/$op.log")
        insns=$(awk '/^insns/{print $2}' "$tmp/$op.log")
        [[ -n "$cycles" && -n "$insns" ]] || { echo "missing result for $elf $op" >&2; return 1; }
        awk -v op="$op" -v c=$((cycles - base_cycles)) -v i=$((insns - base_insns)) \
            'BEGIN { printf "%-14s %14d %14d %7.3f\n", op, c, i, c / i }'
    done
}

for n in ${SIZES:-512 1024}; do
    for bk in ${BACKENDS:-ref armv7em m4f}; do
        dir=build/$bk-$n-cortex-m4
        make -s -C platform/cortex-m END_N="$n" BACKEND="$bk" "$dir/bench.elf" > /dev/null
        {
            echo "### Cortex-M4 model END-$n backend=$bk"
            measure_ops "platform/cortex-m/$dir/bench.elf" 'arg=bench,arg=1,arg=--insn' $OPS
            if [[ "$bk" == m4f ]]; then
                make -s -C platform/cortex-m END_N="$n" BACKEND="$bk" "$dir/bench_ones257.elf" > /dev/null
                echo "### Cortex-M4 standalone prefix product END-$n"
                measure_ops "platform/cortex-m/$dir/bench_ones257.elf" 'arg=bench_ones257' none ref dsp
            fi
        } | tee -a "$results"
    done
done
mv -- "$results" "$OUT"
