#!/bin/bash
# Consolidated dynamic-instruction-count table (QEMU + libinsn plugin) for all backends / parameter sets.
# Requires the plugin-enabled QEMU in /opt/qemu-plugins (see docs/PERFORMANCE.md) and the bench binaries:
#   AArch64 : build/aarch64-<backend>-<N>/bench          (make BACKEND=... ARCH=aarch64 CC=aarch64-linux-gnu-gcc ...)
#   Cortex-M: platform/cortex-m/build/<backend>-<N>-cortex-m4/bench.elf
set -euo pipefail
cd "$(dirname "$0")/.." || exit 1
OUT=${1:-docs/icount_results.txt}
results=$(mktemp "${OUT}.XXXXXX")
trap 'rm -f -- "$results"' EXIT
for n in 512 1024; do
    for bk in ref aarch64; do
        bin=build/aarch64-$bk-$n/bench
        if [ ! -x "$bin" ]; then
            make -s END_N=$n BACKEND=$bk ARCH=aarch64 CC=aarch64-linux-gnu-gcc CFLAGS="-O3" "$bin" > /dev/null 2>&1
        fi
        echo "### AArch64 END-$n backend=$bk (qemu-aarch64 user mode)" | tee -a "$results"
        bash scripts/qemu_icount.sh aarch64 "$bin" | tee -a "$results"
    done
    for bk in ref armv7em; do
        elf=platform/cortex-m/build/$bk-$n-cortex-m4/bench.elf
        [ -f "$elf" ] || continue
        echo "### Cortex-M4 END-$n backend=$bk (qemu-system-arm mps2-an386)" | tee -a "$results"
        bash scripts/qemu_icount.sh cortex-m "$elf" | tee -a "$results"
    done
done
mv -- "$results" "$OUT"
