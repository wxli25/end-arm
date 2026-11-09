#!/bin/bash
# stack high-water marks and code sizes of the Cortex-M4 bench images (QEMU; DWT cycle counts are 0 under QEMU)
cd "$(dirname "$0")/.." || exit 1
for cfg in ref-512 armv7em-512 ref-1024 armv7em-1024; do
    elf=platform/cortex-m/build/$cfg-cortex-m4/bench.elf
    [ -f "$elf" ] || { echo "missing $elf"; continue; }
    echo "=== Cortex-M4 $cfg ==="
    bash scripts/qemu_m_run.sh "$elf" cortex-m4 mps2-an386 bench 2 | grep -a -i -E "stack|mismatch"
    arm-none-eabi-size "$elf" | tail -1
done
