#!/bin/bash
# Run the Cortex-M test images under QEMU (semihosting).  usage: qemu_m_tests.sh <END_N> <BACKEND> [CPU] [MACHINE]
cd "$(dirname "$0")/.." || exit 1
N=${1:-512}; BK=${2:-ref}; CPU=${3:-cortex-m4}; MACH=${4:-mps2-an386}
B=platform/cortex-m/build/$BK-$N-$CPU
Q="qemu-system-arm -M $MACH -cpu $CPU -nographic -monitor none -serial none"
run() { # run <elf> <args...>
    local elf=$1; shift
    local args=""
    for a in "$@"; do args="$args,arg=$a"; done
    timeout 1800 $Q -semihosting-config "enable=on,target=native$args" -kernel "$elf" 2>&1 | grep -v "^\s*\[PASS\]" | tail -3
}
echo "=== END-$N $BK $CPU: test_kernels ==="; run $B/test_kernels.elf
echo "=== test_kat ===";       run $B/test_kat.elf x kat/kat_end${N}_ref.txt ${KAT_MAX:-100}
echo "=== test_roundtrip ===";  run $B/test_roundtrip.elf x ${RT_N:-10}
echo "=== test_wraps ===";      run $B/test_wraps.elf x kat/kat_end${N}_ref.txt kat/wraps_end${N}.txt
