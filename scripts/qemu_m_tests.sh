#!/bin/bash
# Run the Cortex-M test images under QEMU (semihosting).  usage: qemu_m_tests.sh <END_N> <BACKEND> [CPU] [MACHINE]
set -euo pipefail
cd "$(dirname "$0")/.." || exit 1
N=${1:-512}; BK=${2:-ref}; CPU=${3:-cortex-m4}; MACH=${4:-mps2-an386}
B=platform/cortex-m/build/$BK-$N-$CPU
run() {
    local name=$1
    shift
    echo "=== END-$N $BK $CPU: $name ==="
    bash scripts/qemu_m_run.sh "$B/$name.elf" "$CPU" "$MACH" "$name" "$@" |
        sed '/^[[:space:]]*\[PASS\]/d'
}
kat="kat/kat_end${N}_ref.txt"
run test_kernels
if [[ "$N" == 512 ]]; then
    run test_kat "$kat" "${KAT_MAX:-100}"
else
    run test_kat "$kat" "${KAT_MAX:-60}"
fi
run test_roundtrip "${RT_N:-10}"
run test_wraps "$kat" "kat/wraps_end${N}.txt"
