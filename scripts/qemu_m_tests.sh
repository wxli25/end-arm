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
# The test vectors (kat/) are not stored in the repository: a replay whose vector file is absent is skipped.
have() {
    local f
    for f in "$@"; do
        [[ -f "$f" ]] || { echo "=== END-$N $BK $CPU: SKIPPED, $f not found ==="; return 1; }
    done
}
kat="kat/kat_end${N}_ref.txt"
run test_kernels
if [[ "$N" == 512 ]]; then
    if have kat/kat_end512_shake256.txt; then run test_kat kat/kat_end512_shake256.txt "${KAT_MAX:-100}"; fi
    if have "$kat"; then run test_kat "$kat" "${KAT_MAX:-100}" --import; fi
else
    if have "$kat"; then run test_kat "$kat" "${KAT_MAX:-60}"; fi
fi
run test_roundtrip "${RT_N:-10}"
if have "$kat" "kat/wraps_end${N}.txt"; then run test_wraps "$kat" "kat/wraps_end${N}.txt"; fi
