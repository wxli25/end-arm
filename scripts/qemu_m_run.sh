#!/bin/bash
# Run one Cortex-M ELF under QEMU with semihosting.  usage: qemu_m_run.sh <elf> [cpu] [machine] [args...]
set -eu
cd "$(dirname "$0")/.." || exit 1
ELF=${1:?usage: qemu_m_run.sh <elf> [cpu] [machine] [args...]}
CPU=${2:-cortex-m4}
MACH=${3:-mps2-an386}
shift "$(( $# < 3 ? $# : 3 ))"
args=""
for a in "$@"; do args="$args,arg=$a"; done
timeout "${QEMU_TIMEOUT:-1800}" qemu-system-arm -M "$MACH" -cpu "$CPU" -nographic -monitor none -serial none \
    -semihosting-config "enable=on,target=native$args" -kernel "$ELF" 2>&1
