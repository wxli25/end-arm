#!/bin/bash
# Cortex-M4 cycle model + instruction counts for KEMs inside pqm4 (QEMU mps2-an386 + m4cycles plugin).
#
#   bash pqm4/m4model/run_m4model.sh <pqm4-dir> <out-dir> <impl-path>...
#        e.g. crypto_kem/end512/m4 mupq/pqclean/crypto_kem/ml-kem-512/clean
#
# Needs: pqm4/export_pqm4.py <pqm4> --icount (installs mupq/crypto_kem/icount_*.c) and a QEMU built
# with --enable-plugins.  For every implementation it builds elf/<name>_icount_<op>.elf for
# op in $OPS and writes <out-dir>/<name>_<op>.log (plugin output with a per-function profile).
# m4model_report.py turns the logs into tables.  Environment:
#   QEMU          qemu-system-arm with plugin support   (default /opt/qemu-plugins/bin/qemu-system-arm)
#   QEMU_INC      directory holding qemu-plugin.h       (default /opt/qemu-src/qemu-9.1.3/include/qemu)
#   PLUGIN_ARGS   extra model options, e.g. ",pdir=1,mla=1"
#   OPS           default "none kg10 keypair enc dec";  JOBS  parallel QEMU runs (default nproc)
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
P=$(cd "${1:?usage: run_m4model.sh <pqm4-dir> <out-dir> <impl-path>...}" && pwd)
OUT=${2:?usage: run_m4model.sh <pqm4-dir> <out-dir> <impl-path>...}
shift 2
mkdir -p "$OUT"
OUT=$(cd "$OUT" && pwd)
QEMU=${QEMU:-/opt/qemu-plugins/bin/qemu-system-arm}
QEMU_INC=${QEMU_INC:-/opt/qemu-src/qemu-9.1.3/include/qemu}
OPS=${OPS:-none kg10 keypair enc dec}
JOBS=${JOBS:-$(nproc)}
export QEMU PLUGIN_ARGS=${PLUGIN_ARGS:-} QEMU_TIMEOUT=${QEMU_TIMEOUT:-7200}
export SO="$OUT/m4cycles.so"

gcc -O2 -shared -fPIC -Wall -I"$QEMU_INC" $(pkg-config --cflags glib-2.0) -o "$SO" "$HERE/m4cycles.c"

cd "$P"
jobs_file="$OUT/jobs.txt"
: > "$jobs_file"
for impl in "$@"; do
    name=${impl//\//_}
    for op in $OPS; do
        make -s PLATFORM=mps2-an386 IMPLEMENTATION_PATH="$impl" "elf/${name}_icount_${op}.elf" \
            > /dev/null 2> "$OUT/${name}_build.err" || { cat "$OUT/${name}_build.err" >&2; exit 1; }
        echo "elf/${name}_icount_${op}.elf $OUT/${name}_${op}.log" >> "$jobs_file"
    done
done

run_one() {  # run_one <elf> <log>
    timeout "$QEMU_TIMEOUT" "$QEMU" -M mps2-an386 -nographic -semihosting \
        -plugin "$SO,profile=on$PLUGIN_ARGS" -d plugin -D "$2" -kernel "$1" < /dev/null > "$2.stdout" 2>&1 || true
    if grep -q "icount OK" "$2.stdout"; then echo "ok   $1"; else echo "FAIL $1" >&2; return 1; fi
}
export -f run_one
xargs -P "$JOBS" -L 1 bash -c 'run_one "$0" "$1"' < "$jobs_file"
