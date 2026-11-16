#!/bin/bash
# Host x86-64 baseline for the paper's Table 5: on one machine and one pinned core, time
#   - the paper's reference END code (END-{512,1024}-257-C and -AVX2, its own main.c and Makefile),
#   - ML-KEM from pq-crystals/kyber (ref and avx2, test/test_speed*),
#   - this repository's portable C implementation (bench/speed_host.c).
# Each is also built a second way so that compiler auto-vectorisation can be separated out: the
# paper builds END with -march=native and Kyber's ref with plain -O3.
#
#   bash scripts/host_baseline.sh <work-dir> [cpu]
#
# REF       reference implementation directory "END KEM" (default ../third_party/Explore_NTRU_Decryption/END KEM)
# NTEST_REF iterations for the reference main.c (default 2000; the paper uses 100000, which
#           takes hours for END-1024 KeyGen) -- the copy in <work-dir> is patched, not the original
# cpu       core to pin to (default: the fastest core found by a short calibration loop)
# The cycle counter is RDTSC (fixed frequency on current CPUs), as in the paper's and Kyber's code.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
W=${1:?usage: host_baseline.sh <work-dir> [cpu]}
REF=${REF:-$ROOT/../third_party/Explore_NTRU_Decryption/END KEM}
NTEST_REF=${NTEST_REF:-2000}
mkdir -p "$W"
W=$(cd "$W" && pwd)
cd "$W"

if [[ $# -ge 2 ]]; then
    CPU=$2
else  # pick the fastest core (hybrid CPUs: performance vs efficiency cores)
    cat > spin.c <<'EOF'
#include <stdio.h>
#include <stdint.h>
#include <x86intrin.h>
int main(void) { volatile uint64_t x = 1; uint64_t t0 = __rdtsc(); for (long i = 0; i < 200000000; i++) x = x * 3 + 1; printf("%llu\n", (unsigned long long)(__rdtsc() - t0)); return 0; }
EOF
    gcc -O2 spin.c -o spin
    best=; CPU=0
    for c in $(seq 0 $(($(nproc) - 1))); do
        t=$(taskset -c "$c" ./spin)
        if [[ -z "$best" || $t -lt $best ]]; then best=$t; CPU=$c; fi
    done
fi
PIN="taskset -c $CPU"
echo "# host: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //'), pinned to cpu $CPU"
echo "# compilers: $(clang --version | head -1); $(gcc --version | head -1)"

# 1. the paper's reference code
for n in 512 1024; do
    for v in C AVX2; do
        d=END-$n-$v
        rm -rf "$d"
        cp -r "$REF/END-$n/END-$n-257-$v" "$d"
        sed -i "s/^#define NTEST .*/#define NTEST $NTEST_REF/" "$d/main.c"
        make -s -C "$d" main > /dev/null 2>&1
        echo "## paper-ref END-$n $v (mean of $NTEST_REF)"
        $PIN "./$d/main" | grep -E 'cycles'
    done
    d=END-$n-C   # same sources without -march=native (no AVX2 auto-vectorisation)
    (cd "$d" && clang -O3 -DNDEBUG -fomit-frame-pointer -w -o main_nonative big_integer.c fixed_point.c keygen.c \
        code.c module_p.c mq_ntt.c ntru.c pke.c poly.c randombytes.c fips202.c kem.c main.c)
    echo "## paper-ref END-$n C without -march=native (mean of $NTEST_REF)"
    $PIN "./$d/main_nonative" | grep -E 'cycles'
done

# 2. ML-KEM (pq-crystals/kyber)
[[ -d kyber ]] || git clone -q https://github.com/pq-crystals/kyber.git
echo "# kyber commit $(git -C kyber log --oneline -1)"
rm -rf kyber/ref-native && cp -r kyber/ref kyber/ref-native
for v in ref ref-native avx2; do
    if [[ $v == ref-native ]]; then
        CFLAGS="-march=native -mtune=native" make -s -C "kyber/$v" speed > /dev/null 2>&1
    else
        make -s -C "kyber/$v" speed > /dev/null 2>&1
    fi
    for k in 512 768 1024; do
        echo "## ML-KEM-$k $v"
        $PIN "kyber/$v/test/test_speed$k" | grep -A2 -E '^kyber_(keypair|encaps|decaps): ' | grep -v '^--'
    done
done

# 3. this repository (portable C, clang as the paper's build), with and without -march=native
for n in 512 1024; do
    for flags in native nonative; do
        cf="-O3 -fomit-frame-pointer"
        [[ $flags == native ]] && cf="$cf -march=native -mtune=native"
        make -s -C "$ROOT" END_N=$n BACKEND=ref CC=clang CFLAGS="$cf" ARCH=hostbase-$flags \
            build/hostbase-$flags-ref-$n/speed_host > /dev/null
        echo "## end-arm END-$n ref (clang, $flags)"
        $PIN "$ROOT/build/hostbase-$flags-ref-$n/speed_host" 10000 200
    done
done
