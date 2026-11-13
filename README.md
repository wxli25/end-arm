# END KEM on ARM (Cortex-A / Cortex-M)

**English** | [简体中文](README.zh-CN.md)

This project is an ARM implementation of the **END KEM** (END-512 / END-1024) from the paper *How Compact Can NTRU Encryption Be? Heuristic Frontiers and Practical Schemes* (Liu, Zhang, Lu, Cheng, Yin). It contains three kernel backends and one shared scheme layer:

| Directory | Contents |
|---|---|
| `src/common/` | Backend-independent scheme layer: parameters, encoding/decoding, key generation (including the ntrugen solver), DecodeTrap decoding, ACWC/FO transform, API |
| `src/ref/` | Portable C kernels (baseline, correctness reference, fallback path) |
| `src/aarch64/` | Cortex-A: AArch64 NEON intrinsics kernels |
| `src/armv7em/` | Cortex-M4/M33: ARMv7E-M DSP kernels (packed 16-bit SIMD + SMULxx Montgomery) |
| `platform/cortex-m/` | Bare-metal startup code, linker scripts, semihosting, QEMU run scripts |
| `golden/` | Exact-arithmetic Python golden model (executable specification) and wrap-injection tests |
| `test/`, `bench/`, `kat/` | Test programs, benchmarks, KAT vectors |
| `docs/` | Design, code guide, test report, performance report (in Chinese) |

END-512 KeyGen samples from a continuous SHAKE256 stream. Its deterministic regression vectors are in `kat/kat_end512_shake256.txt`. The original END-512 vectors are still used to check Encaps and Decaps on their existing keys. END-1024 uses the KeyGen, Encaps and Decaps vectors of the reference implementation (commit `dde0a77`). See the [test report](docs/TEST_REPORT.md) for details.

## Quick start

All build and test commands run on Linux (WSL Ubuntu 24.04 recommended). Dependencies: `gcc`, `make`, `python3` (+numpy), `gcc-aarch64-linux-gnu`, `qemu-user`, `gcc-arm-none-eabi`, `qemu-system-arm`.

```bash
# 1. Portable C backend (host)
make END_N=512 BACKEND=ref CC=gcc test_all
make END_N=1024 BACKEND=ref CC=gcc test_all

# 2. Cortex-A: NEON backend, cross-compiled, tests run under qemu-aarch64
make END_N=512 BACKEND=aarch64 ARCH=aarch64 CC=aarch64-linux-gnu-gcc CFLAGS=-O3 \
     RUN="qemu-aarch64 -L /usr/aarch64-linux-gnu" test_all

# 3. Cortex-M4: bare-metal build + QEMU (mps2-an386)
cd platform/cortex-m && make END_N=512 BACKEND=armv7em all && cd ../..
bash scripts/qemu_m_tests.sh 512 armv7em

# 4. Golden model
python3 golden/end_golden.py kat/kat_end512_shake256.txt --n 512 --max 100
python3 golden/test_wraps.py
```

On Windows, call the scripts in `scripts/` through WSL, for example `MSYS_NO_PATHCONV=1 wsl.exe -d Ubuntu-24.04 -- bash scripts/wsl_make.sh 512 ref CC=gcc test_all`.

## API (`src/common/end_api.h`)

```c
int end512_crypto_kem_keypair(uint8_t *pk, uint8_t *sk);                 /* needs randombytes() */
int end512_crypto_kem_keypair_derand(uint8_t *pk, uint8_t *sk, const uint8_t seed[32], void *tmp, size_t tmplen);
int end512_crypto_kem_enc(uint8_t *ct, uint8_t *ss, const uint8_t *pk);
int end512_crypto_kem_enc_derand(uint8_t *ct, uint8_t *ss, const uint8_t *pk, const uint8_t m[32]);
int end512_crypto_kem_dec(uint8_t *ss, const uint8_t *ct, const uint8_t *sk);
int end512_expkey_from_sk(end_expkey *ek, const uint8_t *sk);           /* pre-expanded secret key (saves ~15 NTTs) */
int end512_crypto_kem_dec_exp(uint8_t *ss, const uint8_t *ct, const end_expkey *ek);
```

The END-1024 functions use the `end1024_` prefix.

Sizes (bytes): END-512 pk 514 / ct 384 / sk 3618 / ss 32; END-1024 pk 1027 / ct 832 / sk 7235 / ss 32.

The randomized KeyGen and Encaps return `-1` if the entropy source fails, and leave the output buffers unchanged. An external `randombytes` must return `int`: `0` after filling the whole buffer, `-1` on failure. The Cortex-M QEMU configuration uses a fixed-seed test source. `REAL_BOARD=1` builds into separate `*-board` directories. There, the randomized APIs need a board-provided `platform_entropy` that returns `1` after filling the whole buffer. Board startup and linker configuration must be adapted to the target board.

## Documentation

The documents below are in Chinese.

- `docs/DESIGN.md`: design document (decisions D0–D4, mathematical justification V1–V7, module and kernel design, value ranges, constant-time behavior, DR records)
- `docs/CODE_GUIDE.md`: code reading guide (function, inputs/outputs and call relations of each module, suggested reading order)
- `docs/TEST_REPORT.md`: test and acceptance report
- `docs/PERFORMANCE.md`: performance numbers (host cycle counts, QEMU dynamic instruction counts, stack and code size, comparison with ML-KEM in the pqm4 framework)
- `pqm4/README.md`: how to integrate with the pqm4 benchmarking framework (export, test, benchmark, real-board measurement)

## License and provenance

`src/common/ntrugen/` comes from Thomas Pornin's ntrugen (MIT, as modified by the paper's reference implementation). `src/common/fips202.c` comes from pq-crystals (CC0). The encoding format in `src/common/end_codec.c` matches the reference implementation. All other code was written for this project.
