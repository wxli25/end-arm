# END KEM on ARM (Cortex-A / Cortex-M)

本工程是论文 *How Compact Can NTRU Encryption Be? Heuristic Frontiers and Practical Schemes*（Liu, Zhang, Lu, Cheng, Yin）中 **END KEM**（END-512 / END-1024）的 ARM 平台实现，包含三套内核后端和一个共享的方案层：

| 目录 | 内容 |
|---|---|
| `src/common/` | 与后端无关的方案层：参数、编解码、密钥生成（含 ntrugen 求解器）、DecodeTrap 解码、ACWC/FO 变换、API |
| `src/ref/` | 可移植 C 内核（基线、正确性参照、回退路径） |
| `src/aarch64/` | Cortex-A：AArch64 NEON intrinsics 内核 |
| `src/armv7em/` | Cortex-M4/M33：ARMv7E-M DSP（打包 16 位 SIMD + SMULxx Montgomery）内核 |
| `platform/cortex-m/` | 裸机启动代码、链接脚本、半主机、QEMU 运行脚本 |
| `golden/` | Python 精确算术 golden model（可执行规范）与 wrap 注入测试 |
| `test/`, `bench/`, `kat/` | 测试程序、基准程序、参考实现 KAT |
| `docs/` | 设计文档、代码阅读文档、测试报告、性能报告 |

END-512 的 KeyGen 使用连续 SHAKE256 流，确定性回归使用 `kat/kat_end512_shake256.txt`。原 END-512 向量继续验证已有密钥的封装和解封装。END-1024 仍使用参考实现 `dde0a77` 的 KeyGen、封装和解封装向量。详见 [测试报告](docs/TEST_REPORT.md)。

## 快速开始

所有构建/测试命令在 Linux（推荐 WSL Ubuntu 24.04）下执行。依赖：`gcc`、`make`、`python3`(+numpy)、`gcc-aarch64-linux-gnu`、`qemu-user`、`gcc-arm-none-eabi`、`qemu-system-arm`。

```bash
# 1. 可移植 C 后端（主机）
make END_N=512 BACKEND=ref CC=gcc test_all
make END_N=1024 BACKEND=ref CC=gcc test_all

# 2. Cortex-A：NEON 后端，交叉编译并在 qemu-aarch64 下运行测试
make END_N=512 BACKEND=aarch64 ARCH=aarch64 CC=aarch64-linux-gnu-gcc CFLAGS=-O3 \
     RUN="qemu-aarch64 -L /usr/aarch64-linux-gnu" test_all

# 3. Cortex-M4：裸机构建 + QEMU (mps2-an386) 运行
cd platform/cortex-m && make END_N=512 BACKEND=armv7em all && cd ../..
bash scripts/qemu_m_tests.sh 512 armv7em

# 4. golden model
python3 golden/end_golden.py kat/kat_end512_shake256.txt --n 512 --max 100
python3 golden/test_wraps.py
```

Windows 下可通过 `MSYS_NO_PATHCONV=1 wsl.exe -d Ubuntu-24.04 -- bash scripts/wsl_make.sh 512 ref CC=gcc test_all` 调用（见 `scripts/`）。

## API（`src/common/end_api.h`）

```c
int end512_crypto_kem_keypair(uint8_t *pk, uint8_t *sk);                 /* 需要 randombytes() */
int end512_crypto_kem_keypair_derand(uint8_t *pk, uint8_t *sk, const uint8_t seed[32], void *tmp, size_t tmplen);
int end512_crypto_kem_enc(uint8_t *ct, uint8_t *ss, const uint8_t *pk);
int end512_crypto_kem_enc_derand(uint8_t *ct, uint8_t *ss, const uint8_t *pk, const uint8_t m[32]);
int end512_crypto_kem_dec(uint8_t *ss, const uint8_t *ct, const uint8_t *sk);
int end512_expkey_from_sk(end_expkey *ek, const uint8_t *sk);           /* 预展开私钥（省 ~15 次 NTT） */
int end512_crypto_kem_dec_exp(uint8_t *ss, const uint8_t *ct, const end_expkey *ek);
```

尺寸（字节）：END-512 pk 514 / ct 384 / sk 3618 / ss 32；END-1024 pk 1027 / ct 832 / sk 7235 / ss 32。

随机版 KeyGen 和 Encaps 在熵源失败时返回 `-1`，输出缓冲区保持原值。外部 `randombytes` 实现须返回 `int`，成功填满缓冲区后返回 `0`，失败返回 `-1`。Cortex-M 的 QEMU 配置使用固定种子测试源。`REAL_BOARD=1` 使用独立的 `*-board` 构建目录，随机版 API 要求真机提供 `platform_entropy`，成功填满缓冲区后返回 `1`。板级启动和链接配置按目标板接入。

## 文档

- `docs/DESIGN.md`：设计文档（决策 D0–D4、数学依据 V1–V7、模块与内核设计、数值范围、常数时间、DR 记录）
- `docs/CODE_GUIDE.md`：代码阅读理解文档（每个模块的功能、入口/出口、调用关系、阅读顺序）
- `docs/TEST_REPORT.md`：测试验收报告
- `docs/PERFORMANCE.md`：性能统计（主机周期数、QEMU 动态指令数、栈/代码尺寸）

## 许可与来源

`src/common/ntrugen/` 来自 Thomas Pornin 的 ntrugen（MIT，经论文参考实现修改）；`src/common/fips202.c` 来自 pq-crystals（CC0）；`src/common/end_codec.c` 的编码格式与参考实现一致。其余代码为本工程编写。
