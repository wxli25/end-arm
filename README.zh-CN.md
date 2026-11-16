# END KEM on ARM (Cortex-A / Cortex-M)

[English](README.md) | **简体中文**

本工程是论文 *How Compact Can NTRU Encryption Be? Heuristic Frontiers and Practical Schemes*（Liu, Zhang, Lu, Cheng, Yin）中 **END KEM**（END-512 / END-1024）的 ARM 平台实现，包含三套内核后端和一个共享的方案层：

| 目录 | 内容 |
|---|---|
| `src/common/` | 与后端无关的方案层：参数、编解码、密钥生成（含 ntrugen 求解器）、DecodeTrap 解码、ACWC/FO 变换、API |
| `src/ref/` | 可移植 C 内核（基线、正确性参照、回退路径） |
| `src/aarch64/` | Cortex-A：AArch64 NEON intrinsics 内核 |
| `src/armv7em/` | Cortex-M4/M33：C + DSP 内核（`BACKEND=armv7em`）和手写 Thumb-2 汇编内核（`BACKEND=m4f`：Plantard 算术、NTT 层合并、FPU 寄存器作暂存） |
| `platform/cortex-m/` | 裸机启动代码、链接脚本、半主机、QEMU 运行脚本 |
| `golden/` | Python 精确算术 golden model（可执行规范）与 wrap 注入测试 |
| `test/`, `bench/` | 测试程序、基准程序 |

END-512 的 KeyGen 使用连续 SHAKE256 流，密钥因此与参考实现不同；封装和解封装仍用参考实现的向量在已有密钥上验证。END-1024 复现参考实现 `dde0a77` 的 KeyGen、封装和解封装向量。

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
cd platform/cortex-m && make END_N=512 BACKEND=m4f all && cd ../..      # 汇编内核；BACKEND=armv7em 为 C + DSP
bash scripts/qemu_m_tests.sh 512 m4f

# 4. golden model（需要测试向量，见下文）
python3 golden/end_golden.py kat/kat_end512_shake256.txt --n 512 --max 100
python3 golden/test_wraps.py
```

Windows 下在 WSL 里执行这些命令。

仓库只存放源码：

- `src/common/end_tables.h` 和 `src/armv7em/end_tables_plant.h`（twiddle 等常数表）在首次构建时由 `scripts/gen_tables.py` 生成。
- 测试向量不随仓库发布。`test_kat`、`test_wraps` 和 golden model 的命令从 `kat/` 读取向量；某个向量文件不存在时，`make test_all` 和 `scripts/qemu_m_tests.sh` 跳过需要它的回放并给出提示。向量文件有三类：`kat_end{512,1024}_ref.txt`（由参考实现 `dde0a77` 生成）、`kat_end512_shake256.txt`（`test_kat kat/kat_end512_ref.txt 100 --generate`）、`wraps_end{512,1024}.txt`（`python3 golden/test_wraps.py --dump`）。
- 设计文档、测试报告和性能报告不随仓库发布。

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

`pqm4/README.md` 说明接入 pqm4 基准框架的步骤：导出、测试、基准、真板测量。

## 许可与来源

`src/common/ntrugen/` 来自 Thomas Pornin 的 ntrugen（MIT，经论文参考实现修改）；`src/common/fips202.c` 来自 pq-crystals（CC0）；`platform/cortex-m/keccakf1600.S` 是 Keccak 团队（XKCP）的 ARMv7-M Keccak-f[1600] 汇编，含 Alexandre Adomnicai 的优化，取自 pqm4（CC0）；`src/common/end_codec.c` 的编码格式与参考实现一致。其余代码为本工程编写，其中 Cortex-M4 汇编内核的 Plantard 运算序列和合并层 NTT 的安排参照 pqm4 中的 ML-KEM 代码。
