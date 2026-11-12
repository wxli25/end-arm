# END KEM ARM 实现 · 代码阅读理解文档

本文档按"先整体后局部"的顺序介绍工程中的每个模块：它做什么、入口/出口（函数签名、输入输出及取值范围）、调用谁、被谁调用，以及阅读时需要留意的约定。配合 `DESIGN.md`（为什么这样设计）阅读。

---

## 0. 建议阅读顺序

1. `src/common/end_params.h` —— 所有常量与尺寸（10 分钟）
2. `golden/end_golden.py` —— 用 Python 读一遍整个方案：编解码、NTT 排列、解码器（30 分钟）
3. `src/common/end_kem.c` —— DPKE/PKE/KEM 与展开密钥（15 分钟）
4. `src/common/end_decode.c` —— 解码器主体（30 分钟）
5. `src/common/poly.h` + `src/ref/*.c` —— 内核契约与可移植实现（30 分钟）
6. `src/aarch64/*.c`、`src/armv7em/*.c` —— 两套 ARM 内核（各 30 分钟）
7. `src/common/end_keygen.c` + `src/common/ntrugen/` —— 密钥生成（求解器可当黑盒）
8. `test/`、`bench/`、`platform/cortex-m/`、`scripts/` —— 测试、基准、裸机平台与自动化

---

## 1. 目录总览

```
end-arm/
├── Makefile                  主机 / AArch64 Linux 构建（BACKEND=ref|aarch64, ARCH=host|aarch64, END_N=512|1024）
├── src/common/               与后端无关的方案层 + 共享头文件
├── src/ref/                  可移植 C 内核（基线）
├── src/aarch64/              NEON 内核
├── src/armv7em/              Cortex-M4/M33 DSP 内核
├── platform/cortex-m/        裸机：Makefile、link.ld、startup.c、randombytes_m.c
├── golden/                   Python golden model、wrap 注入测试（可导出 C 测试用例）
├── test/                     C 测试：内核自检、KAT、往返、wrap 注入、KeyGen scratch 测量
├── bench/                    基准程序与计时器
├── kat/                      参考实现 KAT（100 × END-512，60 × END-1024）与 wrap 实例
├── scripts/                  表生成、WSL/QEMU 辅助脚本、指令计数
└── docs/                     文档
```

### 1.1 三类符号前缀

- 公开 API：`end512_crypto_kem_*` / `end1024_crypto_kem_*`（`END_NS()` 宏生成）。
- 方案层内部函数：`pk_encode` 等短名在 `end_internal.h` 中 `#define` 为带前缀的符号。
- 内核：`poly_*`、`fcl_select` 同理在 `poly.h` 中加前缀。

---

## 2. `src/common/` 方案层

### 2.1 `end_params.h` —— 参数、尺寸、开关

- `END_N`（编译期）决定一切派生常量：`END_LOGN`、`END_LEAF`（4/8）、`END_DE`（64/90）、`END_DELTA`（1/0）、`END_MSGBYTES`（32/64）、`END_PKBYTES`、`END_CTBYTES`、`END_GAMMA2`（6/3）、`END_HASH1_BYTES`（32/64）。
- 模数：`END_Q=257`、`END_QP=12289`、`END_QQP`、`END_2QQP`、`END_P3=3329`。
- FCL：`END_M0_GF`、`END_M0_GF2`（默认 2），`END_NCAND_GF*`（候选数 9）。
- sk 布局偏移 `END_SK_OFF_*`；`END_KEYGEN_TMPBYTES = 22n+16`；`END_ALIGN`（16 字节对齐属性，所有多项式缓冲区都带）。

### 2.2 `end_api.h` —— 公开接口与展开密钥

```c
typedef struct END_ALIGN {
    int16_t h_ntt[N];        // pk 的 NTT-257 形式 [0,257)
    uint8_t k[MSGBYTES];     // 隐式拒绝密钥
    int16_t f16[N], F16[N];  // f（三元）、F（|.|≤31）作为 int16（旋转输入）
    int16_t w[N];
    int16_t fn257[N], Fn257[N], wfn257[N];   // NTT257(f), NTT257(F), NTT257(w·f)
    int16_t const1[N], const2[N], wconst[N]; // (g+δf)·𝟙, (G+δF)·𝟙, w(g+δf)·𝟙  (mod 257)
    int16_t wn12289[N];                      // NTT12289(w)
    int16_t fn3329[N], FQn3329[N];           // NTT3329(f), NTT3329(QF − wf)
    uint16_t p1, p2, pwc1;                   // 奇偶性：c1 的、G s̄+F μ̄ 的、w·c1 的
} end_expkey;
```
入口函数见 README。`crypto_kem_dec` 内部先 `expkey_from_sk` 再 `crypto_kem_dec_exp`；在 Cortex-M 上 `END_STATIC_EXPKEY` 让展开密钥成为静态单实例。

### 2.3 `end_internal.h` —— 内部接口

声明编解码（`pk_encode/decode`, `ct_encode/decode`, `s_from_mr`, `mr_from_s`）、压缩/解压内联函数、`end_hash1`（SHA3-256/512）、`decode_trap`、`expkey_from_parts`、`dpke_encrypt`、`pke_encrypt`、`pke_decrypt`、`ternary_sample`。

### 2.4 `end_arith.h` —— L1 标量模算术与常数时间原语

每个函数的输入/输出区间写在文件头注释中（也见 DESIGN §6 L1）。重点：
- `mont_reduce16(a, p, pinv_neg)`：通用 Montgomery（R=2^16），低 16 位用无符号乘法避免 UB。
- `red257(x)`：任意 int16 → [0,257)，Fermat 移位 + 两次修正。
- `round_div_2qqp(v)`、`center_2qqp(v)`：解码器用的舍入与中心化。
- `ct_lt_u32/u64`、`ct_mask_lt_u32`、`ct_is_zero_u32`、`ct_sel_i32`。

### 2.5 `end_tables.h`（生成文件）

由 `scripts/gen_tables.py` 生成：`end_zetas257[128]`、`end_zetas257_inv[128]`（及 `_mont` 形式）、`end_zetas3329*`、`end_zetas12289_{512,1024}`、`end_zetas12289_inv_*`（GS 消费顺序）、`END_INV128_*`、`END_NINV_MONT12289_*`、`END_MONT_R2_*`、`END_QINV16_*`、`END_INV257_MOD12289`、`END_INV257_MONT_MOD12289`、`END_INV_QQP_MONT_MOD3329`。

### 2.6 `end_codec.c` —— 编解码

| 函数 | 输入 | 输出 |
|---|---|---|
| `pk_encode(out[PK], h[N])` | h：NTT-257 域、任意 int16（按 mod 258 取） | 参考格式字节串 |
| `pk_decode(h[N], in[PK])` | 任意字节串 | h ∈ [0,257)（257 折到 0） |
| `ct_encode(out[CT], c[N])` | c ∈ [0, DE) | 字节串 |
| `ct_decode(c[N], in[CT])` | 任意字节串 | c ∈ [0, DE) |
| `s_from_mr(s, m′, r)` / `mr_from_s(m′, r, s)` | 字节 / 二元向量 | 二元向量 / 字节 |
| `end_compress(x)` / `end_decompress(c)`（内联） | x∈[0,256] / c∈[0,DE) | c / x′∈[0,256] |

内部：`ct_divmod_small_u64(x, d, ⌊2^32/d⌋)` 16 位 limb 倒数除法（常数时间、无 128 位）。格式细节见 DESIGN §8。

### 2.7 `end_keygen.c` —— 密钥生成

- `ternary_sample(a, state)`：512 版从 SHAKE256 连续读取 384 字节，每个 6 位组映射为一个三元系数。1024 版沿用参考实现的 5 个 SHAKE256 块和位运算构造。
- `compute_w(logn, f, g, F, G, tmp)`：参考实现原样（RNS NTT mod 31 位素数 + 32.32 定点 FFT 除法），输出 `w·20478` 于 `tmp[0..n)`；返回 0 表示溢出 → 拒绝密钥。
- `crypto_kem_keypair_derand(pk, sk, seed, tmp, tmplen)`：流程见 DESIGN §6 L4；`tmp` 需 ≥ `END_KEYGEN_TMPBYTES`。
- `crypto_kem_keypair(pk, sk)`：`randombytes(seed)` + 上式；`END_STATIC_KEYGEN_TMP` 时 scratch 为静态。
- 调试：`-DEND_KEYGEN_TRACE` 打印各阶段。

### 2.8 `end_kem.c` —— DPKE / PKE / KEM / 展开密钥

- `expkey_from_parts(ek, k, h_ntt, f, g, F, G, w)`：计算 2.2 中所有预计算量（8 次 NTT-257 前向等）。
- `expkey_from_sk(ek, sk)`：解析紧凑 sk 后调用上式。
- `dpke_encrypt(c, h_ntt, s)`：c = compress(INTT(NTT(h)∘NTT(s)))。
- `pke_encrypt(ct, h_ntt, m, r)`：ACWC 封装。
- `pke_decrypt(m, ek, ct)`：解码 + `mr_from_s` + 去掉 H1(r)；返回 `decode_trap` 的 found 标志。
- `crypto_kem_enc_derand / enc / dec_exp / dec`：FO 变换；`ct_memneq` 常数时间比较；失败时 `m ← k` 掩码选择。

### 2.9 `end_decode.c` —— DecodeTrap（核心）

```c
int decode_trap(int16_t sbar[N], const end_expkey *ek, const int16_t cprime[N], uint32_t diag[8]);
```
- 输入：`cprime` ∈ [0,256]（解压后的密文）。
- 输出：`sbar`（成功时每项 ±1），返回 1=GF 阶段某候选通过二元检查，0=失败（上层隐式拒绝）。`diag`（可为 NULL）：[0] gf 候选序号、[1] GF 候选序号、[2..3] gf 的 FCL 下标、[4..5] GF 的 FCL 下标、[6] found。
- 内部步骤与内核调用见 DESIGN §6 L6。关键局部函数：`lift_center_parity(r, parity)`（[0,257) → [−257,256]）、`residue_lift(...)`、`enumerate(cand, m0)`（候选族，顺序与 golden 一致）。
- 栈用量：约 14 个 `int16_t[N]` + 3 个 `int32_t[N]`（512：~20 KB；1024：~40 KB）。

### 2.10 `fips202.[ch]`、`randombytes.c`

pq-crystals Keccak；主机随机数（`/dev/urandom` / `CryptGenRandom`）。裸机用 `platform/cortex-m/randombytes_m.c`。

### 2.11 `ntrugen/`

Thomas Pornin 的 ntrugen（经参考实现修改）：`ntru.c`（`solve_NTRU`、`recover_G`）、`module_p.c`（31 位素数 RNS NTT）、`big_integer.c`（`zint_*`）、`fixed_point.c`（`fxr` 32.32 定点与 FFT）、`poly_big.c`（从参考 `poly.c` 抽出的大多项式辅助）。`param.h`/`ntrugen_poly.h` 为本工程的适配头。调用契约：`solve_NTRU(prof, logn, g, f, tmp)` 输出 `tmp` 前 2n 字节 = G‖F（int8），满足 gF − fG = q；scratch 22n 字节。

---

## 3. `src/common/poly.h` —— 内核接口（L2）

契约见 DESIGN §6 L2（输入/输出范围）。三套后端必须实现同一语义；NTT-257 的数据排列固定为参考排列（公钥线格式）；基乘法的叶子 blk 的模为 x^L ∓ `end_zetas257[64 + blk/2]`（偶 blk 取 +，奇 blk 取 −）。

---

## 4. `src/ref/` 可移植内核

| 文件 | 函数 | 要点 |
|---|---|---|
| `ntt257.c` | `poly_ntt257/intt257/basemul257` | 前向：输入 red257，7 层 CT 懒加法（<2056）；逆向：GS，和每层 red257；基乘法 schoolbook + `red257_32` |
| `baseinv257.c` | `poly_baseinv257` | 范数下降求逆，`inv257` = a^255；返回非零表示不可逆 |
| `ntt3329.c` | `poly_ntt3329/intt3329/basemul3329` | Montgomery 乘 twiddle（`_mont` 表）；基乘法 hi 先 `mont3329` 再 ×ζR²（`fqmul3329(ζ_mont, R2)`） |
| `ntt12289.c` | `poly_ntt12289/intt12289/pointwise12289` | [0,q) 条件加减；逐点乘 = 两次 Montgomery（第二次乘 R²） |
| `polyvec.c` | `poly_rot16/sqnorm16/cmov16/cmov32/fcl_select` | 旋转：log₂n 级桶形（两个临时数组）；FCL：有序表掩码插入 |

---

## 5. `src/aarch64/` NEON 内核

- `neon_common.h`：`neon_montmul`（4 条指令 Montgomery）、`neon_red257`、`neon_barrett`、`neon_center_exact`、`neon_canon`、`neon_montred32`（int32x4 → int16x4）、64/32 位 lane 的 zip/uzp 包装。常量 `PINV16_*`、`BARRETT_V_*`。
- `ntt257.c`：`ct_layer`（len ≥ 8 的整向量层）；len = 4 层用 `neon_zip64_lo/hi` 拆合；基乘法：`TRANSPOSE4`/`UNTRANSPOSE4` 宏（UZP/ZIP 两轮）把 8 个 4 次叶子转成系数主序，`vmull/vmlal` 累加，`red257_from32` 约简；1024 用 `transpose8`（TRN 三级）。
- `ntt3329.c`：同结构 + `leafcoef3329`（hi 先 `neon_montred32`，再 `vmlal(h, ζR)`，再 Montgomery + R² 修正，Barrett 规范化）。
- `ntt12289.c`：`ct_tail`/`gs_tail` 处理 len = 4/2/1（64/32/16 位 lane 的 ZIP/UZP），twiddle 向量按层构造；`reduce_all` 每两层。
- `polyvec.c`：`rot_stage_small`（`vextq` 跨向量移位 + `vbslq`）与整向量级；`sqnorm` 每 4 积宽化到 64 位；`fcl_select` 与 ref 相同。

NEON 中所有数据排列与 ref 相同（`test_kernels` 的 NTT-257 指纹 `93b98e37`/`e9b280d7` 与 golden 一致）。

---

## 6. `src/armv7em/` Cortex-M DSP 内核

- `dsp.h`：`dsp_smulbb/bt/tb/tt`、`dsp_smlabb/..`、`dsp_smuad/smuadx/smlad/smladx`、`dsp_smlald`、`dsp_sadd16/ssub16`、`dsp_pkhbt/pkhtb`（内联汇编），`dsp_montmul_pair(y, z, zq, −p)`（两系数 Montgomery 乘 twiddle）。
- `ntt257.c`：以 `int32_t *w = (int32_t *)a` 视图处理（需 `END_ALIGN`）；`len` 以系数计、`hl = len/2` 以字计；基乘法（L=4）用 SMUAD/SMUADX 组合。
- `ntt3329.c`：同结构；逆向每 3 层一次 Barrett（`barrett_pair`）。
- `ntt12289.c`：len ≥ 2 打包；len = 1 在字内做（`SMULTB` 取高半字）；逆向每层约简。
- `polyvec.c`：旋转第 0 级 `PKHBT` 重打包，其余级整字搬移，回绕字 `SSUB16(0, x)`；`sqnorm` 用 `SMLALD`。

---

## 7. `platform/cortex-m/`

- `link.ld`：FLASH 0x00000000/4 MB、RAM 0x20000000/4 MB（QEMU mps2）；可用 `-Wl,--defsym` 覆盖；`STACK_SIZE` 默认 128 KB；提供 `_sidata/_sdata/_edata/_sbss/_ebss/_estack/_sstack/_heap_limit/end`。
- `startup.c`：向量表（`.vectors`）、`Reset_Handler`（拷 .data、清 .bss、开 FPU、栈涂色、`initialise_monitor_handles`、半主机 `SYS_GET_CMDLINE` 取 argv、`main`、`exit`）、`platform_stack_used()`、`_sbrk`。
- `randombytes_m.c`：QEMU 配置启用测试用 DRBG（SHAKE256(seed‖ctr)）。板级构建调用 `platform_entropy()`，熵源失败时随机版 KEM API 返回 `-1`。
- `Makefile`：`END_N`、`BACKEND`（ref/armv7em）、`CPU`（cortex-m4/cortex-m33 → 自动选 FPU 标志与 QEMU 机型）、`OPT`；目标 `all` 生成 `test_kernels/test_kat/test_roundtrip/test_wraps/bench` 五个 ELF；`run-<name>` 在 QEMU 下运行。

---

## 8. `golden/`

- `end_golden.py`：精确算术规范。主要对象：`Params`（P512/P1024）、`pmul`（负循环乘，numpy int64 或纯 Python）、`ntt257/intt257/ntt257_mul/ntt257_inv`（参考排列）、`encode_pk/decode_pk/encode_ct/decode_ct`、`expanded_key`、`bnp_initial`（步骤 1–4）、`decode_trap`（含 FCL、候选、两阶段，返回 diag）、`kem_encaps/kem_decaps`、`parse_kat/check_entry`（对 KAT 逐项验证：gF−fG=q、h=g/f、pk 往返、ct/ss/解封装）。
- `test_wraps.py`：构造 wrap 实例（对齐 μ̄ 与 f、对齐 s̄ 与 G_Q），验证 V5 恒等式与解码；`--dump` 导出到 `kat/wraps_end*.txt` 供 `test/test_wraps.c` 回放；包含 S1/S6 的反例。

---

## 9. `test/`、`bench/`、`scripts/`

| 文件 | 用途 | 运行 |
|---|---|---|
| `test/test_kernels.c` | L1 穷举/随机 + L2 往返/教科书比对/求逆/旋转/FCL（指纹与 golden 比对） | `make test_kernels` |
| `test/test_kat.c` | KAT 回放：keypair_derand/enc_derand/dec/dec_exp 逐位比对 | `make test_kat` |
| `test/test_roundtrip.c` | 随机往返 + 篡改密文的隐式拒绝（K = H(k‖ct′)） | `make test_roundtrip` |
| `test/test_wraps.c` | 回放 golden 导出的 wrap 实例 | `make test_wraps` |
| `test/test_keygen_tmp.c` | KeyGen scratch 用量 canary | 手动编译 |
| `bench/bench.c` | 各操作与内核的中位数周期/时间；`--insn <op>` 单次模式供 QEMU 插件计数 | `make bench` |
| `bench/cycles.[ch]` | 平台计时器 | |
| `scripts/gen_tables.py` | 生成 `end_tables.h` | |
| `scripts/wsl_make.sh` | Windows → WSL 调用 make | |
| `scripts/qemu_m_tests.sh` / `qemu_m_run.sh` / `m_stack.sh` | QEMU 上运行 Cortex-M 测试/单个镜像/栈用量 | |
| `scripts/qemu_icount.sh` / `icount_all.sh` | 插件版 QEMU 的动态指令数统计 | |

---

## 10. 构建矩阵与一致性

| 平台 | 命令 | 验证内容 |
|---|---|---|
| 主机 ref | `make END_N=512 BACKEND=ref CC=gcc test_all` | 全部测试 |
| AArch64 ref/NEON | `make ... ARCH=aarch64 CC=aarch64-linux-gnu-gcc RUN="qemu-aarch64 -L /usr/aarch64-linux-gnu" test_all` | 全部测试 |
| Cortex-M4 ref/armv7em | `cd platform/cortex-m && make END_N=512 BACKEND=armv7em all`；`bash scripts/qemu_m_tests.sh 512 armv7em` | 全部测试 |
| Cortex-M33 | `make CPU=cortex-m33 ...`；`qemu_m_run.sh ... cortex-m33 mps2-an505` | 抽检 |

三个后端对同一 KAT 文件给出逐位相同的 pk/sk/ct/ss，即为跨后端一致性测试；NTT 域中间值因排列相同可直接 memcmp（`test_kernels` 的指纹即为示例）。
