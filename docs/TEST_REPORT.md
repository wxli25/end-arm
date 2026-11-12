# END KEM ARM 实现 · 测试验收报告

## 2026-11-12 接入 pqm4

- 新增 `pqm4/`：导出脚本、胶水代码、指令数 harness。KeyGen 的 SHAKE256 改经 `src/common/end_xof.h` 适配，pqm4 构建使用 pqm4 的 Keccak 汇编，本仓库构建仍用自带 `fips202.c`，输出逐字节相同。
- 修复：ARMv7E-M 内核的系数对读写改用 `__builtin_memcpy`。pqm4 平台带 `-ffreestanding`，普通 `memcpy` 在 DSP 内核里变成了函数调用，使 m4 版解封装慢到 7.9 M 条指令；修复后为 2.73 M。

| 验证 | 结果 |
|---|---|
| pqm4 `test.py -p mps2-an386`：end512/end1024 的 m4、ref | 4 个实现全部 SUCCESSFUL（含无效密文的隐式拒绝检查） |
| pqm4 `testvectors.py -p mps2-an386`：板上实现与主机 ref 逐字节比对 | 4 个实现全部 SUCCESSFUL |
| pqm4 `benchmarks.py -p mps2-an386`：speed、stack、hashing、size | 全部运行成功，数据见 PERFORMANCE §6 |
| 本仓库全套回归（主机、qemu-aarch64 NEON、QEMU Cortex-M4/M33、KAT 重新生成、golden model） | 见下一节同日记录，修改后重新运行并通过 |

## 2026-11-12 修复回归

| 环境 | 本轮结果 |
|---|---|
| WSL 主机、QEMU AArch64 ref/NEON、QEMU Cortex-M4 ref/DSP | END-512 通过 100 组新 KAT 和 100 组参考密钥导入回归。END-1024 通过 60 组参考 KAT。各配置的内核、往返、隐式拒绝及 9 个 wrap 用例均通过 |
| QEMU Cortex-M33 DSP | END-512 通过内核测试、新旧各 5 组 KAT、5 次往返及 9 个 wrap 用例 |
| WSL 主机扩展测试 | 两套参数各 2,000 次随机往返及篡改密文检查，零失败。FCL 两阶段取 4 和 8 时，两套参数的 KAT、往返及熵源测试通过 ASan/UBSan。取 0 或 9 时编译报错 |
| Python golden model | 检查 100 组新 END-512 向量、60 组 END-1024 向量和两套 wrap 注入用例 |
| 失败路径 | 缺失 ELF、KAT 数量不足、模拟器失败和基准文件缺失均返回非零。板级熵源模拟测试确认失败时 API 返回 `-1` 且输出缓冲区不变 |

Cortex-M4 板级配置编译通过，随机数对象中没有固定种子测试源。DSP NTT 的打包访存仍编译为内联加载和存储。真机和板级 TRNG 接入待验证。

使用原 END-512 向量的种子和消息重新生成回归向量：

```bash
make END_N=512 build/host-ref-512/test_kat
build/host-ref-512/test_kat kat/kat_end512_ref.txt 100 --generate > kat/kat_end512_shake256.txt
python3 golden/end_golden.py kat/kat_end512_shake256.txt --n 512 --max 100
```

## 2026-11-12 工程加固回归（提交 3f4f568）

该提交只含不改变计算结果的加固：`randombytes` 返回错误码；Cortex-M 测试 DRBG 只编进 QEMU 配置；ARMv7E-M 内核改用 `memcpy` 打包访存；`END_M0_GF*` 限定在 1–8；测试与脚本对失败路径返回非零。当时 END-512 KeyGen 仍与参考实现逐位一致，主机、qemu-aarch64 NEON、QEMU Cortex-M4（DSP/ref）、Cortex-M33 及 golden model 的 16 项回归全部通过。END-512 采样随后按上节修改。

## 2026-11-06 初版记录

日期：2026-11-06。环境：WSL Ubuntu 24.04（gcc 13.3、aarch64-linux-gnu-gcc 13.3、arm-none-eabi-gcc 13.2、QEMU 8.2 / 自建 9.1）。
参考实现：`Explore_NTRU_Decryption` @ `dde0a77`，用确定性 DRBG（`seed‖ctr` 的 SHAKE256）在主机上生成 KAT：`kat/kat_end512_ref.txt`（100 组）、`kat/kat_end1024_ref.txt`（60 组）；每组含 `keygen_seed, pk, f, g, F, G, w, k, m, ct, ss`。

## 1. 验收矩阵

| 测试 | 内容 | 主机 ref ×2 | AArch64 ref ×2 | AArch64 NEON ×2 | Cortex-M4 ref ×2 | Cortex-M4 DSP ×2 | Cortex-M33 DSP (512) |
|---|---|---|---|---|---|---|---|
| test_kernels | L1 约简 int16 全域穷举 + 2×10⁶ 随机 int32；NTT-257/3329/12289 往返；基乘法/逐点乘 vs 教科书乘法；NTT-257 求逆；旋转；范数；FCL；NTT-257 排列指纹（与 golden 一致：`93b98e37` / `e9b280d7`） | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| test_kat | keypair_derand → pk,k,f,g,F,G,w 逐位；enc_derand → ct,ss 逐位；dec / dec_exp → ss 逐位 | ✅ 100/60 | ✅ 100 (512) | ✅ 100/60 | ✅ 100/60 | ✅ 100/60 | ✅ 5 |
| test_roundtrip | 随机 keypair/enc/dec；篡改密文 → ss 不同且等于 SHA3-256(k‖ct′) | ✅ 20 | ✅ 20 | ✅ 20 | ✅ 10 | ✅ 10 | — |
| test_wraps | golden 导出的 9 个 wrap 实例（KAT 密钥 0）：无 wrap、单 gf-wrap ×3（坐标 0/17/n−1）、双 gf-wrap {3,200}、单 GF-wrap ×2（1.3qQ/1.15qQ）、gf+GF、3 个 gf-wrap（须失败） | ✅ 9/9 | ✅ | ✅ | ✅ | ✅ | — |
| golden（Python） | KAT 20/10 组全链路（含 gF−fG=q、h=g/f、pk 往返）；wrap 注入 13 项 ×2（含 V5 恒等式、S1、S6/V2） | ✅ | — | — | — | — | — |

"×2" 表示 END-512 与 END-1024 都通过。Cortex-M 上通过 QEMU（mps2-an386 / mps2-an505 + 半主机）运行；M33 做了 test_kernels 与 5 组 KAT 的抽检。

## 2. 关键用例说明

### 2.1 wrap 注入（`golden/test_wraps.py` → `kat/wraps_end*.txt` → `test/test_wraps.c`）

构造方法：取 KAT 密钥 0；随机 s；把 μ̄ 的取值在 f 的非零位置上与 (fμ̄)_t 的符号对齐，使 g s̄ + f μ̄ 在坐标 t 上略超过 q（gf-wrap）；或把 s̄ 与 G_Q 的大系数对齐使 G_Q s̄ + F_Q μ̄ 在坐标 0 上达到 1.15–1.3 qQ（GF-wrap）。golden 对每个实例验证 V5 恒等式（真实 e₁、e₂），再由 C 解码器回放：

```
[PASS] single gf wrap at 0 (gf_cand=2 GF_cand=0)      ← 选中候选 +x^0
[PASS] double gf wrap at 3,200 (gf_cand=8 GF_cand=0)   ← 选中候选 x^3 + x^200
[PASS] single GF wrap at 0 target 1.15 (gf_cand=0 GF_cand=2)
[PASS] gf wrap at 5 + GF wrap (gf_cand=2 GF_cand=2)
[PASS] 3 gf wraps beyond B0 (found=0 same=0)           ← 正确地报告失败 → 隐式拒绝
```

golden 中另有两组反例：用参考实现的 gf 阶段打分（`128 − |c₁|`，按绝对值排序）**选不到** wrap 坐标（S1）；用参考的修正符号（`s′ − rf − ẽF`）对真候选**得不到** s（S6）；而修正后的解码器恢复 s（V2）。

FCL 的前提是 wrap 坐标的未约简值只略超边界：实验中把 GF 目标设为 2.4 qQ 时（残差落在区间中部）FCL 不再命中，这是算法的固有假设，不是实现缺陷（见 DESIGN §14）。

### 2.2 对参考实现的实验

- END-1024 参考实现 400 次解封装失败 14 次（3.5%）；把 `poly_norm2` 的返回类型从 int16 改为 uint32 后 0 次失败（`int16` 溢出使候选范数比较出错）。END-512 参考实现 2000 次 0 失败（其 `poly_norm2` 为 uint32）。
- 2000 组 END-512 / 400 组 END-1024 密钥统计：Σf、Σg 100% 为奇数（ntrugen `zint_bezout` 需要奇结数），ΣF 与 ΣG 恰一奇一偶，Σw 奇偶各半——解释了参考实现 CRT 奇偶处理"恰好正确"的原因（DESIGN DR-04）。

### 2.3 KeyGen scratch

`test/test_keygen_tmp.c` 以 0xA5 涂色 64·n 字节后运行 keypair_derand：两种参数集、含基准种子在内的 12 个种子均恰好触及 22·n 字节（`END_KEYGEN_TMPBYTES = 22n + 16`）。

## 3. 构建警告与已知事项

- 测试程序在 32 位目标上的 `%d/%u` 格式告警已通过显式转换消除；库代码在 `-Wall -Wextra` 下无告警。
- Cortex-M 的 `randombytes_m.c` 为测试 DRBG，固定种子（真机须接 TRNG）。
- QEMU 不模拟周期，也不实现 DWT：Cortex-M 构建默认 `-DEND_NO_DWT`（`cycles_now()` 返回 0），真机用 `make REAL_BOARD=1` 启用 DWT CYCCNT；性能见 `PERFORMANCE.md` 的指令数。
- newlib-nano 的 `printf` 不支持 `%llu`/`%f`：基准程序只使用 `%lu` 输出（早期版本在 Cortex-M 上因此打印乱码，已修复）。
- Cortex-M33 在 QEMU `mps2-an505` 上需使用安全别名（代码 0x10000000、RAM 0x38000000），`platform/cortex-m/Makefile` 对 `CPU=cortex-m33` 已默认设置。
- 调试期间发现并修复的问题（均已被自动测试覆盖）：`red257` 对 x=−257 输出 257（穷举测试）；NTT-257 曾错误地依赖 int16 回绕（往返测试）；Montgomery 低 16 位乘法的有符号溢出 UB；`F_Q` 预计算的 int16 溢出（KAT 5/100 失败→0）；NEON Barrett 的 ±1 边界（往返测试）；NEON `sqnorm` 的 int32 累加溢出（全域随机测试）。

## 4. 复现

```bash
make END_N=512 BACKEND=ref CC=gcc test_all
make END_N=1024 BACKEND=ref CC=gcc test_all
make END_N=512 BACKEND=aarch64 ARCH=aarch64 CC=aarch64-linux-gnu-gcc CFLAGS=-O3 RUN="qemu-aarch64 -L /usr/aarch64-linux-gnu" test_all
make END_N=1024 BACKEND=aarch64 ARCH=aarch64 CC=aarch64-linux-gnu-gcc CFLAGS=-O3 RUN="qemu-aarch64 -L /usr/aarch64-linux-gnu" test_all
(cd platform/cortex-m && make END_N=512 BACKEND=armv7em all && make END_N=1024 BACKEND=armv7em all && make END_N=512 BACKEND=ref all && make END_N=1024 BACKEND=ref all)
bash scripts/qemu_m_tests.sh 512 armv7em; bash scripts/qemu_m_tests.sh 1024 armv7em
bash scripts/qemu_m_tests.sh 512 ref;     bash scripts/qemu_m_tests.sh 1024 ref
(cd platform/cortex-m && make END_N=512 BACKEND=armv7em CPU=cortex-m33 all) && KAT_MAX=5 RT_N=5 bash scripts/qemu_m_tests.sh 512 armv7em cortex-m33 mps2-an505
python3 golden/end_golden.py kat/kat_end512_shake256.txt --n 512 --max 100 && python3 golden/test_wraps.py
```
