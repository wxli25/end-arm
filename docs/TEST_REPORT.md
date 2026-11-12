# END KEM ARM 实现 · 测试验收报告

## 2026-11-12 工程加固回归

本轮改动不改变任何计算结果：`randombytes` 返回错误码，随机版 KeyGen/Encaps 在熵源失败时返回 `-1` 且不写输出；Cortex-M 的测试 DRBG 只编进 QEMU 配置，`REAL_BOARD=1` 使用板级熵源；ARMv7E-M 内核改用 `memcpy` 打包访存（消除别名规则违例）；`END_M0_GF*` 限定在 1–8 并在编译期检查；测试与脚本对失败路径返回非零。

| 配置 | 结果 |
|---|---|
| 主机 ref、qemu-aarch64 NEON（END-512 / END-1024） | 内核自检、参考 KAT（100 / 60 组，含 KeyGen 逐位比对）、20 次往返与篡改、9 个 wrap 用例、熵源失败/成功路径，全部通过 |
| QEMU Cortex-M4 armv7em 与 ref（END-512 / END-1024） | 内核自检、参考 KAT（100 / 60 组）、10 次往返、9 个 wrap 用例，全部通过 |
| QEMU Cortex-M33 armv7em（END-512） | 内核自检、5 组 KAT、5 次往返、9 个 wrap 用例，通过 |
| Cortex-M4 板级配置（`REAL_BOARD=1`） | 编译通过 |
| Python golden model | 参考 KAT 100 / 60 组、wrap 注入用例，通过 |

END-512 KeyGen 的 XOF 读取方式（参考实现会重复 48 个采样字节，见 DESIGN DR-13）是否修改，待与开发人员确认；本轮保持与参考实现兼容。

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
(cd platform/cortex-m && make END_N=512 BACKEND=armv7em CPU=cortex-m33 all) && bash scripts/qemu_m_run.sh platform/cortex-m/build/armv7em-512-cortex-m33/test_kat.elf cortex-m33 mps2-an505 x kat/kat_end512_ref.txt 5
python3 golden/end_golden.py kat/kat_end512_ref.txt --n 512 --max 20 && python3 golden/test_wraps.py
```
