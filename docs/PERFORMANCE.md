# END KEM ARM 实现 · 性能统计

> 本轮没有 ARM 真机。给出的 ARM 数据是 **QEMU TCG 插件统计的动态指令数**（与硬件无关、可复现，能反映算法/内核层面的差异，但不等于周期数），以及 **x86-64 主机的周期数**（可移植 C 后端，用于校准数量级）。真机周期数的测量方法与脚本已就绪（§5）。
> 表中数据来自初版提交 `9caa833`。END-512 采样和 Cortex-M DSP 访存修改后，可用各节命令重新测量。
> 2026-11-12 的固定种子抽测中，END-512 KeyGen 在 AArch64 NEON 和 Cortex-M4 DSP 上分别执行约 44.3 M 和 68.2 M 条指令。KeyGen 重试次数随种子及采样流变化，这两个数值对应单个基准种子。

## 1. 结论

| 指标（END-512） | 可移植 C | 优化后端 | 提升 |
|---|---|---|---|
| AArch64 Decaps（展开密钥）动态指令数 | 839 k | **438 k**（NEON） | 1.92× |
| AArch64 Encaps 动态指令数 | 122 k | **92 k** | 1.33× |
| Cortex-M4 Decaps（展开密钥）动态指令数 | 3.02 M | **1.98 M**（DSP） | 1.52× |
| Cortex-M4 Encaps 动态指令数 | 380 k | **277 k** | 1.37× |
| x86-64 周期数（参考数量级） | KeyGen 3.67 M / Encaps 35 k / Decaps 158 k（展开）240 k（紧凑） | — | — |

论文 Table 5（i9-11900K，参考 C 实现）：END-512 KeyGen 10 408 k / Encaps 23 k / Decaps 267 k cycles。本实现主机 C 版本（不同机器、`-O3 -march=native`）的 Decaps 与之同量级，Encaps 略高（本实现在 Encaps 内多做一次 `pk_decode`；参考在 KEM 层之外解码）。

## 2. AArch64（qemu-aarch64 用户态，`-O3`）动态指令数

单位：指令。`dec` = 紧凑密钥（含展开），`dec_exp` = 展开密钥；`expkey` = `expkey_from_sk`。

| 操作 | END-512 ref | END-512 NEON | END-1024 ref | END-1024 NEON |
|---|---:|---:|---:|---:|
| keypair | 17 361 146 | 17 333 447 | 51 503 820 | 51 448 886 |
| enc | 122 205 | 92 147 | 241 316 | 184 630 |
| dec（紧凑） | 1 241 236 | 616 702 | 2 675 939 | 1 424 397 |
| dec_exp | 839 428 | 438 052 | 1 872 414 | 1 028 597 |
| expkey | 401 554 | 178 556 | 803 435 | 395 581 |
| decode_trap | 763 897 | 392 521 | 1 672 385 | 885 379 |
| poly_ntt257 | 14 760 | 5 692 | 18 820 | 10 813 |
| poly_intt257 | 19 920 | 8 569 | 26 564 | 16 606 |
| poly_basemul257 | 12 153 | 2 669 | 52 170 | 13 256 |
| poly_ntt3329 | 12 715 | 5 053 | 15 542 | 9 429 |
| poly_intt3329 | 15 332 | 6 462 | 19 078 | 12 312 |
| poly_basemul3329 | 12 191 | 3 039 | 55 899 | 14 553 |
| poly_ntt12289 | 34 742 | 7 828 | 71 226 | 17 217 |
| poly_intt12289 | 34 741 | 7 962 | 71 012 | 16 769 |
| poly_rot16 | 5 489 | 5 959 | 11 385 | 13 016 |
| sha3_256（416 B） | 15 556 | 15 628 | 15 651 | 15 560 |

观察：
- NEON 对 NTT-12289（完全 NTT）收益最大（4.4×），对 7 层不完全 NTT 约 2.4–2.6×，对基乘法 4.5×；`decode_trap` 中剩余的标量部分（候选循环、CRT、FCL、旋转）限制了整体收益到 ~1.9×。
- `rot16` 的 NEON 版本未优于标量（桶形移位各级都要整段读写并经 `memcpy`），列入后续优化。
- KeyGen 由 ntrugen 求解器主导（后端无关）。

## 3. Cortex-M4（QEMU mps2-an386，`-O3`，newlib-nano）动态指令数

| 操作 | END-512 ref | END-512 DSP | END-1024 ref | END-1024 DSP |
|---|---:|---:|---:|---:|
| keypair | 28 768 853 | 28 513 991 | 82 447 100 | 82 394 740 |
| enc | 379 828 | 276 825 | 591 952 | 542 904 |
| dec（紧凑） | 4 238 364 | 2 750 745 | 7 759 575 | 6 042 204 |
| dec_exp | 3 018 418 | 1 982 242 | 5 587 829 | 4 352 614 |
| expkey | 1 220 039 | 768 596 | 2 171 839 | 1 689 683 |
| decode_trap | 2 746 870 | 1 762 505 | 5 073 322 | 3 887 155 |
| poly_ntt257 | 41 666 | 29 054 | 81 858 | 55 678 |
| poly_intt257 | 58 912 | 47 692 | 116 000 | 93 132 |
| poly_basemul257 | 38 963 | 11 068 | 53 813 | 53 813 |
| poly_ntt3329 | 47 421 | 29 432 | 93 245 | 56 056 |
| poly_intt3329 | 65 583 | 32 682 | 129 071 | 62 506 |
| poly_basemul3329 | 47 797 | 20 281 | 68 533 | 68 533 |
| poly_ntt12289 | 74 690 | 61 749 | 163 020 | 137 040 |
| poly_intt12289 | 82 041 | 54 434 | 176 527 | 116 160 |
| poly_rot16 | 26 850 | 26 259 | 58 512 | 57 579 |
| sha3_256（416 B） | 64 104 | 64 104 | 64 104 | 64 104 |

观察：
- 打包 16 位 DSP 蝶形把 NTT-3329 逆向降到 1/2，NTT-257 前向降到 70%；NTT-12289 收益较小（每层 Barrett 用标量半字处理，且 len=1 层在字内）。
- END-1024 的基乘法（8 次叶子）DSP 版本与 ref 相同（未做 SMUAD 组合，后续工作）。
- Keccak 占 Encaps 的 ~17%、Decaps 的 ~7%（XKCP 汇编可再降一半）。
- 以 M4 常见的 CPI≈1.2–1.5 估算：Encaps ≈ 0.35–0.4 M cycles，Decaps（展开）≈ 2.4–3 M cycles（168 MHz 下约 15–18 ms），KeyGen ≈ 35–45 M cycles（约 0.2–0.3 s）。**真机数值待测**。

## 4. 主机 x86-64（WSL，gcc 13，`-O3 -march=native`，可移植 C，rdtsc 中位数）

| 操作 | END-512 | END-1024 |
|---|---:|---:|
| keypair | 3 672 398 | 11 556 351 |
| enc | 34 556 | 69 734 |
| dec（紧凑） | 239 766 | 520 562 |
| expkey_from_sk | 80 636 | 164 856 |
| dec_exp | 158 426 | 352 888 |
| poly_ntt257 / intt257 / basemul257 | 1 980 / 3 282 / 2 264 | 3 278 / 5 258 / 8 470 |
| poly_ntt3329 / intt3329 / basemul3329 | 2 160 / 2 910 / 4 134 | 3 106 / 4 076 / 10 784 |
| poly_ntt12289 / intt12289 | 5 992 / 6 344 | 12 022 / 12 106 |
| poly_rot16 | 924 | 1 676 |
| sha3_256（416 B） | 2 016 | 2 056 |
| decode_trap | 146 224 | 296 608 |

## 5. 存储与栈

| 项 | END-512 | END-1024 |
|---|---|---|
| pk / ct / sk / ss | 514 / 384 / 3618 / 32 B | 1027 / 832 / 7235 / 32 B |
| 展开密钥 `end_expkey` | 约 19 KB | 约 38 KB |
| KeyGen scratch（静态或栈） | 11 280 B（实测使用 11 264） | 22 544 B |
| `decode_trap` 栈帧 | ≈ 20 KB | ≈ 40 KB |
| Cortex-M4 镜像 bench.elf（含 newlib-nano/半主机/测试代码） | text 89 KB（ref）/ 95 KB（DSP），bss 65 KB | text 98 KB / 100 KB，bss 112 KB |
| Cortex-M4 栈峰值（涂色法，bench 全部操作含 KeyGen） | 26 052 B（ref）/ 26 080 B（DSP） | 51 764 B / 51 768 B |

bss 中包含静态展开密钥（19/38 KB）、静态 KeyGen scratch（11/22 KB）与基准程序的静态缓冲；库本身的 .data 可忽略。栈峰值主要来自 `decode_trap`（≈20/40 KB）加 KEM/PKE 层的缓冲，后续可通过复用缓冲区压缩（DESIGN §13）。复现命令：`bash scripts/m_stack.sh`。

## 6. pqm4 框架下与 ML-KEM 的对比（2026-11-12）

END 已接入 pqm4（见 `pqm4/README.md`），以下数据都用 pqm4 自己的 harness、HAL、随机源和 Keccak 汇编测得，可以直接与 pqm4 中其他方案比较。平台为 QEMU mps2-an386（Cortex-M4），pqm4 提交 `0dd10ee`，`arm-none-eabi-gcc` 13.2，`-O3`。

**动态指令数**（QEMU TCG 插件 libinsn，`pqm4/icount_pqm4.sh`，单次运行，pqm4 的确定性随机源）：

| 方案 / 实现 | KeyGen | Encaps | Decaps |
|---|---:|---:|---:|
| END-512 m4（DSP） | 27 954 153 | 241 173 | 2 727 488 |
| END-512 ref | 28 523 996 | 266 039 | 3 591 739 |
| ML-KEM-512 m4fspeed | 353 341 | 350 139 | 376 445 |
| ML-KEM-512 m4fstack | 353 830 | 352 041 | 378 729 |
| ML-KEM-512 clean | 484 775 | 566 409 | 696 340 |
| ML-KEM-768 m4fspeed | 569 590 | 583 332 | 619 022 |
| END-1024 m4（DSP） | 128 739 670 | 490 789 | 6 008 730 |
| END-1024 ref | 130 180 098 | 538 570 | 7 745 877 |
| ML-KEM-1024 m4fspeed | 900 011 | 909 698 | 956 512 |
| ML-KEM-1024 clean | 1 238 928 | 1 378 798 | 1 593 997 |

**栈峰值**（pqm4 stack 基准，字节）：

| 方案 / 实现 | KeyGen | Encaps | Decaps |
|---|---:|---:|---:|
| END-512 m4 / ref | 15 644 | 5 024 | 41 328 / 41 324 |
| ML-KEM-512 m4fspeed / m4fstack | 4 372 / 2 300 | 5 436 / 2 348 | 5 412 / 2 332 |
| END-1024 m4 / ref | 34 024 | 11 324 | 82 416 / 82 412 |
| ML-KEM-1024 m4fspeed / m4fstack | 6 436 / 3 332 | 7 500 / 3 372 | 7 484 / 3 356 |

**代码尺寸**（pqm4 size 基准，`.text` 字节）：END-512 m4 72 218、ref 70 550；END-1024 m4 77 122、ref 75 334；ML-KEM-512 m4fspeed 15 848；ML-KEM-1024 m4fspeed 16 916。

**哈希占比**（pqm4 hashing 基准，QEMU 时钟下的比例，仅作参考）：END-512 m4 KeyGen 6.8%、Encaps 9.2%、Decaps 0.7%；ML-KEM-512 m4fspeed 为 52.6%、13.2%、17.4%。

结论：
- **Encaps**：END-512 比 ML-KEM-512 m4fspeed 少约 31% 指令（241 k 对 350 k），END-1024 比 ML-KEM-1024 少约 46%。END 的封装只有一次 NTRU 乘法加压缩。
- **Decaps**：END-512 约为 ML-KEM-512 m4fspeed 的 7.2 倍，END-1024 约为 6.3 倍。主要开销在陷门解码（9 个候选，每个候选一次 NTT-3329 乘法）和 NTT-12289，是后续手写汇编优化的重点。
- **KeyGen**：END 需要求解 NTRU 方程（ntrugen），约为 ML-KEM 的 80–140 倍，与论文的定位一致（适合公钥长期复用）。
- **栈与代码**：END 的解封装栈约为 ML-KEM m4fspeed 的 7.6（512）/ 11（1024）倍，代码尺寸约为 4.5 倍，主要来自求解器和展开密钥。
- 以上为指令数。真板周期数需在 NUCLEO-L4R5ZI 上用 `benchmarks.py -p nucleo-l4r5zi` 测量（命令见 `pqm4/README.md`）。

## 7. 测量方法与复现

```bash
# 主机周期数
make END_N=512 BACKEND=ref CC=gcc CFLAGS="-O3 -march=native" ARCH=hostO3 build/hostO3-ref-512/bench && ./build/hostO3-ref-512/bench 300
# 动态指令数（需 --enable-plugins 的 QEMU，见下）
bash scripts/icount_all.sh build/icount_results.txt
# Cortex-M4 栈峰值与尺寸
bash scripts/m_stack.sh
```

QEMU 插件：Ubuntu 的 QEMU 包不含插件支持，本仓库在 WSL 中从源码构建了 QEMU 9.1.3（`--enable-plugins --target-list=aarch64-linux-user,arm-softmmu,aarch64-softmmu`，安装到 `/opt/qemu-plugins`，插件 `tests/tcg/plugins/libinsn.so`）。`bench --insn <op>` 单次执行一个操作，脚本用 `总指令数(op) − 总指令数(none)` 得到每操作的指令数。

真机：
- Cortex-A Linux：`bench/cycles.c` 默认读 `CNTVCT_EL0`（通用计时器，换算为 ns）；编译 `-DEND_AARCH64_PMU` 时读 `PMCCNTR_EL0`（需 `echo 1 > /proc/sys/kernel/perf_user_access` 或内核模块开放用户态访问）；固定频率、绑核、取中位数；同平台以 ML-KEM-512 NEON 实现为对照。
- Cortex-M：`bench/cycles.c` 使用 DWT CYCCNT；`platform/cortex-m/startup.c` 栈涂色；按 pqm4 惯例降频运行（如 24 MHz）以消除 flash 等待周期；同板对照 pqm4 的 ML-KEM-512。
