# END KEM 接入 pqm4

本目录把 END-512 / END-1024 接入 [pqm4](https://github.com/mupq/pqm4)（Cortex-M4 后量子密码基准框架），用 pqm4 统一的测试与基准流程和 ML-KEM 等算法做同口径比较。

| 文件 | 作用 |
|---|---|
| `export_pqm4.py` | 把本仓库源码导出为 pqm4 要求的扁平实现目录 |
| `pqm4_glue.c` | 把 `end512_crypto_kem_*` 映射为 pqm4 调用的 `crypto_kem_*`，`randombytes` 接到 pqm4 的随机源 |
| `icount/` | KEM 的累积 harness（`icount_{none,kg10,keypair,enc,dec}.c`），用于指令数和周期模型，可选安装 |
| `icount_sign/` | 签名方案的同类 harness（KeyGen / Sign / Verify），用于用 ML-DSA 的真板数据校验周期模型 |
| `icount_pqm4.sh` | 在 QEMU + TCG 插件下统计任意 pqm4 KEM 的动态指令数 |
| `m4model/m4cycles.c` | QEMU TCG 插件：按 Cortex-M4 TRM 时序估计周期数，并按函数剖析 |
| `m4model/run_m4model.sh`、`m4model/m4model_report.py` | 对任意 pqm4 实现跑周期模型，汇总成表，并与 pqm4 公布的真板数据对比 |

## 导出的实现

| pqm4 路径 | 内核 | 说明 |
|---|---|---|
| `crypto_kem/end512/m4f`、`crypto_kem/end1024/m4f` | `src/armv7em` + `*.S`（`END_M4_ASM`） | Cortex-M4 手写汇编内核，最快的版本；用 FPU 寄存器作暂存，所以按 pqm4 惯例叫 `m4f` |
| `crypto_kem/end512/m4`、`crypto_kem/end1024/m4` | `src/armv7em`（C + DSP） | 不依赖 FPU 的优化版 |
| `crypto_kem/end512/ref`、`crypto_kem/end1024/ref` | `src/ref`（可移植 C） | 参照实现；pqm4 还会在主机上编译它，`testvectors.py` 用它与板上的 m4 版逐字节比对 |

导出时每个实现会生成 `api.h`（pqm4 需要的尺寸常量与函数原型）和 `end_config.h`（`END_N`、`END_PQM4`，`m4f` 还有 `END_M4_ASM`）。`END_PQM4` 让 KeyGen 的 SHAKE256 走 pqm4 自带的 Keccak 汇编（`src/common/end_xof.h`），因此 pqm4 的 hashing 基准能统计到 END 的哈希时间，结果与本仓库自带的 `fips202.c` 逐字节相同。

导出的目录是生成物：修改源码后重新运行导出脚本，不要直接改 pqm4 里的副本。

## 使用步骤

pqm4 仓库含大量符号链接，必须在 Linux 文件系统中克隆（WSL 下用 `/root/...` 这类路径，不要放在 `/mnt/c`、`/mnt/d`）。

```bash
git clone --recurse-submodules https://github.com/mupq/pqm4.git ~/pqm4
python3 pqm4/export_pqm4.py ~/pqm4            # 加 --icount 同时安装指令数 harness；--clean 卸载
cd ~/pqm4
# QEMU（功能、栈、代码尺寸）
python3 test.py        -p mps2-an386 end512 end1024
python3 testvectors.py -p mps2-an386 end512 end1024
python3 benchmarks.py  -p mps2-an386 end512 end1024 ml-kem-512 ml-kem-768 ml-kem-1024
python3 convert_benchmarks.py md
# 真机周期数（NUCLEO-L4R5ZI，需 ST-Link/OpenOCD，参考 pqm4 README）
python3 benchmarks.py  -p nucleo-l4r5zi -u /dev/ttyACM0 end512 end1024 ml-kem-512 ml-kem-1024
```

指令数对比（需要用 `--enable-plugins` 从源码编译的 QEMU，指令计数插件是其中的 `tests/tcg/plugins/libinsn.so`）：

```bash
python3 pqm4/export_pqm4.py ~/pqm4 --icount
bash pqm4/icount_pqm4.sh ~/pqm4 crypto_kem/end512/m4 crypto_kem/ml-kem-512/m4fspeed
```

Cortex-M4 周期模型（同样需要带插件的 QEMU；没有真板时用它代替 QEMU 的周期数）：

```bash
python3 pqm4/export_pqm4.py ~/pqm4 --icount
bash pqm4/m4model/run_m4model.sh ~/pqm4 ~/m4model crypto_kem/end512/m4 crypto_kem/ml-kem-512/m4fspeed \
     mupq/pqclean/crypto_kem/ml-kem-512/clean crypto_sign/ml-dsa-44/m4f
python3 pqm4/m4model/m4model_report.py ~/m4model --board ~/pqm4/results/speed/ml-kem.md \
     --board ~/pqm4/results/speed/ml-dsa.md --profile end512/m4
```

每个实现跑 5 个程序（空程序、10 次 KeyGen、累积到 KeyGen / Encaps / Decaps），由差值得到各操作的周期数和指令数。`--board` 读取 pqm4 的真板结果表，并输出模型误差。

## 注意事项

- **QEMU 下的周期数不能用于比较。** mps2-an386 的 SysTick 只是模拟器时钟，例如 ML-KEM 的 clean 版在 QEMU 中反而比汇编版"快"。在 QEMU 上请看栈、代码尺寸、指令数和周期模型的估计；准确的周期数要在 NUCLEO-L4R5ZI 真板上测。
- pqm4 的平台编译选项带 `-ffreestanding`，普通 `memcpy` 不会被内联。DSP 内核的系数对读写因此使用 `__builtin_memcpy`。
- END 在 pqm4 中按 pqm4 惯例全部使用栈内存（不使用 `END_STATIC_*` 静态缓冲），栈数据可以直接与其他方案比较。
- 以下告警来自 ntrugen 头文件中未使用的 static 原型，无害：`'zint_co_reduce' declared 'static' but never defined`。
