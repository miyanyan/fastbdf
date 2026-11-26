# Python 读取基准

对照入口为已安装的 fastbdf nanobind 扩展和 MNE 的 `read_raw_bdf(preload=False)` / `Raw.get_data()`。两端返回同一通道顺序、形状和单位的 C 连续 float64 矩阵。fastbdf 在计时段内分配数组、批量读取并把每通道的文件单位换算为伏特。比较计入各自的 Python 接口、文件打开和头部解析成本。这个实验测量两个独立 Python 读者；将扩展接入 MNE Raw 后需要重新测量。

## 输入与复现

```sh
pixi run --locked -e mne test-binding
pixi run --locked -e mne bench-binding
```

第二条命令生成 64 通道、512 Hz、120 秒的合成 BDF，数字码有独立公式可校验。需要复现真实文件时，下载 [26 By Biosemi 的 sub-001 BDF](https://s3.amazonaws.com/openneuro.org/ds005284/sub-001/eeg/sub-001_task-26ByBiosemi_eeg.bdf)，然后执行：

```sh
pixi run --locked -e mne python tests/python_bench.py --file sub-001_task-26ByBiosemi_eeg.bdf
```

实测文件为 65,511,936 字节，SHA-256 为 `a3dcbdc69c5cd2b47c19cd2fcbe27ddbca0baa3deca6881b37e8e02fed970cda`。数据来自 Zhao Xiangyue 等人的 OpenNeuro ds005284 v1.0.0，CC0。记录时长 328 秒，采样率 1024 Hz。对照选取 A1–A32 和 B1–B32 共 64 个 EEG 通道，Status 为独立触发通道，未纳入电压比较。[数据集及作者](https://doi.org/10.18112/openneuro.ds005284.v1.0.0)

## 本机结果

Windows 11 x64、AMD Ryzen 7 9700X、MSVC Release、Python 3.12.16、NumPy 2.5.5、nanobind 3.1.0。首屏和完整读取含打开、读取、关闭，预热后交替顺序测量 5 对，报告中位数。窗口测试保持文件打开，使用同一随机种子的 200 个请求，每次随机选 4 个通道及一个起始样本读取 4 秒，报告 P95。各次输出在计时外逐点核对，物理值容差为 `rtol=2e-11, atol=1e-15 V`。文件处于操作系统热缓存。峰值内存按每个读者单独启动新进程测量，含公共模块导入、各自运行时和完整输出矩阵。

| 数据与对照版本 | 读者 | 首屏 10 秒 | 随机窗 P95 | 全量读取 | 全量进程峰值内存 |
|---|---|---:|---:|---:|---:|
| 合成 / MNE 1.13.2 | fastbdf | 0.901 ms | 0.151 ms | 10.686 ms | 76.0 MiB |
| | MNE | 4.478 ms | 0.596 ms | 18.608 ms | 121.6 MiB |
| 真实 / MNE 1.13.2 | fastbdf | 1.612 ms | 0.169 ms | 54.422 ms | 210.5 MiB |
| | MNE | 5.576 ms | 1.037 ms | 94.956 ms | 255.7 MiB |
| 真实 / MNE main | fastbdf | 1.562 ms | 0.170 ms | 49.926 ms | 209.6 MiB |
| | MNE | 5.519 ms | 1.062 ms | 87.303 ms | 254.6 MiB |

MNE main 对照单独运行，将固定提交的源码目录加入同一解释器的 `PYTHONPATH`，已核对实际导入路径。稳定版和 main 表中的 fastbdf 数字来自各自配对的那一轮。单机热缓存结果不涵盖冷盘、不同存储设备或任意通道组合。[MNE 源码提交](https://github.com/mne-tools/mne-python/commit/c128077643f905ebf8b75a6f7664c8c20243e39c)
