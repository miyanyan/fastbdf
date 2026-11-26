# fastbdf

C++20 header-only 库，用于读取 EDF、BDF 波形和 EDF+/BDF+ 注释。可以按样本位置、通道游标或多通道时间窗读取，也可以逐通道或整条记录写入。应用只需 C++20 标准库和平台文件 API。EDFlib 1.27 用作独立的测试参考库。

## 接入应用

将 `include` 加入头文件搜索路径并启用 C++20，或通过 CMake 的 `add_subdirectory` 使用 `fastbdf::fastbdf`。也可以安装后使用包配置。

```sh
cmake -S . -B build/package -DCMAKE_INSTALL_PREFIX=build/install
cmake --install build/package
```

给调用项目的 `CMAKE_PREFIX_PATH` 指定安装目录的绝对路径。

```cmake
find_package(fastbdf CONFIG REQUIRED)
target_link_libraries(my_application PRIVATE fastbdf::fastbdf)
```

## 读取混合采样率的波形

```cpp
#include <fastbdf/fastbdf.hpp>
using namespace std::chrono_literals;

auto file = fastbdf::Reader::open("night.bdf");
auto ids = file.select_unique({"Fp1", "Temp"});
auto values = file.read_physical(ids, {30s, 40s});
auto first = file.sample_index_at_or_after(ids[0], 30s);
```

`values[i]` 对应 `ids[i]`，各通道保留自己的样本数。时间相对于实际记录起点，窗口包含起点并排除终点。超出文件尾的部分返回空或短结果。通道 ID 从 0 开始，仅计入普通信号。标签有重名时可以遍历 `signals()` 选择 ID。反复翻页可用 `read_physical_into(ids, range, values)` 或 `read_digital_into` 复用结果容量。

对于已有的样本，`file.sample_time(id, index)` 返回相对于第一条记录的精确采样时间，包含记录间断。结果由 `std::chrono::nanoseconds` 类型的 `whole_nanoseconds` 和额外的 `numerator / denominator` 纳秒组成。例如，1 秒记录内有 3 个样本时，下标 1 对应 333,333,333 个整纳秒加 1/3 纳秒。

```cpp
std::array<std::int32_t, 4096> chunk;
auto cursor = file.cursor(ids[0]);
cursor.seek(-512, std::ios_base::end);  // 该通道至少有 512 个样本
auto n = cursor.read(chunk);    // 只使用 chunk 的前 n 个元素
```

`file.read_at(id, first, out)` 支持按绝对位置并发读取。`out` 可以是可写数组、`std::vector` 或 `std::span`；`std::int32_t` 接收数字码，`double` 接收物理值。返回值是实际读取的样本数。每个 `Cursor` 独立维护位置，支持 `seek/tell/rewind`。游标借用 Reader，须在 Reader 移动、重新打开或销毁前结束使用。便捷的 `file.read(id, out)` 使用文件对象自己的每通道游标，多线程共享这些隐式游标时由调用方同步。

## 写入并补记事件

```cpp
#include <array>
#include <fastbdf/fastbdf.hpp>
using namespace std::chrono_literals;

fastbdf::RecordingInfo info;
info.format = fastbdf::Format::bdf_plus;
info.record_duration = 1s;
auto& eeg = info.add_signal("Fp1", 4);  // 每条记录 4 个样本
eeg.unit = "uV";
eeg.physical_min = -100;
eeg.physical_max = 100;
info.add_signal("Temp", 1).unit = "C";

auto writer = fastbdf::Writer::create("new.bdf", info);
writer.write_digital_record(std::array<std::int32_t, 5>{-200, 0, 200, 400, 20});
writer.add_annotation({250ms, std::nullopt, "stimulus"});
writer.finish();
```

`create` 会覆盖目标路径。起始时间默认取本地系统时间，导出历史数据时应设置 `info.start` 和 `info.subsecond`。先选择格式再调用 `add_signal`，它会按 16 或 24 位设置数字范围。每次整记录输入按文件通道顺序连续排列。设备已输出小端样本时，`write_packed_record` 可直接接收普通通道字节，注释通道由库插入。

逐通道写入使用 `write_samples(samples)`，也可以用 `write_samples(channel, samples)` 显式核对通道顺序。每次传入该通道一条记录的样本，支持 C 数组、`std::array`、`std::vector` 和 `std::span` 等连续范围。`std::int16_t` 和 `std::int32_t` 按数字码写入，`float` 和 `double` 按物理值标定，换算使用双精度。

`flush()` 刷出已完成的记录。`finish()` 校验是否有未完成的通道、写回晚到事件和记录数，成功后可以重复调用。要确认文件写入成功，必须显式调用它。写入期间内存保存尚待写回的事件。通过 `WriterOptions` 预留每条记录的注释通道数、每通道字节数和位置，容量不足会报错。

## 格式与兼容行为

| 能力 | 行为 |
|---|---|
| 经典 EDF、BDF | 读写 16 位或 24 位小端信号，保留经典 patient 和 recording 字段 |
| 连续 EDF+、BDF+ | 读写结构化元数据、负增益、亚秒起点、混合采样率、UTF-8 事件、多注释通道及仅含注释的文件 |
| 不连续 EDF+D、BDF+D | 可读取有正记录时长的波形，`record_start(r)` 返回考虑间断后的相对时间。时间窗跳过间断，不填充样本。Writer 输出连续记录。零记录时长仅支持纯注释文件 |
| 数值 | digital 输出按头部数字范围截断；physical 依据头部标定转换。physical 写入将物理端点精确映射到对应数字端点，越界值饱和，区间内换算后向零取整。packed 写保留输入字节 |
| 时间与编码 | 时间使用标准 `std::chrono::nanoseconds`，精度为 1 ns。事件缺失时长用 `std::nullopt`，零时长用 `std::chrono::nanoseconds{0}`。事件正文统一使用 `std::string`。调用方负责提供符合 EDF+ 规范的 UTF-8 字节。头部字段由调用方提供可打印 ASCII 字节 |
| 注释读取 | 默认 `lazy` 在首次访问 `annotations()` 时加载。`none` 跳过事件，`all` 在打开时加载全部，`until_recording_ends` 复现 EDFlib 的停止标记规则。延迟加载时，连续记录的首次打开只需读取头部和第一条记录的计时信息；不连续记录要建立时间索引 |
| 错误与文件生命周期 | 捕获 `fastbdf::Error` 并查看 `code()`。应用直接持有 `Reader` 和 `Writer`；读完后释放 Reader，写完后调用 `Writer::finish()` |

连续记录的总时长可以超出 `std::chrono::nanoseconds` 的表示范围。按样本下标仍可读取末尾波形；若请求的记录起点超出范围，`record_start(r)` 会抛出 `fastbdf::Error`。时间窗的两个端点及解析的 EDF+/BDF+ 时间戳仍须能用 `std::chrono::nanoseconds` 表示。

EDFlib 的 `*_hr` 注释写函数接收微秒，旧函数接收 100 微秒单位，迁移时使用相应的 `std::chrono::duration` 转换。EDFlib 用户从原生 C++ 接口迁移，本项目不导出其 C ABI。

### 与 EDFlib 1.27 交换文件

上游 50 个公开函数对应的 C++ 用法、已验证行为和待补测试见 [EDFlib 1.27 功能对照表](docs/edflib-1.27.md)。

双向测试覆盖连续 EDF+ 和 BDF+。测试用 EDFlib 生成文件，由 fastbdf 读取实际头部、样本和事件并重新生成文件，最后由 EDFlib 打开并核对。经典 EDF、BDF 可以由 fastbdf 生成后交给 EDFlib 读取；EDFlib 的写接口只生成 Plus 格式。

- 空描述事件保留起点、时长和顺序。每条记录开头专门用于计时的空描述不会计入事件数。
- 事件时间的小数部分保留至纳秒，超过 9 位的小数向零截断。fastbdf 解析完整的 TAL 时间字段，不采用 EDFlib 持续时间的前 15 字符限制；与 EDFlib 交换文件时应注意其时间精度及解析范围。
- 事件正文按字节处理。写入最多取 40 字节，遇到 NUL 提前结束，其余小于 0x20 的字节及 DEL（0x7F）替换为句点。读取最多返回 512 字节，允许 TAB、LF、CR，拒绝正文中其余 ASCII 控制字符及 DEL。两端均不校验或转换 UTF-8，按字节截断可能切开多字节字符。结构化患者和记录字段中的附加说明保留原有下划线及内部空格。
- 有波形通道的文件必须写完至少一条记录才能成功 `finish()`。可在补齐后重试。

## Python / NumPy 读取

从源码目录运行 `python -m pip install .`，需要 C++20 编译器。Python 扩展使用 nanobind；C++ 用户仍可直接包含头文件。当前 Python 接口提供 Reader，可读取 EDF、BDF 和 Plus 文件。

```python
import fastbdf
import numpy as np

with fastbdf.Reader("recording.bdf") as reader:
    signal = reader.signals[0]
    codes = np.empty(1024, dtype=np.int32)
    n = reader.read_digital_into(0, 0, codes)
    first_second = reader.read_physical_window([0, 1], 0, 1_000_000_000)
```

`read_digital_into` 和 `read_physical_into` 从通道的原始样本下标开始，返回实际写入点数。使用可写、连续、对齐的 NumPy 数组，数字码用 `int32`，物理值用 `float64`，单位见 `signal.unit`。时间窗使用整数纳秒，包含起点并排除终点。各通道样本数相同时，`read_physical_window_into(ids, start_ns, stop_ns, out)` 可直接填充形状为「通道数 × 样本数」的二维数组。混合采样率可用 `read_physical_window`，按请求顺序返回独立的一维数组。I/O 期间释放 GIL，并发请求应使用独立的输出缓冲区。

Python 默认打开时解析事件。`reader.annotations` 的元素是 `(onset_ns, duration_ns_or_None, text_bytes)`，也可在打开时指定 `annotation_mode=fastbdf.AnnotationMode.lazy`。构建、测试及基准命令见 [Python 绑定与 MNE 性能对照](docs/python-benchmark.md)。

## 与 MNE-Python 对照

正确性测试读取 fastbdf 生成的 EDF、BDF、EDF+ 和 BDF+，逐通道核对采样率、样本数、物理量和读取窗口，并比较 Plus 文件的事件文本、起点和持续时间。反向测试用 MNE 和 edfio 导出 EDF，再由 fastbdf 读取，并按 EDF 的量化精度与原始浮点信号比较。

MNE 默认将电压信号换算为伏特，fastbdf 保留通道头部的单位。对照时统一换算为伏特。混合采样率文件按通道单独交给 MNE 打开，以保留各自的原始采样点。[MNE EDF 读取接口](https://mne.tools/stable/generated/mne.io.read_raw_edf.html) · [MNE BDF 读取接口](https://mne.tools/stable/generated/mne.io.read_raw_bdf.html)

```sh
pixi run --locked -e mne test-mne
pixi run --locked -e mne bench-mne
```

Windows 和 Ubuntu CI 运行正确性测试。性能测试在本机单独运行，比较同采样率的连续 EDF+ 与 BDF+。两端都计入打开文件、读取注释、分配并读取全部通道的 float64 样本、换算为伏特和关闭文件。每种格式预热一次，取 5 次的中位数；进程启动、Python 导入、生成文件和逐点核对均在计时外。

## 开发与验证

开发时安装一次提交钩子。提交时会按仓库的 `.clang-format` 格式化 C++，用 Ruff 检查并格式化 Python，同时检查配置文件、空白和合并冲突标记。

```sh
pipx install pre-commit
pre-commit install
pre-commit run --all-files
```

钩子自动修正文件后，检查差异并重新暂存再提交。CI 对整个仓库运行同一套检查。

安装 Pixi、C++20 编译器和已 bootstrap 的 vcpkg，设置 `VCPKG_ROOT`。Windows 使用带 C++ 桌面开发组件的 Visual Studio，Linux 可用 GCC。Pixi 锁定开发工具，vcpkg manifest 的 `interop` feature 提供 EDFlib 1.27。

```sh
pixi install --locked
pixi run test
pixi run bench
```

原生测试覆盖手工构造的经典格式、混合采样率、间断、负增益、截断文件、UTF-8、时长、独立游标与并发定位读取、跨翻译单元链接，以及写入中断后的行为。互操作测试让两个库互读对方生成的文件。也可给 `fastbdf_interop` 可执行文件传入一个或多个 EDF/BDF 路径，分块比对全部 digital、physical 和注释。Windows 开发预设下可执行文件位于 `build/dev/Release/`。

已用 Sleep-EDF Expanded 的一份 22 小时 5 分钟记录验证全部 24,168,000 个数字样本；物理量误差在 `1e-10 × max(1, |参考值|)` 内。配套文件的 154 条事件按文本、起点、持续时间逐条匹配。测试数据独立下载，常规测试不需要联网下载录制数据。[PSG 文件](https://physionet.org/files/sleep-edfx/1.0.0/sleep-cassette/SC4001E0-PSG.edf) · [标注文件](https://physionet.org/files/sleep-edfx/1.0.0/sleep-cassette/SC4001EC-Hypnogram.edf)

基准使用 16 通道、每记录每通道 1024 点的合成 BDF+，比较相同输入的整文件、单通道、多通道时间窗和短随机读取，核对校验和。将具体结果连同编译器、文件大小和缓存状态一起解释。


在本机 Windows x64、MSVC 19.51 Release、1,000 条记录、46.875 MiB 样本的热缓存测试中，得到以下结果。吞吐操作取预热后 5 次的中位数；短读 P95 取 512 对请求，首屏 P95 取 40 对打开并读取 10 秒全通道数据的请求。写入计时包含创建、写记录和关闭，结果在计时外校验。

| 工作负载 | fastbdf | EDFlib 1.27 |
|---|---:|---:|
| 全通道批量读取 | 27.91 ms | 668.42 ms |
| 单通道完整读取 | 5.04 ms | 41.86 ms |
| 整文件数字样本写入 | 29.23 ms | 73.20 ms |
| 打开到首屏 P95 | 0.89 ms | 17.43 ms |
| 随机读取 64 点 P95 | 8.4 µs | 8.8 µs |

以上数字反映这一台机器和这一种缓存状态。真实睡眠记录用于验证兼容性。Windows 和 Linux 已执行双向互操作、安装后消费测试，Linux 另通过 AddressSanitizer 和 UndefinedBehaviorSanitizer 原生测试。
