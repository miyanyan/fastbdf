# EDFlib 1.27 功能对照

基线为 EDFlib v1.27 和 fastbdf 提交 `72e7ff3`。下表收录上游头文件中全部 **50 个不同名称的公开函数**，包括废弃别名；同名的编译器条件声明计一次。比较目标是 C++ 调用方可用的功能、生成的文件和读回结果。C ABI、全局句柄编号和整数错误码不作为替代接口承诺。[上游头文件](https://gitlab.com/Teuniz/EDFlib/-/blob/v1.27/edflib.h) · [上游实现](https://gitlab.com/Teuniz/EDFlib/-/blob/v1.27/edflib.c)

**互验**表示现有测试已用真实 EDFlib 对照该行注明的结果；**原生**表示有 fastbdf 自身断言；**待验**表示有实现或组合用法，但缺少所列行为的直接对照；**差异**表示源码或测试已经能确认行为不同。同一行可以同时有已验证的正常路径和待验证的边界。

## 打开、关闭和定位

| EDFlib 入口 | fastbdf 用法 | 已有证据或差异 |
|---|---|---|
| `edfopen_file_readonly` | `Reader::open(path, mode)`，通过 `header()`、`signals()`、`record_count()` 和 `sample_count(id)` 取头部与计数 | **互验**。连续 EDF+/BDF+ 两个写入方向、经典 EDF/BDF 的 fastbdf 写 → EDFlib 读、注释读取模式。见 `verify`、`verify_writer_formats`、`annotation_modes`。fastbdf 另提供延迟加载。 |
| `edfopen_file_writeonly` | 准备 `RecordingInfo` 及其通道，然后 `Writer::create(path, info)` | **互验**。上游写 EDF+/BDF+，fastbdf 还可写经典格式。 |
| `edfopen_file_writeonly_with_params` | 对每个通道调用 `info.add_signal(label, samples_per_record)`，设置相同单位和对称物理范围后创建 Writer | **待验**。上游会生成从 `chan. 1` 起的标签；C++ 调用方显式给出标签。尚无两种创建方式的头部对照。 |
| `edfclose_file` | Reader 随对象销毁关闭；Writer 显式 `finish()`，或由 `FileSet::close()` 完成写入 | **互验**正常文件收尾；**差异**：单独销毁 Writer 不会代为调用 `finish()`。 |
| `edfread_digital_samples`、`edfread_physical_samples` | `Reader::read_digital`、`read_physical`；绝对位置版本为 `read_digital_at`、`read_physical_at` | **互验**。逐通道比较原始数字值和标定后的物理值，覆盖 16/24 位；跨记录短读另有预期样本断言。见 `verify`、`verify_external`。 |
| `edfseek`、`edftell`、`edfrewind` | `reader.cursor(id)` 的 `seek/tell/rewind`；Reader 自带的每通道位置也可用 `seek/tell` 操作 | **原生**：相对定位、回绕、越界、独立游标。**差异**：上游把超出范围的位置夹到文件首尾，C++ 游标抛出异常。**待验**：对同一混合采样率文件交错执行定位和读取，逐步比较返回位置与样本。 |
| `edflib_is_file_used`、`edflib_get_number_of_open_files`、`edflib_get_handle` | `FileSet::contains/size/handles` | **原生、差异**。每个 FileSet 只跟踪经它打开的对象。上游查询进程内 EDFlib 的全局表。 |
| `edflib_version` | 包版本见 CMake 项目版本；测试调用上游函数确认返回 127 | **差异**。fastbdf 没有返回 EDFlib 版本号的公共 C++ 函数。 |

## 头部、标定与记录布局

下列设置在创建 Writer 前填入。读回时查看 `Reader::header()` 和 `Reader::signals()`。上游部分文本 setter 会处理遗留编码；fastbdf 的头部字段使用可打印 ASCII，可显式调用 `header_text_from_latin1`。注释正文的规则见下一节。

| EDFlib 入口 | fastbdf 字段或选项 | 已有证据或差异 |
|---|---|---|
| `edf_set_samplefrequency` | `SignalSpec::samples_per_record` | **互验**。值是每条数据记录的采样点数。测试使用 0.5 秒记录、每通道分别 7 和 3 点；实际 Hz 由 `sample_rate(id)` 返回。 |
| `edf_set_physical_minimum`、`edf_set_physical_maximum`、`edf_set_digital_minimum`、`edf_set_digital_maximum` | `SignalSpec::physical_min/physical_max/digital_min/digital_max` | **互验**。既对照头部四个端点，也对照样本标定结果。负增益和量化临界值的组合仍可扩充。 |
| `edf_set_label`、`edf_set_physical_dimension` | `SignalSpec::label/unit` | **互验**已有 ASCII 标签和单位样例。**待验**：字段宽度边界及不同输入编码。 |
| `edf_set_prefilter`、`edf_set_transducer` | `SignalSpec::prefilter/transducer` | **原生**有读回传感器文字的断言；**待验**两字段的非空值经两个写入方向逐字段比较。 |
| `edf_set_startdatetime`、`edf_set_subsecond_starttime` | `RecordingInfo::start/subsecond` | **互验**亚秒起点和相对事件时间；**待验**两边显式设置同一日历时间并比较全部日期时间分量。 |
| `edf_set_patientname`、`edf_set_patientcode`、`edf_set_equipment` | `RecordingInfo::patient_name/patient_code/equipment` | **互验**。两个读取器都核对非空患者名、患者代码和设备名。见 `verify`。 |
| `edf_set_sex`、`edf_set_gender`、`edf_set_birthdate` | `RecordingInfo::sex/birthdate`，出生日期采用 `YYYY-MM-DD` | **原生**有出生日期样例；`gender` 是上游废弃别名。**待验**不同性别值、缺失值和出生日期的双向字段断言。 |
| `edf_set_admincode`、`edf_set_technician` | `RecordingInfo::admin_code/technician` | **待验**。支持序列化和解析，尚缺两个库对非空字段的逐项断言。 |
| `edf_set_patient_additional`、`edf_set_recording_additional` | `RecordingInfo::patient_additional/recording_additional` | **互验**。覆盖内部空格、下划线、字母 X 和读出后重写。见 `annotation_edges`。 |
| `edf_set_datarecord_duration`、`edf_set_micro_datarecord_duration` | `RecordingInfo::record_duration` 使用 100 ns 的 `Tick` | **互验**常规记录时长。两个上游 setter 的参数单位分别为 10 µs、1 µs。**待验**微秒级极值、允许范围与 8 字节头部字段的可表示性。 |
| `edf_set_number_of_annotation_signals`、`edf_set_annot_chan_idx_pos` | `WriterOptions::annotation_signals/annotation_position` | **互验**单、双注释通道的读取模式；**原生**有首、中、末位置布局断言。**待验**上游设置位置后，两端读取相同的逻辑通道顺序。 |

## 写样本与事件

| EDFlib 入口 | fastbdf 用法 | 已有证据或差异 |
|---|---|---|
| `edfwrite_digital_samples`、`edf_blockwrite_digital_samples` | `Writer::write_digital_samples/write_digital_record` | **互验**。上游逐通道写与 fastbdf 整记录写生成的样本均经双库核对，覆盖 BDF 有符号端点。**待验**同输入直接对照上游 block 入口。 |
| `edfwrite_physical_samples`、`edf_blockwrite_physical_samples` | `Writer::write_physical_samples/write_physical_record` | **互验**。上游逐通道写与 fastbdf 整记录写后比较数字码，已有饱和、截断样例。**待验**两个 block 入口直接对照，以及负增益下相邻量化码的边界。 |
| `edfwrite_digital_short_samples`、`edf_blockwrite_digital_short_samples` | `Writer::write_short_samples/write_short_record` | **差异**。上游两个 short 入口只接受 EDF 系列；fastbdf 也可将 int16 输入写为 24 位 BDF。**原生**覆盖 BDF 整记录；**待验**共同支持的 EDF+ 逐通道和整记录路径。 |
| `edf_blockwrite_digital_3byte_samples` | `Writer::write_packed_record`，传入一条记录普通通道的小端字节 | **互验**。EDFlib 能从 fastbdf 的 packed 文件读回有符号端点；**待验**将同一字节缓冲区交给上游 3-byte 写入口，对照产物。 |
| `edf_get_annotation` | `Reader::annotations()` 返回事件序列，元素包含 `onset/duration/text` | **互验**。包括空正文、零时长与缺失时长、亚秒起点、原始字节、超 512 字节正文后的下一条事件。见 `annotation_edges`、`verify_writer_formats`。 |
| `edfwrite_annotation_utf8_hr`、`edfwrite_annotation_utf8` | `Writer::add_annotation`，以 `Tick` 表示时间，正文为 `std::string` | **互验**高精度入口的共同时值及正文处理：取前 40 字节、NUL 终止、低位控制字符替换为句点。上游两个入口分别采用 µs、100 µs；**待验**旧入口的单位换算和负起点的接受行为。 |
| `edfwrite_annotation_latin1_hr`、`edfwrite_annotation_latin1` | 调用方先将 Latin-1 转为目标字节，再调用 `add_annotation` | **有意差异**。库内没有自动转换注释的入口。若要复现上游对超长 Latin-1 输入的结果，还要考虑上游转换缓冲区与截断顺序。 |

## 下一批互操作用例

1. **优先验证样本值。** 在相同标定下分别经逐通道、整记录和 packed 接口写入。覆盖 EDF+ short、BDF+ 24 位端点及其相邻值、负增益，以及量化分界线两侧的浮点输入。逐点比较读回的数字码。
2. **验证位置和时间。** 对混合采样率通道交错执行合法的起点、当前位置和末尾定位，覆盖跨记录、短读和 EOF；对越界分别断言两个库已知的行为。再用非整秒记录和亚秒起点核对普通与高精度注释入口。
3. **补全头部字段。** 固定起始日期时间，填入非空的传感器、预滤波、管理员代码、技师和出生日期。两个写入方向都逐字段比较，并校验不同注释通道位置下的信号顺序。

证据入口：双库对照见 [interop.cpp](../tests/interop.cpp)，原生断言见 [reader.cpp](../tests/reader.cpp)、[writer.cpp](../tests/writer.cpp) 和 [native.cpp](../tests/native.cpp)。实现入口为 [Reader](../include/fastbdf/reader.hpp)、[Writer](../include/fastbdf/writer.hpp)、[公共数据类型](../include/fastbdf/types.hpp) 和 [FileSet](../include/fastbdf/file_set.hpp)。
