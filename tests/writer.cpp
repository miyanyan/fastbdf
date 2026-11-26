#include <fastbdf/writer.hpp>
#include <fstream>
#include <iterator>
#include <limits>
#include <stdexcept>

using namespace fastbdf;
namespace
{
void check(bool condition)
{
    if (!condition) throw std::runtime_error("writer edge check failed");
}
template<class F>
void throws(F f)
{
    bool caught = false;
    try {
        f();
    }
    catch (const Error&) {
        caught = true;
    }
    check(caught);
}
std::string contents(const std::filesystem::path& p)
{
    std::ifstream f(p, std::ios::binary);
    return {std::istreambuf_iterator<char>(f), {}};
}
int sample(const std::string& s, std::size_t p)
{
    auto n = static_cast<unsigned char>(s[p]) | (static_cast<unsigned char>(s[p + 1]) << 8);
    return n >= 32768 ? n - 65536 : n;
}
} // namespace
int main()
{
    auto path =
        std::filesystem::temp_directory_path() /
        ("fastbdf-writer-edge-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".edf");
    RecordingInfo info;
    SignalSpec signal;
    signal.label = "x";
    signal.samples_per_record = 3;
    signal.digital_min = -100;
    signal.digital_max = 100;
    info.signals = {signal, signal};
    for (auto position : {AnnotationPosition::start, AnnotationPosition::middle, AnnotationPosition::end}) {
        WriterOptions options;
        options.annotation_position = position;
        auto w = Writer::create(path, info, options);
        std::int32_t digital[] = {-101, 0, 101};
        w.write_samples(digital);
        throws([&] { w.finish(); });
        double bad[] = {0, std::numeric_limits<double>::quiet_NaN(), 1};
        throws([&] { w.write_samples(bad); });
        double physical[] = {-1e300, 0, 1e300};
        w.write_samples(physical);
        w.add_annotation({std::chrono::nanoseconds{1}, {}, "late \xc3\xa9"});
        w.flush();
        w.finish();
        w.finish();
        throws([&] { w.flush(); });
        throws([&] { w.write_digital_record(digital); });
        auto bytes = contents(path);
        check(bytes.substr(236, 8) == "1       ");
        auto first = position == AnnotationPosition::start ? 1264u : 1024u;
        auto second = position == AnnotationPosition::end ? 1030u : 1270u;
        check(sample(bytes, first) == -100 && sample(bytes, first + 4) == 100);
        check(sample(bytes, second) == -100 && sample(bytes, second + 4) == 100);
        check(bytes.find("late \xc3\xa9") != std::string::npos);
    }
    {
        auto w = Writer::create(path, info);
        throws(
            [&] { w.add_annotation({std::chrono::nanoseconds{}, std::chrono::nanoseconds{-1}, "negative duration"}); });
        throws([&] { w.finish(); });
        check(!w.finished());
    }
    {
        for (auto format : {Format::edf, Format::bdf, Format::edf_plus, Format::bdf_plus}) {
            auto empty = info;
            empty.format = format;
            auto w = Writer::create(path, empty);
            throws([&] { w.finish(); });
            std::int32_t data[6]{};
            w.write_digital_record(data);
            w.finish();
        }
    }
    {
        auto reserved = info;
        reserved.signals[0].label = "EDF Annotations ";
        throws([&] { auto w = Writer::create(path, reserved); });
        reserved.signals[0].label = " EDF Annotations";
        throws([&] { auto w = Writer::create(path, reserved); });
        reserved.format = Format::bdf_plus;
        reserved.signals[0].label = " BDF Annotations";
        throws([&] { auto w = Writer::create(path, reserved); });
        WriterOptions opts;
        opts.annotation_bytes_per_signal = 600;
        auto w = Writer::create(path, info, opts);
        std::int32_t data[6]{};
        w.write_digital_record(data);
        w.add_annotation({std::chrono::nanoseconds{}, {}, std::string(512, 'a')});
        w.finish();
        check(contents(path).find(std::string(40, 'a') + '\x14') != std::string::npos);
    }
    {
        auto young = info;
        young.birthdate = "0001-01-01";
        auto w = Writer::create(path, young);
        std::int32_t data[6]{};
        w.write_digital_record(data);
        w.finish();
        check(contents(path).find("01-JAN-0001") != std::string::npos);
    }
    {
        RecordingInfo only;
        auto w = Writer::create(path, only);
        for (int i = 0; i < 20; ++i)
            w.add_annotation({std::chrono::nanoseconds{i}, {}, "latin\xe9"});
        w.finish();
        check(w.records_written() > 0);
        check(contents(path).find("latin\xe9") != std::string::npos);
    }
    {
    // A 32-byte channel has room for 22 description bytes after the
    // five-byte timekeeping TAL and five bytes of event framing.
        WriterOptions opts;
        opts.annotation_bytes_per_signal = 32;
        opts.annotation_signals = 2;
        RecordingInfo only;
        auto w = Writer::create(path, only, opts);
        w.add_annotation({std::chrono::nanoseconds{}, {}, std::string(22, 'a')});
        w.add_annotation({std::chrono::nanoseconds{}, {}, std::string(27, 'b')});
        w.add_annotation({std::chrono::nanoseconds{}, {}, "c"});
        w.finish();
        check(w.records_written() == 2);
        const auto bytes = contents(path);
        check(bytes.substr(768 + 8, 22) == std::string(22, 'a'));
        check(bytes.substr(768 + 32 + 3, 27) == std::string(27, 'b'));
        check(bytes[768 + 64 + 8] == 'c');
    }
    {
        WriterOptions opts;
        opts.annotation_bytes_per_signal = 32;
        auto w = Writer::create(path, info, opts);
        std::int32_t data[6]{};
        w.write_digital_record(data);
        w.add_annotation({std::chrono::nanoseconds{}, {}, std::string(22, 'a')});
        w.add_annotation({std::chrono::nanoseconds{}, {}, "later"});
        throws([&] { w.finish(); });
        check(w.records_written() == 1);
        w.write_digital_record(data);
        w.finish();
        check(w.records_written() == 2);
        check(contents(path).find("later") != std::string::npos);
    }
    {
        WriterOptions opts;
        opts.annotation_bytes_per_signal = 32;
        auto w = Writer::create(path, RecordingInfo{}, opts);
        w.add_annotation({std::chrono::nanoseconds{}, {}, std::string(23, 'a')});
        throws([&] { w.finish(); });
        check(w.records_written() == 0);
    }
    {
        WriterOptions opts;
        opts.annotation_bytes_per_signal = std::numeric_limits<std::size_t>::max() - 1;
        throws([&] { auto w = Writer::create(path, info, opts); });
    }
    {
        RecordingInfo tiny = info;
        tiny.signals[0].physical_min = 1e-9;
        tiny.signals[0].physical_max = 2e-9;
        auto w = Writer::create(path, tiny);
        std::int32_t data[6]{};
        w.write_digital_record(data);
        w.finish();
    }
    {
        RecordingInfo huge = info;
        for (auto& s : huge.signals) {
            s.physical_min = -1e308;
            s.physical_max = 1e308;
        }
        auto w = Writer::create(path, huge);
        double data[] = {-1e308, 0, 1e308, -1e308, 0, 1e308};
        w.write_physical_record(data);
        w.finish();
        auto b = contents(path);
        check(sample(b, 1024) == -100 && sample(b, 1026) == 0 && sample(b, 1028) == 100);
    }
    {
        RecordingInfo only;
        only.record_duration = std::chrono::nanoseconds{0};
        auto w = Writer::create(path, only);
        for (int i = 0; i < 30; ++i)
            w.add_annotation({std::chrono::nanoseconds{i}, std::chrono::nanoseconds{0}, "zero duration"});
        w.finish();
        check(w.records_written() > 1);
        check(contents(path).substr(244, 8) == "0       ");
    }
    {
        auto bdf = info;
        bdf.format = Format::bdf_plus;
        auto w = Writer::create(path, bdf);
        const std::byte packed[] = {std::byte{0xff}, std::byte{0xff}, std::byte{0xff}, std::byte{0}, std::byte{0},
                                    std::byte{0},    std::byte{1},    std::byte{0},    std::byte{0}, std::byte{2},
                                    std::byte{0},    std::byte{0},    std::byte{3},    std::byte{0}, std::byte{0},
                                    std::byte{4},    std::byte{0},    std::byte{0}};
        w.write_packed_record(packed);
        std::int16_t short_data[] = {-101, 0, 101, -101, 0, 101};
        w.write_short_record(short_data);
        w.finish();
        const auto b = contents(path);
        check(b.substr(1024, sizeof(packed)) == std::string(reinterpret_cast<const char*>(packed), sizeof(packed)));
        const auto second = 1024 + 18 + 240;
        check(static_cast<unsigned char>(b[second]) == 156 && static_cast<unsigned char>(b[second + 1]) == 255 &&
              static_cast<unsigned char>(b[second + 2]) == 255);
    }
    std::filesystem::remove(path);
}
