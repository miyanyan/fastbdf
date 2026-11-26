#include <algorithm>
#include <array>
#include <chrono>
#include <fastbdf/fastbdf.hpp>
#include <filesystem>
#include <fstream>
#include <future>
#include <iostream>
#include <list>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

using namespace std::chrono_literals;
using namespace fastbdf;
namespace fs = std::filesystem;
std::int64_t nanoseconds_in_one_second();
namespace {
template <class R>
concept CanReadSamples = requires(Reader &reader, Reader::Cursor &cursor, R &&out) {
  reader.read_at(0, 0, std::forward<R>(out));
  reader.read(0, std::forward<R>(out));
  cursor.read(std::forward<R>(out));
};
static_assert(CanReadSamples<std::int32_t (&)[2]> &&
              CanReadSamples<std::vector<double> &> &&
              CanReadSamples<std::span<double, 2>>);
static_assert(!CanReadSamples<std::span<const double>> &&
              !CanReadSamples<std::vector<double>> &&
              !CanReadSamples<std::array<float, 2> &>);

template <class R>
concept CanWriteSamples = requires(Writer &w, const R &samples) {
  w.write_samples(samples);
  w.write_samples(0, samples);
};

static_assert(CanWriteSamples<std::int16_t[2]> &&
              CanWriteSamples<std::array<std::int32_t, 2>> &&
              CanWriteSamples<std::vector<float>> &&
              CanWriteSamples<std::span<const double>>);
static_assert(!CanWriteSamples<std::array<unsigned, 2>> &&
              !CanWriteSamples<std::array<bool, 2>> &&
              !CanWriteSamples<std::vector<std::int64_t>> &&
              !CanWriteSamples<std::list<float>>);

void check(bool condition, const char *message) {
  if (!condition)
    throw std::runtime_error(message);
}
template <class F> void expect_error(F action) {
  try {
    action();
  } catch (const Error &) {
    return;
  }
  throw std::runtime_error("expected a fastbdf::Error");
}
struct Temp {
  fs::path dir =
      fs::temp_directory_path() /
      ("fastbdf-native-" +
       std::to_string(
           std::chrono::steady_clock::now().time_since_epoch().count()));
  Temp() {
    if (!fs::create_directory(dir))
      throw std::runtime_error("temporary directory unavailable");
  }
  ~Temp() {
    std::error_code ignored;
    fs::remove_all(dir, ignored);
  }
};

// Independent bytes exercise classic EDF/BDF and samples outside the declared
// digital range.
void classic_fixture(const fs::path &path, bool bdf) {
  auto field = [](std::string s, unsigned n) {
    s.resize(n, ' ');
    return s;
  };
  std::string h = bdf ? std::string("\xff"
                                    "BIOSEMI",
                                    8)
                      : "0       ";
  h += field("patient-42", 80) + field("night", 80) + "31.12.24" + "23.59.58";
  h += field("768", 8) + field("", 44) + field("3", 8) + field("0.5", 8) +
       field("2", 4);
  for (auto [a, b, n] :
       std::initializer_list<std::tuple<std::string, std::string, unsigned>>{
           {"EEG", "Temp", 16},
           {"AgCl", "Therm", 80},
           {"uV", "C", 8},
           {"10", "-100", 8},
           {"-10", "100", 8},
           {"-10", "-100", 8},
           {"+10", "+100", 8},
           {"HP:0.1Hz", "", 80},
           {"3", "1", 8},
           {"", "", 32}})
    h += field(a, n) + field(b, n);
  std::ofstream f(path, std::ios::binary);
  f << h;
  for (int value : {-100, -10, 0, 20, 1, 9, 100, 21, -3, -2, -1, 22}) {
    auto bits = static_cast<unsigned>(value);
    for (unsigned byte = 0; byte < (bdf ? 3u : 2u); ++byte)
      f.put(static_cast<char>(bits >> (8 * byte)));
  }
}
void test_classic(const fs::path &p, bool bdf) {
  classic_fixture(p, bdf);
  auto file = Reader::open(p);
  check(file.record_count() == 3 && file.sample_rate(0) == 6 &&
            file.sample_rate(1) == 2,
        "mixed sample rates");
  check(file.header().patient == "patient-42" &&
            file.signals()[0].transducer == "AgCl",
        "classic metadata");
  std::array<std::int32_t, 12> x{};
  check(file.read_at(0, 0, x) == 9 && x[0] == -10 && x[5] == 10,
        "digital clipping and EOF");
  std::array<double, 3> phys{};
  check(file.read_at(0, 1, phys) == 3 &&
            phys == std::array<double, 3>{10, 0, -1},
        "negative gain");
  const auto t = file.sample_time(0, 4); // 0.5 s + 1/6 s
  check(t.whole_nanoseconds == 666'666'666ns && t.numerator * 3 == 2 * t.denominator,
        "exact fractional nanosecond for a sample in the second record");
  check(file.sample_index_at_or_after(0, t.whole_nanoseconds) == 4 &&
            file.sample_index_at_or_after(0, t.whole_nanoseconds + 1ns) == 5,
        "lookup on either side of a fractional sample time");
  expect_error([&] { file.sample_time(0, file.sample_count(0)); });
  auto ids = file.select_unique({"Temp", "EEG"});
  auto data = file.read_digital(ids, {250ms, 1s});
  check(data == std::vector<std::vector<std::int32_t>>{{21}, {0, 1, 9, 10}},
        "half-open mixed-rate window");
  check(file.read_digital(ids, {500ms, 500ms}) ==
            std::vector<std::vector<std::int32_t>>(2),
        "empty window");
  auto one = file.cursor(0), two = file.cursor(0);
  one.seek(8);
  check(one.read(x) == 1 && one.tell() == 9 && two.tell() == 0,
        "independent cursors");
  check(one.read(x) == 0, "cursor EOF");
  std::array<double, 1> physical_sample{};
  check(two.read(physical_sample) == 1 && physical_sample[0] == 10 &&
            two.tell() == 1 && file.tell(0) == 0, "independent physical cursor");
  check(file.read(0, physical_sample) == 1 && physical_sample[0] == 10 &&
            file.tell(0) == 1, "implicit physical cursor");
  expect_error([&] { file.select_unique("absent"); });
  expect_error([&] { file.read_at(2, 0, x); });
  expect_error([&] { file.read_digital(ids, {2s, 1s}); });
  auto background = std::async(std::launch::async, [&] {
    std::array<std::int32_t, 6> out{};
    for (int i = 0; i < 256; i++)
      check(file.read_at(0, 1, out) == 6 && out[5] == -3,
            "concurrent EEG read");
  });
  for (int i = 0; i < 256; i++)
    check(file.read_at(1, 2, x) == 1 && x[0] == 22,
          "concurrent Temp read");
  background.get();
}
void test_malformed(const fs::path &p) {
  classic_fixture(p, false);
  fs::resize_file(p, fs::file_size(p) - 1);
  expect_error([&] { auto bad = Reader::open(p); });
  classic_fixture(p, false);
  {
    std::fstream f(p, std::ios::binary | std::ios::in | std::ios::out);
    f.seekp(236);
    f << "10      ";
  }
  expect_error([&] { auto bad = Reader::open(p); });
}
void test_long_recording(const fs::path &path) {
  RecordingInfo info;
  info.format = Format::edf;
  info.record_duration = 99999999s;
  info.add_signal("x", 1);
  const auto last_in_range = (std::chrono::nanoseconds::max)().count() / info.record_duration.count();
  auto writer = Writer::create(path, info);
  for (std::int64_t r = 0; r <= last_in_range + 1; ++r)
    writer.write_digital_record(std::array{static_cast<std::int32_t>(r)});
  writer.finish();

  auto reader = Reader::open(path);
  auto start = reader.record_start(last_in_range);
  check(start == info.record_duration * last_in_range, "last representable start");
  expect_error([&] { reader.record_start(last_in_range + 1); });
  check(reader.read_digital(std::array<std::size_t, 1>{0}, {start, (std::chrono::nanoseconds::max)()}) ==
            std::vector<std::vector<std::int32_t>>{{static_cast<std::int32_t>(last_in_range)}},
        "window at nanosecond limit");
  std::array<std::int32_t, 1> sample{};
  check(reader.read_at(0, last_in_range + 1, sample) == 1 &&
            sample[0] == last_in_range + 1, "sample beyond nanosecond limit");
}
void test_annotation_nanoseconds(const fs::path &path) {
  using std::chrono::nanoseconds;
  const auto maximum = (nanoseconds::max)();
  const auto minimum = (nanoseconds::min)();
  auto writer = Writer::create(path, RecordingInfo{});
  writer.add_annotation({minimum, maximum, "limits"});
  writer.add_annotation({-1ns, 1ns, "one nanosecond"});
  writer.finish();

  auto reader = Reader::open(path);
  const auto &events = reader.annotations();
  check(events.size() == 2 && events[0].onset == minimum &&
            events[0].duration == maximum && events[1].onset == -1ns &&
            events[1].duration == 1ns, "signed nanosecond annotation round-trip");
}
void test_sample_offset_overflow(const fs::path &path) {
  RecordingInfo info;
  info.format = Format::edf;
  info.record_duration = 90'000'000s;
  constexpr unsigned samples = 25'005;
  info.add_signal("x", samples);
  auto writer = Writer::create(path, info);
  writer.write_short_record(std::vector<std::int16_t>(samples));
  writer.finish();

  auto reader = Reader::open(path);
  // At 14/15 of this record, elapsed nanoseconds times samples exceeds UINT64_MAX.
  const std::chrono::nanoseconds boundary = 84'000'000s;
  constexpr auto index = samples / 15 * 14;
  const auto t = reader.sample_time(0, index);
  check(t.whole_nanoseconds == boundary && t.numerator == 0,
        "sample time when index times duration overflows");
  check(reader.sample_index_at_or_after(0, boundary - 1ns) == index &&
            reader.sample_index_at_or_after(0, boundary) == index &&
            reader.sample_index_at_or_after(0, boundary + 1ns) == index + 1 &&
            reader.sample_index_at_or_after(0, info.record_duration - 1ns) == samples,
        "exact sample boundaries when the intermediate product overflows");
}
void test_write_samples(const fs::path &path) {
  RecordingInfo info;
  for (const char *label : {"A", "B"}) {
    auto &s = info.add_signal(label, 2);
    s.digital_min = -100;
    s.digital_max = 100;
  }
  auto w = Writer::create(path, info);
  const std::int16_t first[] = {-20, 30};
  expect_error([&] { w.write_samples(1, first); });
  w.write_samples(first);
  w.write_samples(1, std::array<std::int32_t, 2>{-200, 200});
  const std::vector<float> physical{-0.5f, 0.5f};
  w.write_samples(0, physical);
  const double second[] = {0.5, -0.5};
  w.write_samples(std::span{second});
  w.finish();
  auto reader = Reader::open(path);
  check(reader.read_digital(reader.select_unique({"A", "B"}), {0s, 2s}) ==
            std::vector<std::vector<std::int32_t>>{{-20, 30, -50, 50},
                                                   {-100, 100, 50, -50}},
        "typed channel writes");
}
} // namespace
int main() {
  try {
    check(nanoseconds_in_one_second() == 1'000'000'000, "header-only cross-TU linkage");
    Temp tmp;
    test_classic(tmp.dir / "classic.edf", false);
    test_classic(tmp.dir / "classic.bdf", true);
    test_malformed(tmp.dir / "bad.edf");
    test_long_recording(tmp.dir / "long.edf");
    test_annotation_nanoseconds(tmp.dir / "annotation-ns.edf");
    test_sample_offset_overflow(tmp.dir / "sample-offset.edf");
    test_write_samples(tmp.dir / "typed.edf");
    RecordingInfo invalid_header;
    invalid_header.add_signal("Andr\xe9", 1);
    expect_error(
        [&] { Writer::create(tmp.dir / "non-ascii.edf", invalid_header); });
    RecordingInfo info;
    info.format = Format::bdf_plus;
    auto &lead = info.add_signal("Fp1", 2);
    lead.physical_min = -1e20;
    lead.physical_max = 1;
    check(lead.digital_max == 8388607, "BDF signal defaults");
    auto path = tmp.dir / "endpoints.bdf";
    const std::array<std::int32_t, 2> endpoints{-8388608, 8388607};
    {
      auto writer = Writer::create(path, info);
      writer.write_digital_record(endpoints);
      writer.add_annotation({std::chrono::nanoseconds{}, {}, "a\x7f" "b"});
      writer.finish();
    }
    auto reader = Reader::open(path);
    std::array<std::int32_t, 2> data{};
    check(reader.read_at(0, 0, data) == data.size() && data == endpoints,
          "BDF digital endpoint round-trip");
    std::array<double, 2> physical{};
    check(reader.read_at(0, 0, physical) == physical.size() &&
              physical == std::array<double, 2>{-1e20, 1},
          "asymmetric physical endpoints");
    check(reader.annotations().at(0).text == "a.b", "DEL annotation round-trip");
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
