#include <algorithm>
#include <array>
#include <chrono>
#include <fastbdf/fastbdf.hpp>
#include <filesystem>
#include <fstream>
#include <future>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

using namespace std::chrono_literals;
using namespace fastbdf;
namespace fs = std::filesystem;
std::int64_t ticks_in_one_second();
namespace {
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
  check(file.read_digital_at(0, 0, x) == 9 && x[0] == -10 && x[5] == 10,
        "digital clipping and EOF");
  std::array<double, 3> phys{};
  check(file.read_physical_at(0, 1, phys) == 3 &&
            phys == std::array<double, 3>{10, 0, -1},
        "negative gain");
  auto ids = file.select_unique({"Temp", "EEG"});
  auto data = file.read_digital(ids, {250ms, 1s});
  check(data == std::vector<std::vector<std::int32_t>>{{21}, {0, 1, 9, 10}},
        "half-open mixed-rate window");
  check(file.read_digital(ids, {500ms, 500ms}) ==
            std::vector<std::vector<std::int32_t>>(2),
        "empty window");
  auto one = file.cursor(0), two = file.cursor(0);
  one.seek(8);
  check(one.read_digital(x) == 1 && one.tell() == 9 && two.tell() == 0,
        "independent cursors");
  check(one.read_digital(x) == 0, "cursor EOF");
  expect_error([&] { file.select_unique("absent"); });
  expect_error([&] { file.read_digital_at(2, 0, x); });
  expect_error([&] { file.read_digital(ids, {2s, 1s}); });
  auto background = std::async(std::launch::async, [&] {
    std::array<std::int32_t, 6> out{};
    for (int i = 0; i < 256; i++)
      check(file.read_digital_at(0, 1, out) == 6 && out[5] == -3,
            "concurrent EEG read");
  });
  for (int i = 0; i < 256; i++)
    check(file.read_digital_at(1, 2, x) == 1 && x[0] == 22,
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
} // namespace
int main() {
  try {
    check(ticks_in_one_second() == 10000000, "header-only cross-TU linkage");
    Temp tmp;
    test_classic(tmp.dir / "classic.edf", false);
    test_classic(tmp.dir / "classic.bdf", true);
    test_malformed(tmp.dir / "bad.edf");
    RecordingInfo invalid_header;
    invalid_header.add_signal("Andr\xe9", 1);
    expect_error(
        [&] { Writer::create(tmp.dir / "non-ascii.edf", invalid_header); });
    RecordingInfo info;
    info.format = Format::bdf_plus;
    auto &lead = info.add_signal("Fp1", 2);
    check(lead.digital_max == 8388607, "BDF signal defaults");
    auto path = tmp.dir / "endpoints.bdf";
    const std::array<std::int32_t, 2> endpoints{-8388608, 8388607};
    {
      auto writer = Writer::create(path, info);
      writer.write_digital_record(endpoints);
      writer.finish();
    }
    auto reader = Reader::open(path);
    std::array<std::int32_t, 2> data{};
    check(reader.read_digital_at(0, 0, data) == data.size() && data == endpoints,
          "BDF digital endpoint round-trip");
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
