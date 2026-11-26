#include <array>
#include <cmath>
#include <fastbdf/reader.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <tuple>
using namespace fastbdf;
using namespace std::chrono_literals;
namespace fs = std::filesystem;
namespace {
void check(bool ok, const char *message) {
  if (!ok)
    throw std::runtime_error(message);
}
template <class F> void invalid(F action) {
  try {
    action();
  } catch (const Error &e) {
    check(e.code() == Errc::invalid_file, "invalid file error category");
    return;
  }
  throw std::runtime_error("accepted malformed file");
}
template <class F> void bad_seek(F action) {
  try {
    action();
  } catch (const Error &e) {
    check(e.code() == Errc::invalid_argument, "seek error category");
    return;
  }
  throw std::runtime_error("accepted invalid seek");
}
std::string field(std::string value, unsigned width) {
  value.resize(width, ' ');
  return value;
}
void fixture(const fs::path &path, bool discontinuous = true,
             bool extreme = false) {
  std::string h = "0       ";
  h += field("P M 01-JAN-2000 N", 80) +
       field("Startdate 01-JAN-2024 X X X", 80) + "01.01.24" + "00.00.00";
  h += field("1024", 8) + field(discontinuous ? "EDF+D" : "EDF+C", 44) +
       field("2", 8) + field("1", 8) + field("3", 4);
  for (auto [a, b, c, n] : std::initializer_list<
           std::tuple<std::string, std::string, std::string, unsigned>>{
           {"A", "B", "EDF Annotations", 16},
           {"", "", "", 80},
           {"uV", "C", "", 8},
           {extreme ? "-1e308" : "-1", "-1", "-1", 8},
           {extreme ? "1e308" : "1", "1", "1", 8},
           {"-10", "-10", "-32768", 8},
           {"10", "10", "32767", 8},
           {"", "", "", 80},
           {"2", "1", "32", 8},
           {"", "", "", 32}})
    h += field(a, n) + field(b, n) + field(c, n);
  std::ofstream out(path, std::ios::binary);
  out << h;
  for (int r = 0; r < 2; ++r) {
    for (int value : (extreme ? std::array<int, 3>{r ? 10 : -10, 0, 0}
                              : std::array<int, 3>{-10, 10, 0})) {
      out.put(static_cast<char>(value));
      out.put(static_cast<char>(static_cast<unsigned>(value) >> 8));
    }
    auto onset = r ? (discontinuous ? "+2.125" : "+1.125") : "+0.125";
    std::string tal = onset;
    tal += '\x14';
    tal += '\x14';
    tal += '\0';
    if (!r) {
      tal += "+0.625";
      tal += '\x15';
      tal += "0.25";
      tal += '\x14';
      tal += "\xc3\xa9";
      tal += '\x14';
      tal += '\0';
    }
    tal.resize(64, '\0');
    out.write(tal.data(), tal.size());
  }
}
void patch(const fs::path &path, std::uint64_t at, std::string_view value) {
  std::fstream f(path, std::ios::binary | std::ios::in | std::ios::out);
  f.seekp(static_cast<std::streamoff>(at));
  f.write(value.data(), static_cast<std::streamsize>(value.size()));
}
} // namespace
int main() {
  auto path =
      fs::temp_directory_path() /
      ("fastbdf-reader-" +
       std::to_string(
           std::chrono::steady_clock::now().time_since_epoch().count()) +
       ".edf");
  struct Cleanup {
    fs::path p;
    ~Cleanup() {
      std::error_code e;
      fs::remove(p, e);
    }
  } cleanup{path};
  try {
    fixture(path);
    auto r = Reader::open(path);
    auto ids = r.select_unique({"A", "B"});
    check(r.discontinuous() && r.record_onsets().size() == 2,
          "discontinuous metadata");
    check(r.header().birthdate == "2000-01-01" && r.header().subsecond == 125ms,
          "structured metadata");
    check(r.read_digital(ids, {500ms, 2500ms}) ==
              std::vector<std::vector<std::int32_t>>{{10, -10}, {0}},
          "discontinuous window");
    const auto after_gap = r.sample_time(ids[0], 2);
    check(after_gap.whole_nanoseconds == 2s && after_gap.numerator == 0 &&
              r.sample_index_at_or_after(ids[0], 1500ms) == 2,
          "sample time includes the discontinuity");
    check(r.read_digital(ids, {1s, 2s}) ==
              std::vector<std::vector<std::int32_t>>(2),
          "gap window");
    const auto &a = r.annotations();
    check(a.size() == 1 && a[0].onset == 500ms && a[0].duration == 250ms &&
              a[0].text == "\xc3\xa9",
          "lazy UTF8 annotations");
    std::vector<std::vector<std::int32_t>> reuse;
    r.read_digital_into(ids, {0s, 3s}, reuse);
    auto capacity = reuse[0].capacity();
    r.read_digital_into(ids, {0s, 500ms}, reuse);
    check(reuse[0] == std::vector<std::int32_t>{-10} &&
              reuse[0].capacity() == capacity,
          "reused window capacity");
    auto cursor = r.cursor(0);
    cursor.seek(-1, std::ios_base::end);
    check(cursor.tell() == 3, "end seek");
    cursor.seek(-2, std::ios_base::cur);
    check(cursor.tell() == 1, "relative seek");
    cursor.rewind();
    check(cursor.tell() == 0, "rewind");
    bad_seek([&] {
      cursor.seek(std::numeric_limits<std::int64_t>::min(), std::ios_base::cur);
    });
    bad_seek([&] { cursor.seek(1, std::ios_base::end); });
    fixture(path, false, true);
    auto extreme = Reader::open(path);
    std::array<double, 4> physical{};
    extreme.read_at(0, 0, physical);
    check(physical == std::array<double, 4>{-1e308, 0, 1e308, 0},
          "stable extreme interpolation");
    fixture(path);
    fs::resize_file(path, fs::file_size(path) - 1);
    invalid([&] { Reader::open(path); });
    fixture(path);
    patch(path, 168, "32.01.24");
    invalid([&] { Reader::open(path); });
    fixture(path, false);
    patch(path, 1024 + 70 + 6, "x");
    auto lazy = Reader::open(path);
    invalid([&] { lazy.annotations(); });
    invalid([&] { Reader::open(path, AnnotationMode::all); });
    fixture(path);
    patch(path, 1024 + 6 + 9 + 6 + 1 + 4 + 1, std::string_view("\xc0\xaf", 2));
    auto utf = Reader::open(path);
    check(utf.annotations().at(0).text == "\xc0\xaf",
          "preserve raw annotation bytes");
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
