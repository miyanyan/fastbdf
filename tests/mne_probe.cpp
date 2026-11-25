#include <fastbdf/fastbdf.hpp>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <vector>
using namespace fastbdf;
namespace {
void generate(const std::filesystem::path &dir) {
  std::filesystem::create_directories(dir);
  const char *names[] = {"edf", "bdf", "edf_plus", "bdf_plus"};
  const Format formats[] = {Format::edf, Format::bdf, Format::edf_plus,
                            Format::bdf_plus};
  for (int f = 0; f < 4; ++f) {
    RecordingInfo h;
    h.format = formats[f];
    h.start = std::chrono::local_days{std::chrono::year{2000} / 1 / 1};
    h.record_duration = Tick{5000000};
    for (auto n : {16u, 7u}) {
      auto &s = h.add_signal("channel" + std::to_string(h.signals.size()), n);
      s.unit = h.signals.size() == 1 ? "uV" : "V";
      s.physical_min = h.signals.size() == 1 ? -100 : -0.1;
      s.physical_max = h.signals.size() == 1 ? 100 : 0.1;
    }
    auto w = Writer::create(
        dir / (std::string(names[f]) + (f % 2 ? ".bdf" : ".edf")), h);
    if (f >= 2) {
      w.add_annotation({Tick{2500000}, Tick{1250000}, "event \xce\xb1"});
      w.add_annotation({Tick{12500000}, {}, "end"});
    }
    for (int r = 0; r < 4; ++r) {
      std::vector<std::int32_t> values;
      for (const auto &s : h.signals)
        for (unsigned i = 0; i < s.samples_per_record; ++i)
          values.push_back(i % 5 == 0   ? s.digital_min
                           : i % 5 == 1 ? s.digital_max
                                        : (r * 31 + int(i) - 40));
      w.write_digital_record(values);
    }
    w.finish();
  }
}
void array(const std::vector<double> &v) {
  std::cout << '[';
  for (std::size_t i = 0; i < v.size(); ++i) {
    if (i)
      std::cout << ',';
    std::cout << v[i];
  }
  std::cout << ']';
}
std::string json_string(std::string_view value) {
  std::string result = "\"";
  const char *hex = "0123456789abcdef";
  for (unsigned char c : value) {
    if (c == '"' || c == '\\') {
      result += '\\';
      result += char(c);
    } else if (c < 32) {
      result += "\\u00";
      result += hex[c >> 4];
      result += hex[c & 15];
    } else
      result += char(c);
  }
  return result + "\"";
}
void inspect(const std::filesystem::path &path, bool bounded) {
  Reader r(path);
  std::cout << std::setprecision(17) << "{\"records\":" << r.record_count()
            << ",\"record_duration\":"
            << r.header().record_duration.count() / 1e7 << ",\"channels\":[";
  for (std::size_t id = 0; id < r.signals().size(); ++id) {
    if (id)
      std::cout << ',';
    const auto &s = r.signals()[id];
    const auto count = r.sample_count(id);
    std::vector<double> full(bounded ? 0 : count);
    if (!bounded && r.read_physical_at(id, 0, full) != full.size())
      throw std::runtime_error("short full read");
    std::cout << "{\"name\":" << json_string(s.label)
              << ",\"unit\":" << json_string(s.unit) << ",\"rate\":"
              << s.samples_per_record * 1e7 / r.header().record_duration.count()
              << ",\"count\":" << count << ",\"values\":";
    array(full);
    std::cout << ",\"windows\":[";
    const auto width =
        (std::min)(count, bounded ? std::uint64_t{256} : std::uint64_t{5});
    const std::uint64_t starts[] = {0, (count - width) / 2, count - width};
    for (int k = 0; k < 3; ++k) {
      std::vector<double> samples(static_cast<std::size_t>(width));
      if (r.read_physical_at(id, starts[k], samples) != samples.size())
        throw std::runtime_error("short window read");
      if (k)
        std::cout << ',';
      std::cout << "{\"first\":" << starts[k] << ",\"values\":";
      array(samples);
      std::cout << '}';
    }
    std::cout << "]}";
  }
  std::cout << "],\"annotations\":[";
  bool comma = false;
  for (const auto &a : r.annotations()) {
    if (comma)
      std::cout << ',';
    comma = true;
    std::cout << "{\"onset\":" << a.onset.count() / 1e7 << ",\"duration\":"
              << (a.duration ? a.duration->count() / 1e7 : 0)
              << ",\"text\":" << json_string(a.text) << '}';
  }
  std::cout << "]}\n";
}
} // namespace
int main(int argc, char **argv) {
  try {
    if (argc != 3)
      throw std::runtime_error("usage: mne_probe generate <directory> | "
                               "inspect <file> | windows <file>");
    if (std::string(argv[1]) == "generate")
      generate(argv[2]);
    else if (std::string(argv[1]) == "inspect" ||
             std::string(argv[1]) == "windows")
      inspect(argv[2], std::string(argv[1]) == "windows");
    else
      throw std::runtime_error("unknown command");
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
