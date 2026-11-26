#pragma once
#include <chrono>
#include <cstdint>
#include <ctime>
#include <optional>
#include <ratio>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace fastbdf {
enum class Format { edf, bdf, edf_plus, bdf_plus };
using Tick = std::chrono::duration<std::int64_t, std::ratio<1, 10'000'000>>;
struct TimeRange {
  Tick begin{}, end{};
};

struct SignalSpec {
  std::string label, unit, transducer, prefilter;
  double physical_min = -1, physical_max = 1;
  std::int32_t digital_min = -32768, digital_max = 32767;
  std::uint32_t samples_per_record = 1;
};

enum class Errc { invalid_argument, invalid_file, io, invalid_state };

class Error : public std::runtime_error {
public:
  explicit Error(std::string message, Errc code = Errc::invalid_argument)
      : std::runtime_error(std::move(message)), code_(code) {}
  Errc code() const noexcept { return code_; }

private:
  Errc code_;
};

namespace detail {
inline std::string_view trim_spaces(std::string_view s) {
  while (!s.empty() && s.front() == ' ')
    s.remove_prefix(1);
  while (!s.empty() && s.back() == ' ')
    s.remove_suffix(1);
  return s;
}

inline std::chrono::local_seconds local_now() {
  const auto now = std::time(nullptr);
  std::tm tm{};
#ifdef _WIN32
  if (localtime_s(&tm, &now) != 0)
#else
  if (!localtime_r(&now, &tm))
#endif
    throw Error("Cannot determine local recording time", Errc::io);
  using namespace std::chrono;
  return local_days{year{tm.tm_year + 1900} /
                    month{static_cast<unsigned>(tm.tm_mon + 1)} /
                    day{static_cast<unsigned>(tm.tm_mday)}} +
         hours{tm.tm_hour} + minutes{tm.tm_min} + seconds{tm.tm_sec};
}
} // namespace detail

struct RecordingInfo {
  Format format = Format::edf_plus;
  std::chrono::local_seconds start = detail::local_now();
  Tick subsecond{}, record_duration = std::chrono::seconds{1};
  std::string patient, recording;
  std::string patient_code, patient_name, sex, birthdate, patient_additional;
  std::string admin_code, technician, equipment, recording_additional;
  std::vector<SignalSpec> signals;

  // Set format first. The full digital range follows the sample width.
  SignalSpec &add_signal(std::string label, std::uint32_t samples_per_record) {
    if (!samples_per_record)
      throw Error("A signal needs at least one sample per record");
    SignalSpec signal;
    signal.label = std::move(label);
    signal.samples_per_record = samples_per_record;
    if (format == Format::bdf || format == Format::bdf_plus) {
      signal.digital_min = -8388608;
      signal.digital_max = 8388607;
    }
    return signals.emplace_back(std::move(signal));
  }
};

struct Annotation {
  Tick onset{};
  std::optional<Tick> duration;
  std::string text;
};

} // namespace fastbdf
