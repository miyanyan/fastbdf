#pragma once

#include "text.hpp"
#include "types.hpp"
#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace fastbdf {

enum class AnnotationPosition { start, middle, end };

struct WriterOptions {
  std::size_t annotation_signals = 1;
  std::size_t annotation_bytes_per_signal = 240;
  AnnotationPosition annotation_position = AnnotationPosition::end;
};

// Metadata is supplied before creation. Samples are signal-major within each
// record. finish() is required: it validates and writes late annotations and
// the record count.
class Writer {
public:
  static Writer create(const std::filesystem::path &path, RecordingInfo info,
                       WriterOptions options = {}) {
    return Writer(path, std::move(info), options);
  }
  Writer(const Writer &) = delete;
  Writer &operator=(const Writer &) = delete;
  Writer(Writer &&) noexcept = default;
  Writer &operator=(Writer &&) noexcept = default;
  ~Writer() = default;

  const RecordingInfo &info() const noexcept { return info_; }
  std::uint64_t records_written() const noexcept { return records_; }
  bool finished() const noexcept { return finished_; }

  void write_digital_record(std::span<const std::int32_t> samples) {
    write_record(samples,
                 [this](auto x, const auto &s) { return digital(x, s); });
  }
  void write_short_record(std::span<const std::int16_t> samples) {
    write_record(samples,
                 [this](auto x, const auto &s) { return digital(x, s); });
  }
  void write_physical_record(std::span<const double> samples) {
    write_record(samples,
                 [this](auto x, const auto &s) { return physical(x, s); });
  }
  void write_digital_samples(std::span<const std::int32_t> samples) {
    write_channel(next_signal_, samples,
                  [this](auto x, const auto &s) { return digital(x, s); });
  }
  void write_short_samples(std::span<const std::int16_t> samples) {
    write_channel(next_signal_, samples,
                  [this](auto x, const auto &s) { return digital(x, s); });
  }
  void write_physical_samples(std::span<const double> samples) {
    write_channel(next_signal_, samples,
                  [this](auto x, const auto &s) { return physical(x, s); });
  }
  void write_digital_samples(std::size_t channel,
                             std::span<const std::int32_t> samples) {
    write_channel(channel, samples,
                  [this](auto x, const auto &s) { return digital(x, s); });
  }
  void write_short_samples(std::size_t channel,
                           std::span<const std::int16_t> samples) {
    write_channel(channel, samples,
                  [this](auto x, const auto &s) { return digital(x, s); });
  }
  void write_physical_samples(std::size_t channel,
                              std::span<const double> samples) {
    write_channel(channel, samples,
                  [this](auto x, const auto &s) { return physical(x, s); });
  }
  void write_packed_record(std::span<const std::byte> bytes) {
    ensure_open();
    if (next_signal_ || bytes.size() != data_bytes_)
      fail("packed record size or channel order mismatch");
    emit_record(reinterpret_cast<const char *>(bytes.data()));
  }
  void add_annotation(Annotation annotation) {
    ensure_open();
    if (!plus_)
      fail("annotations require EDF+ or BDF+");
    if (annotation.duration && annotation.duration->count() < 0)
      fail("negative annotation duration");
    if (annotation.duration &&
        seconds(annotation.duration->count(), false).size() > 15)
      fail("annotation duration exceeds the EDFlib precision limit");
    if (annotation.text.size() > 512)
      fail("annotation text exceeds 512 bytes");
    validate_utf8(annotation.text);
    annotations_.push_back(std::move(annotation));
  }
  void add_annotation_latin1(Tick onset, std::optional<Tick> duration,
                             std::string_view text) {
    std::string utf8;
    for (unsigned char c : text) {
      if (c < 128)
        utf8.push_back(static_cast<char>(c));
      else {
        utf8.push_back(static_cast<char>(0xc0 | (c >> 6)));
        utf8.push_back(static_cast<char>(0x80 | (c & 63)));
      }
    }
    add_annotation({onset, duration, std::move(utf8)});
  }
  void flush() {
    ensure_open();
    output_.flush();
    check_io();
  }
  void finish() {
    if (finished_)
      return;
    ensure_open();
    if (next_signal_)
      fail("cannot finish a partially written data record");
    if (info_.signals.empty() && records_ == 0) {
      std::size_t slot = 0, used = timekeeping(0).size();
      for (const auto &a : annotations_) {
        const auto n = annotation_tal(a).size();
        if (n > options_.annotation_bytes_per_signal)
          fail("annotation exceeds channel capacity");
        while (used + n > options_.annotation_bytes_per_signal) {
          ++slot;
          used = slot % options_.annotation_signals == 0
                     ? timekeeping(slot / options_.annotation_signals).size()
                     : 0;
          if (n + used > options_.annotation_bytes_per_signal &&
              options_.annotation_signals == 1)
            fail("annotation exceeds timekeeping channel capacity");
        }
        used += n;
      }
      const auto count = slot / options_.annotation_signals + 1;
      for (std::size_t i = 0; i < count; ++i)
        emit_record(nullptr);
    }
    if (!records_)
      fail("cannot finish without at least one data record");
    // Plan every annotation block first so a capacity failure cannot partially
    // patch a file.
    if (plus_ && !annotations_.empty()) {
      std::map<std::size_t, std::string> blocks;
      if (records_ > (std::numeric_limits<std::size_t>::max)() /
                         options_.annotation_signals)
        fail("too many records");
      const auto capacity =
          static_cast<std::size_t>(records_) * options_.annotation_signals;
      auto block_for = [&](std::size_t slot) -> std::string & {
        auto [it, inserted] = blocks.try_emplace(slot);
        if (inserted && slot % options_.annotation_signals == 0)
          it->second = timekeeping(slot / options_.annotation_signals);
        return it->second;
      };
      std::size_t slot = 0;
      for (const auto &a : annotations_) {
        auto tal = annotation_tal(a);
        if (tal.size() > options_.annotation_bytes_per_signal)
          fail("annotation exceeds channel capacity");
        while (slot < capacity && block_for(slot).size() + tal.size() >
                                      options_.annotation_bytes_per_signal)
          ++slot;
        if (slot == capacity)
          fail("annotation capacity exhausted; increase annotation channels or "
               "bytes before writing");
        block_for(slot) += tal;
      }
      for (auto &[slot, block] : blocks) {
        const auto record = slot / options_.annotation_signals,
                   channel = slot % options_.annotation_signals;
        output_.seekp(static_cast<std::streamoff>(
            header_bytes_ + record * record_bytes_ + annotation_offset_ +
            channel * options_.annotation_bytes_per_signal));
        block.resize(options_.annotation_bytes_per_signal, '\0');
        output_.write(block.data(), static_cast<std::streamsize>(block.size()));
      }
    }
    output_.seekp(236);
    const auto count = field(std::to_string(records_), 8);
    output_.write(count.data(), 8);
    output_.flush();
    check_io();
    output_.close();
    check_io();
    finished_ = true;
  }

private:
  RecordingInfo info_;
  WriterOptions options_;
  std::ofstream output_;
  std::vector<Annotation> annotations_;
  std::vector<char> pending_, encoded_;
  std::size_t sample_bytes_ = 0, data_bytes_ = 0, record_bytes_ = 0,
              header_bytes_ = 0, next_signal_ = 0, annotation_offset_ = 0,
              annotation_index_ = 0;
  std::uint64_t records_ = 0;
  bool plus_ = false, finished_ = false;
  [[noreturn]] static void fail(std::string message) {
    throw Error(std::move(message));
  }
  void ensure_open() const {
    if (finished_ || !output_.is_open())
      throw Error("Writer is closed", Errc::invalid_state);
    check_io();
  }
  void check_io() const {
    if (!output_)
      throw Error("Recording write failed", Errc::io);
  }
  static std::string field(std::string value, std::size_t width) {
    if (value.size() > width)
      fail("EDF header field exceeds fixed width");
    for (unsigned char c : value)
      if (c < 32 || c > 126)
        fail("EDF header fields must be printable ASCII");
    value.resize(width, ' ');
    return value;
  }
  static std::string decimal(double value) {
    std::array<char, 64> buf{};
    for (int precision = 7; precision >= 0; --precision) {
      auto result = std::to_chars(buf.data(), buf.data() + buf.size(), value,
                                  std::chars_format::fixed, precision);
      if (result.ec == std::errc{} && result.ptr - buf.data() <= 8) {
        std::string out(buf.data(), result.ptr);
        if (out.find('.') != std::string::npos) {
          while (out.back() == '0')
            out.pop_back();
          if (out.back() == '.')
            out.pop_back();
        }
        if ((out == "0" || out == "-0") && value != 0)
          continue;
        return out;
      }
    }
    for (int precision = 7; precision >= 1; --precision) {
      auto result = std::to_chars(buf.data(), buf.data() + buf.size(), value,
                                  std::chars_format::general, precision);
      if (result.ec == std::errc{} && result.ptr - buf.data() <= 8)
        return {buf.data(), result.ptr};
    }
    fail("numeric header value cannot fit eight bytes");
  }
  static std::string seconds(std::int64_t ticks, bool sign) {
    const bool negative = ticks < 0;
    const auto n =
        negative ? std::uint64_t(-(ticks + 1)) + 1 : std::uint64_t(ticks);
    std::string out = negative ? "-" : (sign ? "+" : "");
    out += std::to_string(n / 10000000);
    auto frac = n % 10000000;
    if (frac) {
      auto s = std::to_string(frac + 10000000).substr(1);
      while (s.back() == '0')
        s.pop_back();
      out += '.';
      out += s;
    }
    return out;
  }
  static void validate_utf8(std::string_view text) {
    for (std::size_t i = 0; i < text.size();) {
      const auto c = static_cast<unsigned char>(text[i++]);
      if (c < 128) {
        if (!detail::allowed_annotation_byte(c))
          fail("annotation contains a forbidden control byte");
        continue;
      }
      unsigned n;
      std::uint32_t cp, min;
      if (c >= 0xc2 && c <= 0xdf) {
        n = 1;
        cp = c & 31;
        min = 128;
      } else if (c >= 0xe0 && c <= 0xef) {
        n = 2;
        cp = c & 15;
        min = 2048;
      } else if (c >= 0xf0 && c <= 0xf4) {
        n = 3;
        cp = c & 7;
        min = 65536;
      } else
        fail("invalid annotation UTF-8");
      if (i + n > text.size())
        fail("truncated annotation UTF-8");
      while (n--) {
        auto b = static_cast<unsigned char>(text[i++]);
        if ((b & 192) != 128)
          fail("invalid annotation UTF-8 continuation");
        cp = (cp << 6) | (b & 63);
      }
      if (cp < min || cp > 0x10ffff || (cp >= 0xd800 && cp <= 0xdfff))
        fail("invalid annotation Unicode scalar");
    }
  }
  std::string timekeeping(std::uint64_t record) const {
    const auto duration = info_.record_duration.count();
    if (duration && record > static_cast<std::uint64_t>(
                                 ((std::numeric_limits<std::int64_t>::max)() -
                                  info_.subsecond.count()) /
                                 duration))
      fail("record time overflow");
    auto s = seconds(static_cast<std::int64_t>(record) * duration +
                         info_.subsecond.count(),
                     true);
    s.append("\x14\x14\0", 3);
    return s;
  }
  std::string annotation_tal(const Annotation &a) const {
    if (a.onset.count() >
        (std::numeric_limits<std::int64_t>::max)() - info_.subsecond.count())
      fail("annotation onset overflow");
    auto s = seconds(a.onset.count() + info_.subsecond.count(), true);
    if (a.duration) {
      s += '\x15';
      s += seconds(a.duration->count(), false);
    }
    s += '\x14';
    s += a.text;
    s.append("\x14\0", 2);
    return s;
  }
  std::int32_t digital(std::int32_t x, const SignalSpec &s) const {
    return std::clamp(x, s.digital_min, s.digital_max);
  }
  std::int32_t physical(double x, const SignalSpec &s) const {
    if (!std::isfinite(x))
      fail("non-finite physical sample");
    if (s.physical_min < s.physical_max) {
      if (x <= s.physical_min)
        return s.digital_min;
      if (x >= s.physical_max)
        return s.digital_max;
    } else {
      if (x >= s.physical_min)
        return s.digital_min;
      if (x <= s.physical_max)
        return s.digital_max;
    }
    const auto bitvalue = (s.physical_max - s.physical_min) /
                          (double(s.digital_max) - s.digital_min);
    const auto offset = s.physical_max / bitvalue - s.digital_max;
    auto value = x / bitvalue - offset;
    if (!std::isfinite(value) || !std::isfinite(bitvalue) || bitvalue == 0) {
      const auto scale =
          (std::max)(std::abs(s.physical_min), std::abs(s.physical_max));
      const auto fraction = (x / scale - s.physical_min / scale) /
                            (s.physical_max / scale - s.physical_min / scale);
      value =
          fraction * (double(s.digital_max) - s.digital_min) + s.digital_min;
    }
    return static_cast<std::int32_t>(
        std::clamp(value, double(s.digital_min), double(s.digital_max)));
  }
  template <std::size_t Width, class T, class Convert>
  static void encode(std::span<const T> samples, const SignalSpec &signal,
                     char *dest, Convert convert) {
    for (auto sample : samples) {
      const auto bits = static_cast<std::uint32_t>(convert(sample, signal));
      dest[0] = static_cast<char>(bits);
      dest[1] = static_cast<char>(bits >> 8);
      if constexpr (Width == 3)
        dest[2] = static_cast<char>(bits >> 16);
      dest += Width;
    }
  }
  template <class T, class Convert>
  void encode_channel(std::span<const T> samples, const SignalSpec &signal,
                      char *dest, Convert convert) const {
    if (sample_bytes_ == 3)
      encode<3>(samples, signal, dest, convert);
    else
      encode<2>(samples, signal, dest, convert);
  }
  template <class T, class Convert>
  void write_record(std::span<const T> samples, Convert convert) {
    ensure_open();
    if (next_signal_ || samples.size() != data_bytes_ / sample_bytes_)
      fail("record sample count or channel order mismatch");
    encoded_.resize(data_bytes_);
    std::size_t pos = 0;
    for (const auto &signal : info_.signals) {
      encode_channel(samples.subspan(pos, signal.samples_per_record), signal,
                     encoded_.data() + pos * sample_bytes_, convert);
      pos += signal.samples_per_record;
    }
    emit_record(encoded_.data());
  }
  template <class T, class Convert>
  void write_channel(std::size_t channel, std::span<const T> samples,
                     Convert convert) {
    ensure_open();
    if (channel != next_signal_ || channel >= info_.signals.size())
      fail("signals must be written in header order");
    const auto &s = info_.signals[channel];
    if (samples.size() != s.samples_per_record)
      fail("channel sample count mismatch");
    encoded_.resize(samples.size() * sample_bytes_);
    encode_channel(samples, s, encoded_.data(), convert);
    pending_.insert(pending_.end(), encoded_.begin(), encoded_.end());
    if (++next_signal_ == info_.signals.size()) {
      emit_record(pending_.data());
      pending_.clear();
      next_signal_ = 0;
    }
  }
  void emit_record(const char *data) {
    if (records_ >= 99999999)
      fail("EDF record count exceeds eight digits");
    auto block = plus_ ? timekeeping(records_) : std::string{};
    if (plus_)
      block.resize(options_.annotation_bytes_per_signal *
                       options_.annotation_signals,
                   '\0');
    if (annotation_offset_)
      output_.write(data, static_cast<std::streamsize>(annotation_offset_));
    if (plus_)
      output_.write(block.data(), static_cast<std::streamsize>(block.size()));
    if (data_bytes_ > annotation_offset_)
      output_.write(
          data + annotation_offset_,
          static_cast<std::streamsize>(data_bytes_ - annotation_offset_));
    check_io();
    ++records_;
  }
  static std::string token(std::string value) {
    if (value.empty())
      return "X";
    std::replace(value.begin(), value.end(), ' ', '_');
    return value;
  }
  static std::string date(std::chrono::year_month_day d, bool extended) {
    const auto y = int(d.year());
    const auto m = unsigned(d.month());
    const auto day = unsigned(d.day());
    auto pad = [](unsigned n) {
      auto s = std::to_string(n);
      if (s.size() < 2)
        s = '0' + s;
      return s;
    };
    if (extended) {
      static constexpr const char *months[] = {"JAN", "FEB", "MAR", "APR",
                                               "MAY", "JUN", "JUL", "AUG",
                                               "SEP", "OCT", "NOV", "DEC"};
      auto year = std::to_string(y);
      year.insert(0, 4 - year.size(), '0');
      return pad(day) + "-" + months[m - 1] + "-" + year;
    }
    return pad(day) + '.' + pad(m) + '.' + pad(unsigned(y % 100));
  }
  static std::string birthdate(std::string value) {
    if (value.empty() || value == "X")
      return "X";
    if (value.size() != 10 || value[4] != '-' || value[7] != '-')
      fail("birthdate must be YYYY-MM-DD");
    auto number = [](std::string_view s) {
      unsigned n = 0;
      auto r = std::from_chars(s.data(), s.data() + s.size(), n);
      if (r.ec != std::errc{} || r.ptr != s.data() + s.size())
        fail("invalid birthdate");
      return n;
    };
    std::chrono::year_month_day d{
        std::chrono::year(int(number(std::string_view(value).substr(0, 4)))),
        std::chrono::month(number(std::string_view(value).substr(5, 2))),
        std::chrono::day(number(std::string_view(value).substr(8, 2)))};
    if (!d.ok())
      fail("invalid birthdate");
    return date(d, true);
  }
  Writer(const std::filesystem::path &path, RecordingInfo info,
         WriterOptions options)
      : info_(std::move(info)), options_(options) {
    sample_bytes_ =
        (info_.format == Format::bdf || info_.format == Format::bdf_plus) ? 3
                                                                          : 2;
    plus_ =
        info_.format == Format::edf_plus || info_.format == Format::bdf_plus;
    const auto max_record_bytes = sample_bytes_ * (5u * 1024u * 1024u);
    if ((info_.signals.empty() && !plus_) || info_.signals.size() > 4096)
      fail("invalid signal count");
    if (info_.record_duration.count() < 0 ||
        (info_.record_duration.count() == 0 &&
         (!plus_ || !info_.signals.empty())) ||
        seconds(info_.record_duration.count(), false).size() > 8)
      fail("invalid data record duration");
    if (info_.subsecond.count() < 0 || info_.subsecond.count() >= 10000000 ||
        (!plus_ && info_.subsecond.count() != 0))
      fail("invalid subsecond start");
    if (plus_ &&
        (!options_.annotation_signals || options_.annotation_signals > 64 ||
         options_.annotation_bytes_per_signal < 32 ||
         options_.annotation_bytes_per_signal > max_record_bytes ||
         options_.annotation_bytes_per_signal % sample_bytes_))
      fail("invalid annotation channel allocation");
    if (!plus_)
      options_.annotation_signals = 0;
    for (auto &s : info_.signals) {
      if (!s.samples_per_record || s.samples_per_record > 99999999)
        fail("invalid samples per record");
      const auto lo = sample_bytes_ == 2 ? -32768 : -8388608;
      const auto hi = sample_bytes_ == 2 ? 32767 : 8388607;
      if (s.digital_min < lo || s.digital_max > hi ||
          s.digital_min >= s.digital_max)
        fail("invalid digital limits");
      if (!std::isfinite(s.physical_min) || !std::isfinite(s.physical_max) ||
          s.physical_min == s.physical_max)
        fail("invalid physical limits");
      if (decimal(s.physical_min) == decimal(s.physical_max))
        fail("physical limits collapse in header precision");
      auto parse = [](const std::string &str) {
        double x{};
        auto r = std::from_chars(str.data(), str.data() + str.size(), x);
        if (r.ec != std::errc{})
          fail("invalid calibration");
        return x;
      };
      s.physical_min = parse(decimal(s.physical_min));
      s.physical_max = parse(decimal(s.physical_max));
      auto label = std::string_view(s.label);
      while (!label.empty() && label.back() == ' ')
        label.remove_suffix(1);
      if (label == "EDF Annotations" || label == "BDF Annotations")
        fail("reserved annotation label");
      if (s.samples_per_record >
          (max_record_bytes - data_bytes_) / sample_bytes_)
        fail("data record exceeds the EDFlib size limit");
      data_bytes_ += std::size_t(s.samples_per_record) * sample_bytes_;
    }
    annotation_index_ =
        options_.annotation_position == AnnotationPosition::start
            ? 0
            : (options_.annotation_position == AnnotationPosition::middle
                   ? info_.signals.size() / 2
                   : info_.signals.size());
    for (std::size_t i = 0; i < annotation_index_; ++i)
      annotation_offset_ +=
          std::size_t(info_.signals[i].samples_per_record) * sample_bytes_;
    if (options_.annotation_signals &&
        options_.annotation_bytes_per_signal >
            (max_record_bytes - data_bytes_) / options_.annotation_signals)
      fail("data record exceeds the EDFlib size limit");
    record_bytes_ = data_bytes_ + options_.annotation_signals *
                                      options_.annotation_bytes_per_signal;
    if (record_bytes_ > max_record_bytes)
      fail("data record exceeds the EDFlib size limit");
    const auto ns = info_.signals.size() + options_.annotation_signals;
    if (ns > 4096)
      fail("More than 4096 signals including annotation channels");
    header_bytes_ = 256 * (ns + 1);
    const auto days = std::chrono::floor<std::chrono::days>(info_.start);
    const std::chrono::year_month_day ymd{days};
    if (!ymd.ok() || int(ymd.year()) < 1985 || int(ymd.year()) > 2084)
      fail("EDF start year must be 1985 through 2084");
    const auto tod = std::chrono::hh_mm_ss(info_.start - days);
    auto pad = [](auto n) {
      auto s = std::to_string(n);
      if (s.size() < 2)
        s = '0' + s;
      return s;
    };
    std::string h;
    if (sample_bytes_ == 3)
      h = std::string("\xff"
                      "BIOSEMI",
                      8);
    else
      h = "0       ";
    std::string patient = info_.patient, recording = info_.recording;
    if (plus_) {
      if (!info_.sex.empty() && info_.sex != "M" && info_.sex != "F" &&
          info_.sex != "X")
        fail("sex must be M, F or X");
      // EDFlib requires the separator after the patient code at byte 72 or
      // earlier.
      if (info_.patient_code.size() > 72)
        fail("patient code exceeds the EDFlib header limit");
      patient = token(info_.patient_code) + ' ' + token(info_.sex) + ' ' +
                birthdate(info_.birthdate) + ' ' + token(info_.patient_name);
      if (!info_.patient_additional.empty())
        patient += ' ' + info_.patient_additional;
      recording = "Startdate " + date(ymd, true) + ' ' +
                  token(info_.admin_code) + ' ' + token(info_.technician) +
                  ' ' + token(info_.equipment);
      if (!info_.recording_additional.empty())
        recording += ' ' + info_.recording_additional;
    }
    h += field(patient, 80) + field(recording, 80) + field(date(ymd, false), 8);
    h += field(pad(tod.hours().count()) + '.' + pad(tod.minutes().count()) +
                   '.' + pad(tod.seconds().count()),
               8);
    h += field(std::to_string(header_bytes_), 8) +
         field(plus_ ? (sample_bytes_ == 3 ? "BDF+C" : "EDF+C") : "", 44);
    h += field("-1", 8) +
         field(seconds(info_.record_duration.count(), false), 8) +
         field(std::to_string(ns), 4);
    std::vector<SignalSpec> specs = info_.signals;
    for (std::size_t i = 0; i < options_.annotation_signals; ++i) {
      SignalSpec s{};
      s.label = sample_bytes_ == 3 ? "BDF Annotations" : "EDF Annotations";
      s.physical_min = -1;
      s.physical_max = 1;
      s.digital_min = sample_bytes_ == 3 ? -8388608 : -32768;
      s.digital_max = sample_bytes_ == 3 ? 8388607 : 32767;
      s.samples_per_record = static_cast<std::uint32_t>(
          options_.annotation_bytes_per_signal / sample_bytes_);
      specs.insert(specs.begin() +
                       static_cast<std::ptrdiff_t>(annotation_index_ + i),
                   std::move(s));
    }
    auto fields = [&](std::size_t width, auto get) {
      for (const auto &s : specs)
        h += field(get(s), width);
    };
    fields(16, [](const auto &s) { return s.label; });
    fields(80, [](const auto &s) { return s.transducer; });
    fields(8, [](const auto &s) { return s.unit; });
    fields(8, [](const auto &s) { return decimal(s.physical_min); });
    fields(8, [](const auto &s) { return decimal(s.physical_max); });
    fields(8, [](const auto &s) { return std::to_string(s.digital_min); });
    fields(8, [](const auto &s) { return std::to_string(s.digital_max); });
    fields(80, [](const auto &s) { return s.prefilter; });
    fields(8,
           [](const auto &s) { return std::to_string(s.samples_per_record); });
    fields(32, [](const auto &) { return std::string{}; });
    output_.open(path, std::ios::binary | std::ios::trunc);
    if (!output_)
      throw Error("Cannot create recording file", Errc::io);
    output_.write(h.data(), static_cast<std::streamsize>(h.size()));
    check_io();
  }
};
} // namespace fastbdf
