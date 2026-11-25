#pragma once
#include "detail/file_source.hpp"
#include "types.hpp"
#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <filesystem>
#include <initializer_list>
#include <ios>
#include <limits>
#include <memory>
#include <mutex>
#include <span>
#include <string_view>
#include <type_traits>
#include <vector>

namespace fastbdf {
enum class AnnotationMode { none, lazy, all, until_recording_ends };

// Position-independent reads and independent Cursor objects may run
// concurrently. Reopening the reader and its legacy implicit cursors require
// synchronization.
class Reader {
  struct Channel {
    SignalSpec spec;
    std::uint64_t offset{};
    bool annotation{};
  };
  RecordingInfo info_{};
  std::optional<detail::FileSource> file_;
  std::vector<Channel> channels_;
  std::vector<std::size_t> data_channels_;
  std::vector<std::uint64_t> cursors_;
  mutable std::vector<Annotation> annotations_;
  std::vector<Tick> record_onsets_;
  std::uint64_t records_{}, header_bytes_{}, record_bytes_{};
  unsigned width_{};
  bool discontinuous_{};
  mutable Tick current_record_onset_{};
  mutable bool annotations_loaded_ = true;
  std::unique_ptr<std::mutex> annotation_mutex_ =
      std::make_unique<std::mutex>();
  std::string raw_patient_, raw_recording_, reserved_;
  static constexpr std::int64_t ticks_per_second = 10000000;

  static std::string trim(std::string_view s) {
    while (!s.empty() && s.back() == ' ')
      s.remove_suffix(1);
    while (!s.empty() && s.front() == ' ')
      s.remove_prefix(1);
    return std::string(s);
  }
  static std::int64_t integer(std::string_view s, const char *field) {
    auto t = trim(s);
    if (t.size() >= 2 && t[0] == '+' && t[1] >= '0' && t[1] <= '9')
      t.erase(t.begin());
    std::int64_t n{};
    auto [p, e] = std::from_chars(t.data(), t.data() + t.size(), n);
    if (e != std::errc{} || p != t.data() + t.size() || t.empty())
      throw Error(std::string("Invalid ") + field);
    return n;
  }
  static double number(std::string_view s, const char *field) {
    auto t = trim(s);
    double n{};
    if (!t.empty() && t.front() == '+')
      t.erase(t.begin());
    auto [p, e] = std::from_chars(t.data(), t.data() + t.size(), n);
    if (e != std::errc{} || p != t.data() + t.size() || t.empty() ||
        !std::isfinite(n))
      throw Error(std::string("Invalid ") + field);
    return n;
  }
  static Tick ticks(std::string_view s, const char *field) {
    auto t = trim(s);
    std::size_t p = 0;
    bool negative = false;
    if (p < t.size() && (t[p] == '+' || t[p] == '-'))
      negative = t[p++] == '-';
    std::uint64_t whole = 0, fraction = 0;
    unsigned digits = 0;
    bool any = false;
    const auto limit =
        static_cast<std::uint64_t>((std::numeric_limits<std::int64_t>::max)()) +
        (negative ? 1u : 0u);
    while (p < t.size() && t[p] >= '0' && t[p] <= '9') {
      any = true;
      auto d = static_cast<unsigned>(t[p++] - '0');
      if (whole > (limit / ticks_per_second - d) / 10)
        throw Error("Time value overflows");
      whole = whole * 10 + d;
    }
    if (p < t.size() && t[p] == '.') {
      ++p;
      while (p < t.size() && t[p] >= '0' && t[p] <= '9') {
        any = true;
        auto d = static_cast<unsigned>(t[p++] - '0');
        if (digits < 7) {
          fraction = fraction * 10 + d;
          ++digits;
        } else if (d)
          throw Error("Time precision exceeds 100 ns");
      }
    }
    if (!any || p != t.size())
      throw Error(std::string("Invalid ") + field);
    while (digits++ < 7)
      fraction *= 10;
    if (whole > limit / ticks_per_second ||
        whole * ticks_per_second > limit - fraction)
      throw Error("Time value overflows");
    const auto magnitude = whole * ticks_per_second + fraction;
    if (negative && magnitude == limit)
      return (Tick::min)();
    auto value = static_cast<std::int64_t>(magnitude);
    return Tick{negative ? -value : value};
  }
  static std::uint64_t product(std::uint64_t a, std::uint64_t b) {
    if (b && a > (std::numeric_limits<std::uint64_t>::max)() / b)
      throw Error("File size overflows");
    return a * b;
  }
  void normalize_annotations() const {
    for (auto &annotation : annotations_) {
      if (annotation.onset.count() <
          (std::numeric_limits<std::int64_t>::min)() + info_.subsecond.count())
        throw Error("Annotation onset overflows");
      annotation.onset -= info_.subsecond;
    }
  }
  // Exact ceil(a*b/c), with a<c and b<=UINT32_MAX, without nonstandard
  // 128-bit integers or floating-point loss at sample boundaries.
  static std::uint64_t sample_fraction(std::uint64_t a, std::uint32_t b,
                                       std::uint64_t c) {
    if (!a)
      return 0;
    if (a <= (std::numeric_limits<std::uint64_t>::max)() / b) {
      const auto p = a * b;
      return p / c + (p % c != 0);
    }
    std::uint64_t quotient = 0, remainder = 0;
    for (int bit = 31; bit >= 0; --bit) {
      quotient *= 2;
      remainder *= 2;
      if (remainder >= c) {
        remainder -= c;
        ++quotient;
      }
      if ((b >> bit) & 1) {
        remainder += a;
        if (remainder >= c) {
          remainder -= c;
          ++quotient;
        }
      }
    }
    return quotient + (remainder != 0);
  }
  void read_bytes(std::uint64_t offset, std::span<char> out) const {
    if (!file_)
      throw Error("Reader is not open");
    file_->read_at(offset, out);
  }
  static bool utf8(std::string_view s) {
    std::size_t i = 0;
    while (i < s.size()) {
      auto c = static_cast<unsigned char>(s[i++]);
      if (c < 0x80)
        continue;
      unsigned n;
      std::uint32_t cp;
      if (c >= 0xc2 && c <= 0xdf) {
        n = 1;
        cp = c & 31;
      } else if (c >= 0xe0 && c <= 0xef) {
        n = 2;
        cp = c & 15;
      } else if (c >= 0xf0 && c <= 0xf4) {
        n = 3;
        cp = c & 7;
      } else
        return false;
      if (i + n > s.size())
        return false;
      for (unsigned j = 0; j < n; ++j) {
        auto d = static_cast<unsigned char>(s[i++]);
        if ((d & 0xc0) != 0x80)
          return false;
        cp = (cp << 6) | (d & 63);
      }
      if ((n == 1 && cp < 0x80) || (n == 2 && cp < 0x800) ||
          (n == 3 && cp < 0x10000) || cp > 0x10ffff ||
          (cp >= 0xd800 && cp <= 0xdfff))
        return false;
    }
    return true;
  }
  static std::vector<std::string> words(std::string_view s) {
    std::vector<std::string> result;
    std::size_t p = 0;
    while (p < s.size()) {
      while (p < s.size() && s[p] == ' ')
        ++p;
      if (p == s.size())
        break;
      auto e = s.find(' ', p);
      if (e == s.npos)
        e = s.size();
      result.emplace_back(s.substr(p, e - p));
      p = e;
    }
    return result;
  }
  static std::string decoded(std::string s) {
    if (s == "X")
      return {};
    std::replace(s.begin(), s.end(), '_', ' ');
    return s;
  }
  void parse_identification() {
    info_.patient = raw_patient_;
    info_.recording = raw_recording_;
    if (info_.format != Format::edf_plus && info_.format != Format::bdf_plus)
      return;
    auto p = words(raw_patient_);
    auto r = words(raw_recording_);
    if (p.size() < 4 || r.size() < 5 || r[0] != "Startdate")
      throw Error("Invalid EDF+ identification fields");
    info_.patient_code = decoded(p[0]);
    info_.sex = decoded(p[1]);
    if (!info_.sex.empty() && info_.sex != "M" && info_.sex != "F")
      throw Error("Invalid patient sex");
    info_.birthdate = decoded(p[2]);
    info_.patient_name = decoded(p[3]);
    if (!info_.birthdate.empty()) {
      static constexpr std::array<std::string_view, 12> months = {
          "JAN", "FEB", "MAR", "APR", "MAY", "JUN",
          "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"};
      const auto &b = info_.birthdate;
      if (b.size() != 11 || b[2] != '-' || b[6] != '-')
        throw Error("Invalid patient birthdate");
      auto m = std::find(months.begin(), months.end(),
                         std::string_view(b).substr(3, 3));
      if (m == months.end())
        throw Error("Invalid patient birthdate");
      auto day = integer(std::string_view(b).substr(0, 2), "birth day"),
           year = integer(std::string_view(b).substr(7), "birth year");
      auto date =
          std::chrono::year{static_cast<int>(year)} /
          std::chrono::month{static_cast<unsigned>(m - months.begin() + 1)} /
          std::chrono::day{static_cast<unsigned>(day)};
      if (!date.ok())
        throw Error("Invalid patient birthdate");
      const auto month = static_cast<unsigned>(m - months.begin() + 1);
      info_.birthdate = b.substr(7) + "-" + (month < 10 ? "0" : "") +
                        std::to_string(month) + "-" + b.substr(0, 2);
    }
    for (std::size_t i = 4; i < p.size(); ++i) {
      if (i > 4)
        info_.patient_additional += ' ';
      info_.patient_additional += decoded(p[i]);
    }
    info_.admin_code = decoded(r[2]);
    info_.technician = decoded(r[3]);
    info_.equipment = decoded(r[4]);
    for (std::size_t i = 5; i < r.size(); ++i) {
      if (i > 5)
        info_.recording_additional += ' ';
      info_.recording_additional += decoded(r[i]);
    }
    // EDF+ uses the unambiguous four-digit date in the recording field.
    if (r[1] != "X") {
      static constexpr std::array<std::string_view, 12> months = {
          "JAN", "FEB", "MAR", "APR", "MAY", "JUN",
          "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"};
      const auto &d = r[1];
      if (d.size() != 11 || d[2] != '-' || d[6] != '-')
        throw Error("Invalid EDF+ start date");
      auto it = std::find(months.begin(), months.end(),
                          std::string_view(d).substr(3, 3));
      if (it == months.end())
        throw Error("Invalid EDF+ month");
      auto date =
          std::chrono::year{static_cast<int>(
              integer(std::string_view(d).substr(7), "year"))} /
          std::chrono::month{static_cast<unsigned>(it - months.begin() + 1)} /
          std::chrono::day{static_cast<unsigned>(
              integer(std::string_view(d).substr(0, 2), "day"))};
      if (!date.ok())
        throw Error("Invalid EDF+ date");
      auto original = std::chrono::year_month_day{
          std::chrono::floor<std::chrono::days>(info_.start)};
      if (original.month() != date.month() || original.day() != date.day() ||
          int(original.year()) % 100 != int(date.year()) % 100)
        throw Error("EDF+ start dates disagree");
      auto tod =
          info_.start - std::chrono::floor<std::chrono::days>(info_.start);
      info_.start = std::chrono::local_days{date} + tod;
    }
  }
  // TAL bytes are stored consecutively, including all bytes of 16/24 bit words.
  bool parse_tals(std::span<const char> bytes, std::uint64_t,
                  bool first_channel, AnnotationMode mode) const {
    std::size_t p = 0;
    if (!first_channel && mode == AnnotationMode::none)
      return false;
    bool first = true;
    while (p < bytes.size()) {
      if (bytes[p] == '\0') {
        ++p;
        continue;
      }
      const auto begin = p;
      while (p < bytes.size() && bytes[p] != '\0')
        ++p;
      if (p == bytes.size())
        throw Error("Unterminated annotation TAL");
      std::string_view tal(bytes.data() + begin, p - begin);
      ++p;
      auto sep = tal.find_first_of("\x14\x15");
      if (sep == tal.npos || tal.empty() || (tal[0] != '+' && tal[0] != '-'))
        throw Error("Invalid annotation onset");
      auto onset = ticks(tal.substr(0, sep), "annotation onset");
      std::optional<Tick> duration;
      auto q = sep;
      if (tal[q] == '\x15') {
        auto e = tal.find('\x14', q + 1);
        if (e == tal.npos)
          throw Error("Invalid annotation duration");
        duration = ticks(tal.substr(q + 1, e - q - 1), "annotation duration");
        if (duration->count() < 0)
          throw Error("Negative annotation duration");
        q = e;
      }
      if (first && first_channel) {
        if (q + 1 >= tal.size() || tal[q + 1] != '\x14' || duration)
          throw Error("Missing record timekeeping TAL");
        current_record_onset_ = onset;
        if (mode == AnnotationMode::none)
          return false;
      }
      first = false;
      while (q < tal.size()) {
        if (tal[q] != '\x14')
          throw Error("Invalid annotation separator");
        auto e = tal.find('\x14', q + 1);
        if (e == tal.npos)
          throw Error("Unterminated annotation text");
        auto text = tal.substr(q + 1, e - q - 1);
        if (!text.empty() && mode != AnnotationMode::none) {
          if (!utf8(text))
            throw Error("Invalid annotation UTF-8");
          annotations_.push_back({onset, duration, std::string(text)});
          if (mode == AnnotationMode::until_recording_ends &&
              text.starts_with("Recording ends"))
            return true;
        }
        q = e;
        if (q + 1 == tal.size())
          break;
      }
    }
    if (first && first_channel)
      throw Error("Missing record timekeeping TAL");
    return false;
  }
  const Channel &channel(std::size_t signal) const {
    if (signal >= data_channels_.size())
      throw Error("Signal index out of range");
    return channels_[data_channels_[signal]];
  }
  template <class T>
  void decode(const unsigned char *ptr, std::span<T> out,
              const SignalSpec &s) const {
    const double gain = (s.physical_max - s.physical_min) /
                        (static_cast<double>(s.digital_max) - s.digital_min);
    for (auto &sample : out) {
      std::uint32_t v = std::uint32_t(ptr[0]) | (std::uint32_t(ptr[1]) << 8);
      if (width_ == 3)
        v |= std::uint32_t(ptr[2]) << 16;
      const auto raw = width_ == 2
                           ? static_cast<std::int32_t>(v ^ 0x8000) - 0x8000
                           : static_cast<std::int32_t>(v ^ 0x800000) - 0x800000;
      const auto value = std::clamp(raw, s.digital_min, s.digital_max);
      if constexpr (std::is_same_v<T, std::int32_t>)
        sample = value;
      else if (std::isfinite(gain))
        sample = s.physical_min +
                 (static_cast<double>(value) - s.digital_min) * gain;
      else
        sample =
            std::lerp(s.physical_min, s.physical_max,
                      (static_cast<double>(value) - s.digital_min) /
                          (static_cast<double>(s.digital_max) - s.digital_min));
      ptr += width_;
    }
  }
  template <class T>
  std::size_t read_at(std::size_t signal, std::uint64_t first,
                      std::span<T> out) const {
    const auto &ch = channel(signal);
    const auto total = sample_count(signal);
    if (first > total)
      throw Error("Sample offset out of range");
    auto count = static_cast<std::size_t>(
        std::min<std::uint64_t>(out.size(), total - first));
    if (!count)
      return 0;
    const auto spr = ch.spec.samples_per_record;
    const auto first_record = first / spr,
               last_record = (first + count - 1) / spr;
    if (first_record == last_record && count <= 1024 / width_) {
      std::array<char, 1024> bytes;
      read_bytes(header_bytes_ + first_record * record_bytes_ + ch.offset +
                     (first % spr) * width_,
                 std::span<char>(bytes).first(count * width_));
      decode(reinterpret_cast<const unsigned char *>(bytes.data()),
             out.first(count), ch.spec);
      return count;
    }
    if (first_record == last_record ||
        product(spr, width_) * 4 < record_bytes_) {
      std::vector<char> bytes;
      std::size_t written = 0;
      for (auto r = first_record; r <= last_record; ++r) {
        const auto lo = r == first_record ? first % spr : 0;
        const auto hi = r == last_record ? (first + count - 1) % spr + 1 : spr;
        bytes.resize(static_cast<std::size_t>(hi - lo) * width_);
        read_bytes(header_bytes_ + r * record_bytes_ + ch.offset + lo * width_,
                   bytes);
        decode(reinterpret_cast<const unsigned char *>(bytes.data()),
               out.subspan(written, static_cast<std::size_t>(hi - lo)),
               ch.spec);
        written += static_cast<std::size_t>(hi - lo);
      }
      return written;
    }
    // Read contiguous slabs rather than issuing one seek per sample or record.
    // Bound temporary storage so enormous recordings remain streamable.
    constexpr std::uint64_t target = 8 * 1024 * 1024;
    const auto batch_records =
        std::max<std::uint64_t>(1, target / record_bytes_);
    std::vector<char> buffer;
    std::size_t written = 0;
    for (auto r = first_record; r <= last_record;) {
      const auto n = (std::min)(batch_records, last_record - r + 1);
      const auto size = product(n, record_bytes_);
      if (size > buffer.max_size())
        throw Error("Data record too large");
      buffer.resize(static_cast<std::size_t>(size));
      read_bytes(header_bytes_ + product(r, record_bytes_), buffer);
      for (std::uint64_t j = 0; j < n; ++j) {
        const auto record = r + j;
        const auto lo = record == first_record ? first % spr : 0;
        const auto hi =
            record == last_record ? (first + count - 1) % spr + 1 : spr;
        auto *ptr = reinterpret_cast<const unsigned char *>(
            buffer.data() + j * record_bytes_ + ch.offset + lo * width_);
        decode(ptr, out.subspan(written, static_cast<std::size_t>(hi - lo)),
               ch.spec);
        written += static_cast<std::size_t>(hi - lo);
      }
      r += n;
    }
    return written;
  }
  template <class T>
  void window_into(std::span<const std::size_t> signals, TimeRange range,
                   std::vector<std::vector<T>> &result) const {
    if (range.begin.count() < 0 || range.end < range.begin)
      throw Error("Invalid time range");
    std::vector<std::pair<std::uint64_t, std::uint64_t>> bounds;
    bounds.reserve(signals.size());
    std::uint64_t selected_bytes = 0, first_record = records_, end_record = 0;
    for (auto id : signals) {
      const auto spr = channel(id).spec.samples_per_record;
      const auto lo = sample_index_at_or_after(id, range.begin);
      const auto hi = sample_index_at_or_after(id, range.end);
      if (hi - lo > std::vector<T>{}.max_size())
        throw Error("Requested window too large");
      bounds.emplace_back(lo, hi);
      selected_bytes += product(spr, width_);
      if (lo < hi) {
        first_record = (std::min)(first_record, lo / spr);
        end_record = (std::max)(end_record, (hi - 1) / spr + 1);
      }
    }
    result.resize(signals.size());
    for (std::size_t i = 0; i < signals.size(); ++i)
      result[i].resize(
          static_cast<std::size_t>(bounds[i].second - bounds[i].first));
    if (first_record >= end_record)
      return;
    if (signals.size() == 1 || selected_bytes < record_bytes_ / 4) {
      for (std::size_t i = 0; i < signals.size(); ++i)
        read_at(signals[i], bounds[i].first, std::span<T>(result[i]));
      return;
    }
    const auto batch =
        std::max<std::uint64_t>(1, 8 * 1024 * 1024 / record_bytes_);
    std::vector<char> buffer;
    for (auto r = first_record; r < end_record;) {
      const auto n = (std::min)(batch, end_record - r);
      buffer.resize(static_cast<std::size_t>(product(n, record_bytes_)));
      read_bytes(header_bytes_ + r * record_bytes_, buffer);
      for (std::size_t i = 0; i < signals.size(); ++i) {
        const auto &ch = channel(signals[i]);
        const auto spr = ch.spec.samples_per_record;
        const auto [lo, hi] = bounds[i];
        for (auto rec = r; rec < r + n; ++rec) {
          const auto begin = (std::max)(lo, rec * spr),
                     end = (std::min)(hi, (rec + 1) * spr);
          if (end <= begin)
            continue;
          auto *ptr = reinterpret_cast<const unsigned char *>(
              buffer.data() + (rec - r) * record_bytes_ + ch.offset +
              (begin - rec * spr) * width_);
          decode(ptr,
                 std::span<T>(result[i]).subspan(
                     static_cast<std::size_t>(begin - lo),
                     static_cast<std::size_t>(end - begin)),
                 ch.spec);
        }
      }
      r += n;
    }
  }
  template <class T>
  std::vector<std::vector<T>> window(std::span<const std::size_t> signals,
                                     TimeRange range) const {
    std::vector<std::vector<T>> result;
    window_into(signals, range, result);
    return result;
  }

public:
  class Cursor {
    const Reader *reader_;
    std::size_t id_;
    std::uint64_t position_{};

  public:
    Cursor(const Reader &reader, std::size_t id) : reader_(&reader), id_(id) {
      reader.channel(id);
    }
    std::uint64_t tell() const noexcept { return position_; }
    void seek(std::uint64_t sample) {
      if (sample > reader_->sample_count(id_))
        throw Error("Sample offset out of range");
      position_ = sample;
    }
    void rewind() noexcept { position_ = 0; }
    void seek(std::int64_t offset, std::ios_base::seekdir direction) {
      std::uint64_t base;
      if (direction == std::ios_base::beg)
        base = 0;
      else if (direction == std::ios_base::cur)
        base = position_;
      else if (direction == std::ios_base::end)
        base = reader_->sample_count(id_);
      else
        throw Error("Invalid seek direction");
      if (offset < 0) {
        auto magnitude = static_cast<std::uint64_t>(-(offset + 1)) + 1;
        if (magnitude > base)
          throw Error("Seek before first sample");
        seek(base - magnitude);
      } else {
        auto positive = static_cast<std::uint64_t>(offset);
        if (positive > (std::numeric_limits<std::uint64_t>::max)() - base)
          throw Error("Seek overflows");
        seek(base + positive);
      }
    }
    std::size_t read_digital(std::span<std::int32_t> out) {
      auto n = reader_->read_digital_at(id_, position_, out);
      position_ += n;
      return n;
    }
    std::size_t read_physical(std::span<double> out) {
      auto n = reader_->read_physical_at(id_, position_, out);
      position_ += n;
      return n;
    }
  };
  Reader() = default;
  explicit Reader(const std::filesystem::path &path,
                  AnnotationMode mode = AnnotationMode::lazy) {
    reopen(path, mode);
  }
  Reader(const Reader &) = delete;
  Reader &operator=(const Reader &) = delete;
  Reader(Reader &&) = default;
  Reader &operator=(Reader &&) = default;
  static Reader open(const std::filesystem::path &path,
                     AnnotationMode mode = AnnotationMode::lazy) {
    return Reader(path, mode);
  }
  void reopen(const std::filesystem::path &path,
              AnnotationMode mode = AnnotationMode::lazy) {
    // Parse into a temporary object so failed reopen preserves the previous
    // reader.
    Reader next;
    try {
      next.open_impl(path, mode);
    } catch (const Error &error) {
      if (error.code() == Errc::io)
        throw;
      throw Error(error.what(), Errc::invalid_file);
    }
    *this = std::move(next);
  }

private:
  void open_impl(const std::filesystem::path &path, AnnotationMode mode) {
    file_.emplace(path);
    std::array<char, 256> fixed{};
    read_bytes(0, fixed);
    std::string_view f(fixed.data(), fixed.size());
    const bool bdf =
        static_cast<unsigned char>(f[0]) == 0xff && f.substr(1, 7) == "BIOSEMI";
    if (!bdf && f.substr(0, 8) != "0       ")
      throw Error("Unsupported EDF/BDF version");
    width_ = bdf ? 3 : 2;
    for (std::size_t i = 8; i < f.size(); ++i)
      if (static_cast<unsigned char>(f[i]) < 32 ||
          static_cast<unsigned char>(f[i]) > 126)
        throw Error("Non-ASCII EDF/BDF header");
    raw_patient_ = trim(f.substr(8, 80));
    raw_recording_ = trim(f.substr(88, 80));
    reserved_ = trim(f.substr(192, 44));
    const auto plus = reserved_.starts_with(bdf ? "BDF+" : "EDF+");
    if (plus && !(reserved_.starts_with(bdf ? "BDF+C" : "EDF+C") ||
                  reserved_.starts_with(bdf ? "BDF+D" : "EDF+D")))
      throw Error("Invalid EDF+ continuity marker");
    discontinuous_ = plus && reserved_[4] == 'D';
    info_.format = bdf ? (plus ? Format::bdf_plus : Format::bdf)
                       : (plus ? Format::edf_plus : Format::edf);
    auto date = f.substr(168, 8), time = f.substr(176, 8);
    if (date[2] != '.' || date[5] != '.' || time[2] != '.' || time[5] != '.')
      throw Error("Invalid start date/time");
    auto y = integer(date.substr(6, 2), "year");
    y += y >= 85 ? 1900 : 2000;
    auto d = std::chrono::year{static_cast<int>(y)} /
             std::chrono::month{
                 static_cast<unsigned>(integer(date.substr(3, 2), "month"))} /
             std::chrono::day{
                 static_cast<unsigned>(integer(date.substr(0, 2), "day"))};
    auto hh = integer(time.substr(0, 2), "hour"),
         mm = integer(time.substr(3, 2), "minute"),
         ss = integer(time.substr(6, 2), "second");
    if (!d.ok() || hh < 0 || hh > 23 || mm < 0 || mm > 59 || ss < 0 || ss > 59)
      throw Error("Invalid start date/time");
    info_.start = std::chrono::local_days{d} + std::chrono::hours{hh} +
                  std::chrono::minutes{mm} + std::chrono::seconds{ss};
    auto ns = integer(f.substr(252, 4), "signal count");
    if (ns <= 0 || ns > 4096)
      throw Error("Invalid signal count");
    auto hb = integer(f.substr(184, 8), "header length");
    if (hb != 256 + 256 * ns)
      throw Error("Invalid header length");
    header_bytes_ = static_cast<std::uint64_t>(hb);
    const auto declared = integer(f.substr(236, 8), "record count");
    if (declared < -1)
      throw Error("Invalid record count");
    info_.record_duration = ticks(f.substr(244, 8), "record duration");
    if (info_.record_duration.count() < 0)
      throw Error("Record duration must be nonnegative");
    std::vector<char> raw(static_cast<std::size_t>(ns) * 256);
    read_bytes(256, raw);
    for (char c : raw)
      if (static_cast<unsigned char>(c) < 32 ||
          static_cast<unsigned char>(c) > 126)
        throw Error("Non-ASCII signal header");
    auto field = [&](std::size_t offset, std::size_t length, std::size_t i) {
      return std::string_view(raw.data() + offset * ns + length * i, length);
    };
    for (std::size_t i = 0; i < static_cast<std::size_t>(ns); ++i) {
      Channel c;
      auto &s = c.spec;
      s.label = trim(field(0, 16, i));
      s.transducer = trim(field(16, 80, i));
      s.unit = trim(field(96, 8, i));
      s.physical_min = number(field(104, 8, i), "physical minimum");
      s.physical_max = number(field(112, 8, i), "physical maximum");
      auto mn = integer(field(120, 8, i), "digital minimum"),
           mx = integer(field(128, 8, i), "digital maximum");
      auto limit = bdf ? 8388608 : 32768;
      if (mn < -limit || mx >= limit || mn >= mx)
        throw Error("Invalid digital range");
      s.digital_min = static_cast<std::int32_t>(mn);
      s.digital_max = static_cast<std::int32_t>(mx);
      if (s.physical_min == s.physical_max)
        throw Error("Invalid physical range");
      s.prefilter = trim(field(136, 80, i));
      auto samples = integer(field(216, 8, i), "samples per record");
      if (samples <= 0 || samples > (std::numeric_limits<std::uint32_t>::max)())
        throw Error("Invalid samples per record");
      s.samples_per_record = static_cast<std::uint32_t>(samples);
      c.offset = record_bytes_;
      record_bytes_ += product(static_cast<std::uint64_t>(samples), width_);
      c.annotation =
          plus && s.label == (bdf ? "BDF Annotations" : "EDF Annotations");
      if (!c.annotation) {
        data_channels_.push_back(channels_.size());
        info_.signals.push_back(s);
      }
      channels_.push_back(std::move(c));
    }
    if (info_.record_duration.count() == 0 &&
        (!plus || !data_channels_.empty()))
      throw Error("Zero duration requires annotation-only EDF+");
    if (record_bytes_ > width_ * (5u * 1024u * 1024u))
      throw Error("Data record exceeds the EDFlib size limit");
    const auto size = file_->size();
    if (size < header_bytes_ || (size - header_bytes_) % record_bytes_)
      throw Error("Incomplete data record");
    records_ = (size - header_bytes_) / record_bytes_;
    if (declared >= 0 && static_cast<std::uint64_t>(declared) != records_)
      throw Error("Record count does not match file size");
    if (info_.record_duration.count() &&
        records_ > static_cast<std::uint64_t>(
                       (std::numeric_limits<std::int64_t>::max)() /
                       info_.record_duration.count()))
      throw Error("Recording duration overflows");
    if (discontinuous_) {
      if (records_ > record_onsets_.max_size())
        throw Error("Too many data records");
      record_onsets_.resize(static_cast<std::size_t>(records_));
    }
    parse_identification();
    cursors_.resize(data_channels_.size());
    if (plus) {
      if (std::none_of(channels_.begin(), channels_.end(),
                       [](const Channel &c) { return c.annotation; }))
        throw Error("EDF+ requires an annotation signal");
      // EDFlib only honors the stop marker for a single annotation channel.
      if (mode == AnnotationMode::until_recording_ends &&
          channels_.size() - data_channels_.size() > 1)
        mode = AnnotationMode::all;
      annotations_loaded_ = mode != AnnotationMode::lazy;
      bool stopped = false;
      std::vector<char> bytes;
      Tick previous{};
      const auto scan_records =
          (!discontinuous_ &&
           (mode == AnnotationMode::none || mode == AnnotationMode::lazy))
              ? std::min<std::uint64_t>(1, records_)
              : records_;
      for (std::uint64_t r = 0; r < scan_records; ++r) {
        bool first = true;
        for (const auto &c : channels_)
          if (c.annotation) {
            bytes.resize(static_cast<std::size_t>(c.spec.samples_per_record) *
                         width_);
            read_bytes(header_bytes_ + r * record_bytes_ + c.offset, bytes);
            const bool stop = parse_tals(
                bytes, r, first,
                (stopped || mode == AnnotationMode::lazy) ? AnnotationMode::none
                                                          : mode);
            first = false;
            stopped |= stop;
            if (stopped || mode == AnnotationMode::none ||
                mode == AnnotationMode::lazy)
              break;
          }
        if (r == 0) {
          info_.subsecond = current_record_onset_;
          if (info_.subsecond.count() < 0 ||
              info_.subsecond.count() >= ticks_per_second)
            throw Error("Invalid first record onset");
        }
        if (r && previous.count() > (std::numeric_limits<std::int64_t>::max)() -
                                        info_.record_duration.count())
          throw Error("Record onset overflows");
        if (r &&
            (discontinuous_
                 ? current_record_onset_ < previous + info_.record_duration
                 : current_record_onset_ != previous + info_.record_duration))
          throw Error("Invalid record onset sequence");
        if (discontinuous_)
          record_onsets_[r] = current_record_onset_;
        previous = current_record_onset_;
        if (stopped && !discontinuous_)
          break;
      }
      normalize_annotations();
    }
  }

public:
  const RecordingInfo &header() const noexcept { return info_; }
  const std::vector<SignalSpec> &signals() const noexcept {
    return info_.signals;
  }
  const std::vector<Annotation> &annotations() const {
    std::lock_guard lock(*annotation_mutex_);
    if (!annotations_loaded_) {
      std::uint64_t record = 0;
      try {
        std::vector<char> bytes;
        Tick previous{};
        for (std::uint64_t r = 0; r < records_; ++r) {
          record = r;
          bool first = true;
          for (const auto &c : channels_)
            if (c.annotation) {
              bytes.resize(static_cast<std::size_t>(c.spec.samples_per_record) *
                           width_);
              read_bytes(header_bytes_ + r * record_bytes_ + c.offset, bytes);
              parse_tals(bytes, r, first, AnnotationMode::all);
              first = false;
            }
          if (r &&
              previous.count() > (std::numeric_limits<std::int64_t>::max)() -
                                     info_.record_duration.count())
            throw Error("Record onset overflows");
          if (r &&
              (discontinuous_
                   ? current_record_onset_ < previous + info_.record_duration
                   : current_record_onset_ != previous + info_.record_duration))
            throw Error("Invalid record onset sequence");
          previous = current_record_onset_;
        }
        normalize_annotations();
        annotations_loaded_ = true;
      } catch (const Error &error) {
        annotations_.clear();
        if (error.code() == Errc::io)
          throw;
        throw Error("Annotation record " + std::to_string(record) + ": " +
                        error.what(),
                    Errc::invalid_file);
      } catch (...) {
        annotations_.clear();
        throw;
      }
    }
    return annotations_;
  }
  const std::vector<Tick> &record_onsets() const noexcept {
    return record_onsets_;
  }
  std::uint64_t record_count() const noexcept { return records_; }
  bool discontinuous() const noexcept { return discontinuous_; }
  std::uint64_t sample_count(std::size_t id) const {
    return product(records_, channel(id).spec.samples_per_record);
  }
  // Elapsed from the first record, including any discontinuities.
  Tick record_start(std::uint64_t record) const {
    if (record >= records_)
      throw Error("Record index out of range");
    if (discontinuous_)
      return record_onsets_[static_cast<std::size_t>(record)] - info_.subsecond;
    return Tick{static_cast<std::int64_t>(record) *
                info_.record_duration.count()};
  }

  // A time inside a gap maps to the next sample. Beyond EOF returns
  // sample_count(id).
  std::uint64_t sample_index_at_or_after(std::size_t id, Tick time) const {
    const auto spr = channel(id).spec.samples_per_record;
    if (time < Tick::zero())
      throw Error("Negative sample lookup time");
    if (!records_)
      return 0;
    const auto duration = info_.record_duration.count();
    std::uint64_t record;
    if (discontinuous_) {
      auto next = std::upper_bound(
          record_onsets_.begin(), record_onsets_.end(), time,
          [&](Tick t, Tick onset) { return t < onset - info_.subsecond; });
      record = static_cast<std::uint64_t>(next - record_onsets_.begin() - 1);
    } else {
      record = static_cast<std::uint64_t>(time.count() / duration);
    }
    if (record >= records_)
      return sample_count(id);
    const auto delta = (time - record_start(record)).count();
    if (delta >= duration)
      return (record + 1) * spr;
    return record * spr + sample_fraction(static_cast<std::uint64_t>(delta),
                                          spr,
                                          static_cast<std::uint64_t>(duration));
  }

  double sample_rate(std::size_t id) const {
    return static_cast<double>(channel(id).spec.samples_per_record) *
           ticks_per_second / info_.record_duration.count();
  }
  std::size_t select_unique(std::string_view label) const {
    std::size_t found = info_.signals.size();
    for (std::size_t i = 0; i < info_.signals.size(); ++i)
      if (info_.signals[i].label == label) {
        if (found != info_.signals.size())
          throw Error("Signal label is ambiguous");
        found = i;
      }
    if (found == info_.signals.size())
      throw Error("Signal label not found");
    return found;
  }
  std::vector<std::size_t>
  select_unique(std::span<const std::string_view> labels) const {
    std::vector<std::size_t> result;
    result.reserve(labels.size());
    for (auto label : labels)
      result.push_back(select_unique(label));
    return result;
  }
  std::vector<std::size_t>
  select_unique(std::initializer_list<std::string_view> labels) const {
    return select_unique(
        std::span<const std::string_view>(labels.begin(), labels.size()));
  }
  Cursor cursor(std::size_t id) const { return Cursor(*this, id); }
  void seek(std::size_t id, std::uint64_t sample) {
    if (sample > sample_count(id))
      throw Error("Sample offset out of range");
    cursors_[id] = sample;
  }
  std::uint64_t tell(std::size_t id) const {
    channel(id);
    return cursors_[id];
  }
  std::size_t read_digital_at(std::size_t id, std::uint64_t first,
                              std::span<std::int32_t> out) const {
    return read_at(id, first, out);
  }
  std::size_t read_physical_at(std::size_t id, std::uint64_t first,
                               std::span<double> out) const {
    return read_at(id, first, out);
  }
  std::size_t read_digital(std::size_t id, std::span<std::int32_t> out) {
    auto n = read_at(id, tell(id), out);
    cursors_[id] += n;
    return n;
  }
  std::size_t read_physical(std::size_t id, std::span<double> out) {
    auto n = read_at(id, tell(id), out);
    cursors_[id] += n;
    return n;
  }
  std::vector<std::vector<std::int32_t>>
  read_digital(std::span<const std::size_t> ids, TimeRange range) const {
    return window<std::int32_t>(ids, range);
  }
  std::vector<std::vector<double>>
  read_physical(std::span<const std::size_t> ids, TimeRange range) const {
    return window<double>(ids, range);
  }
  // Results follow the requested signal order. Mixed-rate signals have separate
  // lengths; discontinuous gaps contribute no samples. The into variants reuse
  // the capacity of each output channel across repeated windows.
  void read_digital_into(std::span<const std::size_t> ids, TimeRange range,
                         std::vector<std::vector<std::int32_t>> &out) const {
    window_into(ids, range, out);
  }
  void read_physical_into(std::span<const std::size_t> ids, TimeRange range,
                          std::vector<std::vector<double>> &out) const {
    window_into(ids, range, out);
  }
};
} // namespace fastbdf
