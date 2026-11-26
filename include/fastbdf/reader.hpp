#pragma once
#include "detail/read_header.hpp"
#include "detail/sample_codec.hpp"
#include "detail/tal.hpp"
#include "types.hpp"
#include <algorithm>
#include <array>
#include <concepts>
#include <filesystem>
#include <initializer_list>
#include <ios>
#include <limits>
#include <memory>
#include <mutex>
#include <span>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace fastbdf {
template <class R>
concept SampleOutput = std::convertible_to<R, std::span<std::int32_t>> ||
                       std::convertible_to<R, std::span<double>>;

// Position-independent reads and independent Cursor objects may run
// concurrently. Reopening the reader and its legacy implicit cursors require
// synchronization.
class Reader {
  detail::ReadHeader header_;
  std::optional<detail::FileSource> file_;
  std::vector<std::uint64_t> cursors_;
  mutable std::vector<Annotation> annotations_;
  std::vector<Tick> record_onsets_;
  mutable bool annotations_loaded_ = true;
  std::unique_ptr<std::mutex> annotation_mutex_ =
      std::make_unique<std::mutex>();
  void normalize_annotations() const {
    for (auto &annotation : annotations_) {
      if (annotation.onset.count() <
          (std::numeric_limits<std::int64_t>::min)() +
              header_.info.subsecond.count())
        throw Error("Annotation onset overflows");
      annotation.onset -= header_.info.subsecond;
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
  const detail::Channel &channel(std::size_t signal) const {
    if (signal >= header_.data_channels.size())
      throw Error("Signal index out of range");
    return header_.channels[header_.data_channels[signal]];
  }
  template <class T>
  std::size_t read_at_impl(std::size_t signal, std::uint64_t first,
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
    if (first_record == last_record && count <= 1024 / header_.width) {
      std::array<char, 1024> bytes;
      read_bytes(header_.header_bytes + first_record * header_.record_bytes +
                     ch.offset + (first % spr) * header_.width,
                 std::span<char>(bytes).first(count * header_.width));
      detail::decode(header_.width,
                     reinterpret_cast<const unsigned char *>(bytes.data()),
                     out.first(count), ch.spec);
      return count;
    }
    if (first_record == last_record ||
        detail::product(spr, header_.width) * 4 < header_.record_bytes) {
      std::vector<char> bytes;
      std::size_t written = 0;
      for (auto r = first_record; r <= last_record; ++r) {
        const auto lo = r == first_record ? first % spr : 0;
        const auto hi = r == last_record ? (first + count - 1) % spr + 1 : spr;
        bytes.resize(static_cast<std::size_t>(hi - lo) * header_.width);
        read_bytes(header_.header_bytes + r * header_.record_bytes + ch.offset +
                       lo * header_.width,
                   bytes);
        detail::decode(header_.width,
                       reinterpret_cast<const unsigned char *>(bytes.data()),
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
        std::max<std::uint64_t>(1, target / header_.record_bytes);
    std::vector<char> buffer;
    std::size_t written = 0;
    for (auto r = first_record; r <= last_record;) {
      const auto n = (std::min)(batch_records, last_record - r + 1);
      const auto size = detail::product(n, header_.record_bytes);
      if (size > buffer.max_size())
        throw Error("Data record too large");
      buffer.resize(static_cast<std::size_t>(size));
      read_bytes(header_.header_bytes +
                     detail::product(r, header_.record_bytes),
                 buffer);
      for (std::uint64_t j = 0; j < n; ++j) {
        const auto record = r + j;
        const auto lo = record == first_record ? first % spr : 0;
        const auto hi =
            record == last_record ? (first + count - 1) % spr + 1 : spr;
        auto *ptr = reinterpret_cast<const unsigned char *>(
            buffer.data() + j * header_.record_bytes + ch.offset +
            lo * header_.width);
        detail::decode(header_.width, ptr,
                       out.subspan(written, static_cast<std::size_t>(hi - lo)),
                       ch.spec);
        written += static_cast<std::size_t>(hi - lo);
      }
      r += n;
    }
    return written;
  }
  using SampleBounds = std::pair<std::uint64_t, std::uint64_t>;
  template <class T>
  static void prepare_window(std::vector<std::vector<T>> &out,
                             const std::vector<SampleBounds> &bounds) {
    out.resize(bounds.size());
    for (std::size_t i = 0; i < bounds.size(); ++i)
      out[i].resize(
          static_cast<std::size_t>(bounds[i].second - bounds[i].first));
  }
  template <class T>
  static void prepare_window(std::span<std::span<T>> out,
                             const std::vector<SampleBounds> &bounds) {
    if (out.size() != bounds.size())
      throw Error("Window output must have one row per requested signal");
    for (std::size_t i = 0; i < bounds.size(); ++i)
      if (out[i].size() != bounds[i].second - bounds[i].first)
        throw Error("Window output length differs from native sample count");
  }
  template <class Rows>
  void window_into(std::span<const std::size_t> signals, TimeRange range,
                   Rows &result) const {
    using T = typename Rows::value_type::value_type;
    if (range.begin.count() < 0 || range.end < range.begin)
      throw Error("Invalid time range");
    std::vector<SampleBounds> bounds;
    bounds.reserve(signals.size());
    std::uint64_t selected_bytes = 0, first_record = header_.records,
                  end_record = 0;
    for (auto id : signals) {
      const auto spr = channel(id).spec.samples_per_record;
      const auto lo = sample_index_at_or_after(id, range.begin);
      const auto hi = sample_index_at_or_after(id, range.end);
      if (hi - lo > std::vector<T>{}.max_size())
        throw Error("Requested window too large");
      bounds.emplace_back(lo, hi);
      selected_bytes += detail::product(spr, header_.width);
      if (lo < hi) {
        first_record = (std::min)(first_record, lo / spr);
        end_record = (std::max)(end_record, (hi - 1) / spr + 1);
      }
    }
    prepare_window(result, bounds);
    if (first_record >= end_record)
      return;
    if (signals.size() == 1 || selected_bytes < header_.record_bytes / 4) {
      for (std::size_t i = 0; i < signals.size(); ++i)
        read_at(signals[i], bounds[i].first, std::span<T>(result[i]));
      return;
    }
    const auto batch =
        std::max<std::uint64_t>(1, 8 * 1024 * 1024 / header_.record_bytes);
    std::vector<char> buffer;
    for (auto r = first_record; r < end_record;) {
      const auto n = (std::min)(batch, end_record - r);
      buffer.resize(
          static_cast<std::size_t>(detail::product(n, header_.record_bytes)));
      read_bytes(header_.header_bytes + r * header_.record_bytes, buffer);
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
              buffer.data() + (rec - r) * header_.record_bytes + ch.offset +
              (begin - rec * spr) * header_.width);
          detail::decode(header_.width, ptr,
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
    template <SampleOutput R> std::size_t read(R &&out) {
      auto n = reader_->read_at(id_, position_, std::forward<R>(out));
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
  bool read_annotation_record(std::uint64_t record, AnnotationMode mode,
                              std::vector<char> &bytes, Tick &onset) const {
    bool first = true;
    for (const auto &channel : header_.channels) {
      if (!channel.annotation)
        continue;
      bytes.resize(static_cast<std::size_t>(channel.spec.samples_per_record) *
                   header_.width);
      read_bytes(header_.header_bytes + record * header_.record_bytes +
                     channel.offset,
                 bytes);
      if (detail::parse_tals(bytes, first, mode, onset, annotations_))
        return true;
      if (mode == AnnotationMode::none)
        break;
      first = false;
    }
    return false;
  }
  void validate_next_onset(Tick previous, Tick onset) const {
    const Tick duration = header_.info.record_duration;
    if (previous > (Tick::max)() - duration)
      throw Error("Record onset overflows");
    const Tick expected = previous + duration;
    if (header_.discontinuous ? onset < expected : onset != expected)
      throw Error("Invalid record onset sequence");
  }
  void open_impl(const std::filesystem::path &path, AnnotationMode mode) {
    file_.emplace(path);
    header_ = detail::read_header(*file_);
    if (header_.discontinuous) {
      if (header_.records > record_onsets_.max_size())
        throw Error("Too many data records");
      record_onsets_.resize(static_cast<std::size_t>(header_.records));
    }
    cursors_.resize(header_.data_channels.size());
    if (header_.info.format == Format::edf_plus ||
        header_.info.format == Format::bdf_plus) {
      if (std::none_of(header_.channels.begin(), header_.channels.end(),
                       [](const detail::Channel &c) { return c.annotation; }))
        throw Error("EDF+ requires an annotation signal");
      // EDFlib only honors the stop marker for a single annotation channel.
      if (mode == AnnotationMode::until_recording_ends &&
          header_.channels.size() - header_.data_channels.size() > 1)
        mode = AnnotationMode::all;
      annotations_loaded_ = mode != AnnotationMode::lazy;
      bool stopped = false;
      std::vector<char> bytes;
      Tick previous{};
      const auto scan_records =
          (!header_.discontinuous &&
           (mode == AnnotationMode::none || mode == AnnotationMode::lazy))
              ? std::min<std::uint64_t>(1, header_.records)
              : header_.records;
      for (std::uint64_t r = 0; r < scan_records; ++r) {
        Tick onset{};
        const auto scan_mode = stopped || mode == AnnotationMode::lazy
                                   ? AnnotationMode::none
                                   : mode;
        if (read_annotation_record(r, scan_mode, bytes, onset))
          stopped = true;
        if (r == 0) {
          header_.info.subsecond = onset;
          if (header_.info.subsecond.count() < 0 ||
              header_.info.subsecond.count() >= detail::ticks_per_second)
            throw Error("Invalid first record onset");
        }
        if (r)
          validate_next_onset(previous, onset);
        if (header_.discontinuous)
          record_onsets_[r] = onset;
        previous = onset;
        if (stopped && !header_.discontinuous)
          break;
      }
      normalize_annotations();
    }
  }

public:
  const RecordingInfo &header() const noexcept { return header_.info; }
  const std::vector<SignalSpec> &signals() const noexcept {
    return header_.info.signals;
  }
  const std::vector<Annotation> &annotations() const {
    std::lock_guard lock(*annotation_mutex_);
    if (!annotations_loaded_) {
      std::uint64_t record = 0;
      try {
        std::vector<char> bytes;
        Tick previous{};
        for (std::uint64_t r = 0; r < header_.records; ++r) {
          record = r;
          Tick onset{};
          read_annotation_record(r, AnnotationMode::all, bytes, onset);
          if (r)
            validate_next_onset(previous, onset);
          previous = onset;
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
  std::uint64_t record_count() const noexcept { return header_.records; }
  bool discontinuous() const noexcept { return header_.discontinuous; }
  std::uint64_t sample_count(std::size_t id) const {
    return detail::product(header_.records,
                           channel(id).spec.samples_per_record);
  }
  // Elapsed from the first record, including any discontinuities.
  Tick record_start(std::uint64_t record) const {
    if (record >= header_.records)
      throw Error("Record index out of range");
    if (header_.discontinuous)
      return record_onsets_[static_cast<std::size_t>(record)] -
             header_.info.subsecond;
    const auto duration = header_.info.record_duration.count();
    if (duration &&
        record > static_cast<std::uint64_t>((Tick::max)().count() / duration))
      throw Error("Record start exceeds Tick range");
    return Tick{static_cast<std::int64_t>(record) * duration};
  }

  // A time inside a gap maps to the next sample. Beyond EOF returns
  // sample_count(id).
  std::uint64_t sample_index_at_or_after(std::size_t id, Tick time) const {
    const auto spr = channel(id).spec.samples_per_record;
    if (time < Tick::zero())
      throw Error("Negative sample lookup time");
    if (!header_.records)
      return 0;
    const auto duration = header_.info.record_duration.count();
    std::uint64_t record;
    if (header_.discontinuous) {
      auto next = std::upper_bound(record_onsets_.begin(), record_onsets_.end(),
                                   time, [&](Tick t, Tick onset) {
                                     return t < onset - header_.info.subsecond;
                                   });
      record = static_cast<std::uint64_t>(next - record_onsets_.begin() - 1);
    } else {
      record = static_cast<std::uint64_t>(time.count() / duration);
    }
    if (record >= header_.records)
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
           detail::ticks_per_second / header_.info.record_duration.count();
  }
  std::size_t select_unique(std::string_view label) const {
    std::size_t found = header_.info.signals.size();
    for (std::size_t i = 0; i < header_.info.signals.size(); ++i)
      if (header_.info.signals[i].label == label) {
        if (found != header_.info.signals.size())
          throw Error("Signal label is ambiguous");
        found = i;
      }
    if (found == header_.info.signals.size())
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
  std::size_t read_at(std::size_t id, std::uint64_t first,
                      std::span<std::int32_t> out) const {
    return read_at_impl(id, first, out);
  }
  std::size_t read_at(std::size_t id, std::uint64_t first,
                      std::span<double> out) const {
    return read_at_impl(id, first, out);
  }
  template <SampleOutput R> std::size_t read(std::size_t id, R &&out) {
    auto n = read_at(id, tell(id), std::forward<R>(out));
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
  // Read directly into caller-owned rows. Validate all lengths before decoding.
  void read_physical_into(std::span<const std::size_t> ids, TimeRange range,
                          std::span<std::span<double>> out) const {
    window_into(ids, range, out);
  }
};
} // namespace fastbdf
