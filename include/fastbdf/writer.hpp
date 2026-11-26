#pragma once

#include "detail/sample_codec.hpp"
#include "detail/tal.hpp"
#include "detail/write_header.hpp"
#include "types.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <ranges>
#include <span>
#include <string>
#include <vector>

namespace fastbdf {

template <class R>
concept SampleRange = std::ranges::contiguous_range<const R> &&
                      std::ranges::sized_range<const R> &&
                      (DigitalSample<std::ranges::range_value_t<const R>> ||
                       PhysicalSample<std::ranges::range_value_t<const R>>) &&
                      requires(const R &r) {
                        std::span<const std::ranges::range_value_t<const R>>{r};
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

  template <SampleRange R> void write_samples(const R &samples) {
    write_samples(next_signal_, samples);
  }
  template <SampleRange R>
  void write_samples(std::size_t channel, const R &samples) {
    using T = std::ranges::range_value_t<const R>;
    write_channel(channel, std::span<const T>{samples});
  }
  void write_digital_record(std::span<const std::int32_t> samples) {
    write_record(samples);
  }
  void write_short_record(std::span<const std::int16_t> samples) {
    write_record(samples);
  }
  void write_physical_record(std::span<const double> samples) {
    write_record(samples);
  }
  void write_packed_record(std::span<const std::byte> bytes) {
    ensure_open();
    if (next_signal_ || bytes.size() != layout_.data_bytes)
      fail("packed record size or channel order mismatch");
    emit_record(reinterpret_cast<const char *>(bytes.data()));
  }
  void add_annotation(Annotation annotation) {
    ensure_open();
    if (!layout_.plus)
      fail("annotations require EDF+ or BDF+");
    if (annotation.duration && annotation.duration->count() < 0)
      fail("negative annotation duration");
    // Match EDFlib 1.27's byte-oriented description buffer.
    auto &text = annotation.text;
    text.resize(std::min(text.size(),
                         detail::edflib_limits::written_description_bytes));
    if (auto end = text.find('\0'); end != text.npos)
      text.resize(end);
    for (char &c : text)
      if (static_cast<unsigned char>(c) < 32 || c == '\x7f')
        c = '.';
    annotations_.push_back(std::move(annotation));
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
    const bool generate_records = info_.signals.empty() && records_ == 0;
    if (!records_ && !generate_records)
      fail("cannot finish without at least one data record");
    auto blocks = plan_annotation_blocks(generate_records);
    if (generate_records) {
      const auto count =
          blocks.empty()
              ? 1
              : (blocks.size() - 1) / options_.annotation_signals + 1;
      while (records_ < count)
        emit_record(nullptr);
    }
    for (std::size_t slot = 0; slot < blocks.size(); ++slot) {
      const auto record = slot / options_.annotation_signals,
                 channel = slot % options_.annotation_signals;
      output_.seekp(static_cast<std::streamoff>(
          layout_.header_bytes +
          static_cast<std::uint64_t>(record) * layout_.record_bytes +
          layout_.annotation_offset +
          channel * options_.annotation_bytes_per_signal));
      auto &block = blocks[slot];
      block.resize(options_.annotation_bytes_per_signal, '\0');
      output_.write(block.data(), static_cast<std::streamsize>(block.size()));
    }
    constexpr auto record_count = detail::fixed_header::record_count;
    output_.seekp(record_count.offset);
    const auto count =
        detail::field(std::to_string(records_), record_count.width);
    output_.write(count.data(), static_cast<std::streamsize>(count.size()));
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
  std::vector<char> record_buffer_;
  std::size_t channel_offset_{};
  detail::WriteLayout layout_;
  std::size_t next_signal_{};
  std::uint64_t records_ = 0;
  bool finished_ = false;
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
  // Slots are filled in file order. Plan them before emitting or patching
  // records.
  std::vector<std::string> plan_annotation_blocks(bool generate_records) const {
    std::vector<std::string> blocks;
    const auto block_bytes = options_.annotation_bytes_per_signal;
    for (const auto &annotation : annotations_) {
      auto tal = detail::annotation_tal(info_, annotation);
      if (tal.size() > block_bytes)
        fail("annotation exceeds channel capacity");
      while (blocks.empty() ||
             blocks.back().size() + tal.size() > block_bytes) {
        const auto slot = blocks.size();
        const auto record = slot / options_.annotation_signals;
        if (!generate_records && record >= records_)
          fail("annotation capacity exhausted; increase annotation channels or "
               "bytes before writing");
        if (record >=
            detail::max_decimal_count(detail::fixed_header::record_count))
          fail("EDF record count exceeds eight digits");
        auto &block = blocks.emplace_back();
        if (slot % options_.annotation_signals == 0) {
          block = detail::timekeeping(info_, record);
          if (generate_records && options_.annotation_signals == 1 &&
              record > 0 && block.size() + tal.size() > block_bytes)
            fail("annotation exceeds timekeeping channel capacity");
        }
      }
      blocks.back() += tal;
    }
    return blocks;
  }
  template <class T> void write_record(std::span<const T> samples) {
    ensure_open();
    if (next_signal_ ||
        samples.size() != layout_.data_bytes / layout_.sample_bytes)
      fail("record sample count or channel order mismatch");
    record_buffer_.resize(layout_.data_bytes);
    std::size_t pos = 0;
    for (const auto &signal : info_.signals) {
      detail::encode_channel(
          layout_.sample_bytes, samples.subspan(pos, signal.samples_per_record),
          signal, record_buffer_.data() + pos * layout_.sample_bytes);
      pos += signal.samples_per_record;
    }
    emit_record(record_buffer_.data());
  }
  template <class T>
  void write_channel(std::size_t channel, std::span<const T> samples) {
    ensure_open();
    if (channel != next_signal_ || channel >= info_.signals.size())
      fail("signals must be written in header order");
    const auto &s = info_.signals[channel];
    if (samples.size() != s.samples_per_record)
      fail("channel sample count mismatch");
    record_buffer_.resize(layout_.data_bytes);
    detail::encode_channel(layout_.sample_bytes, samples, s,
                           record_buffer_.data() + channel_offset_);
    if (channel + 1 == info_.signals.size()) {
      emit_record(record_buffer_.data());
      next_signal_ = 0;
      channel_offset_ = 0;
    } else {
      ++next_signal_;
      channel_offset_ += samples.size() * layout_.sample_bytes;
    }
  }
  void emit_record(const char *data) {
    if (records_ >=
        detail::max_decimal_count(detail::fixed_header::record_count))
      fail("EDF record count exceeds eight digits");
    auto block =
        layout_.plus ? detail::timekeeping(info_, records_) : std::string{};
    if (layout_.plus)
      block.resize(options_.annotation_bytes_per_signal *
                       options_.annotation_signals,
                   '\0');
    if (layout_.annotation_offset)
      output_.write(data,
                    static_cast<std::streamsize>(layout_.annotation_offset));
    if (layout_.plus)
      output_.write(block.data(), static_cast<std::streamsize>(block.size()));
    if (layout_.data_bytes > layout_.annotation_offset)
      output_.write(data + layout_.annotation_offset,
                    static_cast<std::streamsize>(layout_.data_bytes -
                                                 layout_.annotation_offset));
    check_io();
    ++records_;
  }
  Writer(const std::filesystem::path &path, RecordingInfo info,
         WriterOptions options)
      : info_(std::move(info)), options_(options),
        layout_(detail::prepare_layout(info_, options_)) {
    const auto header = detail::encode_header(info_, options_, layout_);
    output_.open(path, std::ios::binary | std::ios::trunc);
    if (!output_)
      throw Error("Cannot create recording file", Errc::io);
    output_.write(header.data(), static_cast<std::streamsize>(header.size()));
    check_io();
  }
};
} // namespace fastbdf
