#pragma once

#include "edflib_limits.hpp"
#include "time.hpp"
#include <algorithm>
#include <limits>
#include <span>

namespace fastbdf::detail {

inline Tick tal_time(std::string_view text, bool is_duration) {
  auto digits = text;
  if (!is_duration) {
    if (digits.empty() || (digits.front() != '+' && digits.front() != '-'))
      throw Error("Invalid annotation onset");
    digits.remove_prefix(1);
  }
  const auto dot = digits.find('.');
  if (digits.empty() || digits.front() == '.' || digits.back() == '.' ||
      (dot != digits.npos && digits.find('.', dot + 1) != digits.npos) ||
      digits.find_first_not_of("0123456789.") != digits.npos)
    throw Error("Invalid TAL time");
  // EDFlib 1.27 validates the full duration, then copies 15 characters
  // before converting it to ticks. Onsets use the complete spelling.
  return ticks(is_duration ? text.substr(0, edflib_limits::duration_text_chars)
                           : text,
               "TAL time");
}

// TAL bytes are stored consecutively, including all bytes of 16/24 bit words.
inline bool parse_tals(std::span<const char> bytes, bool first_channel,
                       AnnotationMode mode, Tick &record_onset,
                       std::vector<Annotation> &annotations) {
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
    if (sep == tal.npos)
      throw Error("Invalid annotation onset");
    auto onset = tal_time(tal.substr(0, sep), false);
    std::optional<Tick> duration;
    auto q = sep;
    if (tal[q] == '\x15') {
      auto e = tal.find('\x14', q + 1);
      if (e == tal.npos)
        throw Error("Invalid annotation duration");
      duration = tal_time(tal.substr(q + 1, e - q - 1), true);
      q = e;
    }
    if (first && first_channel) {
      if (q + 1 >= tal.size() || tal[q + 1] != '\x14' || duration)
        throw Error("Missing record timekeeping TAL");
      record_onset = onset;
      if (mode == AnnotationMode::none)
        return false;
      // Skip just the first empty description reserved for record
      // timekeeping. Subsequent descriptions in the same TAL may themselves
      // be events.
      ++q;
    }
    first = false;
    while (q + 1 < tal.size()) {
      if (tal[q] != '\x14')
        throw Error("Invalid annotation separator");
      auto e = tal.find('\x14', q + 1);
      if (e == tal.npos)
        throw Error("Unterminated annotation text");
      if (mode != AnnotationMode::none) {
        // EDFlib exposes and checks the first 512 description bytes.
        auto text = tal.substr(
            q + 1, std::min(e - q - 1, edflib_limits::read_description_bytes));
        if (!std::all_of(text.begin(), text.end(), [](unsigned char c) {
              return (c >= 32 && c != 127) || c == '\t' || c == '\n' ||
                     c == '\r';
            }))
          throw Error("Invalid annotation text");
        annotations.push_back({onset, duration, std::string(text)});
        if (mode == AnnotationMode::until_recording_ends &&
            text.starts_with("Recording ends"))
          return true;
      }
      q = e;
    }
  }
  if (first && first_channel)
    throw Error("Missing record timekeeping TAL");
  return false;
}

inline std::string timekeeping(const RecordingInfo &info,
                               std::uint64_t record) {
  const auto duration = info.record_duration.count();
  if (duration && record > static_cast<std::uint64_t>(
                               ((std::numeric_limits<std::int64_t>::max)() -
                                info.subsecond.count()) /
                               duration))
    throw Error("record time overflow");
  auto s = seconds(static_cast<std::int64_t>(record) * duration +
                       info.subsecond.count(),
                   true);
  s.append("\x14\x14\0", 3);
  return s;
}
inline std::string annotation_tal(const RecordingInfo &info,
                                  const Annotation &a) {
  if (a.onset.count() >
      (std::numeric_limits<std::int64_t>::max)() - info.subsecond.count())
    throw Error("annotation onset overflow");
  auto s = seconds(a.onset.count() + info.subsecond.count(), true);
  if (a.duration) {
    s += '\x15';
    s += seconds(a.duration->count(), false);
  }
  s += '\x14';
  s += a.text;
  s.append("\x14\0", 2);
  return s;
}

} // namespace fastbdf::detail
