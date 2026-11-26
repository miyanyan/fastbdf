#pragma once

#include "../types.hpp"
#include <limits>

namespace fastbdf::detail {

inline constexpr std::int64_t ticks_per_second = Tick::period::den;

inline Tick ticks(std::string_view s, const char *field) {
  auto t = trim_spaces(s);
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
      if (digits < std::chrono::hh_mm_ss<Tick>::fractional_width) {
        fraction = fraction * 10 + d;
        ++digits;
      } // EDFlib truncates sub-tick fractional digits toward zero.
    }
  }
  if (!any || p != t.size())
    throw Error(std::string("Invalid ") + field);
  while (digits++ < std::chrono::hh_mm_ss<Tick>::fractional_width)
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

inline std::string seconds(std::int64_t ticks, bool sign) {
  const bool negative = ticks < 0;
  const auto n =
      negative ? std::uint64_t(-(ticks + 1)) + 1 : std::uint64_t(ticks);
  std::string out = negative ? "-" : (sign ? "+" : "");
  out += std::to_string(n / ticks_per_second);
  auto frac = n % ticks_per_second;
  if (frac) {
    auto s = std::to_string(frac + ticks_per_second).substr(1);
    while (s.back() == '0')
      s.pop_back();
    out += '.';
    out += s;
  }
  return out;
}

} // namespace fastbdf::detail
