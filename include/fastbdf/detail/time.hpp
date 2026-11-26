#pragma once

#include <array>
#include <limits>
#include "../types.hpp"

namespace fastbdf::detail
{

inline constexpr std::array<std::string_view, 12> month_names = {"JAN", "FEB", "MAR", "APR", "MAY", "JUN",
                                                                 "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"};

inline constexpr std::int64_t nanoseconds_per_second = std::chrono::nanoseconds::period::den;

inline std::chrono::nanoseconds parse_seconds(std::string_view s, const char* field)
{
    auto t = trim_spaces(s);
    std::size_t p = 0;
    bool negative = false;
    if (p < t.size() && (t[p] == '+' || t[p] == '-')) negative = t[p++] == '-';
    std::uint64_t whole = 0;
    std::uint64_t fraction = 0;
    unsigned digits = 0;
    constexpr auto fractional_width = std::chrono::hh_mm_ss<std::chrono::nanoseconds>::fractional_width;
    bool any = false;
    const auto limit = static_cast<std::uint64_t>((std::numeric_limits<std::int64_t>::max)()) + (negative ? 1u : 0u);
    while (p < t.size() && t[p] >= '0' && t[p] <= '9') {
        any = true;
        auto d = static_cast<unsigned>(t[p++] - '0');
        if (whole > (limit / nanoseconds_per_second - d) / 10) throw Error("Time value overflows");
        whole = whole * 10 + d;
    }
    if (p < t.size() && t[p] == '.') {
        ++p;
        while (p < t.size() && t[p] >= '0' && t[p] <= '9') {
            any = true;
            auto d = static_cast<unsigned>(t[p++] - '0');
            if (digits < fractional_width) {
                fraction = fraction * 10 + d;
                ++digits;
            } // Truncate sub-nanosecond digits toward zero.
        }
    }
    if (!any || p != t.size()) throw Error(std::string("Invalid ") + field);
    while (digits++ < fractional_width)
        fraction *= 10;
    if (whole > limit / nanoseconds_per_second || whole * nanoseconds_per_second > limit - fraction)
        throw Error("Time value overflows");
    const auto magnitude = whole * nanoseconds_per_second + fraction;
    if (negative && magnitude == limit) return (std::chrono::nanoseconds::min)();
    auto value = static_cast<std::int64_t>(magnitude);
    return std::chrono::nanoseconds{negative ? -value : value};
}

inline std::string format_seconds(std::chrono::nanoseconds duration, bool sign)
{
    const auto ns = duration.count();
    const bool negative = ns < 0;
    const auto n = negative ? std::uint64_t(-(ns + 1)) + 1 : std::uint64_t(ns);
    std::string out = negative ? "-" : (sign ? "+" : "");
    out += std::to_string(n / nanoseconds_per_second);
    auto frac = n % nanoseconds_per_second;
    if (frac) {
        auto s = std::to_string(frac + nanoseconds_per_second).substr(1);
        while (s.back() == '0')
            s.pop_back();
        out += '.';
        out += s;
    }
    return out;
}

} // namespace fastbdf::detail
