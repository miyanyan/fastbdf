#pragma once

#include <algorithm>
#include <cmath>
#include <span>
#include <type_traits>
#include "../types.hpp"

namespace fastbdf::detail
{

template<class T>
inline void decode(unsigned width, const unsigned char* ptr, std::span<T> out, const SignalSpec& s)
{
    const double gain = (s.physical_max - s.physical_min) / (static_cast<double>(s.digital_max) - s.digital_min);
    for (auto& sample : out) {
        std::uint32_t v = std::uint32_t(ptr[0]) | (std::uint32_t(ptr[1]) << 8);
        if (width == 3) v |= std::uint32_t(ptr[2]) << 16;
        const auto raw = width == 2 ? static_cast<std::int32_t>(v ^ 0x8000) - 0x8000
                                    : static_cast<std::int32_t>(v ^ 0x800000) - 0x800000;
        const auto value = std::clamp(raw, s.digital_min, s.digital_max);
        if constexpr (std::is_same_v<T, std::int32_t>)
            sample = value;
        else if (value == s.digital_min)
            sample = s.physical_min;
        else if (value == s.digital_max)
            sample = s.physical_max;
        else if (std::isfinite(gain))
            sample = s.physical_min + (static_cast<double>(value) - s.digital_min) * gain;
        else
            sample = std::lerp(
                s.physical_min, s.physical_max,
                (static_cast<double>(value) - s.digital_min) / (static_cast<double>(s.digital_max) - s.digital_min));
        ptr += width;
    }
}

inline std::int32_t to_digital(DigitalSample auto x, const SignalSpec& s)
{
    return std::clamp<std::int32_t>(x, s.digital_min, s.digital_max);
}
inline std::int32_t to_digital(PhysicalSample auto sample, const SignalSpec& s)
{
    const double x = sample;
    if (!std::isfinite(x)) throw Error("non-finite physical sample");
    if (s.physical_min < s.physical_max) {
        if (x <= s.physical_min) return s.digital_min;
        if (x >= s.physical_max) return s.digital_max;
    }
    else {
        if (x >= s.physical_min) return s.digital_min;
        if (x <= s.physical_max) return s.digital_max;
    }
    const auto bitvalue = (s.physical_max - s.physical_min) / (double(s.digital_max) - s.digital_min);
    const auto offset = s.physical_max / bitvalue - s.digital_max;
    auto value = x / bitvalue - offset;
    if (!std::isfinite(value) || !std::isfinite(bitvalue) || bitvalue == 0) {
        const auto scale = (std::max)(std::abs(s.physical_min), std::abs(s.physical_max));
        const auto fraction = (x / scale - s.physical_min / scale) / (s.physical_max / scale - s.physical_min / scale);
        value = fraction * (double(s.digital_max) - s.digital_min) + s.digital_min;
    }
    return static_cast<std::int32_t>(std::clamp(value, double(s.digital_min), double(s.digital_max)));
}
template<std::size_t Width, class T>
inline void encode(std::span<const T> samples, const SignalSpec& signal, char* dest)
{
    for (auto sample : samples) {
        const auto bits = static_cast<std::uint32_t>(to_digital(sample, signal));
        dest[0] = static_cast<char>(bits);
        dest[1] = static_cast<char>(bits >> 8);
        if constexpr (Width == 3) dest[2] = static_cast<char>(bits >> 16);
        dest += Width;
    }
}
template<class T>
inline void encode_channel(std::size_t width, std::span<const T> samples, const SignalSpec& signal, char* dest)
{
    if (width == 3)
        encode<3>(samples, signal, dest);
    else
        encode<2>(samples, signal, dest);
}

} // namespace fastbdf::detail
