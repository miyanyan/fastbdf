#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>

namespace fastbdf::detail
{

struct HeaderField
{
    std::size_t offset;
    std::size_t width;

    constexpr std::size_t end() const { return offset + width; }
    constexpr std::string_view read(std::string_view header) const { return header.substr(offset, width); }

  // Signal fields are stored column by column, with one entry per signal.
    constexpr std::size_t column_offset(std::size_t signals, std::size_t index) const
    {
        return offset * signals + width * index;
    }
};

namespace fixed_header
{
inline constexpr HeaderField version{0, 8};
inline constexpr HeaderField patient{version.end(), 80};
inline constexpr HeaderField recording{patient.end(), 80};
inline constexpr HeaderField start_date{recording.end(), 8};
inline constexpr HeaderField start_time{start_date.end(), 8};
inline constexpr HeaderField header_bytes{start_time.end(), 8};
inline constexpr HeaderField reserved{header_bytes.end(), 44};
inline constexpr HeaderField record_count{reserved.end(), 8};
inline constexpr HeaderField record_duration{record_count.end(), 8};
inline constexpr HeaderField signal_count{record_duration.end(), 4};
inline constexpr auto size = signal_count.end();
static_assert(size == 256);
} // namespace fixed_header

namespace signal_header
{
inline constexpr HeaderField label{0, 16};
inline constexpr HeaderField transducer{label.end(), 80};
inline constexpr HeaderField unit{transducer.end(), 8};
inline constexpr HeaderField physical_min{unit.end(), 8};
inline constexpr HeaderField physical_max{physical_min.end(), 8};
inline constexpr HeaderField digital_min{physical_max.end(), 8};
inline constexpr HeaderField digital_max{digital_min.end(), 8};
inline constexpr HeaderField prefilter{digital_max.end(), 80};
inline constexpr HeaderField samples_per_record{prefilter.end(), 8};
inline constexpr HeaderField reserved{samples_per_record.end(), 32};
inline constexpr auto size = reserved.end();
static_assert(size == 256);
} // namespace signal_header

constexpr std::size_t header_size(std::size_t signals)
{
    return fixed_header::size + signal_header::size * signals;
}

constexpr std::uint64_t max_decimal_count(HeaderField field)
{
    std::uint64_t value = 0;
    for (std::size_t digit = 0; digit < field.width; ++digit)
        value = value * 10 + 9;
    return value;
}

inline constexpr std::string_view bdf_version =
    "\xff"
    "BIOSEMI";
inline constexpr std::string_view edf_version = "0       ";
static_assert(bdf_version.size() == fixed_header::version.width);
static_assert(edf_version.size() == fixed_header::version.width);

} // namespace fastbdf::detail
