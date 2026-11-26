#pragma once

#include <cstddef>

namespace fastbdf::detail::edflib_limits {

// EDFlib 1.27 counts description bytes, including for UTF-8 text.
inline constexpr std::size_t written_description_bytes = 40;
inline constexpr std::size_t read_description_bytes = 512;
inline constexpr std::size_t duration_text_chars = 15;

// The byte limit depends on the sample width (EDF 16-bit, BDF 24-bit).
inline constexpr std::size_t samples_per_record_limit = 5 * 1024 * 1024;

} // namespace fastbdf::detail::edflib_limits
