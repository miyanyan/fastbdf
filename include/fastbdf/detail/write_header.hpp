#pragma once

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <format>
#include "edflib_limits.hpp"
#include "header_layout.hpp"
#include "time.hpp"

namespace fastbdf::detail
{

struct WriteLayout
{
    std::size_t sample_bytes{};
    std::size_t data_bytes{};
    std::size_t record_bytes{};
    std::size_t header_bytes{};
    std::size_t annotation_offset{};
    std::size_t annotation_index{};
    bool plus{};
};

inline std::string field(std::string value, std::size_t width)
{
    if (value.size() > width) throw Error("EDF header field exceeds fixed width");
    for (unsigned char c : value)
        if (c < 32 || c > 126) throw Error("EDF header fields must be printable ASCII");
    value.resize(width, ' ');
    return value;
}
inline std::string decimal(double value)
{
    std::array<char, 64> buf{};
    for (int precision = 7; precision >= 0; --precision) {
        auto result = std::to_chars(buf.data(), buf.data() + buf.size(), value, std::chars_format::fixed, precision);
        if (result.ec == std::errc{} && result.ptr - buf.data() <= 8) {
            std::string out(buf.data(), result.ptr);
            if (out.find('.') != std::string::npos) {
                while (out.back() == '0')
                    out.pop_back();
                if (out.back() == '.') out.pop_back();
            }
            if ((out == "0" || out == "-0") && value != 0) continue;
            return out;
        }
    }
    for (int precision = 7; precision >= 1; --precision) {
        auto result = std::to_chars(buf.data(), buf.data() + buf.size(), value, std::chars_format::general, precision);
        if (result.ec == std::errc{} && result.ptr - buf.data() <= 8) return {buf.data(), result.ptr};
    }
    throw Error("numeric header value cannot fit eight bytes");
}
inline std::string token(std::string value)
{
    if (value.empty()) return "X";
    std::replace(value.begin(), value.end(), ' ', '_');
    return value;
}
inline std::string date(std::chrono::year_month_day d, bool extended)
{
    const auto y = int(d.year());
    const auto m = unsigned(d.month());
    const auto day = unsigned(d.day());
    if (extended) return std::format("{:02}-{}-{:04}", day, month_names[m - 1], y);
    return std::format("{:02}.{:02}.{:02}", day, m, y % 100);
}
inline std::string birthdate(std::string value)
{
    if (value.empty() || value == "X") return "X";
    if (value.size() != 10 || value[4] != '-' || value[7] != '-') throw Error("birthdate must be YYYY-MM-DD");
    auto number = [](std::string_view s) {
        unsigned n = 0;
        auto r = std::from_chars(s.data(), s.data() + s.size(), n);
        if (r.ec != std::errc{} || r.ptr != s.data() + s.size()) throw Error("invalid birthdate");
        return n;
    };
    std::chrono::year_month_day d{std::chrono::year(int(number(std::string_view(value).substr(0, 4)))),
                                  std::chrono::month(number(std::string_view(value).substr(5, 2))),
                                  std::chrono::day(number(std::string_view(value).substr(8, 2)))};
    if (!d.ok()) throw Error("invalid birthdate");
    return date(d, true);
}

inline WriteLayout prepare_layout(RecordingInfo& info, WriterOptions& options)
{
    WriteLayout layout;
    layout.sample_bytes = (info.format == Format::bdf || info.format == Format::bdf_plus) ? 3 : 2;
    layout.plus = info.format == Format::edf_plus || info.format == Format::bdf_plus;
    const auto max_record_bytes = layout.sample_bytes * edflib_limits::samples_per_record_limit;
    if ((info.signals.empty() && !layout.plus) || info.signals.size() > 4096) throw Error("invalid signal count");
    if (info.record_duration.count() < 0 ||
        (info.record_duration.count() == 0 && (!layout.plus || !info.signals.empty())) ||
        format_seconds(info.record_duration, false).size() > fixed_header::record_duration.width)
        throw Error("invalid data record duration");
    if (info.subsecond.count() < 0 || info.subsecond >= std::chrono::seconds{1} ||
        (!layout.plus && info.subsecond.count() != 0))
        throw Error("invalid subsecond start");
    if (layout.plus &&
        (!options.annotation_signals || options.annotation_signals > 64 || options.annotation_bytes_per_signal < 32 ||
         options.annotation_bytes_per_signal > max_record_bytes ||
         options.annotation_bytes_per_signal % layout.sample_bytes))
        throw Error("invalid annotation channel allocation");
    if (!layout.plus) options.annotation_signals = 0;
    for (auto& s : info.signals) {
        if (!s.samples_per_record || s.samples_per_record > max_decimal_count(signal_header::samples_per_record))
            throw Error("invalid samples per record");
        const auto lo = layout.sample_bytes == 2 ? -32768 : -8388608;
        const auto hi = layout.sample_bytes == 2 ? 32767 : 8388607;
        if (s.digital_min < lo || s.digital_max > hi || s.digital_min >= s.digital_max)
            throw Error("invalid digital limits");
        if (!std::isfinite(s.physical_min) || !std::isfinite(s.physical_max) || s.physical_min == s.physical_max)
            throw Error("invalid physical limits");
        const auto physical_min = decimal(s.physical_min);
        const auto physical_max = decimal(s.physical_max);
        if (physical_min == physical_max) throw Error("physical limits collapse in header precision");
        auto parse = [](const std::string& str) {
            double x{};
            auto r = std::from_chars(str.data(), str.data() + str.size(), x);
            if (r.ec != std::errc{}) throw Error("invalid calibration");
            return x;
        };
        s.physical_min = parse(physical_min);
        s.physical_max = parse(physical_max);
        const auto label = detail::trim_spaces(s.label);
        if (label == "EDF Annotations" || label == "BDF Annotations") throw Error("reserved annotation label");
        if (s.samples_per_record > (max_record_bytes - layout.data_bytes) / layout.sample_bytes)
            throw Error("data record exceeds the EDFlib size limit");
        layout.data_bytes += std::size_t(s.samples_per_record) * layout.sample_bytes;
    }
    layout.annotation_index = options.annotation_position == AnnotationPosition::start
                                  ? 0
                                  : (options.annotation_position == AnnotationPosition::middle ? info.signals.size() / 2
                                                                                               : info.signals.size());
    for (std::size_t i = 0; i < layout.annotation_index; ++i)
        layout.annotation_offset += std::size_t(info.signals[i].samples_per_record) * layout.sample_bytes;
    if (options.annotation_signals &&
        options.annotation_bytes_per_signal > (max_record_bytes - layout.data_bytes) / options.annotation_signals)
        throw Error("data record exceeds the EDFlib size limit");
    layout.record_bytes = layout.data_bytes + options.annotation_signals * options.annotation_bytes_per_signal;
    if (layout.record_bytes > max_record_bytes) throw Error("data record exceeds the EDFlib size limit");
    const auto ns = info.signals.size() + options.annotation_signals;
    if (ns > 4096) throw Error("More than 4096 signals including annotation channels");
    layout.header_bytes = header_size(ns);
    return layout;
}

inline std::string encode_header(const RecordingInfo& info, const WriterOptions& options, const WriteLayout& layout)
{
    const auto ns = info.signals.size() + options.annotation_signals;
    const auto days = std::chrono::floor<std::chrono::days>(info.start);
    const std::chrono::year_month_day ymd{days};
    if (!ymd.ok() || int(ymd.year()) < 1985 || int(ymd.year()) > 2084)
        throw Error("EDF start year must be 1985 through 2084");
    const auto tod = std::chrono::hh_mm_ss(info.start - days);
    std::string h(layout.header_bytes, ' ');
    h.replace(fixed_header::version.offset, fixed_header::version.width,
              layout.sample_bytes == 3 ? bdf_version : edf_version);
    auto put = [&](HeaderField spec, std::string value) {
        h.replace(spec.offset, spec.width, field(std::move(value), spec.width));
    };
    std::string patient = info.patient;
    std::string recording = info.recording;
    if (layout.plus) {
        if (!info.sex.empty() && info.sex != "M" && info.sex != "F" && info.sex != "X")
            throw Error("sex must be M, F or X");
    // EDFlib requires the separator after the patient code at byte 72 or
    // earlier.
        if (info.patient_code.size() > 72) throw Error("patient code exceeds the EDFlib header limit");
        patient = token(info.patient_code) + ' ' + token(info.sex) + ' ' + birthdate(info.birthdate) + ' ' +
                  token(info.patient_name);
        if (!info.patient_additional.empty()) patient += ' ' + info.patient_additional;
        recording = "Startdate " + date(ymd, true) + ' ' + token(info.admin_code) + ' ' + token(info.technician) + ' ' +
                    token(info.equipment);
        if (!info.recording_additional.empty()) recording += ' ' + info.recording_additional;
    }
    put(fixed_header::patient, patient);
    put(fixed_header::recording, recording);
    put(fixed_header::start_date, date(ymd, false));
    put(fixed_header::start_time,
        std::format("{:02}.{:02}.{:02}", tod.hours().count(), tod.minutes().count(), tod.seconds().count()));
    put(fixed_header::header_bytes, std::to_string(layout.header_bytes));
    put(fixed_header::reserved, layout.plus ? (layout.sample_bytes == 3 ? "BDF+C" : "EDF+C") : "");
    put(fixed_header::record_count, "-1");
    put(fixed_header::record_duration, format_seconds(info.record_duration, false));
    put(fixed_header::signal_count, std::to_string(ns));
    std::vector<SignalSpec> specs = info.signals;
    for (std::size_t i = 0; i < options.annotation_signals; ++i) {
        SignalSpec s{};
        s.label = layout.sample_bytes == 3 ? "BDF Annotations" : "EDF Annotations";
        s.physical_min = -1;
        s.physical_max = 1;
        s.digital_min = layout.sample_bytes == 3 ? -8388608 : -32768;
        s.digital_max = layout.sample_bytes == 3 ? 8388607 : 32767;
        s.samples_per_record = static_cast<std::uint32_t>(options.annotation_bytes_per_signal / layout.sample_bytes);
        specs.insert(specs.begin() + static_cast<std::ptrdiff_t>(layout.annotation_index + i), std::move(s));
    }
    auto fields = [&](HeaderField spec, auto get) {
        for (std::size_t i = 0; i < specs.size(); ++i)
            h.replace(fixed_header::size + spec.column_offset(ns, i), spec.width, field(get(specs[i]), spec.width));
    };
    fields(signal_header::label, [](const auto& s) { return s.label; });
    fields(signal_header::transducer, [](const auto& s) { return s.transducer; });
    fields(signal_header::unit, [](const auto& s) { return s.unit; });
    fields(signal_header::physical_min, [](const auto& s) { return decimal(s.physical_min); });
    fields(signal_header::physical_max, [](const auto& s) { return decimal(s.physical_max); });
    fields(signal_header::digital_min, [](const auto& s) { return std::to_string(s.digital_min); });
    fields(signal_header::digital_max, [](const auto& s) { return std::to_string(s.digital_max); });
    fields(signal_header::prefilter, [](const auto& s) { return s.prefilter; });
    fields(signal_header::samples_per_record, [](const auto& s) { return std::to_string(s.samples_per_record); });
    return h;
}

} // namespace fastbdf::detail
