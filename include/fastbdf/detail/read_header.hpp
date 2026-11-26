#pragma once

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <format>
#include <limits>
#include "edflib_limits.hpp"
#include "file_source.hpp"
#include "header_layout.hpp"
#include "time.hpp"

namespace fastbdf::detail
{

struct Channel
{
    SignalSpec spec;
    std::uint64_t offset{};
    bool annotation{};
};

struct ReadHeader
{
    RecordingInfo info;
    std::vector<Channel> channels;
    std::vector<std::size_t> data_channels;
    std::uint64_t records{};
    std::uint64_t header_bytes{};
    std::uint64_t record_bytes{};
    unsigned width{};
    bool discontinuous{};
};

inline std::string trim(std::string_view s)
{
    return std::string(detail::trim_spaces(s));
}
inline std::int64_t integer(std::string_view s, const char* field)
{
    auto t = trim(s);
    if (t.size() >= 2 && t[0] == '+' && t[1] >= '0' && t[1] <= '9') t.erase(t.begin());
    std::int64_t n{};
    auto [p, e] = std::from_chars(t.data(), t.data() + t.size(), n);
    if (e != std::errc{} || p != t.data() + t.size() || t.empty()) throw Error(std::string("Invalid ") + field);
    return n;
}
inline double number(std::string_view s, const char* field)
{
    auto t = trim(s);
    double n{};
    if (!t.empty() && t.front() == '+') t.erase(t.begin());
    auto [p, e] = std::from_chars(t.data(), t.data() + t.size(), n);
    if (e != std::errc{} || p != t.data() + t.size() || t.empty() || !std::isfinite(n))
        throw Error(std::string("Invalid ") + field);
    return n;
}

inline std::uint64_t product(std::uint64_t a, std::uint64_t b)
{
    if (b && a > (std::numeric_limits<std::uint64_t>::max)() / b) throw Error("File size overflows");
    return a * b;
}

// Only the required prefix fields use underscores and space delimiters.
// The remainder is free text; consume exactly one separator before it.
template<std::size_t Count>
inline std::array<std::string, Count> identification_fields(std::string_view& tail)
{
    std::array<std::string, Count> fields;
    for (auto& field : fields) {
        if (tail.empty() || tail.front() == ' ') throw Error("Missing EDF+ identification field");
        const auto end = tail.find(' ');
        field = tail.substr(0, end);
        tail.remove_prefix(end == tail.npos ? tail.size() : end + 1);
    }
    return fields;
}
inline std::string decoded(std::string s)
{
    if (s == "X") return {};
    std::replace(s.begin(), s.end(), '_', ' ');
    return s;
}
inline std::chrono::year_month_day parse_extended_date(std::string_view text, const char* field)
{
    if (text.size() != 11 || text[2] != '-' || text[6] != '-') throw Error(std::string("Invalid ") + field);
    const auto month = std::find(month_names.begin(), month_names.end(), text.substr(3, 3));
    if (month == month_names.end()) throw Error(std::string("Invalid ") + field);
    const auto date = std::chrono::year{static_cast<int>(integer(text.substr(7), field))} /
                      std::chrono::month{static_cast<unsigned>(month - month_names.begin() + 1)} /
                      std::chrono::day{static_cast<unsigned>(integer(text.substr(0, 2), field))};
    if (!date.ok()) throw Error(std::string("Invalid ") + field);
    return date;
}

inline void parse_identification(RecordingInfo& info)
{
    if (info.format != Format::edf_plus && info.format != Format::bdf_plus) return;
    std::string_view patient_tail = info.patient;
    std::string_view recording_tail = info.recording;
    auto p = identification_fields<4>(patient_tail);
    auto r = identification_fields<5>(recording_tail);
    if (r[0] != "Startdate") throw Error("Invalid EDF+ identification fields");
    info.patient_code = decoded(p[0]);
    info.sex = decoded(p[1]);
    if (!info.sex.empty() && info.sex != "M" && info.sex != "F") throw Error("Invalid patient sex");
    info.birthdate = decoded(p[2]);
    info.patient_name = decoded(p[3]);
    if (!info.birthdate.empty()) {
        const auto date = parse_extended_date(info.birthdate, "patient birthdate");
        const std::string_view text = info.birthdate;
        info.birthdate = std::format("{}-{:02}-{}", text.substr(7), unsigned(date.month()), text.substr(0, 2));
    }
    info.patient_additional = patient_tail;
    info.admin_code = decoded(r[2]);
    info.technician = decoded(r[3]);
    info.equipment = decoded(r[4]);
    info.recording_additional = recording_tail;
  // EDF+ uses the unambiguous four-digit date in the recording field.
    if (r[1] != "X") {
        const auto date = parse_extended_date(r[1], "EDF+ start date");
        auto original = std::chrono::year_month_day{std::chrono::floor<std::chrono::days>(info.start)};
        if (original.month() != date.month() || original.day() != date.day() ||
            int(original.year()) % 100 != int(date.year()) % 100)
            throw Error("EDF+ start dates disagree");
        auto tod = info.start - std::chrono::floor<std::chrono::days>(info.start);
        info.start = std::chrono::local_days{date} + tod;
    }
}

inline ReadHeader read_header(const FileSource& file)
{
    ReadHeader header;
    std::array<char, fixed_header::size> fixed{};
    file.read_at(0, fixed);
    std::string_view f(fixed.data(), fixed.size());
    const bool bdf = fixed_header::version.read(f) == bdf_version;
    if (!bdf && fixed_header::version.read(f) != edf_version) throw Error("Unsupported EDF/BDF version");
    header.width = bdf ? 3 : 2;
    for (std::size_t i = fixed_header::version.end(); i < f.size(); ++i)
        if (static_cast<unsigned char>(f[i]) < 32 || static_cast<unsigned char>(f[i]) > 126)
            throw Error("Non-ASCII EDF/BDF header");
    header.info.patient = trim(fixed_header::patient.read(f));
    header.info.recording = trim(fixed_header::recording.read(f));
    const auto reserved = trim(fixed_header::reserved.read(f));
    const auto plus = reserved.starts_with(bdf ? "BDF+" : "EDF+");
    if (plus && !(reserved.starts_with(bdf ? "BDF+C" : "EDF+C") || reserved.starts_with(bdf ? "BDF+D" : "EDF+D")))
        throw Error("Invalid EDF+ continuity marker");
    header.discontinuous = plus && reserved[4] == 'D';
    header.info.format = bdf ? (plus ? Format::bdf_plus : Format::bdf) : (plus ? Format::edf_plus : Format::edf);
    auto date = fixed_header::start_date.read(f);
    auto time = fixed_header::start_time.read(f);
    if (date[2] != '.' || date[5] != '.' || time[2] != '.' || time[5] != '.') throw Error("Invalid start date/time");
    auto y = integer(date.substr(6, 2), "year");
    y += y >= 85 ? 1900 : 2000;
    auto d = std::chrono::year{static_cast<int>(y)} /
             std::chrono::month{static_cast<unsigned>(integer(date.substr(3, 2), "month"))} /
             std::chrono::day{static_cast<unsigned>(integer(date.substr(0, 2), "day"))};
    auto hh = integer(time.substr(0, 2), "hour");
    auto mm = integer(time.substr(3, 2), "minute");
    auto ss = integer(time.substr(6, 2), "second");
    if (!d.ok() || hh < 0 || hh > 23 || mm < 0 || mm > 59 || ss < 0 || ss > 59) throw Error("Invalid start date/time");
    header.info.start =
        std::chrono::local_days{d} + std::chrono::hours{hh} + std::chrono::minutes{mm} + std::chrono::seconds{ss};
    auto ns = integer(fixed_header::signal_count.read(f), "signal count");
    if (ns <= 0 || ns > 4096) throw Error("Invalid signal count");
    auto hb = integer(fixed_header::header_bytes.read(f), "header length");
    if (hb < 0 || static_cast<std::uint64_t>(hb) != header_size(static_cast<std::size_t>(ns)))
        throw Error("Invalid header length");
    header.header_bytes = static_cast<std::uint64_t>(hb);
    const auto declared = integer(fixed_header::record_count.read(f), "record count");
    if (declared < -1) throw Error("Invalid record count");
    header.info.record_duration = parse_seconds(fixed_header::record_duration.read(f), "record duration");
    if (header.info.record_duration.count() < 0) throw Error("Record duration must be nonnegative");
    std::vector<char> raw(static_cast<std::size_t>(ns) * signal_header::size);
    file.read_at(fixed_header::size, raw);
    for (char c : raw)
        if (static_cast<unsigned char>(c) < 32 || static_cast<unsigned char>(c) > 126)
            throw Error("Non-ASCII signal header");
    auto field = [&](HeaderField spec, std::size_t i) {
        return std::string_view(raw.data() + spec.column_offset(static_cast<std::size_t>(ns), i), spec.width);
    };
    for (std::size_t i = 0; i < static_cast<std::size_t>(ns); ++i) {
        Channel c;
        auto& s = c.spec;
        s.label = trim(field(signal_header::label, i));
        s.transducer = trim(field(signal_header::transducer, i));
        s.unit = trim(field(signal_header::unit, i));
        s.physical_min = number(field(signal_header::physical_min, i), "physical minimum");
        s.physical_max = number(field(signal_header::physical_max, i), "physical maximum");
        auto mn = integer(field(signal_header::digital_min, i), "digital minimum");
        auto mx = integer(field(signal_header::digital_max, i), "digital maximum");
        auto limit = bdf ? 8388608 : 32768;
        if (mn < -limit || mx >= limit || mn >= mx) throw Error("Invalid digital range");
        s.digital_min = static_cast<std::int32_t>(mn);
        s.digital_max = static_cast<std::int32_t>(mx);
        if (s.physical_min == s.physical_max) throw Error("Invalid physical range");
        s.prefilter = trim(field(signal_header::prefilter, i));
        auto samples = integer(field(signal_header::samples_per_record, i), "samples per record");
        if (samples <= 0 || samples > (std::numeric_limits<std::uint32_t>::max)())
            throw Error("Invalid samples per record");
        s.samples_per_record = static_cast<std::uint32_t>(samples);
        c.offset = header.record_bytes;
        header.record_bytes += product(static_cast<std::uint64_t>(samples), header.width);
        c.annotation = plus && s.label == (bdf ? "BDF Annotations" : "EDF Annotations");
        if (!c.annotation) {
            header.data_channels.push_back(header.channels.size());
            header.info.signals.push_back(s);
        }
        header.channels.push_back(std::move(c));
    }
    if (header.info.record_duration.count() == 0 && (!plus || !header.data_channels.empty()))
        throw Error("Zero duration requires annotation-only EDF+");
    if (header.record_bytes > header.width * edflib_limits::samples_per_record_limit)
        throw Error("Data record exceeds the EDFlib size limit");
    const auto size = file.size();
    if (size < header.header_bytes || (size - header.header_bytes) % header.record_bytes)
        throw Error("Incomplete data record");
    header.records = (size - header.header_bytes) / header.record_bytes;
    if (declared >= 0 && static_cast<std::uint64_t>(declared) != header.records)
        throw Error("Record count does not match file size");
    parse_identification(header.info);
    return header;
}

} // namespace fastbdf::detail
