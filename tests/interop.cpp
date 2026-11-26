#include "edflib.h"
#include <array>
#include <cmath>
#include <fastbdf/fastbdf.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <vector>
using namespace fastbdf;
static void check(bool v, const char *message) {
  if (!v)
    throw std::runtime_error(message);
}
static void ok(int v) { check(v == 0, "upstream call failed"); }
static RecordingInfo info(bool bdf) {
  RecordingInfo h;
  h.format = bdf ? Format::bdf_plus : Format::edf_plus;
  h.record_duration = Tick{5000000};
  h.subsecond = Tick{1234500};
  h.patient_code = "P123";
  h.patient_name = "Test Person";
  h.sex = "F";
  h.equipment = "fastbdf";
  for (unsigned n : {7u, 3u}) {
    SignalSpec s;
    s.label = "channel" + std::to_string(h.signals.size());
    s.unit = "uV";
    s.samples_per_record = n;
    s.physical_min = -100;
    s.physical_max = 100;
    s.digital_min = bdf ? -8388608 : -32768;
    s.digital_max = bdf ? 8388607 : 32767;
    h.signals.push_back(s);
  }
  return h;
}
static std::vector<std::int32_t> values(const SignalSpec &s, unsigned r) {
  std::vector<std::int32_t> v(s.samples_per_record);
  for (unsigned i = 0; i < v.size(); ++i)
    v[i] = i % 3 == 0   ? s.digital_min
           : i % 3 == 1 ? s.digital_max
                        : int(r * 31 + i) - 100;
  return v;
}
static void create_fast(const std::filesystem::path &p, const RecordingInfo &h,
                        const std::vector<Annotation> &annotations = {
                            {Tick{2500000}, Tick{1250000}, "event \xCE\xB1"},
                            {Tick{9000000}, std::nullopt, "end"}}) {
  auto w = Writer::create(p, h);
  for (const auto &a : annotations)
    w.add_annotation(a);
  for (unsigned r = 0; r < 4; ++r) {
    std::vector<std::int32_t> all;
    for (auto &s : h.signals) {
      auto v = values(s, r);
      all.insert(all.end(), v.begin(), v.end());
    }
    w.write_digital_record(all);
  }
  w.finish();
}
static void
create_upstream(const std::filesystem::path &p, const RecordingInfo &h,
                bool bdf,
                const std::vector<Annotation> &annotations = {
                    {Tick{2500000}, Tick{1250000}, "event \xCE\xB1"},
                    {Tick{9000000}, std::nullopt, "end"}}) {
  int f = edfopen_file_writeonly(
      p.string().c_str(),
      bdf ? EDFLIB_FILETYPE_BDFPLUS : EDFLIB_FILETYPE_EDFPLUS, 2);
  check(f >= 0, "upstream create");
  try {
    ok(edf_set_startdatetime(f, 2000, 1, 1, 0, 0, 0));
    ok(edf_set_datarecord_duration(f, 50000));
    ok(edf_set_subsecond_starttime(f, 1234500));
    ok(edf_set_patientcode(f, "P123"));
    ok(edf_set_patientname(f, "Test Person"));
    ok(edf_set_sex(f, 0));
    ok(edf_set_equipment(f, "fastbdf"));
    ok(edf_set_patient_additional(f, h.patient_additional.c_str()));
    ok(edf_set_recording_additional(f, h.recording_additional.c_str()));
    for (int i = 0; i < 2; ++i) {
      auto &s = h.signals[i];
      ok(edf_set_samplefrequency(f, i, s.samples_per_record));
      ok(edf_set_label(f, i, s.label.c_str()));
      ok(edf_set_physical_dimension(f, i, s.unit.c_str()));
      ok(edf_set_physical_minimum(f, i, s.physical_min));
      ok(edf_set_physical_maximum(f, i, s.physical_max));
      ok(edf_set_digital_minimum(f, i, s.digital_min));
      ok(edf_set_digital_maximum(f, i, s.digital_max));
    }
    for (const auto &a : annotations)
      ok(edfwrite_annotation_utf8_hr(f, a.onset.count() / 10,
                                     a.duration ? a.duration->count() / 10 : -1,
                                     a.text.c_str()));
    for (unsigned r = 0; r < 4; ++r)
      for (auto &s : h.signals) {
        auto v = values(s, r);
        ok(edfwrite_digital_samples(f, v.data()));
      }
    ok(edfclose_file(f));
  } catch (...) {
    edfclose_file(f);
    throw;
  }
}
static void verify(const std::filesystem::path &p, const RecordingInfo &h) {
  Reader reader(p);
  auto hdr = std::make_unique<edflib_hdr_t>();
  ok(edfopen_file_readonly(p.string().c_str(), hdr.get(),
                           EDFLIB_READ_ALL_ANNOTATIONS));
  int f = hdr->handle;
  try {
    check(hdr->edfsignals == 2 && reader.signals().size() == 2, "signal count");
    check(reader.record_count() == 4 && hdr->datarecords_in_file == 4,
          "record count");
    check(reader.header().subsecond.count() == hdr->starttime_subsecond,
          "subsecond");
    check(reader.header().record_duration.count() == hdr->datarecord_duration,
          "record duration");
    check(reader.header().patient_code == h.patient_code &&
              reader.header().patient_name == h.patient_name &&
              reader.header().equipment == h.equipment,
          "structured patient metadata");
    check(std::string(hdr->patientcode) == h.patient_code &&
              std::string(hdr->patient_name) == h.patient_name &&
              std::string(hdr->equipment) == h.equipment,
          "upstream structured metadata");
    for (unsigned id = 0; id < 2; ++id) {
      auto &spec = reader.signals()[id];
      auto &ref = hdr->signalparam[id];
      check(spec.label == h.signals[id].label &&
                spec.unit == h.signals[id].unit,
            "signal metadata");
      check(spec.samples_per_record == ref.smp_in_datarecord &&
                spec.digital_min == ref.dig_min &&
                spec.digital_max == ref.dig_max &&
                spec.physical_min == ref.phys_min &&
                spec.physical_max == ref.phys_max,
            "signal scaling differential");
      auto &s = h.signals[id];
      std::vector<std::int32_t> expected;
      for (unsigned r = 0; r < 4; ++r) {
        auto v = values(s, r);
        expected.insert(expected.end(), v.begin(), v.end());
      }
      std::vector<std::int32_t> a(expected.size());
      std::vector<int> b(expected.size());
      check(reader.read_at(id, 0, a) == a.size(), "fast read count");
      check(edfread_digital_samples(f, id, int(b.size()), b.data()) ==
                int(b.size()),
            "upstream read count");
      check(a == expected && std::equal(a.begin(), a.end(), b.begin()),
            "digital differential");
      std::vector<double> pa(a.size()), pb(a.size());
      reader.read_at(id, 0, pa);
      check(edfseek(f, id, 0, EDFSEEK_SET) == 0, "seek");
      check(edfread_physical_samples(f, id, int(pb.size()), pb.data()) ==
                int(pb.size()),
            "physical count");
      for (unsigned k = 0; k < pa.size(); ++k)
        check(std::abs(pa[k] - pb[k]) < 1e-10, "physical differential");
      std::array<std::int32_t, 5> slice{};
      reader.read_at(id, s.samples_per_record - 1, slice);
      check(std::equal(slice.begin(), slice.end(),
                       expected.begin() + s.samples_per_record - 1),
            "cross-record partial read");
    }
    std::array<std::size_t, 2> ids{1, 0};
    auto window = reader.read_digital(ids, {Tick{1250000}, Tick{8750000}});
    for (unsigned i = 0; i < ids.size(); ++i) {
      auto &s = h.signals[ids[i]];
      std::vector<std::int32_t> expected;
      for (unsigned r = 0; r < 4; ++r) {
        auto v = values(s, r);
        expected.insert(expected.end(), v.begin(), v.end());
      }
      auto lo =
               static_cast<std::size_t>(std::ceil(0.25 * s.samples_per_record)),
           hi =
               static_cast<std::size_t>(std::ceil(1.75 * s.samples_per_record));
      check(window[i] == std::vector<std::int32_t>(expected.begin() + lo,
                                                   expected.begin() + hi),
            "mixed-rate half-open time window");
    }
    check(reader.annotations().size() == 2 && hdr->annotations_in_file == 2,
          "annotation count");
    for (int i = 0; i < 2; ++i) {
      edflib_annotation_t a{};
      ok(edf_get_annotation(f, i, &a));
      auto &b = reader.annotations()[i];
      check(b.text == a.annotation, "annotation text");
      check(b.onset.count() == a.onset, "annotation onset");
      check(
          (b.duration ? b.duration->count() == a.duration_l : a.duration_l < 0),
          "annotation duration");
    }
    ok(edfclose_file(f));
  } catch (...) {
    edfclose_file(f);
    throw;
  }
}
static void conversion_tests(const std::filesystem::path &dir, bool bdf) {
  auto h = info(bdf);
  h.subsecond = Tick{};
  h.record_duration = Tick{10000000};
  for (auto &s : h.signals) {
    s.digital_min = -100;
    s.digital_max = 200;
  }
  auto f = dir / (bdf ? "physical-fast.bdf" : "physical-fast.edf"),
       u = dir / (bdf ? "physical-upstream.bdf" : "physical-upstream.edf");
  std::vector<double> samples{-200, -100,  -99.99, -1,  0,
                              1,    99.99, 100,    200, 0.0001};
  {
    auto w = Writer::create(f, h);
    w.write_physical_record(samples);
    w.finish();
  }
  int handle = edfopen_file_writeonly(
      u.string().c_str(),
      bdf ? EDFLIB_FILETYPE_BDFPLUS : EDFLIB_FILETYPE_EDFPLUS, 2);
  check(handle >= 0, "physical upstream create");
  for (unsigned c = 0; c < 2; ++c) {
    auto &s = h.signals[c];
    ok(edf_set_samplefrequency(handle, c, s.samples_per_record));
    ok(edf_set_physical_minimum(handle, c, s.physical_min));
    ok(edf_set_physical_maximum(handle, c, s.physical_max));
    ok(edf_set_digital_minimum(handle, c, s.digital_min));
    ok(edf_set_digital_maximum(handle, c, s.digital_max));
  }
  ok(edfwrite_physical_samples(handle, samples.data()));
  ok(edfwrite_physical_samples(handle, samples.data() + 7));
  ok(edfclose_file(handle));
  Reader a(f), b(u);
  for (unsigned c = 0; c < 2; ++c) {
    std::vector<std::int32_t> x(h.signals[c].samples_per_record), y(x.size());
    a.read_at(c, 0, x);
    b.read_at(c, 0, y);
    check(x == y, "physical write clipping/truncation differential");
  }
  // Packed words preserve signed endpoints in their little-endian 16/24-bit
  // representation.
  h = info(bdf);
  auto packed = dir / (bdf ? "packed.bdf" : "packed.edf");
  std::vector<std::byte> bytes;
  unsigned width = bdf ? 3 : 2;
  for (auto &s : h.signals)
    for (auto x : values(s, 0))
      for (unsigned k = 0; k < width; ++k)
        bytes.push_back(static_cast<std::byte>(
            (static_cast<std::uint32_t>(x) >> (8 * k)) & 255));
  {
    auto w = Writer::create(packed, h);
    for (unsigned r = 0; r < 4; ++r)
      w.write_packed_record(bytes);
    w.finish();
  }
  auto hdr = std::make_unique<edflib_hdr_t>();
  ok(edfopen_file_readonly(packed.string().c_str(), hdr.get(),
                           EDFLIB_READ_ALL_ANNOTATIONS));
  for (unsigned c = 0; c < 2; ++c) {
    auto expected = values(h.signals[c], 0);
    std::vector<int> got(expected.size());
    check(edfread_digital_samples(hdr->handle, c, int(got.size()),
                                  got.data()) == int(got.size()),
          "packed upstream count");
    check(std::equal(got.begin(), got.end(), expected.begin()),
          "packed signed endpoints");
  }
  ok(edfclose_file(hdr->handle));
}

static void annotation_modes(const std::filesystem::path &dir, bool bdf,
                             std::size_t tracks) {
  auto h = info(bdf);
  h.subsecond = Tick{};
  auto path = dir / (bdf ? "modes.bdf" : "modes.edf");
  {
    auto w =
        Writer::create(path, h, WriterOptions{.annotation_signals = tracks});
    w.add_annotation({Tick{0}, std::nullopt, "before"});
    w.add_annotation({Tick{1000000}, Tick{0}, "Recording ends"});
    w.add_annotation({Tick{9000000}, Tick{2500000}, "after"});
    std::vector<std::int32_t> record(10);
    for (unsigned i = 0; i < 4; ++i)
      w.write_digital_record(record);
    w.finish();
  }
  for (auto mode :
       {AnnotationMode::none, AnnotationMode::lazy, AnnotationMode::all,
        AnnotationMode::until_recording_ends}) {
    Reader reader(path, mode);
    int upstream_mode = mode == AnnotationMode::none
                            ? EDFLIB_DO_NOT_READ_ANNOTATIONS
                        : mode == AnnotationMode::until_recording_ends
                            ? EDFLIB_READ_ANNOTATIONS
                            : EDFLIB_READ_ALL_ANNOTATIONS;
    auto hdr = std::make_unique<edflib_hdr_t>();
    ok(edfopen_file_readonly(path.string().c_str(), hdr.get(), upstream_mode));
    auto &annotations = reader.annotations();
    auto count = mode == AnnotationMode::none ? 0
                 : mode == AnnotationMode::until_recording_ends && tracks == 1
                     ? 2
                     : 3;
    check(annotations.size() == static_cast<std::size_t>(count) &&
              hdr->annotations_in_file == count,
          "annotation mode count");
    for (int i = 0; i < count; ++i) {
      edflib_annotation_t a{};
      ok(edf_get_annotation(hdr->handle, i, &a));
      auto &b = annotations[i];
      check(b.text == a.annotation && b.onset.count() == a.onset,
            "annotation mode content");
      check(b.duration ? b.duration->count() == a.duration_l : a.duration_l < 0,
            "annotation mode duration");
    }
    if (count) {
      check(!annotations[0].duration, "absent duration");
      if (count > 1)
        check(annotations[1].duration && annotations[1].duration->count() == 0,
              "explicit zero duration");
      if (count > 2)
        check(annotations[2].duration &&
                  annotations[2].duration->count() == 2500000,
              "positive duration");
    }
    ok(edfclose_file(hdr->handle));
  }
}
static std::string unpad(std::string s) {
  while (!s.empty() && s.back() == ' ')
    s.pop_back();
  return s;
}
static void annotation_edges(const std::filesystem::path &dir, bool bdf) {
  auto h = info(bdf);
  h.patient_additional = "X  patient_under_score   X";
  h.recording_additional = "X  recording_under_score   X";
  for (const auto &annotations : std::vector<std::vector<Annotation>>{
           {},
           {{Tick{0}, std::nullopt, ""}},
           {{Tick{2500000}, Tick{1250000}, "first\tline\nnext\rlast"},
            {Tick{9000000}, std::nullopt, std::string(39, 'a') + "\xc3\xa9"},
            {Tick{12000000}, Tick{0}, std::string("hi\0end", 6)}}}) {
    for (bool upstream : {false, true}) {
      auto path = dir / (std::string(upstream ? "edges-up" : "edges-fast") +
                         (bdf ? ".bdf" : ".edf"));
      if (upstream)
        create_upstream(path, h, bdf, annotations);
      else
        create_fast(path, h, annotations);
      auto expected = annotations;
      if (expected.size() == 3) {
        expected[0].text = "first.line.next.last";
        expected[1].text = std::string(39, 'a') + "\xc3";
        expected[2].text = "hi";
      }
      Reader reader(path);
      auto hdr = std::make_unique<edflib_hdr_t>();
      ok(edfopen_file_readonly(path.string().c_str(), hdr.get(),
                               EDFLIB_READ_ALL_ANNOTATIONS));
      try {
        check(reader.header().patient_additional == h.patient_additional &&
                  unpad(hdr->patient_additional) == h.patient_additional,
              "patient free text preserves spaces, underscores and X");
        check(reader.header().recording_additional == h.recording_additional &&
                  unpad(hdr->recording_additional) == h.recording_additional,
              "recording free text preserves spaces, underscores and X");
        check(reader.annotations().size() == expected.size() &&
                  hdr->annotations_in_file == expected.size(),
              "genuine writer annotation edge count");
        for (std::size_t i = 0; i < expected.size(); ++i) {
          edflib_annotation_t a{};
          ok(edf_get_annotation(hdr->handle, int(i), &a));
          const auto &got = reader.annotations()[i];
          check(got.text == expected[i].text && got.text == a.annotation &&
                    got.onset == expected[i].onset &&
                    got.onset.count() == a.onset &&
                    got.duration == expected[i].duration &&
                    (got.duration ? got.duration->count() == a.duration_l
                                  : a.duration_l < 0),
                "genuine writer annotation edge content");
        }
        if (upstream) {
          auto chained = dir / (bdf ? "chain.bdf" : "chain.edf");
          {
            auto w = Writer::create(chained, reader.header());
            for (const auto &a : reader.annotations())
              w.add_annotation(a);
            for (std::uint64_t r = 0; r < reader.record_count(); ++r) {
              std::vector<std::int32_t> record;
              for (std::size_t c = 0; c < reader.signals().size(); ++c) {
                auto n = reader.signals()[c].samples_per_record;
                std::vector<std::int32_t> samples(n);
                check(reader.read_at(c, r * n, samples) == n,
                      "chain source samples");
                record.insert(record.end(), samples.begin(), samples.end());
              }
              w.write_digital_record(record);
            }
            w.finish();
          }
          Reader copy(chained);
          auto ref = std::make_unique<edflib_hdr_t>();
          ok(edfopen_file_readonly(chained.string().c_str(), ref.get(),
                                   EDFLIB_READ_ALL_ANNOTATIONS));
          try {
            check(copy.header().patient_additional ==
                          reader.header().patient_additional &&
                      unpad(ref->patient_additional) ==
                          reader.header().patient_additional &&
                      copy.header().recording_additional ==
                          reader.header().recording_additional &&
                      unpad(ref->recording_additional) ==
                          reader.header().recording_additional &&
                      copy.record_count() == reader.record_count() &&
                      copy.header().record_duration ==
                          reader.header().record_duration &&
                      copy.header().subsecond == reader.header().subsecond,
                  "chain preserves extracted header");
            for (std::size_t c = 0; c < reader.signals().size(); ++c) {
              std::vector<std::int32_t> source(reader.sample_count(c));
              std::vector<int> target(source.size());
              reader.read_at(c, 0, source);
              check(
                  edfread_digital_samples(ref->handle, int(c),
                                          int(target.size()), target.data()) ==
                          int(target.size()) &&
                      std::equal(source.begin(), source.end(), target.begin()),
                  "chain preserves extracted samples");
            }
            check(ref->annotations_in_file == reader.annotations().size(),
                  "chain preserves extracted annotation count");
            for (std::size_t i = 0; i < reader.annotations().size(); ++i) {
              edflib_annotation_t a{};
              ok(edf_get_annotation(ref->handle, int(i), &a));
              const auto &original = reader.annotations()[i];
              check(original.text == a.annotation &&
                        original.onset.count() == a.onset &&
                        (original.duration
                             ? original.duration->count() == a.duration_l
                             : a.duration_l < 0),
                    "chain preserves extracted annotations");
            }
            ok(edfclose_file(ref->handle));
          } catch (...) {
            edfclose_file(ref->handle);
            throw;
          }
        }
        ok(edfclose_file(hdr->handle));
      } catch (...) {
        edfclose_file(hdr->handle);
        throw;
      }
    }
  }
  // Keep a genuine upstream header/data layout and replace only the first TAL
  // payload. EDFlib accepts further fraction digits and truncates to ticks.
  auto fractional = dir / (bdf ? "fractional.bdf" : "fractional.edf");
  create_upstream(fractional, h, bdf, {});
  {
    std::fstream file(fractional,
                      std::ios::in | std::ios::out | std::ios::binary);
    const auto width = bdf ? 3u : 2u;
    const auto record_bytes =
        (std::filesystem::file_size(fractional) - 1024) / 4;
    std::string payload(record_bytes - 10 * width, '\0');
    std::string tal = "+0.1234500\x14\x14";
    tal.push_back('\0');
    tal += "+0.3734500999\x15"
           "0.1250000999\x14"
           "fraction\x14";
    tal.push_back('\0');
    tal += "+0.5\x15"
           "123456789.1234567\x14"
           "long duration\x14";
    tal.push_back('\0');
    check(tal.size() <= payload.size(), "fraction fixture capacity");
    payload.replace(0, tal.size(), tal);
    file.seekp(1024 + 10 * width);
    file.write(payload.data(), payload.size());
    check(bool(file), "fraction fixture patch");
  }
  {
    Reader reader(fractional);
    auto hdr = std::make_unique<edflib_hdr_t>();
    ok(edfopen_file_readonly(fractional.string().c_str(), hdr.get(),
                             EDFLIB_READ_ALL_ANNOTATIONS));
    edflib_annotation_t a{};
    ok(edf_get_annotation(hdr->handle, 0, &a));
    edflib_annotation_t long_duration{};
    ok(edf_get_annotation(hdr->handle, 1, &long_duration));
    ok(edfclose_file(hdr->handle));
    check(reader.annotations().size() == 2 && hdr->annotations_in_file == 2,
          "long fractional TAL count");
    const auto &got = reader.annotations()[0];
    check(got.text == "fraction" && got.text == a.annotation &&
              got.onset.count() == 2500000 && got.onset.count() == a.onset &&
              got.duration && got.duration->count() == 1250000 &&
              got.duration->count() == a.duration_l,
          "long fractional TAL truncation differential");
    check(reader.annotations()[1].duration &&
              reader.annotations()[1].duration->count() == 1234567891234500LL &&
              reader.annotations()[1].duration->count() ==
                  long_duration.duration_l,
          "long TAL duration matches upstream read buffer precision");
  }
}

static void verify_external(const std::filesystem::path &path) {
  Reader reader(path, AnnotationMode::all);
  auto hdr = std::make_unique<edflib_hdr_t>();
  ok(edfopen_file_readonly(path.string().c_str(), hdr.get(),
                           EDFLIB_READ_ALL_ANNOTATIONS));
  auto handle = hdr->handle;
  try {
    check(reader.signals().size() == static_cast<std::size_t>(hdr->edfsignals),
          "external signal count");
    check(reader.record_count() ==
              static_cast<std::uint64_t>(hdr->datarecords_in_file),
          "external record count");
    check(reader.header().record_duration.count() == hdr->datarecord_duration,
          "external duration");
    check(reader.header().subsecond.count() == hdr->starttime_subsecond,
          "external subsecond");
    std::uint64_t total = 0;
    std::vector<std::int32_t> fast(65536);
    std::vector<int> up(fast.size());
    std::vector<double> pf(fast.size()), pu(fast.size());
    for (std::size_t c = 0; c < reader.signals().size(); ++c) {
      auto n = reader.sample_count(c);
      check(n == static_cast<std::uint64_t>(hdr->signalparam[c].smp_in_file),
            "external sample count");
      total += n;
      for (std::uint64_t offset = 0; offset < n;) {
        auto count = static_cast<std::size_t>(
            std::min<std::uint64_t>(fast.size(), n - offset));
        check(reader.read_at(c, offset, std::span(fast).first(count)) ==
                  count,
              "external digital count");
        check(edfread_digital_samples(handle, c, int(count), up.data()) ==
                  int(count),
              "external upstream count");
        check(std::equal(fast.begin(), fast.begin() + count, up.begin()),
              "external digital mismatch");
        offset += count;
      }
      check(edfseek(handle, c, 0, EDFSEEK_SET) == 0, "external seek");
      for (std::uint64_t offset = 0; offset < n;) {
        auto count = static_cast<std::size_t>(
            std::min<std::uint64_t>(pf.size(), n - offset));
        check(reader.read_at(c, offset, std::span(pf).first(count)) ==
                  count,
              "external physical count");
        check(edfread_physical_samples(handle, c, int(count), pu.data()) ==
                  int(count),
              "external upstream physical count");
        for (std::size_t k = 0; k < count; ++k)
          check(std::abs(pf[k] - pu[k]) <=
                    1e-10 * std::max(1.0, std::abs(pu[k])),
                "external physical mismatch");
        offset += count;
      }
      std::cout << "signal " << c << ": " << n << " samples, "
                << reader.sample_rate(c) << " Hz\n";
    }
    auto &annotations = reader.annotations();
    check(annotations.size() ==
              static_cast<std::size_t>(hdr->annotations_in_file),
          "external annotation count");
    for (std::size_t i = 0; i < annotations.size(); ++i) {
      edflib_annotation_t a{};
      ok(edf_get_annotation(handle, int(i), &a));
      auto &b = annotations[i];
      check(b.text == a.annotation && b.onset.count() == a.onset,
            "external annotation mismatch");
      check(b.duration ? b.duration->count() == a.duration_l : a.duration_l < 0,
            "external annotation duration");
    }
    std::cout << "External differential passed: bytes="
              << std::filesystem::file_size(path)
              << ", records=" << reader.record_count()
              << ", record_duration_ticks="
              << reader.header().record_duration.count()
              << ", total_samples=" << total
              << ", annotations=" << annotations.size() << "\n";
    ok(edfclose_file(handle));
  } catch (...) {
    edfclose_file(handle);
    throw;
  }
}

static void verify_writer_formats(const std::filesystem::path &dir, bool bdf) {
  auto h = info(bdf);
  h.format = bdf ? Format::bdf : Format::edf;
  h.subsecond = Tick{};
  auto classic = dir / (bdf ? "classic.bdf" : "classic.edf");
  create_fast(classic, h, {});
  verify_external(classic);

  h = info(bdf);
  h.signals.clear();
  h.record_duration = Tick{};
  auto events = dir / (bdf ? "events.bdf" : "events.edf");
  WriterOptions options;
  options.annotation_bytes_per_signal = 600;
  auto w = Writer::create(events, h, options);
  w.add_annotation({Tick{}, std::nullopt, ""});
  w.add_annotation({Tick{1}, Tick{}, std::string(512, 'a')});
  w.finish();
  verify_external(events);
  // A third-party file can have descriptions larger than either writer emits.
  // The reader must consume the entire TAL but return only its first 512 bytes.
  {
    std::fstream file(events, std::ios::binary | std::ios::in | std::ios::out);
    std::string bytes((std::istreambuf_iterator<char>(file)), {});
    auto offset = bytes.find(std::string(40, 'a'));
    check(offset != bytes.npos, "find annotation bytes");
    auto tail = std::string(512, 'x') + "extra\x14";
    tail.push_back('\0');
    tail += "+0.5\x14next\x14";
    check(offset + tail.size() < bytes.size(), "annotation fixture capacity");
    file.seekp(static_cast<std::streamoff>(offset));
    file.write(tail.data(), tail.size());
    check(bool(file), "write long annotation fixture");
  }
  verify_external(events);
  Reader long_reader(events);
  const auto &annotations = long_reader.annotations();
  check(annotations.size() == 3 &&
            annotations[1].text == std::string(512, 'x') &&
            annotations[2].text == "next",
        "long description and following TAL");
}

void write_code_tests(const std::filesystem::path &dir);

int main(int argc, char **argv) {
  auto dir = std::filesystem::temp_directory_path() /
             ("fastbdf-interop-" +
              std::to_string(
                  std::chrono::steady_clock::now().time_since_epoch().count()));
  try {
    check(edflib_version() == 127, "requires upstream EDFlib 1.27");
    if (argc > 1) {
      for (int i = 1; i < argc; ++i)
        verify_external(argv[i]);
      return 0;
    }
    std::filesystem::create_directory(dir);
    write_code_tests(dir);
    for (bool bdf : {false, true}) {
      verify_writer_formats(dir, bdf);
      annotation_edges(dir, bdf);
      annotation_modes(dir, bdf, 1);
      annotation_modes(dir, bdf, 2);
      conversion_tests(dir, bdf);
      auto h = info(bdf);
      auto f = dir / (bdf ? "fast.bdf" : "fast.edf"),
           u = dir / (bdf ? "upstream.bdf" : "upstream.edf");
      create_fast(f, h);
      verify(f, h);
      create_upstream(u, h, bdf);
      verify(u, h);
    }
    std::filesystem::remove_all(dir);
    std::cout
        << "EDFlib 1.27 bidirectional EDF+/BDF+ interoperability passed\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << "; fixtures: " << dir << '\n';
    return 1;
  }
}
