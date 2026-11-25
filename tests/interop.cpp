#include "edflib.h"
#include <array>
#include <cmath>
#include <fastbdf/fastbdf.hpp>
#include <filesystem>
#include <iostream>
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
static void create_fast(const std::filesystem::path &p,
                        const RecordingInfo &h) {
  auto w = Writer::create(p, h);
  w.add_annotation({Tick{2500000}, Tick{1250000}, "event \xCE\xB1"});
  w.add_annotation({Tick{9000000}, std::nullopt, "end"});
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
static void create_upstream(const std::filesystem::path &p,
                            const RecordingInfo &h, bool bdf) {
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
    ok(edfwrite_annotation_utf8_hr(f, 2500000, 1250000, "event \xCE\xB1"));
    ok(edfwrite_annotation_utf8_hr(f, 9000000, -1, "end"));
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
      check(reader.read_digital_at(id, 0, a) == a.size(), "fast read count");
      check(edfread_digital_samples(f, id, int(b.size()), b.data()) ==
                int(b.size()),
            "upstream read count");
      check(a == expected && std::equal(a.begin(), a.end(), b.begin()),
            "digital differential");
      std::vector<double> pa(a.size()), pb(a.size());
      reader.read_physical_at(id, 0, pa);
      check(edfseek(f, id, 0, EDFSEEK_SET) == 0, "seek");
      check(edfread_physical_samples(f, id, int(pb.size()), pb.data()) ==
                int(pb.size()),
            "physical count");
      for (unsigned k = 0; k < pa.size(); ++k)
        check(std::abs(pa[k] - pb[k]) < 1e-10, "physical differential");
      std::array<std::int32_t, 5> slice{};
      reader.read_digital_at(id, s.samples_per_record - 1, slice);
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
    a.read_digital_at(c, 0, x);
    b.read_digital_at(c, 0, y);
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
        check(reader.read_digital_at(c, offset, std::span(fast).first(count)) ==
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
        check(reader.read_physical_at(c, offset, std::span(pf).first(count)) ==
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
    for (bool bdf : {false, true}) {
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
