#include "edflib.h"
#include <algorithm>
#include <cmath>
#include <fastbdf/fastbdf.hpp>
#include <filesystem>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <vector>

using namespace fastbdf;

namespace {
void require(bool condition, const char *what) {
  if (!condition)
    throw std::runtime_error(what);
}

using Codes = std::vector<int>;

// A record contains consecutive channel blocks. Read each file through both
// libraries and compare to explicit codes, including the returned sample count.
void check_codes(const std::filesystem::path &path, const RecordingInfo &info,
                 const std::vector<Codes> &records) {
  Reader reader(path);
  auto hdr = std::make_unique<edflib_hdr_t>();
  require(edfopen_file_readonly(path.string().c_str(), hdr.get(),
                                EDFLIB_READ_ALL_ANNOTATIONS) == 0,
          "open result");
  try {
    require(reader.record_count() == records.size() &&
                hdr->datarecords_in_file ==
                    static_cast<long long>(records.size()),
            "record count");
    for (size_t r = 0; r < records.size(); ++r) {
      size_t offset = 0;
      for (size_t c = 0; c < info.signals.size(); ++c) {
        int n = static_cast<int>(info.signals[c].samples_per_record);
        std::vector<std::int32_t> fast(n);
        Codes ref(n);
        require(reader.read(c, fast) == fast.size() &&
                    edfread_digital_samples(hdr->handle, static_cast<int>(c), n,
                                            ref.data()) == n,
                "read count");
        for (int i = 0; i < n; ++i) {
          int expected = records[r].at(offset++);
          if (fast[i] != expected || ref[i] != expected) {
            std::ostringstream msg;
            msg << path.filename() << " record=" << r << " channel=" << c
                << " sample=" << i << " expected=" << expected
                << " fastbdf reader=" << fast[i] << " EDFlib reader=" << ref[i];
            throw std::runtime_error(msg.str());
          }
        }
      }
      require(offset == records[r].size(), "fixture record size");
    }
  } catch (...) {
    edfclose_file(hdr->handle);
    throw;
  }
  require(edfclose_file(hdr->handle) == 0, "close result");
}

template <class Write>
void compare_writers(const std::filesystem::path &dir, const char *name,
                     const RecordingInfo &info,
                     const std::vector<Codes> &expected, Write write,
                     const std::vector<Codes> &reference_expected = {}) {
  bool bdf = info.format == Format::bdf_plus;
  auto prefix = std::string(name);
  auto fast = dir / (prefix + "-fast" + (bdf ? ".bdf" : ".edf"));
  auto ref = dir / (prefix + "-ref" + (bdf ? ".bdf" : ".edf"));
  auto writer = Writer::create(fast, info);
  int h = edfopen_file_writeonly(ref.string().c_str(),
                                 bdf ? EDFLIB_FILETYPE_BDFPLUS
                                     : EDFLIB_FILETYPE_EDFPLUS,
                                 static_cast<int>(info.signals.size()));
  require(h >= 0, "create reference");
  try {
    for (size_t c = 0; c < info.signals.size(); ++c) {
      const auto &s = info.signals[c];
      int id = static_cast<int>(c);
      require(edf_set_samplefrequency(h, id, s.samples_per_record) == 0 &&
                  edf_set_digital_minimum(h, id, s.digital_min) == 0 &&
                  edf_set_digital_maximum(h, id, s.digital_max) == 0 &&
                  edf_set_physical_minimum(h, id, s.physical_min) == 0 &&
                  edf_set_physical_maximum(h, id, s.physical_max) == 0,
              "configure reference");
    }
    write(writer, h);
    writer.finish();
  } catch (...) {
    edfclose_file(h);
    throw;
  }
  require(edfclose_file(h) == 0, "finish reference");
  check_codes(ref, info,
              reference_expected.empty() ? expected : reference_expected);
  check_codes(fast, info, expected);
}

void short_edf(const std::filesystem::path &dir) {
  RecordingInfo info;
  info.add_signal("full", 5);
  auto &s = info.add_signal("limited", 4);
  s.digital_min = -10;
  s.digital_max = 10;
  std::vector<short> first{-32768, -32767, 0, 32766, 32767, -11, -10, 10, 11};
  std::vector<short> second{32767, 1, -1, 0, -32768, 9, -9, 32767, -32768};
  compare_writers(
      dir, "short", info,
      {{-32768, -32767, 0, 32766, 32767, -10, -10, 10, 10},
       {32767, 1, -1, 0, -32768, 9, -9, 10, -10}},
      [&](Writer &w, int h) {
        w.write_samples(std::span(first).first(5));
        w.write_samples(std::span(first).subspan(5));
        require(edfwrite_digital_short_samples(h, first.data()) == 0 &&
                    edfwrite_digital_short_samples(h, first.data() + 5) == 0,
                "write short channels");
        w.write_short_record(second);
        require(edf_blockwrite_digital_short_samples(h, second.data()) == 0,
                "write short block");
      });
}

void digital_bdf(const std::filesystem::path &dir) {
  RecordingInfo info;
  info.format = Format::bdf_plus;
  info.add_signal("A", 7);
  info.add_signal("B", 4);
  constexpr int lo = -8388608, hi = 8388607;
  Codes block{lo - 1, lo, lo + 1, -1,       0,        1,
              hi - 1, hi, hi + 1, 0x123456, -0x123456};
  Codes clipped = block;
  clipped.front() = lo;
  clipped[8] = hi;
  Codes packed_codes{hi, hi - 1, 1,     0,   -1,  lo + 1,
                     lo, -65537, 65536, 255, -256};
  std::vector<std::byte> packed;
  for (int code : packed_codes)
    for (unsigned shift : {0u, 8u, 16u})
      packed.push_back(static_cast<std::byte>(
          (static_cast<std::uint32_t>(code) >> shift) & 255));
  compare_writers(
      dir, "digital24", info, {clipped, packed_codes}, [&](Writer &w, int h) {
        w.write_digital_record(block);
        require(edf_blockwrite_digital_samples(h, block.data()) == 0,
                "write int block");
        w.write_packed_record(packed);
        require(edf_blockwrite_digital_3byte_samples(h, packed.data()) == 0,
                "write packed block");
      });
}

void negative_gain(const std::filesystem::path &dir, Format fmt) {
  RecordingInfo info;
  info.format = fmt;
  auto &exact = info.add_signal("exact", 11);
  exact.physical_min = 4;
  exact.physical_max = -4;
  exact.digital_min = -2;
  exact.digital_max = 2;
  auto &asym = info.add_signal("asymmetric", 5);
  asym.physical_min = 19;
  asym.physical_max = -7;
  asym.digital_min = -100;
  asym.digital_max = 200;
  std::vector<double> samples{5,
                              4,
                              std::nextafter(2., 4.),
                              2.,
                              std::nextafter(2., 0.),
                              0.,
                              std::nextafter(-2., 0.),
                              -2.,
                              std::nextafter(-2., -4.),
                              -4,
                              -5,
                              20,
                              19,
                              6,
                              -7,
                              -8};
  Codes expected{-2, -2, -1, -1, 0, 0, 0, 1, 1, 2, 2, -100, -100, 50, 200, 200};
  Codes reference = expected;
  // At physical_min = 19, EDFlib truncates -99.99999999999999 to -99.
  reference[12] = -99;
  compare_writers(
      dir, "negative", info, {expected, expected},
      [&](Writer &w, int h) {
        w.write_samples(std::span(samples).first(11));
        w.write_samples(std::span(samples).subspan(11));
        require(edfwrite_physical_samples(h, samples.data()) == 0 &&
                    edfwrite_physical_samples(h, samples.data() + 11) == 0,
                "write physical channels");
        w.write_physical_record(samples);
        require(edf_blockwrite_physical_samples(h, samples.data()) == 0,
                "write physical block");
      },
      {reference, reference});
}
} // namespace

void write_code_tests(const std::filesystem::path &dir) {
  short_edf(dir);
  digital_bdf(dir);
  negative_gain(dir, Format::edf_plus);
  negative_gain(dir, Format::bdf_plus);
}
