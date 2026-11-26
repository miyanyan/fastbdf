#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>

#include <psapi.h>
#endif
#include "edflib.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fastbdf/fastbdf.hpp>
#include <filesystem>
#include <iostream>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <vector>
using Clock = std::chrono::steady_clock;
static void check(bool x) {
  if (!x)
    throw std::runtime_error("benchmark I/O or checksum failed");
}
template <class F> double median(F f) {
  f();
  std::array<double, 5> t{};
  for (auto &x : t) {
    auto start = Clock::now();
    f();
    x = std::chrono::duration<double>(Clock::now() - start).count();
  }
  std::sort(t.begin(), t.end());
  return t[2];
}
static double quantile(std::vector<double> values, double q) {
  std::sort(values.begin(), values.end());
  return values[static_cast<std::size_t>(std::ceil(q * values.size())) - 1];
}
static void writer_benchmark(const std::filesystem::path &path,
                             const fastbdf::RecordingInfo &info,
                             const std::vector<std::int32_t> &record,
                             unsigned records) {
  const unsigned channels = static_cast<unsigned>(info.signals.size()),
                 spr = info.signals[0].samples_per_record;
  auto verify = [&] {
    auto header = std::make_unique<edflib_hdr_t>();
    check(edfopen_file_readonly(path.string().c_str(), header.get(),
                                EDFLIB_DO_NOT_READ_ANNOTATIONS) == 0);
    check(header->datarecords_in_file == records &&
          header->edfsignals == int(channels));
    std::vector<int> out(spr);
    for (unsigned c = 0; c < channels; ++c)
      for (unsigned r = 0; r < records; ++r) {
        check(edfread_digital_samples(header->handle, c, spr, out.data()) ==
              int(spr));
        check(std::equal(out.begin(), out.end(), record.begin() + c * spr));
      }
    check(edfclose_file(header->handle) == 0);
  };
  auto run = [&](bool fast) {
    std::vector<double> times;
    for (unsigned i = 0; i < 6; ++i) {
      auto begin = Clock::now();
      if (fast) {
        auto writer = fastbdf::Writer::create(path, info);
        for (unsigned r = 0; r < records; ++r)
          writer.write_digital_record(record);
        writer.finish();
      } else {
        auto handle = edfopen_file_writeonly(path.string().c_str(),
                                             EDFLIB_FILETYPE_BDFPLUS, channels);
        check(handle >= 0);
        for (unsigned c = 0; c < channels; ++c) {
          auto &s = info.signals[c];
          check(edf_set_samplefrequency(handle, c, s.samples_per_record) == 0);
          check(edf_set_digital_minimum(handle, c, s.digital_min) == 0);
          check(edf_set_digital_maximum(handle, c, s.digital_max) == 0);
          check(edf_set_physical_minimum(handle, c, s.physical_min) == 0);
          check(edf_set_physical_maximum(handle, c, s.physical_max) == 0);
          check(edf_set_label(handle, c, s.label.c_str()) == 0);
        }
        for (unsigned r = 0; r < records; ++r)
          check(edf_blockwrite_digital_samples(
                    handle, const_cast<std::int32_t *>(record.data())) == 0);
        check(edfclose_file(handle) == 0);
      }
      auto elapsed =
          std::chrono::duration<double>(Clock::now() - begin).count();
      verify();
      if (i)
        times.push_back(elapsed);
    }
    return quantile(times, .5);
  };
  double fast = run(true), up = run(false),
         mib = double(records) * channels * spr * 3 / 1048576;
  std::cout << "writer create/configure/records/finalize median5: fastbdf "
            << fast * 1000 << " ms (" << mib / fast << " MiB/s), upstream "
            << up * 1000 << " ms (" << mib / up << " MiB/s), ratio "
            << up / fast << "; every output fully verified outside timing\n";
}
int main(int argc, char **argv) {
  try {
    const unsigned records = argc > 1 ? std::stoul(argv[1]) : 1000;
    check(records >= 4);
    constexpr unsigned channels = 16, spr = 1024;
    auto path =
        std::filesystem::temp_directory_path() /
        ("fastbdf-bench-" +
         std::to_string(Clock::now().time_since_epoch().count()) + ".bdf");
    fastbdf::RecordingInfo h;
    h.format = fastbdf::Format::bdf_plus;
    for (unsigned c = 0; c < channels; ++c) {
      fastbdf::SignalSpec s;
      s.label = "ch" + std::to_string(c);
      s.samples_per_record = spr;
      s.digital_min = -8388608;
      s.digital_max = 8388607;
      h.signals.push_back(s);
    }
    std::vector<std::int32_t> record(channels * spr);
    for (unsigned i = 0; i < record.size(); ++i)
      record[i] = int(i) - 8192;
    writer_benchmark(path, h, record, records);
    auto start = Clock::now();
    {
      auto w = fastbdf::Writer::create(path, h);
      for (unsigned r = 0; r < records; ++r)
        w.write_digital_record(record);
      w.finish();
    }
    double write = std::chrono::duration<double>(Clock::now() - start).count();
    fastbdf::Reader reader(path, fastbdf::AnnotationMode::none);
    auto hdr = std::make_unique<edflib_hdr_t>();
    check(edfopen_file_readonly(path.string().c_str(), hdr.get(),
                                EDFLIB_DO_NOT_READ_ANNOTATIONS) == 0);
    auto handle = hdr->handle;
    const auto total = static_cast<std::size_t>(records) * spr;
    std::vector<std::int32_t> data(total);
    std::int64_t fastsum{}, upsum{};
    auto sum = [](const auto &v) {
      return std::accumulate(v.begin(), v.end(), std::int64_t{});
    };
    auto report = [&](const char *name, double fast, double up, double bytes) {
      check(fastsum == upsum);
      std::cout << name << ": fastbdf " << fast * 1000 << " ms, upstream "
                << up * 1000 << " ms; " << bytes / 1048576 / fast << " vs "
                << bytes / 1048576 / up << " MiB/s; ratio " << up / fast
                << "; checksum " << fastsum << '\n';
    };
    auto fast_scan = [&](unsigned n) {
      fastsum = 0;
      for (unsigned c = 0; c < n; ++c) {
        check(reader.read_at(c, 0, data) == data.size());
        fastsum += sum(data);
      }
    };
    auto upstream_scan = [&](unsigned n) {
      upsum = 0;
      for (unsigned c = 0; c < n; ++c) {
        check(edfseek(handle, c, 0, EDFSEEK_SET) == 0);
        check(edfread_digital_samples(handle, c, int(data.size()),
                                      data.data()) == int(data.size()));
        upsum += sum(data);
      }
    };
    const double payload = double(records) * channels * spr * 3;
    std::cout << "BDF+ " << channels << " x " << spr << ", " << records
              << " records; writer " << payload / 1048576 / write
              << " MiB/s\nWarm cache, median of 5 after warmup; readers opened "
                 "before timing, checksum included\n";
    double f = median([&] { fast_scan(channels); }),
           u = median([&] { upstream_scan(channels); });
    report("all channels sequential", f, u, payload);
    f = median([&] { fast_scan(1); });
    u = median([&] { upstream_scan(1); });
    report("single channel", f, u, payload / channels);
    std::array<std::size_t, channels> ids{};
    std::iota(ids.begin(), ids.end(), 0);
    fastbdf::TimeRange range{fastbdf::Tick{},
                             fastbdf::Tick{std::int64_t(records) * 10000000}};
    f = median([&] {
      auto all = reader.read_digital(ids, range);
      fastsum = 0;
      for (auto &v : all)
        fastsum += sum(v);
    });
    u = median([&] { upstream_scan(channels); });
    report("bulk full time window (result allocation included)", f, u, payload);
    std::array<std::int32_t, 64> shortdata{};
    constexpr unsigned windows = 512;
    auto offset = [&](unsigned i) {
      return std::uint64_t(i * 2654435761u) % (total - shortdata.size());
    };
    f = median([&] {
      fastsum = 0;
      for (unsigned i = 0; i < windows; ++i) {
        check(reader.read_at(i % channels, offset(i), shortdata) ==
              shortdata.size());
        fastsum += sum(shortdata);
      }
    });
    u = median([&] {
      upsum = 0;
      for (unsigned i = 0; i < windows; ++i) {
        check(edfseek(handle, i % channels, offset(i), EDFSEEK_SET) ==
              static_cast<long long>(offset(i)));
        check(
            edfread_digital_samples(handle, i % channels, int(shortdata.size()),
                                    shortdata.data()) == int(shortdata.size()));
        upsum += sum(shortdata);
      }
    });
    report("random 64-sample windows", f, u, windows * shortdata.size() * 3);
    std::vector<double> fast_latency, up_latency;
    std::array<std::int32_t, 64> other{};
    for (unsigned i = 0; i < windows; ++i) {
      auto begin = Clock::now();
      check(reader.read_at(i % channels, offset(i), shortdata) ==
            shortdata.size());
      fast_latency.push_back(
          std::chrono::duration<double>(Clock::now() - begin).count());
      begin = Clock::now();
      check(edfseek(handle, i % channels, offset(i), EDFSEEK_SET) ==
            static_cast<long long>(offset(i)));
      check(edfread_digital_samples(handle, i % channels, int(other.size()),
                                    other.data()) == int(other.size()));
      up_latency.push_back(
          std::chrono::duration<double>(Clock::now() - begin).count());
      check(shortdata == other);
    }
    std::cout << "random 64-sample request P95 (512 paired requests, verify "
                 "outside timing): fastbdf "
              << quantile(fast_latency, .95) * 1e6 << " us, upstream "
              << quantile(up_latency, .95) * 1e6 << " us\n";
    check(edfclose_file(handle) == 0);
    fast_latency.clear();
    up_latency.clear();
    const auto screen_samples = std::min(total, std::size_t(spr) * 10);
    // Each request owns a new reader and output allocation, reading 10 seconds
    // of every channel.
    for (unsigned request = 0; request < 41; ++request) {
      auto begin = Clock::now();
      std::vector<std::int32_t> screen_fast(channels * screen_samples);
      {
        fastbdf::Reader first(path);
        for (unsigned c = 0; c < channels; ++c)
          check(first.read_at(
                    c, 0,
                    std::span(screen_fast)
                        .subspan(c * screen_samples, screen_samples)) ==
                screen_samples);
      }
      auto fast_elapsed =
          std::chrono::duration<double>(Clock::now() - begin).count();
      begin = Clock::now();
      std::vector<int> screen_up(channels * screen_samples);
      auto header = std::make_unique<edflib_hdr_t>();
      check(edfopen_file_readonly(path.string().c_str(), header.get(),
                                  EDFLIB_READ_ALL_ANNOTATIONS) == 0);
      for (unsigned c = 0; c < channels; ++c)
        check(edfread_digital_samples(header->handle, c, int(screen_samples),
                                      screen_up.data() + c * screen_samples) ==
              int(screen_samples));
      check(edfclose_file(header->handle) == 0);
      auto up_elapsed =
          std::chrono::duration<double>(Clock::now() - begin).count();
      check(std::equal(screen_fast.begin(), screen_fast.end(),
                       screen_up.begin()));
      if (request) {
        fast_latency.push_back(fast_elapsed);
        up_latency.push_back(up_elapsed);
      }
    }
    std::cout << "open-to-first-screen P95 (40 warm-cache paired opens, "
                 "10s/all channels, allocation/read/close included): fastbdf "
              << quantile(fast_latency, .95) * 1000 << " ms, upstream "
              << quantile(up_latency, .95) * 1000 << " ms\n";
#if defined(_WIN32)
    PROCESS_MEMORY_COUNTERS memory{};
    memory.cb = sizeof(memory);
    if (GetProcessMemoryInfo(GetCurrentProcess(), &memory, sizeof(memory)))
      std::cout << "Whole benchmark process peak working set: "
                << double(memory.PeakWorkingSetSize) / 1048576
                << " MiB (shared process; not per-implementation memory "
                   "comparison)\n";
#endif
    std::filesystem::remove(path);
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
