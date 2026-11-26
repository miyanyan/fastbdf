#include <fastbdf/reader.hpp>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/filesystem.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace fastbdf;

namespace {
template <class T>
using Array =
    nb::ndarray<nb::numpy, T, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

class PythonReader {
  std::shared_ptr<const Reader> reader_;

public:
  PythonReader(const std::filesystem::path &path, AnnotationMode mode) {
    nb::gil_scoped_release release;
    reader_ = std::make_shared<Reader>(path, mode);
  }
  std::shared_ptr<const Reader> get() const {
    if (!reader_)
      throw nb::value_error("Reader is closed");
    return reader_;
  }
  void close() { reader_.reset(); }

  template <class T>
  size_t read_into(size_t channel, uint64_t start, Array<T> out) const {
    // Keep the file alive if another Python thread calls close() during I/O.
    auto reader = get();
    if (reinterpret_cast<uintptr_t>(out.data()) % alignof(T))
      throw nb::value_error("Output array must be aligned");
    auto buffer = std::span<T>(out.data(), out.size());
    nb::gil_scoped_release release;
    return reader->read_at(channel, start, buffer);
  }

  void window_into(
      const std::vector<size_t> &channels, int64_t start_tick,
      int64_t stop_tick,
      nb::ndarray<nb::numpy, double, nb::ndim<2>, nb::c_contig, nb::device::cpu>
          out) const {
    auto reader = get();
    if (channels.size() != out.shape(0))
      throw nb::value_error("Expected one output row per channel");
    if (reinterpret_cast<uintptr_t>(out.data()) % alignof(double))
      throw nb::value_error("Output array must be aligned");
    std::vector<std::span<double>> rows;
    rows.reserve(channels.size());
    for (size_t i = 0; i < channels.size(); ++i)
      rows.emplace_back(out.shape(1) ? out.data() + i * out.shape(1)
                                     : out.data(),
                        out.shape(1));
    nb::gil_scoped_release release;
    reader->read_physical_into(channels, {Tick{start_tick}, Tick{stop_tick}},
                               rows);
  }

  nb::list window(const std::vector<size_t> &channels, int64_t start_tick,
                  int64_t stop_tick) const {
    auto reader = get();
    std::vector<std::vector<double>> values;
    {
      nb::gil_scoped_release release;
      reader->read_physical_into(channels, {Tick{start_tick}, Tick{stop_tick}},
                                 values);
    }
    nb::list result;
    for (auto &channel : values) {
      auto storage = std::make_unique<std::vector<double>>(std::move(channel));
      if (storage->empty())
        storage->reserve(1); // provide a valid array address for empty windows
      nb::capsule owner(storage.get(), [](void *p) noexcept {
        delete static_cast<std::vector<double> *>(p);
      });
      auto *v = storage.release();
      result.append(Array<double>(v->data(), {v->size()}, owner).cast());
    }
    return result;
  }
};
} // namespace

NB_MODULE(fastbdf, m) {
  nb::enum_<AnnotationMode>(m, "AnnotationMode")
      .value("none", AnnotationMode::none)
      .value("lazy", AnnotationMode::lazy)
      .value("all", AnnotationMode::all)
      .value("until_recording_ends", AnnotationMode::until_recording_ends);
  nb::class_<SignalSpec>(m, "Signal")
      .def_ro("label", &SignalSpec::label)
      .def_ro("unit", &SignalSpec::unit)
      .def_ro("samples_per_record", &SignalSpec::samples_per_record)
      .def_ro("physical_min", &SignalSpec::physical_min)
      .def_ro("physical_max", &SignalSpec::physical_max)
      .def_ro("digital_min", &SignalSpec::digital_min)
      .def_ro("digital_max", &SignalSpec::digital_max);
  nb::class_<PythonReader>(m, "Reader")
      .def(nb::init<const std::filesystem::path &, AnnotationMode>(), "path"_a,
           "annotation_mode"_a = AnnotationMode::all)
      .def("close", &PythonReader::close)
      .def(
          "__enter__",
          [](PythonReader &self) -> PythonReader & {
            self.get();
            return self;
          },
          nb::rv_policy::reference)
      .def("__exit__", [](PythonReader &self, nb::args) { self.close(); })
      .def_prop_ro(
          "signals",
          [](const PythonReader &self) { return self.get()->signals(); })
      .def_prop_ro(
          "annotations",
          [](const PythonReader &self) {
            auto reader = self.get();
            nb::list events;
            for (const auto &event : reader->annotations()) {
              auto duration = event.duration ? nb::cast(event.duration->count())
                                             : nb::none();
              events.append(nb::make_tuple(
                  event.onset.count(), duration,
                  nb::bytes(event.text.data(), event.text.size())));
            }
            return events;
          },
          "List of (onset_ticks, duration_ticks_or_None, text_bytes).")
      .def(
          "sample_rate",
          [](const PythonReader &self, size_t channel) {
            return self.get()->sample_rate(channel);
          },
          "channel"_a)
      .def(
          "sample_count",
          [](const PythonReader &self, size_t channel) {
            return self.get()->sample_count(channel);
          },
          "channel"_a)
      .def("read_physical_into", &PythonReader::read_into<double>, "channel"_a,
           "start"_a, "out"_a.noconvert())
      .def("read_digital_into", &PythonReader::read_into<std::int32_t>,
           "channel"_a, "start"_a, "out"_a.noconvert())
      .def("read_physical_window_into", &PythonReader::window_into,
           "channels"_a, "start_tick"_a, "stop_tick"_a, "out"_a.noconvert(),
           "Fill a writable C-contiguous float64 matrix. Each selected channel "
           "must have exactly out.shape[1] native samples in the window.")
      .def("read_physical_window", &PythonReader::window, "channels"_a,
           "start_tick"_a, "stop_tick"_a,
           "Return one NumPy array per channel, in file units, for "
           "[start_tick, stop_tick). A tick is 100 ns.");
}
