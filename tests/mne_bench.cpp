// Optional helper for mne_bench.py; no Python or EDFlib build dependency.
#include <algorithm>
#include <chrono>
#include <fastbdf/fastbdf.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>

int main(int argc, char** argv)
{
    try {
        if (argc != 5) throw std::runtime_error("usage: fastbdf_mne_bench edf|bdf fixture output records");
        const bool bdf = std::string(argv[1]) == "bdf";
        if (!bdf && std::string(argv[1]) != "edf") throw std::runtime_error("format must be edf or bdf");
        constexpr unsigned channels = 16;
        constexpr unsigned rate = 512;
        constexpr unsigned repeats = 5;
        const auto record_text = std::string(argv[4]);
        std::size_t consumed = 0;
        const auto parsed_records = std::stoul(record_text, &consumed);
        if (consumed != record_text.size() || parsed_records < 1 || parsed_records > 4096)
            throw std::runtime_error("records must be in [1, 4096]");
        const auto records = static_cast<unsigned>(parsed_records);
        fastbdf::RecordingInfo info;
        info.format = bdf ? fastbdf::Format::bdf_plus : fastbdf::Format::edf_plus;
        info.start = std::chrono::local_days{std::chrono::year{2024} / 1 / 1};
        for (unsigned c = 0; c < channels; ++c) {
            auto& s = info.add_signal("EEG " + std::to_string(c), rate);
            s.unit = "uV";
            s.physical_min = -200;
            s.physical_max = 300;
        }
        {
            auto writer = fastbdf::Writer::create(argv[2], info);
            std::vector<std::int32_t> data(channels * rate);
            const std::uint64_t span = bdf ? 16777216 : 65536;
            const std::int32_t low = bdf ? -8388608 : -32768;
            for (unsigned r = 0; r < records; ++r) {
                for (unsigned c = 0; c < channels; ++c)
                    for (unsigned i = 0; i < rate; ++i)
                        data[c * rate + i] =
                            low + static_cast<std::int32_t>(((std::uint64_t(r) * rate + i) * 7919 + c * 104729) % span);
                writer.write_digital_record(data);
            }
            writer.finish();
        }
        std::vector<std::size_t> ids(channels);
        std::iota(ids.begin(), ids.end(), 0);
        std::vector<std::vector<double>> output;
        std::vector<double> seconds;
        for (unsigned run = 0; run <= repeats; ++run) {
            output.clear();
            const auto begin = std::chrono::steady_clock::now();
            {
                auto reader = fastbdf::Reader::open(argv[2], fastbdf::AnnotationMode::all);
                output = reader.read_physical(ids, {{}, std::chrono::seconds{records}});
        // MNE returns SI volts; fastbdf preserves header units.
                for (auto& channel : output)
                    for (auto& value : channel)
                        value *= 1e-6;
            }
            const auto elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
            if (run) seconds.push_back(elapsed);
        }
    // Export every value outside the timed region for Python validation.
        std::ofstream binary(argv[3], std::ios::binary);
        for (const auto& channel : output) {
            if (channel.size() != rate * records) throw std::runtime_error("unexpected output length");
            binary.write(reinterpret_cast<const char*>(channel.data()),
                         static_cast<std::streamsize>(channel.size() * sizeof(double)));
        }
        binary.close();
        if (!binary) throw std::runtime_error("output write failed");
        std::sort(seconds.begin(), seconds.end());
        std::cout << std::setprecision(17) << "{\"seconds\":" << seconds[repeats / 2] << ",\"channels\":" << channels
                  << ",\"samples\":" << rate * records << "}\n";
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
