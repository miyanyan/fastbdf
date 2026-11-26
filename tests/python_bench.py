"""Warm-cache Python Reader/MNE benchmark (requires fastbdf, numpy and MNE).

Run: python tests/python_bench.py [--file recording.bdf] [--records 120]
Synthetic data: 64 EEG channels, 512 Hz, one-second BDF records. Timings include
allocated float64 output and conversion to SI volts. Imports, fixture creation,
verification and output disposal are excluded. This is not a cold-disk benchmark.
"""

import argparse
import ctypes
import json
import os
import platform
import statistics
import subprocess
import sys
import tempfile
import time
from importlib.metadata import version
from pathlib import Path

import fastbdf
import mne
import numpy as np

SCALES = {"V": 1.0, "mV": 1e-3, "uV": 1e-6, "µV": 1e-6, "μV": 1e-6, "nV": 1e-9}
NANOSECONDS_PER_SECOND = 1_000_000_000


def synthetic(path, records, rate=512, channels=64):
    """Independent BDF fixture; avoids an optional BDF export dependency."""

    def field(value, width):
        encoded = str(value).encode("ascii")
        assert len(encoded) <= width
        return encoded.ljust(width, b" ")

    header = b"\xffBIOSEMI" + field("benchmark", 80) + field("benchmark", 80)
    header += field("01.01.24", 8) + field("00.00.00", 8)
    header += field(256 * (channels + 1), 8) + field("", 44)
    header += field(records, 8) + field(1, 8) + field(channels, 4)
    for values, width in (
        ([f"EEG {i + 1:03}" for i in range(channels)], 16),
        ([""] * channels, 80),
        (["uV"] * channels, 8),
        ([-200] * channels, 8),
        ([300] * channels, 8),
        ([-8388608] * channels, 8),
        ([8388607] * channels, 8),
        ([""] * channels, 80),
        ([rate] * channels, 8),
        ([""] * channels, 32),
    ):
        header += b"".join(field(v, width) for v in values)
    with path.open("wb") as output:
        output.write(header)
        channel = np.arange(channels, dtype=np.int64)[:, None]
        for record in range(records):
            sample = np.arange(record * rate, (record + 1) * rate, dtype=np.int64)
            code = (sample * 7919 + channel * 104729) % 16777216 - 8388608
            packed = np.empty((channels, rate, 3), dtype=np.uint8)
            for byte in range(3):
                packed[:, :, byte] = (code >> (byte * 8)) & 255
            output.write(packed.tobytes())


def open_mne(path, names):
    raw = mne.io.read_raw_bdf(str(path), include=names, preload=False, verbose="ERROR")
    raw.reorder_channels(names)
    return raw


def selection(path, requested):
    with fastbdf.Reader(str(path)) as reader:
        candidates = [
            (i, s.label, s.unit)
            for i, s in enumerate(reader.signals)
            if s.unit in SCALES and s.label.casefold() != "status"
        ]
        if requested:
            candidates = [item for item in candidates if item[1] in requested]
            if set(requested) != {item[1] for item in candidates}:
                raise ValueError("requested channels must exist and have voltage units")
        with open_mne(path, [c[1] for c in candidates]) as raw:
            eeg = {
                name
                for name, kind in zip(raw.ch_names, raw.get_channel_types())
                if kind == "eeg"
            }
        candidates = [c for c in candidates if c[1] in eeg]
        if not candidates:
            raise ValueError("no EEG channels with voltage units")
        rate = reader.sample_rate(candidates[0][0])
        count = reader.sample_count(candidates[0][0])
        candidates = [
            c
            for c in candidates
            if reader.sample_rate(c[0]) == rate and reader.sample_count(c[0]) == count
        ]
    if requested and set(requested) != {c[1] for c in candidates}:
        raise ValueError(
            "requested channels must all be EEG with matching sample rate/count"
        )
    if len(candidates) < 4:
        raise ValueError("need at least four equal-rate EEG channels")
    if len({c[1] for c in candidates}) != len(candidates):
        raise ValueError("duplicate channel labels are unsupported")
    with open_mne(path, [c[1] for c in candidates]) as raw:
        assert raw.ch_names == [c[1] for c in candidates]
        assert raw.n_times == count and raw.info["sfreq"] == rate
    return candidates, float(rate), int(count)


def fast_read(reader, channels, start_ns, stop_ns, count):
    result = np.empty((len(channels), count), dtype=np.float64)
    reader.read_physical_window_into(
        [c[0] for c in channels], start_ns, stop_ns, result
    )
    for row, (_, _, unit) in zip(result, channels):
        row *= SCALES[unit]
    return result


def peak_rss():
    if os.name == "nt":
        from ctypes import wintypes

        class Counters(ctypes.Structure):
            _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD)] + [
                (name, ctypes.c_size_t)
                for name in (
                    "PeakWorkingSetSize",
                    "WorkingSetSize",
                    "QuotaPeakPagedPoolUsage",
                    "QuotaPagedPoolUsage",
                    "QuotaPeakNonPagedPoolUsage",
                    "QuotaNonPagedPoolUsage",
                    "PagefileUsage",
                    "PeakPagefileUsage",
                )
            ]

        counters = Counters()
        counters.cb = ctypes.sizeof(counters)
        kernel = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel.GetCurrentProcess.restype = wintypes.HANDLE
        psapi = ctypes.WinDLL("psapi", use_last_error=True)
        psapi.GetProcessMemoryInfo.argtypes = [
            wintypes.HANDLE,
            ctypes.POINTER(Counters),
            wintypes.DWORD,
        ]
        if not psapi.GetProcessMemoryInfo(
            kernel.GetCurrentProcess(), ctypes.byref(counters), counters.cb
        ):
            raise ctypes.WinError(ctypes.get_last_error())
        return counters.PeakWorkingSetSize
    import resource

    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def memory_child(payload):
    channels = payload["channels"]
    count, rate = payload["count"], payload["rate"]
    if payload["engine"] == "fastbdf":
        with fastbdf.Reader(payload["path"]) as reader:
            result = fast_read(reader, channels, 0, count * NANOSECONDS_PER_SECOND // int(rate), count)
    else:
        with open_mne(payload["path"], [c[1] for c in channels]) as raw:
            result = raw.get_data(start=0, stop=count)
    assert result.shape == (len(channels), count) and result.dtype == np.float64
    touch = float(result[0, 0] + result[-1, -1])
    print(json.dumps({"peak_rss": peak_rss(), "touch": touch}))


def report_memory(path, channels, rate, count):
    for engine in ("fastbdf", "MNE"):
        payload = {
            "path": str(path),
            "channels": channels,
            "rate": rate,
            "count": count,
            "engine": engine,
        }
        child = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--memory-child"],
            input=json.dumps(payload),
            check=True,
            capture_output=True,
            text=True,
        )
        measured = json.loads(child.stdout)
        print(
            f"open+full peak RSS {engine}: {measured['peak_rss'] / (1024**2):.1f} MiB"
        )


def check(actual, expected):
    assert (
        actual.dtype == expected.dtype == np.float64 and actual.shape == expected.shape
    )
    assert actual.flags.c_contiguous and expected.flags.c_contiguous
    np.testing.assert_allclose(actual, expected, rtol=2e-11, atol=1e-15)


def benchmark(path, args, is_synthetic=False):
    channels, rate, count = selection(path, args.channels)
    if count < round(10 * rate):
        raise ValueError("recording must contain at least ten seconds")
    if not rate.is_integer():
        raise ValueError(
            "benchmark requires integer-Hz EEG sampling for nanosecond window alignment"
        )
    names = [c[1] for c in channels]
    print(f"\n{path}: {len(channels)} EEG channels x {count} samples at {rate:g} Hz")
    print("Selected: " + ", ".join(names))
    times = {
        case: {engine: [] for engine in ("fastbdf", "MNE")}
        for case in ("open+first10s", "random4ch4s", "open+full")
    }
    rng = np.random.default_rng(args.seed)
    length = round(4 * rate)
    windows = [
        (
            rng.choice(len(channels), 4, replace=False).tolist(),
            int(rng.integers(0, count - length + 1)),
        )
        for _ in range(args.windows)
    ]
    # Alternate implementation order per repetition; first pair is warm-up.
    for case in ("open+first10s", "open+full"):
        stop = min(count, round(10 * rate)) if case == "open+first10s" else count
        for run in range(args.repeats + 1):
            results = {}
            for engine in ("fastbdf", "MNE") if run % 2 == 0 else ("MNE", "fastbdf"):
                begin = time.perf_counter()
                if engine == "fastbdf":
                    with fastbdf.Reader(str(path)) as reader:
                        result = fast_read(
                            reader, channels, 0, stop * NANOSECONDS_PER_SECOND // int(rate), stop
                        )
                else:
                    with open_mne(path, names) as raw:
                        result = raw.get_data(start=0, stop=stop)
                elapsed = time.perf_counter() - begin
                results[engine] = result
                if run:
                    times[case][engine].append(elapsed)
            check(results["fastbdf"], results["MNE"])
            if is_synthetic:
                for row, (index, _, _) in zip(results["fastbdf"], channels):
                    code = (
                        np.arange(stop, dtype=np.int64) * 7919 + index * 104729
                    ) % 16777216
                    check(row, (code * (500.0 / 16777215) - 200) * 1e-6)
            del results, result
    with fastbdf.Reader(str(path)) as reader, open_mne(path, names) as raw:
        # Same requests and outputs, resident readers; warm one request each.
        for run, (picks, start) in enumerate([windows[0]] + windows):
            # Floor to the nanosecond before each sample: the reader selects the first
            # native sample at or after that nanosecond. Windows may cross any record.
            start_ns = start * NANOSECONDS_PER_SECOND // int(rate)
            stop_ns = (start + length) * NANOSECONDS_PER_SECOND // int(rate)
            results = {}
            selected = [channels[i] for i in picks]
            for engine in ("fastbdf", "MNE") if run % 2 == 0 else ("MNE", "fastbdf"):
                begin = time.perf_counter()
                if engine == "fastbdf":
                    result = fast_read(reader, selected, start_ns, stop_ns, length)
                else:
                    result = raw.get_data(picks=picks, start=start, stop=start + length)
                elapsed = time.perf_counter() - begin
                results[engine] = result
                if run:
                    times["random4ch4s"][engine].append(elapsed)
            check(results["fastbdf"], results["MNE"])
            del results, result
    for case, engines in times.items():
        values = {
            engine: statistics.median(samples) for engine, samples in engines.items()
        }
        print(
            f"{case}: median fastbdf {values['fastbdf'] * 1000:.3f} ms; "
            f"MNE {values['MNE'] * 1000:.3f} ms; MNE/fastbdf {values['MNE'] / values['fastbdf']:.2f}x"
        )
        if case == "random4ch4s":
            print(
                "  P95: "
                + "; ".join(
                    f"{engine} {np.percentile(samples, 95) * 1000:.3f} ms"
                    for engine, samples in engines.items()
                )
            )
    print("All measured outputs verified in SI volts outside timing.")
    report_memory(path, channels, rate, count)


def main():
    if len(sys.argv) == 2 and sys.argv[1] == "--memory-child":
        memory_child(json.load(sys.stdin))
        return
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--file", type=Path, action="append", default=[])
    parser.add_argument("--channels", nargs="+", help="explicit real-file EEG labels")
    parser.add_argument("--records", type=int, default=120)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--windows", type=int, default=200)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args()
    if args.records < 10 or args.repeats < 1 or args.windows < 1:
        parser.error("records >= 10, repeats >= 1 and windows >= 1 required")
    print(
        f"Python {platform.python_version()} | fastbdf {version('fastbdf')} | "
        f"MNE {mne.__version__} | NumPy {np.__version__} | {platform.platform()}"
    )
    print(
        "Warm OS cache; float64 SI volts; both parse annotations on open; MNE preload=False; close included."
    )
    print(
        "Peak RSS: fresh process per engine/full read, includes imports and native allocations."
    )
    with tempfile.TemporaryDirectory(prefix="fastbdf-python-") as directory:
        path = Path(directory) / "synthetic.bdf"
        synthetic(path, args.records)
        requested = args.channels
        args.channels = None
        benchmark(path, args, is_synthetic=True)
        args.channels = requested
    for path in args.file:
        benchmark(path.resolve(strict=True), args)


if __name__ == "__main__":
    main()
