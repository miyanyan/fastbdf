"""Python binding correctness checks: python python_binding.py <mne_probe>."""

import gc
import subprocess
import sys
import tempfile
from pathlib import Path

import fastbdf
import mne
import numpy as np


def rejects(call, *args):
    try:
        call(*args)
    except (TypeError, ValueError, IndexError, OverflowError, RuntimeError):
        return
    raise AssertionError("invalid argument was accepted")


def check_file(path):
    reader = fastbdf.Reader(str(path))
    signals = reader.signals
    assert isinstance(signals, list)
    assert [s.label for s in signals] == ["channel0", "channel1"]
    assert [s.unit for s in signals] == ["uV", "V"]
    annotations = reader.annotations
    assert annotations == (
        [(2500000, 1250000, "event α".encode()), (12500000, None, b"end")]
        if "plus" in path.stem
        else []
    )
    minimum = -8388608 if path.suffix == ".bdf" else -32768
    maximum = -minimum - 1
    expected_physical = []
    for channel, (per_record, bound, rate) in enumerate([(16, 100, 32), (7, 0.1, 14)]):
        signal = signals[channel]
        assert signal.digital_min == minimum
        assert signal.digital_max == maximum
        count = 4 * per_record
        assert reader.sample_count(channel) == count
        assert reader.sample_rate(channel) == rate
        digital = np.array(
            [
                minimum
                if i % 5 == 0
                else maximum
                if i % 5 == 1
                else record * 31 + i - 40
                for record in range(4)
                for i in range(per_record)
            ],
            dtype=np.int32,
        )
        physical = (digital.astype(np.float64) - minimum) * (
            2 * bound / (maximum - minimum)
        ) - bound
        expected_physical.append(physical)
        for method, dtype, expected in (
            (reader.read_digital_into, np.int32, digital),
            (reader.read_physical_into, np.float64, physical),
        ):
            output = np.full(count, 123, dtype=dtype)
            assert method(channel, 0, output) == count
            np.testing.assert_allclose(output, expected, rtol=2e-12, atol=1e-13)
            # Cross a record boundary and preserve the unwritten EOF tail.
            output = np.full(9, 123, dtype=dtype)
            assert method(channel, per_record - 2, output) == 9
            np.testing.assert_allclose(
                output,
                expected[per_record - 2 : per_record + 7],
                rtol=2e-12,
                atol=1e-13,
            )
            output.fill(123)
            assert method(channel, count - 2, output) == 2
            np.testing.assert_allclose(
                output[:2], expected[-2:], rtol=2e-12, atol=1e-13
            )
            np.testing.assert_array_equal(output[2:], 123)
            assert method(channel, count, output) == 0
            np.testing.assert_array_equal(output[2:], 123)
            assert method(channel, 0, np.empty(0, dtype=dtype)) == 0
            rejects(method, channel, -1, output)
            rejects(method, 2, 0, output)
            rejects(method, channel, 0, np.empty(8, dtype=np.float32))
            rejects(method, channel, 0, np.empty(16, dtype=dtype)[::2])
            rejects(method, channel, 0, np.empty((2, 4), dtype=dtype))
            readonly = np.empty(8, dtype=dtype)
            readonly.flags.writeable = False
            rejects(method, channel, 0, readonly)
        read = mne.io.read_raw_bdf if path.suffix == ".bdf" else mne.io.read_raw_edf
        with read(path, include=[signal.label], preload=False, verbose="ERROR") as raw:
            assert raw.info["sfreq"] == rate
            assert raw.n_times == count
            scale = 1e-6 if signal.unit == "uV" else 1
            np.testing.assert_allclose(
                raw.get_data()[0], physical * scale, rtol=2e-12, atol=1e-15
            )

    ids = [1, 0, 1]
    for start, stop in [
        (0, 20000000),
        (2500000, 12500000),
        (1, 5000000),
        (5000000, 5000000),
        (20000000, 30000000),
        (19000000, 30000000),
    ]:
        result = reader.read_physical_window(ids, start, stop)
        assert isinstance(result, list) and len(result) == len(ids)
        for channel, output in zip(ids, result):
            rate = [32, 14][channel]
            first = (start * rate + 9999999) // 10000000
            last = (stop * rate + 9999999) // 10000000
            assert isinstance(output, np.ndarray)
            assert output.dtype == np.float64 and output.ndim == 1
            assert output.flags.writeable
            np.testing.assert_allclose(
                output, expected_physical[channel][first:last], rtol=2e-12, atol=1e-13
            )
    assert reader.read_physical_window([], 0, 20000000) == []
    for channel in [-1, 2]:
        rejects(reader.sample_rate, channel)
        rejects(reader.sample_count, channel)
        rejects(reader.read_physical_window, [channel], 0, 1)
    rejects(reader.read_physical_window, [0], -1, 1)
    rejects(reader.read_physical_window, [0], 2, 1)
    # Rectangular output is valid when every selected native grid has the
    # requested number of samples; repeated channels retain selection order.
    for channel in [0, 1]:
        rate = [32, 14][channel]
        rectangular = np.empty((2, rate), dtype=np.float64)
        assert (
            reader.read_physical_window_into(
                [channel, channel], 0, 10000000, rectangular
            )
            is None
        )
        np.testing.assert_allclose(
            rectangular,
            np.tile(expected_physical[channel][:rate], (2, 1)),
            rtol=2e-12,
            atol=1e-13,
        )
    assert (
        reader.read_physical_window_into(
            [], 0, 10000000, np.empty((0, 0), dtype=np.float64)
        )
        is None
    )
    mixed = np.empty((3, 1), dtype=np.float64)
    assert reader.read_physical_window_into([1, 0, 1], 0, 1, mixed) is None
    np.testing.assert_allclose(
        mixed[:, 0],
        [expected_physical[i][0] for i in [1, 0, 1]],
        rtol=2e-12,
        atol=1e-13,
    )
    assert (
        reader.read_physical_window_into(
            [1, 0, 1], 5000000, 5000000, np.empty((3, 0), dtype=np.float64)
        )
        is None
    )
    valid = np.empty((1, 32), dtype=np.float64)
    readonly = valid.copy()
    readonly.flags.writeable = False
    unaligned = np.ndarray(
        (1, 32), dtype=np.float64, buffer=bytearray(32 * 8 + 1), offset=1
    )
    for channels, start, stop, output in [
        ([0, 1], 0, 10000000, np.empty((2, 32))),
        ([0], 0, 10000000, np.empty((2, 32))),
        ([0], 0, 10000000, np.empty((1, 31))),
        ([0], 0, 10000000, np.empty(32)),
        ([0], 0, 10000000, np.empty((1, 32), dtype=np.float32)),
        ([0], 0, 10000000, np.empty((1, 64))[:, ::2]),
        ([0], 0, 10000000, readonly),
        ([0], 0, 10000000, unaligned),
        ([2], 0, 10000000, valid),
        ([-1], 0, 10000000, valid),
        ([0], -1, 10000000, valid),
        ([0], 10000000, 0, valid),
    ]:
        rejects(reader.read_physical_window_into, channels, start, stop, output)
    retained = reader.read_physical_window(ids, 0, 20000000)
    reader.close()
    reader.close()  # Closing twice is harmless.
    rejects(reader.sample_count, 0)
    rejects(reader.sample_rate, 0)
    rejects(getattr, reader, "signals")
    rejects(getattr, reader, "annotations")
    rejects(reader.read_physical_into, 0, 0, np.empty(1))
    rejects(reader.read_physical_window, [0], 0, 1)
    rejects(reader.read_physical_window_into, [0], 0, 10000000, valid)
    del method
    del reader
    gc.collect()
    assert [s.label for s in signals] == ["channel0", "channel1"]
    assert signals[0].digital_min == minimum
    assert annotations == (
        [(2500000, 1250000, "event α".encode()), (12500000, None, b"end")]
        if "plus" in path.stem
        else []
    )
    for channel, output in zip(ids, retained):
        np.testing.assert_allclose(
            output, expected_physical[channel], rtol=2e-12, atol=1e-13
        )
    assert not np.shares_memory(retained[0], retained[2])
    retained[0][0] = 999
    assert retained[2][0] == expected_physical[1][0]
    print(f"PASS Python binding: {path.name}")


def main():
    executable = Path(sys.argv[1]).resolve()
    with tempfile.TemporaryDirectory(prefix="fastbdf-python-") as directory:
        subprocess.run([str(executable), "generate", directory], check=True)
        for stem in ["edf", "bdf", "edf_plus", "bdf_plus"]:
            suffix = ".bdf" if stem.startswith("bdf") else ".edf"
            check_file(Path(directory) / (stem + suffix))
        # Also consume an independently produced EDF+ fixture (edfio via MNE).
        exported = Path(directory) / "mne_export.edf"
        data = np.array([20e-6 * np.sin(np.arange(64) / 32 * 2 * np.pi)])
        raw = mne.io.RawArray(
            data, mne.create_info(["EEG"], 32, "eeg"), verbose="ERROR"
        )
        mne.export.export_raw(
            exported, raw, fmt="edf", physical_range=(-100, 100), verbose="ERROR"
        )
        reader = fastbdf.Reader(str(exported))
        output = np.empty(64, dtype=np.float64)
        assert reader.read_physical_into(0, 0, output) == 64
        np.testing.assert_allclose(
            output * 1e-6, data[0], rtol=0, atol=100e-6 / 32767 + 1e-12
        )
        print("PASS Python binding: MNE/edfio export")
        rejects(fastbdf.Reader, str(Path(directory) / "missing.edf"))


if __name__ == "__main__":
    main()
