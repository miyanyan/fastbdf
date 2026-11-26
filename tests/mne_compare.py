"""Differential EDF/BDF checks. Usage: python tests/mne_compare.py <mne_probe>."""

import argparse
import json
import subprocess
import tempfile
from pathlib import Path

import mne
import numpy as np


def probe(executable, command, path):
    result = subprocess.run(
        [str(executable), command, str(path)],
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    return json.loads(result.stdout) if command in ("inspect", "windows") else None


def compare_annotations(actual, expected):
    assert list(actual.description) == [a["text"] for a in expected]
    np.testing.assert_allclose(
        actual.onset, [a["onset"] for a in expected], rtol=0, atol=1e-7
    )
    np.testing.assert_allclose(
        actual.duration, [a["duration"] for a in expected], rtol=0, atol=1e-7
    )


def compare(executable, path, external=False):
    native = probe(executable, "windows" if external else "inspect", path)
    read = mne.io.read_raw_bdf if path.suffix.lower() == ".bdf" else mne.io.read_raw_edf
    if not native["channels"]:
        compare_annotations(mne.read_annotations(path), native["annotations"])
        print(
            f"PASS {path}: annotation-only, {len(native['annotations'])} events; 0 channels/0 samples"
        )
        return native
    with read(path, preload=False, verbose="ERROR") as raw:
        assert raw.ch_names == [c["name"] for c in native["channels"]]
        assert len(raw.ch_names) == len(native["channels"])
        compare_annotations(raw.annotations, native["annotations"])
    compared = samples = 0
    skipped = []
    scales = {"uV": 1e-6, "µV": 1e-6, "μV": 1e-6, "mV": 1e-3, "nV": 1e-9, "V": 1.0}
    # Opening only one channel preserves its native grid instead of MNE's
    # mixed-rate upsampling. Each read is bounded for external recordings.
    for channel in native["channels"]:
        if channel["unit"] not in scales:
            if not external:
                raise AssertionError(f"unexpected fixture unit {channel['unit']!r}")
            skipped.append(f"{channel['name']} ({channel['unit']!r})")
            continue
        scale = scales[channel["unit"]]
        with read(
            path, include=[channel["name"]], preload=False, verbose="ERROR"
        ) as raw:
            assert raw.ch_names == [channel["name"]]
            assert raw.n_times == channel["count"]
            assert raw.info["sfreq"] == channel["rate"]
            if not external:
                expected = np.asarray(channel["values"]) * scale
                np.testing.assert_allclose(
                    raw.get_data()[0], expected, rtol=2e-12, atol=1e-15
                )
            for window in channel["windows"]:
                first, values = window["first"], window["values"]
                np.testing.assert_allclose(
                    raw.get_data(start=first, stop=first + len(values))[0],
                    np.asarray(values) * scale,
                    rtol=2e-12,
                    atol=1e-15,
                )
                samples += len(values)
            compared += 1
    if external:
        print(
            f"PASS {path}: compared {compared} channels/{samples} window samples, "
            f"{len(native['annotations'])} events; skipped {len(skipped)}: "
            + (", ".join(skipped) if skipped else "none")
        )
    return native


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("probe", type=Path)
    parser.add_argument("--file", action="append", default=[], type=Path)
    args = parser.parse_args()
    executable = args.probe.resolve()
    with tempfile.TemporaryDirectory(prefix="fastbdf-mne-") as directory:
        directory = Path(directory)
        probe(executable, "generate", directory)
        for stem in ("edf", "bdf", "edf_plus", "bdf_plus"):
            path = directory / (stem + (".bdf" if stem.startswith("bdf") else ".edf"))
            result = compare(executable, path)
            assert result["records"] == 4
            assert result["record_duration"] == 0.5
            assert [c["count"] for c in result["channels"]] == [64, 28]
            # Independent fixture oracle catches matching reader/writer mistakes.
            digital_min = -8388608 if stem.startswith("bdf") else -32768
            digital_max = -digital_min - 1
            for channel, per_record, bound in zip(
                result["channels"], [16, 7], [100, 0.1]
            ):
                digital = np.array(
                    [
                        digital_min
                        if i % 5 == 0
                        else digital_max
                        if i % 5 == 1
                        else r * 31 + i - 40
                        for r in range(4)
                        for i in range(per_record)
                    ]
                )
                expected = (digital - digital_min) * (
                    2 * bound / (digital_max - digital_min)
                ) - bound
                np.testing.assert_allclose(
                    channel["values"], expected, rtol=2e-12, atol=1e-13
                )
            assert len(result["annotations"]) == (2 if "plus" in stem else 0)
            if "plus" in stem:
                assert result["annotations"][0]["text"] == "event α"
                np.testing.assert_allclose(
                    [a["onset"] for a in result["annotations"]], [0.25, 1.25]
                )
                np.testing.assert_allclose(
                    [a["duration"] for a in result["annotations"]], [0.125, 0]
                )
            print(f"PASS fastbdf -> MNE: {stem}")
        info = mne.create_info(["EEG 001", "EEG 002"], sfreq=32, ch_types="eeg")
        time = np.arange(64) / 32
        data = np.array(
            [20e-6 * np.sin(2 * np.pi * time), 30e-6 * np.cos(4 * np.pi * time)]
        )
        raw = mne.io.RawArray(data, info, verbose="ERROR")
        raw.set_annotations(mne.Annotations([0.25, 1.25], [0.125, 0], ["event", "end"]))
        exported = directory / "mne_export.edf"
        mne.export.export_raw(
            exported, raw, fmt="edf", physical_range=(-100, 100), verbose="ERROR"
        )
        result = compare(executable, exported)
        for channel, expected in zip(result["channels"], data):
            np.testing.assert_allclose(
                np.asarray(channel["values"]) * 1e-6,
                expected,
                rtol=0,
                atol=100e-6 / 32767 + 1e-12,
            )
        assert [a["text"] for a in result["annotations"]] == ["event", "end"]
        print("PASS MNE/edfio -> fastbdf: EDF+ (quantization bounded)")

    for path in args.file:
        compare(executable, path.resolve(), external=True)


if __name__ == "__main__":
    main()
