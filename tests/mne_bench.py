"""Optional warm-cache physical full-window benchmark.

Install numpy and mne, build fastbdf_mne_bench in Release mode, then run:
    python tests/mne_bench.py <path-to-fastbdf_mne_bench>

Continuous EDF+/BDF+, 16 equal-rate EEG channels, 512 Hz, 512 seconds by default.
Both timings include open/header parsing, full allocated float64 output, unit
conversion to volts, and reader close. Imports, process startup, fixture writing,
output export, correctness checks, and output disposal are excluded. One warm-up
and five measurements per implementation; reports medians, not cold-disk speed.
This does not benchmark preprocessing, lazy slices, mixed rates or annotations.
"""

import argparse
import json
import platform
import statistics
import subprocess
import tempfile
import time
from pathlib import Path

import mne
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument(
        "--records",
        type=int,
        default=512,
        help="one-second records per format (1..4096; default 512)",
    )
    args = parser.parse_args()
    if not 1 <= args.records <= 4096:
        parser.error("--records must be in [1, 4096]")
    executable = str(args.executable.resolve(strict=True))
    print(
        f"Python {platform.python_version()} | MNE {mne.__version__} | NumPy {np.__version__} | {platform.platform()}"
    )
    print("Open + full float64 read + SI volts + close; warm-up + median of 5")
    with tempfile.TemporaryDirectory(prefix="fastbdf-mne-") as directory:
        for fmt in ("edf", "bdf"):
            fixture = Path(directory) / ("fixture." + fmt)
            exported = Path(directory) / "native.bin"
            native = json.loads(
                subprocess.check_output(
                    [executable, fmt, str(fixture), str(exported), str(args.records)],
                    text=True,
                )
            )
            shape = (native["channels"], native["samples"])
            fast = np.fromfile(str(exported), dtype=np.float64).reshape(shape)
            reader = mne.io.read_raw_edf if fmt == "edf" else mne.io.read_raw_bdf
            durations = []
            data = None
            for run in range(6):
                data = None  # Dispose previous result outside timing.
                begin = time.perf_counter()
                with reader(str(fixture), preload=False, verbose="ERROR") as raw:
                    data = raw.get_data()
                elapsed = time.perf_counter() - begin
                if run:
                    durations.append(elapsed)
                # Check every channel/value outside timing, including warm-up.
                assert data.dtype == np.float64 and data.shape == shape
                np.testing.assert_allclose(data, fast, rtol=1e-12, atol=1e-15)
            # Independent deterministic digital fixture oracle in SI volts.
            span = 16777216 if fmt == "bdf" else 65536
            for channel in range(shape[0]):
                code = (
                    np.arange(shape[1], dtype=np.int64) * 7919 + channel * 104729
                ) % span
                expected = (
                    code.astype(np.float64) * (500.0 / (span - 1)) - 200.0
                ) * 1e-6
                np.testing.assert_allclose(
                    fast[channel], expected, rtol=1e-12, atol=1e-15
                )
            mne_seconds = statistics.median(durations)
            print(
                "{}+: {} x {} values; fastbdf {:.3f} ms; MNE {:.3f} ms; "
                "MNE/fastbdf {:.2f}x; all values verified".format(
                    fmt.upper(),
                    *shape,
                    native["seconds"] * 1000,
                    mne_seconds * 1000,
                    mne_seconds / native["seconds"],
                )
            )


if __name__ == "__main__":
    main()
