"""Light cross-check: Guppy path vs Stage 1 reference sampler.

**Not** the full Stage 4 comparison. This compares one scalar -- the mean
detection-event rate per round -- at small distance and a handful of rounds, to
catch a gross disagreement between the two implementations early.

Both paths are driven from the same ``(p, lam, dt, b_read)`` and the same
geometry, but the noise probabilities are computed by two independent
implementations of Eqs. 1-3 (``noise/params.py`` and ``reference/noise.py``),
so this exercises the whole chain rather than a shared helper.

Only the ``n_rounds`` round blocks are compared; the final data-readout block
has different statistics (it involves the data readout rather than an ancilla
measurement) and is excluded on both sides.

Usage:
    python scripts/stage2_crosscheck.py [--shots N]
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from qec_timing.noise import NoiseParams
from qec_timing.reference.circuit import build_memory_circuit
from qec_timing.layout import RotatedSurfaceCodeZSector
from qec_timing.reference.sampler import clopper_pearson
from qec_timing.stub_circuit import build_memory_program

OUT_DIR = Path(__file__).resolve().parents[1] / "results" / "stage2"
BASE_SEED = 20260927
CONFIDENCE = 0.999

CASES = [
    # (d, p, lam, dt, b_read, n_rounds)
    (3, 0.010, 1.0, 0.5, 1.0, 6),
    (3, 0.015, 1.0, 2.0, 1.0, 6),
    (5, 0.010, 1.0, 0.5, 1.0, 8),
    (5, 0.015, 0.0, 1.0, 1.0, 8),
    (5, 0.005, 3.0, 1.0, 2.0, 8),
]


def guppy_rate(layout, params, n_rounds, shots):
    """Mean detection-event rate per check per round on the Guppy path."""
    program = build_memory_program(layout, emit="packed")
    t0 = time.perf_counter()
    result = program.run(
        p_data=params.p_data,
        p_read=params.p_read,
        n_rounds=n_rounds,
        base_seed=BASE_SEED,
        shots=shots,
    )
    events = np.stack(program.detection_events(result, n_rounds))
    elapsed = time.perf_counter() - t0
    rounds_only = events[:, :n_rounds, :]  # drop the final readout block
    fired = int(rounds_only.sum())
    trials = int(rounds_only.size)
    return fired, trials, elapsed


def reference_rate(layout, params, n_rounds, shots):
    """Same quantity from the Stage 1 Stim reference circuit."""
    circuit = build_memory_circuit(
        layout,
        [params.p_data] * n_rounds,
        params.p_read,
        with_coords=False,
    )
    t0 = time.perf_counter()
    detectors = circuit.compile_detector_sampler(seed=BASE_SEED).sample(shots=shots)
    elapsed = time.perf_counter() - t0
    rounds_only = detectors[:, : n_rounds * layout.n_checks]
    fired = int(np.count_nonzero(rounds_only))
    trials = int(rounds_only.size)
    return fired, trials, elapsed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shots", type=int, default=400)
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    print(
        f"{'d':>2} {'p':>7} {'lam':>4} {'dt':>5} {'b':>4} "
        f"{'guppy rate':>22} {'reference rate':>22}  verdict"
    )
    for d, p, lam, dt, b_read, n_rounds in CASES:
        layout = RotatedSurfaceCodeZSector(d)
        params = NoiseParams(p=p, lam=lam, dt=dt, b_read=b_read)

        g_fired, g_trials, g_s = guppy_rate(layout, params, n_rounds, args.shots)
        r_fired, r_trials, r_s = reference_rate(
            layout, params, n_rounds, args.shots
        )

        g_lo, g_hi = clopper_pearson(g_fired, g_trials, CONFIDENCE)
        r_lo, r_hi = clopper_pearson(r_fired, r_trials, CONFIDENCE)
        g_rate, r_rate = g_fired / g_trials, r_fired / r_trials
        # Agreement = the two exact intervals overlap.
        agree = (g_lo <= r_hi) and (r_lo <= g_hi)

        rows.append(
            {
                "d": d, "p": p, "lam": lam, "dt": dt, "b_read": b_read,
                "n_rounds": n_rounds, "shots": args.shots,
                "p_data": params.p_data, "p_read": params.p_read,
                "guppy_fired": g_fired, "guppy_trials": g_trials,
                "guppy_rate": g_rate, "guppy_ci": [g_lo, g_hi],
                "reference_fired": r_fired, "reference_trials": r_trials,
                "reference_rate": r_rate, "reference_ci": [r_lo, r_hi],
                "relative_difference": (g_rate - r_rate) / r_rate if r_rate else 0.0,
                "intervals_overlap": bool(agree),
                "guppy_seconds": g_s, "reference_seconds": r_s,
            }
        )
        print(
            f"{d:>2} {p:>7} {lam:>4} {dt:>5} {b_read:>4} "
            f"{g_rate:.5f} [{g_lo:.5f},{g_hi:.5f}] "
            f"{r_rate:.5f} [{r_lo:.5f},{r_hi:.5f}]  "
            f"{'OK' if agree else 'MISMATCH'} "
            f"({100 * rows[-1]['relative_difference']:+.2f}%)"
        )

    (OUT_DIR / "crosscheck.json").write_text(
        json.dumps(rows, indent=2), encoding="utf-8"
    )
    n_bad = sum(not r["intervals_overlap"] for r in rows)
    print(f"\n{len(rows) - n_bad}/{len(rows)} cases agree at {CONFIDENCE:.1%}")
    if n_bad:
        print("Disagreements reported, not debugged (Stage 4 owns the full sweep).")


if __name__ == "__main__":
    main()
