"""Benchmark the three detection-event wire formats.

Stage 0 measured ~40 us per emitted bit when dumping 1500 loose bits per shot
and concluded result streaming dominates run time. This benchmark tests whether
that still holds for the memory experiment. It does not: against a no-emission
control, bit-packed emission costs only 3-4% of run time, and the cost is
dominated by RNG draws and gate operations instead. The "none" format is that
control.

The choice is a crossover, not a single measurement:

* **sparse** costs one word per *fired detector*
* **packed** costs ``ceil(n_checks / 64)`` words per *round*, whatever fires
* **dense** costs ``n_checks`` booleans per round

so sparse wins only when the number of fired detectors per round falls below
``ceil(n_checks / 64)`` -- **1** at d = 11 (60 checks) and **6** at d = 27
(364 checks). This script measures both the realised detection rate and the
wall-clock cost across the Fig. 2 parameter range (p = 0.015, lam = 1,
dt spanning 0.01 to 10), at both distances.

`n_rounds` is held fixed and small: the question is cost *per round*, and using
the physical `round(T/dt)` would confound that with a 1000x change in round
count. This is a benchmark, not a sweep.

Usage:
    python scripts/stage2_format_benchmark.py [--shots N] [--rounds N] [--quick]
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from qec_timing.noise import NoiseParams
from qec_timing.layout import RotatedSurfaceCodeZSector
from qec_timing.stub_circuit import WORD_BITS, build_memory_program

OUT_DIR = Path(__file__).resolve().parents[1] / "results" / "stage2"
BASE_SEED = 20260927

# Fig. 2 / Fig. 9 parameters.
P = 0.015
LAM = 1.0
B_READ = 1.0
FORMATS = ("none", "packed", "sparse", "dense")

# Per-distance budgets: d=27 is ~4x the cost of d=11 per shot.
BUDGET = {11: dict(shots=40, rounds=10), 27: dict(shots=12, rounds=10)}


@dataclass
class Row:
    d: int
    n_checks: int
    emit: str
    dt: float
    p_data: float
    fired_per_round: float
    detection_rate: float
    words_per_round: float
    shots: int
    rounds: int
    build_s: float
    run_s: float
    ms_per_shot: float
    us_per_round: float


def words_per_round(emit: str, n_checks: int, fired: float) -> float:
    """Result-stream entries per round, in 64-bit words."""
    if emit == "packed":
        return math.ceil(n_checks / WORD_BITS)
    if emit == "sparse":
        return fired
    if emit == "dense":
        return float(n_checks)  # one entry per check
    return 0.0  # "none" control emits nothing


def benchmark(distances, dts, shots_override=None, rounds_override=None) -> list[Row]:
    rows: list[Row] = []
    for d in distances:
        budget = BUDGET.get(d, dict(shots=20, rounds=10))
        shots = shots_override or budget["shots"]
        rounds = rounds_override or budget["rounds"]
        layout = RotatedSurfaceCodeZSector(d)
        for emit in FORMATS:
            t0 = time.perf_counter()
            program = build_memory_program(layout, emit=emit)
            build_s = time.perf_counter() - t0

            for dt in dts:
                params = NoiseParams(p=P, lam=LAM, dt=dt, b_read=B_READ)
                t1 = time.perf_counter()
                result = program.run(
                    p_data=params.p_data,
                    p_read=params.p_read,
                    n_rounds=rounds,
                    base_seed=BASE_SEED,
                    shots=shots,
                )
                if emit == "none":
                    # Control: nothing to unpack; force the stream to drain.
                    _ = [shot.collate_tags() for shot in result.results]
                    run_s = time.perf_counter() - t1
                    fired, rate = float("nan"), float("nan")
                else:
                    events = program.detection_events(result, rounds)
                    run_s = time.perf_counter() - t1
                    stacked = np.stack(events)  # (shots, rounds+1, n_checks)
                    fired = float(stacked.sum(axis=2).mean())
                    rate = float(stacked.mean())

                rows.append(
                    Row(
                        d=d,
                        n_checks=layout.n_checks,
                        emit=emit,
                        dt=dt,
                        p_data=params.p_data,
                        fired_per_round=fired,
                        detection_rate=rate,
                        words_per_round=words_per_round(emit, layout.n_checks, fired),
                        shots=shots,
                        rounds=rounds,
                        build_s=build_s,
                        run_s=run_s,
                        ms_per_shot=run_s / shots * 1e3,
                        us_per_round=run_s / shots / (rounds + 1) * 1e6,
                    )
                )
                print(
                    f"  d={d:<3} {emit:<7} dt={dt:<7} "
                    f"fired/round={fired:7.2f} "
                    f"words/round={rows[-1].words_per_round:7.2f} "
                    f"{rows[-1].ms_per_shot:8.2f} ms/shot",
                    flush=True,
                )
    return rows


def summarise(rows: list[Row]) -> dict:
    out: dict = {"crossover": {}, "winner": {}}
    for d in sorted({r.d for r in rows}):
        n_checks = next(r.n_checks for r in rows if r.d == d)
        packed_words = math.ceil(n_checks / WORD_BITS)
        fired = [r.fired_per_round for r in rows if r.d == d and r.emit == "packed"]
        baseline = {r.dt: r.ms_per_shot for r in rows if r.d == d and r.emit == "none"}
        out["crossover"][str(d)] = {
            "n_checks": n_checks,
            "packed_words_per_round": packed_words,
            "sparse_wins_below_fired_per_round": packed_words,
            "observed_fired_per_round_min": min(fired),
            "observed_fired_per_round_max": max(fired),
            "sparse_ever_below_crossover": min(fired) < packed_words,
        }
        for dt in sorted({r.dt for r in rows if r.d == d}):
            cell = {r.emit: r.ms_per_shot for r in rows if r.d == d and r.dt == dt}
            emitting = {k: v for k, v in cell.items() if k != "none"}
            best = min(emitting, key=emitting.get)
            base = baseline.get(dt, float("nan"))
            out["winner"][f"d={d},dt={dt}"] = {
                "fastest": best,
                "ms_per_shot": cell,
                "no_emission_baseline_ms": base,
                "emission_cost_ms": {k: v - base for k, v in emitting.items()},
                "emission_share_of_runtime": {
                    k: (v - base) / v for k, v in emitting.items()
                },
                "packed_speedup_vs_sparse": cell["sparse"] / cell["packed"],
                "packed_speedup_vs_dense": cell["dense"] / cell["packed"],
            }
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shots", type=int, default=None)
    parser.add_argument("--rounds", type=int, default=None)
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()
    args.shots_set = args.shots is not None
    args.rounds_set = args.rounds is not None

    distances = (11,) if args.quick else (11, 27)
    dts = (0.01, 10.0) if args.quick else (0.01, 1.0, 10.0)
    shots = 6 if args.quick else args.shots
    rounds = 5 if args.quick else args.rounds

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(
        f"Format benchmark: d={distances}, dt={dts}, p={P}, lam={LAM}; "
        f"budgets {BUDGET}"
    )
    rows = benchmark(
        distances,
        dts,
        shots_override=shots if args.quick or args.shots_set else None,
        rounds_override=rounds if args.quick or args.rounds_set else None,
    )
    summary = summarise(rows)

    suffix = "_quick" if args.quick else ""
    with (OUT_DIR / f"format_benchmark{suffix}.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0])))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))
    (OUT_DIR / f"format_benchmark{suffix}_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    print("\nCrossover analysis:")
    for d, info in summary["crossover"].items():
        print(
            f"  d={d}: {info['n_checks']} checks -> packed costs "
            f"{info['packed_words_per_round']} word(s)/round; sparse would need "
            f"< {info['sparse_wins_below_fired_per_round']} fired/round. "
            f"Observed {info['observed_fired_per_round_min']:.1f}-"
            f"{info['observed_fired_per_round_max']:.1f}. "
            f"Sparse ever competitive: {info['sparse_ever_below_crossover']}"
        )
    print("\nFastest format per configuration:")
    for key, info in summary["winner"].items():
        print(
            f"  {key:<16} -> {info['fastest']:<7} "
            f"(packed is {info['packed_speedup_vs_sparse']:.2f}x sparse, "
            f"{info['packed_speedup_vs_dense']:.2f}x dense)"
        )


if __name__ == "__main__":
    main()
