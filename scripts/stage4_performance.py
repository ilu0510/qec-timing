"""Measure per-round cost on both paths and extrapolate the Fig. 2 sweep.

Three measurements, then a cost model:

1. **Guppy/Selene marginal cost per round** at every distance d = 11..27. The
   marginal cost is cheap to measure even at d = 27 -- it is the *shot count*
   that makes the sweep expensive, not the round count -- so this is measured
   directly rather than extrapolated from a fitted scaling law.
2. **Reference (Stim + pymatching) cost per detector**, same distances.
3. **RNG cost per draw**, via the Stage 2 ``record_flip_probe``, which is a bare
   Bernoulli draw with no quantum operations. Combined with the per-round
   operation count this splits the Guppy round cost into RNG vs everything else,
   which is the evidence for whether geometric-gap sampling would pay.

The cost model then asks: how long would one dt point per d cost on each path,
at the Fig. 2 targets, with enough shots to resolve the minimum?

Usage:
    python scripts/stage4_performance.py [--quick]
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from qec_timing import ansatz
from qec_timing.layout import RotatedSurfaceCodeZSector
from qec_timing.noise import NoiseParams
from qec_timing.noise.probes import record_flip_probe
from qec_timing.reference.circuit import build_memory_circuit
from qec_timing.reference.decoder import build_decoder
from qec_timing.stub_circuit import build_memory_program

OUT_DIR = Path(__file__).resolve().parents[1] / "results" / "stage4a"
BASE_SEED = 20260927

# Fig. 2 / Fig. 9 targets.
P, LAM, B_READ, T = 0.015, 1.0, 1.0, 200.0
P_TH, BETA, G = 0.029, 2.0, 0.8

DISTANCES = (11, 13, 15, 17, 19, 21, 23, 25, 27)

# Paper-fitted A where quoted (Fig. 5), linearly interpolated between.
# 4a fits A independently at d = 11, 13, 15 and reports the comparison.
A_BY_D = {11: 1.4, 13: 1.45, 15: 1.5, 17: 1.6, 19: 1.75,
          21: 1.9, 23: 1.95, 25: 2.0, 27: 2.05}

# Shots per point are sized to this many logical failures.
TARGET_FAILURES = 400


def measure_guppy(distances, rounds_pair=(50, 200), shots=12):
    """Marginal ms/round and fixed ms/shot for the Guppy path."""
    out = {}
    params = NoiseParams(p=P, lam=LAM, dt=1.0, b_read=B_READ)
    for d in distances:
        layout = RotatedSurfaceCodeZSector(d)
        program = build_memory_program(layout, emit="packed")
        timings = []
        for n_rounds in rounds_pair:
            t0 = time.perf_counter()
            result = program.run(
                p_data=params.p_data,
                p_read=params.p_read,
                dt=1.0,
                n_rounds=n_rounds,
                base_seed=BASE_SEED,
                shots=shots,
            )
            program.detection_events(result, n_rounds)
            timings.append((time.perf_counter() - t0) / shots * 1e3)
        (lo_r, hi_r), (lo_t, hi_t) = rounds_pair, timings
        per_round = (hi_t - lo_t) / (hi_r - lo_r)
        fixed = lo_t - per_round * lo_r
        ops = layout.n_data + layout.n_checks  # RNG draws per round
        out[d] = {
            "n_data": layout.n_data,
            "n_checks": layout.n_checks,
            "rng_draws_per_round": ops,
            "ms_per_round": per_round,
            "fixed_ms_per_shot": fixed,
            "us_per_draw_equivalent": per_round * 1e3 / ops,
            "raw": dict(zip(map(str, rounds_pair), timings)),
        }
        print(
            f"  guppy d={d:<3} {per_round:8.4f} ms/round  "
            f"(fixed {fixed:6.1f} ms/shot, {ops} draws/round)",
            flush=True,
        )
    return out


def measure_reference(distances, n_rounds=200, shots=200):
    """ms per detector, and the implied ms/round, for the reference path."""
    out = {}
    params = NoiseParams(p=P, lam=LAM, dt=1.0, b_read=B_READ)
    for d in distances:
        layout = RotatedSurfaceCodeZSector(d)
        circuit = build_memory_circuit(
            layout, [params.p_data] * n_rounds, params.p_read, with_coords=False
        )
        decoder = build_decoder(circuit)
        t0 = time.perf_counter()
        detectors, _ = circuit.compile_detector_sampler(
            seed=BASE_SEED
        ).sample(shots=shots, separate_observables=True)
        decoder.decode_batch(detectors)
        ms_per_shot = (time.perf_counter() - t0) / shots * 1e3
        out[d] = {
            "n_detectors": circuit.num_detectors,
            "ms_per_shot_at_200_rounds": ms_per_shot,
            "ns_per_detector": ms_per_shot * 1e6 / circuit.num_detectors,
            "ms_per_round": ms_per_shot / n_rounds,
        }
        print(
            f"  ref   d={d:<3} {out[d]['ms_per_round']:8.4f} ms/round  "
            f"({out[d]['ns_per_detector']:6.1f} ns/detector)",
            flush=True,
        )
    return out


def measure_rng(n_reps=200_000, shots=8):
    """Cost of one bare Bernoulli draw, with no quantum operations."""
    t0 = time.perf_counter()
    result = (
        record_flip_probe.emulator(n_qubits=1)
        .stabilizer_sim()
        .with_shots(shots)
        .run(p_read=0.015, n_reps=n_reps, base_seed=BASE_SEED)
    )
    total = time.perf_counter() - t0
    draws = sum(s.as_dict()["trials"] for s in result.results)
    us = total / draws * 1e6
    print(f"  rng   {us:.4f} us/draw  ({draws:,} draws in {total:.1f} s)")
    return {"us_per_draw": us, "draws": draws, "seconds": total}


def cost_model(guppy, reference, distances):
    """Cost of one dt point per d -- the minimum -- on each path."""
    rows = []
    for d in distances:
        A = A_BY_D[d]
        dt_star = ansatz.optimal_interval(d=d, lam=LAM, g=G)
        n_rounds = max(round(T / dt_star), 1)
        p_L = ansatz.p_L_constant_interval(
            dt_star, n_rounds=n_rounds, d=d, p=P, lam=LAM, A=A,
            p_th=P_TH, beta=BETA, g=G,
        )
        shots = TARGET_FAILURES / p_L
        g, r = guppy[d], reference[d]
        guppy_s = shots * (g["ms_per_round"] * n_rounds + g["fixed_ms_per_shot"]) / 1e3
        ref_s = shots * r["ms_per_round"] * n_rounds / 1e3
        rows.append(
            {
                "d": d, "A_assumed": A, "dt_star": dt_star, "n_rounds": n_rounds,
                "p_L_at_dt_star": p_L, "shots_for_target_failures": shots,
                "guppy_seconds": guppy_s, "reference_seconds": ref_s,
                "ratio": guppy_s / ref_s if ref_s else float("nan"),
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()
    distances = (11, 15, 27) if args.quick else DISTANCES

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("Measuring per-round cost on both paths...")
    guppy = measure_guppy(distances, shots=6 if args.quick else 12)
    reference = measure_reference(distances, shots=50 if args.quick else 200)
    rng = measure_rng(n_reps=20_000 if args.quick else 200_000,
                      shots=4 if args.quick else 8)

    # RNG share of the Guppy round cost.
    for d in distances:
        draws = guppy[d]["rng_draws_per_round"]
        rng_ms = draws * rng["us_per_draw"] / 1e3
        guppy[d]["rng_ms_per_round"] = rng_ms
        guppy[d]["rng_share_of_round"] = rng_ms / guppy[d]["ms_per_round"]

    rows = cost_model(guppy, reference, distances)

    payload = {
        "targets": {"p": P, "lam": LAM, "b_read": B_READ, "T": T,
                    "p_th": P_TH, "beta": BETA, "g": G},
        "target_failures_per_point": TARGET_FAILURES,
        "guppy": guppy, "reference": reference, "rng": rng,
        "cost_model": rows,
    }
    suffix = "_quick" if args.quick else ""
    (OUT_DIR / f"performance{suffix}.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )

    print("\nRNG share of the Guppy per-round cost:")
    for d in distances:
        g = guppy[d]
        print(
            f"  d={d:<3} {g['ms_per_round']:7.4f} ms/round, "
            f"{g['rng_ms_per_round']:7.4f} ms is RNG "
            f"({100 * g['rng_share_of_round']:4.1f}%)"
        )

    print(f"\nCost of ONE dt point per d (the minimum, {TARGET_FAILURES} failures):")
    print(f"{'d':>3} {'dt*':>7} {'p_L*':>9} {'rounds':>7} {'shots':>10} "
          f"{'reference':>12} {'guppy':>14} {'ratio':>8}")
    tot_g = tot_r = 0.0
    for row in rows:
        tot_g += row["guppy_seconds"]
        tot_r += row["reference_seconds"]
        print(
            f"{row['d']:>3} {row['dt_star']:>7.4f} {row['p_L_at_dt_star']:>9.2e} "
            f"{row['n_rounds']:>7} {row['shots_for_target_failures']:>10.3g} "
            f"{_fmt(row['reference_seconds']):>12} {_fmt(row['guppy_seconds']):>14} "
            f"{row['ratio']:>7.0f}x"
        )
    print(f"\nTotal, one point per d: reference {_fmt(tot_r)}, guppy {_fmt(tot_g)}")


def _fmt(seconds: float) -> str:
    if seconds < 90:
        return f"{seconds:.1f} s"
    if seconds < 5400:
        return f"{seconds / 60:.1f} min"
    if seconds < 86400 * 2:
        return f"{seconds / 3600:.1f} h"
    return f"{seconds / 86400:.1f} d"


if __name__ == "__main__":
    main()
