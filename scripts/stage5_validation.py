"""Test syndrome timing with full X/Z extraction and independent confirmation.

Run from the repository root (Python >=3.12):
    python scripts/stage5_validation.py --quick
    python scripts/stage5_validation.py --q-d 2

All configurations have EXACTLY the same T: dt = T / n_rounds. Both check
types are extracted on every patch, but this is Z memory with the paper's
X-only data faults, noisy Z records, ideal complementary X records and gates.
The unchanged Z-sector reference supplies the matching decoder and an
independent statistical cross-check. This is local Selene, not H2 deployment.

The scan selects the lowest mean patch failure probability. Fresh, fixed-size
samples then compare that interval against a baseline chosen before scanning.
Confirmation ratio bounds use Bonferroni-adjusted Clopper-Pearson intervals
across both probabilities and every distance/patch comparison. No ansatz is
used to generate, select or fit the measured logical-error probabilities.

Each invocation creates results/stage5/<UTC timestamp>/ with points.csv,
summary.json, logical_error_vs_interval.png and run.log. --quick tests the
entire workflow with low statistics; it is not a scientific validation run.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

# Direct script execution puts scripts/, not src/, on Python's import path.
# Resolve from this file so execution also works from an IDE or another cwd.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from qec_timing.circuit import build_surface_code_program
from qec_timing.noise import NoiseParams
from qec_timing.reference.circuit import build_memory_circuit
from qec_timing.reference.decoder import build_decoder
from qec_timing.reference.sampler import clopper_pearson, estimate_logical_error
from qec_timing.wire import unpack_packed_words

OUT_DIR = Path(__file__).resolve().parents[1] / "results" / "stage5"


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--distances", nargs="+", type=int, default=[5, 7])
    parser.add_argument("--q-d", type=int, default=1, help="independent logical patches")
    parser.add_argument("--rounds", nargs="+", type=int, default=[100, 50, 25, 20, 10, 5])
    parser.add_argument("--T", type=float, default=10.0, help="fixed storage duration")
    parser.add_argument("--p", type=float, default=0.015)
    parser.add_argument("--lam", type=float, default=1.0)
    parser.add_argument("--b-read", type=float, default=1.0)
    parser.add_argument("--shots", type=int, default=1000, help="full-circuit shots per scan point")
    parser.add_argument("--confirmation-shots", type=int, default=2000)
    parser.add_argument("--reference-shots", type=int, default=20_000)
    parser.add_argument("--baseline-rounds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20261002)
    parser.add_argument("--confidence", type=float, default=0.95)
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args(argv)
    if args.quick:
        args.distances, args.rounds, args.T = [3], [20, 10, 5], 2.0
        args.shots, args.confirmation_shots, args.reference_shots = 32, 64, 512
        args.baseline_rounds = 5
    args.distances = sorted(set(args.distances))
    args.rounds = sorted(set(args.rounds), reverse=True)
    if any(d < 3 or d % 2 == 0 for d in args.distances):
        parser.error("distances must be odd integers >=3")
    if len(args.rounds) < 3 or min(args.rounds) < 1:
        parser.error("provide at least three distinct positive round counts")
    if min(args.q_d, args.shots, args.confirmation_shots, args.reference_shots,
           args.baseline_rounds) < 1:
        parser.error("patch count, shot counts and baseline rounds must be positive")
    if not math.isfinite(args.T) or args.T <= 0:
        parser.error("T must be finite and positive")
    if not 0 < args.confidence < 1:
        parser.error("confidence must lie strictly between 0 and 1")
    if args.seed < 0:
        parser.error("seed must be non-negative")
    try:
        params = NoiseParams(args.p, args.lam, args.T / min(args.rounds), args.b_read)
        _ = params.p_data, params.p_read
        if not all(math.isfinite(x) for x in (args.p, args.lam, args.b_read)):
            raise ValueError("noise parameters must be finite")
    except ValueError as exc:
        parser.error(str(exc))
    return args


def ratio_summary(baseline, selected, confidence):
    """Conservative ratio interval; None means an unbounded upper endpoint."""
    b = baseline["failures"] / baseline["shots"]
    s = selected["failures"] / selected["shots"]
    bl, bh = clopper_pearson(baseline["failures"], baseline["shots"], confidence)
    sl, sh = clopper_pearson(selected["failures"], selected["shots"], confidence)
    lower, upper = bl / sh, bh / sl if sl > 0 else None
    return {
        "gamma": b / s if s > 0 else None,
        "gamma_low": lower, "gamma_high": upper,
        "status": "improved" if lower > 1 else
                  "worse" if upper is not None and upper < 1 else "inconclusive",
    }


def _row(d, patch, phase, path, n_rounds, args, shots, failures, seed,
         p_data, p_read, seconds, confidence):
    low, high = clopper_pearson(failures, shots, confidence)
    return dict(d=d, patch=patch, phase=phase, path=path, dt=args.T / n_rounds,
                n_rounds=n_rounds, T_actual=args.T, shots=shots, failures=failures,
                p_L=failures / shots, ci_low=low, ci_high=high,
                logical_error_per_time=failures / shots / args.T,
                p_data=p_data, p_read=p_read, base_seed=seed, seconds=seconds,
                interval_confidence=confidence)


def run_point(program, n_rounds, phase, shots, seed, ref_seed, args):
    dt = args.T / n_rounds
    params = NoiseParams(args.p, args.lam, dt, args.b_read)
    circuit = build_memory_circuit(program.layout, [params.p_data] * n_rounds,
                                   params.p_read, with_coords=False)
    decoder = build_decoder(circuit)
    start = time.perf_counter()
    result = program.run(p_data=params.p_data, p_read=params.p_read, dt=dt,
                         n_rounds=n_rounds, base_seed=seed, shots=shots)
    predictions = program.decode(result, n_rounds, decoder)
    observed = program.observables(result)
    if predictions.shape != (shots, args.q_d) or observed.shape != predictions.shape:
        raise RuntimeError("unexpected logical output shape")
    for shot in result.results:
        tags = shot.collate_tags()
        if len(tags["complementary_syn"]) != n_rounds * args.q_d:
            raise RuntimeError("missing complementary check measurements")
        if unpack_packed_words(tags["complementary_syn"], program.layout.n_checks).any():
            raise RuntimeError("ideal complementary stabilizers changed sign")
    elapsed = time.perf_counter() - start
    full = [_row(program.layout.d, patch, phase, "full_circuit", n_rounds, args,
                 shots, int(np.count_nonzero(predictions[:, patch] != observed[:, patch])),
                 seed, params.p_data, params.p_read, elapsed, args.confidence)
            for patch in range(args.q_d)]
    start = time.perf_counter()
    ref = estimate_logical_error(circuit, decoder, shots=args.reference_shots,
                                 seed=ref_seed, confidence=args.confidence)
    reference = _row(program.layout.d, -1, phase, "reference", n_rounds, args,
                     ref.shots, ref.failures, ref_seed, params.p_data, params.p_read,
                     time.perf_counter() - start, args.confidence)
    # As in Stage 4a: 99.9% interval overlap is a consistency diagnostic,
    # not a formal equivalence test. Store counts so stronger tests are possible.
    rl, rh = clopper_pearson(ref.failures, ref.shots, 0.999)
    for row in full:
        fl, fh = clopper_pearson(row["failures"], shots, 0.999)
        row["reference_intervals_overlap_999"] = bool(fl <= rh and rl <= fh)
    reference["reference_intervals_overlap_999"] = ""
    return full, reference


def save_points(rows, directory):
    with (directory / "points.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot_points(rows, args, directory):
    fig, axes = plt.subplots(1, len(args.distances), squeeze=False,
                             figsize=(6 * len(args.distances), 4.5))
    for ax, d in zip(axes[0], args.distances, strict=True):
        for path, patch in [("full_circuit", i) for i in range(args.q_d)] + [("reference", -1)]:
            points = sorted((r for r in rows if r["phase"] == "sweep" and
                             r["d"] == d and r["path"] == path and r["patch"] == patch),
                            key=lambda r: r["dt"])
            probabilities = np.array([r["p_L"] for r in points])
            errors = np.array([[r["p_L"] - r["ci_low"] for r in points],
                               [r["ci_high"] - r["p_L"] for r in points]])
            ax.errorbar([r["dt"] for r in points], probabilities / args.T,
                        yerr=errors / args.T, fmt="o-" if patch >= 0 else "s--",
                        capsize=3, label=f"full patch {patch}" if patch >= 0 else "reference")
        ax.axvline(args.T / args.baseline_rounds, color="grey", linestyle=":", label="baseline")
        ax.set_xscale("log")
        ax.set_xlabel("Syndrome interval Δt")
        ax.set_ylabel("Decoded logical failure probability / T")
        ax.set_title(f"d={d}, fixed T={args.T:g}")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle("Full X/Z extraction: timing sweep" + (" (smoke test)" if args.quick else ""))
    fig.tight_layout()
    fig.savefig(directory / "logical_error_vs_interval.png", dpi=150)
    plt.close(fig)


def main(argv=None):
    args = parse_args(argv)
    directory = OUT_DIR / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")
    directory.mkdir(parents=True, exist_ok=False)
    rows, confirmations = [], []
    comparisons = len(args.distances) * args.q_d
    marginal_confidence = 1 - (1 - args.confidence) / (2 * comparisons)
    stride = max(args.shots, args.confirmation_shots, args.reference_shots) + 1000
    seed_counter = 0
    started = time.perf_counter()
    payload = dict(configuration=vars(args), circuit="full X/Z extraction, Z memory",
                   noise_model="X-only phenomenological; ideal gates and complementary X records",
                   packages={p: version(p) for p in ["guppylang", "selene-sim", "stim", "pymatching"]},
                   confirmation_family_confidence=args.confidence,
                   confirmation_marginal_confidence=marginal_confidence,
                   status="running", confirmations=confirmations)

    def save_summary():
        payload["elapsed_seconds"] = time.perf_counter() - started
        (directory / "summary.json").write_text(
            json.dumps(payload, indent=2, allow_nan=False), encoding="utf-8")

    with (directory / "run.log").open("w", encoding="utf-8") as log:
        def report(message):
            print(message, flush=True)
            log.write(message + "\n")
            log.flush()

        def point(program, rounds, phase, shots):
            nonlocal seed_counter
            seed = args.seed + seed_counter * stride
            ref_seed = seed + stride
            seed_counter += 2
            full, ref = run_point(program, rounds, phase, shots, seed, ref_seed, args)
            rows.extend([*full, ref])
            save_points(rows, directory)
            report(f"d={program.layout.d} {phase} dt={args.T / rounds:.6g} rounds={rounds} "
                   f"full failures={[r['failures'] for r in full]}/{shots}, "
                   f"reference={ref['failures']}/{ref['shots']}")
            return full

        report(f"Results: {directory}")
        save_summary()
        try:
            for d in args.distances:
                program = build_surface_code_program(d=d, q_d=args.q_d, basis="Z", emit="packed")
                sweep = [(n, point(program, n, "sweep", args.shots)) for n in args.rounds]
                selected_rounds, _ = min(sweep, key=lambda entry: np.mean([r["p_L"] for r in entry[1]]))
                baseline = point(program, args.baseline_rounds, "confirmation_baseline",
                                 args.confirmation_shots)
                selected = point(program, selected_rounds, "confirmation_selected",
                                 args.confirmation_shots)
                for patch, (b, s) in enumerate(zip(baseline, selected, strict=True)):
                    conclusion = dict(d=d, patch=patch, baseline_dt=b["dt"], selected_dt=s["dt"],
                                      selected_is_scan_boundary=selected_rounds in
                                      (min(args.rounds), max(args.rounds)),
                                      baseline_failures=b["failures"], selected_failures=s["failures"],
                                      shots_each=args.confirmation_shots,
                                      **ratio_summary(b, s, marginal_confidence))
                    if selected_rounds == args.baseline_rounds:
                        conclusion["status"] = "baseline_selected"
                    confirmations.append(conclusion)
                    report(f"confirmation d={d} patch={patch}: {conclusion['status']}; "
                           f"gamma={conclusion['gamma']} interval="
                           f"[{conclusion['gamma_low']}, {conclusion['gamma_high']}]")
                save_summary()
            full_rows = [r for r in rows if r["path"] == "full_circuit"]
            payload["all_reference_intervals_overlap_999"] = all(
                r["reference_intervals_overlap_999"] for r in full_rows)
            payload["saturated_points"] = sum(r["p_L"] >= 0.45 for r in full_rows)
            improved = sum(c["status"] == "improved" for c in confirmations)
            payload["premise_result"] = "improved_all_patches" if improved == comparisons else \
                                        "improved_some_patches" if improved else "not_demonstrated"
            payload["interpretation"] = (
                "Low-statistics smoke test only; no scientific validation claim."
                if args.quick else
                "Fresh-shot confirmation for the specified phenomenological model; "
                "reference interval overlap is a diagnostic, not proof of equivalence."
            )
            payload["status"] = "completed"
            plot_points(rows, args, directory)
            save_summary()
            report(f"Premise: {payload['premise_result']}; "
                   f"reference consistency: {payload['all_reference_intervals_overlap_999']}")
        except Exception as exc:
            payload.update(status="failed", error=f"{type(exc).__name__}: {exc}")
            save_summary()
            report(payload["error"])
            raise
    return directory


if __name__ == "__main__":
    main()
