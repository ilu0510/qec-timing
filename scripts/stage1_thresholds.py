"""Validation run 1: reproduce the Appendix C code thresholds.

Protocol (Appendix C): set ``lam = 0``, run ``3d`` rounds of syndrome
measurement, and sweep the physical error rate ``p`` around the threshold. The
``p_L(p)`` curves for different ``d`` cross at ``p_th``. The crossing is located
by fitting the finite-size ansatz

    p_L = a0 + a1 x + a2 x**2,    x = (p - p_th) d**(1/nu)

jointly over all distances for each read-out coefficient ``b_read``.

Expected: p_th ~ 0.036 (b_read=0.5), 0.029 (b_read=1), 0.023 (b_read=2).

Usage:
    python scripts/stage1_thresholds.py [--shots N] [--quick]
"""

from __future__ import annotations

import argparse
import csv
import json
import time
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from qec_timing.reference import (
    RotatedSurfaceCodeZSector,
    build_decoder,
    build_memory_circuit_from_schedule,
    constant_schedule,
    estimate_logical_error,
)

OUT_DIR = Path(__file__).resolve().parents[1] / "results" / "stage1"
SEED = 20260927

# Sweep windows centred on the expected thresholds.
SWEEPS: dict[float, tuple[float, float]] = {
    0.5: (0.030, 0.042),
    1.0: (0.024, 0.034),
    2.0: (0.019, 0.027),
}
# Precise values from the Fig. 8 caption in Appendix C; docs/paper_equations.md
# quotes them rounded to 0.036 / 0.029 / 0.023.
EXPECTED: dict[float, float] = {0.5: 0.03603, 1.0: 0.02878, 2.0: 0.0231}
DISTANCES = (7, 11, 15, 19)
N_P_POINTS = 9

# The paper fits d = 7..35 (d = 15..45 for b_read = 2); we stop at 19, so the
# fit retains finite-size contamination from the smallest distances. Refitting
# on progressively larger-d subsets exposes that drift.
D_SUBSETS = ((7, 11, 15, 19), (11, 15, 19), (15, 19))


@dataclass
class Row:
    b_read: float
    d: int
    p: float
    n_rounds: int
    shots: int
    failures: int
    p_L: float
    ci_low: float
    ci_high: float
    stderr: float
    seconds: float


@dataclass
class FitResult:
    b_read: float
    p_th: float
    p_th_err: float
    nu: float
    nu_err: float
    expected: float
    coeffs: list[float] = field(default_factory=list)

    @property
    def deviation_sigma(self) -> float:
        """How far the fit sits from the paper value, in fitted sigmas."""
        if self.p_th_err <= 0:
            return float("inf")
        return abs(self.p_th - self.expected) / self.p_th_err


def collect(shots: int, distances: tuple[int, ...], n_points: int) -> list[Row]:
    rows: list[Row] = []
    for b_read, (p_lo, p_hi) in SWEEPS.items():
        p_values = np.linspace(p_lo, p_hi, n_points)
        for d in distances:
            layout = RotatedSurfaceCodeZSector(d)
            # lam = 0, so p_data = p exactly and dt is immaterial; 3d rounds.
            schedule = constant_schedule(dt=1.0, T=float(3 * d), lam=0.0)
            assert schedule.n_rounds == 3 * d
            for p in p_values:
                t0 = time.perf_counter()
                circuit = build_memory_circuit_from_schedule(
                    layout, schedule, p=float(p), b_read=b_read, with_coords=False
                )
                decoder = build_decoder(circuit)
                est = estimate_logical_error(
                    circuit,
                    decoder,
                    shots=shots,
                    seed=SEED + d * 1000 + int(p * 1e6),
                )
                elapsed = time.perf_counter() - t0
                rows.append(
                    Row(
                        b_read=b_read,
                        d=d,
                        p=float(p),
                        n_rounds=schedule.n_rounds,
                        shots=est.shots,
                        failures=est.failures,
                        p_L=est.p_L,
                        ci_low=est.ci_low,
                        ci_high=est.ci_high,
                        stderr=est.stderr,
                        seconds=elapsed,
                    )
                )
                print(
                    f"  b_read={b_read:<4} d={d:<3} p={p:.5f} "
                    f"p_L={est.p_L:.4f} +/- {est.stderr:.4f}  ({elapsed:.1f}s)",
                    flush=True,
                )
    return rows


def _model(stacked: np.ndarray, p_th: float, nu: float, a0: float, a1: float, a2: float):
    """Appendix C ansatz, evaluated on stacked (p, d) columns."""
    p, d = stacked[0], stacked[1]
    x = (p - p_th) * d ** (1.0 / nu)
    return a0 + a1 * x + a2 * x * x


def read_csv(path: Path) -> list[Row]:
    """Reload a completed sweep so fits can be redone without resampling."""
    rows: list[Row] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for rec in csv.DictReader(handle):
            rows.append(
                Row(
                    b_read=float(rec["b_read"]),
                    d=int(rec["d"]),
                    p=float(rec["p"]),
                    n_rounds=int(rec["n_rounds"]),
                    shots=int(rec["shots"]),
                    failures=int(rec["failures"]),
                    p_L=float(rec["p_L"]),
                    ci_low=float(rec["ci_low"]),
                    ci_high=float(rec["ci_high"]),
                    stderr=float(rec["stderr"]),
                    seconds=float(rec["seconds"]),
                )
            )
    return rows


def fit_threshold(
    rows: list[Row], b_read: float, distances: tuple[int, ...] | None = None
) -> FitResult:
    subset = [
        r
        for r in rows
        if r.b_read == b_read and (distances is None or r.d in distances)
    ]
    p = np.array([r.p for r in subset])
    d = np.array([float(r.d) for r in subset])
    y = np.array([r.p_L for r in subset])
    sigma = np.array([max(r.stderr, 1e-6) for r in subset])

    guess = [EXPECTED[b_read], 1.0, float(np.mean(y)), 1.0, 0.0]
    popt, pcov = curve_fit(
        _model,
        np.vstack([p, d]),
        y,
        p0=guess,
        sigma=sigma,
        absolute_sigma=True,
        maxfev=200_000,
    )
    err = np.sqrt(np.diag(pcov))
    return FitResult(
        b_read=b_read,
        p_th=float(popt[0]),
        p_th_err=float(err[0]),
        nu=float(popt[1]),
        nu_err=float(err[1]),
        expected=EXPECTED[b_read],
        coeffs=[float(v) for v in popt[2:]],
    )


def write_csv(rows: list[Row], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "b_read", "d", "p", "n_rounds", "shots", "failures",
                "p_L", "ci_low", "ci_high", "stderr", "seconds",
            ]
        )
        for r in rows:
            writer.writerow(
                [r.b_read, r.d, f"{r.p:.8g}", r.n_rounds, r.shots, r.failures,
                 f"{r.p_L:.8g}", f"{r.ci_low:.8g}", f"{r.ci_high:.8g}",
                 f"{r.stderr:.8g}", f"{r.seconds:.3f}"]
            )


def plot(rows: list[Row], fits: list[FitResult], path: Path) -> None:
    fig, axes = plt.subplots(1, len(fits), figsize=(5 * len(fits), 4.2), squeeze=False)
    for ax, fit in zip(axes[0], fits, strict=True):
        subset = [r for r in rows if r.b_read == fit.b_read]
        for d in sorted({r.d for r in subset}):
            pts = [r for r in subset if r.d == d]
            xs = np.array([r.p for r in pts])
            ys = np.array([r.p_L for r in pts])
            errs = np.array([r.stderr for r in pts])
            line = ax.errorbar(xs, ys, yerr=errs, marker="o", ms=4, lw=0,
                               elinewidth=1, capsize=2, label=f"d = {d}")
            dense = np.linspace(xs.min(), xs.max(), 200)
            ax.plot(
                dense,
                _model(np.vstack([dense, np.full_like(dense, float(d))]),
                       fit.p_th, fit.nu, *fit.coeffs),
                ls="--", lw=1, color=line[0].get_color(),
            )
        ax.axvline(fit.p_th, color="k", ls=":", lw=1)
        ax.set_title(
            f"$b_{{read}}$ = {fit.b_read}\n"
            f"fitted $p_{{th}}$ = {fit.p_th:.5f}, paper {fit.expected}"
        )
        ax.set_xlabel("physical error rate $p$")
        ax.set_ylabel("logical failure rate $p_L$")
        ax.legend(fontsize=8)
    fig.suptitle("Appendix C threshold fit ($\\lambda = 0$, $3d$ rounds)")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shots", type=int, default=40_000)
    parser.add_argument(
        "--quick", action="store_true", help="small run for smoke-testing"
    )
    parser.add_argument(
        "--refit-only",
        action="store_true",
        help="reload results/stage1/thresholds.csv and redo the fits only",
    )
    args = parser.parse_args()

    shots = 2_000 if args.quick else args.shots
    distances = (7, 11) if args.quick else DISTANCES
    n_points = 5 if args.quick else N_P_POINTS

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    suffix_in = "_quick" if args.quick else ""
    if args.refit_only:
        rows = read_csv(OUT_DIR / f"thresholds{suffix_in}.csv")
        elapsed = sum(r.seconds for r in rows)
        distances = tuple(sorted({r.d for r in rows}))
        print(f"Refitting {len(rows)} existing points, d={distances}")
    else:
        print(
            f"Appendix C thresholds: shots={shots}, d={distances}, "
            f"{n_points} p-values"
        )
        t0 = time.perf_counter()
        rows = collect(shots, distances, n_points)
        elapsed = time.perf_counter() - t0

    fits = [fit_threshold(rows, b) for b in SWEEPS]

    # Finite-size drift: refit on progressively larger-d subsets.
    drift: dict[str, dict[str, dict[str, float]]] = {}
    for b_read in SWEEPS:
        per_subset: dict[str, dict[str, float]] = {}
        for subset in D_SUBSETS:
            if not set(subset) <= set(distances):
                continue
            fit = fit_threshold(rows, b_read, subset)
            per_subset[str(subset)] = {
                "p_th": fit.p_th,
                "stderr": fit.p_th_err,
                "relative_deviation_percent": 100.0
                * (fit.p_th - fit.expected)
                / fit.expected,
            }
        drift[str(b_read)] = per_subset

    suffix = "_quick" if args.quick else ""
    write_csv(rows, OUT_DIR / f"thresholds{suffix}.csv")
    plot(rows, fits, OUT_DIR / f"thresholds{suffix}.png")

    summary = {
        "shots_per_point": shots,
        "distances": list(distances),
        "n_p_points": n_points,
        "total_seconds": elapsed,
        "fits": [
            {
                "b_read": f.b_read,
                "p_th_fitted": f.p_th,
                "p_th_stderr": f.p_th_err,
                "p_th_paper": f.expected,
                "deviation_sigma": f.deviation_sigma,
                "nu": f.nu,
                "nu_stderr": f.nu_err,
                "relative_deviation_percent": 100.0
                * (f.p_th - f.expected)
                / f.expected,
            }
            for f in fits
        ],
        "finite_size_drift": drift,
    }
    (OUT_DIR / f"thresholds{suffix}_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    print(f"\nTotal sampling+decoding time: {elapsed / 60:.1f} min")
    print(
        f"{'b_read':>7} {'fitted p_th':>20} {'paper':>9} {'rel dev':>9} "
        f"{'sigma':>7} {'nu':>16}"
    )
    for f in fits:
        rel = 100.0 * (f.p_th - f.expected) / f.expected
        print(
            f"{f.b_read:>7} {f.p_th:>13.5f}+/-{f.p_th_err:.5f} {f.expected:>9} "
            f"{rel:>+8.2f}% {f.deviation_sigma:>6.1f} "
            f"{f.nu:>11.3f}+/-{f.nu_err:.3f}"
        )

    print("\nFinite-size drift (the paper fits up to d=35; we stop at 19):")
    for b_read, per_subset in drift.items():
        if not per_subset:
            continue
        print(f"  b_read = {b_read}")
        for subset, values in per_subset.items():
            print(
                f"    d in {subset:<18} p_th = {values['p_th']:.5f}"
                f" +/- {values['stderr']:.5f}"
                f"  ({values['relative_deviation_percent']:+.2f}%)"
            )


if __name__ == "__main__":
    main()
