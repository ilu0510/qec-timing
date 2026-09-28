"""Validation run 2: reproduce Fig. 6 and fit the ansatz prefactor A.

Fig. 6 of the paper: total logical error rate ``p_L`` after time ``T`` against
the syndrome interval ``dt``, at ``d = 5``, ``b_read = 1``, for several physical
error rates. The paper fixes ``beta = 2``, ``g = 0.8``, ``p_th = 0.029`` and
fits only ``A``, reporting ``A = 0.75``.

The paper does **not** state T for Fig. 6. We use ``T = 100``, the value given
for Fig. 5a which is also ``d = 5``. Override with ``--T``.

Caveat worth reading before trusting the fit: Fig. 6's vertical axis spans only
1e-3 to 1e-1, but at ``T = 100`` the ansatz puts the ``p = 0.01`` curve at
``p_L -> 0.5`` for small ``dt`` (many noisy rounds), i.e. saturated against
Eq. A3's ceiling and off the published axis. Saturated points carry no
information about ``A``, so they are excluded from the fit and reported
separately.

Usage:
    python scripts/stage1_fig6.py [--T 100] [--target-failures 200] [--quick]
"""

from __future__ import annotations

import argparse
import csv
import json
import time
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from qec_timing import ansatz
from qec_timing.reference import (
    RotatedSurfaceCodeZSector,
    build_decoder,
    build_memory_circuit_from_schedule,
    constant_schedule,
    estimate_logical_error_adaptive,
)

OUT_DIR = Path(__file__).resolve().parents[1] / "results" / "stage1"
SEED = 20260927

D = 5
B_READ = 1.0
LAM = 1.0
P_VALUES = (0.001, 0.0035, 0.01)
DT_RANGE = (0.01, 10.0)
N_DT = 13
TICK = 1e-5  # fine enough that snapping dt to the tick grid is <0.02%

BETA, G, P_TH = 2.0, 0.8, 0.029
A_PAPER = 0.75
# Above this p_L, Eq. A3 is effectively at its 0.5 ceiling and tells us nothing
# about A.
SATURATION_CUTOFF = 0.40


@dataclass
class Point:
    p: float
    dt_requested: float
    dt_actual: float
    n_rounds: int
    T_actual: float
    shots: int
    failures: int
    p_L: float
    ci_low: float
    ci_high: float
    stderr: float
    seconds: float

    @property
    def saturated(self) -> bool:
        return self.p_L > SATURATION_CUTOFF

    @property
    def usable(self) -> bool:
        return self.failures > 0 and not self.saturated


def snap_to_tick(dt: float, tick: float = TICK) -> float:
    """Snap an interval onto the tick grid (CLAUDE.md: exact integer time)."""
    return round(dt / tick) * tick


def collect(T: float, target_failures: int, n_dt: int, max_shots: int) -> list[Point]:
    layout = RotatedSurfaceCodeZSector(D)
    dt_grid = np.geomspace(*DT_RANGE, n_dt)
    points: list[Point] = []

    for p in P_VALUES:
        for dt_req in dt_grid:
            dt = snap_to_tick(float(dt_req))
            schedule = constant_schedule(dt=dt, T=T, lam=LAM, tick=TICK)
            t0 = time.perf_counter()
            circuit = build_memory_circuit_from_schedule(
                layout, schedule, p=p, b_read=B_READ, with_coords=False
            )
            decoder = build_decoder(circuit)
            est = estimate_logical_error_adaptive(
                circuit,
                decoder,
                seed=SEED + int(p * 1e6) * 100 + int(dt * 1e5),
                target_failures=target_failures,
                max_shots=max_shots,
            )
            elapsed = time.perf_counter() - t0
            points.append(
                Point(
                    p=p,
                    dt_requested=float(dt_req),
                    dt_actual=dt,
                    n_rounds=schedule.n_rounds,
                    T_actual=schedule.total_time,
                    shots=est.shots,
                    failures=est.failures,
                    p_L=est.p_L,
                    ci_low=est.ci_low,
                    ci_high=est.ci_high,
                    stderr=est.stderr,
                    seconds=elapsed,
                )
            )
            flag = " SATURATED" if points[-1].saturated else ""
            print(
                f"  p={p:<7} dt={dt:<9.5f} rounds={schedule.n_rounds:<6} "
                f"p_L={est.p_L:.3e} +/-{est.stderr:.1e} "
                f"({est.failures}/{est.shots}, {elapsed:.1f}s){flag}",
                flush=True,
            )
    return points


def fit_A(points: list[Point]) -> tuple[float, float, int]:
    """Fit the single prefactor A on log p_L over unsaturated points."""
    usable = [pt for pt in points if pt.usable]
    if len(usable) < 2:
        raise RuntimeError("not enough unsaturated points to fit A")

    dt = np.array([pt.dt_actual for pt in usable])
    p = np.array([pt.p for pt in usable])
    n_rounds = np.array([pt.n_rounds for pt in usable])
    y = np.log(np.array([pt.p_L for pt in usable]))
    # d(log p_L) = stderr / p_L
    sigma = np.array([pt.stderr / pt.p_L for pt in usable])

    def model(_x, A):
        return np.log(
            [
                ansatz.p_L_constant_interval(
                    float(dt_i), n_rounds=int(n_i), d=D, p=float(p_i),
                    lam=LAM, A=A, p_th=P_TH, beta=BETA, g=G,
                )
                for dt_i, n_i, p_i in zip(dt, n_rounds, p, strict=True)
            ]
        )

    popt, pcov = curve_fit(
        model, np.zeros(len(usable)), y, p0=[A_PAPER],
        sigma=sigma, absolute_sigma=True, maxfev=100_000,
    )
    return float(popt[0]), float(np.sqrt(pcov[0, 0])), len(usable)


def fit_A_per_p(points: list[Point]) -> dict[float, tuple[float, float, int]]:
    out: dict[float, tuple[float, float, int]] = {}
    for p in P_VALUES:
        subset = [pt for pt in points if pt.p == p]
        try:
            out[p] = fit_A(subset)
        except RuntimeError:
            out[p] = (float("nan"), float("nan"), 0)
    return out


def observed_minimum(points: list[Point], p: float) -> tuple[float, float]:
    """dt of the lowest measured p_L for one p, and that p_L."""
    subset = [pt for pt in points if pt.p == p and pt.failures > 0]
    best = min(subset, key=lambda pt: pt.p_L)
    return best.dt_actual, best.p_L


def write_csv(points: list[Point], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["p", "dt_requested", "dt_actual", "n_rounds", "T_actual", "shots",
             "failures", "p_L", "ci_low", "ci_high", "stderr", "saturated", "seconds"]
        )
        for pt in points:
            writer.writerow(
                [pt.p, f"{pt.dt_requested:.8g}", f"{pt.dt_actual:.8g}", pt.n_rounds,
                 f"{pt.T_actual:.8g}", pt.shots, pt.failures, f"{pt.p_L:.8g}",
                 f"{pt.ci_low:.8g}", f"{pt.ci_high:.8g}", f"{pt.stderr:.8g}",
                 int(pt.saturated), f"{pt.seconds:.3f}"]
            )


def plot(points: list[Point], A: float, T: float, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    dense = np.geomspace(*DT_RANGE, 300)

    for p, colour in zip(P_VALUES, ("tab:blue", "tab:orange", "tab:green"), strict=True):
        subset = [pt for pt in points if pt.p == p]
        xs = np.array([pt.dt_actual for pt in subset])
        ys = np.array([pt.p_L for pt in subset])
        errs = np.array([pt.stderr for pt in subset])
        sat = np.array([pt.saturated for pt in subset])

        ax.errorbar(xs[~sat], ys[~sat], yerr=errs[~sat], fmt="o", ms=5, color=colour,
                    capsize=2, label=f"p = {p}")
        if sat.any():
            ax.errorbar(xs[sat], ys[sat], yerr=errs[sat], fmt="s", ms=5, mfc="none",
                        color=colour, capsize=2)

        curve = [
            ansatz.p_L_constant_interval(
                float(dt), n_rounds=max(round(T / float(dt)), 1), d=D, p=p,
                lam=LAM, A=A, p_th=P_TH, beta=BETA, g=G,
            )
            for dt in dense
        ]
        ax.plot(dense, curve, ls="--", lw=1.2, color=colour)

    dt_star = ansatz.optimal_interval(d=D, lam=LAM, g=G)
    ax.axvline(dt_star, color="k", ls=":", lw=1,
               label=f"$\\Delta t^*$ = {dt_star:.3f} (Eq. 6)")
    ax.axhline(SATURATION_CUTOFF, color="grey", lw=0.8, alpha=0.5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("syndrome interval $\\Delta t$")
    ax.set_ylabel("logical error probability $p_L$")
    ax.set_title(
        f"Fig. 6 reproduction: d = {D}, $b_{{read}}$ = {B_READ}, "
        f"$\\lambda$ = {LAM}, T = {T}\n"
        f"fitted A = {A:.3f} (paper 0.75); open squares are saturated, excluded"
    )
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--T", type=float, default=100.0)
    parser.add_argument("--target-failures", type=int, default=200)
    parser.add_argument("--max-shots", type=int, default=400_000)
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()

    n_dt = 5 if args.quick else N_DT
    target = 40 if args.quick else args.target_failures
    max_shots = 20_000 if args.quick else args.max_shots

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(
        f"Fig. 6: d={D}, b_read={B_READ}, lam={LAM}, T={args.T}, "
        f"{n_dt} dt points, target {target} failures"
    )
    t0 = time.perf_counter()
    points = collect(args.T, target, n_dt, max_shots)
    elapsed = time.perf_counter() - t0

    A, A_err, n_used = fit_A(points)
    per_p = fit_A_per_p(points)
    dt_star = ansatz.optimal_interval(d=D, lam=LAM, g=G)

    suffix = "_quick" if args.quick else ""
    write_csv(points, OUT_DIR / f"fig6{suffix}.csv")
    plot(points, A, args.T, OUT_DIR / f"fig6{suffix}.png")

    minima = {p: observed_minimum(points, p) for p in P_VALUES}
    summary = {
        "d": D, "b_read": B_READ, "lam": LAM, "T_nominal": args.T,
        "beta": BETA, "g": G, "p_th": P_TH,
        "A_fitted": A, "A_stderr": A_err, "A_paper": A_PAPER,
        "n_points_used": n_used,
        "n_points_saturated": sum(pt.saturated for pt in points),
        "A_per_p": {str(p): {"A": v[0], "stderr": v[1], "n_points": v[2]}
                     for p, v in per_p.items()},
        "dt_star_eq6": dt_star,
        "observed_minima": {str(p): {"dt": v[0], "p_L": v[1]}
                             for p, v in minima.items()},
        "total_seconds": elapsed,
        "total_shots": sum(pt.shots for pt in points),
    }
    (OUT_DIR / f"fig6{suffix}_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    print(f"\nTotal time: {elapsed / 60:.1f} min, {summary['total_shots']:,} shots")
    print(f"Fitted A = {A:.4f} +/- {A_err:.4f}   (paper 0.75), "
          f"using {n_used}/{len(points)} points "
          f"({summary['n_points_saturated']} saturated, excluded)")
    for p, (a, err, n) in per_p.items():
        print(f"   p={p:<7} A = {a:.4f} +/- {err:.4f}  ({n} points)")
    print(f"\nEq. 6 optimum dt* = {dt_star:.4f}")
    for p, (dt_min, p_L_min) in minima.items():
        print(f"   p={p:<7} observed minimum at dt = {dt_min:.4f} "
              f"(p_L = {p_L_min:.3e})")


if __name__ == "__main__":
    main()
