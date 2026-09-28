"""Stage 1 carry-over: locate the Fig. 6 minimum properly.

Stage 1 sampled dt on a coarse log grid (spacing factor 1.78) and could only say
that the observed minimum sat on one of the two grid points bracketing
``dt* = 0.714``. Here we refine the grid *inside* that bracket and fit a
parabola in ``log(dt)``, which uses every point rather than just the lowest one,
and gives the vertex an actual uncertainty.

Reference model only (no Selene). d = 5, b_read = 1, lam = 1, T = 100 -- the
same T assumption as Stage 1, documented there.

Usage:
    python scripts/stage2_refine_minimum.py [--budget-seconds N]
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
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

OUT_DIR = Path(__file__).resolve().parents[1] / "results" / "stage2"
SEED = 20260927

D, B_READ, LAM, T = 5, 1.0, 1.0, 100.0
TICK = 1e-5
G = 0.8

# Window brackets Stage 1's adjacent grid points (0.562 and 1.0) with room
# either side; spacing factor ~1.25 instead of Stage 1's 1.78.
DT_WINDOW = (0.3, 1.8)

# p -> (n_points, target_failures, max_shots). p = 0.001 is ~40x more expensive
# per failure than p = 0.0035, so it gets a smaller, capped budget.
PLAN = {
    0.0035: (9, 800, 400_000),
    0.001: (9, 800, 1_200_000),
}


def snap(dt: float) -> float:
    return round(dt / TICK) * TICK


def parabola(u, a, b, c):
    return a + b * u + c * u * u


def fit_vertex(dt: np.ndarray, p_L: np.ndarray, rel_err: np.ndarray):
    """Fit log(p_L) vs log(dt) with a parabola; return vertex and its error."""
    u = np.log(dt)
    y = np.log(p_L)
    popt, pcov = curve_fit(parabola, u, y, sigma=rel_err, absolute_sigma=True)
    a, b, c = popt
    if c <= 0:
        raise RuntimeError(f"fit is not convex (c = {c:.4g}); no interior minimum")

    u_star = -b / (2.0 * c)
    # Propagate covariance through u* = -b / (2c).
    grad = np.array([0.0, -1.0 / (2.0 * c), b / (2.0 * c * c)])
    var_u = float(grad @ pcov @ grad)
    sigma_u = math.sqrt(max(var_u, 0.0))

    dt_star = math.exp(u_star)
    return dt_star, dt_star * sigma_u, popt, pcov


def collect(p: float, n_points: int, target: int, max_shots: int):
    layout = RotatedSurfaceCodeZSector(D)
    grid = [snap(v) for v in np.geomspace(*DT_WINDOW, n_points)]
    rows = []
    for dt in grid:
        schedule = constant_schedule(dt=dt, T=T, lam=LAM, tick=TICK)
        circuit = build_memory_circuit_from_schedule(
            layout, schedule, p=p, b_read=B_READ, with_coords=False
        )
        decoder = build_decoder(circuit)
        t0 = time.perf_counter()
        est = estimate_logical_error_adaptive(
            circuit,
            decoder,
            seed=SEED + int(dt * 1e5),
            target_failures=target,
            max_shots=max_shots,
        )
        elapsed = time.perf_counter() - t0
        rows.append(
            {
                "p": p, "dt": dt, "n_rounds": schedule.n_rounds,
                "T_actual": schedule.total_time, "shots": est.shots,
                "failures": est.failures, "p_L": est.p_L,
                "stderr": est.stderr, "ci_low": est.ci_low,
                "ci_high": est.ci_high, "seconds": elapsed,
            }
        )
        print(
            f"  p={p:<7} dt={dt:<8.5f} rounds={schedule.n_rounds:<5} "
            f"p_L={est.p_L:.4e} +/-{100 * est.relative_stderr:.1f}%  "
            f"({est.failures}/{est.shots}, {elapsed:.1f}s)",
            flush=True,
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    dt_star_eq6 = ansatz.optimal_interval(d=D, lam=LAM, g=G)
    print(f"Refining the Fig. 6 minimum. Eq. 6 predicts dt* = {dt_star_eq6:.4f}")

    all_rows, results = [], {}
    t_start = time.perf_counter()
    for p, (n_points, target, max_shots) in PLAN.items():
        if args.quick:
            n_points, target, max_shots = 5, 60, 20_000
        rows = collect(p, n_points, target, max_shots)
        all_rows.extend(rows)

        dt = np.array([r["dt"] for r in rows])
        p_L = np.array([r["p_L"] for r in rows])
        rel = np.array([r["stderr"] / r["p_L"] for r in rows])
        try:
            vertex, vertex_err, popt, _ = fit_vertex(dt, p_L, rel)
            deviation = (vertex - dt_star_eq6) / dt_star_eq6
            sigmas = abs(vertex - dt_star_eq6) / vertex_err if vertex_err else float("inf")
        except RuntimeError as exc:
            vertex = vertex_err = deviation = sigmas = float("nan")
            popt = [float("nan")] * 3
            print(f"  fit failed for p={p}: {exc}")
        results[str(p)] = {
            "dt_min_fitted": vertex,
            "dt_min_stderr": vertex_err,
            "dt_star_eq6": dt_star_eq6,
            "relative_deviation": deviation,
            "deviation_sigma": sigmas,
            "parabola_coefficients": [float(v) for v in popt],
            "n_points": len(rows),
            "total_shots": sum(r["shots"] for r in rows),
        }
    elapsed = time.perf_counter() - t_start

    suffix = "_quick" if args.quick else ""
    with (OUT_DIR / f"refine_minimum{suffix}.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(all_rows[0]))
        writer.writeheader()
        writer.writerows(all_rows)
    (OUT_DIR / f"refine_minimum{suffix}_summary.json").write_text(
        json.dumps(
            {"dt_star_eq6": dt_star_eq6, "T": T, "d": D, "lam": LAM,
             "total_seconds": elapsed, "fits": results},
            indent=2,
        ),
        encoding="utf-8",
    )

    # Plot
    fig, ax = plt.subplots(figsize=(7, 4.6))
    dense = np.geomspace(*DT_WINDOW, 300)
    for (p, _), colour in zip(PLAN.items(), ("tab:blue", "tab:orange")):
        rows = [r for r in all_rows if r["p"] == p]
        dt = np.array([r["dt"] for r in rows])
        p_L = np.array([r["p_L"] for r in rows])
        err = np.array([r["stderr"] for r in rows])
        ax.errorbar(dt, p_L, yerr=err, fmt="o", ms=5, color=colour,
                    capsize=2, label=f"p = {p}")
        coeffs = results[str(p)]["parabola_coefficients"]
        if not any(math.isnan(c) for c in coeffs):
            ax.plot(dense, np.exp(parabola(np.log(dense), *coeffs)),
                    ls="--", lw=1.2, color=colour)
            v = results[str(p)]["dt_min_fitted"]
            ax.axvline(v, color=colour, ls="-", lw=0.9, alpha=0.6)
    ax.axvline(dt_star_eq6, color="k", ls=":", lw=1.4,
               label=f"$\\Delta t^*$ = {dt_star_eq6:.3f} (Eq. 6)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("syndrome interval $\\Delta t$")
    ax.set_ylabel("logical error probability $p_L$")
    ax.set_title(f"Refined minimum, d = {D}, T = {T} (parabola fit in log $\\Delta t$)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(OUT_DIR / f"refine_minimum{suffix}.png", dpi=150)
    plt.close(fig)

    print(f"\nTotal time {elapsed / 60:.1f} min")
    print(f"Eq. 6: dt* = {dt_star_eq6:.4f}")
    for p, info in results.items():
        print(
            f"  p={p:<8} fitted minimum dt = {info['dt_min_fitted']:.4f}"
            f" +/- {info['dt_min_stderr']:.4f}"
            f"  ({100 * info['relative_deviation']:+.1f}%, "
            f"{info['deviation_sigma']:.1f} sigma)"
        )


if __name__ == "__main__":
    main()
