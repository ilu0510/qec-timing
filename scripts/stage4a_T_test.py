"""Discriminate two explanations for dt* sitting above Eq. 6.

At T = 10 the fitted minima sat 17-31% above Eq. 6's dt*. Two candidate causes:

(a) **Effective g < 0.8 at small d.** Eq. 6 is `dt* = 1/(lambda(alpha - 1))`
    with `alpha = g(d+1)/2`; the paper fits g = 0.8 against d = 11-27. If g is
    really ~0.70 at d = 5-9, dt* is genuinely larger and the location is a
    property of the bulk physics.

(b) **A finite-rounds boundary effect.** At T = 10 the large-dt points have very
    few rounds -- d = 9 at dt = 1.33 runs only 8 rounds, fewer than d -- so the
    final-readout detector block is a large fraction of the history and the
    paper's T/dt = O(d) condition is violated. That could push the apparent
    minimum to larger dt.

**The discriminator:** in the ansatz, dt* does not depend on T at all. R(dt) is
T-independent, and `p_L = 1/2(1 - exp(-2 R T))` is monotone in `R T`, so
`argmin_dt p_L = argmin_dt R` for any fixed T. So tripling T changes the round
count at every dt by 3x while leaving the predicted optimum untouched:

* if the minima **move toward Eq. 6**, the T = 10 offset was cause (b);
* if they **stay put**, it is cause (a).

Cost is roughly T-independent: shots scale as 1/p_L ~ 1/(R T) while rounds scale
as T/dt, so shots x rounds ~ 1/(R dt). Reference path only.

Usage:
    python scripts/stage4a_T_test.py
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

import stage4a_validate as S
from qec_timing import ansatz
from qec_timing.layout import RotatedSurfaceCodeZSector

OUT_DIR = Path(__file__).resolve().parents[1] / "results" / "stage4a"

T_NEW = 30.0
DISTANCES = (7, 9)
N_POINTS = 9
TARGET_FAILURES = 800
MAX_SHOTS = 800_000

# Fitted dt* from the T = 10 run, used only to place the scan window.
DT_AT_T10 = {7: 0.5546, 9: 0.3888}
DT_ERR_AT_T10 = {7: 0.0114, 9: 0.0109}


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    S.T = T_NEW  # reference_point reads this global
    t_start = time.perf_counter()
    results = {}

    print(f"Repeating the dt* scan at T = {T_NEW} (reference path only)")
    for d in DISTANCES:
        layout = RotatedSurfaceCodeZSector(d)
        dt_eq6 = ansatz.optimal_interval(d=d, lam=S.LAM, g=S.G)
        centre = DT_AT_T10[d]
        # Window is wide enough to contain both the T=10 minimum and Eq. 6's.
        grid = np.geomspace(centre / 2.5, centre * 2.5, N_POINTS)
        assert grid[0] < dt_eq6 < grid[-1], "window must bracket Eq. 6's dt*"

        points = []
        for dt in grid:
            pt = S.reference_point(
                layout, float(dt), target_failures=TARGET_FAILURES,
                max_shots=MAX_SHOTS, seed_offset=31_000,
            )
            points.append(pt)
            print(
                f"  d={d} dt={pt['dt']:<8.4f} rounds={pt['n_rounds']:<5} "
                f"p_L={pt['p_L']:.4f} "
                f"({pt['failures']}/{pt['shots']}, {pt['seconds']:.0f}s)",
                flush=True,
            )

        dt_fit, dt_err, coeffs, n_used = S.fit_vertex(points)
        moved = dt_fit - DT_AT_T10[d]
        # Did it move toward Eq. 6, and by how much of the original gap?
        gap_before = DT_AT_T10[d] - dt_eq6
        closed = (gap_before - (dt_fit - dt_eq6)) / gap_before
        sigma_move = abs(moved) / np.hypot(dt_err, DT_ERR_AT_T10[d])

        results[str(d)] = {
            "T": T_NEW, "dt_star_fitted": dt_fit, "dt_star_stderr": dt_err,
            "dt_star_at_T10": DT_AT_T10[d], "dt_star_eq6": dt_eq6,
            "deviation_at_T10_percent": 100 * gap_before / dt_eq6,
            "deviation_at_T30_percent": 100 * (dt_fit - dt_eq6) / dt_eq6,
            "fraction_of_gap_closed": closed,
            "shift_sigma": float(sigma_move),
            "n_rounds_at_min": max(round(T_NEW / dt_fit), 1),
            "points_used": n_used, "parabola_coefficients": coeffs,
            "points": points,
        }
        print(
            f"  -> d={d}: dt* = {dt_fit:.4f} +/- {dt_err:.4f}  "
            f"(T=10 gave {DT_AT_T10[d]:.4f}, Eq. 6 {dt_eq6:.4f})\n"
            f"     deviation {100 * gap_before / dt_eq6:+.1f}% -> "
            f"{100 * (dt_fit - dt_eq6) / dt_eq6:+.1f}%, "
            f"gap closed {100 * closed:+.0f}%, shift {sigma_move:.1f} sigma\n",
            flush=True,
        )

    elapsed = time.perf_counter() - t_start
    closed = [r["fraction_of_gap_closed"] for r in results.values()]
    shifts = [r["shift_sigma"] for r in results.values()]
    verdict = (
        "(b) finite-rounds boundary effect"
        if min(closed) > 0.5 and min(shifts) > 3
        else "(a) effective g < 0.8 at small d"
        if max(abs(c) for c in closed) < 0.25 and max(shifts) < 3
        else "inconclusive / mixed"
    )
    payload = {"T": T_NEW, "verdict": verdict, "per_distance": results,
               "total_seconds": elapsed}
    (OUT_DIR / "T_test.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    print(f"=== T-test done in {elapsed / 60:.1f} min ===")
    print(f"{'d':>3} {'dt*(T=10)':>10} {'dt*(T=30)':>18} {'Eq. 6':>8} "
          f"{'gap closed':>11} {'shift':>8}")
    for d, r in results.items():
        print(f"{d:>3} {r['dt_star_at_T10']:>10.4f} "
              f"{r['dt_star_fitted']:>11.4f}+/-{r['dt_star_stderr']:.4f} "
              f"{r['dt_star_eq6']:>8.4f} {100 * r['fraction_of_gap_closed']:>10.0f}% "
              f"{r['shift_sigma']:>6.1f}s")
    print(f"\nVERDICT: {verdict}")


if __name__ == "__main__":
    main()
