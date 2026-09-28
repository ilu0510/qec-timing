"""Stage 4a: validate the noise model. Cheap by design.

The goal is a **validated noise model**, not a reproduction of Fig. 2. Success,
in priority order:

1. Guppy and the reference agree within error bars across the dt range.
2. Each d shows a U-shaped curve with a minimum near Eq. 6's dt*.
3. 1/dt* is roughly linear in d, with real error bars on the slope.

Design choices that keep it cheap:

* **d = 5, 7, 9.** d = 11 was costed at ~74 min on the Guppy path on its own and
  dropped; d = 5..9 together run in ~25 min. See docs/performance.md.
* **T = 10 for every d.** Chosen so p_L stays unsaturated (~0.03-0.3) at every
  sampled dt, and so the round count at the optimum is ~3d -- matching Appendix
  C's convention and the paper's requirement that T/dt = O(d). Saturation would
  compress the U and bias the fit.
* **~100 logical failures per point.** ~10% error bars, which is ample to see
  agreement and a U-shape.
* **Shots are sized up front from the ansatz, and each (d, dt) point is a single
  ``run()`` call.** The Guppy path pays ~9.3 s *per call* regardless of shot
  count (docs/performance.md), so an adaptive sample-until-N-failures loop would
  waste ~9 s per extra batch. Adaptive sampling is used on the reference path,
  where per-call cost is negligible.
* **Guppy runs the coarse grid only** (criterion 1). The refined grid used to
  locate dt* runs on the reference path, which is ~200x faster -- the two-path
  division of labour in CLAUDE.md.

Minima are located by a **parabola fit in log(dt)**, never by taking the lowest
sampled point (Stage 2 lesson). The parabola is fitted to **log R**, not
log p_L: Eq. A3 is inverted first, ``R = -ln(1 - 2 p_L) / (2 T_actual)``, which
removes the ceiling distortion and makes the fitted quantity the one Eq. 6
actually minimises.

Usage:
    python scripts/stage4a_validate.py [--quick] [--skip-guppy]
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
from qec_timing.layout import RotatedSurfaceCodeZSector
from qec_timing.noise import NoiseParams
from qec_timing.reference.circuit import build_memory_circuit
from qec_timing.reference.decoder import build_decoder
from qec_timing.reference.sampler import (
    clopper_pearson,
    estimate_logical_error_adaptive,
)
from qec_timing.schedule import constant_schedule
from qec_timing.stub_circuit import build_memory_program

OUT_DIR = Path(__file__).resolve().parents[1] / "results" / "stage4a"
BASE_SEED = 20260927

P, LAM, B_READ = 0.015, 1.0, 1.0
T = 10.0
P_TH, BETA, G = 0.029, 2.0, 0.8
TICK = 1e-5

DISTANCES = (5, 7, 9)
A_GUESS = {5: 0.75, 7: 1.0, 9: 1.2, 11: 1.4}

N_COARSE = 6
N_REFINED = 9
TARGET_FAILURES = 100
MIN_SHOTS, MAX_SHOTS = 150, 20_000

# Guppy cost is unstable to ~3x between runs (docs/performance.md), so the shot
# count is additionally capped by a per-point TIME budget. That bounds the total
# runtime by construction; points where the cap binds simply get wider error
# bars, which is acceptable here.
GUPPY_SECONDS_PER_POINT = 150.0
GUPPY_FIXED_S = 9.3
# Effective ms/round/shot, calibrated from an observed d=5 run and scaled by
# (n_data + n_checks). Deliberately pessimistic.
GUPPY_MS_PER_ROUND = {d: 1.27 * (d * d + (d * d - 1) // 2) / 37.0
                      for d in (5, 7, 9, 11)}
CONFIDENCE = 0.999
P_L_SATURATED = 0.45


def snap(dt: float) -> float:
    return round(dt / TICK) * TICK


def coarse_grid(d: int) -> list[float]:
    """Log-spaced dt centred on Eq. 6's optimum for this d."""
    dt_star = ansatz.optimal_interval(d=d, lam=LAM, g=G)
    return [snap(x) for x in np.geomspace(dt_star / 5.0, dt_star * 4.0, N_COARSE)]


def predicted_p_L(d: int, dt: float) -> float:
    n_rounds = max(round(T / dt), 1)
    return ansatz.p_L_constant_interval(
        dt, n_rounds=n_rounds, d=d, p=P, lam=LAM, A=A_GUESS[d],
        p_th=P_TH, beta=BETA, g=G,
    )


def planned_shots(d: int, dt: float) -> tuple[int, bool]:
    """Shots sized up front, so each point is exactly one run() call.

    Returns (shots, capped_by_time).
    """
    p_L = predicted_p_L(d, dt)
    wanted = TARGET_FAILURES / max(p_L, 1e-9)
    n_rounds = max(round(T / dt), 1)
    budget = (GUPPY_SECONDS_PER_POINT - GUPPY_FIXED_S) / (
        n_rounds * GUPPY_MS_PER_ROUND[d] / 1000.0
    )
    shots = int(min(max(min(wanted, budget), MIN_SHOTS), MAX_SHOTS))
    return shots, budget < wanted


def rate_from_p_L(p_L: float, stderr: float, T_actual: float):
    """Invert Eq. A3: p_L = 1/2 (1 - exp(-2 R T)) -> R, with error on log R."""
    if not 0.0 < p_L < 0.5:
        return float("nan"), float("nan")
    R = -math.log1p(-2.0 * p_L) / (2.0 * T_actual)
    sigma_log_R = stderr / ((1.0 - 2.0 * p_L) * T_actual * R)
    return R, sigma_log_R


# ---------------------------------------------------------------------------
# the two paths
# ---------------------------------------------------------------------------


def reference_point(layout, dt, *, target_failures, max_shots, seed_offset=0):
    schedule = constant_schedule(dt=snap(dt), T=T, lam=LAM, tick=TICK)
    params = NoiseParams(p=P, lam=LAM, dt=schedule.dt, b_read=B_READ)
    circuit = build_memory_circuit(
        layout, [params.p_data] * schedule.n_rounds, params.p_read,
        with_coords=False,
    )
    decoder = build_decoder(circuit)
    t0 = time.perf_counter()
    est = estimate_logical_error_adaptive(
        circuit, decoder,
        seed=BASE_SEED + seed_offset + int(schedule.dt * 1e5),
        target_failures=target_failures, max_shots=max_shots,
    )
    return {
        "dt": schedule.dt, "n_rounds": schedule.n_rounds,
        "T_actual": schedule.total_time, "shots": est.shots,
        "failures": est.failures, "p_L": est.p_L, "stderr": est.stderr,
        "ci_low": est.ci_low, "ci_high": est.ci_high,
        "seconds": time.perf_counter() - t0,
    }


def guppy_point(program, layout, dt, shots):
    """One (d, dt) point in a single run() call -- see docs/performance.md."""
    schedule = constant_schedule(dt=snap(dt), T=T, lam=LAM, tick=TICK)
    params = NoiseParams(p=P, lam=LAM, dt=schedule.dt, b_read=B_READ)
    circuit = build_memory_circuit(
        layout, [params.p_data] * schedule.n_rounds, params.p_read,
        with_coords=False,
    )
    decoder = build_decoder(circuit)

    t0 = time.perf_counter()
    result = program.run(
        p_data=params.p_data, p_read=params.p_read, dt=schedule.dt,
        n_rounds=schedule.n_rounds, base_seed=BASE_SEED, shots=shots,
    )
    tags = [s.collate_tags() for s in result.results]
    events = program.detection_events(result, schedule.n_rounds)
    predictions = decoder.decode_batch(np.stack([e.ravel() for e in events]))
    observed = np.array([[bool(t["obs"][0])] for t in tags])
    failures = int(np.count_nonzero(predictions != observed))
    elapsed = time.perf_counter() - t0

    p_L = failures / shots
    low, high = clopper_pearson(failures, shots, CONFIDENCE)
    return {
        "dt": schedule.dt, "n_rounds": schedule.n_rounds,
        "T_actual": schedule.total_time, "shots": shots, "failures": failures,
        "p_L": p_L,
        "stderr": math.sqrt(max(p_L * (1 - p_L), 0.0) / shots),
        "ci_low": low, "ci_high": high, "seconds": elapsed,
    }


# ---------------------------------------------------------------------------
# fits
# ---------------------------------------------------------------------------


def parabola(u, a, b, c):
    return a + b * u + c * u * u


def fit_vertex(points):
    u, y, sigma = [], [], []
    for p in points:
        if not (0.0 < p["p_L"] < P_L_SATURATED) or p["failures"] == 0:
            continue
        R, s = rate_from_p_L(p["p_L"], p["stderr"], p["T_actual"])
        if math.isfinite(R) and math.isfinite(s) and s > 0:
            u.append(math.log(p["dt"]))
            y.append(math.log(R))
            sigma.append(s)
    if len(u) < 4:
        raise RuntimeError(f"only {len(u)} usable points")
    popt, pcov = curve_fit(
        parabola, np.array(u), np.array(y), sigma=np.array(sigma),
        absolute_sigma=True,
    )
    a, b, c = popt
    if c <= 0:
        raise RuntimeError(f"not convex (c={c:.4g})")
    u_star = -b / (2.0 * c)
    grad = np.array([0.0, -1.0 / (2.0 * c), b / (2.0 * c * c)])
    sigma_u = math.sqrt(max(float(grad @ pcov @ grad), 0.0))
    dt_star = math.exp(u_star)
    return dt_star, dt_star * sigma_u, [float(v) for v in popt], len(u)


def fit_line(ds, ys, errs):
    """1/dt* = eta d + z, weighted, with errors on both coefficients."""
    def line(x, eta, z):
        return eta * x + z
    popt, pcov = curve_fit(
        line, np.array(ds, float), np.array(ys), sigma=np.array(errs),
        absolute_sigma=True,
    )
    err = np.sqrt(np.diag(pcov))
    return float(popt[0]), float(err[0]), float(popt[1]), float(err[1])


# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--skip-guppy", action="store_true")
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    distances = (5,) if args.quick else DISTANCES
    t_start = time.perf_counter()
    cross, per_d = [], {}

    for d in distances:
        layout = RotatedSurfaceCodeZSector(d)
        grid = coarse_grid(d)[: 3 if args.quick else None]
        dt_eq6 = ansatz.optimal_interval(d=d, lam=LAM, g=G)
        print(f"\nd = {d}  (Eq. 6 dt* = {dt_eq6:.4f}, T = {T})", flush=True)

        program = None if args.skip_guppy else build_memory_program(layout, emit="packed")
        coarse_ref, coarse_guppy = [], []

        for dt in grid:
            ref = reference_point(
                layout, dt, target_failures=TARGET_FAILURES * 4,
                max_shots=MAX_SHOTS * 20,
            )
            coarse_ref.append(ref)
            if program is None:
                print(f"  dt={dt:<8.4f} rounds={ref['n_rounds']:<5} "
                      f"ref p_L={ref['p_L']:.4f}", flush=True)
                continue

            shots, capped = (400, False) if args.quick else planned_shots(d, dt)
            gup = guppy_point(program, layout, dt, shots)
            coarse_guppy.append(gup)
            overlap = (gup["ci_low"] <= ref["ci_high"]) and (
                ref["ci_low"] <= gup["ci_high"]
            )
            rel = (gup["p_L"] - ref["p_L"]) / ref["p_L"] if ref["p_L"] else float("nan")
            cross.append({
                "d": d, "dt": dt, "n_rounds": gup["n_rounds"],
                "guppy_p_L": gup["p_L"], "guppy_ci_low": gup["ci_low"],
                "guppy_ci_high": gup["ci_high"], "guppy_shots": gup["shots"],
                "guppy_failures": gup["failures"], "guppy_seconds": gup["seconds"],
                "reference_p_L": ref["p_L"], "reference_ci_low": ref["ci_low"],
                "reference_ci_high": ref["ci_high"], "reference_shots": ref["shots"],
                "reference_failures": ref["failures"],
                "relative_difference": rel, "intervals_overlap": bool(overlap),
                "shots_time_capped": bool(capped),
            })
            print(
                f"  dt={dt:<8.4f} rounds={gup['n_rounds']:<5} "
                f"guppy {gup['p_L']:.4f} ({gup['failures']}/{gup['shots']}) "
                f"ref {ref['p_L']:.4f} ({ref['failures']}/{ref['shots']})  "
                f"{'OK' if overlap else 'MISMATCH'} {100 * rel:+.1f}%  "
                f"[{gup['seconds']:.0f}s{'*' if capped else ''}]",
                flush=True,
            )

        # Refinement on the reference path only.
        best = min((p for p in coarse_ref if p["failures"] > 0), key=lambda p: p["p_L"])
        lo, hi = best["dt"] / 2.5, best["dt"] * 2.5
        refined = []
        for dt in np.geomspace(lo, hi, 4 if args.quick else N_REFINED):
            pt = reference_point(
                layout, float(dt), target_failures=TARGET_FAILURES * 8,
                max_shots=MAX_SHOTS * 40, seed_offset=7_000,
            )
            refined.append(pt)
        print(f"  refined {len(refined)} reference points around dt={best['dt']:.4f}",
              flush=True)

        try:
            dt_fit, dt_err, coeffs, n_used = fit_vertex(refined)
        except RuntimeError as exc:
            print(f"  vertex fit failed: {exc}")
            dt_fit = dt_err = float("nan")
            coeffs, n_used = [float("nan")] * 3, 0

        per_d[str(d)] = {
            "dt_star_fitted": dt_fit, "dt_star_stderr": dt_err,
            "dt_star_eq6": dt_eq6,
            "relative_deviation": (dt_fit - dt_eq6) / dt_eq6,
            "inv_dt_star": 1.0 / dt_fit if dt_fit else float("nan"),
            "inv_dt_star_err": dt_err / dt_fit**2 if dt_fit else float("nan"),
            "parabola_coefficients": coeffs, "vertex_points_used": n_used,
            "coarse_reference": coarse_ref, "coarse_guppy": coarse_guppy,
            "refined_reference": refined,
        }
        print(f"  dt* = {dt_fit:.4f} +/- {dt_err:.4f}  "
              f"(Eq. 6 {dt_eq6:.4f}, {100 * per_d[str(d)]['relative_deviation']:+.1f}%)",
              flush=True)

    # 1/dt* vs d
    usable = [(int(d), i["inv_dt_star"], i["inv_dt_star_err"])
              for d, i in per_d.items() if math.isfinite(i["inv_dt_star"])]
    slope = None
    if len(usable) >= 2:
        eta, eta_err, z, z_err = fit_line(*zip(*usable, strict=True))
        slope = {"eta": eta, "eta_stderr": eta_err, "z": z, "z_stderr": z_err,
                 "paper_eta": 0.402}

    elapsed = time.perf_counter() - t_start
    n_bad = sum(not c["intervals_overlap"] for c in cross)
    payload = {
        "targets": {"p": P, "lam": LAM, "b_read": B_READ, "T": T,
                    "p_th": P_TH, "beta": BETA, "g": G},
        "distances": list(distances), "target_failures": TARGET_FAILURES,
        "cross_check_points": len(cross), "cross_check_disagreements": n_bad,
        "per_distance": per_d, "inv_dt_star_fit": slope,
        "total_seconds": elapsed,
    }
    suffix = "_quick" if args.quick else ""
    (OUT_DIR / f"stage4a{suffix}.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    if cross:
        with (OUT_DIR / f"stage4a_crosscheck{suffix}.csv").open(
            "w", newline="", encoding="utf-8"
        ) as h:
            w = csv.DictWriter(h, fieldnames=list(cross[0]))
            w.writeheader()
            w.writerows(cross)
    flat = [{"d": d, "grid": k, "path": p, **pt}
            for d, i in per_d.items()
            for k, p, key in (("coarse", "reference", "coarse_reference"),
                              ("coarse", "guppy", "coarse_guppy"),
                              ("refined", "reference", "refined_reference"))
            for pt in i[key]]
    if flat:
        with (OUT_DIR / f"stage4a_points{suffix}.csv").open(
            "w", newline="", encoding="utf-8"
        ) as h:
            w = csv.DictWriter(h, fieldnames=list(flat[0]))
            w.writeheader()
            w.writerows(flat)

    _plot(per_d, slope, OUT_DIR / f"stage4a{suffix}.png")

    print(f"\n=== Stage 4a, {elapsed / 60:.1f} min ===")
    print(f"1. Agreement: {len(cross) - n_bad}/{len(cross)} points overlap at "
          f"{CONFIDENCE:.1%}")
    print(f"{'d':>3} {'dt* fitted':>20} {'Eq. 6':>8} {'dev':>8} {'1/dt*':>8}")
    for d, i in per_d.items():
        print(f"{d:>3} {i['dt_star_fitted']:>13.4f}+/-{i['dt_star_stderr']:.4f} "
              f"{i['dt_star_eq6']:>8.4f} {100 * i['relative_deviation']:>+7.1f}% "
              f"{i['inv_dt_star']:>8.3f}")
    if slope:
        print(f"3. 1/dt* = ({slope['eta']:.3f} +/- {slope['eta_stderr']:.3f}) d "
              f"+ ({slope['z']:.3f} +/- {slope['z_stderr']:.3f})   "
              f"[paper eta = 0.402]")


def _plot(per_d, slope, path):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    colours = ("tab:blue", "tab:orange", "tab:green", "tab:red")
    for (d, info), colour in zip(per_d.items(), colours, strict=False):
        ref = info["coarse_reference"] + info["refined_reference"]
        ref = sorted(ref, key=lambda p: p["dt"])
        axes[0].errorbar(
            [p["dt"] for p in ref], [p["p_L"] for p in ref],
            yerr=[p["stderr"] for p in ref], fmt="o-", ms=3, lw=0.8,
            color=colour, capsize=2, label=f"d = {d} (reference)",
        )
        if info["coarse_guppy"]:
            g = info["coarse_guppy"]
            axes[0].errorbar(
                [p["dt"] for p in g], [p["p_L"] for p in g],
                yerr=[p["stderr"] for p in g], fmt="s", ms=7, mfc="none",
                lw=0, color=colour, capsize=3, label=f"d = {d} (Guppy)",
            )
        if math.isfinite(info["dt_star_fitted"]):
            axes[0].axvline(info["dt_star_fitted"], color=colour, ls="-", lw=0.8, alpha=0.5)
        axes[0].axvline(info["dt_star_eq6"], color=colour, ls=":", lw=1.2)
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[0].set_xlabel("syndrome interval $\\Delta t$")
    axes[0].set_ylabel("logical error probability $p_L$")
    axes[0].set_title(f"Guppy vs reference, T = {T}\n"
                      "dotted = Eq. 6 $\\Delta t^*$, solid = fitted")
    axes[0].legend(fontsize=7)
    axes[0].grid(alpha=0.3, which="both")

    ds = [int(d) for d, i in per_d.items() if math.isfinite(i["inv_dt_star"])]
    ys = [per_d[str(d)]["inv_dt_star"] for d in ds]
    es = [per_d[str(d)]["inv_dt_star_err"] for d in ds]
    axes[1].errorbar(ds, ys, yerr=es, fmt="o", ms=6, capsize=3, color="k",
                     label="fitted $1/\\Delta t^*$")
    axes[1].plot(ds, [1.0 / ansatz.optimal_interval(d=d, lam=LAM, g=G) for d in ds],
                 "^--", ms=6, color="tab:purple", label="Eq. 6")
    if slope:
        xs = np.linspace(min(ds) - 1, max(ds) + 1, 50)
        axes[1].plot(xs, slope["eta"] * xs + slope["z"], "-", lw=1, color="tab:red",
                     label=f"fit: $\\eta$ = {slope['eta']:.3f} $\\pm$ "
                           f"{slope['eta_stderr']:.3f}")
    axes[1].set_xlabel("code distance $d$")
    axes[1].set_ylabel("$1/\\Delta t^*$")
    axes[1].set_title("Optimal rate vs distance")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
