"""Stage 6c: dt sweeps on the validated Stim twin -- how does dt* scale with d?

Twin only (no Selene), in the configuration Stage 6b validated against Selene:
``order="hook_safe"``, ``native_charging=True``.

Variants
* A ``depolarizing_kops`` (Z memory): k charged DEPOLARIZE1(p_1q) ops per unit
  time on every data qubit -- the paper's linear idle law. Ansatz (Eq. 5, D6):
  R = (1/dt) C_d (1 + lam dt)^alpha, dt* = 1 / (lam (alpha - 1)).
* B ``dephasing_h2`` (X memory, f_lin = 0): p_z = sin^2(0.5 f_quad dt), i.e.
  quadratic idle. R = (1/dt) C_d (1 + kap dt^2)^alpha, dt* = 1/sqrt(kap (2 alpha - 1)).
with alpha = g (d + 1) / 2 and C_d = A d^-beta (p/p_th)^((d+1)/2); at circuit
level there is no single p, so (p/p_th) is one fitted constant, exp(c).

Figure of merit: R from Eq. A3 inverted, R = -ln(1 - 2 p_L) / (2 T_actual), as
in Stage 4a. Each point picks its own round count n (>= 2d, <= NMAX_ROUNDS) so
p_L per shot sits near P_L_TARGET: cost per failure does not depend on n, but
more rounds shrink the T-independent boundary error (preparation, final
readout) that would otherwise bias R by ~1/n. T = n dt. A 2T spot check per
variant verifies R does not depend on T. Points with p_L > P_L_SATURATED are
flagged and excluded from fits (Eq. A3 is ill-conditioned near 1/2).

Grid per (variant, d): 12 log-spaced dt over [0.1, 10] x the coarse dt*
(COARSE_DTSTAR, from a coarse scan done before this script), plus the
d-independent baseline dt (DT_BASE); then 6 refinement points log-spaced over
[1/1.6, 1.6] x the measured minimum. Shots: sized for >= TARGET_FAILURES.

Seeds (6b discipline): every point owns a disjoint block of SEED_STRIDE seeds,
recorded in points.csv. Pilot batches use block offsets [0, PILOT_SEEDS), the
main run one seed per batch from PILOT_SEEDS on. A run needing more batches
than its block holds is refused.

Budget: the whole run (pilots included) is capped by --budget-min (45 by
default); a point whose predicted cost no longer fits is skipped and recorded.

Usage:
    python scripts/stage6c_sweep.py --plan-only     # pilots + time plan, no main runs
    python scripts/stage6c_sweep.py --smoke         # one d=3 point end to end
    python scripts/stage6c_sweep.py                 # full sweep + analysis
    python scripts/stage6c_sweep.py --analyse DIR   # re-run the analysis only
Results: results/stage6c/<UTC timestamp>[-smoke|-plan]/
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import importlib.metadata as md
import json
import math
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_timing import ansatz  # noqa: E402
from qec_timing.circuit_level import (  # noqa: E402
    CircuitNoise,
    ReadoutNoise,
    build_constant_interval_memory,
    rate_per_unit_time,
)
from qec_timing.reference.sampler import clopper_pearson  # noqa: E402

# ---------------------------------------------------------------------------
# Configuration (chosen before this script; see docs/stage6c_report.md)
# ---------------------------------------------------------------------------

# Gate noise: the 6b ratio (p_2q/p_1q = 50, p_meas = 0.4 p_2q), level chosen so
# Lambda >= 2 per d+2 at every d (measured 2.6-3.2 in both memories, idle off).
P2, P1, PM, PINIT = 4e-3, 8e-5, 1.6e-3, 0.0
K_A = 100.0                   # variant A: idle ops per data qubit per unit time (p_op = p_1q)
FQ_B, FL_B = 0.1, 0.0         # variant B: H2-style quadratic dephasing
DISTANCES = (3, 5, 7, 9, 11)  # planned ~30 min of the 45 min budget; d = 11 runs last

# Coarse dt* per (variant, d) from the pre-sweep scan (2d rounds, <= 200 failures).
# d = 11 values are extrapolated from d = 3-9 (not scanned).
COARSE_DTSTAR = {
    "A": {3: 2.0, 5: 0.9, 7: 0.8, 9: 0.5, 11: 0.4},
    "B": {3: 1.1, 5: 1.0, 7: 0.75, 9: 0.65, 11: 0.55},
}
DT_BASE = {"A": 2.0, "B": 1.1}  # d-independent baseline: coarse dt* at d = 3

GRID_N, GRID_LO, GRID_HI = 12, 0.1, 10.0
REFINE_N, REFINE_SPAN = 6, 1.6
TARGET_FAILURES = 1000
P_L_TARGET = 0.2
P_L_SATURATED = 0.4
NMAX_ROUNDS = 60
TICK = 1e-3
BATCH = 20_000
MAX_SHOT_FACTOR = 3

SEED_STRIDE = 100_000
PILOT_SEEDS = 1_000
VARIANTS = ("A", "B")
SLOT_BASE, SLOT_REFINE0, SLOT_SPOT = 12, 20, 40

PALETTE_D = {3: "#86b6ef", 5: "#5598e7", 7: "#2a78d6", 9: "#184f95", 11: "#0d366b"}
INK, INK2, MUTED, GRID_C, AXIS_C, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
FIT_POWER, FIT_LINEAR, FIT_ANSATZ = "#2a78d6", "#eb6834", "#1baf7a"


@dataclass(frozen=True)
class Point:
    variant: str
    d: int
    dt: float
    kind: str        # grid | base | refine | spot
    slot: int
    rounds_factor: int = 1

    @property
    def basis(self) -> str:
        return "Z" if self.variant == "A" else "X"

    @property
    def label(self) -> str:
        return f"{self.variant}_d{self.d}_{self.kind}{self.slot}_dt{self.dt:g}"


def noise_for(variant: str) -> CircuitNoise:
    common = dict(p1=P1, p2=P2, p_init=PINIT, readout=ReadoutNoise.symmetric(PM))
    if variant == "A":
        return CircuitNoise.depolarizing_kops(k=K_A, p_op=P1, **common)
    return CircuitNoise.dephasing_h2(f_quad=FQ_B, f_lin=FL_B, **common)


def snap(dt: float) -> float:
    return max(TICK, round(dt / TICK) * TICK)


def seed_block(seed: int, pt: Point) -> int:
    index = (VARIANTS.index(pt.variant) * 16 + (pt.d - 3) // 2) * 128 + pt.slot
    return seed + index * SEED_STRIDE


def build(pt: Point, n_rounds: int):
    return build_constant_interval_memory(
        pt.d, pt.dt, n_rounds * pt.dt, noise_for(pt.variant), basis=pt.basis,
        order="hook_safe", native_charging=True)


def grid_points(variant: str, d: int) -> list[Point]:
    c = COARSE_DTSTAR[variant][d]
    dts = [snap(x) for x in c * np.geomspace(GRID_LO, GRID_HI, GRID_N)]
    pts = [Point(variant, d, dt, "grid", i) for i, dt in enumerate(dts)]
    base = snap(DT_BASE[variant])
    if all(abs(base / p.dt - 1) > 0.02 for p in pts):
        pts.append(Point(variant, d, base, "base", SLOT_BASE))
    return pts


def p_L_from_eps(eps: float, n: int) -> float:
    return 0.5 * (1.0 - math.exp(-2.0 * eps * n))


def rounds_for(eps: float, d: int, factor: int) -> int:
    n = math.ceil(-math.log(1 - 2 * P_L_TARGET) / (2 * max(eps, 1e-12)))
    return factor * int(min(max(n, 2 * d), NMAX_ROUNDS))


def run_batches(mem, decoder, shots: int, first_seed: int,
                target: int = TARGET_FAILURES) -> tuple[int, int]:
    """Planned shots, then onward (one seed per batch) until `target` failures,
    up to MAX_SHOT_FACTOR x planned: the pilot's rate is only good to ~20%."""
    failures = done = 0
    b = 0
    cap = MAX_SHOT_FACTOR * shots
    while done < shots or (failures < target and done < cap):
        n = BATCH if done >= shots else min(BATCH, shots - done)
        dets, obs = mem.sample(n, seed=first_seed + b)
        failures += int(np.count_nonzero(np.any(decoder.decode_batch(dets) != obs, axis=1)))
        done += n
        b += 1
    return failures, done


def pilot(pt: Point, block: int) -> dict:
    """Cheap estimate of the failure rate per round and the cost per shot-round."""
    n = 2 * pt.d
    t0 = time.perf_counter()
    mem = build(pt, n)
    dec = mem.decoder()
    t_build = time.perf_counter() - t0
    t1 = time.perf_counter()
    failures = shots = b = 0
    while failures < 30 and shots < 60_000 and time.perf_counter() - t1 < 2.0:
        dets, obs = mem.sample(2_000, seed=block + b)
        failures += int(np.count_nonzero(np.any(dec.decode_batch(dets) != obs, axis=1)))
        shots += 2_000
        b += 1
    us_round = (time.perf_counter() - t1) / (shots * n) * 1e6
    pl = max(failures, 1) / shots
    eps = rate_per_unit_time(min(pl, 0.45), mem.schedule.total_time) * pt.dt
    return dict(eps=eps, us_round=us_round, build_s=t_build, pilot_failures=failures,
                pilot_shots=shots)


def plan(pt: Point, pil: dict) -> dict:
    n = rounds_for(pil["eps"], pt.d, pt.rounds_factor)
    pl = p_L_from_eps(pil["eps"], n)
    shots = max(2_000, math.ceil(1.15 * TARGET_FAILURES / max(pl, 1e-9)))
    # Build scales ~linearly with rounds; sampling+decoding per shot-round.
    cost = pil["build_s"] * n / (2 * pt.d) + shots * n * pil["us_round"] * 1e-6
    return dict(n_rounds=n, p_L_pred=pl, shots=shots, cost_s=cost)


def measure_point(pt: Point, pl_plan: dict, block: int) -> dict:
    n_batches = math.ceil(MAX_SHOT_FACTOR * pl_plan["shots"] / BATCH) + 1
    if PILOT_SEEDS + n_batches > SEED_STRIDE:
        raise ValueError(f"{pt.label}: {n_batches} batches overflow the seed block")
    t0 = time.perf_counter()
    mem = build(pt, pl_plan["n_rounds"])
    dec = mem.decoder()
    failures, shots = run_batches(mem, dec, pl_plan["shots"], block + PILOT_SEEDS)
    T = mem.schedule.total_time
    pl = failures / shots
    lo, hi = clopper_pearson(failures, shots, 0.95)
    R = rate_per_unit_time(pl, T)
    return dict(n_rounds=mem.n_rounds, T_actual=T, shots=shots, failures=failures, p_L=pl,
                p_L_lo=lo, p_L_hi=hi, R=R, R_lo=rate_per_unit_time(lo, T),
                R_hi=rate_per_unit_time(min(hi, 0.4999), T),
                saturated=pl > P_L_SATURATED, seconds=time.perf_counter() - t0,
                seed_block_start=block, seed_block_end=block + SEED_STRIDE,
                batches=math.ceil(shots / BATCH))


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------


def R_of(p, T):
    return -np.log1p(-2.0 * np.asarray(p, float)) / (2.0 * np.asarray(T, float))


def shape_min(dts, Rs, sds, variant, x0=None):
    """Per-d fit of the ansatz SHAPE, log R = a - log dt + alpha_d log(1 + lam_d u),
    u = dt (A) or dt^2 (B), with a, alpha_d, lam_d free at this d.

    Primary dt* estimator. A symmetric parabola in log dt is biased toward small
    dt (R rises like 1/dt on the left but like dt^(alpha-1) or dt^(2 alpha-1) on
    the right; -5 to -13% on synthetic ansatz data). Returns
    (dt*, R*, params) or None if the fit has no interior minimum.
    """
    from scipy.optimize import least_squares
    x, y, s = np.asarray(dts), np.log(Rs), np.maximum(np.asarray(sds), 1e-9)
    if len(x) < 4:
        return None
    u = x if variant == "A" else x ** 2

    def resid(th):
        a, al, ll = th
        return (a - np.log(x) + al * np.log1p(np.exp(ll) * u) - y) / s

    i0 = int(np.argmin(y))
    th0 = x0 if x0 is not None else np.array([y[i0] + math.log(x[i0]), 2.0,
                                               -math.log(u[i0])])
    try:
        th = least_squares(resid, th0, bounds=([-50, 0.01, -20], [50, 50, 20])).x
    except ValueError:
        return None
    a, al, ll = th
    lam = math.exp(ll)
    if variant == "A":
        if al <= 1:
            return None
        ts = 1 / (lam * (al - 1))
    else:
        if al <= 0.5:
            return None
        ts = 1 / math.sqrt(lam * (2 * al - 1))
    us = ts if variant == "A" else ts ** 2
    return float(ts), float(math.exp(a - math.log(ts) + al * math.log1p(lam * us))), th


def parabola_min(dts, Rs, sds, centre):
    """Cross-check: weighted parabola in log dt within a factor 2 of `centre`."""
    x, y, s = np.log(dts), np.log(Rs), np.asarray(sds)
    sel = np.abs(x - math.log(centre)) <= math.log(2.0)
    if sel.sum() < 4:
        return None
    w = 1.0 / np.maximum(s[sel], 1e-9) ** 2
    c2, c1, c0 = np.polyfit(x[sel], y[sel], 2, w=np.sqrt(w))
    if c2 <= 0:
        return None
    xs = -c1 / (2 * c2)
    return float(np.exp(xs)), float(np.exp(c0 - c1 ** 2 / (4 * c2))), (c2, c1, c0)


def ansatz_logR(theta, d, dt, variant):
    logA, beta, c, g, loglam = theta
    alpha = g * (d + 1) / 2
    idle = np.log1p(np.exp(loglam) * (dt if variant == "A" else dt ** 2))
    return -np.log(dt) + logA - beta * np.log(d) + (d + 1) / 2 * c + alpha * idle


def ansatz_dtstar(theta, d, variant):
    _, _, _, g, loglam = theta
    alpha = g * (d + 1) / 2
    lam = math.exp(loglam)
    if variant == "A":
        return 1 / (lam * (alpha - 1)) if alpha > 1 else float("inf")
    return 1 / math.sqrt(lam * (2 * alpha - 1)) if alpha > 0.5 else float("inf")


def fit_ansatz(d, dt, y, s, variant, theta0=None):
    from scipy.optimize import least_squares
    theta0 = theta0 if theta0 is not None else np.array([0.0, 1.0, -1.0, 0.8, 0.0])
    res = least_squares(lambda th: (ansatz_logR(th, d, dt, variant) - y) / s, theta0,
                        bounds=([-30, -10, -10, 0.05, -15], [30, 10, 5, 5, 10]))
    chi2 = float(np.sum(res.fun ** 2))
    return res.x, chi2, len(y) - len(theta0)


def linfit(x, y, s):
    w = 1 / np.asarray(s) ** 2
    X = np.vstack([np.ones_like(x), x]).T
    cov = np.linalg.inv(X.T @ (w[:, None] * X))
    b = cov @ X.T @ (w * y)
    chi2 = float(np.sum(w * (y - X @ b) ** 2))
    return b, cov, chi2


def analyse(out: Path, n_boot: int = 400, seed: int = 6) -> dict:
    rows = [r for r in csv.DictReader(open(out / "points.csv")) if r["status"] == "done"]
    rng = np.random.default_rng(seed)
    result = {"n_boot": n_boot, "bootstrap_seed": seed, "variants": {}}
    per_variant = {}
    for v in VARIANTS:
        vr = [r for r in rows if r["variant"] == v and r["kind"] != "spot"
              and r["saturated"] == "False" and int(r["failures"]) > 0]
        if not vr:
            continue
        ds = sorted({int(r["d"]) for r in vr})
        data = {}
        for d in ds:
            rr = sorted((r for r in vr if int(r["d"]) == d), key=lambda r: float(r["dt"]))
            data[d] = dict(dt=np.array([float(r["dt"]) for r in rr]),
                           n=np.array([int(r["shots"]) for r in rr]),
                           k=np.array([int(r["failures"]) for r in rr]),
                           T=np.array([float(r["T_actual"]) for r in rr]),
                           kind=[r["kind"] for r in rr])
        per_variant[v] = data

        def estimates(kmat):
            """dt*, R* per d and Gamma per d for one (possibly resampled) dataset."""
            out_d = {}
            for d in ds:
                D = data[d]
                p = kmat[d] / D["n"]
                ok = (p > 0) & (p < P_L_SATURATED)
                R = np.where(ok, R_of(np.clip(p, 1e-12, 0.4999), D["T"]), np.nan)
                sd = np.sqrt(np.clip(p * (1 - p), 1e-12, None) / D["n"]) / (
                    D["T"] * (1 - 2 * np.clip(p, 0, 0.4999)) * np.where(ok, R, 1.0))
                good = ok & np.isfinite(R)
                if good.sum() < 4:
                    continue
                fit = shape_min(D["dt"][good], R[good], sd[good], v, start.get(d))
                if fit is None:
                    continue
                start.setdefault(d, fit[2])
                centre = D["dt"][good][np.argmin(R[good])]
                par = parabola_min(D["dt"][good], R[good], sd[good], centre)
                base = [i for i, kd in enumerate(D["kind"]) if kd == "base"
                        or abs(D["dt"][i] / snap(DT_BASE[v]) - 1) < 0.02]
                gamma = float(R[base[0]] / fit[1]) if base and good[base[0]] else float("nan")
                out_d[d] = (fit[0], fit[1], gamma, par[0] if par else float("nan"),
                            float(fit[2][1]), float(math.exp(fit[2][2])))
            return out_d

        start: dict = {}

        k0 = {d: data[d]["k"] for d in ds}
        point = estimates(k0)
        boots = []
        for _ in range(n_boot):
            kb = {d: rng.binomial(data[d]["n"], data[d]["k"] / data[d]["n"]) for d in ds}
            boots.append(estimates(kb))
        dres = {}
        for d in point:
            bd = np.array([b[d][0] for b in boots if d in b])
            bg = np.array([b[d][2] for b in boots if d in b and np.isfinite(b[d][2])])
            bp = np.array([b[d][3] for b in boots if d in b and np.isfinite(b[d][3])])
            dres[d] = dict(
                dtstar=point[d][0], R_star=point[d][1], gamma=point[d][2],
                alpha_d=point[d][4], idle_coeff_d=point[d][5],
                dtstar_parabola=point[d][3],
                dtstar_parabola_ci95=([float(np.percentile(bp, 2.5)), float(np.percentile(bp, 97.5))]
                                      if len(bp) else None),
                dtstar_ci95=[float(np.percentile(bd, 2.5)), float(np.percentile(bd, 97.5))],
                sd_log_dtstar=float(np.std(np.log(bd))),
                gamma_ci95=([float(np.percentile(bg, 2.5)), float(np.percentile(bg, 97.5))]
                            if len(bg) else None),
                boot_ok=int(len(bd)))
        vres = {"per_d": dres}

        # Scaling fits on log dt* vs log d, and the two linearised forms.
        dd = np.array(sorted(dres), float)
        if len(dd) >= 3:
            ts = np.array([dres[int(d)]["dtstar"] for d in dd])
            sl = np.array([dres[int(d)]["sd_log_dtstar"] for d in dd])
            b, cov, chi2 = linfit(np.log(dd), np.log(ts), sl)
            slopes = []
            for bt in boots:
                if all(int(d) in bt for d in dd):
                    tb = np.array([bt[int(d)][0] for d in dd])
                    slopes.append(np.polyfit(np.log(dd), np.log(tb), 1, w=1 / sl)[0])
            vres["power_law"] = dict(
                slope=float(b[1]), intercept=float(b[0]), slope_sd=float(math.sqrt(cov[1, 1])),
                slope_ci95=[float(np.percentile(slopes, 2.5)), float(np.percentile(slopes, 97.5))],
                chi2=chi2, dof=len(dd) - 2, expected_slope=-1.0 if v == "A" else -0.5)
            for name, power in (("inv_dtstar_linear", 1), ("inv_dtstar2_linear", 2)):
                y = ts ** -power
                sy = power * y * sl
                bl, covl, chi2l = linfit(dd, y, sy)
                vres[name] = dict(intercept=float(bl[0]), slope=float(bl[1]),
                                  slope_sd=float(math.sqrt(covl[1, 1])), chi2=chi2l,
                                  dof=len(dd) - 2, d_intercept=float(-bl[0] / bl[1]))
            vres["expected_linear_form"] = "inv_dtstar_linear" if v == "A" else "inv_dtstar2_linear"

        # Global ansatz fit.
        d_all = np.concatenate([np.full(len(data[d]["dt"]), d, float) for d in ds])
        dt_all = np.concatenate([data[d]["dt"] for d in ds])
        n_all = np.concatenate([data[d]["n"] for d in ds])
        T_all = np.concatenate([data[d]["T"] for d in ds])
        k_all = np.concatenate([data[d]["k"] for d in ds])
        p_all = k_all / n_all
        R_all = R_of(p_all, T_all)
        s_all = np.sqrt(p_all * (1 - p_all) / n_all) / (T_all * (1 - 2 * p_all) * R_all)
        theta, chi2, dof = fit_ansatz(d_all, dt_all, np.log(R_all), s_all, v)
        th_boot = []
        for _ in range(min(n_boot, 200)):
            kb = rng.binomial(n_all, p_all)
            pb = np.clip(kb / n_all, 1e-9, 0.4999)
            yb = np.log(R_of(pb, T_all))
            th_boot.append(fit_ansatz(d_all, dt_all, yb, s_all, v, theta)[0])
        th_boot = np.array(th_boot)
        names = ["logA", "beta", "c", "g", "log_idle_coeff"]
        ci = {nm: [float(np.percentile(th_boot[:, i], 2.5)), float(np.percentile(th_boot[:, i], 97.5))]
              for i, nm in enumerate(names)}
        vres["ansatz"] = dict(
            A=float(math.exp(theta[0])), beta=float(theta[1]), p_over_pth=float(math.exp(theta[2])),
            g=float(theta[3]), idle_coeff=float(math.exp(theta[4])),
            idle_coeff_name="lambda" if v == "A" else "kappa",
            chi2=chi2, dof=dof, ci95=ci,
            dtstar_pred={int(d): ansatz_dtstar(theta, d, v) for d in ds},
            theta=[float(x) for x in theta])
        preds = np.array([ansatz_dtstar(theta, d, v) for d in ds])
        if np.all(np.isfinite(preds)) and len(ds) >= 2:
            vres["ansatz"]["implied_loglog_slope"] = float(
                np.polyfit(np.log(ds), np.log(preds), 1)[0])
        if v == "A":
            lam, g = math.exp(theta[4]), theta[3]
            vres["ansatz"]["gamma_pred"] = {
                int(d): float(ansatz.gamma_opt(snap(DT_BASE[v]), lam=lam, d=d, g=g))
                if ansatz.alpha(d, g) > 1 else None for d in ds}
        else:
            kap, g = math.exp(theta[4]), theta[3]
            gp = {}
            for d in ds:
                a = g * (d + 1) / 2
                ts_ = 1 / math.sqrt(kap * (2 * a - 1))
                tb = snap(DT_BASE[v])
                gp[int(d)] = float((ts_ / tb) * ((1 + kap * tb ** 2) / (1 + kap * ts_ ** 2)) ** a)
            vres["ansatz"]["gamma_pred"] = gp
        result["variants"][v] = vres

    (out / "fits.json").write_text(json.dumps(result, indent=2))
    try:
        plot(out, rows, result)
    except Exception as e:  # noqa: BLE001  plots must never lose the numbers
        print(f"plotting failed: {e!r}")
    return result


def plot(out: Path, rows, result) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FormatStrFormatter, LogLocator, NullFormatter
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "axes.edgecolor": AXIS_C,
        "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
        "axes.spines.top": False, "axes.spines.right": False, "grid.color": GRID_C,
        "font.size": 9, "axes.titlesize": 10, "legend.frameon": False, "lines.linewidth": 2,
    })

    # 1. R vs dt per d, minima marked.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    for ax, v in zip(axes, VARIANTS):
        vres = result["variants"].get(v, {}).get("per_d", {})
        for d in sorted({int(r["d"]) for r in rows if r["variant"] == v}):
            rr = sorted((r for r in rows if r["variant"] == v and int(r["d"]) == d
                         and r["kind"] != "spot" and int(r["failures"]) > 0),
                        key=lambda r: float(r["dt"]))
            dt = np.array([float(r["dt"]) for r in rr])
            R = np.array([float(r["R"]) for r in rr])
            lo = np.array([float(r["R_lo"]) for r in rr])
            hi = np.array([float(r["R_hi"]) for r in rr])
            sat = np.array([r["saturated"] == "True" for r in rr])
            col = PALETTE_D[d]
            ax.plot(dt[~sat], R[~sat], "-", color=col, lw=1.5, label=f"d = {d}")
            ax.errorbar(dt[~sat], R[~sat], yerr=[np.clip(R[~sat] - lo[~sat], 0, None),
                                                 np.clip(hi[~sat] - R[~sat], 0, None)],
                        fmt="o", ms=4, color=col, elinewidth=1, capsize=0)
            if sat.any():
                ax.plot(dt[sat], R[sat], "o", ms=5, mfc=SURFACE, mec=col, mew=1.2)
            if str(d) in map(str, vres):
                m = vres[d] if d in vres else vres[str(d)]
                ax.plot([m["dtstar"]], [m["R_star"]], marker="D", ms=8, color=col,
                        mec=SURFACE, mew=1.5, zorder=5)
        ax.axvline(snap(DT_BASE[v]), color=AXIS_C, lw=1, ls="--", zorder=0)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("syndrome interval dt"); ax.set_ylabel("logical error rate R per unit time")
        ax.set_title("A: depolarizing k-ops, Z memory" if v == "A"
                     else "B: quadratic dephasing, X memory", loc="left")
        ax.grid(True, which="major", lw=0.6)
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.yaxis.set_minor_formatter(NullFormatter())
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=len(labels), bbox_to_anchor=(0.5, 1.0))
    fig.text(0.5, 0.01, "diamond = fitted minimum dt*;  dashed line = d-independent baseline dt;  "
             "hollow = saturated p_L (excluded from fits);  bars = 95% binomial",
             ha="center", color=INK2, fontsize=8)
    fig.tight_layout(rect=(0, 0.04, 1, 0.93))
    fig.savefig(out / "R_vs_dt.png", dpi=150)
    plt.close(fig)

    # 2. dt* vs d, log-log, with power-law, linearised and ansatz fits.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, v in zip(axes, VARIANTS):
        vres = result["variants"].get(v)
        if not vres or "power_law" not in vres:
            continue
        per = {int(k): x for k, x in vres["per_d"].items()}
        dd = np.array(sorted(per), float)
        ts = np.array([per[int(d)]["dtstar"] for d in dd])
        lo = np.array([per[int(d)]["dtstar_ci95"][0] for d in dd])
        hi = np.array([per[int(d)]["dtstar_ci95"][1] for d in dd])
        xs = np.linspace(dd.min() * 0.9, dd.max() * 1.1, 100)
        pl = vres["power_law"]
        ax.plot(xs, np.exp(pl["intercept"]) * xs ** pl["slope"], color=FIT_POWER, lw=1.5,
                label=f"power law, slope {pl['slope']:.2f} (expect {pl['expected_slope']:g})")
        key = vres["expected_linear_form"]
        lin = vres[key]
        power = 1 if key == "inv_dtstar_linear" else 2
        yl = lin["intercept"] + lin["slope"] * xs
        ax.plot(xs[yl > 0], yl[yl > 0] ** (-1 / power), color=FIT_LINEAR, lw=1.5, ls="--",
                label=("1/dt* linear in d" if power == 1 else "1/dt*^2 linear in d"))
        an = vres["ansatz"]
        pred = [ansatz_dtstar(an["theta"], x, v) for x in xs]
        ax.plot(xs, pred, color=FIT_ANSATZ, lw=1.5, ls=":", label="ansatz (Eq. 5) fit")
        ax.errorbar(dd, ts, yerr=[np.clip(ts - lo, 0, None), np.clip(hi - ts, 0, None)],
                    fmt="o", ms=6, color=INK,
                    elinewidth=1, capsize=0, label="measured dt* (95% bootstrap)", zorder=5)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xticks(dd); ax.set_xticklabels([str(int(d)) for d in dd])
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 3.0, 5.0)))
        ax.yaxis.set_major_formatter(FormatStrFormatter("%g"))
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.set_xlabel("code distance d"); ax.set_ylabel("optimal interval dt*")
        ax.set_title("A: depolarizing k-ops" if v == "A" else "B: quadratic dephasing", loc="left")
        ax.grid(True, which="major", lw=0.6)
        ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "dtstar_vs_d.png", dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def git_head() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    except Exception:  # noqa: BLE001
        return "unknown"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--plan-only", action="store_true", help="pilots and time plan only")
    ap.add_argument("--smoke", action="store_true", help="one d=3 point end to end")
    ap.add_argument("--analyse", type=Path, default=None, help="re-run analysis on a results dir")
    ap.add_argument("--distances", type=int, nargs="+", default=list(DISTANCES))
    ap.add_argument("--budget-min", type=float, default=45.0)
    ap.add_argument("--seed", type=int, default=60_000_000)
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "stage6c")
    args = ap.parse_args(argv)

    if args.analyse:
        res = analyse(args.analyse)
        print(json.dumps(res, indent=2)[:4000])
        return 0

    stamp = _dt.datetime.now(_dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    suffix = "-smoke" if args.smoke else ("-plan" if args.plan_only else "")
    out = args.out / (stamp + suffix)
    out.mkdir(parents=True, exist_ok=False)
    log_file = open(out / "run.log", "w", encoding="utf-8")

    def log(msg: str) -> None:
        print(msg, flush=True)
        log_file.write(msg + "\n")
        log_file.flush()

    budget = args.budget_min * 60.0
    job_start = time.perf_counter()
    log(f"stage6c sweep -> {out}")
    log(f"noise p_2q={P2} p_1q={P1} p_meas={PM} p_init={PINIT}; A: k={K_A}; "
        f"B: f_quad={FQ_B} f_lin={FL_B}; d={args.distances}; budget {args.budget_min} min")

    rows: list[dict] = []

    def record(pt: Point, block: int, pil: dict, pl_plan: dict, res: dict | None, status: str):
        rows.append({**asdict(pt), "label": pt.label, "basis": pt.basis, "status": status,
                     "seed_block_start": block, "seed_block_end": block + SEED_STRIDE,
                     "pilot_eps": pil["eps"], "planned_rounds": pl_plan["n_rounds"],
                     "planned_shots": pl_plan["shots"], "planned_cost_s": pl_plan["cost_s"],
                     **(res or {})})

    def run_point(pt: Point, pil: dict | None = None) -> dict | None:
        block = seed_block(args.seed, pt)
        pil = pil or pilot(pt, block)
        pl_plan = plan(pt, pil)
        spent = time.perf_counter() - job_start
        if spent + pl_plan["cost_s"] > budget:
            log(f"[skip] {pt.label}: predicted {pl_plan['cost_s']:.0f}s, {budget - spent:.0f}s left")
            record(pt, block, pil, pl_plan, None, "skipped_budget")
            return None
        res = measure_point(pt, pl_plan, block)
        record(pt, block, pil, pl_plan, res, "done")
        log(f"[pt] {pt.label}: n={res['n_rounds']} T={res['T_actual']:.3f} "
            f"p_L={res['p_L']:.4f} ({res['failures']}/{res['shots']}) R={res['R']:.3e} "
            f"{'SATURATED ' if res['saturated'] else ''}{res['seconds']:.1f}s "
            f"(planned {pl_plan['cost_s']:.1f}s) seeds [{block}, {block + SEED_STRIDE})")
        return res

    if args.smoke:
        pt = Point("A", 3, snap(COARSE_DTSTAR["A"][3]), "grid", 99)
        run_point(pt)
        write(out, rows, args, job_start)
        log_file.close()
        return 0

    # Stage 1: grid + baseline, with pilots (also the plan).
    # Ordered by d (both variants per d) so a budget shortfall skips the largest d first.
    stage1 = [(v, d, pt) for d in args.distances for v in VARIANTS for pt in grid_points(v, d)]
    pilots = {}
    plan_cost = 0.0
    for v, d, pt in stage1:
        pil = pilot(pt, seed_block(args.seed, pt))
        pilots[pt] = pil
        pp = plan(pt, pil)
        plan_cost += pp["cost_s"]
        log(f"[pilot] {pt.label}: eps/round={pil['eps']:.2e} ({pil['pilot_failures']}/"
            f"{pil['pilot_shots']}) -> n={pp['n_rounds']} shots={pp['shots']} "
            f"~{pp['cost_s']:.1f}s")
    # Refinement and spot checks, costed from the stage-1 point nearest the coarse dt*.
    extra = 0.0
    for v in VARIANTS:
        for d in args.distances:
            near = min((pt for vv, dd, pt in stage1 if vv == v and dd == d),
                       key=lambda pt: abs(math.log(pt.dt / COARSE_DTSTAR[v][d])))
            c = plan(near, pilots[near])["cost_s"]
            extra += REFINE_N * c + (2 * c if d == 5 else 0)
    pilot_time = time.perf_counter() - job_start
    log(f"plan: pilots took {pilot_time:.0f}s; stage 1 ~{plan_cost:.0f}s; refinement + "
        f"spot checks ~{extra:.0f}s; total ~{pilot_time + plan_cost + extra:.0f}s "
        f"against {budget:.0f}s budget")
    if args.plan_only:
        log_file.close()
        return 0

    for v, d, pt in stage1:
        run_point(pt, pilots[pt])

    # Stage 2: refinement around each measured minimum.
    for d in args.distances:
        for v in VARIANTS:
            done = [r for r in rows if r["variant"] == v and r["d"] == d and r["status"] == "done"
                    and not r["saturated"] and r["failures"] > 0]
            if len(done) < 3:
                continue
            best = min(done, key=lambda r: r["R"])["dt"]
            have = [r["dt"] for r in rows if r["variant"] == v and r["d"] == d]
            for j, x in enumerate(np.geomspace(best / REFINE_SPAN, best * REFINE_SPAN, REFINE_N)):
                dt = snap(x)
                if any(abs(dt / h - 1) < 0.02 for h in have):
                    continue
                run_point(Point(v, d, dt, "refine", SLOT_REFINE0 + j))

    # Spot check: R must not depend on T. d = 5, at 2T, on the grid point nearest the
    # COARSE dt* -- chosen before the data. (Choosing the lowest measured R instead
    # biases the reference low: a winner's-curse z of about +2 in a d = 3, 5 trial.)
    for v in VARIANTS:
        done = [r for r in rows if r["variant"] == v and r["d"] == 5 and r["kind"] == "grid"
                and r["status"] == "done" and not r["saturated"] and r["failures"] > 0]
        if not done:
            continue
        ref = min(done, key=lambda r: abs(math.log(r["dt"] / COARSE_DTSTAR[v][5])))
        pt = Point(v, 5, ref["dt"], "spot", SLOT_SPOT, rounds_factor=2)
        res = run_point(pt)
        if res:
            z = (math.log(res["R"]) - math.log(ref["R"])) / math.hypot(
                (math.log(res["R_hi"]) - math.log(res["R_lo"])) / 3.92,
                (math.log(ref["R_hi"]) - math.log(ref["R_lo"])) / 3.92)
            log(f"[spot] {v} d=5 dt={ref['dt']:g}: R(T={ref['T_actual']:.3f})={ref['R']:.4e}, "
                f"R(2T={res['T_actual']:.3f})={res['R']:.4e}, z(log R)={z:+.2f}")
            rows[-1]["spot_ref_label"] = ref["label"]
            rows[-1]["spot_z_logR"] = z

    write(out, rows, args, job_start)
    res = analyse(out)
    for v, vres in res["variants"].items():
        log(f"[fit] {v}: " + ", ".join(f"d={d}: dt*={x['dtstar']:.3f} "
                                       f"[{x['dtstar_ci95'][0]:.3f},{x['dtstar_ci95'][1]:.3f}]"
                                       for d, x in vres["per_d"].items()))
        if "power_law" in vres:
            pl = vres["power_law"]
            log(f"[fit] {v}: log-log slope {pl['slope']:.3f} 95% [{pl['slope_ci95'][0]:.3f}, "
                f"{pl['slope_ci95'][1]:.3f}] (expect {pl['expected_slope']})")
    log(f"done in {time.perf_counter() - job_start:.0f}s")
    log_file.close()
    return 0


def write(out: Path, rows, args, job_start) -> None:
    keys = sorted({k for r in rows for k in r})
    with open(out / "points.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    summary = dict(
        stage="6c", git_head=git_head(), args={k: str(v) for k, v in vars(args).items()},
        noise=dict(p_2q=P2, p_1q=P1, p_meas=PM, p_init=PINIT), k_A=K_A, f_quad_B=FQ_B,
        f_lin_B=FL_B, coarse_dtstar=COARSE_DTSTAR, dt_base=DT_BASE,
        target_failures=TARGET_FAILURES, p_L_target=P_L_TARGET, n_max_rounds=NMAX_ROUNDS,
        seed_stride=SEED_STRIDE, versions={p: md.version(p) for p in ("stim", "PyMatching", "numpy")},
        seconds=time.perf_counter() - job_start,
        n_done=sum(r["status"] == "done" for r in rows),
        n_skipped=sum(r["status"] != "done" for r in rows))
    (out / "summary.json").write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    sys.exit(main())
