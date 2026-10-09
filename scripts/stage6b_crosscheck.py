"""Stage 6b: cross-check the Stage 6a Stim twin against Selene's native noise.

Selene runs ``circuit/program.py`` under ``DepolarizingErrorModel`` (native
p_1q / p_2q / p_meas, charged per compiled op) on the Stim backend. The twin is
``qec_timing.circuit_level`` with ``native_charging=True`` and
``mirror_guppy_program=True``. Both are decoded by the SAME decoder (the twin's
detector error model), so only the noise processes are being compared.

Idle noise:
* variant A (depolarizing_kops, Z memory): ``idle_ops = k * dt`` charged
  ``phased_x(q, 0, 0)`` per data qubit per round in Selene; the twin's
  ``DEPOLARIZE1(compose(p_1q, k * dt))``.
* variant B (dephasing_h2, X memory): program.py's injection with
  ``p_data = p_z(dt)`` (``h``-wrapped, so a Z fault); the twin mirrors the
  wrapper and the charge on each injected ``x``.

Pass criterion per point (99%):
* p_L: two-proportion z, |z| < 2.576 (Selene and twin 99% Clopper-Pearson CIs
  are also reported).
* detection events: for every round r, the mean number of detection events per
  shot agrees, Welch z over shots, Bonferroni-corrected to 99% across the
  n_rounds + 1 rounds of that point. Catches compensating errors that leave p_L
  unchanged.

Guppy rules (CLAUDE.md): ONE run() per configuration with shots sized up front
(~``--target-failures`` expected failures from the twin's p_L), base_seed per
shot via the program, and the whole job capped by ``--budget-min``. Shots are
scaled down uniformly if the cost model says the budget would be exceeded, and
a point whose predicted cost no longer fits is skipped (recorded), never
truncated mid-run.

Usage:
    python scripts/stage6b_crosscheck.py                # full cross-check
    python scripts/stage6b_crosscheck.py --smoke        # one small point
Results: results/stage6b/<UTC timestamp>[-smoke]/
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
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_timing.circuit import build_surface_code_program  # noqa: E402
from qec_timing.circuit_level import (  # noqa: E402
    CircuitNoise,
    ReadoutNoise,
    build_constant_interval_memory,
    dephasing_h2_probability,
)
from qec_timing.reference.sampler import clopper_pearson  # noqa: E402

# Approved noise levels (p_2q / p_1q = 50).
P1, P2, PM, PINIT = 1e-4, 5e-3, 2e-3, 0.0
T = 12.0
K_A = 100                    # variant A: idle ops per data qubit per unit time
FQ_B, FL_B = 0.1, 0.0        # variant B: p_z(dt) = sin^2(0.05 dt)
CONF = 0.99
SECONDS_PER_RUN = 9.3        # fixed run() cost (docs/performance.md)


@dataclass(frozen=True)
class Point:
    variant: str   # "A" or "B"
    d: int
    order: str
    dt: float

    @property
    def basis(self) -> str:
        return "Z" if self.variant == "A" else "X"

    @property
    def k(self) -> int | None:
        return K_A if self.variant == "A" else None

    @property
    def idle_ops(self) -> int:
        if self.variant != "A":
            return 0
        n = K_A * self.dt
        if abs(n - round(n)) > 1e-9:
            raise ValueError(f"k*dt = {n} is not an integer")
        return int(round(n))

    @property
    def p_data(self) -> float:
        return dephasing_h2_probability(self.dt, FQ_B, FL_B) if self.variant == "B" else 0.0

    @property
    def label(self) -> str:
        return f"{self.variant}_d{self.d}_{self.order}_dt{self.dt:g}"


POINTS = [
    *(Point("A", d, "hook_safe", dt) for d in (3, 5) for dt in (1.0, 4.0, 12.0)),
    *(Point("B", d, "hook_safe", dt) for d in (3, 5) for dt in (1.0, 4.0, 6.0)),
]
SMOKE = Point("B", 3, "hook_safe", 4.0)


def twin_noise(pt: Point) -> CircuitNoise:
    common = dict(p1=P1, p2=P2, p_init=PINIT, readout=ReadoutNoise.symmetric(PM))
    if pt.variant == "A":
        return CircuitNoise.depolarizing_kops(k=K_A, p_op=P1, **common)
    return CircuitNoise.dephasing_h2(f_quad=FQ_B, f_lin=FL_B, **common)


def build_twin(pt: Point):
    return build_constant_interval_memory(
        pt.d, pt.dt, T, twin_noise(pt), basis=pt.basis, order=pt.order,
        native_charging=True, mirror_guppy_program=True)


def native_op_count(pt: Point, twin) -> int:
    """Charged native ops per Selene shot, from the twin's gate content."""
    n = 0
    for inst in twin.circuit.flattened():
        targets = inst.targets_copy()
        if not targets or targets[0].is_measurement_record_target:
            continue
        if inst.name == "CX":
            n += 5 * (len(targets) // 2)    # Rzz + 4 single-qubit ops
        elif inst.name == "H":
            n += 2 * len(targets)
        elif inst.name in ("R", "M"):
            n += len(targets)
    return n + pt.idle_ops * pt.d * pt.d * twin.n_rounds


def per_round_counts(dets: np.ndarray, n_rounds: int, n_checks: int) -> np.ndarray:
    """(shots, n_rounds + 1): detection events per shot in each round."""
    return np.asarray(dets, dtype=np.int64).reshape(-1, n_rounds + 1, n_checks).sum(axis=2)


def welch_z(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    se = np.sqrt(a.var(axis=0, ddof=1) / len(a) + b.var(axis=0, ddof=1) / len(b))
    diff = a.mean(axis=0) - b.mean(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(se > 0, diff / se, np.where(diff == 0, 0.0, np.inf))


def two_proportion_z(k1: int, n1: int, k2: int, n2: int) -> float:
    p = (k1 + k2) / (n1 + n2)
    se = math.sqrt(p * (1 - p) * (1 / n1 + 1 / n2))
    return (k1 / n1 - k2 / n2) / se if se > 0 else 0.0


def run_twin(pt: Point, twin, decoder, seed: int, target_failures: int, max_shots: int):
    """Twin p_L (adaptive, cheap) and per-round detection counts."""
    sampler = twin.sampler(seed)
    n_checks = twin.layout.n_checks
    failures = shots = 0
    counts, first = [], None
    while shots < max_shots and (failures < target_failures or shots < 20_000):
        n = min(20_000, max_shots - shots)
        dets, obs = sampler.sample(n)
        failures += int(np.count_nonzero(np.any(decoder.decode_batch(dets) != obs, axis=1)))
        counts.append(per_round_counts(dets, twin.n_rounds, n_checks))
        if first is None:
            first = np.asarray(dets, dtype=np.bool_)
        shots += n
    return failures, shots, np.concatenate(counts), first


def busy_duplicates(events: np.ndarray) -> int:
    """Shots with >= 1 detection whose full pattern repeats an earlier shot."""
    flat = np.asarray(events, dtype=np.bool_).reshape(len(events), -1)
    busy = flat[flat.any(axis=1)]
    return len(busy) - len({row.tobytes() for row in busy})


def run_selene(pt: Point, twin, decoder, shots: int, seed: int, n_processes: int):
    prog = build_surface_code_program(pt.d, basis=pt.basis, order=pt.order,
                                      idle_ops=pt.idle_ops)
    from selene_sim import DepolarizingErrorModel
    em = DepolarizingErrorModel(p_1q=P1, p_2q=P2, p_meas=PM, p_init=PINIT)
    inst = prog.emulator().with_error_model(em).with_shots(shots).with_seed(seed)
    if n_processes > 1:
        inst = inst.with_n_processes(n_processes)
    t0 = time.perf_counter()
    result = inst.run(p_data=pt.p_data, p_read=0.0, dt=pt.dt, n_rounds=twin.n_rounds,
                      base_seed=seed + 1, inject_round=-1, inject_qubit=-1, inject_patch=0)
    elapsed = time.perf_counter() - t0
    pred = prog.decode(result, twin.n_rounds, decoder)[:, 0]
    obs = prog.observables(result)[:, 0]
    failures = int(np.count_nonzero(pred != obs))
    events = np.asarray(prog.detection_events(result, twin.n_rounds))
    counts = events.reshape(len(events), twin.n_rounds + 1, -1).sum(axis=2).astype(np.int64)
    return failures, len(events), counts, elapsed, events


def git_head() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    except Exception:  # noqa: BLE001
        return "unknown"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--smoke", action="store_true", help="one small point end to end")
    ap.add_argument("--plan-only", action="store_true",
                    help="twin predictions and the shot/time plan; no Selene runs")
    ap.add_argument("--budget-min", type=float, default=30.0, help="cap on Selene wall time")
    ap.add_argument("--target-failures", type=int, default=100)
    ap.add_argument("--min-shots", type=int, default=200)
    ap.add_argument("--smoke-shots", type=int, default=300)
    ap.add_argument("--n-processes", type=int, default=4)
    ap.add_argument("--us-per-op", type=float, default=None,
                    help="Selene cost model, microseconds per native op per shot "
                         "(single process); default from the smoke calibration")
    ap.add_argument("--twin-target-failures", type=int, default=2_000)
    ap.add_argument("--twin-max-shots", type=int, default=400_000)
    ap.add_argument("--seed", type=int, default=20261009)
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "stage6b")
    args = ap.parse_args(argv)

    stamp = _dt.datetime.now(_dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    suffix = "-smoke" if args.smoke else ("-plan" if args.plan_only else "")
    out = args.out / (stamp + suffix)
    out.mkdir(parents=True, exist_ok=False)
    log_file = open(out / "run.log", "w", encoding="utf-8")

    def log(msg: str) -> None:
        print(msg, flush=True)
        log_file.write(msg + "\n")
        log_file.flush()

    points = [SMOKE] if args.smoke else POINTS
    us_per_op = args.us_per_op if args.us_per_op is not None else US_PER_OP_DEFAULT
    log(f"stage6b cross-check -> {out}")
    log(f"noise: p_1q={P1} p_2q={P2} p_meas={PM} p_init={PINIT}; T={T}; "
        f"A: k={K_A}; B: f_quad={FQ_B} f_lin={FL_B}")
    log(f"budget {args.budget_min} min, target {args.target_failures} failures/point, "
        f"n_processes={args.n_processes}, cost model {us_per_op} us/op/shot")

    # 1. Twin first (cheap): p_L prediction, shot sizing, reference statistics.
    plan = []
    for i, pt in enumerate(points):
        seed = args.seed + 1000 * i
        twin = build_twin(pt)
        dec = twin.decoder()
        t0 = time.perf_counter()
        tf, ts, tcounts, tfirst = run_twin(pt, twin, dec, seed + 2,
                                           args.twin_target_failures, args.twin_max_shots)
        p_twin = tf / ts
        ops = native_op_count(pt, twin)
        wanted = args.smoke_shots if args.smoke else max(
            args.min_shots, math.ceil(args.target_failures / max(p_twin, 1e-4)))
        plan.append(dict(pt=pt, seed=seed, twin=twin, dec=dec, tf=tf, ts=ts,
                         tcounts=tcounts, tfirst=tfirst, ops=ops, wanted=wanted))
        log(f"[twin] {pt.label}: rounds={twin.n_rounds} p_L={p_twin:.4f} ({tf}/{ts}) "
            f"ops/shot={ops} -> wanted {wanted} shots, "
            f"~{wanted * sec_per_shot(ops, us_per_op, args.n_processes):.0f}s "
            f"[{time.perf_counter() - t0:.1f}s]")

    # 2. Budget. Shots for each run are fixed just before it, from the remaining
    #    budget and the worst per-op cost seen so far (cost is unstable ~3x).
    budget = args.budget_min * 60.0
    usable = USABLE_FRACTION * budget

    def sized(remaining_plan, spent, us):
        need = sum(SECONDS_PER_RUN + q["wanted"] * sec_per_shot(q["ops"], us, args.n_processes)
                   for q in remaining_plan)
        scale = min(1.0, max(0.0, usable - spent) / need) if need else 1.0
        return scale

    scale0 = sized(plan, 0.0, us_per_op)
    log(f"plan: wanted {sum(SECONDS_PER_RUN + q['wanted'] * sec_per_shot(q['ops'], us_per_op, args.n_processes) for q in plan):.0f}s "
        f"against {usable:.0f}s usable -> initial shot scale {scale0:.3f}")
    for q in plan:
        log(f"   {q['pt'].label}: ~{max(args.min_shots, int(q['wanted'] * scale0))} shots")
    if args.plan_only:
        log_file.close()
        return 0

    # 3. Selene: one run() per configuration.
    rows, det_rows = [], []
    job_start = time.perf_counter()
    z_crit = norm.ppf(1 - (1 - CONF) / 2)
    us_seen = us_per_op
    for idx, p in enumerate(plan):
        pt, twin = p["pt"], p["twin"]
        spent = time.perf_counter() - job_start
        if args.smoke:
            shots = p["wanted"]
        else:
            scale = sized(plan[idx:], spent, us_seen)
            shots = max(args.min_shots, int(p["wanted"] * scale))
        predicted = SECONDS_PER_RUN + shots * sec_per_shot(p["ops"], us_seen, args.n_processes)
        if not args.smoke and spent + predicted > budget:
            log(f"[skip] {pt.label}: predicted {predicted:.0f}s, {budget - spent:.0f}s left")
            rows.append(dict(asdict(pt), status="skipped_budget"))
            continue
        sf, ss, scounts, elapsed, sevents = run_selene(pt, twin, p["dec"], shots,
                                                       p["seed"], args.n_processes)
        tf, ts, tcounts = p["tf"], p["ts"], p["tcounts"]
        dups = busy_duplicates(sevents)
        twin_dups = busy_duplicates(p["tfirst"][:ss])
        z_pl = two_proportion_z(sf, ss, tf, ts)
        s_lo, s_hi = clopper_pearson(sf, ss, CONF)
        t_lo, t_hi = clopper_pearson(tf, ts, CONF)
        zr = welch_z(scounts, tcounts)
        z_crit_rounds = norm.ppf(1 - (1 - CONF) / 2 / zr.size)
        det_ok = bool(np.all(np.abs(zr) < z_crit_rounds))
        pl_ok = abs(z_pl) < z_crit
        measured_us = elapsed / ss / p["ops"] * 1e6 * parallel_speedup(args.n_processes)
        us_seen = max(us_seen, measured_us)
        rows.append(dict(
            asdict(pt), basis=pt.basis, k=pt.k, idle_ops=pt.idle_ops, p_data=pt.p_data,
            n_rounds=twin.n_rounds, T_actual=twin.schedule.total_time,
            selene_failures=sf, selene_shots=ss, selene_p_L=sf / ss,
            selene_ci_low=s_lo, selene_ci_high=s_hi,
            twin_failures=tf, twin_shots=ts, twin_p_L=tf / ts,
            twin_ci_low=t_lo, twin_ci_high=t_hi, z_p_L=z_pl,
            det_max_abs_z=float(np.max(np.abs(zr))), det_z_crit=float(z_crit_rounds),
            pass_p_L=pl_ok, pass_detection=det_ok, status="pass" if pl_ok and det_ok else "FAIL",
            selene_seconds=elapsed, ops_per_shot=p["ops"], measured_us_per_op=measured_us,
            duplicate_busy_shots_selene=dups, duplicate_busy_shots_twin_same_n=twin_dups))
        for r in range(zr.size):
            det_rows.append(dict(label=pt.label, round=r,
                                 selene_mean=float(scounts[:, r].mean()),
                                 twin_mean=float(tcounts[:, r].mean()), z=float(zr[r])))
        log(f"[selene] {pt.label}: p_L={sf / ss:.4f} [{s_lo:.4f},{s_hi:.4f}] ({sf}/{ss}) "
            f"twin {tf / ts:.4f} [{t_lo:.4f},{t_hi:.4f}] z={z_pl:+.2f} | "
            f"det max|z|={np.max(np.abs(zr)):.2f} (crit {z_crit_rounds:.2f}) | "
            f"{elapsed:.0f}s, {measured_us:.2f} us/op, dup selene/twin={dups}/{twin_dups} "
            f"-> {rows[-1]['status']}")

    # 4. Write results.
    if rows:
        keys = sorted({k for r in rows for k in r})
        with open(out / "points.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(rows)
    if det_rows:
        with open(out / "detection_rates.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(det_rows[0]))
            w.writeheader()
            w.writerows(det_rows)
    summary = dict(
        stage="6b", smoke=args.smoke, git_head=git_head(), args={k: str(v) for k, v in vars(args).items()},
        noise=dict(p_1q=P1, p_2q=P2, p_meas=PM, p_init=PINIT), T=T, k_A=K_A,
        f_quad_B=FQ_B, f_lin_B=FL_B, confidence=CONF,
        versions={pkg: md.version(pkg) for pkg in
                  ("stim", "PyMatching", "selene-sim", "guppylang", "numpy")},
        selene_seconds_total=time.perf_counter() - job_start,
        n_points=len(rows), n_pass=sum(r.get("status") == "pass" for r in rows),
        n_fail=sum(r.get("status") == "FAIL" for r in rows),
        n_skipped=sum(r.get("status") == "skipped_budget" for r in rows),
    )
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    log(f"done: {summary['n_pass']} pass, {summary['n_fail']} fail, "
        f"{summary['n_skipped']} skipped; Selene {summary['selene_seconds_total']:.0f}s")
    log_file.close()
    return 0


# Cost model. Calibrated by the smoke run (B, d=3, 300 shots, 4 processes:
# 6.3 us/op single-process equivalent); 8 leaves room for the ~3x instability.
US_PER_OP_DEFAULT = 8.0
USABLE_FRACTION = 0.85


def parallel_speedup(n_processes: int) -> float:
    """interfaces.md: ~3x at 4 processes, flat beyond."""
    return max(1.0, min(n_processes, 4) * 0.75)


def sec_per_shot(ops: int, us_per_op: float, n_processes: int) -> float:
    return ops * us_per_op * 1e-6 / parallel_speedup(n_processes)

if __name__ == "__main__":
    sys.exit(main())
