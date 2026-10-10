"""Stage 7a: dt* under H2-1E noise, estimated with the Stage 6c twin (main venv).

The H2 idle model comes from ``src/qec_timing/h2/noise_params.py``, loaded BY
FILE PATH so the main venv never imports the ``qec_timing.h2`` package (which
also holds Guppy-0.x code). Twin: variant B (dephasing, X memory),
``order="hook_safe"``, ``native_charging=True``, H2-1E p1 / p2 / asymmetric
p_meas / p_init. Time unit = one round (t_round); dt is the syndrome PERIOD,
so dt >= 1 round.

Approximations, beyond the twin's own (see docs/stage7a_report.md):
* H2 charges no error on Rz (virtual); the twin's native charging charges p1 on
  every Rz. With p1 = 1.89e-5 this changes very little next to p2 = 1.05e-3.
* Not modelled: emission, crosstalk, leakage, ancilla dephasing during
  extraction.

This is an estimate (a few hundred failures per point), not a sweep: it fits
dt* with the Stage 6c per-d shape fit and a binomial bootstrap.

Usage: python scripts/stage7a_dtstar.py [--t-round 2e-3] [--distances 3 5]
"""

from __future__ import annotations

import argparse
import datetime as _dt
import importlib.util
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_timing.circuit_level import (  # noqa: E402
    CircuitNoise,
    DephasingH2,
    ReadoutNoise,
    build_constant_interval_memory,
    rate_per_unit_time,
)


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


h2n = _load("h2_noise_params", ROOT / "src" / "qec_timing" / "h2" / "noise_params.py")
s6c = _load("stage6c_sweep", ROOT / "scripts" / "stage6c_sweep.py")

TARGET_FAILURES = 400
NMAX = 60
SEED = 70_000_000
STRIDE = 100_000


def twin_noise(f_quad: float, f_lin: float) -> CircuitNoise:
    n = h2n.H2_1E
    return CircuitNoise(DephasingH2(f_quad, f_lin), p1=n.p1, p2=n.p2, p_init=n.p_init,
                        readout=ReadoutNoise(n.p_meas_0, n.p_meas_1))


def measure(d: int, dt: float, noise: CircuitNoise, block: int) -> dict:
    pilot = build_constant_interval_memory(d, dt, 2 * d * dt, noise, basis="X",
                                           order="hook_safe", native_charging=True)
    dec = pilot.decoder()
    f = sh = b = 0
    t0 = time.perf_counter()
    while f < 20 and sh < 100_000 and time.perf_counter() - t0 < 3:
        dets, obs = pilot.sample(5_000, seed=block + b)
        f += int(np.count_nonzero(np.any(dec.decode_batch(dets) != obs, axis=1)))
        sh += 5_000
        b += 1
    eps = rate_per_unit_time(min(max(f, 1) / sh, 0.45), pilot.schedule.total_time) * dt
    n = int(min(max(math.ceil(-math.log(0.6) / (2 * eps)), 2 * d), NMAX))
    mem = build_constant_interval_memory(d, dt, n * dt, noise, basis="X",
                                         order="hook_safe", native_charging=True)
    dec = mem.decoder()
    p_pred = 0.5 * (1 - math.exp(-2 * eps * n))
    shots = max(2_000, math.ceil(1.15 * TARGET_FAILURES / max(p_pred, 1e-9)))
    fails, done = s6c.run_batches(mem, dec, shots, block + 1_000, target=TARGET_FAILURES)
    T = mem.schedule.total_time
    p = fails / done
    return dict(d=d, dt_rounds=dt, n_rounds=mem.n_rounds, shots=done, failures=fails,
                p_L=p, T_rounds=T, R_per_round=rate_per_unit_time(min(p, 0.4999), T),
                saturated=p > s6c.P_L_SATURATED, seed_block=[block, block + STRIDE])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--t-round", type=float, default=h2n.DEFAULT_T_ROUND, help="seconds")
    ap.add_argument("--distances", type=int, nargs="+", default=[3, 5])
    ap.add_argument("--conventions", nargs="+", default=list(h2n.CONVENTIONS))
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "stage7a")
    args = ap.parse_args(argv)

    out = args.out / (_dt.datetime.now(_dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-dtstar")
    out.mkdir(parents=True, exist_ok=False)
    summary = {"t_round_s": args.t_round, "noise": h2n.H2_1E.__dict__, "runs": {}}
    rng = np.random.default_rng(7)
    block_i = 0
    for conv in args.conventions:
        f_quad, f_lin = h2n.twin_dephasing(args.t_round, convention=conv)
        noise = twin_noise(f_quad, f_lin)
        # dt* ~ 1 / f_quad-ish scale; span it generously, period >= 1 round.
        hi = min(20.0 / f_quad, 1e5)
        grid = sorted({max(1.0, round(x, 3)) for x in np.geomspace(1.0, hi, 16)})
        print(f"[{conv}] f_quad={f_quad:.4e} rad/round, f_lin={f_lin:.3e} /round; "
              f"p_z(1 round)={h2n.idle_pz(args.t_round, convention=conv):.3e}", flush=True)
        for d in args.distances:
            pts = []
            for dt in grid:
                pt = measure(d, dt, noise, SEED + block_i * STRIDE)
                block_i += 1
                pts.append(pt)
                print(f"   d={d} dt={dt:g} rounds: R/round={pt['R_per_round']:.3e} "
                      f"({pt['failures']}/{pt['shots']}, n={pt['n_rounds']})"
                      f"{' SAT' if pt['saturated'] else ''}", flush=True)
            ok = [p for p in pts if not p["saturated"] and p["failures"] > 0]
            x = np.array([p["dt_rounds"] for p in ok])
            k = np.array([p["failures"] for p in ok]); n = np.array([p["shots"] for p in ok])
            T = np.array([p["T_rounds"] for p in ok])

            def fit(kk):
                pp = np.clip(kk / n, 1e-12, 0.4999)
                R = s6c.R_of(pp, T)
                sd = np.sqrt(pp * (1 - pp) / n) / (T * (1 - 2 * pp) * R)
                return s6c.shape_min(x, R, sd, "B")

            best = fit(k)
            boots = [fit(rng.binomial(n, k / n)) for _ in range(300)]
            bt = np.array([b[0] for b in boots if b is not None])
            res = dict(points=pts)
            if best:
                res.update(dtstar_rounds=best[0],
                           dtstar_rounds_ci95=[float(np.percentile(bt, 2.5)),
                                               float(np.percentile(bt, 97.5))],
                           dtstar_seconds=best[0] * args.t_round,
                           R_star_per_round=best[1],
                           R_at_1_round=next(p["R_per_round"] for p in pts if p["dt_rounds"] == 1.0))
                res["gain_vs_every_round"] = res["R_at_1_round"] / best[1]
                print(f"   -> d={d} [{conv}]: dt* = {best[0]:.1f} rounds "
                      f"[{res['dtstar_rounds_ci95'][0]:.1f}, {res['dtstar_rounds_ci95'][1]:.1f}] "
                      f"= {res['dtstar_seconds'] * 1e3:.0f} ms; R(1 round)/R(dt*) = "
                      f"{res['gain_vs_every_round']:.2f}", flush=True)
            summary["runs"][f"{conv}_d{d}"] = dict(convention=conv, d=d, f_quad=f_quad,
                                                    f_lin=f_lin, **res)
    (out / "dtstar.json").write_text(json.dumps(summary, indent=2))
    print(f"wrote {out / 'dtstar.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
