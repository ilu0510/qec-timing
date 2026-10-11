#!/usr/bin/env python
"""Generate, store and decode Selene memory shots for the 06 interval study. Cluster-ready.

Each block is one Selene run (``surface_memory.Geometry``) of ``--shots`` shots at a fixed
(tau, rounds, noise, seed, shot_offset). Blocks are stored as raw detection events by
``memory_shots`` (never overwritten: a complete block is loaded, not regenerated) and decoded by
both decoders into ``selene_memory_decoded``. Everything lands under
``results/circuit_level/cache/``, keyed by every parameter including the distance.

Work is split by ``--worker i --workers n``: worker i takes every n-th job of a fixed,
cost-sorted list, so n independent processes (or a cluster job array) never share a block.
Store directories from different machines merge by copying, since keys never collide.

    # check the Stim copy at this distance first (exit code 1 if it fails)
    python scripts/circuit_level_generate.py --distance 5 --verify-only
    # one worker of 16
    python scripts/circuit_level_generate.py --distance 5 --blocks 46 --worker 3 --workers 16
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time

REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import numpy as np                                             # noqa: E402

from qec_timing.circuit_level import selene_mirror as sm       # noqa: E402
from qec_timing.circuit_level.surface_memory import Geometry   # noqa: E402


def parse():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--code", default="rotated")
    ap.add_argument("--distance", type=int, required=True)
    ap.add_argument("--total-idle", type=int, default=96, help="R = total / tau")
    ap.add_argument("--taus", type=int, nargs="*", default=None,
                    help="default: every divisor of total with R >= d")
    ap.add_argument("--p1q", type=float, default=0.0005)
    ap.add_argument("--pmeas", type=float, default=0.002)
    ap.add_argument("--shots", type=int, default=10_000, help="shots per block")
    ap.add_argument("--blocks", type=int, default=5, help="blocks per tau (shot offsets k * shots)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--worker", type=int, default=0)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--verify-only", action="store_true",
                    help="run the 06a checks of the Stim copy at this distance and exit")
    return ap.parse_args()


def verify(G, out_path):
    """06a's checks at this distance: noiseless determinism and statistics against Selene."""
    res = {"code": G.code, "distance": G.d}
    for tau, R in ((1, 3), (8, 12)):
        c0 = sm.mirror_circuit(G, tau, R, 0.0, 0.0)
        d0, o0 = c0.compile_detector_sampler().sample(2000, separate_observables=True)
        res[f"noiseless tau={tau} R={R}"] = {"bits": int(d0.sum()), "observable": int(o0.sum())}
    tau, R, n, p1, pm = 8, 12, 20_000, 0.002, 0.004
    cm = sm.mirror_circuit(G, tau, R, p1, pm)

    def stats(A, B):
        pA, pB = A.mean(0), B.mean(0)
        z = (pA - pB) / (np.sqrt((pA * (1 - pA) + pB * (1 - pB)) / n) + 1e-12)
        iu = np.triu_indices(A.shape[1], 1)
        cA = (A.T.astype(np.int64) @ A.astype(np.int64))[iu] / n
        cB = (B.T.astype(np.int64) @ B.astype(np.int64))[iu] / n
        k = (cA + cB) * n / 2 >= 20
        zc = (cA - cB)[k] / np.sqrt((cA * (1 - cA) + cB * (1 - cB))[k] / n)
        return {"per-detector chi2/dof": float(np.mean(z ** 2)), "pairwise chi2/dof": float(np.mean(zc ** 2)),
                "pairs": int(k.sum())}

    blocks = []
    for off in (0, n):
        s, d = G.run(tau, R, n, offset=off, seed=11, p_1q=p1, p_meas=pm)
        blocks.append(G.detectors(s, d).astype(bool))
    M = [cm.compile_detector_sampler(seed=sd).sample(n, separate_observables=True)[0] for sd in (21, 22)]
    res["statistics"] = {"Selene vs Selene": stats(blocks[0], blocks[1]),
                         "copy vs copy": stats(M[0], M[1]),
                         "Selene vs copy 1": stats(blocks[0], M[0]),
                         "Selene vs copy 2": stats(blocks[1], M[1])}
    st = res["statistics"]
    refs, tests = ("Selene vs Selene", "copy vs copy"), ("Selene vs copy 1", "Selene vs copy 2")
    ref = max(st[k]["per-detector chi2/dof"] for k in refs)
    refp = max(st[k]["pairwise chi2/dof"] for k in refs)
    worst = max(st[k]["per-detector chi2/dof"] for k in tests)
    worstp = max(st[k]["pairwise chi2/dof"] for k in tests)
    noiseless_ok = all(v["bits"] == 0 and v["observable"] == 0 for k, v in res.items() if k.startswith("noiseless"))
    # pass: noiseless exact, and Selene-vs-copy within 1.5x of the larger reference
    res["pass"] = bool(noiseless_ok and worst <= 1.5 * ref and worstp <= 1.5 * refp)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1), flush=True)
    return res["pass"]


def main():
    a = parse()
    G = Geometry(a.code, a.distance)
    if a.verify_only:
        ok = verify(G, REPO / "results" / "circuit_level" / f"verify_mirror_{a.code}_d{a.distance}.json")
        sys.exit(0 if ok else 1)
    taus = a.taus or [t for t in range(1, a.total_idle + 1)
                      if a.total_idle % t == 0 and a.total_idle // t >= a.distance]
    jobs = [(t, k) for t in taus for k in range(a.blocks)]
    jobs.sort(key=lambda j: (-(a.total_idle // j[0]), j[0], j[1]))      # costliest first, fixed order
    mine = jobs[a.worker::a.workers]
    print(f"worker {a.worker}/{a.workers}: {len(mine)} of {len(jobs)} blocks, d={a.distance}, "
          f"taus={taus}, {a.shots} shots/block", flush=True)
    for tau, k in mine:
        t0 = time.time()
        r = sm.decoded(G, tau, a.total_idle // tau, a.shots, offset=k * a.shots, seed=a.seed,
                       p_1q=a.p1q, p_meas=a.pmeas)
        print(f"d={a.distance} tau={tau} R={a.total_idle // tau} block={k}: hand {r['hand_failures']} "
              f"mirror {r['mirror_failures']} / {r['shots']} [{time.time() - t0:.0f}s]", flush=True)
    print(f"worker {a.worker} done", flush=True)


if __name__ == "__main__":
    main()
