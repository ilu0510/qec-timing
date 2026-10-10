"""Stage 7a: build the H2 test program and convert it to QIR. Run in .venv-h2 ONLY.

Guppy 0.21 -> HUGR (``.compile()``) -> QIR bitcode (``hugr_qir.hugr_to_qir``,
``validate_qir=True``). Writes ``<name>.bc`` (what Nexus takes), ``<name>.ll``
(readable IR) and ``<name>.json`` (parameters, sizes, timings, op counts).

Usage (from the repo root):
    .venv-h2/Scripts/python.exe scripts/stage7a_build_qir.py                 # d=3, 2 rounds, Z
    .venv-h2/Scripts/python.exe scripts/stage7a_build_qir.py --checks X
    .venv-h2/Scripts/python.exe scripts/stage7a_build_qir.py --scaling       # size/time vs d, rounds
"""

from __future__ import annotations

import argparse
import datetime as _dt
import importlib.metadata as md
import json
import re
import sys
import threading
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from hugr_qir.hugr_to_qir import hugr_to_qir  # noqa: E402
from hugr_qir.output import OutputFormat  # noqa: E402

from qec_timing.h2.program import build_program, effective_p_inject  # noqa: E402

H2_1E_STABILIZER_QUBITS = 56  # docs: H2-1E "Qubits (Stabilizer) 56"

# quantinuum_qircheck (run by validate_qir=True) checks for cycles by recursing
# once per basic block. Every RNG-conditioned Z adds blocks, so beyond ~600
# injections (e.g. d=5, 32 rounds) it hits Python's default recursion limit.
# Workaround with validation still ON: a high recursion limit, inside a worker
# thread with a large stack (the Windows default thread stack would overflow).
RECURSION_LIMIT = 1_000_000
THREAD_STACK = 1 << 26  # 64 MiB: Windows rejects >= 256 MiB; Python 3.11 frames need little C stack


def _qir(pkg, fmt):
    box = {}

    def work():
        sys.setrecursionlimit(RECURSION_LIMIT)
        try:
            box["out"] = hugr_to_qir(pkg, validate_qir=True, output_format=fmt)
        except BaseException as e:  # noqa: BLE001  re-raised in the caller
            box["err"] = e

    old = threading.stack_size(THREAD_STACK)
    try:
        t = threading.Thread(target=work)
        t.start()
        t.join()
    finally:
        threading.stack_size(old)
    if "err" in box:
        raise box["err"]
    return box["out"]


def build(d: int, rounds: int, checks: str, p_inject: float, base_seed: int) -> dict:
    t0 = time.perf_counter()
    pkg = build_program(d, rounds, p_inject, base_seed, checks).compile()
    t_hugr = time.perf_counter() - t0
    t1 = time.perf_counter()
    bc = _qir(pkg, OutputFormat.BITCODE)
    t_qir = time.perf_counter() - t1
    ll = _qir(pkg, OutputFormat.LLVM_IR)
    calls = Counter(re.findall(r"call [^@]*@([\w.]+)", ll))
    return dict(d=d, rounds=rounds, checks=checks, p_inject=p_inject,
                p_inject_effective=effective_p_inject(p_inject), base_seed=base_seed,
                n_qubits=d * d + 1, hugr_bytes=len(pkg.to_bytes()), qir_bitcode_bytes=len(bc),
                qir_ir_chars=len(ll), seconds_hugr=t_hugr, seconds_qir=t_qir,
                fits_h2_1e_stabilizer=d * d + 1 <= H2_1E_STABILIZER_QUBITS,
                qir_calls=dict(calls), _bc=bc, _ll=ll)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--d", type=int, default=3)
    ap.add_argument("--rounds", type=int, default=2)
    ap.add_argument("--checks", choices=("Z", "X"), default="Z")
    ap.add_argument("--p-inject", type=float, default=0.05)
    ap.add_argument("--base-seed", type=int, default=7_000_000)
    ap.add_argument("--scaling", action="store_true")
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "stage7a")
    args = ap.parse_args(argv)

    versions = {p: md.version(p) for p in ("guppylang", "guppylang-internals", "hugr",
                                           "hugr-qir", "pyqir", "quantinuum-qircheck")}
    if args.scaling:
        out = args.out / (_dt.datetime.now(_dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
                          + "-qir-scaling")
        out.mkdir(parents=True)
        rows = []
        for d, rounds in [(3, 2), (3, 8), (3, 32), (3, 64), (5, 2), (5, 8), (5, 32),
                          (7, 2), (7, 8), (7, 16)]:
            r = build(d, rounds, "X", args.p_inject, args.base_seed)
            row = {k: v for k, v in r.items() if not k.startswith("_") and k != "qir_calls"}
            row["n_calls"] = sum(r["qir_calls"].values())
            rows.append(row)
            print(f"d={d} rounds={rounds}: qubits={row['n_qubits']} hugr={row['hugr_bytes']} B "
                  f"bc={row['qir_bitcode_bytes']} B calls={row['n_calls']} "
                  f"t_hugr={row['seconds_hugr']:.1f}s t_qir={row['seconds_qir']:.1f}s", flush=True)
        (out / "scaling.json").write_text(json.dumps({"versions": versions, "rows": rows}, indent=2))
        print(f"wrote {out / 'scaling.json'}")
        return 0

    r = build(args.d, args.rounds, args.checks, args.p_inject, args.base_seed)
    name = f"d{args.d}_r{args.rounds}_{args.checks}"
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / f"{name}.bc").write_bytes(r["_bc"])
    (args.out / f"{name}.ll").write_text(r["_ll"])
    meta = {k: v for k, v in r.items() if not k.startswith("_")}
    meta["versions"] = versions
    (args.out / f"{name}.json").write_text(json.dumps(meta, indent=2))
    print(f"{name}: validated QIR, {meta['qir_bitcode_bytes']} B bitcode -> "
          f"{args.out / (name + '.bc')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
