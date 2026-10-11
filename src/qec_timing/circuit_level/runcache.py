"""A content-addressed cache for simulation results.

Every result is keyed by the full parameter set that produced it, so a point is
computed once and read back thereafter. Change any parameter and you get a new
key and a fresh run; change nothing and the sweep is instant.

    results/circuit_level/cache/<experiment>/<key>/meta.json     what was run, and when
    results/circuit_level/cache/<experiment>/<key>/result.json   what came out

`DATA` defaults to the repo's `results/circuit_level/cache` (gitignored). Assign
`runcache.DATA` to keep runs elsewhere.

    import runcache
    hit = runcache.load("06_selene_interval", params)
    if hit is None:
        hit = compute(...)
        runcache.save("06_selene_interval", params, hit)

The key is a hash of the parameters rather than a human-readable name on purpose:
a name that has to encode ten parameters either grows unbounded or silently
collides. `meta.json` carries the readable version, and `index()` reads them back.
"""

from __future__ import annotations

import datetime
import hashlib
import json
import pathlib


def _repo_root() -> pathlib.Path:
    for base in [pathlib.Path.cwd(), *pathlib.Path.cwd().parents]:
        if (base / ".git").exists():
            return base
    return pathlib.Path.cwd()


DATA = pathlib.Path(__file__).resolve().parents[3] / "results" / "circuit_level" / "cache"


def key(params: dict) -> str:
    """Stable 16-hex-char digest of a parameter dict.

    Floats are formatted rather than repr'd so 0.002 and 2e-3 land on one key.
    """
    norm = {k: (f"{v:.12g}" if isinstance(v, float) else v) for k, v in params.items()}
    blob = json.dumps(norm, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha1(blob.encode()).hexdigest()[:16]


def path(experiment: str, params: dict) -> pathlib.Path:
    return DATA / experiment / key(params)


def load(experiment: str, params: dict):
    """The stored result, or None."""
    p = path(experiment, params) / "result.json"
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text())
    except json.JSONDecodeError:
        return None          # a half-written file is a miss, not a crash


_ENV: dict | None = None


def environment() -> dict:
    """Provenance: what code produced a result, on what.

    Captured once per process. Without this a cached number is unfalsifiable --
    you cannot tell whether it predates a bug fix, and a stale entry looks
    exactly like a fresh one.
    """
    global _ENV
    if _ENV is not None:
        return _ENV
    import platform
    import subprocess

    def _git(*args):
        try:
            return subprocess.run(["git", "-C", str(_repo_root()), *args],
                                  capture_output=True, text=True,
                                  timeout=10).stdout.strip() or None
        except Exception:
            return None

    from importlib.metadata import PackageNotFoundError, version
    versions = {}
    for dist in ("selene-sim", "selene-hugr-qis-compiler", "guppylang", "stim",
                 "pymatching", "numpy"):
        try:
            versions[dist] = version(dist)
        except PackageNotFoundError:
            pass

    _ENV = {
        "git_commit": _git("rev-parse", "--short", "HEAD"),
        "git_dirty": bool(_git("status", "--porcelain")),
        "python": platform.python_version(),
        "host": platform.node(),
        "packages": versions,
    }
    return _ENV


def save(experiment: str, params: dict, result: dict,
         extra: dict | None = None) -> pathlib.Path:
    d = path(experiment, params)
    d.mkdir(parents=True, exist_ok=True)
    (d / "meta.json").write_text(json.dumps(
        {"experiment": experiment,
         "key": key(params),
         "params": params,
         "saved": datetime.datetime.now().isoformat(timespec="seconds"),
         "environment": environment(),
         **({"extra": extra} if extra else {})},
        indent=1, default=str))
    # write result last: meta without result reads as a miss, so an interrupted
    # save never looks like a complete one
    (d / "result.json").write_text(json.dumps(result, indent=1, default=str))
    return d


def save_arrays(experiment: str, params: dict, **arrays) -> pathlib.Path:
    """Store raw arrays (e.g. detection events) next to a run's meta/result.

    Call BEFORE ``save``: ``save`` writes ``result.json`` last, so arrays without a result
    read as a miss and an interrupted run is redone rather than half-trusted.
    """
    import os
    import numpy as np
    d = path(experiment, params)
    d.mkdir(parents=True, exist_ok=True)
    final = d / "arrays.npz"
    if final.exists() and (d / "result.json").exists():
        raise FileExistsError(f"refusing to overwrite a complete stored run: {d}")
    # write to a temporary name, then rename: a crash mid-write never leaves a truncated
    # arrays.npz behind under the real name
    tmp = d / f"arrays.tmp.{os.getpid()}.npz"
    np.savez_compressed(tmp, **arrays)
    os.replace(tmp, final)
    return final


def load_arrays(experiment: str, params: dict):
    """The stored arrays as a dict, or None if the run is absent or incomplete."""
    import numpy as np
    d = path(experiment, params)
    if not ((d / "result.json").exists() and (d / "arrays.npz").exists()):
        return None
    with np.load(d / "arrays.npz") as z:
        return {k: z[k] for k in z.files}


def index(experiment: str | None = None) -> list[dict]:
    """Every cached run, newest first, as flat dicts of params + result."""
    roots = [DATA / experiment] if experiment else (
        [p for p in DATA.iterdir() if p.is_dir()] if DATA.exists() else [])
    rows = []
    for root in roots:
        if not root.exists():
            continue
        for d in root.iterdir():
            meta_p, res_p = d / "meta.json", d / "result.json"
            if not (meta_p.exists() and res_p.exists()):
                continue
            try:
                meta, res = json.loads(meta_p.read_text()), json.loads(res_p.read_text())
            except json.JSONDecodeError:
                continue
            rows.append({"experiment": meta["experiment"], "key": d.name,
                         "saved": meta["saved"], **meta["params"], **res})
    return sorted(rows, key=lambda r: r["saved"], reverse=True)


def summary() -> str:
    rows = index()
    if not rows:
        return f"{DATA}: empty"
    by = {}
    for r in rows:
        by[r["experiment"]] = by.get(r["experiment"], 0) + 1
    parts = ", ".join(f"{k}: {v}" for k, v in sorted(by.items()))
    return f"{DATA}  ({len(rows)} cached runs -- {parts})"
