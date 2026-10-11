"""The effective noise model of one surface-code extraction round, measured in Selene.

The procedure of `04_extraction_cost_surface` (p_1q counts) and `04c_extraction_cost_p2q`
(p_2q counts), packaged so a script can run it for any distance:

    eff = measure("rotated", 5, with_p2q=True)
    eff["ops_per_qubit"], eff["ops_anc"]        # charged 1q ops per round, data / Z-ancilla
    eff["n2_per_qubit"],  eff["n2_anc"]         # Rzz events per round reaching each

The proxy is `|0>^n` with both sectors extracted and no `h` on the X-ancillas, so the round
is logically inert while every CNOT is still compiled and charged. A per-qubit flip rate
inverts to an operation count, and the count per round is the slope in R.

Runs are cached under the SAME experiment names and parameter sets the notebooks use, so a
distance the notebooks already measured costs nothing here and gives identical numbers.
"""

from __future__ import annotations

import importlib.util
import json
import pathlib
import sys
import tempfile

import numpy as np

from . import runcache
from .surface_memory import PREAMBLE, SPEC, layout

ROUNDS = (1, 2, 4, 6, 8)
SHOTS = 12000
SEED = 1
P1Q_SLOPE, P1Q_ANC = 0.001, 0.006      # 04_extraction_cost_surface
P2Q_SLOPE, P2Q_ANC = 0.001, 0.006      # 04c_extraction_cost_p2q
ANC_SHOTS = 6000
PROXY = "zero_state_both_sectors_no_h"


def f1(p_1q):
    """Z-basis flip per charged 1q op: X and Y flip |0>, Z does not."""
    return 2 * p_1q / 3


def f2(p_2q):
    """Flip of one qubit per Rzz: 8 of the 15 two-qubit Paulis carry X or Y on it."""
    return 8 * p_2q / 15


def invert(rate, f):
    """Flips compose mod 2: rate = (1 - (1-2f)^n)/2, solved for n."""
    r = np.clip(np.asarray(rate, float), 0, 0.4999)
    return np.log(1 - 2 * r) / np.log(1 - 2 * f)


class Proxy:
    """The inert extraction circuit for one code, built once per round count."""

    def __init__(self, code: str, d: int):
        self.code, self.d = code, d
        self.L = L = layout(SPEC[code], d)
        self.didx = {q: i for i, q in enumerate(L["data"])}
        self.aidx = {q: i for i, q in enumerate(L["zanc"] + L["xanc"])}
        self._gen = pathlib.Path(tempfile.mkdtemp(prefix=f"proxy_{code}_d{d}_"))
        self._built = {}

    def source(self, rounds: int) -> str:
        L, didx, aidx = self.L, self.didx, self.aidx
        o = ["@guppy", "def prog() -> None:",
             f"    qs = array(qubit() for _ in range(py({len(didx)})))",
             f"    for _ in range(py({rounds})):",
             f"        anc = array(qubit() for _ in range(py({len(aidx)})))"]
        for ctrl, tgt in L["cx"]:
            cs = f"anc[{aidx[ctrl]}]" if ctrl in aidx else f"qs[{didx[ctrl]}]"
            ts = f"anc[{aidx[tgt]}]" if tgt in aidx else f"qs[{didx[tgt]}]"
            o.append(f"        cx({cs}, {ts})")
        o.append('        output("s", collect_measurements(measure_array(anc)))')
        o.append('    output("data", collect_measurements(measure_array(qs)))')
        return "\n".join(o)

    def _runner(self, rounds: int):
        if rounds not in self._built:
            from selene_sim import build
            name = f"proxy_{self.code}_d{self.d}_r{rounds}"
            path = self._gen / f"{name}.py"
            # Guppy reads the decorated function's source from disk; exec() of a string fails.
            path.write_text(PREAMBLE + self.source(rounds) + "\n")
            spec = importlib.util.spec_from_file_location(name, path)
            mod = importlib.util.module_from_spec(spec)
            sys.modules[name] = mod
            spec.loader.exec_module(mod)
            self._built[rounds] = build(mod.prog.compile())
        return self._built[rounds]

    def run(self, rounds: int, shots: int, **ch):
        """(syndromes [shots, rounds, n_anc], data [shots, n_data])."""
        from selene_sim import DepolarizingErrorModel, Stim
        syn, dat = [], []
        # Each yielded shot is a lazy generator inside the `with TCPStream` block: consume
        # it as it arrives, or the child dies with ConnectionRefused.
        for shot in self._runner(rounds).run_shots(
                simulator=Stim(), n_qubits=len(self.didx) + len(self.aidx), n_shots=shots,
                error_model=DepolarizingErrorModel(random_seed=SEED, **ch), random_seed=SEED):
            pairs = list(shot)
            syn.append([[int(b) for b in v] for k, v in pairs if k == "s"])
            dat.append([int(b) for b in [v for k, v in pairs if k == "data"][0]])
        return np.array(syn), np.array(dat)

    def rates(self, rounds: int, shots: int, **ch) -> dict:
        """Per-data-qubit flip rates, and round-0 Z-ancilla rates conditioned on clean data."""
        syn, dat = self.run(rounds, shots, **ch)
        raw, local = [], []
        for a in self.L["zanc"]:
            col = syn[:, 0, self.aidx[a]]
            clean = dat[:, [self.didx[q] for q in self.L["zsup"][a]]].sum(axis=1) == 0
            raw.append(float(col.mean()))
            local.append(float(col[clean].mean()))
        return {"data_rates": dat.mean(axis=0).tolist(), "z_raw": raw, "z_local": local,
                "n_clean": int(clean.sum())}


def _sample_p1q(px: Proxy, rounds: int, p_1q: float, shots: int) -> dict:
    params = dict(code=px.code, d=px.d, rounds=rounds, p_1q=p_1q, p_2q=0.0,
                  shots=shots, seed=SEED, proxy=PROXY)
    hit = runcache.load("04_surface_extraction", params)
    if hit is None:
        hit = px.rates(rounds, shots, p_1q=p_1q)
        runcache.save("04_surface_extraction", params, hit)
    return hit


def _sample_p2q(px: Proxy, rounds: int, p_2q: float, shots: int) -> dict:
    params = dict(code=px.code, d=px.d, rounds=rounds, tau=0, p_1q=0.0, p_2q=p_2q,
                  shots=shots, seed=SEED, proxy=PROXY)
    hit = runcache.load("04c_extraction_cost_p2q", params)
    if hit is None:
        hit = px.rates(rounds, shots, p_2q=p_2q)
        hit.pop("n_clean")
        runcache.save("04c_extraction_cost_p2q", params, hit)
    return hit


def _slope(px: Proxy, sampler, p: float, f) -> dict:
    """Per-qubit count per round (slope in R), with the checks that make it trustworthy."""
    R = np.array(ROUNDS, float)
    hits = {r: sampler(px, r, p, SHOTS) for r in ROUNDS}
    n = np.array([invert(hits[r]["data_rates"], f(p)) for r in ROUNDS])
    slope, icept = np.array([np.polyfit(R, n[:, i], 1) for i in range(n.shape[1])]).T
    se = icept.std() / np.sqrt(len(icept))
    return {"per_qubit": slope,
            "intercept_t": float(abs(icept.mean()) / se) if se > 0 else 0.0,
            "max_rate": float(max(hits[ROUNDS[-1]]["data_rates"]))}


def measure(code: str, d: int, with_p2q: bool = False, log=print) -> dict:
    """Measure (or read from cache) the effective noise model for one code."""
    px = Proxy(code, d)
    log(f"  effective noise model, {code} d={d}: p_1q counts")
    s1 = _slope(px, _sample_p1q, P1Q_SLOPE, f1)
    anc1 = _sample_p1q(px, 1, P1Q_ANC, ANC_SHOTS)
    eff = {"code": code, "d": d,
           "ops_per_qubit": s1["per_qubit"],
           "ops_anc": float(invert(np.mean(anc1["z_local"]), f1(P1Q_ANC))),
           "checks": {"p1q_intercept_in_se": round(s1["intercept_t"], 2),
                      "p1q_max_rate": round(s1["max_rate"], 4)}}
    if with_p2q:
        log(f"  effective noise model, {code} d={d}: p_2q counts")
        s2 = _slope(px, _sample_p2q, P2Q_SLOPE, f2)
        anc2 = _sample_p2q(px, 1, P2Q_ANC, ANC_SHOTS)
        eff["n2_per_qubit"] = s2["per_qubit"]
        eff["n2_anc"] = float(invert(np.mean(anc2["z_local"]), f2(P2Q_ANC)))
        eff["checks"].update({"p2q_intercept_in_se": round(s2["intercept_t"], 2),
                              "p2q_max_rate": round(s2["max_rate"], 4)})
    for k, t in eff["checks"].items():
        if k.endswith("intercept_in_se") and t > 4:
            log(f"  WARNING {code} d={d}: {k} = {t} -- something outside the round contributes")
    return eff


def save(eff: dict, directory: pathlib.Path) -> pathlib.Path:
    """Write the model as JSON, with lambda = 1/n_1 and the b_read floor alongside."""
    directory.mkdir(parents=True, exist_ok=True)
    ops = np.asarray(eff["ops_per_qubit"])
    out = {k: (np.round(v, 4).tolist() if isinstance(v, np.ndarray) else v)
           for k, v in eff.items()}
    out["ops_mean"] = round(float(ops.mean()), 4)
    out["lambda_per_idle_op_at_p2q_0"] = round(1.0 / float(ops.mean()), 4)
    out["b_read_floor_at_p2q_0"] = round(eff["ops_anc"] / float(ops.mean()), 4)
    if "n2_per_qubit" in eff:
        out["n2_mean"] = round(float(np.mean(eff["n2_per_qubit"])), 4)
        out["lambda_formula"] = "k / (ops_mean + 0.8 * (p_2q/p_1q) * n2_mean)"
    out["units"] = ("ops_* in charged 1q ops (flip 2p_1q/3 each); n2_* in Rzz events "
                    "(flip 8p_2q/15 each); all per extraction round")
    path = directory / f"{eff['code']}_d{eff['d']}.json"
    path.write_text(json.dumps(out, indent=1))
    return path
