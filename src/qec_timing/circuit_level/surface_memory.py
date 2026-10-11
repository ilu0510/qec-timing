"""Shared machinery for surface-code memory experiments in Selene.

Factored out of `06_selene_interval.ipynb` so that `07_memory_error_rate.ipynb` runs the
*same* circuit, the same CNOT order and the same decoder wiring. Two notebooks that differ
only in how they sweep must not differ in what they simulate.

    L = layout("surface_code:rotated_memory_z", 3)
    G = Geometry(L)
    syn, dat = G.run(tau=8, rounds=12, shots=1000, p_1q=0.001)
    pL = G.decode(syn, dat, tau=8, rounds=12, p_1q=0.001, p_meas=0.002)

`06` retains its own inline copy; this module is not retrofitted to it, so a change here
does not silently alter a notebook that has already been run and recorded.
"""

from __future__ import annotations

import importlib.util
import pathlib
import sys
import tempfile
from collections import defaultdict

import numpy as np
import pymatching
import stim
import yaml

SPEC = {"rotated": "surface_code:rotated_memory_z",
        "unrotated": "surface_code:unrotated_memory_z"}
# X-basis memory (|+_L>, decoded on X-checks): the only basis that sees Z (dephasing) noise.
SPEC_X = {"rotated": "surface_code:rotated_memory_x",
          "unrotated": "surface_code:unrotated_memory_x"}

PREAMBLE = """from guppylang import guppy
from guppylang.std.angles import angle
from guppylang.std.builtins import array, output, py
from guppylang.std.quantum import collect_measurements, cx, h, measure, measure_array, qubit
from guppylang.std.qsystem import phased_x


"""


CONFIG = pathlib.Path(__file__).resolve().parents[3] / "config"   # the repo's config/


def find_config(name: str):
    p = CONFIG / name
    return p if p.exists() else None


def _config_constants(code: str, d: int):
    """Measured charged-op counts for this code, from 04_extraction_cost_surface.

    Returns (ops_per_data_qubit array, mean, ancilla_z). Raises if the measured file is
    absent rather than falling back to the derived one: the derived surface rows price the
    CONTROL side of every CNOT, while a surface data qubit is the TARGET of each X-check,
    which understates the data side by ~2.6x and lands straight on the decoder priors.
    """
    p = find_config("surface_extraction_measured.yaml")
    if p is None:
        raise FileNotFoundError(
            "config/surface_extraction_measured.yaml missing -- "
            "run 04_extraction_cost_surface.ipynb first")
    m = yaml.safe_load(p.read_text())["surface_measured"][f"{code}_d{d}"]
    ops = np.array(m["ops_per_round"]["data_per_qubit"], float)
    return ops, float(m["ops_per_round"]["data_mean"]), float(m["ops_per_round"]["ancilla_z"])


def layout(spec: str, d: int) -> dict:
    """Data qubits, Z/X stabiliser supports in CNOT order, and the logical support.

    The CNOT *order* is taken from Stim rather than written by hand: the order decides how
    an ancilla fault propagates back onto the data.
    """
    c = stim.Circuit.generated(spec, distance=d, rounds=1)
    anc, data_order, cx_pairs, obs = set(), [], [], []
    for inst in c.flattened():
        t = [g.value for g in inst.targets_copy() if g.is_qubit_target]
        if inst.name in ("CX", "CNOT", "ZCX"):
            cx_pairs += [(t[i], t[i + 1]) for i in range(0, len(t), 2)]
        elif inst.name in ("MR", "MRZ"):
            anc.update(t)
        elif inst.name in ("M", "MZ", "MX"):
            data_order = t
        elif inst.name == "OBSERVABLE_INCLUDE":
            recs = [g.value for g in inst.targets_copy() if g.is_measurement_record_target]
            obs = [data_order[len(data_order) + r] for r in recs]
    data = [q for q in data_order if q not in anc]
    zsup, xsup = defaultdict(list), defaultdict(list)
    for ctrl, tgt in cx_pairs:
        if tgt in anc:
            zsup[tgt].append(ctrl)
        if ctrl in anc:
            xsup[ctrl].append(tgt)
    return dict(data=data, zanc=sorted(zsup), xanc=sorted(xsup),
                zsup=dict(zsup), xsup=dict(xsup), obs=obs, cx=cx_pairs)


def accumulate(n, f):
    """P(odd number of flips) after n charged operations."""
    return (1 - (1 - 2 * f) ** n) / 2


class Geometry:
    """One surface-code memory: program generation, execution and decoding."""

    def __init__(self, code: str, d: int, constants: dict | None = None, basis: str = "Z"):
        """`constants` overrides the config file, for distances it does not cover:

            {"ops_per_qubit": array, "ops_anc": float,            # p_1q counts
             "n2_per_qubit": array, "n2_anc": float}              # p_2q counts, optional

        Without it, the p_1q counts come from surface_extraction_measured.yaml and
        p_2q priors are unavailable.

        `basis="X"` prepares |+_L> (h on every data qubit after allocation), reads the data out
        in the X basis (h before the final measurement) and decodes on the X-checks. Its
        program, CNOT order and logical support come from Stim's ``*_memory_x`` circuit. The
        hand-built decoder (`priors`, `matching`) is calibrated for the Z sector only and is
        not available in this basis; decode with ``selene_mirror`` instead.
        """
        if basis not in ("Z", "X"):
            raise ValueError(f"basis must be 'Z' or 'X', not {basis!r}")
        self.code, self.d, self.basis = code, d, basis
        self.L = L = layout((SPEC if basis == "Z" else SPEC_X)[code], d)
        self.n_data = len(L["data"])
        self.n_anc = len(L["zanc"]) + len(L["xanc"])
        self.didx = {q: i for i, q in enumerate(L["data"])}
        self.aidx = {q: i for i, q in enumerate(L["zanc"] + L["xanc"])}
        self.zcols = [self.aidx[a] for a in L["zanc"]]
        self.zsup = [[self.didx[q] for q in L["zsup"][a]] for a in L["zanc"]]
        self.xcols = [self.aidx[a] for a in L["xanc"]]
        self.xsup = [[self.didx[q] for q in L["xsup"][a]] for a in L["xanc"]]
        # the checks this memory is decoded on: detectors, final syndrome, mirror layout
        self.cols, self.sup = (self.zcols, self.zsup) if basis == "Z" else (self.xcols, self.xsup)
        self.obs = [self.didx[q] for q in L["obs"]]
        self.n2_per_qubit = self.n2_anc = None
        if basis == "X":
            self.ops_per_qubit = self.ops_mean = self.ops_anc = None
        elif constants is None:
            self.ops_per_qubit, self.ops_mean, self.ops_anc = _config_constants(code, d)
        else:
            self.ops_per_qubit = np.asarray(constants["ops_per_qubit"], float)
            self.ops_mean = float(self.ops_per_qubit.mean())
            self.ops_anc = float(constants["ops_anc"])
            if constants.get("n2_per_qubit") is not None:
                self.n2_per_qubit = np.asarray(constants["n2_per_qubit"], float)
                self.n2_anc = float(constants["n2_anc"])
        assert self.ops_per_qubit is None or len(self.ops_per_qubit) == self.n_data, \
            "measured op counts do not match this layout's data-qubit count"
        self._gen = pathlib.Path(tempfile.mkdtemp(prefix=f"guppy_{code}_d{d}_"))
        self._built = {}

    # ---- program ---------------------------------------------------------------
    def source(self, tau: int, rounds: int) -> str:
        """One round is: tau idle ops on every data qubit, then BOTH check sectors."""
        L, didx, aidx = self.L, self.didx, self.aidx
        o = ["@guppy", "def prog() -> None:",
             f"    qs = array(qubit() for _ in range(py({self.n_data})))"]
        if self.basis == "X":
            o += [f"    for i in range(py({self.n_data})):", "        h(qs[i])"]
        o.append(f"    for _ in range(py({rounds})):")
        if tau:
            o += [f"        for _ in range(py({tau})):",
                  f"            for i in range(py({self.n_data})):",
                   "                phased_x(qs[i], angle(0.0), angle(0.0))"]
        o.append(f"        anc = array(qubit() for _ in range(py({self.n_anc})))")
        o += [f"        h(anc[{aidx[a]}])" for a in L["xanc"]]
        for ctrl, tgt in L["cx"]:
            cs = f"anc[{aidx[ctrl]}]" if ctrl in aidx else f"qs[{didx[ctrl]}]"
            ts = f"anc[{aidx[tgt]}]" if tgt in aidx else f"qs[{didx[tgt]}]"
            o.append(f"        cx({cs}, {ts})")
        o += [f"        h(anc[{aidx[a]}])" for a in L["xanc"]]
        # measure_array consumes the whole array: Guppy's linear types forbid consuming a
        # qubit out of an array by subscript.
        o.append('        output("s", collect_measurements(measure_array(anc)))')
        if self.basis == "X":
            o += [f"    for i in range(py({self.n_data})):", "        h(qs[i])"]
        o.append('    output("data", collect_measurements(measure_array(qs)))')
        return "\n".join(o)

    def _runner(self, tau: int, rounds: int):
        key = (tau, rounds)
        if key not in self._built:
            from selene_sim import build
            name = f"mem_{self.code}_d{self.d}_t{tau}_r{rounds}" + ("_x" if self.basis == "X" else "")
            path = self._gen / f"{name}.py"
            # Guppy parses the decorated function with inspect.getsourcelines, so the
            # program must exist as real source on disk; exec()ing a string fails.
            path.write_text(PREAMBLE + self.source(tau, rounds) + "\n")
            spec = importlib.util.spec_from_file_location(name, path)
            mod = importlib.util.module_from_spec(spec)
            sys.modules[name] = mod
            spec.loader.exec_module(mod)
            self._built[key] = build(mod.prog.compile())
        return self._built[key]

    def run(self, tau: int, rounds: int, shots: int, offset: int = 0, seed: int = 1, **ch):
        """(syndromes [shots, rounds, n_anc], data [shots, n_data])."""
        from selene_sim import DepolarizingErrorModel, Stim
        syn, dat = [], []
        # Each yielded shot is a lazy generator whose `with TCPStream` block spans the whole
        # loop: consume as it arrives, or the child dies with ConnectionRefused.
        for shot in self._runner(tau, rounds).run_shots(
                simulator=Stim(), n_qubits=self.n_data + self.n_anc, n_shots=shots,
                error_model=DepolarizingErrorModel(random_seed=seed, **ch),
                random_seed=seed, shot_offset=offset):
            pairs = list(shot)
            syn.append([[int(b) for b in v] for k, v in pairs if k == "s"])
            dat.append([int(b) for b in [v for k, v in pairs if k == "data"][0]])
        return np.array(syn), np.array(dat)

    # ---- decoding --------------------------------------------------------------
    def observable(self, dat):
        return dat[:, self.obs].sum(axis=1) % 2

    def final_syndrome(self, dat):
        return np.array([dat[:, s].sum(axis=1) % 2 for s in self.sup]).T

    def detectors(self, syn, dat):
        """[shots, (rounds+1) * n_checks] detector array. D[j,c] = s[j,c] ^ s[j-1,c].

        Over the checks of the memory basis (`cols`): Z-checks for basis Z, X-checks for X.
        """
        s = syn[:, :, self.cols]
        s = np.concatenate([s, self.final_syndrome(dat)[:, None, :]], axis=1)
        prev = np.concatenate([np.zeros_like(s[:, :1, :]), s[:, :-1, :]], axis=1)
        return (s ^ prev).reshape(s.shape[0], -1)

    def priors(self, tau: int, p_1q: float, p_meas: float, p_2q: float = 0.0):
        """Edge weights from the measured per-qubit op counts.

        The data prior is accumulate(ops + tau): the extraction term dominates at small
        tau and the idle term at large tau, so an error in `ops` is tau-dependent.
        With p_2q on, the Rzz term composes mod 2 with it, n_2 events at 8*p_2q/15 each.
        """
        if self.basis != "Z":
            raise ValueError("the hand-built priors are calibrated for the Z sector only")
        f = 2 * p_1q / 3
        p_d = accumulate(self.ops_per_qubit + tau, f)
        p_anc = accumulate(self.ops_anc, f)
        if p_2q:
            if self.n2_per_qubit is None:
                raise ValueError("p_2q > 0 needs n_2 counts: pass them in `constants`")
            f2 = 8 * p_2q / 15
            q_d, q_anc = accumulate(self.n2_per_qubit, f2), accumulate(self.n2_anc, f2)
            p_d = p_d + q_d - 2 * p_d * q_d
            p_anc = p_anc + q_anc - 2 * p_anc * q_anc
        return np.maximum(p_d, 1e-5), max(1 - (1 - p_meas) * (1 - p_anc), 1e-5)

    def matching(self, rounds: int, p_d, p_m):
        """`p_d` is per data qubit, or a scalar broadcast over all of them."""
        n_z = len(self.zsup)
        p_d = np.broadcast_to(np.asarray(p_d, float), (self.n_data,))
        m = pymatching.Matching()
        w = lambda p: float(np.log((1 - p) / max(p, 1e-9)))
        for j in range(rounds + 1):
            for qi in range(self.n_data):
                ch = [c for c, s in enumerate(self.zsup) if qi in s]
                fset = {0} if qi in self.obs else set()
                if len(ch) == 2:
                    m.add_edge(j * n_z + ch[0], j * n_z + ch[1], fault_ids=fset,
                               weight=w(p_d[qi]), merge_strategy="independent")
                elif len(ch) == 1:
                    m.add_boundary_edge(j * n_z + ch[0], fault_ids=fset,
                                        weight=w(p_d[qi]), merge_strategy="independent")
        for j in range(rounds):
            for c in range(n_z):
                m.add_edge(j * n_z + c, (j + 1) * n_z + c, fault_ids=set(),
                           weight=w(p_m), merge_strategy="independent")
        return m

    def failures(self, tau, rounds, shots, offset=0, p_1q=0.0, p_2q=0.0, p_meas=0.0):
        """Number of shots the decoder got wrong."""
        syn, dat = self.run(tau, rounds, shots, offset=offset,
                            p_1q=p_1q, p_2q=p_2q, p_meas=p_meas)
        p_d, p_m = self.priors(tau, p_1q, p_meas, p_2q)
        det = self.detectors(syn, dat).astype(np.uint8)
        pred = self.matching(rounds, p_d, p_m).decode_batch(det)[:, 0]
        return int(np.sum(pred != self.observable(dat)))
