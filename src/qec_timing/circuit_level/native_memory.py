"""Surface-code memory on Selene's NATIVE gates, in Stim, with emulator-style idle noise.

``stim_memory.StimMemory`` runs Stim's textbook circuit: a CX is one gate, charged ``p2`` once.
Selene compiles each CX into ``Rzz`` plus single-qubit ``Rxy``/``Rz`` gates (06a), so on hardware
a round carries more, and differently placed, gate faults. Selene itself cannot express the idle
noise trapped ions have (Z only, linear and quadratic in time: no duration term, notebook 03).
This class takes the native gate sequence from Selene and the idle channel from ``stim_memory``:

    gates     Selene's compiled program for ``surface_memory.Geometry`` (tau = 0: no idle ops),
              read off one noiseless shot and rebuilt gate for gate by ``selene_mirror``
    gate noise  DEPOLARIZE1(p1) on every Rxy (and Rz if ``rz_noise``), DEPOLARIZE2(p2) on every
              Rzz, X_ERROR(p_init) after every reset, the record flipped with p_meas
    idle      ``stim_memory.idle_channel(tau, ...)`` on every data qubit once per round, where
              Stim's generated circuit puts it, so the two classes differ ONLY in the gates

``p2`` here is a modelling choice (H2's two-qubit error on each native Rzz), not a copy of
Selene's ``p_2q``, whose placement ``selene_mirror`` has not verified; a depolarizing channel
commutes with its own Clifford gate, so before/after placement is the same channel either way.

    M = NativeMemory("rotated", 5, basis="X", p1=9e-5, p2=5e-3, p_meas=4.7e-3, p_init=1.9e-4,
                     dephase_linear=1.5e-3)
    pL, sigma, fails = M.logical_error(tau=4, rounds=48, shots=100_000)
"""

from __future__ import annotations

import numpy as np
import pymatching
import stim

from . import runcache
from . import selene_mirror as sm
from .stim_memory import idle_channel
from .surface_memory import Geometry

OPS_EXPERIMENT = "native_memory_ops"
_KINDS = ("Reset", "Rxy", "Rz", "Rzz", "MeasureRequest", "QAlloc", "QFree", "BatchStart")
_FIELDS = {"qubit", "qubit0", "qubit1", "theta", "phi"}


def _pack(ops: list[dict]) -> dict:
    """Native ops as flat arrays (a d = 7, 192-round program is ~2e5 ops: too big for JSON)."""
    for o in ops:
        if o["op"] not in _KINDS or set(o) - {"op"} - _FIELDS:
            raise ValueError(f"cannot pack native op: {o}")
    nan = float("nan")
    return dict(kind=np.array([_KINDS.index(o["op"]) for o in ops], np.int8),
                q0=np.array([o.get("qubit", o.get("qubit0", -1)) for o in ops], np.int32),
                q1=np.array([o.get("qubit1", -1) for o in ops], np.int32),
                theta=np.array([o.get("theta", nan) for o in ops], float),
                phi=np.array([o.get("phi", nan) for o in ops], float))


def _unpack(a: dict) -> list[dict]:
    ops = []
    for k, q0, q1, th, ph in zip(a["kind"], a["q0"], a["q1"], a["theta"], a["phi"]):
        kind = _KINDS[k]
        if kind == "Rzz":
            o = {"op": kind, "qubit0": int(q0), "qubit1": int(q1), "theta": float(th)}
        else:
            o = {"op": kind}
            if q0 >= 0:
                o["qubit"] = int(q0)
            if not np.isnan(th):
                o["theta"] = float(th)
            if not np.isnan(ph):
                o["phi"] = float(ph)
        ops.append(o)
    return ops


class NativeMemory:
    """One code, one basis, one set of rates; circuits and decoders are built per (tau, R)."""

    _geometries: dict[tuple, Geometry] = {}

    def __init__(self, code: str, d: int, basis: str = "X", p1: float = 0.0,
                 p2: float = 0.0, p_meas: float = 0.0, p_init: float = 0.0,
                 idle_depol: float = 0.0, dephase_linear: float = 0.0,
                 dephase_quad: float = 0.0, rz_noise: bool = False):
        key = (code, d, basis)
        if key not in self._geometries:
            self._geometries[key] = Geometry(code, d, basis=basis)
        self.G = self._geometries[key]
        self.code, self.d, self.basis = code, d, basis
        self.p1, self.p2, self.p_meas, self.p_init = p1, p2, p_meas, p_init
        self.rz_noise = rz_noise
        self.idle = dict(idle_depol=idle_depol, dephase_linear=dephase_linear,
                         dephase_quad=dephase_quad)

    def ops(self, rounds: int) -> list[dict]:
        """Selene's native operations for `rounds` rounds, no idle ops; cached on disk."""
        G = self.G
        prm = dict(code=self.code, distance=self.d, basis=self.basis, rounds=int(rounds),
                   circuit="surface_memory.Geometry, tau = 0")
        arr = runcache.load_arrays(OPS_EXPERIMENT, prm)
        if arr is None:
            ops = sm.native_ops(G._runner(0, rounds), G.n_data + G.n_anc)
            runcache.save_arrays(OPS_EXPERIMENT, prm, **_pack(ops))
            runcache.save(OPS_EXPERIMENT, prm, {"n_ops": len(ops)})
            return ops
        return _unpack(arr)

    def circuit(self, tau: float, rounds: int) -> stim.Circuit:
        return sm.to_stim(self.ops(rounds), sm.MemoryLayout.from_geometry(self.G, rounds),
                          p_1q=self.p1, p_2q=self.p2, p_meas=self.p_meas,
                          allow_unverified_p2q=True, idle=idle_channel(tau, **self.idle),
                          p_init=self.p_init, rz_noise=self.rz_noise)

    def logical_error(self, tau: float, rounds: int, shots: int, seed: int = 1):
        """(pL, sigma, failures), decoded with PyMatching on the circuit's detector error model."""
        c = self.circuit(tau, rounds)
        dem = c.detector_error_model(decompose_errors=True)
        matcher = pymatching.Matching.from_detector_error_model(dem)
        det, obs = c.compile_detector_sampler(seed=seed).sample(
            shots, separate_observables=True)
        fails = int(np.sum(matcher.decode_batch(det)[:, 0] != obs[:, 0]))
        e = fails / shots
        return e, float(np.sqrt(max(e * (1 - e), 1e-12) / shots)), fails


# ---- time-resolved idle noise on the native circuit -------------------------------------------

def segments(ops: list[dict], n_data: int, n_anc: int) -> tuple[list, list[list], list]:
    """Split a native op list into (prep, rounds, final).

    The data qubits are the first `n_data` reset. A round starts at the first reset of a
    non-data qubit and ends when `n_anc` non-data qubits have been measured; what follows the
    last round (the data's final basis change and readout) is `final`.
    """
    ops = [o for o in ops if o["op"] not in ("QAlloc", "QFree", "BatchStart")]
    data, i = set(), 0
    while len(data) < n_data:
        if ops[i]["op"] == "Reset":
            data.add(ops[i]["qubit"])
        i += 1
    prep, rounds, final, cur, n_meas = ops[:i], [], [], None, 0
    for o in ops[i:]:
        q = o.get("qubit")
        if cur is None:
            if o["op"] == "Reset" and q not in data:
                cur, n_meas, final = [], 0, []      # ops between rounds were not final after all
            else:
                (final if rounds else prep).append(o)
                continue
        cur.append(o)
        if o["op"] == "MeasureRequest" and q not in data:
            n_meas += 1
            if n_meas == n_anc:
                rounds.append(cur)
                cur = None
    assert cur is None, "op list ends inside a round"
    return prep, rounds, final


def _qubits(o: dict) -> list[int]:
    return [o["qubit0"], o["qubit1"]] if o["op"] == "Rzz" else [o["qubit"]]


def layers(round_ops: list[dict]) -> list[list[dict]]:
    """ASAP schedule: each op in the first layer after every earlier op on its qubits.

    Order on each qubit is preserved, so the circuit is unchanged; only commuting operations on
    disjoint qubits are moved together.
    """
    free, out = {}, []
    for o in round_ops:
        qs = _qubits(o)
        k = max(free.get(q, 0) for q in qs)
        for q in qs:
            free[q] = k + 1
        while len(out) <= k:
            out.append([])
        out[k].append(o)
    return out


class TimedNativeMemory(NativeMemory):
    """``NativeMemory`` with time running through the round, and the same idle noise on every
    qubit (data and ancilla alike).

    Each round's native ops are scheduled ASAP into layers; the round lasts ``t_R = 1`` and each
    of its L layers ``1 / L``. After every layer, every qubit gets ``idle_channel(1 / L)``. At the
    start of each round every qubit gets ``idle_channel(tau - 1)`` for the wait (``tau`` is the full
    round period, >= 1). The data preparation before the first round and the readout after the
    last are not given idle time. Gate noise is ``NativeMemory``'s.
    """

    def schedule(self, rounds: int) -> tuple[list, list[list[list[dict]]], list]:
        prep, rnds, final = segments(self.ops(rounds), self.G.n_data, self.G.n_anc)
        assert len(rnds) == rounds, (len(rnds), rounds)
        return prep, [layers(r) for r in rnds], final

    def circuit(self, tau: float, rounds: int) -> stim.Circuit:
        if tau < 1:
            raise ValueError("tau is the full round period and cannot be shorter than a round")
        prep, rnds, final = self.schedule(rounds)
        allq = sorted({q for o in self.ops(rounds) if "qubit" in o or "qubit0" in o
                       for q in _qubits(o)})
        wait = idle_channel(tau - 1.0, **self.idle)
        has = lambda ch: any(v > 0 for v in ch)          # noqa: E731
        ops = list(prep)
        for lay in rnds:
            if has(wait):
                ops.append({"op": "Idle", "qubits": allq, "channel": wait})
            step = idle_channel(1.0 / len(lay), **self.idle)
            for layer_ops in lay:
                ops += layer_ops
                if has(step):
                    ops.append({"op": "Idle", "qubits": allq, "channel": step})
        ops += final
        return sm.to_stim(ops, sm.MemoryLayout.from_geometry(self.G, rounds),
                          p_1q=self.p1, p_2q=self.p2, p_meas=self.p_meas,
                          allow_unverified_p2q=True, p_init=self.p_init, rz_noise=self.rz_noise)
