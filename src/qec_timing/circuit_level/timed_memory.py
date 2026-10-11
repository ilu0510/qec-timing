"""Circuit-level memory with time-resolved dephasing on EVERY qubit, and the measurement of the
phenomenological parameters from its detection events.

``stim_memory.StimMemory`` dephases the data once per interval and never the ancillas. Here time
runs through the round as well: Stim's generated memory circuit has 7 ``TICK``-separated layers
per round (H, four CX, H, MR), and each layer is taken to last ``t_R / 7``. **Every qubit has the
same noise model, ancillas included**: after every ``TICK`` every qubit gets the idle channel for
that layer's duration, and over the wait between rounds, ``tau - 1`` (``tau`` is the full round
period, >= 1), every qubit gets it again (on a freshly reset ancilla a Z does nothing, which is a
consequence, not a rule). Gate noise is ``StimMemory``'s and does not care which qubit is which:
``p1`` after each H, ``p2`` (DEPOLARIZE2) after every CX, ``p_meas`` before each measurement,
``p_init`` after each reset. ``ancilla_idle=False`` exempts the ancillas from idle noise; it exists
only as a diagnostic, to show what ancilla dephasing contributes. Nothing here sets a readout multiplier: whatever ``b_read`` the circuit has comes out
of it, and is measured by ``edge_probabilities``.

    M = TimedMemory("rotated", 5, p1=9e-5, p2=5e-3, p_meas=4.7e-3, p_init=1.9e-4,
                    dephase_linear=1.5e-3)
    pL, sigma, fails = M.logical_error(tau=4, rounds=48, shots=100_000)

``phenomenological_circuit`` is the paper's model on the same code (data Z flips once per round,
syndrome-record flips, ideal gates), for predictions from measured parameters.
"""

from __future__ import annotations

from collections import defaultdict

import numpy as np
import pymatching
import stim

from .stim_memory import SPEC, _GATE, _IDLE, _MEAS, _RESET, _is, idle_channel

LAYERS = 7          # TICKs per round in Stim's generated rotated memory circuit


def _decode(c: stim.Circuit, shots: int, seed: int):
    dem = c.detector_error_model(decompose_errors=True)
    matcher = pymatching.Matching.from_detector_error_model(dem)
    det, obs = c.compile_detector_sampler(seed=seed).sample(shots, separate_observables=True)
    fails = int(np.sum(matcher.decode_batch(det)[:, 0] != obs[:, 0]))
    e = fails / shots
    return e, float(np.sqrt(max(e * (1 - e), 1e-12) / shots)), fails


def _data_qubits(c: stim.Circuit) -> list[int]:
    """Data qubits: those measured by the final (non-resetting) data measurement."""
    last = [i for i in c.flattened() if i.name in ("M", "MX", "MZ")][-1]
    return [t.value for t in last.targets_copy()]


class TimedMemory:
    """One code, one basis, one set of rates; circuits and decoders are built per (tau, R)."""

    def __init__(self, code: str, d: int, basis: str = "X", p1: float = 0.0,
                 p2: float = 0.0, p_meas: float = 0.0, p_init: float = 0.0,
                 dephase_linear: float = 0.0, dephase_quad: float = 0.0,
                 ancilla_idle: bool = True):
        if (code, basis) not in SPEC:
            raise ValueError(f"unknown code/basis {code}/{basis}")
        self.code, self.d, self.basis = code, d, basis
        self.p1, self.p2, self.p_meas, self.p_init = p1, p2, p_meas, p_init
        self.gamma, self.kappa = dephase_linear, dephase_quad
        self.ancilla_idle = ancilla_idle

    def idle(self, duration: float) -> tuple[float, float, float]:
        return idle_channel(duration, 0.0, self.gamma, self.kappa)

    def circuit(self, tau: float, rounds: int) -> stim.Circuit:
        if tau < 1:
            raise ValueError("tau is the full round period and cannot be shorter than a round")
        tagged = stim.Circuit.generated(
            SPEC[(self.code, self.basis)], distance=self.d, rounds=rounds,
            before_round_data_depolarization=_IDLE, after_clifford_depolarization=_GATE,
            before_measure_flip_probability=_MEAS, after_reset_flip_probability=_RESET)
        data = _data_qubits(tagged)
        used = sorted({t.value for i in tagged.flattened() for t in i.targets_copy()
                       if t.is_qubit_target})
        idle_q = used if self.ancilla_idle else data
        layer, wait = self.idle(1.0 / LAYERS), self.idle(tau - 1.0)
        has = lambda ch: any(v > 0 for v in ch)          # noqa: E731

        def rewrite(block: stim.Circuit) -> stim.Circuit:
            out = stim.Circuit()
            for inst in block:
                if isinstance(inst, stim.CircuitRepeatBlock):
                    out.append(stim.CircuitRepeatBlock(inst.repeat_count,
                                                       rewrite(inst.body_copy())))
                    continue
                name, args, tg = inst.name, inst.gate_args_copy(), inst.targets_copy()
                if name == "TICK":
                    out.append(inst)
                    if has(layer):
                        out.append("PAULI_CHANNEL_1", idle_q, list(layer))
                elif name == "DEPOLARIZE1" and _is(args[0], _IDLE):
                    if has(wait):
                        out.append("PAULI_CHANNEL_1", idle_q, list(wait))
                elif name == "DEPOLARIZE1" and _is(args[0], _GATE):
                    out.append("DEPOLARIZE1", tg, [self.p1])
                elif name == "DEPOLARIZE2" and _is(args[0], _GATE):
                    out.append("DEPOLARIZE2", tg, [self.p2])
                elif name in ("X_ERROR", "Z_ERROR") and _is(args[0], _MEAS):
                    out.append(name, tg, [self.p_meas])
                elif name in ("X_ERROR", "Z_ERROR") and _is(args[0], _RESET):
                    out.append(name, tg, [self.p_init])
                else:
                    out.append(inst)
            return out

        c = rewrite(tagged)
        n_tick = sum(1 for i in c.flattened() if i.name == "TICK")
        assert n_tick == LAYERS * rounds, (n_tick, rounds)
        left = [v for v in (_IDLE, _GATE, _MEAS, _RESET) if f"{v}" in str(c)]
        assert not left, f"untranslated noise tags remain: {left}"
        return c

    def logical_error(self, tau: float, rounds: int, shots: int, seed: int = 1):
        """(pL, sigma, failures), decoded with PyMatching on the circuit's detector error model."""
        return _decode(self.circuit(tau, rounds), shots, seed)


def phenomenological_circuit(code: str, d: int, rounds: int, p_data: float, p_read: float,
                             basis: str = "X") -> stim.Circuit:
    """The paper's model: ideal gates; once per round a flip on every data qubit with
    probability `p_data` (Z for an X-basis memory), and every recorded bit -- syndrome and final
    data readout -- flipped with probability `p_read`."""
    tagged = stim.Circuit.generated(SPEC[(code, basis)], distance=d, rounds=rounds,
                                    before_round_data_depolarization=_IDLE,
                                    before_measure_flip_probability=_MEAS)
    flip = "Z_ERROR" if basis == "X" else "X_ERROR"

    def rewrite(block: stim.Circuit) -> stim.Circuit:
        out = stim.Circuit()
        for inst in block:
            if isinstance(inst, stim.CircuitRepeatBlock):
                out.append(stim.CircuitRepeatBlock(inst.repeat_count, rewrite(inst.body_copy())))
                continue
            name, args, tg = inst.name, inst.gate_args_copy(), inst.targets_copy()
            if name == "DEPOLARIZE1" and _is(args[0], _IDLE):
                out.append(flip, tg, [p_data])
            elif name in ("X_ERROR", "Z_ERROR") and _is(args[0], _MEAS):
                out.append(name, tg, [p_read])
            else:
                out.append(inst)
        return out

    return rewrite(tagged)


def phenomenological_error(code, d, rounds, p_data, p_read, shots, seed=1, basis="X"):
    return _decode(phenomenological_circuit(code, d, rounds, p_data, p_read, basis), shots, seed)


# ---- measuring the phenomenological parameters from detection events ---------------------------

def check_edges(c: stim.Circuit) -> dict:
    """Detector pairs of the decoded checks, classified, from a generated memory circuit.

    The decoded checks are the ones with a detector in the first round (X-checks for an X-basis
    memory). Returns lists of detector-index pairs, bulk rounds only (neither the first round nor
    the final data-readout layer):

        time      (c, r) - (c, r+1)                       a flipped syndrome record
        space     (a, r) - (b, r),  a, b share a data qubit
        diag      (a, r) - (b, r+1) and (a, r+1) - (b, r)  a fault between a's and b's CNOTs
    with ``pair[k]`` naming the check pair (a, b) behind space[k] and diag[k].
    """
    coords = c.get_detector_coordinates()
    qcoord = {q: tuple(v) for q, v in c.get_final_qubit_coordinates().items()}
    first = {tuple(v[:2]) for v in coords.values() if v[2] == 0}
    det_at = {(tuple(v[:2]), int(v[2])): k for k, v in coords.items() if tuple(v[:2]) in first}
    rounds = max(r for _, r in det_at)                     # syndrome rounds
    # support of each decoded check: CX partners of its ancilla
    sup = defaultdict(set)
    for inst in c.flattened():
        if inst.name == "CX":
            t = [g.value for g in inst.targets_copy()]
            for a, b in zip(t[::2], t[1::2]):
                for anc, other in ((a, b), (b, a)):
                    xy = qcoord[anc][:2]
                    if xy in first:
                        sup[xy].add(other)
    return _classify(det_at, sorted(first), sup, rounds)


def layout_edges(sup: list[list[int]], rounds: int) -> dict:
    """``check_edges`` for detectors laid out round-major over checks, as ``selene_mirror``
    writes them: detector ``j * len(sup) + i`` compares check ``i`` in rounds ``j`` and ``j-1``."""
    n = len(sup)
    det_at = {(i, j): j * n + i for j in range(rounds + 1) for i in range(n)}
    return _classify(det_at, list(range(n)), {i: set(s) for i, s in enumerate(sup)}, rounds)


def _classify(det_at: dict, checks: list, sup: dict, rounds: int) -> dict:
    pairs = [(a, b) for i, a in enumerate(checks) for b in checks[i + 1:] if sup[a] & sup[b]]
    assert all(len(sup[a] & sup[b]) == 1 for a, b in pairs), "checks sharing >1 data qubit"
    out = {"time": [], "space": [], "diag": [], "pair": [], "checks": checks, "pairs": pairs}
    for r in range(1, rounds - 1):       # r and r+1 both in 1 .. rounds-1 (syndrome rounds)
        for ch in checks:
            out["time"].append((det_at[(ch, r)], det_at[(ch, r + 1)]))
        for a, b in pairs:
            out["space"].append((det_at[(a, r)], det_at[(b, r)]))
            out["diag"].append([(det_at[(a, r)], det_at[(b, r + 1)]),
                                (det_at[(a, r + 1)], det_at[(b, r)])])
            out["pair"].append((a, b))
    return out


def pij(det: np.ndarray, pairs) -> np.ndarray:
    """Edge probability between detectors i and j, from their firing statistics.

    p_ij = 1/2 - 1/2 sqrt(1 - 4 (<x_i x_j> - <x_i><x_j>) / (1 - 2<x_i> - 2<x_j> + 4<x_i x_j>)),
    the standard estimator for independent graph-like errors (Google, Nature 595, 383 (2021),
    supplement). It is exact when every error mechanism flips at most two detectors.
    """
    p = np.asarray(pairs)
    x = det.mean(0)
    return _pij(x[p[:, 0]], x[p[:, 1]], (det[:, p[:, 0]] & det[:, p[:, 1]]).mean(0))


def _pij(xi, xj, xij):
    num = 4 * (xij - xi * xj)
    den = 1 - 2 * xi - 2 * xj + 4 * xij
    return 0.5 - 0.5 * np.sqrt(np.clip(1 - num / den, 0.0, None))


def xor(*ps):
    """Probability of an odd number of independent flips."""
    q = 1.0
    for p in ps:
        q *= 1 - 2 * np.asarray(p)
    return (1 - q) / 2


def edge_probabilities(c: stim.Circuit, shots: int, seed: int = 1, chunk: int = 50_000,
                       edges: dict | None = None) -> dict:
    """The phenomenological parameters of a circuit, measured from its detection events.

    p_read   mean time-like edge: a check's record flipping, with nothing on the data
    p_data   mean over bulk check pairs of space-like XOR both diagonal edges: everything that
             flips the shared data qubit once per round, whenever in the round it happens
    Also the components, so the share of mid-round (diagonal) faults is visible. Moments are
    accumulated over chunks of `chunk` shots, so memory does not grow with `shots`. `edges`
    defaults to ``check_edges(c)``; pass ``layout_edges`` for a circuit without coordinates.
    """
    E = edges if edges is not None else check_edges(c)
    groups = {"time": np.asarray(E["time"]), "space": np.asarray(E["space"]),
              "d1": np.asarray([x[0] for x in E["diag"]]),
              "d2": np.asarray([x[1] for x in E["diag"]])}
    sampler = c.compile_detector_sampler(seed=seed)
    n, x = 0, 0.0
    joint = {k: 0.0 for k in groups}
    while n < shots:
        det = sampler.sample(min(chunk, shots - n))
        x = x + det.sum(0, dtype=np.int64)
        for k, p in groups.items():
            joint[k] = joint[k] + (det[:, p[:, 0]] & det[:, p[:, 1]]).sum(0, dtype=np.int64)
        n += len(det)
    x = x / n
    e = {k: _pij(x[p[:, 0]], x[p[:, 1]], joint[k] / n) for k, p in groups.items()}
    pd_ = xor(e["space"], e["d1"], e["d2"])
    se = lambda a: float(a.std(ddof=1) / np.sqrt(len(a)))          # noqa: E731
    return dict(p_read=float(e["time"].mean()), p_read_se=se(e["time"]),
                p_data=float(pd_.mean()), p_data_se=se(pd_),
                p_space=float(e["space"].mean()), p_diag=float(xor(e["d1"], e["d2"]).mean()),
                n_time=len(e["time"]), n_pairs=len(e["space"]))
