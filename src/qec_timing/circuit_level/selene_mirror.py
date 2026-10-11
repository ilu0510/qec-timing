"""An exact Stim copy of what Selene simulates, so Selene's shots can be decoded properly.

Why: ``surface_memory.Geometry.matching`` is a hand-built phenomenological graph -- data edges
within a round, measurement edges between rounds. A fault *during* extraction (after some of a
qubit's CNOTs, before others) or one propagated from an ancilla produces detection patterns that
graph can only explain as several independent faults. This module rebuilds Selene's circuit gate
for gate in Stim, so Stim can derive the full detector error model, and Selene's own detection
events can be decoded against it.

How the copy is made (each step checked against the installed selene-sim 0.3.0):

* The native gate sequence is read from one noiseless Selene shot with ``CircuitExtractor``:
  ``get_error_model_output()`` under ``IdealErrorModel`` is exactly what the error model passes
  to the simulator -- ``Reset``, ``Rxy``, ``Rz``, ``Rzz``, ``MeasureRequest``.
* Each native gate is converted through its unitary (``stim.Tableau.from_unitary_matrix``), so no
  angle is hand-translated. Selene's 5-gate CNOT decomposition converts to exactly ``CX``. Stim
  does not reject a non-Clifford unitary (it silently returns a wrong Clifford), so every
  conversion is checked by rebuilding the unitary from the tableau; a mismatch raises.
* Noise placement was read off Selene's own instruction log at ``p_1q = 1``: one uniformly random
  Pauli is inserted **before** every native single-qubit gate (the idle ``phased_x(0, 0)`` arrives
  as ``Rxy(0, 0)`` and is charged like any other). ``p_meas`` flips the reported bit only, so it
  becomes ``M(p_meas)``. Resets are noiseless at ``p_init = 0``.
* ``p_2q`` is placed as ``DEPOLARIZE2`` before each ``Rzz`` by analogy; that placement has NOT been
  read off the log, so ``to_stim`` refuses ``p_2q > 0`` unless ``allow_unverified_p2q`` is set.

Detectors follow ``surface_memory.Geometry.detectors`` exactly: ``D[r, c] = s[r, c] xor s[r-1, c]``
over the ancillas of the memory basis (Z-checks, or X-checks for an X-basis ``Geometry``) in
layout order, round-major, with a final layer from the data readout, and the observable is the
parity of the logical support.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import stim

from . import runcache

_X = np.array([[0, 1], [1, 0]], dtype=complex)
_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
_I = np.eye(2, dtype=complex)


def _rxy(theta, phi):
    return np.cos(theta / 2) * _I - 1j * np.sin(theta / 2) * (np.cos(phi) * _X + np.sin(phi) * _Y)


def _rz(theta):
    return np.diag([np.exp(-1j * theta / 2), np.exp(1j * theta / 2)])


def _rzz(theta):
    return np.diag(np.exp(-1j * theta / 2 * np.array([1, -1, -1, 1])))


_GATE_CACHE: dict[tuple, stim.Circuit] = {}


def _clifford(op: dict) -> tuple[stim.Circuit, list[int]]:
    """Stim circuit on local qubits 0..k-1 for one native gate, and the global qubits."""
    if op["op"] == "Rxy":
        key, qubits = ("Rxy", round(op["theta"], 9), round(op["phi"], 9)), [op["qubit"]]
        make = lambda: _rxy(op["theta"], op["phi"])
    elif op["op"] == "Rz":
        key, qubits = ("Rz", round(op["theta"], 9)), [op["qubit"]]
        make = lambda: _rz(op["theta"])
    elif op["op"] == "Rzz":
        key, qubits = ("Rzz", round(op["theta"], 9)), [op["qubit0"], op["qubit1"]]
        make = lambda: _rzz(op["theta"])
    else:
        raise ValueError(f"not a gate: {op}")
    if key not in _GATE_CACHE:
        U = make()
        tab = stim.Tableau.from_unitary_matrix(U, endian="little")
        # from_unitary_matrix does NOT reject a non-Clifford unitary: it silently returns some
        # Clifford (Rz(0.3) comes back as the identity). Check the round trip up to global phase.
        V = tab.to_unitary_matrix(endian="little")
        if abs(abs(np.trace(V.conj().T @ U)) - U.shape[0]) > 1e-6:
            raise ValueError(f"not a Clifford gate, cannot be mirrored: {op}")
        _GATE_CACHE[key] = tab.to_circuit(method="elimination")
    return _GATE_CACHE[key], qubits


def _append_remapped(out: stim.Circuit, local: stim.Circuit, qubits: list[int]) -> None:
    for inst in local:
        out.append(inst.name, [qubits[t.value] for t in inst.targets_copy()],
                   inst.gate_args_copy())


def native_ops(instance, n_qubits: int) -> list[dict]:
    """The native operations of one noiseless shot of a built Selene program."""
    from selene_sim import IdealErrorModel, Stim
    from selene_sim.event_hooks.instruction_log import CircuitExtractor

    hook = CircuitExtractor()
    for shot in instance.run_shots(simulator=Stim(), n_qubits=n_qubits, n_shots=1,
                                   error_model=IdealErrorModel(), random_seed=1,
                                   event_hook=hook):
        list(shot)
    return hook.shots[0].get_error_model_output()


@dataclass
class MemoryLayout:
    """What the detector definition needs; matches ``surface_memory.Geometry``."""
    n_data: int
    n_anc: int
    rounds: int
    cols: list[int]           # positions of the decoded checks in the per-round ancilla array
    sup: list[list[int]]      # data indices of each decoded check, in cols order
    obs: list[int]            # data indices of the logical observable

    @classmethod
    def from_geometry(cls, G, rounds: int) -> "MemoryLayout":
        return cls(G.n_data, G.n_anc, rounds, list(G.cols), [list(s) for s in G.sup],
                   list(G.obs))


def to_stim(ops: list[dict], layout: MemoryLayout, p_1q: float = 0.0, p_2q: float = 0.0,
            p_meas: float = 0.0, allow_unverified_p2q: bool = False, *,
            idle: tuple[float, float, float] | None = None, p_init: float = 0.0,
            rz_noise: bool = True) -> stim.Circuit:
    """Stim circuit of the same native gates, noise placed as Selene places it, with detectors.

    The defaults reproduce Selene's ``DepolarizingErrorModel``. Three extensions go beyond
    what Selene can express, for modelling other hardware on the same native circuit:

    idle      (pX, pY, pZ): one ``PAULI_CHANNEL_1`` on every data qubit at the start of each
              round, before its ancillas are reset -- the placement of Stim's generated
              ``before_round_data_depolarization``, so one channel covers the whole interval.
    p_init    ``X_ERROR`` after every reset (Stim's ``after_reset_flip_probability``).
    rz_noise  False leaves ``Rz`` noiseless (``p_1q`` then charges ``Rxy`` only).

    ``ops`` may also hold ``{"op": "Idle", "qubits": [...], "channel": (pX, pY, pZ)}`` pseudo-ops,
    written as ``PAULI_CHANNEL_1`` where they stand (``native_memory.TimedNativeMemory``).
    """
    if p_2q and not allow_unverified_p2q:
        raise ValueError("p_2q placement has not been verified against Selene's log")
    L = layout
    c = stim.Circuit()
    resets: list[int] = []
    data_ids: list[int] | None = None
    anc_ids: list[int] = []                   # current round's ancilla qubit ids, array order
    rec: dict[tuple, int] = {}                # (round, anc index) or ('data', index) -> record
    n_meas, rnd = 0, -1
    for op in ops:
        kind = op["op"]
        if kind == "Reset":
            if data_ids is not None and (len(anc_ids) == L.n_anc or rnd < 0):
                anc_ids, rnd = [], rnd + 1                # first ancilla reset of a new round
                if idle is not None:
                    c.append("PAULI_CHANNEL_1", data_ids, list(idle))
            c.append("R", [op["qubit"]])
            if p_init:
                c.append("X_ERROR", [op["qubit"]], p_init)
            resets.append(op["qubit"])
            if data_ids is None:
                if len(resets) == L.n_data:
                    data_ids = list(resets)
            else:
                anc_ids.append(op["qubit"])
        elif kind in ("Rxy", "Rz"):
            if p_1q and (rz_noise or kind == "Rxy"):
                c.append("DEPOLARIZE1", [op["qubit"]], p_1q)
            _append_remapped(c, *_clifford(op))
        elif kind == "Rzz":
            if p_2q:
                c.append("DEPOLARIZE2", [op["qubit0"], op["qubit1"]], p_2q)
            _append_remapped(c, *_clifford(op))
        elif kind == "MeasureRequest":
            q = op["qubit"]
            c.append("M", [q], p_meas) if p_meas else c.append("M", [q])
            if q in anc_ids and (rnd, anc_ids.index(q)) not in rec:
                rec[(rnd, anc_ids.index(q))] = n_meas
            elif data_ids is not None and q in data_ids:
                rec[("data", data_ids.index(q))] = n_meas
            else:
                raise ValueError(f"measurement of unexpected qubit {q}")
            n_meas += 1
        elif kind == "Idle":                     # pseudo-op from native_memory's scheduler
            c.append("PAULI_CHANNEL_1", op["qubits"], list(op["channel"]))
        elif kind in ("QAlloc", "QFree", "BatchStart"):
            continue
        else:
            raise ValueError(f"operation not mirrored: {op}")
    assert data_ids is not None and rnd + 1 == L.rounds, (rnd + 1, L.rounds)
    assert len(resets) == L.n_data + L.rounds * L.n_anc, len(resets)
    assert len(rec) == L.rounds * L.n_anc + L.n_data, len(rec)

    def r(key):                               # absolute record index -> stim relative target
        return stim.target_rec(rec[key] - n_meas)

    for j in range(L.rounds + 1):
        for ci, col in enumerate(L.cols):
            t = []
            if j < L.rounds:
                t.append(r((j, col)))
            else:
                t += [r(("data", q)) for q in L.sup[ci]]
            if j > 0:
                t.append(r((j - 1, col)))
            c.append("DETECTOR", t)
    c.append("OBSERVABLE_INCLUDE", [r(("data", q)) for q in L.obs], 0)
    return c


def matcher(circuit: stim.Circuit):
    """PyMatching decoder built from the circuit's full detector error model."""
    import pymatching
    dem = circuit.detector_error_model(decompose_errors=True)
    return pymatching.Matching.from_detector_error_model(dem)


_OPS_CACHE: dict[tuple, list[dict]] = {}


def mirror_circuit(G, tau: int, rounds: int, p_1q: float, p_meas: float) -> stim.Circuit:
    """The Stim copy of ``G``'s memory circuit at (tau, rounds), with Selene's noise."""
    key = (G.code, G.d, tau, rounds)
    if key not in _OPS_CACHE:
        _OPS_CACHE[key] = native_ops(G._runner(tau, rounds), G.n_data + G.n_anc)
    return to_stim(_OPS_CACHE[key], MemoryLayout.from_geometry(G, rounds), p_1q=p_1q, p_meas=p_meas)


def decode_hand(G, tau: int, rounds: int, det: np.ndarray, obs: np.ndarray,
                p_1q: float, p_meas: float) -> int:
    """Failures of the hand-built graph (``Geometry.matching``, ``Geometry.priors``)."""
    p_d, p_m = G.priors(tau, p_1q, p_meas)
    pred = G.matching(rounds, p_d, p_m).decode_batch(det)[:, 0]
    return int(np.sum(pred != obs))


def decode_mirror(G, tau: int, rounds: int, det: np.ndarray, obs: np.ndarray,
                  p_1q: float, p_meas: float) -> int:
    """Failures of PyMatching on the full detector error model of the Stim copy."""
    pred = matcher(mirror_circuit(G, tau, rounds, p_1q, p_meas)).decode_batch(det)[:, 0]
    return int(np.sum(pred != obs))


def compare_decoders(G, tau: int, rounds: int, shots: int, offset: int = 0, seed: int = 1,
                     p_1q: float = 0.0, p_meas: float = 0.0) -> dict:
    """One block of Selene shots (stored via ``memory_shots``), decoded by both decoders."""
    from . import memory_shots
    det, obs = memory_shots.get(G, tau, rounds, shots, offset=offset, seed=seed,
                                p_1q=p_1q, p_meas=p_meas)
    return {"shots": int(shots),
            "hand_failures": decode_hand(G, tau, rounds, det, obs, p_1q, p_meas),
            "mirror_failures": decode_mirror(G, tau, rounds, det, obs, p_1q, p_meas),
            "raw_observable_flips": int(obs.sum())}


DECODED = "selene_memory_decoded"


def decoded(G, tau: int, rounds: int, shots: int, offset: int = 0, seed: int = 1,
            p_1q: float = 0.0, p_meas: float = 0.0) -> dict:
    """Failure counts of both decoders on one stored block of shots, cached.

    The shots come from ``memory_shots.get`` (stored once, never overwritten); the counts are
    cached under ``selene_memory_decoded`` with the same parameters plus the decoder versions,
    so a notebook can read them back without decoding again.
    """
    from . import memory_shots
    prm = {**memory_shots.params(G, tau, rounds, shots, offset, seed, p_1q, p_meas),
           "decoders": "hand:Geometry.priors|mirror:full-dem-pymatching"}
    hit = runcache.load(DECODED, prm)
    if hit is not None:
        return hit
    det, obs = memory_shots.get(G, tau, rounds, shots, offset=offset, seed=seed,
                                p_1q=p_1q, p_meas=p_meas)
    hit = {"shots": int(shots),
           "hand_failures": decode_hand(G, tau, rounds, det, obs, p_1q, p_meas),
           "mirror_failures": decode_mirror(G, tau, rounds, det, obs, p_1q, p_meas),
           "raw_observable_flips": int(obs.sum())}
    runcache.save(DECODED, prm, hit)
    return hit
