"""Circuit-level Stim twin of the full surface-code memory (Stage 6a).

Gate-level syndrome extraction (no MPP) that follows
``qec_timing.circuit.program`` qubit for qubit, so it can later be compared with
the Selene run of the same circuit:

* Qubits: data ``0 .. d^2 - 1`` in ``layout.data_index`` (row-major) order and
  ONE ancilla, index ``d^2``, reset before every check, as the Guppy program
  allocates a fresh ``qubit()`` per check.
* One check: reset ancilla; for an X-check ``H``; four (or two) CNOTs over the
  support in the order of ``program._supports`` (data -> ancilla for a Z-check,
  ancilla -> data for an X-check); for an X-check ``H``; measure the ancilla.
  The support order comes from the same ``circuit.ordering.order_supports``
  helper and ``program._supports`` table the Guppy program uses, with the same
  ``order`` argument, so the two cannot disagree. The default ``"row_major"``
  is NW, NE, SW, SE for BOTH sectors; it costs the X memory circuit-level
  distance (Z-check hooks run along the logical Z row). ``"hook_safe"`` uses the
  standard Z / N shapes and keeps full distance in both memories. Checks are
  extracted strictly one after another, as in ``program.py``.
* Preparation: data reset (``H`` on all for an X memory); the complementary
  checks are measured once and their random signs corrected with
  ``program._preparation_corrections`` (as Stim feedback).
* Each round: idle channel on every data qubit for that round's ``dt``; then for
  each check index ``ci``: the active check ``ci``, then the complementary check
  ``ci``. Only active checks get detectors (round 0 compares with the
  deterministic +1 value, later rounds with the previous round).
* Final: (``H`` for an X memory) and a noisy readout of every data qubit, final
  detectors from the support parities, logical observable on
  ``logical_z_support`` (Z memory) or ``logical_x_support`` (X memory).

Detector order is round -> check, ``(n_rounds + 1) * n_checks`` in total,
identical to the Guppy program's ``det`` stream for one patch and to the
Stage-1 reference.

By default, two deliberate differences from ``program.py``:

* ``program.py`` conjugates its per-round X-noise by ``H`` in an X memory, to
  turn the phenomenological X channel into Z faults. Those ``H`` gates are a
  noise-injection device, not part of extraction, and are not emitted here: the
  idle channel acts directly.
* The preparation corrections are conditional gates; Stim cannot charge noise
  on a gate that only sometimes runs, so they are noiseless here.

``mirror_guppy_program=True`` (requires ``native_charging=True``; used for the
Stage 6b cross-check against Selene) instead charges what Selene charges when
it runs ``program.py`` with its RNG noise injection standing in for the idle
channel:

* X memory: the per-round ``H`` wrapper on every data qubit, each ``H`` charged
  2 native ops, emitted literally (``H``, noise, ``H``).
* The injected flip is an ``x`` gate that only runs when the RNG fires, and is
  then charged 1 native op. Per data qubit that is exactly the Pauli channel
  "flip w.p. p, then DEPOLARIZE1(p1)", i.e. ``PAULI_CHANNEL_1`` with the flip
  Pauli (X, in the frame inside the ``H`` wrapper) at ``p (1 - p1)`` and each
  other Pauli at ``p p1 / 3``. The idle channel must be one the injection can
  produce: ``DephasingH2`` in an X memory (X inside the wrapper == Z) or
  ``BitFlip`` in a Z memory. ``DepolarizingKOps`` is not injected; its ops are
  real native ops and are charged as such already.
* Each preparation correction runs with probability 1/2 (the projected sign is
  uniformly random) and is then charged 1 native op on each qubit it touches.
  Modelled by its marginal: ``DEPOLARIZE1(compose(p1 / 2, c_i))`` on data qubit
  ``i``, where ``c_i`` is the number of corrections touching it. The joint
  correlation within one correction is dropped; that is O(p1^2).

``program.py`` is imported for its geometry helpers only; building a circuit
here never compiles or runs a Guppy program.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import cached_property

import numpy as np
import stim

from ..circuit.ordering import Order, order_supports
from ..circuit.program import _preparation_corrections, _supports
from ..layout import RotatedSurfaceCode
from ..reference.decoder import PyMatchingDecoder
from ..reference.sampler import LogicalErrorEstimate, clopper_pearson
from ..schedule import RoundSchedule, constant_schedule
from .noise import (
    CX_CONTROL_1Q_OPS,
    CX_TARGET_1Q_OPS,
    H_NATIVE_1Q_OPS,
    PAULI_GATE_1Q_OPS,
    Basis,
    BitFlip,
    CircuitNoise,
    DephasingH2,
    IdleChannel,
    ReadoutNoise,
    compose_depolarize1,
)

__all__ = [
    "CircuitLevelMemory",
    "build_circuit_level_memory",
    "build_constant_interval_memory",
    "idle_probe_circuit",
    "rate_per_unit_time",
    "readout_probe_circuit",
    "sample_measurements",
]

_BATCH = 20_000


# ---------------------------------------------------------------------------
# Gate emission with noise
# ---------------------------------------------------------------------------


class _Emitter:
    """Appends gates with their noise; tracks the measurement count."""

    def __init__(self, noise: CircuitNoise, native_charging: bool,
                 p_meas_in_circuit: float) -> None:
        self.c = stim.Circuit()
        self.noise = noise
        self.native = native_charging
        self.p_meas = p_meas_in_circuit
        self.n_meas = 0
        p1 = noise.p1
        self.p_h = compose_depolarize1(p1, H_NATIVE_1Q_OPS) if native_charging else p1
        def ops(n: int) -> float:
            return compose_depolarize1(p1, n) if native_charging else 0.0

        self.p_target_before, self.p_target_after = map(ops, CX_TARGET_1Q_OPS)
        self.p_control_before, self.p_control_after = map(ops, CX_CONTROL_1Q_OPS)

    def _noise(self, name: str, targets: list[int], p: float) -> None:
        if p > 0.0:
            self.c.append(name, targets, p)

    def reset(self, qs: list[int]) -> None:
        self.c.append("R", qs)
        self._noise("X_ERROR", qs, self.noise.p_init)

    def h(self, qs: list[int]) -> None:
        self.c.append("H", qs)
        self._noise("DEPOLARIZE1", qs, self.p_h)

    def cx(self, control: int, target: int) -> None:
        self._noise("DEPOLARIZE1", [control], self.p_control_before)
        self._noise("DEPOLARIZE1", [target], self.p_target_before)
        self.c.append("CX", [control, target])
        self._noise("DEPOLARIZE2", [control, target], self.noise.p2)
        self._noise("DEPOLARIZE1", [control], self.p_control_after)
        self._noise("DEPOLARIZE1", [target], self.p_target_after)

    def measure(self, qs: list[int]) -> list[int]:
        """Measure; returns the absolute record indices."""
        if self.p_meas > 0.0:
            self.c.append("M", qs, self.p_meas)
        else:
            self.c.append("M", qs)
        first = self.n_meas
        self.n_meas += len(qs)
        return list(range(first, self.n_meas))

    def rec(self, index: int) -> stim.GateTarget:
        return stim.target_rec(index - self.n_meas)


# ---------------------------------------------------------------------------
# The memory experiment
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CircuitLevelMemory:
    """A built circuit-level memory experiment and how to sample/decode it."""

    d: int
    basis: Basis
    schedule: RoundSchedule
    noise: CircuitNoise
    native_charging: bool
    order: Order
    mirror_guppy_program: bool
    circuit: stim.Circuit
    """Decoding circuit: readout as ``M(p)`` (``p_avg`` if asymmetric)."""
    sampling_circuit: stim.Circuit
    """Same circuit; noiseless ``M`` when readout is asymmetric (flips in numpy)."""

    @property
    def layout(self) -> RotatedSurfaceCode:
        return RotatedSurfaceCode(self.d)

    @property
    def n_rounds(self) -> int:
        return self.schedule.n_rounds

    @property
    def num_detectors(self) -> int:
        return self.circuit.num_detectors

    @cached_property
    def detector_error_model(self) -> stim.DetectorErrorModel:
        return self.circuit.detector_error_model(decompose_errors=True)

    def decoder(self) -> PyMatchingDecoder:
        """The existing ``Decoder`` implementation, built from this circuit's DEM."""
        return PyMatchingDecoder(self.detector_error_model)

    def round_metadata(self) -> list[dict[str, float]]:
        """Per round: ``dt``, start time ``t``, idle and readout probabilities."""
        r = self.noise.readout
        return [
            {"dt": dt, "t": t, "p_idle": self.noise.idle.probability(dt),
             "p_read_01": r.p01, "p_read_10": r.p10}
            for dt, t in zip(self.schedule.intervals, self.schedule.start_times)
        ]

    def sampler(self, seed: int) -> "_Sampler":
        return _Sampler(self, seed)

    def sample(self, shots: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
        """``(detectors, observables)`` for ``shots`` shots."""
        return self.sampler(seed).sample(shots)

    def estimate_logical_error(self, shots: int, seed: int, *, decoder=None,
                               confidence: float = 0.999) -> LogicalErrorEstimate:
        if shots <= 0:
            raise ValueError("shots must be positive")
        decoder = decoder if decoder is not None else self.decoder()
        sampler = self.sampler(seed)
        failures, remaining = 0, shots
        while remaining:
            n = min(_BATCH, remaining)
            dets, obs = sampler.sample(n)
            failures += int(np.count_nonzero(
                np.any(decoder.decode_batch(dets) != obs, axis=1)))
            remaining -= n
        low, high = clopper_pearson(failures, shots, confidence)
        return LogicalErrorEstimate(failures, shots, confidence, low, high)


class _Sampler:
    """Stateful sampler: successive calls continue one seeded stream."""

    def __init__(self, memory: CircuitLevelMemory, seed: int) -> None:
        self._readout = memory.noise.readout
        if self._readout.is_symmetric:
            self._det = memory.circuit.compile_detector_sampler(seed=seed)
        else:
            self._det = None
            self._meas = memory.sampling_circuit.compile_sampler(seed=seed)
            self._m2d = memory.sampling_circuit.compile_m2d_converter()
            self._rng = _readout_rng(seed)

    def sample(self, shots: int) -> tuple[np.ndarray, np.ndarray]:
        if self._det is not None:
            return self._det.sample(shots, separate_observables=True)
        m = _flip_records(self._meas.sample(shots), self._readout, self._rng)
        return self._m2d.convert(measurements=m, separate_observables=True)


def _readout_rng(seed: int) -> np.random.Generator:
    # Separate stream from Stim's: same user seed, distinct spawn key.
    return np.random.default_rng(np.random.SeedSequence(seed, spawn_key=(6,)))


def _flip_records(m: np.ndarray, readout: ReadoutNoise,
                  rng: np.random.Generator) -> np.ndarray:
    """Exact asymmetric classical flip of each recorded bit."""
    m = np.asarray(m, dtype=np.bool_)
    p = np.where(m, readout.p10, readout.p01)
    return m ^ (rng.random(m.shape) < p)


def sample_measurements(circuit: stim.Circuit, readout: ReadoutNoise,
                        shots: int, seed: int) -> np.ndarray:
    """Raw records of a circuit with noiseless ``M``, then the readout channel.

    The same flip code path the memory's asymmetric sampler uses.
    """
    m = circuit.compile_sampler(seed=seed).sample(shots)
    return _flip_records(m, readout, _readout_rng(seed))


# ---------------------------------------------------------------------------
# Builders
# ---------------------------------------------------------------------------


def build_circuit_level_memory(
    d: int,
    schedule: RoundSchedule,
    noise: CircuitNoise,
    *,
    basis: Basis | None = None,
    native_charging: bool = False,
    order: Order = "row_major",
    mirror_guppy_program: bool = False,
    with_coords: bool = True,
) -> CircuitLevelMemory:
    """Build the memory for any schedule (constant or variable ``dt``).

    ``schedule.lam`` is ignored: the idle strength lives in ``noise.idle``.
    ``basis`` defaults to the idle channel's natural basis (Z for depolarizing
    and bit-flip, X for dephasing, which a Z memory cannot see). ``order`` is
    the CNOT order within each check, as in ``build_surface_code_program``.
    ``mirror_guppy_program`` adds the charged noise of program.py's injection
    machinery (module docstring); it needs ``native_charging=True``.
    """
    if basis is None:
        basis = noise.idle.default_basis
    if basis not in ("Z", "X"):
        raise ValueError("basis must be 'Z' or 'X'")
    if mirror_guppy_program:
        _check_mirror(noise, basis, native_charging)
    layout = RotatedSurfaceCode(d)
    readout = noise.readout
    decode_p = readout.p01 if readout.is_symmetric else readout.average
    args = (layout, schedule, noise, basis, native_charging, order,
            mirror_guppy_program, with_coords)
    circuit = _build(*args, decode_p)
    sampling = circuit if readout.is_symmetric else _build(*args, 0.0)
    return CircuitLevelMemory(d, basis, schedule, noise, native_charging, order,
                              mirror_guppy_program, circuit, sampling)


def _check_mirror(noise: CircuitNoise, basis: Basis, native_charging: bool) -> None:
    if not native_charging:
        raise ValueError("mirror_guppy_program requires native_charging=True")
    idle = noise.idle
    if isinstance(idle, DephasingH2) and basis != "X":
        raise ValueError("program.py injects Z faults only in an X memory")
    if isinstance(idle, BitFlip) and basis != "Z":
        raise ValueError("program.py injects X faults only in a Z memory")


def build_constant_interval_memory(d: int, dt: float, T: float, noise: CircuitNoise,
                                   **kwargs) -> CircuitLevelMemory:
    """Fixed ``dt``: ``n_rounds = round(T / dt)`` via ``qec_timing.schedule``."""
    # lam is a placeholder: the schedule's idling rate is unused at circuit level.
    return build_circuit_level_memory(d, constant_schedule(dt, T, 0.0), noise, **kwargs)


def _check_centre(support) -> tuple[float, float]:
    return (sum(r for r, _ in support) / len(support),
            sum(c for _, c in support) / len(support))


def _build(layout: RotatedSurfaceCode, schedule: RoundSchedule, noise: CircuitNoise,
           basis: Basis, native_charging: bool, order: Order, mirror: bool,
           with_coords: bool, p_meas: float) -> stim.Circuit:
    x_memory = basis == "X"
    n_data, n_checks = layout.n_data, layout.n_checks
    data = list(range(n_data))
    anc = n_data

    active = layout.x_checks if x_memory else layout.z_checks
    complementary = layout.z_checks if x_memory else layout.x_checks
    # The CNOT order: the same helper and table the Guppy program uses.
    active_type, other_type = ("X", "Z") if x_memory else ("Z", "X")
    sup_active = _supports(layout, order_supports(layout, active, active_type, order))
    sup_other = _supports(layout, order_supports(layout, complementary, other_type, order))
    corrections = _preparation_corrections(layout, complementary)
    logical = [layout.data_index[q] for q in
               (layout.logical_x_support if x_memory else layout.logical_z_support)]
    centres = [_check_centre(s) for s in active]

    e = _Emitter(noise, native_charging, p_meas)
    if with_coords:
        for q, (r, c) in enumerate(layout.data_qubits):
            e.c.append("QUBIT_COORDS", [q], [r, c])

    def extract(ci: int, is_x: bool, table: list[int]) -> int:
        e.reset([anc])
        if is_x:
            e.h([anc])
        for k in range(4):
            idx = table[ci * 4 + k]
            if idx >= 0:
                if is_x:
                    e.cx(anc, idx)
                else:
                    e.cx(idx, anc)
        if is_x:
            e.h([anc])
        return e.measure([anc])[0]

    # Preparation.
    e.reset(data)
    if x_memory:
        e.h(data)
    prep = [extract(ci, not x_memory, sup_other) for ci in range(n_checks)]
    fix = "CX" if x_memory else "CZ"
    for ci in range(n_checks):
        for i in range(n_data):
            if corrections[ci * n_data + i]:
                e.c.append(fix, [e.rec(prep[ci]), i])
    if mirror:
        # Each correction runs w.p. 1/2, then costs PAULI_GATE_1Q_OPS native ops.
        half = compose_depolarize1(noise.p1, PAULI_GATE_1Q_OPS) / 2.0
        by_p: dict[float, list[int]] = {}
        for i in range(n_data):
            c_i = sum(corrections[ci * n_data + i] for ci in range(n_checks))
            if c_i:
                by_p.setdefault(compose_depolarize1(half, c_i), []).append(i)
        for p, qs in sorted(by_p.items()):
            e._noise("DEPOLARIZE1", qs, p)
    e.c.append("TICK")

    # Rounds.
    prev: list[int] | None = None
    injected = mirror and isinstance(noise.idle, (DephasingH2, BitFlip))
    for j, dt in enumerate(schedule.intervals):
        if injected:
            if x_memory:
                e.h(data)
            _append_injection(e, data, noise.idle.probability(dt), noise.p1)
            if x_memory:
                e.h(data)
        else:
            noise.idle.append(e.c, data, dt)
        e.c.append("TICK")
        syn = []
        for ci in range(n_checks):
            syn.append(extract(ci, x_memory, sup_active))
            extract(ci, not x_memory, sup_other)
        for ci in range(n_checks):
            recs = [e.rec(syn[ci])]
            if prev is not None:
                recs.append(e.rec(prev[ci]))
            e.c.append("DETECTOR", recs, [*centres[ci], j] if with_coords else [])
        prev = syn
        e.c.append("TICK")

    # Final readout.
    if x_memory:
        e.h(data)
    final = e.measure(data)
    for ci in range(n_checks):
        recs = [e.rec(final[i]) for i in sup_active[ci * 4:ci * 4 + 4] if i >= 0]
        recs.append(e.rec(prev[ci]))
        e.c.append("DETECTOR", recs,
                   [*centres[ci], schedule.n_rounds] if with_coords else [])
    e.c.append("OBSERVABLE_INCLUDE", [e.rec(final[i]) for i in logical], 0)
    return e.c


def _append_injection(e: _Emitter, data: list[int], p: float, p1: float) -> None:
    """program.py's injected flip (X in its frame) plus the charge on the x gate."""
    if p <= 0.0:
        return
    q = compose_depolarize1(p1, PAULI_GATE_1Q_OPS)
    px, other = p * (1.0 - q), p * q / 3.0
    if other > 0.0:
        e.c.append("PAULI_CHANNEL_1", data, [px, other, other])
    else:
        e.c.append("X_ERROR", data, px)


# ---------------------------------------------------------------------------
# Probes and figure of merit
# ---------------------------------------------------------------------------


def idle_probe_circuit(idle: IdleChannel, dt: float, basis: Basis,
                       n_qubits: int = 100) -> stim.Circuit:
    """Prepare ``n_qubits`` in the basis's +1 state, idle once, measure.

    Noiseless apart from the idle channel, and appended by the same
    ``idle.append`` the memory uses. A measured 1 is exactly one flip in that
    basis.
    """
    qs = list(range(n_qubits))
    c = stim.Circuit()
    c.append("R", qs)
    if basis == "X":
        c.append("H", qs)
    idle.append(c, qs, dt)
    if basis == "X":
        c.append("H", qs)
    c.append("M", qs)
    return c


def readout_probe_circuit(n_zeros: int, n_ones: int) -> stim.Circuit:
    """Noiseless ``|0>^n_zeros |1>^n_ones`` then noiseless ``M`` (zeros first)."""
    c = stim.Circuit()
    qs = list(range(n_zeros + n_ones))
    c.append("R", qs)
    if n_ones:
        c.append("X", qs[n_zeros:])
    c.append("M", qs)
    return c


def rate_per_unit_time(p_L: float, total_time: float) -> float:
    """Logical error rate per unit time from Eq. A3, as in Stage 4a.

    Inverts ``p_L = (1 - exp(-2 R T)) / 2``; use the schedule's realised
    ``total_time``, not the nominal ``T``.
    """
    if not 0.0 <= p_L < 0.5:
        return math.nan
    return -math.log1p(-2.0 * p_L) / (2.0 * total_time)
