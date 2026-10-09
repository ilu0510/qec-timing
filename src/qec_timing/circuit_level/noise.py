"""Circuit-level noise configuration for the Stage 6a Stim twin.

This deliberately departs from the Stage-1 physics (which has no gate noise).
Nothing here is shared with ``qec_timing.noise`` or ``qec_timing.reference``.

Channels
--------
* ``p1``: ``DEPOLARIZE1(p1)`` after every single-qubit gate.
* ``p2``: ``DEPOLARIZE2(p2)`` after every CNOT.
* ``p_init``: ``X_ERROR(p_init)`` after every Z-basis reset (qubit allocation).
* ``readout``: classical flip of every measurement record, possibly asymmetric
  (see ``ReadoutNoise``).
* ``idle``: one of the idle channels below, applied to every data qubit once
  per round for that round's interval ``dt``, before extraction.

Idle channels
-------------
``DepolarizingKOps`` ("depolarizing_kops", the colleague's Selene route):
``k`` charged single-qubit depolarizing ops per unit time, each of strength
``p_op``. ``n`` compositions of ``DEPOLARIZE1(p)`` are again depolarizing,
because the channel is ``rho -> (1 - 4p/3) rho + (4p/3) I/2`` and the
coefficient multiplies::

    p_n = 3/4 * (1 - (1 - 4p/3) ** n)

``n = k * dt`` is applied as ONE ``DEPOLARIZE1(p_n)``. For integer ``n`` this is
exact; for non-integer ``n`` it is the unique continuous-time depolarizing
semigroup that agrees with it at every integer, which is what "k ops per unit
time" means between whole ops. Selene itself can only charge a whole number of
ops; a Selene run with non-integer ``k * dt`` has to round, and that rounding
belongs to 6b, not here.

``DephasingH2`` ("dephasing_h2", the H2 emulator's idle model)::

    p_z(dt) = sin^2(0.5 * f_quad * dt) + f_lin * dt,   clipped to [0, 0.5]

applied as ``Z_ERROR(p_z)``. ``f_quad`` and ``f_lin`` are rates in units of
1 / (the time unit of ``dt``); ``f_quad`` is an angular rate (the argument of
``sin`` is in radians). No hardware values are built in.

``BitFlip`` ("bitflip"): ``X_ERROR(p)`` per round, independent of ``dt``. Exists
only to reduce the circuit-level model to the phenomenological one (the bridge
test to Stage 4a); the caller passes the per-round probability.

Native charging
---------------
Selene charges ``p_1q``/``p_2q`` per COMPILED native op, not per source gate.
The counts below were MEASURED in Stage 6b from the error model's own input
(``CircuitExtractor.get_error_model_input()``, selene-sim 0.3.2,
``DepolarizingErrorModel`` on the Stim backend, one noiseless traced shot of
``program.py``), and confirmed with the model's ``gates_1q`` / ``gates_2q``
counters: every ``Rxy`` and every ``Rz`` is charged ``p_1q`` (``DEPOLARIZE1``),
every ``Rzz`` ``p_2q`` (``DEPOLARIZE2``). See ``docs/stage6b_report.md``.

* ``h``: ``Rxy`` + ``Rz`` = **2** ops (agrees with notebook 03, section 2).
* ``cx(control, target)``: target ``Rxy``; ``Rzz``; control ``Rz``; target
  ``Rxy``, ``Rz``. So the TARGET gets 1 op before and 2 after, the CONTROL 0
  before and 1 after -- regardless of which of them is the ancilla.
  (Stage 6a assumed 1 before + 1 after on both, from notebook 04c's "2 per data
  qubit", which is the average over the control and target roles. Wrong per
  qubit; corrected here.)
* ``x``, ``z`` (preparation corrections, injected flips): **1** op each.
* ``phased_x(q, 0, 0)`` (variant A idle op): **1** op (notebook 03, section 1).

With ``native_charging=True`` the builder therefore emits:

* ``H``: one ``DEPOLARIZE1(compose(p1, 2))``. Exact, because the single-qubit
  depolarizing channel is invariant under conjugation by any single-qubit
  Clifford, so a native op's noise can be moved past the native rotations and
  merged with its neighbours.
* ``CX``: ``DEPOLARIZE1(p1)`` on the target; ``CX``; ``DEPOLARIZE2(p2)``;
  ``DEPOLARIZE1(p1)`` on the control and ``DEPOLARIZE1(compose(p1, 2))`` on the
  target. Both depolarizing channels are invariant under local Cliffords, so
  only the before/after split relative to the entangling ``Rzz`` matters.
* idle (variant A): unchanged, each idle op is already one native op.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Literal, Union

import stim

__all__ = [
    "BitFlip",
    "CircuitNoise",
    "CX_CONTROL_1Q_OPS",
    "CX_TARGET_1Q_OPS",
    "PAULI_GATE_1Q_OPS",
    "DephasingH2",
    "DepolarizingKOps",
    "H_NATIVE_1Q_OPS",
    "IDLE_OP_NATIVE_1Q_OPS",
    "IdleChannel",
    "ReadoutNoise",
    "compose_depolarize1",
    "dephasing_h2_probability",
]

Basis = Literal["Z", "X"]

# Native op counts, measured from Selene's error-model input (module docstring).
H_NATIVE_1Q_OPS = 2          # h -> Rxy, Rz
IDLE_OP_NATIVE_1Q_OPS = 1    # phased_x(q, 0, 0) -> Rxy
PAULI_GATE_1Q_OPS = 1        # x -> Rxy, z -> Rz
CX_TARGET_1Q_OPS = (1, 2)    # (before, after) the Rzz
CX_CONTROL_1Q_OPS = (0, 1)   # (before, after) the Rzz

_MAX_DEPOLARIZE1 = 0.75


def compose_depolarize1(p: float, n: float) -> float:
    """Single ``DEPOLARIZE1`` equivalent to ``n`` compositions of ``DEPOLARIZE1(p)``.

    ``p_n = 3/4 * (1 - (1 - 4p/3) ** n)``; ``n`` may be any non-negative real.
    """
    _check_prob(p, "p", _MAX_DEPOLARIZE1)
    if not math.isfinite(n) or n < 0.0:
        raise ValueError(f"composition count must be finite and >= 0, got {n}")
    if n == 0.0:
        return 0.0
    if n == 1.0:
        return p  # exact, rather than p up to rounding
    return 0.75 * (1.0 - (1.0 - 4.0 * p / 3.0) ** n)


def dephasing_h2_probability(dt: float, f_quad: float, f_lin: float) -> float:
    """``p_z(dt) = sin^2(0.5 f_quad dt) + f_lin dt``, clipped to [0, 0.5]."""
    if dt < 0.0:
        raise ValueError(f"dt must be non-negative, got {dt}")
    value = math.sin(0.5 * f_quad * dt) ** 2 + f_lin * dt
    return min(max(value, 0.0), 0.5)


# ---------------------------------------------------------------------------
# Idle channels
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DepolarizingKOps:
    """Variant A: ``k`` charged ``DEPOLARIZE1(p_op)`` per unit time."""

    k: float
    p_op: float

    name = "depolarizing_kops"
    default_basis = "Z"

    def __post_init__(self) -> None:
        if not math.isfinite(self.k) or self.k < 0.0:
            raise ValueError(f"k must be finite and >= 0, got {self.k}")
        _check_prob(self.p_op, "p_op", _MAX_DEPOLARIZE1)

    def probability(self, dt: float) -> float:
        """Composed ``DEPOLARIZE1`` probability for one interval ``dt``."""
        if dt < 0.0:
            raise ValueError(f"dt must be non-negative, got {dt}")
        return compose_depolarize1(self.p_op, self.k * dt)

    def append(self, circuit: stim.Circuit, targets: list[int], dt: float) -> None:
        p = self.probability(dt)
        if p > 0.0:
            circuit.append("DEPOLARIZE1", targets, p)


@dataclass(frozen=True)
class DephasingH2:
    """Variant B: ``Z_ERROR(p_z(dt))`` with the H2 emulator's idle formula."""

    f_quad: float
    f_lin: float

    name = "dephasing_h2"
    default_basis = "X"

    def __post_init__(self) -> None:
        for value, label in ((self.f_quad, "f_quad"), (self.f_lin, "f_lin")):
            if not math.isfinite(value):
                raise ValueError(f"{label} must be finite, got {value}")
        if self.f_lin < 0.0:
            raise ValueError(f"f_lin must be >= 0, got {self.f_lin}")

    def probability(self, dt: float) -> float:
        return dephasing_h2_probability(dt, self.f_quad, self.f_lin)

    def append(self, circuit: stim.Circuit, targets: list[int], dt: float) -> None:
        p = self.probability(dt)
        if p > 0.0:
            circuit.append("Z_ERROR", targets, p)


@dataclass(frozen=True)
class BitFlip:
    """Phenomenological bridge: ``X_ERROR(p)`` once per round, any ``dt``."""

    p: float

    name = "bitflip"
    default_basis = "Z"

    def __post_init__(self) -> None:
        _check_prob(self.p, "p", 1.0)

    def probability(self, dt: float) -> float:
        return self.p

    def append(self, circuit: stim.Circuit, targets: list[int], dt: float) -> None:
        if self.p > 0.0:
            circuit.append("X_ERROR", targets, self.p)


IdleChannel = Union[DepolarizingKOps, DephasingH2, BitFlip]


# ---------------------------------------------------------------------------
# Readout
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ReadoutNoise:
    """Classical flip of the measurement record, possibly asymmetric.

    ``p01`` = P(report 1 | true outcome 0), ``p10`` = P(report 0 | true 1).

    Stim's ``M(p)`` is symmetric and its records are combined only by XOR, so an
    outcome-dependent flip ``m -> m XOR (e AND m)`` cannot be written inside a
    Stim circuit. When ``p01 != p10`` the memory therefore uses two circuits
    that differ only in their ``M`` arguments:

    * SAMPLING: noiseless ``M``; measurement records are sampled and then
      flipped in numpy, each bit with ``p01`` or ``p10`` according to the bit
      actually recorded. This is EXACT for the classical asymmetric channel.
      Detectors and observables are then recomputed with Stim's
      measurements-to-detections converter.
    * DECODING: ``M(p_avg)`` with ``p_avg = (p01 + p10) / 2`` builds the
      detector error model. Only the decoder's edge weights are approximate,
      as a Pauli-twirled (symmetric) model of the true channel.

    One side effect, documented rather than modelled: Stim feedback (the
    preparation corrections) sees the record BEFORE the classical flip. Those
    corrections only fix the signs of the complementary checks, which are never
    decoded, so neither detectors nor the observable are affected.
    """

    p01: float = 0.0
    p10: float = 0.0

    def __post_init__(self) -> None:
        _check_prob(self.p01, "p01", 1.0)
        _check_prob(self.p10, "p10", 1.0)

    @classmethod
    def symmetric(cls, p: float) -> "ReadoutNoise":
        return cls(p, p)

    @property
    def is_symmetric(self) -> bool:
        return self.p01 == self.p10

    @property
    def average(self) -> float:
        return 0.5 * (self.p01 + self.p10)


@dataclass(frozen=True)
class CircuitNoise:
    """All circuit-level noise parameters for one memory experiment."""

    idle: IdleChannel
    p1: float = 0.0
    p2: float = 0.0
    p_init: float = 0.0
    readout: ReadoutNoise = field(default_factory=ReadoutNoise)

    def __post_init__(self) -> None:
        _check_prob(self.p1, "p1", _MAX_DEPOLARIZE1)
        _check_prob(self.p2, "p2", 15.0 / 16.0)
        _check_prob(self.p_init, "p_init", 1.0)
        if not isinstance(self.idle, (DepolarizingKOps, DephasingH2, BitFlip)):
            raise TypeError(f"unknown idle channel {self.idle!r}")
        if not isinstance(self.readout, ReadoutNoise):
            raise TypeError("readout must be a ReadoutNoise")

    @classmethod
    def depolarizing_kops(cls, *, k: float, p_op: float | None = None,
                          p1: float = 0.0, **kwargs) -> "CircuitNoise":
        """Variant A. ``p_op`` defaults to ``p1``: Selene charges an idle op ``p_1q``."""
        return cls(DepolarizingKOps(k, p1 if p_op is None else p_op), p1=p1, **kwargs)

    @classmethod
    def dephasing_h2(cls, *, f_quad: float, f_lin: float, **kwargs) -> "CircuitNoise":
        """Variant B."""
        return cls(DephasingH2(f_quad, f_lin), **kwargs)


def _check_prob(value: float, name: str, upper: float) -> None:
    if not (isinstance(value, (int, float)) and 0.0 <= value <= upper):
        raise ValueError(f"{name} must lie in [0, {upper}], got {value!r}")
