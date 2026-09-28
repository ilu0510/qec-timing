"""Phenomenological noise parameters for the reference model -- Eqs. 1-3.

Equations 1-3 of ``docs/paper_equations.md``:

* Eq. 1  ``p_idle = 1 - exp(-p * lam * dt)``
* Eq. 2  ``p_stab = p``,  ``p_data = 1 - (1 - p_stab)(1 - p_idle)``
* Eq. 3  ``p_read = b_read * p``

``p_data`` is applied as a *single* Bernoulli draw per data qubit per round, not
as two independent flips -- see CLAUDE.md.

This is the **reference** implementation of Eqs. 1-3, deliberately independent
of ``qec_timing.noise.params`` so the cross-check has something to check. Round
scheduling is *not* duplicated: it is pure bookkeeping and lives in the shared
``qec_timing.schedule``, re-exported here for convenience.
"""

from __future__ import annotations

import math

from ..schedule import (
    DEFAULT_TICK,
    RoundSchedule,
    constant_schedule,
    schedule_from_intervals,
    ticks_for,
)

__all__ = [
    "DEFAULT_TICK",
    "p_idle",
    "p_data",
    "p_read",
    "data_probabilities",
    "RoundSchedule",
    "constant_schedule",
    "schedule_from_intervals",
    "ticks_for",
]


def p_idle(p: float, lam: float, dt: float) -> float:
    """Eq. 1: idling fault probability over an interval ``dt``."""
    _check_prob(p, "p")
    if lam < 0.0:
        raise ValueError(f"lam must be non-negative, got {lam}")
    if dt < 0.0:
        raise ValueError(f"dt must be non-negative, got {dt}")
    return -math.expm1(-p * lam * dt)


def p_data(p: float, lam: float, dt: float) -> float:
    """Eq. 2: total data-qubit fault probability for one round.

    ``p_data = 1 - (1 - p_stab)(1 - p_idle)`` with ``p_stab = p``. The ``p_stab``
    term survives at ``lam = 0`` and must never be dropped: it is what creates
    the trade-off in ``dt``.
    """
    return 1.0 - (1.0 - p) * (1.0 - p_idle(p, lam, dt))


def p_read(p: float, b_read: float) -> float:
    """Eq. 3: syndrome-readout flip probability ``p_read = b_read * p``."""
    _check_prob(p, "p")
    if b_read < 0.0:
        raise ValueError(f"b_read must be non-negative, got {b_read}")
    value = b_read * p
    if value > 1.0:
        raise ValueError(f"p_read = b_read * p = {value} exceeds 1")
    return value



def data_probabilities(schedule: RoundSchedule, p: float) -> tuple[float, ...]:
    """``p_data`` for each round of a schedule, via Eqs. 1-2.

    A free function rather than a method on ``RoundSchedule``: the schedule is
    shared between the reference and product paths, and must not carry either
    side's implementation of the noise equations.
    """
    return tuple(
        p_data(p, lam_j, dt_j)
        for lam_j, dt_j in zip(schedule.lam, schedule.intervals, strict=True)
    )


def _check_prob(value: float, name: str) -> None:
    if not 0.0 <= value <= 1.0:
        raise ValueError(f"{name} must lie in [0, 1], got {value}")
