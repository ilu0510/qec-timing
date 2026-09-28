"""Noise parameters for the product (Guppy) path -- Eqs. 1-3.

Deliberately an **independent** implementation of the same equations as
``qec_timing.reference.noise``. The reference model exists to cross-check the
Selene implementation (CLAUDE.md); if the product path imported its physics from
the reference, a mistake in Eqs. 1-2 would appear identically on both sides and
the cross-check could never see it. ``test_stage2_guppy_noise.py`` asserts the
two implementations agree, which catches drift without creating the coupling.

``exp()`` is not available inside Guppy, so these run in Python and the results
are passed into the Guppy program as runtime float arguments. See
``docs/noise_spec.md``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ..schedule import RoundSchedule

__all__ = [
    "NoiseParams",
    "idle_probability",
    "data_probability",
    "read_probability",
    "schedule_data_probabilities",
]


def idle_probability(p: float, lam: float, dt: float) -> float:
    """Eq. 1: ``p_idle = 1 - exp(-p * lam * dt)``.

    Uses ``expm1`` for accuracy when the exponent is small, which is the regime
    the short-interval end of a dt sweep lives in.
    """
    _require_probability(p, "p")
    _require_non_negative(lam, "lam")
    _require_non_negative(dt, "dt")
    return -math.expm1(-p * lam * dt)


def data_probability(p: float, lam: float, dt: float) -> float:
    """Eq. 2: ``p_data = 1 - (1 - p)(1 - p_idle)``, with ``p_stab = p``.

    The ``p_stab = p`` term survives at ``lam = 0`` and must never be dropped.
    """
    return 1.0 - (1.0 - p) * (1.0 - idle_probability(p, lam, dt))


def read_probability(p: float, b_read: float) -> float:
    """Eq. 3: ``p_read = b_read * p``, independent of ``dt``."""
    _require_probability(p, "p")
    _require_non_negative(b_read, "b_read")
    value = b_read * p
    if value > 1.0:
        raise ValueError(f"p_read = b_read * p = {value} exceeds 1")
    return value


@dataclass(frozen=True)
class NoiseParams:
    """The four numbers the Guppy program needs for a static-lambda run."""

    p: float
    lam: float
    dt: float
    b_read: float

    @property
    def p_data(self) -> float:
        return data_probability(self.p, self.lam, self.dt)

    @property
    def p_read(self) -> float:
        return read_probability(self.p, self.b_read)

    def n_rounds(self, T: float) -> int:
        """``round(T / dt)`` -- never a float accumulator (docs/noise_spec.md)."""
        if T <= 0.0:
            raise ValueError(f"T must be positive, got {T}")
        n = round(T / self.dt)
        if n < 1:
            raise ValueError(f"dt={self.dt} exceeds T={T}: round(T/dt) = {n}")
        return n

    def realised_T(self, T: float) -> float:
        """Actual elapsed time ``n_rounds * dt``, which may differ from ``T``."""
        return self.n_rounds(T) * self.dt


def _require_probability(value: float, name: str) -> None:
    if not 0.0 <= value <= 1.0:
        raise ValueError(f"{name} must lie in [0, 1], got {value}")


def _require_non_negative(value: float, name: str) -> None:
    if value < 0.0:
        raise ValueError(f"{name} must be non-negative, got {value}")


def schedule_data_probabilities(
    schedule: RoundSchedule, p: float
) -> tuple[float, ...]:
    """``p_data`` for each round of a schedule, via Eqs. 1-2 (product path).

    Mirrors ``reference.noise.data_probabilities`` but uses this module's
    independent implementation; a test asserts the two agree.
    """
    return tuple(
        data_probability(p, lam_j, dt_j)
        for lam_j, dt_j in zip(schedule.lam, schedule.intervals, strict=True)
    )
