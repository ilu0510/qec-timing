"""Round scheduling and time bookkeeping.

Shared by the reference model and the Guppy path, like ``layout.py`` and
``ansatz.py``.

**This module contains no physics.** It carries round durations, idling rates
and the clock; it does not know Eqs. 1-3. That separation is deliberate: the
noise equations are implemented twice on purpose (``noise/params.py`` and
``reference/noise.py``) so the reference stays an independent cross-check, and
promoting the schedule to a shared module must not quietly reunify them. Each
side turns a schedule into probabilities with its own implementation.

The rule this module exists to enforce (CLAUDE.md, docs/noise_spec.md section 5):

> **Never loop on a float accumulator.** Accumulating ``t += dt`` while ``t < T``
> silently changes the round count -- with ``dt = 0.1, T = 1.0`` the sum after
> ten steps is ``0.9999999999999999 < 1.0``, giving 11 rounds instead of 10.

So round counts come from ``round(T / dt)`` and the clock is an **integer tick
counter**, converted to physical time only for output.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

__all__ = [
    "DEFAULT_TICK",
    "RoundSchedule",
    "constant_schedule",
    "schedule_from_intervals",
    "ticks_for",
]

DEFAULT_TICK = 1e-3


def ticks_for(dt: float, tick: float = DEFAULT_TICK) -> int:
    """Express ``dt`` as a whole number of ticks, refusing to misrepresent it.

    Raises rather than silently snapping, because a snapped interval would make
    the emitted ``dt`` disagree with the one actually simulated.
    """
    if dt <= 0.0:
        raise ValueError(f"dt must be positive, got {dt}")
    if tick <= 0.0:
        raise ValueError(f"tick must be positive, got {tick}")
    n = round(dt / tick)
    if n == 0:
        raise ValueError(
            f"dt={dt!r} is smaller than one tick ({tick!r}); "
            "use a finer tick to represent it"
        )
    if not math.isclose(n * tick, dt, rel_tol=1e-9, abs_tol=tick * 1e-9):
        raise ValueError(
            f"dt={dt!r} is not an integer multiple of tick={tick!r} "
            f"(nearest is {n} ticks = {n * tick!r}); "
            "choose a compatible tick to keep time bookkeeping exact"
        )
    return n


@dataclass(frozen=True)
class RoundSchedule:
    """Per-round durations and idling rates, with an exact integer clock.

    ``ticks[j]`` is the duration of round ``j`` in ticks and ``lam[j]`` the
    idling rate during it. ``start_ticks[j]`` is the exact time at the start of
    round ``j``, reproducing Eq. A3's ``t_j = sum_{i<j} dt_i`` without drift.
    """

    ticks: tuple[int, ...]
    lam: tuple[float, ...]
    tick: float = DEFAULT_TICK

    def __post_init__(self) -> None:
        if len(self.ticks) != len(self.lam):
            raise ValueError(
                f"ticks and lam must have the same length, got "
                f"{len(self.ticks)} and {len(self.lam)}"
            )
        if not self.ticks:
            raise ValueError("schedule must contain at least one round")
        if any(n <= 0 for n in self.ticks):
            raise ValueError("every round must span a positive number of ticks")
        if self.tick <= 0.0:
            raise ValueError(f"tick must be positive, got {self.tick}")
        if any(x < 0.0 for x in self.lam):
            raise ValueError("idling rates must be non-negative")

    # -- shape ------------------------------------------------------------

    @property
    def n_rounds(self) -> int:
        return len(self.ticks)

    @property
    def is_constant(self) -> bool:
        """True when every round has the same duration.

        The ``t_j = j * dt`` reconstruction used by the emitted metadata is only
        valid in this case.
        """
        return len(set(self.ticks)) == 1

    @property
    def is_static_lambda(self) -> bool:
        """True when the idling rate never changes."""
        return len(set(self.lam)) == 1

    # -- time -------------------------------------------------------------

    @property
    def intervals(self) -> tuple[float, ...]:
        """Round durations ``dt_j`` in physical time units."""
        return tuple(n * self.tick for n in self.ticks)

    @property
    def dt(self) -> float:
        """The single interval, for a constant schedule."""
        if not self.is_constant:
            raise ValueError(
                "schedule has varying intervals; use .intervals instead of .dt"
            )
        return self.ticks[0] * self.tick

    @property
    def start_ticks(self) -> tuple[int, ...]:
        """Exact integer start time of each round."""
        out: list[int] = []
        acc = 0
        for n in self.ticks:
            out.append(acc)
            acc += n
        return tuple(out)

    @property
    def start_times(self) -> tuple[float, ...]:
        """``t_j`` for each round, in physical time units."""
        return tuple(n * self.tick for n in self.start_ticks)

    @property
    def total_ticks(self) -> int:
        return sum(self.ticks)

    @property
    def total_time(self) -> float:
        """Realised total time ``sum_j dt_j``.

        May differ from a nominal ``T``, because ``n_rounds`` is rounded. Use
        this, not the nominal value, whenever a rate per unit time is computed.
        """
        return self.total_ticks * self.tick


def constant_schedule(
    dt: float, T: float, lam: float, *, tick: float = DEFAULT_TICK
) -> RoundSchedule:
    """Fixed-interval schedule with ``n_rounds = round(T / dt)``."""
    if T <= 0.0:
        raise ValueError(f"T must be positive, got {T}")
    if dt > T:
        raise ValueError(
            f"dt={dt} exceeds T={T}: a single round would overrun the total "
            "time budget"
        )
    n_ticks = ticks_for(dt, tick)
    n_rounds = round(T / dt)
    if n_rounds < 1:
        # Unreachable while dt <= T, kept so the invariant is explicit.
        raise ValueError(f"round(T/dt) = {n_rounds} rounds for dt={dt}, T={T}")
    return RoundSchedule(
        ticks=(n_ticks,) * n_rounds, lam=(float(lam),) * n_rounds, tick=tick
    )


def schedule_from_intervals(
    intervals: list[float] | tuple[float, ...],
    lams: list[float] | tuple[float, ...],
    *,
    tick: float = DEFAULT_TICK,
) -> RoundSchedule:
    """Schedule from explicit per-round ``dt_j`` and ``lam_j``.

    Built now so variable-interval work is not blocked later. Nothing in Stages
    0-4 uses it beyond a round-trip test; the burst and adaptive stages are out
    of scope (CLAUDE.md roadmap).
    """
    if len(intervals) != len(lams):
        raise ValueError(
            f"intervals and lams must have the same length, got "
            f"{len(intervals)} and {len(lams)}"
        )
    return RoundSchedule(
        ticks=tuple(ticks_for(dt, tick) for dt in intervals),
        lam=tuple(float(x) for x in lams),
        tick=tick,
    )
