"""Analytic logical-error-rate ansatz from Haug, Bharti & Aolita (2026).

Pure analytic formulas. Equation numbers refer to ``docs/paper_equations.md``.

This module is deliberately self-contained: it must not import from
``qec_timing.reference`` or from any Selene/Guppy code, so that the analytic
model stays an independent yardstick for both the reference sampler and the
Selene implementation.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

__all__ = [
    "alpha",
    "K_alpha",
    "logical_rate",
    "optimal_interval",
    "optimal_rate",
    "gamma_opt",
    "p_L_from_rates",
    "p_L_constant_interval",
]

# Fitted constants from the paper (docs/paper_equations.md, "Fitted constants").
BETA_DEFAULT = 2.0
G_DEFAULT = 0.8
P_TH_BREAD_1 = 0.029


def alpha(d: int | ArrayLike, g: float = G_DEFAULT) -> NDArray[np.float64]:
    """Leading idling-sensitive exponent, Eq. D2: ``alpha = g (d + 1) / 2``."""
    return np.asarray(g * (np.asarray(d, dtype=float) + 1.0) / 2.0, dtype=float)


def K_alpha(a: ArrayLike) -> NDArray[np.float64]:
    """Eq. D8: ``K_alpha = alpha**alpha / (alpha - 1)**(alpha - 1)``.

    Computed in log space; ``K_alpha -> e * alpha`` for large alpha.
    """
    a = np.asarray(a, dtype=float)
    if np.any(a <= 1.0):
        raise ValueError("K_alpha requires alpha > 1 (i.e. g (d + 1) > 2)")
    return np.exp(a * np.log(a) - (a - 1.0) * np.log(a - 1.0))


def logical_rate(
    dt: ArrayLike,
    *,
    d: int,
    p: float,
    lam: float,
    A: float,
    p_th: float = P_TH_BREAD_1,
    beta: float = BETA_DEFAULT,
    g: float = G_DEFAULT,
) -> NDArray[np.float64]:
    """Eq. 5 / D1: logical error rate per unit physical time.

    ``R = (1/dt) (A / d**beta) (p / p_th)**((d+1)/2) (1 + lam dt)**(g (d+1)/2)``
    """
    dt = np.asarray(dt, dtype=float)
    if np.any(dt <= 0.0):
        raise ValueError("dt must be positive")
    a = alpha(d, g)
    prefactor = (A / d**beta) * (p / p_th) ** ((d + 1) / 2.0)
    return prefactor * (1.0 + lam * dt) ** a / dt


def optimal_interval(
    *, d: int, lam: float, g: float = G_DEFAULT
) -> float:
    """Eq. 6 / D6: ``dt* = 2 / (lam (g(d+1) - 2)) = 1 / (lam (alpha - 1))``.

    Requires ``alpha > 1``; otherwise Eq. 5 has no interior minimum.
    """
    if lam <= 0.0:
        raise ValueError("optimal_interval requires lam > 0")
    a = float(alpha(d, g))
    if a <= 1.0:
        raise ValueError(
            f"no interior optimum: alpha = {a:.4f} <= 1 for d={d}, g={g}"
        )
    return 1.0 / (lam * (a - 1.0))


def optimal_rate(
    *,
    d: int,
    p: float,
    lam: float,
    A: float,
    p_th: float = P_TH_BREAD_1,
    beta: float = BETA_DEFAULT,
    g: float = G_DEFAULT,
) -> float:
    """Eq. D7: ``R*(lam) = C K_alpha lam`` with ``C`` from Eq. D4."""
    a = alpha(d, g)
    C = (A / d**beta) * (p / p_th) ** ((d + 1) / 2.0)
    return float(C * K_alpha(a) * lam)


def gamma_opt(
    dt: ArrayLike, *, d: int, lam: float, g: float = G_DEFAULT
) -> NDArray[np.float64]:
    """Eq. 7 / D9, exact and without approximation:

    ``Gamma_opt = R(dt) / R(dt*) = (1 + lam dt)**alpha / (lam dt K_alpha)``

    Independent of ``A``, ``p``, ``p_th`` and ``beta``, which cancel in the
    ratio. Equals 1 at ``dt = dt*`` and is > 1 everywhere else.
    """
    dt = np.asarray(dt, dtype=float)
    if np.any(dt <= 0.0):
        raise ValueError("dt must be positive")
    if lam <= 0.0:
        raise ValueError("gamma_opt requires lam > 0")
    a = alpha(d, g)
    return (1.0 + lam * dt) ** a / (lam * dt * K_alpha(a))


def p_L_from_rates(rates: ArrayLike, intervals: ArrayLike) -> float:
    """Eq. A3: map a schedule of rates to a total logical error probability.

    ``p_L = 1/2 (1 - exp(-2 sum_j R(dt_j, t_j) dt_j))``

    Saturates at 1/2 for large accumulated rate and reduces to
    ``sum_j R_j dt_j`` when that sum is small.
    """
    rates = np.asarray(rates, dtype=float)
    intervals = np.asarray(intervals, dtype=float)
    if rates.shape != intervals.shape:
        raise ValueError("rates and intervals must have the same shape")
    total = float(np.sum(rates * intervals))
    return 0.5 * (1.0 - np.exp(-2.0 * total))


def p_L_constant_interval(
    dt: float,
    *,
    n_rounds: int,
    d: int,
    p: float,
    lam: float,
    A: float,
    p_th: float = P_TH_BREAD_1,
    beta: float = BETA_DEFAULT,
    g: float = G_DEFAULT,
) -> float:
    """Eq. 5 mapped through Eq. A3 for a constant interval.

    Uses ``T = n_rounds * dt`` rather than a nominal T, because the round count
    is ``round(T_nominal / dt)`` and need not divide T exactly.
    """
    R = float(
        logical_rate(dt, d=d, p=p, lam=lam, A=A, p_th=p_th, beta=beta, g=g)
    )
    return p_L_from_rates(np.full(n_rounds, R), np.full(n_rounds, dt))
