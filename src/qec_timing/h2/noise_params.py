"""H2-Emulator noise parameters and the idle (dephasing) conversion. Stdlib only.

Source: Quantinuum emulator docs, read 2026-10-10:
* H2-1E defaults ("H2-1E emulators from 2025-07-16"),
  https://docs.quantinuum.com/systems/user_guide/emulator_user_guide/emulators/h2_emulators.html
* noise model,
  https://docs.quantinuum.com/systems/user_guide/emulator_user_guide/noise_model.html

Unit convention (quoted from the noise-model page):

* Linear: "Pauli Z is applied ... according to the probability of
  linear_dephasing_rate multiplied by duration (fd) where linear_dephasing_rate
  is per second (s^-1), and duration is in units of seconds."
  -> p_lin = r * d, r = 0.0028 s^-1, d in seconds.
* Incoherent quadratic (the stabilizer-simulator default, coherent_dephasing
  False): "Pauli-Z is applied ... according to the probability, sin(fd/2)^2,
  where frequency is equal to quadratic_dephasing_rate multiplied by the
  coherent_to_incoherent_factor (all in units of 2 pi radians per second)."
  -> the factor scales the FREQUENCY, not the probability.

The one ambiguity: "units of 2 pi radians per second" read literally makes the
numeric rate a frequency in cycles/s, so the angle is 2 pi * 0.043 * 2 * d
(convention "cycles", the default here). Read as a loose label, the printed
formula sin(fd/2)^2 takes f = 0.086 rad/s directly (convention "radians").
The two differ by (2 pi)^2 ~ 39 in the quadratic term. "cycles" is the
default because it follows the units statement literally; it is also the
conservative (larger-noise) choice. Every function takes `convention` so the
other reading can be compared, and the report gives both.

How the two terms combine is not stated beyond "we apply both sorts of noise
simultaneously". Here they are composed as two independent Z channels,
p = q + l - 2 q l. The twin's DephasingH2 adds them (and clips), so it is
higher by exactly 2 q l: about 3e-7 at the d=3 dt* (44 rounds of 2 ms, q ~ 5.6e-4,
l ~ 2.5e-4), but 0.5% of p_z at the far end of the dt grid (p_z ~ 1e-2).

Timing convention: dt is the syndrome PERIOD (start of one round to the next),
so dt >= t_round. Data qubits idle and get transported for the whole period,
and the emulator charges memory noise during both, so the idle duration per
round is dt. In twin units (time unit = t_round), dt_rounds = dt / t_round.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

__all__ = ["H2_1E", "H2Noise", "idle_pz", "twin_dephasing", "DEFAULT_T_ROUND"]

DEFAULT_T_ROUND = 2e-3  # seconds per syndrome round (a parameter, not an H2 figure)
CONVENTIONS = ("cycles", "radians")


@dataclass(frozen=True)
class H2Noise:
    p1: float
    p2: float
    p_meas_0: float          # bit flip when the true outcome is 0
    p_meas_1: float          # bit flip when the true outcome is 1
    p_init: float
    quadratic_dephasing_rate: float
    linear_dephasing_rate: float     # s^-1
    coherent_to_incoherent_factor: float


H2_1E = H2Noise(p1=1.89e-5, p2=1.05e-3, p_meas_0=6e-4, p_meas_1=1.39e-3, p_init=4e-5,
                quadratic_dephasing_rate=0.043, linear_dephasing_rate=0.0028,
                coherent_to_incoherent_factor=2.0)


def _omega(noise: H2Noise, convention: str) -> float:
    """Incoherent-quadratic angular frequency in rad/s (factor included)."""
    if convention not in CONVENTIONS:
        raise ValueError(f"convention must be one of {CONVENTIONS}")
    f = noise.quadratic_dephasing_rate * noise.coherent_to_incoherent_factor
    return 2 * math.pi * f if convention == "cycles" else f


def idle_pz(duration_s: float, noise: H2Noise = H2_1E, convention: str = "cycles") -> float:
    """P(Z) on one idling qubit over `duration_s` seconds (stabilizer default model)."""
    if duration_s < 0:
        raise ValueError("duration must be non-negative")
    q = math.sin(0.5 * _omega(noise, convention) * duration_s) ** 2
    lin = min(noise.linear_dephasing_rate * duration_s, 1.0)
    return q + lin - 2 * q * lin


def twin_dephasing(t_round: float = DEFAULT_T_ROUND, noise: H2Noise = H2_1E,
                   convention: str = "cycles") -> tuple[float, float]:
    """(f_quad, f_lin) for the Stage 6 twin's DephasingH2 with time unit = t_round.

    The twin applies p_z(dt) = sin^2(0.5 f_quad dt) + f_lin dt with dt in its
    own time unit. With dt measured in rounds: f_quad = omega * t_round
    (radians per round), f_lin = r * t_round (per round).
    """
    if t_round <= 0:
        raise ValueError("t_round must be positive")
    return _omega(noise, convention) * t_round, noise.linear_dephasing_rate * t_round
