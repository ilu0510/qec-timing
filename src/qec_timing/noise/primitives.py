"""Guppy noise primitives: Bernoulli draws, data faults, record flips.

All of these take the RNG as a **borrowed** parameter (`rng: RNG`, not
`rng: RNG @owned`), so a single instance threads through the whole program
without being consumed. That is required, not stylistic: `RNG(seed)` maps onto
Selene's global runtime PRNG, so a second instance would reseed the first.
See `docs/noise_spec.md` section 7.
"""

from __future__ import annotations

from guppylang import guppy
from guppylang.std.array import array
from guppylang.std.qsystem.random import RNG
from guppylang.std.quantum import qubit, x

__all__ = ["bernoulli", "apply_data_noise", "flip_record"]

N = guppy.nat_var("N")


@guppy
def bernoulli(rng: RNG, p: float) -> bool:
    """Draw ``Bernoulli(p)``.

    ``random_float()`` is uniform on ``[0, 1)``, so ``< p`` is true with
    probability exactly ``p``; ``p = 0`` can never fire and ``p = 1`` always
    does.
    """
    return rng.random_float() < p


@guppy
def apply_data_noise(rng: RNG, data: array[qubit, N], p_data: float) -> None:
    """Apply X to each data qubit with probability ``p_data``.

    Exactly **one** Bernoulli draw per qubit per round. ``p_data`` already
    combines the idling and measurement-induced mechanisms via Eq. 2; drawing
    once for each and applying X twice would let the flips cancel and give a
    different channel (see `docs/noise_spec.md` section 3).
    """
    for i in range(N):
        if bernoulli(rng, p_data):
            x(data[i])


@guppy
def flip_record(rng: RNG, bit: bool, p_read: float) -> bool:
    """Flip a recorded classical bit with probability ``p_read``.

    Operates on the **record only**; the qubit is untouched, so a read-out fault
    cannot propagate into the next round's syndrome.
    """
    if bernoulli(rng, p_read):
        return not bit
    return bit
