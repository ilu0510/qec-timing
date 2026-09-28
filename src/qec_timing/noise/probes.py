"""Probe programs that isolate one noise channel so its rate can be measured.

Each probe exercises a single primitive and reports a raw success count, so the
empirical rate can be compared against the target with a binomial confidence
interval. Both repeat the draw ``n_reps`` times *inside* one shot, which is far
cheaper than paying per-shot overhead for every Bernoulli trial.
"""

from __future__ import annotations

from guppylang import guppy
from guppylang.std.array import array
from guppylang.std.builtins import comptime, output
from guppylang.std.qsystem.random import RNG
from guppylang.std.qsystem.utils import get_current_shot
from guppylang.std.quantum import measure_array, qubit

from .primitives import apply_data_noise, flip_record

__all__ = ["build_data_flip_probe", "record_flip_probe"]


def build_data_flip_probe(n_data: int):
    """Apply ``apply_data_noise`` to ``n_data`` fresh qubits, ``n_reps`` times.

    Every qubit starts in |0>, so a measured 1 is exactly an applied X. Emits
    the total flip count and the number of trials.
    """

    @guppy
    def data_flip_probe(p_data: float, n_reps: int, base_seed: int) -> None:
        rng = RNG(base_seed + int(get_current_shot()))
        flips = 0
        rep = 0
        while rep < n_reps:
            data = array(qubit() for _ in range(comptime(n_data)))
            apply_data_noise(rng, data, p_data)
            measured = measure_array(data)
            for i in range(comptime(n_data)):
                if measured[i].read():
                    flips += 1
            rep += 1
        output("flips", flips)
        output("trials", n_reps * comptime(n_data))
        rng.discard()

    return data_flip_probe


@guppy
def record_flip_probe(p_read: float, n_reps: int, base_seed: int) -> None:
    """Apply ``flip_record`` to a False bit ``n_reps`` times.

    ``flip_record(rng, False, p)`` returns True exactly when it flipped, so the
    count of True results estimates ``p_read``.
    """
    rng = RNG(base_seed + int(get_current_shot()))
    flips = 0
    rep = 0
    while rep < n_reps:
        if flip_record(rng, False, p_read):
            flips += 1
        rep += 1
    output("flips", flips)
    output("trials", n_reps)
    rng.discard()
