"""Minimal surface-code test program for H2-Emulator via QIR. Guppy 0.21.x ONLY.

Run under ``.venv-h2`` (guppylang 0.21.11, hugr-qir 0.1.2). Shaped for
hugr-qir's limits:

* no unbounded loops: the entrypoint is ``@guppy.comptime`` (traced), so every
  Python loop is unrolled at compile time into straight-line code;
* no functions taking or returning arrays, no non-comptime arrays: qubits live
  in a Python list inside the traced function, never in a Guppy array;
* no dynamic qubit allocation: all data qubits and the single ancilla are
  allocated once at the top, and the ancilla is reused via ``project_z`` +
  ``reset``;
* outputs via ``result()``; the entrypoint takes no arguments, so every
  parameter (distance, rounds, injection probability, base_seed) is
  compile-time.

Noise injection: ONE ``RNG`` per program, seeded ``base_seed + get_current_shot()``
(CLAUDE.md); each round, before extraction, each data qubit gets Z with
compile-time probability ``p_inject``, drawn inside a small non-traced
``@guppy`` helper, because traced code cannot branch on a runtime value. This
mimics dephasing idle noise.

The Bernoulli is ``random_int_bounded(BOUND) < round(p * BOUND)``, NOT
``random_float() < p``: hugr-qir 0.1.2 lowers RandomInt, RandomIntBounded,
NewRNGContext and GetCurrentShot to QIR, but has no lowering for RandomFloat
("Unknown op: RandomFloat"). This is the integer fallback CLAUDE.md names. The
effective probability is ``effective_p_inject(p)`` = round(p * BOUND) / BOUND,
a quantisation of at most 1 / (2 BOUND) ~ 5e-7.

``checks="Z"`` (the Stage 7a spec): data in |0>, Z-check extraction. Note that
Z-checks CANNOT see injected Z faults (Stage 6a test 4): the syndromes carry
only the emulator's own gate/measurement noise. ``checks="X"``: data in |+>,
X-check extraction, which does detect the injected Z faults.

Result tags: ``s{round}_{check}`` (syndrome bit), ``m{i}`` (final data
measurement, in the check basis).
"""

from __future__ import annotations

from guppylang import guppy
from guppylang.std.builtins import comptime, result
from guppylang.std.qsystem.random import RNG
from guppylang.std.qsystem.utils import get_current_shot
from guppylang.std.quantum import cx, discard, h, measure, project_z, qubit, reset, z

from .geometry import x_checks_hook_safe, z_checks_hook_safe

__all__ = ["BOUND", "build_program", "effective_p_inject"]

BOUND = 1 << 20


def effective_p_inject(p: float) -> float:
    return round(p * BOUND) / BOUND


def build_program(d: int = 3, n_rounds: int = 2, p_inject: float = 0.05,
                  base_seed: int = 7_000_000, checks: str = "Z"):
    """Return the Guppy entrypoint (a traced function with no arguments)."""
    if checks not in ("Z", "X"):
        raise ValueError("checks must be 'Z' or 'X'")
    if not 0.0 <= p_inject <= 1.0:
        raise ValueError("p_inject must be a probability")
    supports = z_checks_hook_safe(d) if checks == "Z" else x_checks_hook_safe(d)
    threshold = round(p_inject * BOUND)
    n_data = d * d
    x_basis = checks == "X"

    @guppy
    def maybe_z(rng: RNG, q: qubit) -> None:
        if rng.random_int_bounded(comptime(BOUND)) < comptime(threshold):
            z(q)

    @guppy.comptime
    def main() -> None:
        rng = RNG(base_seed + get_current_shot())
        data = [qubit() for _ in range(n_data)]
        anc = qubit()
        if x_basis:
            for q in data:
                h(q)
        for r in range(n_rounds):
            for q in data:
                maybe_z(rng, q)
            for ci, support in enumerate(supports):
                if x_basis:
                    h(anc)
                    for i in support:
                        cx(anc, data[i])
                    h(anc)
                else:
                    for i in support:
                        cx(data[i], anc)
                result(f"s{r}_{ci}", project_z(anc))
                reset(anc)
        for i, q in enumerate(data):
            if x_basis:
                h(q)
            result(f"m{i}", measure(q))
        discard(anc)
        rng.discard()

    return main
