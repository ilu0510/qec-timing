"""PLACEHOLDER extraction circuit -- stands in for teammate #1's circuit.

This is **not** the production syndrome-extraction circuit. It is the simplest
thing that is provably equivalent to the reference model's ideal MPP: for each
Z-check, allocate an ancilla in |0>, CNOT each data qubit of the check's support
onto it, and measure it. All gates and the ancilla measurement are **ideal**;
the only noise is what `qec_timing.noise` injects explicitly.

The check geometry is imported from `qec_timing.layout`, shared with the
reference model so the two paths cannot drift apart. (Through Stage 2 this lived
in `reference/layout.py`, which made product code depend on the cross-check
package; Stage 3 promoted it to a shared module alongside `ansatz.py`.)

Everything about the geometry is **compile-time**: array sizes and loop bounds
must be comptime in Guppy, so a given distance produces one compiled binary.
`p_data`, `p_read`, `n_rounds` and `base_seed` are **runtime** arguments, so a
whole sweep reuses a single binary (Stage 0, docs/interfaces.md section (c)).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from guppylang import guppy
from guppylang.std.array import array
from guppylang.std.builtins import comptime, nat, output
from guppylang.std.qsystem.random import RNG
from guppylang.std.qsystem.utils import get_current_shot
from guppylang.std.quantum import cx, measure, measure_array, qubit, x

from ..noise.primitives import apply_data_noise, flip_record
from ..layout import RotatedSurfaceCodeZSector
from ..wire import WORD_BITS, EmitFormat

__all__ = ["EmitFormat", "MemoryProgram", "build_memory_program", "WORD_BITS"]

_N = guppy.nat_var("_STAGE2_N")


@dataclass
class MemoryProgram:
    """A compiled-per-distance Guppy memory experiment."""

    layout: RotatedSurfaceCodeZSector
    emit: EmitFormat
    program: object
    n_words: int

    @property
    def n_qubits(self) -> int:
        """Data qubits plus one live ancilla (ancillas are used one at a time)."""
        return self.layout.n_data + 1

    def emulator(self, **kwargs):
        """Stim-backed emulator instance for this program."""
        return self.program.emulator(n_qubits=self.n_qubits, **kwargs).stabilizer_sim()

    def run(
        self,
        *,
        p_data: float,
        p_read: float,
        n_rounds: int,
        base_seed: int,
        dt: float = 1.0,
        shots: int = 1,
        inject_round: int = -1,
        inject_qubit: int = -1,
        n_processes: int = 1,
    ):
        """Run the experiment. Injection is off unless ``inject_qubit >= 0``."""
        emulator = self.emulator().with_shots(shots)
        if n_processes > 1:
            emulator = emulator.with_n_processes(n_processes)
        return emulator.run(
            p_data=p_data,
            p_read=p_read,
            dt=dt,
            n_rounds=n_rounds,
            base_seed=base_seed,
            inject_round=inject_round,
            inject_qubit=inject_qubit,
        )

    def detection_events(self, result, n_rounds: int):
        """Unpack every shot into a ``(n_rounds + 1, n_checks)`` bool array."""
        from .reader import unpack_detection_events

        return [
            unpack_detection_events(
                shot.collate_tags(),
                n_checks=self.layout.n_checks,
                n_rounds=n_rounds,
                n_words=self.n_words,
                emit=self.emit,
            )
            for shot in result.results
        ]


def _padded_supports(layout: RotatedSurfaceCodeZSector) -> tuple[list[int], int]:
    """Flatten check supports into a fixed-width table padded with -1."""
    max_weight = max(len(s) for s in layout.z_checks)
    flat: list[int] = []
    for support in layout.z_checks:
        indices = [layout.data_index[q] for q in support]
        flat.extend(indices + [-1] * (max_weight - len(indices)))
    return flat, max_weight


def build_memory_program(
    layout: RotatedSurfaceCodeZSector, emit: EmitFormat = "packed"
) -> MemoryProgram:
    """Build the memory-experiment program for one distance and wire format.

    The returned program takes ``(p_data, p_read, dt, n_rounds, base_seed)`` as
    runtime arguments and emits, per shot:

    * ``det``  -- detection events, ``n_rounds + 1`` blocks (format-dependent)
    * ``obs``  -- logical Z parity of the final data readout
    * ``dt``, ``p_data``, ``p_read`` -- per-round metadata; constant under a
      static schedule, so emitted once per shot, with ``t_j = j * dt``
    * ``n_rounds``, ``n_checks``, ``n_words`` -- so the stream can be reshaped

    The full schema is in ``docs/interfaces.md``.
    """
    if emit not in ("packed", "sparse", "dense", "none"):
        raise ValueError(f"unknown emit format {emit!r}")

    n_data = layout.n_data
    n_checks = layout.n_checks
    supports, max_weight = _padded_supports(layout)
    n_words = math.ceil(n_checks / WORD_BITS)
    logical = [layout.data_index[q] for q in layout.logical_z_support]
    n_logical = len(logical)

    # ---- emission helpers, one per wire format -------------------------

    @guppy
    def _emit_packed(det: array[bool, _N], round_index: int) -> None:
        words = array(nat(0) for _ in range(comptime(n_words)))
        for ci in range(comptime(n_checks)):
            if det[ci]:
                w = ci // comptime(WORD_BITS)
                b = ci % comptime(WORD_BITS)
                words[w] = words[w] | (nat(1) << nat(b))
        output("det", words)

    @guppy
    def _emit_sparse(det: array[bool, _N], round_index: int) -> None:
        for ci in range(comptime(n_checks)):
            if det[ci]:
                output("det", round_index * comptime(n_checks) + ci)

    @guppy
    def _emit_dense(det: array[bool, _N], round_index: int) -> None:
        output("det", det)

    @guppy
    def _emit_none(det: array[bool, _N], round_index: int) -> None:
        """Control: compute detection events but emit nothing.

        Used only by the benchmark, to separate emission cost from the cost of
        simulation and RNG draws. Emits a single reduced parity so the work
        cannot be optimised away entirely.
        """
        parity = False
        for ci in range(comptime(n_checks)):
            if det[ci]:
                parity = not parity
        if round_index < 0:
            output("det", parity)

    emitter = {
        "packed": _emit_packed,
        "sparse": _emit_sparse,
        "dense": _emit_dense,
        "none": _emit_none,
    }[emit]

    # ---- the experiment ------------------------------------------------

    @guppy
    def program(
        p_data: float,
        p_read: float,
        dt: float,
        n_rounds: int,
        base_seed: int,
        inject_round: int,
        inject_qubit: int,
    ) -> None:
        # Exactly one RNG, seeded per shot; see docs/noise_spec.md section 7.
        rng = RNG(base_seed + int(get_current_shot()))

        data = array(qubit() for _ in range(comptime(n_data)))
        prev = array(False for _ in range(comptime(n_checks)))

        j = 0
        while j < n_rounds:
            # Data faults BEFORE extraction, one Bernoulli draw per qubit.
            apply_data_noise(rng, data, p_data)

            # Deterministic fault injection for tests; inject_qubit < 0 is off.
            if j == inject_round:
                if inject_qubit >= 0:
                    x(data[inject_qubit])

            syn = array(False for _ in range(comptime(n_checks)))
            for ci in range(comptime(n_checks)):
                anc = qubit()
                for k in range(comptime(max_weight)):
                    idx = comptime(supports)[ci * comptime(max_weight) + k]
                    if idx >= 0:
                        cx(data[idx], anc)
                # Ideal measurement; only the classical record is flipped.
                syn[ci] = flip_record(rng, measure(anc).read(), p_read)

            # Eq. A1. prev starts all-False, so round 0 gives D = s.
            det = array(False for _ in range(comptime(n_checks)))
            for ci in range(comptime(n_checks)):
                det[ci] = syn[ci] != prev[ci]
                prev[ci] = syn[ci]
            emitter(det, j)
            j += 1

        # Final data readout: ideal measurement, classical flips only.
        measured = measure_array(data)
        final = array(False for _ in range(comptime(n_data)))
        for i in range(comptime(n_data)):
            final[i] = flip_record(rng, measured[i].read(), p_read)

        # Final detector block: support parity of the readout vs last syndrome.
        last = array(False for _ in range(comptime(n_checks)))
        for ci in range(comptime(n_checks)):
            parity = False
            for k in range(comptime(max_weight)):
                idx = comptime(supports)[ci * comptime(max_weight) + k]
                if idx >= 0:
                    if final[idx]:
                        parity = not parity
            last[ci] = parity != prev[ci]
        emitter(last, n_rounds)

        # Logical Z: parity of the final readout over row 0.
        obs = False
        for t in range(comptime(n_logical)):
            if final[comptime(logical)[t]]:
                obs = not obs
        output("obs", obs)

        # Per-round metadata. Under a static schedule dt, p_data and p_read are
        # the same every round, so they are emitted once per shot; t_j = j * dt.
        # See the output schema in docs/interfaces.md.
        output("dt", dt)
        output("p_data", p_data)
        output("p_read", p_read)
        output("n_rounds", n_rounds)
        output("n_checks", comptime(n_checks))
        output("n_words", comptime(n_words))

        rng.discard()

    return MemoryProgram(
        layout=layout, emit=emit, program=program, n_words=n_words
    )
