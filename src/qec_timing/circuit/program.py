"""Full, ideal CSS extraction; phenomenological noise in one Pauli sector.

Each patch is prepared in the +1 eigenspace of every stabilizer. Z memory
stores |0_L> and uses the existing X-only model. X memory stores |+_L> and
conjugates that channel to Z faults. The complementary checks are measured
ideally, without additional record noise. No gate noise or physical gate
duration is modelled. All data stay live; one ancilla is reused sequentially.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

import numpy as np
from guppylang import guppy
from guppylang.std.array import array
from guppylang.std.builtins import comptime, nat, output
from guppylang.std.qsystem.random import RNG
from guppylang.std.qsystem.utils import get_current_shot
from guppylang.std.quantum import cx, h, measure, measure_array, qubit, x, z

from ..layout import RotatedSurfaceCode
from ..noise.primitives import apply_data_noise, flip_record
from ..wire import WORD_BITS, EmitFormat, unpack_detection_events

_N = guppy.nat_var("_SURFACE_N")
Basis = Literal["Z", "X"]


def _supports(layout, checks):
    """Fixed-width support table with -1 padding for boundary checks."""
    flat = []
    for support in checks:
        flat.extend([layout.data_index[q] for q in support] + [-1] * (4 - len(support)))
    return flat


def _preparation_corrections(layout, checks):
    """GF(2) right inverse: correction i flips only preparation syndrome i.

    Row reduction tracks the syndrome transformation, then assigns pivot
    qubits. Z corrections for X checks preserve logical Z (and vice versa).
    """
    rows = [sum(1 << layout.data_index[q] for q in s) for s in checks]
    transforms = [1 << i for i in range(len(rows))]
    pivots = []
    rank = 0
    for col in range(layout.n_data):
        pivot = next((r for r in range(rank, len(rows)) if rows[r] >> col & 1), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        transforms[rank], transforms[pivot] = transforms[pivot], transforms[rank]
        for r in range(len(rows)):
            if r != rank and rows[r] >> col & 1:
                rows[r] ^= rows[rank]
                transforms[r] ^= transforms[rank]
        pivots.append(col)
        rank += 1
    if rank != len(checks):
        raise ValueError("preparation stabilizers must be independent")
    corrections = [False] * (len(checks) * layout.n_data)
    for ci in range(len(checks)):
        for r, q in enumerate(pivots):
            corrections[ci * layout.n_data + q] = bool(transforms[r] >> ci & 1)
    return corrections


@dataclass
class SurfaceCodeProgram:
    """Specialized Guppy circuit and its Selene result/decoder adapter."""

    layout: RotatedSurfaceCode
    q_d: int
    basis: Basis
    emit: EmitFormat
    program: object
    n_words: int

    @property
    def n_qubits(self):
        return self.q_d * self.layout.n_data + 1

    def emulator(self, **kwargs):
        return self.program.emulator(n_qubits=self.n_qubits, **kwargs).stabilizer_sim()

    def run(self, *, p_data, p_read, n_rounds, base_seed, dt=1.0, shots=1,
            inject_round=-1, inject_qubit=-1, inject_patch=0, n_processes=1):
        """Fault indices are local to inject_patch; negative qubit disables it."""
        if not (0 <= p_data <= 1 and 0 <= p_read <= 1):
            raise ValueError("noise probabilities must lie in [0, 1]")
        if not isinstance(n_rounds, int) or n_rounds < 1:
            raise ValueError("n_rounds must be a positive integer")
        if not math.isfinite(dt) or dt <= 0:
            raise ValueError("dt must be finite and positive")
        if not isinstance(shots, int) or shots < 1:
            raise ValueError("shots must be a positive integer")
        if inject_qubit != -1 and not (
            0 <= inject_patch < self.q_d and 0 <= inject_qubit < self.layout.n_data
            and 0 <= inject_round < n_rounds
        ):
            raise ValueError("injection must identify a valid patch, qubit and round")
        emulator = self.emulator().with_shots(shots)
        if n_processes > 1:
            emulator = emulator.with_n_processes(n_processes)
        return emulator.run(p_data=p_data, p_read=p_read, dt=dt,
                            n_rounds=n_rounds, base_seed=base_seed,
                            inject_round=inject_round, inject_qubit=inject_qubit,
                            inject_patch=inject_patch)

    def detection_events(self, result, n_rounds):
        """Per-shot (patch, round, check); one patch preserves legacy 2D shape."""
        events = []
        for shot in result.results:
            blocks = unpack_detection_events(
                shot.collate_tags(), n_checks=self.layout.n_checks,
                n_rounds=(n_rounds + 1) * self.q_d - 1,
                n_words=self.n_words, emit=self.emit,
            )
            patches = blocks.reshape(n_rounds + 1, self.q_d, -1).transpose(1, 0, 2)
            events.append(patches[0] if self.q_d == 1 else patches)
        return events

    def observables(self, result):
        """Logical readout parities, shape (shots, patches)."""
        return np.asarray([s.collate_tags()["obs"] for s in result.results], dtype=bool)

    def decode(self, result, n_rounds, decoder):
        """Use one existing single-sector Decoder independently on each patch."""
        expected = (n_rounds + 1) * self.layout.n_checks
        if decoder.num_detectors != expected or decoder.num_observables != 1:
            raise ValueError("decoder must match one patch and this round count")
        events = np.asarray(self.detection_events(result, n_rounds)).reshape(
            -1, self.q_d, n_rounds + 1, self.layout.n_checks)
        return np.column_stack([
            decoder.decode_batch(events[:, p].reshape(-1, expected))[:, 0]
            for p in range(self.q_d)
        ])


def build_surface_code_program(d: int, q_d: int = 1, *, basis: Basis = "Z",
                               emit: EmitFormat = "packed") -> SurfaceCodeProgram:
    """Compile-time distance, patch count and memory basis; runtime noise/timing.

    det blocks are round -> patch -> check, including the final closure.
    obs is repeated in patch order. complementary_syn packs ideal opposite
    stabilizer measurements (round -> patch), without a final closure.
    """
    if not isinstance(d, int) or isinstance(d, bool):
        raise ValueError("distance must be an odd integer >= 3")
    layout = RotatedSurfaceCode(d)
    if not isinstance(q_d, int) or isinstance(q_d, bool) or q_d < 1:
        raise ValueError("q_d must be a positive integer")
    if basis not in ("Z", "X"):
        raise ValueError("basis must be Z or X")
    if emit not in ("packed", "sparse", "dense", "none"):
        raise ValueError("unknown emit format")
    x_memory = basis == "X"
    n_data, n_checks = layout.n_data, layout.n_checks
    total_data, total_checks = q_d * n_data, q_d * n_checks
    active = layout.x_checks if x_memory else layout.z_checks
    complementary = layout.z_checks if x_memory else layout.x_checks
    supports, other_supports = _supports(layout, active), _supports(layout, complementary)
    corrections = _preparation_corrections(layout, complementary)
    logical = [layout.data_index[q] for q in
               (layout.logical_x_support if x_memory else layout.logical_z_support)]
    n_words = math.ceil(n_checks / WORD_BITS)

    @guppy
    def extract(data: array[qubit, _N], patch: int, ci: int, other: bool) -> bool:
        anc = qubit()
        is_x = comptime(x_memory)
        if other:
            is_x = not is_x
        if is_x:
            h(anc)
        for k in range(4):
            idx = comptime(supports)[ci * 4 + k]
            if other:
                idx = comptime(other_supports)[ci * 4 + k]
            if idx >= 0:
                if is_x:
                    cx(anc, data[patch * comptime(n_data) + idx])
                else:
                    cx(data[patch * comptime(n_data) + idx], anc)
        if is_x:
            h(anc)
        return measure(anc).read()

    @guppy
    def emit_block(det: array[bool, _N], block: int, other: bool) -> None:
        if comptime(emit == "packed") or other:
            words = array(nat(0) for _ in range(comptime(n_words)))
            for ci in range(comptime(n_checks)):
                if det[ci]:
                    w, b = ci // comptime(WORD_BITS), ci % comptime(WORD_BITS)
                    words[w] = words[w] | (nat(1) << nat(b))
            if other:
                output("complementary_syn", words)
            else:
                output("det", words)
        elif comptime(emit == "sparse"):
            for ci in range(comptime(n_checks)):
                if det[ci]:
                    output("det", block * comptime(n_checks) + ci)
        elif comptime(emit == "dense"):
            output("det", det)

    @guppy
    def program(p_data: float, p_read: float, dt: float, n_rounds: int,
                base_seed: int, inject_round: int, inject_qubit: int,
                inject_patch: int) -> None:
        rng = RNG(base_seed + int(get_current_shot()))
        data = array(qubit() for _ in range(comptime(total_data)))
        if comptime(x_memory):
            for i in range(comptime(total_data)):
                h(data[i])

        # Project complementary stabilizers, then correct their random signs.
        # The product preparation already fixes active checks and logical parity.
        for patch in range(comptime(q_d)):
            prep = array(False for _ in range(comptime(n_checks)))
            for ci in range(comptime(n_checks)):
                prep[ci] = extract(data, patch, ci, True)
            for ci in range(comptime(n_checks)):
                if prep[ci]:
                    for i in range(comptime(n_data)):
                        if comptime(corrections)[ci * comptime(n_data) + i]:
                            if comptime(x_memory):
                                x(data[patch * comptime(n_data) + i])
                            else:
                                z(data[patch * comptime(n_data) + i])

        prev = array(False for _ in range(comptime(total_checks)))
        j = 0
        while j < n_rounds:
            if comptime(x_memory):
                for i in range(comptime(total_data)):
                    h(data[i])
            apply_data_noise(rng, data, p_data)
            if comptime(x_memory):
                for i in range(comptime(total_data)):
                    h(data[i])
            if j == inject_round and inject_qubit >= 0:
                idx = inject_patch * comptime(n_data) + inject_qubit
                if comptime(x_memory):
                    z(data[idx])
                else:
                    x(data[idx])
            for patch in range(comptime(q_d)):
                det = array(False for _ in range(comptime(n_checks)))
                other = array(False for _ in range(comptime(n_checks)))
                for ci in range(comptime(n_checks)):
                    syn = flip_record(rng, extract(data, patch, ci, False), p_read)
                    idx = patch * comptime(n_checks) + ci
                    det[ci] = syn != prev[idx]
                    prev[idx] = syn
                    other[ci] = extract(data, patch, ci, True)
                emit_block(det, j * comptime(q_d) + patch, False)
                if comptime(emit != "none"):
                    emit_block(other, j * comptime(q_d) + patch, True)
            j += 1

        if comptime(x_memory):
            for i in range(comptime(total_data)):
                h(data[i])
        measured = measure_array(data)
        final = array(False for _ in range(comptime(total_data)))
        for i in range(comptime(total_data)):
            final[i] = flip_record(rng, measured[i].read(), p_read)
        for patch in range(comptime(q_d)):
            last = array(False for _ in range(comptime(n_checks)))
            for ci in range(comptime(n_checks)):
                parity = False
                for k in range(4):
                    idx = comptime(supports)[ci * 4 + k]
                    if idx >= 0:
                        if final[patch * comptime(n_data) + idx]:
                            parity = not parity
                last[ci] = parity != prev[patch * comptime(n_checks) + ci]
            emit_block(last, n_rounds * comptime(q_d) + patch, False)
            obs = False
            for k in range(comptime(len(logical))):
                if final[patch * comptime(n_data) + comptime(logical)[k]]:
                    obs = not obs
            output("obs", obs)
        output("dt", dt)
        output("p_data", p_data)
        output("p_read", p_read)
        output("n_rounds", n_rounds)
        output("n_checks", comptime(n_checks))
        output("n_words", comptime(n_words))
        output("q_d", comptime(q_d))
        output("x_memory", comptime(x_memory))
        rng.discard()

    return SurfaceCodeProgram(layout, q_d, basis, emit, program, n_words)
