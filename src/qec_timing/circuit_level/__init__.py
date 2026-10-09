"""Stage 6a: circuit-level Stim twin with a switchable idle channel.

A NEW model, beyond the Stage-1 phenomenological physics: gate, reset and
readout noise plus one of two idle-noise shapes. It does not modify or replace
``qec_timing.reference`` / ``qec_timing.noise``; it only reuses the shared
layout, schedule, the Guppy program's CNOT order and the ``Decoder`` seam.
"""

from .memory import (
    CircuitLevelMemory,
    build_circuit_level_memory,
    build_constant_interval_memory,
    idle_probe_circuit,
    rate_per_unit_time,
    readout_probe_circuit,
    sample_measurements,
)
from .noise import (
    BitFlip,
    CircuitNoise,
    DephasingH2,
    DepolarizingKOps,
    ReadoutNoise,
    compose_depolarize1,
    dephasing_h2_probability,
)

__all__ = [
    "BitFlip",
    "CircuitLevelMemory",
    "CircuitNoise",
    "DephasingH2",
    "DepolarizingKOps",
    "ReadoutNoise",
    "build_circuit_level_memory",
    "build_constant_interval_memory",
    "compose_depolarize1",
    "dephasing_h2_probability",
    "idle_probe_circuit",
    "rate_per_unit_time",
    "readout_probe_circuit",
    "sample_measurements",
]
