"""PLACEHOLDER syndrome extraction -- stands in for teammate #1's circuit.

Not the production circuit. `program.py` implements the simplest extraction that
is provably equivalent to the reference model's ideal MPP: per Z-check, an
ancilla in |0>, a CNOT from each data qubit of the support, and an ideal
measurement. All gates are ideal; every fault comes from `qec_timing.noise`.

Replace with the real circuit when it lands. The check geometry comes from the
shared `qec_timing.layout`, so this and the reference model cannot drift apart.
"""

from .program import (
    WORD_BITS,
    EmitFormat,
    MemoryProgram,
    build_memory_program,
)
from .reader import unpack_detection_events, unpack_packed_words

__all__ = [
    "WORD_BITS",
    "EmitFormat",
    "MemoryProgram",
    "build_memory_program",
    "unpack_detection_events",
    "unpack_packed_words",
]
