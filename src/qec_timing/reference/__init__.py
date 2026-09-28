"""Independent pure-Python phenomenological reference model.

Per CLAUDE.md this package exists **only** to cross-check the Selene
implementation; it is not the product. It uses numpy, Stim and pymatching, and
imports nothing from Guppy or Selene.
"""

from .circuit import (
    build_memory_circuit,
    build_memory_circuit_from_schedule,
    data_flip_probe,
    readout_flip_probe,
    single_error_circuit,
)
from .decoder import Decoder, PyMatchingDecoder, build_decoder
from ..layout import RotatedSurfaceCodeZSector
from .noise import (
    DEFAULT_TICK,
    RoundSchedule,
    constant_schedule,
    data_probabilities,
    p_data,
    p_idle,
    p_read,
    schedule_from_intervals,
    ticks_for,
)
from .sampler import (
    LogicalErrorEstimate,
    estimate_logical_error,
    estimate_logical_error_adaptive,
)

__all__ = [
    "DEFAULT_TICK",
    "Decoder",
    "LogicalErrorEstimate",
    "PyMatchingDecoder",
    "RotatedSurfaceCodeZSector",
    "RoundSchedule",
    "build_decoder",
    "build_memory_circuit",
    "build_memory_circuit_from_schedule",
    "constant_schedule",
    "data_probabilities",
    "data_flip_probe",
    "estimate_logical_error",
    "estimate_logical_error_adaptive",
    "p_data",
    "p_idle",
    "p_read",
    "readout_flip_probe",
    "schedule_from_intervals",
    "ticks_for",
    "single_error_circuit",
]
