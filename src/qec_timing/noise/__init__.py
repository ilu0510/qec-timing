"""Noise model for the Guppy/Selene path: Eqs. 1-3 and the Guppy primitives.

``params`` runs in Python (``exp`` is unavailable inside Guppy); ``primitives``
holds the Guppy-side noise injection. The contract both implement is written up
in ``docs/noise_spec.md``.
"""

from .params import (
    NoiseParams,
    data_probability,
    idle_probability,
    read_probability,
    schedule_data_probabilities,
)
from .primitives import apply_data_noise, bernoulli, flip_record

__all__ = [
    "NoiseParams",
    "apply_data_noise",
    "bernoulli",
    "data_probability",
    "flip_record",
    "idle_probability",
    "read_probability",
    "schedule_data_probabilities",
]
