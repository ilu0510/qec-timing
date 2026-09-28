"""Sampling and logical-error estimation for the reference model.

Sampling is batched so that memory stays bounded: a sweep point near the optimal
interval can need several hundred thousand shots, and the small-interval points
have >10^5 detectors per shot, so materialising every shot at once is not an
option.

Shot counts may be chosen adaptively by targeting a fixed number of logical
failures, which keeps the *relative* error bar roughly constant as ``p_L``
varies over orders of magnitude across a ``dt`` sweep.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import stim
from scipy.stats import beta as _beta

from .decoder import Decoder

__all__ = [
    "LogicalErrorEstimate",
    "estimate_logical_error",
    "estimate_logical_error_adaptive",
]

# Cap on bytes held for one batch of detector data (Stim returns one byte per
# detector per shot).
_BATCH_BYTE_BUDGET = 64 * 1024 * 1024


@dataclass(frozen=True)
class LogicalErrorEstimate:
    """Logical error probability with an exact binomial confidence interval."""

    failures: int
    shots: int
    confidence: float
    ci_low: float
    ci_high: float

    @property
    def p_L(self) -> float:
        return self.failures / self.shots if self.shots else float("nan")

    @property
    def stderr(self) -> float:
        """Normal-approximation standard error, for plotting error bars."""
        if not self.shots:
            return float("nan")
        p = self.p_L
        return float(np.sqrt(max(p * (1.0 - p), 0.0) / self.shots))

    @property
    def relative_stderr(self) -> float:
        return self.stderr / self.p_L if self.failures else float("inf")

    def __str__(self) -> str:
        return (
            f"p_L = {self.p_L:.6g} "
            f"[{self.ci_low:.6g}, {self.ci_high:.6g}] "
            f"({self.failures}/{self.shots} failures)"
        )


def clopper_pearson(
    failures: int, shots: int, confidence: float = 0.999
) -> tuple[float, float]:
    """Exact (Clopper-Pearson) binomial interval.

    Preferred over a normal approximation because the logical error rates here
    reach 1e-4 and below, where the normal approximation is not trustworthy.
    """
    if shots <= 0:
        raise ValueError("shots must be positive")
    if not 0 <= failures <= shots:
        raise ValueError(f"failures={failures} outside [0, {shots}]")
    tail = (1.0 - confidence) / 2.0
    low = 0.0 if failures == 0 else float(_beta.ppf(tail, failures, shots - failures + 1))
    high = (
        1.0
        if failures == shots
        else float(_beta.ppf(1.0 - tail, failures + 1, shots - failures))
    )
    return low, high


def _batch_size(num_detectors: int, requested: int | None) -> int:
    if requested is not None:
        return max(1, requested)
    per_shot = max(num_detectors, 1)
    return max(1, min(20_000, _BATCH_BYTE_BUDGET // per_shot))


def _count_failures(
    sampler: stim.CompiledDetectorSampler, decoder: Decoder, n: int
) -> int:
    detectors, observables = sampler.sample(shots=n, separate_observables=True)
    predictions = decoder.decode_batch(detectors)
    # A shot fails when the correction differs from the true error by a
    # nontrivial logical operator (Appendix A).
    return int(np.count_nonzero(np.any(predictions != observables, axis=1)))


def estimate_logical_error(
    circuit: stim.Circuit,
    decoder: Decoder,
    *,
    shots: int,
    seed: int,
    confidence: float = 0.999,
    batch_size: int | None = None,
) -> LogicalErrorEstimate:
    """Estimate the logical error probability over a fixed number of shots."""
    if shots <= 0:
        raise ValueError("shots must be positive")
    sampler = circuit.compile_detector_sampler(seed=seed)
    chunk = _batch_size(circuit.num_detectors, batch_size)

    failures = 0
    remaining = shots
    while remaining > 0:
        n = min(chunk, remaining)
        failures += _count_failures(sampler, decoder, n)
        remaining -= n

    low, high = clopper_pearson(failures, shots, confidence)
    return LogicalErrorEstimate(failures, shots, confidence, low, high)


def estimate_logical_error_adaptive(
    circuit: stim.Circuit,
    decoder: Decoder,
    *,
    seed: int,
    target_failures: int = 200,
    min_shots: int = 2_000,
    max_shots: int = 400_000,
    confidence: float = 0.999,
    batch_size: int | None = None,
) -> LogicalErrorEstimate:
    """Sample until ``target_failures`` logical failures accumulate.

    Stops early at ``max_shots``. Targeting a failure count rather than a shot
    count keeps the relative error bar near ``1/sqrt(target_failures)``
    regardless of how small ``p_L`` is; 200 failures corresponds to roughly 7%.
    """
    if min_shots <= 0 or max_shots < min_shots:
        raise ValueError("require 0 < min_shots <= max_shots")
    sampler = circuit.compile_detector_sampler(seed=seed)
    chunk = _batch_size(circuit.num_detectors, batch_size)

    failures = 0
    shots = 0
    while shots < max_shots:
        if shots >= min_shots and failures >= target_failures:
            break
        n = min(chunk, max_shots - shots)
        failures += _count_failures(sampler, decoder, n)
        shots += n

    low, high = clopper_pearson(failures, shots, confidence)
    return LogicalErrorEstimate(failures, shots, confidence, low, high)
