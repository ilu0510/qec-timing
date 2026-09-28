"""Decoder seam.

The production decoder is owned by teammate #3 (CLAUDE.md, "Team split"). This
module defines the contract the reference model decodes against, plus a
pymatching implementation used for cross-checking. Anything that needs a decoder
should depend on the ``Decoder`` protocol, never on pymatching directly.

The contract is recorded in ``docs/interfaces.md``.
"""

from __future__ import annotations

from typing import Protocol, runtime_checkable

import numpy as np
import pymatching
import stim
from numpy.typing import NDArray

__all__ = ["Decoder", "PyMatchingDecoder", "build_decoder"]


@runtime_checkable
class Decoder(Protocol):
    """Batch decoder over a Stim detector error model.

    Implementations are constructed from a ``stim.DetectorErrorModel``, which
    already carries the log-likelihood edge weights of Eq. A2
    (``w_e = log((1 - p_e) / p_e)``) derived from the circuit's error
    probabilities. A decoder that instead follows a *believed* noise state --
    as the adaptive controller of Appendix F requires -- is constructed from a
    DEM built with those believed probabilities rather than the true ones.
    """

    @property
    def num_detectors(self) -> int:
        """Number of detectors this decoder expects per shot."""
        ...

    @property
    def num_observables(self) -> int:
        """Number of logical observables predicted per shot."""
        ...

    def decode_batch(self, detectors: NDArray[np.bool_]) -> NDArray[np.bool_]:
        """Predict observable flips for a batch of shots.

        Args:
            detectors: shape ``(n_shots, num_detectors)``, boolean.

        Returns:
            shape ``(n_shots, num_observables)``, boolean predicted flips.
        """
        ...


class PyMatchingDecoder:
    """Minimum-weight perfect matching decoder (Appendix A, ref. [37]).

    Handles the degenerate all-zero DEM that arises at ``p = 0``: pymatching
    cannot be built from a model with no error mechanisms, so in that case the
    decoder predicts no flips, which is the correct answer when no error can
    occur.
    """

    def __init__(self, dem: stim.DetectorErrorModel) -> None:
        self._num_detectors = dem.num_detectors
        self._num_observables = dem.num_observables
        self._trivial = len(dem) == 0 or all(
            instruction.type != "error" for instruction in dem.flattened()
        )
        self._matching = (
            None if self._trivial else pymatching.Matching.from_detector_error_model(dem)
        )

    @property
    def num_detectors(self) -> int:
        return self._num_detectors

    @property
    def num_observables(self) -> int:
        return self._num_observables

    @property
    def is_trivial(self) -> bool:
        """True when the error model contains no error mechanisms."""
        return self._trivial

    def decode_batch(self, detectors: NDArray[np.bool_]) -> NDArray[np.bool_]:
        detectors = np.asarray(detectors)
        if detectors.ndim != 2 or detectors.shape[1] != self._num_detectors:
            raise ValueError(
                f"expected detectors of shape (n_shots, {self._num_detectors}), "
                f"got {detectors.shape}"
            )
        if self._trivial:
            return np.zeros(
                (detectors.shape[0], self._num_observables), dtype=np.bool_
            )
        return self._matching.decode_batch(detectors).astype(np.bool_)


def build_decoder(
    circuit: stim.Circuit, *, decompose_errors: bool = False
) -> PyMatchingDecoder:
    """Build the default decoder for a circuit.

    ``decompose_errors=False`` is correct here: every X fault flips at most two
    Z-checks, so the detector error model is already graphlike and needs no
    decomposition into graphlike parts.
    """
    return PyMatchingDecoder(
        circuit.detector_error_model(decompose_errors=decompose_errors)
    )
