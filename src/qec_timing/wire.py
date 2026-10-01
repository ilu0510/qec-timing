"""Decode the detection-event wire formats back into a dense boolean array.

The packed format is the one specified in ``docs/interfaces.md``; the sparse and
dense readers exist so the three can be compared and benchmarked.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from typing import Literal

WORD_BITS = 64
EmitFormat = Literal["packed", "sparse", "dense", "none"]

__all__ = ["unpack_detection_events", "unpack_packed_words"]


def unpack_packed_words(
    words: NDArray[np.uint64] | list, n_checks: int
) -> NDArray[np.bool_]:
    """Unpack bit-packed words into a ``(n_blocks, n_checks)`` boolean array.

    Bit ``b`` of word ``w`` is check index ``w * 64 + b``, LSB first; words run
    in ascending order; the final word is zero-padded in its high bits. This is
    the reference implementation of the snippet in ``docs/interfaces.md``.
    """
    words = np.asarray(words, dtype=np.uint64)
    if words.ndim == 1:
        words = words[None, :]
    shifts = np.arange(WORD_BITS, dtype=np.uint64)
    bits = ((words[:, :, None] >> shifts) & np.uint64(1)).astype(np.bool_)
    return bits.reshape(words.shape[0], -1)[:, :n_checks]


def unpack_detection_events(
    shot_tags: dict[str, list],
    *,
    n_checks: int,
    n_rounds: int,
    n_words: int,
    emit: EmitFormat,
) -> NDArray[np.bool_]:
    """Rebuild the ``(n_rounds + 1, n_checks)`` detection-event array.

    Args:
        shot_tags: one shot's ``QsysShot.collate_tags()``.
        n_checks, n_rounds, n_words: geometry, as emitted by the program.
        emit: which wire format the program used.

    ``n_rounds + 1`` blocks: one per round plus the final data-readout block.
    """
    n_blocks = n_rounds + 1
    # A shot in which nothing ever fires emits no "det" entries at all under the
    # sparse format, so a missing tag means "all zero", not an error.
    raw = shot_tags.get("det", [])

    if emit == "sparse":
        events = np.zeros((n_blocks, n_checks), dtype=np.bool_)
        for index in raw:
            block, check = divmod(int(index), n_checks)
            events[block, check] = True
        return events

    if emit == "dense":
        events = np.asarray(raw, dtype=np.bool_).reshape(n_blocks, n_checks)
        return events

    if emit == "packed":
        words = np.asarray(raw, dtype=np.uint64).reshape(n_blocks, n_words)
        return unpack_packed_words(words, n_checks)

    raise ValueError(f"unknown emit format {emit!r}")
