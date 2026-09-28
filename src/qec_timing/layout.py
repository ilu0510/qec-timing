"""Rotated surface-code geometry, Z-check sector only.

Shared by the reference model and the Guppy path so the two cannot drift apart
(the same arrangement as ``ansatz.py``). Geometry is combinatorial, not physics:
unlike Eqs. 1-3, which are deliberately implemented twice, there is nothing here
for an independent implementation to cross-check.

Appendix A of the paper: the code is CSS, so the two Pauli sectors separate in
the phenomenological model. We simulate only the Z-check sector, which detects
X-type data faults; the X-check sector is equivalent after exchanging X and Z.

Data qubits sit on a ``d x d`` grid indexed by ``(row, col)``. Z-checks are the
weight-4 plaquettes of one checkerboard sublattice plus weight-2 checks on the
left and right columns. Verified properties (see tests):

* ``(d**2 - 1) // 2`` Z-checks, as required for one sector of a distance-d code.
* Every data qubit lies in **at most two** Z-checks, so a single X fault flips
  one or two detectors and the detector error model is graphlike -- which is
  what lets pymatching decode it directly.
* The code distance is exactly ``d``.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property

__all__ = ["RotatedSurfaceCodeZSector"]

Coord = tuple[int, int]


@dataclass(frozen=True)
class RotatedSurfaceCodeZSector:
    """Z-check sector of the rotated surface code of odd distance ``d``."""

    d: int

    def __post_init__(self) -> None:
        if self.d < 3 or self.d % 2 == 0:
            raise ValueError(f"distance must be an odd integer >= 3, got {self.d}")

    # -- data qubits ------------------------------------------------------

    @cached_property
    def data_qubits(self) -> tuple[Coord, ...]:
        """Data qubits in row-major order."""
        return tuple((r, c) for r in range(self.d) for c in range(self.d))

    @cached_property
    def data_index(self) -> dict[Coord, int]:
        """Map ``(row, col)`` to the Stim qubit index."""
        return {q: i for i, q in enumerate(self.data_qubits)}

    @property
    def n_data(self) -> int:
        return self.d * self.d

    # -- Z-checks ---------------------------------------------------------

    @cached_property
    def z_checks(self) -> tuple[tuple[Coord, ...], ...]:
        """Z-check supports, each a sorted tuple of data coordinates.

        Interior weight-4 plaquettes are taken on the ``(r + c) even``
        sublattice; the remaining ``d - 1`` checks are weight-2, split evenly
        between the left column (odd rows) and the right column (even rows).
        That parity choice is what makes every data qubit lie in at most two
        Z-checks.
        """
        d = self.d
        checks: list[tuple[Coord, ...]] = []
        for r in range(d - 1):
            for c in range(d - 1):
                if (r + c) % 2 == 0:
                    checks.append(
                        tuple(sorted({(r, c), (r, c + 1), (r + 1, c), (r + 1, c + 1)}))
                    )
        for r in range(1, d - 1, 2):
            checks.append(tuple(sorted({(r, 0), (r + 1, 0)})))
        for r in range(0, d - 1, 2):
            checks.append(tuple(sorted({(r, d - 1), (r + 1, d - 1)})))
        return tuple(checks)

    @property
    def n_checks(self) -> int:
        return (self.d * self.d - 1) // 2

    def check_coords(self, index: int) -> tuple[float, float]:
        """Centre of a Z-check, for DETECTOR coordinate annotations."""
        support = self.z_checks[index]
        return (
            sum(r for r, _ in support) / len(support),
            sum(c for _, c in support) / len(support),
        )

    # -- logical operator -------------------------------------------------

    @cached_property
    def logical_z_support(self) -> tuple[Coord, ...]:
        """Data qubits whose final-readout parity gives the logical Z.

        Row 0. The logical X operator of this layout is a *column* of X (a
        column commutes with every Z-check, a row does not), and a column meets
        row 0 in exactly one qubit, so a logical X fault flips this parity.
        """
        return tuple((0, c) for c in range(self.d))
