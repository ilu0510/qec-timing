"""Selectable CNOT order within each stabilizer check.

Shared by the Guppy circuit (``program.py``) and the Stage 6 Stim twin
(``qec_timing.circuit_level``), so the two cannot disagree about which data
qubit each check touches when.

Orders
------
``"row_major"`` (default): each check's support exactly as the layout gives it,
unchanged. For this layout that is NW, NE, SW, SE for both sectors.

``"hook_safe"``: the standard rotated-surface-code shapes. Checks are extracted
one at a time, so the only damaging hook is an ancilla fault after the second of
four CNOTs, which leaves an error on the LAST TWO data qubits of the check:

* X-check (ancilla is the CNOT control) -> X hook. Logical X of this layout is a
  column (``logical_x_support``), so the hook pair must be horizontal:
  NW, NE, SW, SE ("Z" shape).
* Z-check (ancilla is the CNOT target) -> Z hook. Logical Z is a row
  (``logical_z_support``), so the hook pair must be vertical:
  NW, SW, NE, SE ("N" shape).

The shapes are derived from the layout's logical supports, not hard-coded, and
are independent of the memory basis: each check type only damages one logical.

Boundary checks: a weight-2 check is a 2x2 plaquette with one side missing. Its
two qubits keep the order of the two time slots they occupy in the chosen
shape. Sorting a pair by the shape's key gives exactly that order, for every
missing side, so weight-2 and weight-4 checks share one rule. (A weight-2
check's hook is weight 1, so its order cannot reduce distance either way.) The
order of present qubits is what the sequential circuit executes; the ``-1``
padding added by ``program._supports`` sits at the end, as before.
"""

from __future__ import annotations

from typing import Literal

from ..layout import RotatedSurfaceCode

__all__ = ["CheckType", "ORDERS", "Order", "order_supports"]

Order = Literal["row_major", "hook_safe"]
CheckType = Literal["X", "Z"]
ORDERS: tuple[str, ...] = ("row_major", "hook_safe")

Coord = tuple[int, int]


def _row_major(q: Coord) -> tuple[int, int]:
    return q


def _column_major(q: Coord) -> tuple[int, int]:
    return (q[1], q[0])


def _is_row(support) -> bool:
    return len({r for r, _ in support}) == 1


def _is_column(support) -> bool:
    return len({c for _, c in support}) == 1


def order_supports(layout: RotatedSurfaceCode, checks, check_type: CheckType,
                   order: Order = "row_major"):
    """Return ``checks`` with each support in the requested CNOT order.

    ``"row_major"`` returns ``checks`` itself, untouched, so the default path is
    identical to calling ``program._supports`` on the layout's checks directly.
    """
    if order not in ORDERS:
        raise ValueError(f"order must be one of {ORDERS}, got {order!r}")
    if check_type not in ("X", "Z"):
        raise ValueError(f"check_type must be 'X' or 'Z', got {check_type!r}")
    if order == "row_major":
        return checks

    # The hook of a check type must be perpendicular to the logical it can damage.
    damaged = layout.logical_x_support if check_type == "X" else layout.logical_z_support
    if _is_column(damaged):
        key = _row_major       # last two share a row: horizontal hook
    elif _is_row(damaged):
        key = _column_major    # last two share a column: vertical hook
    else:
        raise ValueError("logical operator is neither a row nor a column")
    return tuple(tuple(sorted(support, key=key)) for support in checks)
