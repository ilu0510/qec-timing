"""Rotated-surface-code Z-checks in hook_safe order. Stdlib only.

Re-implements, without importing, the logic of ``qec_timing.layout``
(``RotatedSurfaceCodeZSector.z_checks``, ``data_index``) and
``qec_timing.circuit.ordering`` (``order_supports(..., "Z", "hook_safe")``).
A main-venv test loads this file by path and checks it gives exactly the
same supports in the same order.

Layout: data (r, c) on a d x d grid, index r * d + c. Z-checks: weight-4
plaquettes at (r + c) even, then weight-2 checks on the left column (odd rows)
and the right column (even rows). Logical Z is row 0, a horizontal line, so a
Z-check's hook (an error on its last two data qubits) must be vertical:
hook_safe orders each support column-major, NW, SW, NE, SE ("N" shape). A
weight-2 check sorted the same way keeps the time-slot order its present
qubits occupy.
"""

from __future__ import annotations

__all__ = ["data_index", "x_checks_hook_safe", "z_checks_hook_safe"]


def data_index(d: int, r: int, c: int) -> int:
    return r * d + c


def z_checks_hook_safe(d: int) -> list[list[int]]:
    """Data-qubit indices of every Z-check, in CNOT order (hook_safe)."""
    if d < 3 or d % 2 == 0:
        raise ValueError("distance must be an odd integer >= 3")
    checks: list[list[tuple[int, int]]] = []
    for r in range(d - 1):
        for c in range(d - 1):
            if (r + c) % 2 == 0:
                checks.append([(r, c), (r, c + 1), (r + 1, c), (r + 1, c + 1)])
    for r in range(1, d - 1, 2):
        checks.append([(r, 0), (r + 1, 0)])
    for r in range(0, d - 1, 2):
        checks.append([(r, d - 1), (r + 1, d - 1)])
    return [[data_index(d, r, c) for r, c in sorted(s, key=lambda q: (q[1], q[0]))]
            for s in checks]


def x_checks_hook_safe(d: int) -> list[list[int]]:
    """Data-qubit indices of every X-check, in CNOT order (hook_safe).

    Mirrors ``RotatedSurfaceCode.x_checks``: weight-4 plaquettes at (r + c) odd,
    then the top row (even c) and the bottom row (odd c). Logical X is column
    0, a vertical line, so an X-check's hook must be horizontal: hook_safe
    orders each support row-major, NW, NE, SW, SE ("Z" shape).
    """
    if d < 3 or d % 2 == 0:
        raise ValueError("distance must be an odd integer >= 3")
    checks: list[list[tuple[int, int]]] = []
    for r in range(d - 1):
        for c in range(d - 1):
            if (r + c) % 2 == 1:
                checks.append([(r, c), (r, c + 1), (r + 1, c), (r + 1, c + 1)])
    for c in range(0, d - 1, 2):
        checks.append([(0, c), (0, c + 1)])
    for c in range(1, d - 1, 2):
        checks.append([(d - 1, c), (d - 1, c + 1)])
    return [[data_index(d, r, c) for r, c in sorted(s)] for s in checks]
