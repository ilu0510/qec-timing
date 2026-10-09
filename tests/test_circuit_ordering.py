"""Selectable CNOT order (circuit/ordering.py): default unchanged, hook_safe shapes."""

from __future__ import annotations

import hashlib
import json

import pytest

from qec_timing.circuit import build_surface_code_program
from qec_timing.circuit.ordering import order_supports
from qec_timing.circuit.program import _supports
from qec_timing.layout import RotatedSurfaceCode

# sha256 of json([_supports(z_checks), _supports(x_checks)]), recorded from the
# code BEFORE ordering.py existed (commit 04db7a8). Guards the teammate's
# circuit and Stage 5 against any change of the default table.
PRE_ORDERING_TABLES = {
    3: "db59aa2f7bb8a56c36c5dc5bc0eede263d2029d248745524c367f24e5e69115e",
    5: "0b3993d8a7d205e0ca991cda777bb9d9b970dd8becb4788dea3950656337817c",
    7: "e08253cc994b28f4535fcedee7a1e13f8902b3f693b200bf4f72af03642fe03e",
    9: "540e291fea3084050baf73bab6e8d2e7fd47a9e08bc50b8bdeff556ce875fdf0",
    27: "2a2d6321ba7c7c728af73d5d392ebd9c84b660e4d166d0b825d83bc409d244aa",
}


@pytest.mark.parametrize("d", sorted(PRE_ORDERING_TABLES))
def test_default_order_is_byte_identical_to_pre_ordering_tables(d):
    layout = RotatedSurfaceCode(d)
    z = order_supports(layout, layout.z_checks, "Z")
    x = order_supports(layout, layout.x_checks, "X", "row_major")
    assert z is layout.z_checks and x is layout.x_checks  # passed through untouched
    tables = [_supports(layout, z), _supports(layout, x)]
    digest = hashlib.sha256(json.dumps(tables).encode()).hexdigest()
    assert digest == PRE_ORDERING_TABLES[d]


@pytest.mark.parametrize("d", [3, 5, 7, 9])
def test_hook_safe_shapes(d):
    """Last two CNOTs: same column for Z-checks, same row for X-checks."""
    layout = RotatedSurfaceCode(d)
    for checks, kind in ((layout.z_checks, "Z"), (layout.x_checks, "X")):
        ordered = order_supports(layout, checks, kind, "hook_safe")
        assert [set(s) for s in ordered] == [set(s) for s in checks]  # same checks, same order
        for s in ordered:
            if len(s) == 4:
                (r0, c0), (r1, c1), (r2, c2), (r3, c3) = s
                assert (r0, c0) == (min(r for r, _ in s), min(c for _, c in s))  # NW first
                if kind == "Z":
                    assert c2 == c3 and r2 != r3      # vertical hook: N shape
                else:
                    assert r2 == r3 and c2 != c3      # horizontal hook: Z shape
            else:
                # Boundary pairs keep their slot order, which matches the layout's.
                assert s == tuple(sorted(s))


def test_rejects_unknown_order():
    layout = RotatedSurfaceCode(3)
    with pytest.raises(ValueError):
        order_supports(layout, layout.z_checks, "Z", "spiral")
    with pytest.raises(ValueError):
        order_supports(layout, layout.z_checks, "Y", "hook_safe")
    with pytest.raises(ValueError):
        build_surface_code_program(3, order="spiral")
