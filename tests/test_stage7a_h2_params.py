"""Stage 7a: the stdlib-only H2 helpers, checked from the MAIN venv.

``qec_timing.h2`` also holds Guppy-0.x code that the main venv must never
import, so these modules are loaded by file path, not as a package.
"""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import pytest

from qec_timing.circuit.ordering import order_supports
from qec_timing.circuit_level import DephasingH2
from qec_timing.layout import RotatedSurfaceCode

H2_DIR = Path(__file__).resolve().parents[1] / "src" / "qec_timing" / "h2"


def _load(name):
    spec = importlib.util.spec_from_file_location(f"h2_{name}", H2_DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


geometry = _load("geometry")
noise_params = _load("noise_params")


@pytest.mark.parametrize("d", [3, 5, 7, 9, 11])
def test_h2_geometry_matches_hook_safe_ordering(d):
    """The re-implementation equals layout + ordering.order_supports('hook_safe')."""
    L = RotatedSurfaceCode(d)
    z = [[L.data_index[q] for q in s] for s in order_supports(L, L.z_checks, "Z", "hook_safe")]
    x = [[L.data_index[q] for q in s] for s in order_supports(L, L.x_checks, "X", "hook_safe")]
    assert geometry.z_checks_hook_safe(d) == z
    assert geometry.x_checks_hook_safe(d) == x


def test_idle_pz_conventions():
    n = noise_params.H2_1E
    d = 0.1  # seconds
    lin = n.linear_dephasing_rate * d
    for conv, omega in (("cycles", 2 * math.pi * 0.043 * 2.0), ("radians", 0.043 * 2.0)):
        q = math.sin(0.5 * omega * d) ** 2
        assert noise_params.idle_pz(d, convention=conv) == pytest.approx(q + lin - 2 * q * lin)
    assert noise_params.idle_pz(0.0) == 0.0
    with pytest.raises(ValueError):
        noise_params.idle_pz(1.0, convention="hertz")


@pytest.mark.parametrize("conv", ["cycles", "radians"])
@pytest.mark.parametrize("t_round", [1e-3, 2e-3, 5e-3])
def test_twin_dephasing_reproduces_h2_pz_per_period(conv, t_round):
    """Same quadratic and linear terms; they differ only in how the two combine.

    The twin adds them (q + l); noise_params composes two independent Z
    channels (q + l - 2 q l). So twin - H2 = 2 q l exactly.
    """
    f_quad, f_lin = noise_params.twin_dephasing(t_round, convention=conv)
    twin = DephasingH2(f_quad, f_lin)
    n = noise_params.H2_1E
    omega = 2 * math.pi * 0.086 if conv == "cycles" else 0.086
    for dt_rounds in (1, 10, 44, 180):
        secs = dt_rounds * t_round
        q = math.sin(0.5 * omega * secs) ** 2
        lin = n.linear_dephasing_rate * secs
        assert twin.probability(dt_rounds) == pytest.approx(q + lin, rel=1e-12)
        h2 = noise_params.idle_pz(secs, convention=conv)
        assert twin.probability(dt_rounds) - h2 == pytest.approx(2 * q * lin, rel=1e-9, abs=1e-15)
