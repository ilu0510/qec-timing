"""Opt-in idle_ops in the Guppy circuit builder (Stage 6b variant A)."""

from __future__ import annotations

import pytest
from selene_sim import DepolarizingErrorModel

from qec_timing.circuit import build_surface_code_program

RUN = dict(p_data=0.05, p_read=0.03, dt=1.0, n_rounds=3, base_seed=7,
           inject_round=-1, inject_qubit=-1, inject_patch=0)


def _shots(program, error_model=None, shots=20):
    em = program.emulator().with_shots(shots).with_seed(11)
    if error_model is not None:
        em = em.with_error_model(error_model)
    return [s.collate_tags() for s in em.run(**RUN).results]


@pytest.mark.slow
@pytest.mark.parametrize("basis", ["Z", "X"])
def test_default_equals_idle_ops_zero_shot_for_shot(basis):
    """Same seed, same outputs: the default program is unchanged by the option."""
    default = build_surface_code_program(3, basis=basis)
    zero = build_surface_code_program(3, basis=basis, idle_ops=0)
    em = DepolarizingErrorModel(p_1q=1e-2, p_2q=2e-2, p_meas=1e-2, p_init=0.0)
    assert _shots(default) == _shots(zero)
    assert _shots(default, em) == _shots(zero, em)


@pytest.mark.slow
def test_idle_ops_are_charged_and_not_optimised_away():
    """Each phased_x(q, 0, 0) is one charged 1q op: count it in Selene's metrics."""
    d, idle_ops = 3, 2
    em = DepolarizingErrorModel(p_1q=0.0, p_2q=0.0, p_meas=0.0, p_init=0.0)
    # No injected flips. The conditional preparation corrections depend on the
    # seeded projection outcomes, identical in both programs, so they cancel.
    run = {**RUN, "p_data": 0.0, "p_read": 0.0}

    def gates_1q(program):
        res = (program.emulator().with_error_model(em).with_shots(1).with_seed(5)
               .with_metrics().run(**run))
        m = res.metrics if not callable(res.metrics) else res.metrics()
        m = m[0] if isinstance(m, list) else m
        return m["error_model"]["gates_1q"]

    base = gates_1q(build_surface_code_program(d))
    idle = gates_1q(build_surface_code_program(d, idle_ops=idle_ops))
    assert idle - base == idle_ops * d * d * RUN["n_rounds"]


@pytest.mark.parametrize("bad", [-1, 1.5, True])
def test_idle_ops_validation(bad):
    with pytest.raises(ValueError):
        build_surface_code_program(3, idle_ops=bad)
