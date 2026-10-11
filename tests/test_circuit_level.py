"""Tests for qec_timing.circuit_level.selene_mirror: the Stim copy of Selene's circuit.

The copy is only useful if it is the same circuit with the same noise. Checked here: the
gate conversion, that idle ops survive compilation, that the noiseless copy is deterministic,
and that its detector statistics match Selene's within a fixed-seed chi-square bound.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
import stim

from qec_timing.circuit_level import selene_mirror as sm

# Selene's native decomposition of cx(q0, q1), read from CircuitExtractor (selene-sim 0.3.0,
# selene-hugr-qis-compiler 0.4.2).
CX_NATIVE = [
    {"op": "Rxy", "qubit": 1, "theta": -math.pi / 2, "phi": math.pi / 2},
    {"op": "Rzz", "qubit0": 0, "qubit1": 1, "theta": math.pi / 2},
    {"op": "Rz", "qubit": 0, "theta": -math.pi / 2},
    {"op": "Rxy", "qubit": 1, "theta": math.pi / 2, "phi": math.pi},
    {"op": "Rz", "qubit": 1, "theta": -math.pi / 2},
]


def test_native_cx_converts_to_cx():
    c = stim.Circuit()
    for op in CX_NATIVE:
        sm._append_remapped(c, *sm._clifford(op))
    assert c.to_tableau() == stim.Tableau.from_named_gate("CX")


def test_non_clifford_angle_is_refused():
    with pytest.raises(ValueError):
        sm._clifford({"op": "Rz", "qubit": 0, "theta": 0.3})


@pytest.fixture(scope="module")
def geometry():
    from qec_timing.circuit_level.surface_memory import Geometry
    return Geometry("rotated", 3)


@pytest.mark.slow
def test_idle_ops_survive_compilation(geometry):
    """tau phased_x(0, 0) per data qubit per round must reach the error model as Rxy(0, 0)."""
    G, tau, rounds = geometry, 5, 2
    ops = sm.native_ops(G._runner(tau, rounds), G.n_data + G.n_anc)
    idle = [o for o in ops if o["op"] == "Rxy" and o["theta"] == 0.0 and o["phi"] == 0.0]
    assert len(idle) == tau * rounds * G.n_data


@pytest.mark.slow
def test_noiseless_copy_is_deterministic(geometry):
    G, rounds = geometry, 3
    c = sm.to_stim(sm.native_ops(G._runner(1, rounds), G.n_data + G.n_anc),
                   sm.MemoryLayout.from_geometry(G, rounds))
    det, obs = c.compile_detector_sampler(seed=1).sample(500, separate_observables=True)
    assert c.num_detectors == (rounds + 1) * len(G.zcols)
    assert det.sum() == 0 and obs.sum() == 0


@pytest.mark.slow
def test_copy_matches_selene_detector_rates(geometry):
    """Per-detector firing rates, Selene vs the copy, at raised noise and fixed seeds."""
    G, tau, rounds, n, p1, pm = geometry, 8, 12, 20_000, 0.002, 0.004
    c = sm.to_stim(sm.native_ops(G._runner(tau, rounds), G.n_data + G.n_anc),
                   sm.MemoryLayout.from_geometry(G, rounds), p_1q=p1, p_meas=pm)
    syn, dat = G.run(tau, rounds, n, seed=11, p_1q=p1, p_meas=pm)
    a = G.detectors(syn, dat).astype(bool).mean(0)
    b = c.compile_detector_sampler(seed=21).sample(n, separate_observables=True)[0].mean(0)
    z = (a - b) / np.sqrt((a * (1 - a) + b * (1 - b)) / n)
    dof = len(z)
    from scipy.stats import chi2
    assert float(np.sum(z ** 2)) < chi2.ppf(0.999, dof), f"chi2 {np.sum(z**2):.1f} over {dof}"


# ---- X-basis memory and native_memory --------------------------------------------------------

@pytest.fixture(scope="module")
def geometry_x():
    from qec_timing.circuit_level.surface_memory import Geometry
    return Geometry("rotated", 3, basis="X")


@pytest.mark.slow
def test_x_basis_noiseless_is_deterministic(geometry_x):
    """|+_L> decoded on X-checks: no events and no logical flip, in Selene and in the copy."""
    G, tau, rounds = geometry_x, 2, 3
    syn, dat = G.run(tau, rounds, 200, seed=1)
    assert G.detectors(syn, dat).sum() == 0 and G.observable(dat).sum() == 0
    c = sm.to_stim(sm.native_ops(G._runner(tau, rounds), G.n_data + G.n_anc),
                   sm.MemoryLayout.from_geometry(G, rounds))
    det, obs = c.compile_detector_sampler(seed=1).sample(500, separate_observables=True)
    assert c.num_detectors == (rounds + 1) * len(G.xcols)
    assert det.sum() == 0 and obs.sum() == 0


@pytest.mark.slow
def test_x_basis_copy_matches_selene_detector_rates(geometry_x):
    G, tau, rounds, n, p1, pm = geometry_x, 4, 8, 20_000, 0.002, 0.004
    c = sm.to_stim(sm.native_ops(G._runner(tau, rounds), G.n_data + G.n_anc),
                   sm.MemoryLayout.from_geometry(G, rounds), p_1q=p1, p_meas=pm)
    syn, dat = G.run(tau, rounds, n, seed=11, p_1q=p1, p_meas=pm)
    a = G.detectors(syn, dat).astype(bool).mean(0)
    b = c.compile_detector_sampler(seed=21).sample(n, separate_observables=True)[0].mean(0)
    z = (a - b) / np.sqrt((a * (1 - a) + b * (1 - b)) / n)
    from scipy.stats import chi2
    assert float(np.sum(z ** 2)) < chi2.ppf(0.999, len(z)), f"chi2 {np.sum(z**2):.1f} over {len(z)}"


def test_native_ops_pack_round_trip():
    from qec_timing.circuit_level.native_memory import _pack, _unpack
    ops = [{"op": "Reset", "qubit": 0}, *CX_NATIVE, {"op": "MeasureRequest", "qubit": 1},
           {"op": "Rxy", "qubit": 2, "theta": 0.0, "phi": 0.0}]
    assert _unpack(_pack(ops)) == ops


@pytest.mark.slow
def test_native_idle_only_matches_textbook():
    """Gates noiseless, idle on: native and textbook circuits are the same noisy circuit."""
    from qec_timing.circuit_level.native_memory import NativeMemory
    from qec_timing.circuit_level.stim_memory import StimMemory
    kw = dict(basis="X", dephase_linear=0.01, dephase_quad=0.05)
    n, tau, rounds = 100_000, 4.0, 6
    a = StimMemory("rotated", 3, **kw).logical_error(tau, rounds, n, seed=1)[0]
    b = NativeMemory("rotated", 3, **kw).logical_error(tau, rounds, n, seed=2)[0]
    assert abs(a - b) < 4 * np.sqrt((a * (1 - a) + b * (1 - b)) / n), (a, b)


# ---- timed_memory: uniform idle noise and measured phenomenological parameters ---------------

def test_timed_memory_noiseless_and_uniform():
    """No noise: no events. With idle noise: 7 layer channels + 1 wait per round, every qubit."""
    from qec_timing.circuit_level.timed_memory import LAYERS, TimedMemory
    c0 = TimedMemory("rotated", 3).circuit(3.0, 4)
    det, obs = c0.compile_detector_sampler(seed=1).sample(500, separate_observables=True)
    assert det.sum() == 0 and obs.sum() == 0
    flat = list(TimedMemory("rotated", 3, dephase_linear=0.01).circuit(3.0, 4).flattened())
    ch = [i for i in flat if i.name == "PAULI_CHANNEL_1"]
    n_q = len({t.value for i in flat for t in i.targets_copy() if t.is_qubit_target})
    assert len(ch) == 4 * (LAYERS + 1)
    assert all(len(i.targets_copy()) == n_q for i in ch)


def test_edge_estimator_recovers_phenomenological_parameters():
    from qec_timing.circuit_level.timed_memory import edge_probabilities, phenomenological_circuit
    e = edge_probabilities(phenomenological_circuit("rotated", 5, 12, 0.01, 0.02), 100_000)
    assert abs(e["p_data"] - 0.01) < 5e-4 and abs(e["p_read"] - 0.02) < 5e-4
    assert abs(e["p_diag"]) < 3e-4


@pytest.mark.slow
def test_native_schedule_preserves_the_circuit():
    """ASAP layering moves only ops on disjoint qubits: the error model must be identical."""
    from qec_timing.circuit_level.native_memory import NativeMemory, TimedNativeMemory
    g = dict(p1=1e-3, p2=5e-3, p_meas=4e-3, p_init=2e-4)
    a = TimedNativeMemory("rotated", 3, basis="X", **g).circuit(1.0, 4).detector_error_model()
    b = NativeMemory("rotated", 3, basis="X", **g).circuit(1.0, 4).detector_error_model()
    assert a.approx_equals(b, atol=1e-12)
    flat = list(TimedNativeMemory("rotated", 3, basis="X", dephase_linear=0.01)
                .circuit(2.0, 4).flattened())
    ch = [i for i in flat if i.name == "PAULI_CHANNEL_1"]
    n_q = len({t.value for i in flat for t in i.targets_copy() if t.is_qubit_target})
    assert all(len(i.targets_copy()) == n_q for i in ch)
