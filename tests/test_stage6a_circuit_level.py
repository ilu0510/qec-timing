"""Stage 6a: circuit-level Stim twin with a switchable idle channel.

Every stochastic check uses fixed seeds and an exact (Clopper-Pearson) binomial
interval, or a two-proportion z-test when two estimates are compared.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
import stim

from qec_timing.circuit import build_surface_code_program
from qec_timing.circuit.ordering import order_supports
from qec_timing.circuit.program import _supports
from qec_timing.circuit_level import (
    BitFlip,
    CircuitNoise,
    DephasingH2,
    DepolarizingKOps,
    ReadoutNoise,
    build_circuit_level_memory,
    build_constant_interval_memory,
    compose_depolarize1,
    dephasing_h2_probability,
    idle_probe_circuit,
    rate_per_unit_time,
    readout_probe_circuit,
    sample_measurements,
)
from qec_timing.circuit_level.noise import H_NATIVE_1Q_OPS
from qec_timing.layout import RotatedSurfaceCode
from qec_timing.reference import (
    build_decoder,
    build_memory_circuit,
    constant_schedule,
    estimate_logical_error,
    p_data,
)
from qec_timing.reference.sampler import clopper_pearson

CONF = 0.999


def _within_ci(k: int, n: int, target: float) -> bool:
    low, high = clopper_pearson(k, n, CONF)
    return low <= target <= high


def _error_set(dem: stim.DetectorErrorModel) -> list[tuple[float, str]]:
    """Error mechanisms only (no coordinates), order-independent."""
    out = []
    for inst in dem.flattened():
        if inst.type == "error":
            targets = " ".join(str(t) for t in inst.targets_copy())
            out.append((round(inst.args_copy()[0], 12), targets))
    return sorted(out)


def _two_proportion_z(k1: int, n1: int, k2: int, n2: int) -> float:
    p = (k1 + k2) / (n1 + n2)
    se = math.sqrt(p * (1 - p) * (1 / n1 + 1 / n2))
    return abs(k1 / n1 - k2 / n2) / se if se else 0.0


FULL_A = CircuitNoise.depolarizing_kops(
    k=2.0, p1=1e-3, p2=5e-3, p_init=1e-3, readout=ReadoutNoise(1e-3, 3e-3))
FULL_B = CircuitNoise.dephasing_h2(
    f_quad=0.1, f_lin=1e-3, p1=1e-3, p2=5e-3, p_init=1e-3,
    readout=ReadoutNoise.symmetric(2e-3))


# ---------------------------------------------------------------------------
# Structure: CNOT order, native charging, distance
# ---------------------------------------------------------------------------


def _twin_cx_pairs(d, basis, order, n_rounds=1):
    m = build_constant_interval_memory(d, 1.0, float(n_rounds), CircuitNoise(BitFlip(0.0)),
                                       basis=basis, order=order, with_coords=False)
    return [pair for inst in m.circuit.flattened()
            if inst.name == "CX" and not t_is_rec(inst) for pair in _cx_pairs(inst)]


@pytest.mark.parametrize("order", ["row_major", "hook_safe"])
@pytest.mark.parametrize("basis", ["Z", "X"])
def test_cnot_order_follows_guppy_program(basis, order):
    """First round's CNOTs: active ci then complementary ci, program.py supports."""
    d = 3
    layout = RotatedSurfaceCode(d)
    m = build_constant_interval_memory(d, 1.0, 1.0, CircuitNoise(BitFlip(0.0)),
                                       basis=basis, order=order, with_coords=False)
    anc = layout.n_data
    x_mem = basis == "X"
    a_checks, o_checks = (layout.x_checks, layout.z_checks) if x_mem else (layout.z_checks, layout.x_checks)
    a_type, o_type = ("X", "Z") if x_mem else ("Z", "X")
    act = _supports(layout, order_supports(layout, a_checks, a_type, order))
    oth = _supports(layout, order_supports(layout, o_checks, o_type, order))

    def pairs(table, ci, is_x):
        return [(anc, i) if is_x else (i, anc) for i in table[ci * 4:ci * 4 + 4] if i >= 0]

    prep = [p for ci in range(layout.n_checks) for p in pairs(oth, ci, not x_mem)]
    round0 = [p for ci in range(layout.n_checks)
              for p in pairs(act, ci, x_mem) + pairs(oth, ci, not x_mem)]
    got = [pair for inst in m.circuit.flattened()
           if inst.name == "CX" and not t_is_rec(inst) for pair in _cx_pairs(inst)]
    assert got == prep + round0


def _cx_pairs(inst: stim.CircuitInstruction) -> list[tuple[int, int]]:
    """Stim fuses adjacent CX instructions; split back into (control, target)."""
    return [(g[0].value, g[1].value) for g in inst.target_groups()]


def t_is_rec(inst: stim.CircuitInstruction) -> bool:
    return inst.targets_copy()[0].is_measurement_record_target


def test_native_charging_annotation():
    """H -> 2 composed ops; CX -> target 1+2, control 0+1 single-qubit ops, + p2."""
    p1, p2 = 1e-3, 4e-3
    noise = CircuitNoise(BitFlip(0.0), p1=p1, p2=p2)
    plain = build_constant_interval_memory(3, 1.0, 1.0, noise, basis="X").circuit
    native = build_constant_interval_memory(3, 1.0, 1.0, noise, basis="X",
                                            native_charging=True).circuit

    def targets(c, name, p=None):
        """Number of qubit targets of gate ``name`` (with argument ``p``)."""
        return sum(len(i.targets_copy()) for i in c.flattened() if i.name == name
                   and (p is None or i.gate_args_copy()[0] == p)
                   and not i.targets_copy()[0].is_measurement_record_target)

    n_h = targets(plain, "H")
    n_cx = targets(plain, "CX") // 2
    # Plain: DEPOLARIZE1(p1) on every H target, DEPOLARIZE2(p2) on every CX pair.
    assert targets(plain, "DEPOLARIZE1") == targets(plain, "DEPOLARIZE1", p1) == n_h
    assert targets(plain, "DEPOLARIZE2", p2) == 2 * n_cx
    # Native (measured in Stage 6b): H is two composed ops; each CX charges its
    # target 1 op before and 2 after, its control 1 after, plus the p2 channel.
    p2ops = compose_depolarize1(p1, 2)
    assert targets(native, "DEPOLARIZE1", p2ops) == n_h + n_cx   # H, and target-after
    assert targets(native, "DEPOLARIZE1", p1) == 2 * n_cx         # target-before, control-after
    assert targets(native, "DEPOLARIZE2", p2) == 2 * n_cx


@pytest.mark.parametrize("d", [3, 5])
def test_program_py_order_reduces_x_memory_circuit_distance(d):
    """Documents a property of the teammate's order, not a bug in the twin.

    program.py uses row-major (NW, NE, SW, SE) supports for BOTH sectors. Under
    gate noise, a Z-check ancilla fault mid-check leaves a horizontal pair of Z
    errors on data, parallel to the X memory's logical Z (a row). The Z memory
    keeps distance d; the X memory -- the one variant B needs -- does not.
    """
    z = build_constant_interval_memory(d, 1.0, 3.0, FULL_A, basis="Z").circuit
    x = build_constant_interval_memory(d, 1.0, 3.0, FULL_A, basis="X").circuit
    assert len(z.shortest_graphlike_error()) == d
    assert len(x.shortest_graphlike_error()) < d


HOOK_SAFE_DISTANCE_CASES = [(d, basis, noise, native)
                            for d in (3, 5, 7) for basis in ("Z", "X")
                            for noise in ("A", "B") for native in (False, True)]


@pytest.mark.parametrize("d,basis,variant,native", HOOK_SAFE_DISTANCE_CASES)
def test_hook_safe_order_keeps_full_circuit_distance(d, basis, variant, native):
    """order="hook_safe": graphlike circuit distance is d in BOTH memories."""
    noise = FULL_A if variant == "A" else FULL_B
    m = build_constant_interval_memory(d, 1.0, 3.0, noise, basis=basis,
                                       order="hook_safe", native_charging=native)
    assert len(m.circuit.shortest_graphlike_error()) == d


@pytest.mark.slow
@pytest.mark.parametrize("order", ["hook_safe", "row_major"])
@pytest.mark.parametrize("basis", ["Z", "X"])
def test_twin_cnot_pairs_match_compiled_guppy_program(basis, order):
    """Qubit for qubit against the Guppy program as Selene actually executes it.

    One noiseless Selene shot with tracing; each CNOT compiles to one ZZPhase on
    (ancilla, data), so the executed ZZPhase sequence is the CNOT sequence. The
    ancilla is qubit n_data (the only free slot); direction is fixed by the
    check type in both circuits, so the data-qubit sequence decides the match.
    """
    d = 3
    n_data = RotatedSurfaceCode(d).n_data
    prog = build_surface_code_program(d, basis=basis, order=order)
    res = prog.emulator().with_shots(1).with_trace().run(
        p_data=0.0, p_read=0.0, dt=1.0, n_rounds=1, base_seed=1,
        inject_round=-1, inject_qubit=-1, inject_patch=0)
    zz = [[q.index[0] for q in cmd.qubits] for cmd in res.circuits()[0].get_commands()
          if str(cmd.op.type) == "OpType.ZZPhase"]
    assert all(n_data in pair for pair in zz)
    guppy_data = [next(q for q in pair if q != n_data) for pair in zz]

    twin = _twin_cx_pairs(d, basis, order)
    assert all(n_data in pair for pair in twin)
    twin_data = [next(q for q in pair if q != n_data) for pair in twin]
    assert guppy_data == twin_data


# ---------------------------------------------------------------------------
# 1. Noiseless
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("native", [False, True])
@pytest.mark.parametrize("basis", ["Z", "X"])
@pytest.mark.parametrize("idle", [DepolarizingKOps(3.0, 0.0), DephasingH2(0.0, 0.0)])
def test_noiseless_has_no_detections_or_logical_errors(idle, basis, native):
    d, shots = 5, 5_000
    m = build_constant_interval_memory(d, 0.5, 3.0, CircuitNoise(idle), basis=basis,
                                       native_charging=native)
    assert m.num_detectors == (m.n_rounds + 1) * RotatedSurfaceCode(d).n_checks
    assert len(_error_set(m.detector_error_model)) == 0
    dets, obs = m.sample(shots, seed=11)
    assert not dets.any() and not obs.any()

    # Preparation corrections work: complementary checks read +1 every round.
    n = RotatedSurfaceCode(d).n_checks
    meas = m.sampling_circuit.compile_sampler(seed=12).sample(200)
    other = [n + j * 2 * n + 2 * ci + 1 for j in range(m.n_rounds) for ci in range(n)]
    assert not meas[:, other].any()


# ---------------------------------------------------------------------------
# 2. Bridge to Stage 4a: reduction to the phenomenological reference
# ---------------------------------------------------------------------------

BRIDGE_POINTS = [  # (d, p, lam, dt, T, b_read)
    (3, 0.02, 1.0, 0.5, 5.0, 0.0),
    (5, 0.03, 1.0, 1.0, 10.0, 0.0),
    (3, 0.03, 1.0, 2.0, 10.0, 0.0),
    (5, 0.02, 1.0, 1.0, 10.0, 1.0),  # readout on as well: also identical
]


@pytest.mark.parametrize("d,p,lam,dt,T,b_read", BRIDGE_POINTS)
def test_reduces_to_phenomenological_reference_stage4a_bridge(d, p, lam, dt, T, b_read):
    """Gate/init noise off, idle = X_ERROR(p_data): equals qec_timing.reference.

    The idle X_ERROR probability is the reference's own Eq. 1-2 ``p_data`` for
    (p, lam, dt), so dt enters exactly as in Stage 4a. Two checks:
    (a) the detector error models contain identical error mechanisms, and
    (b) independent samples of each give logical error rates within 4 sigma.
    """
    pd_ = p_data(p, lam, dt)
    p_read = b_read * p
    sched = constant_schedule(dt, T, lam)
    twin = build_circuit_level_memory(
        d, sched, CircuitNoise(BitFlip(pd_), readout=ReadoutNoise.symmetric(p_read)),
        basis="Z")
    ref = build_memory_circuit(RotatedSurfaceCode(d), [pd_] * sched.n_rounds, p_read)

    assert twin.num_detectors == ref.num_detectors
    assert _error_set(twin.detector_error_model) == _error_set(
        ref.detector_error_model(decompose_errors=True))

    shots = 20_000
    a = twin.estimate_logical_error(shots, seed=101)
    b = estimate_logical_error(ref, build_decoder(ref), shots=shots, seed=202)
    assert 0.01 < a.p_L < 0.4, "bridge point should be unsaturated and resolvable"
    assert _two_proportion_z(a.failures, a.shots, b.failures, b.shots) < 4.0, (a, b)


# ---------------------------------------------------------------------------
# 3. Idle channel rates
# ---------------------------------------------------------------------------


def test_depolarize1_composition_formula_is_exact():
    """n sequential DEPOLARIZE1(p) == one DEPOLARIZE1(compose(p, n)) in Stim's DEM."""
    p = 0.03
    for n in (1, 2, 5):
        seq = stim.Circuit()
        seq.append("R", [0])
        for _ in range(n):
            seq.append("DEPOLARIZE1", [0], p)
        seq += stim.Circuit("M 0\nDETECTOR rec[-1]")
        flip = seq.detector_error_model(approximate_disjoint_errors=True)
        total = 0.0
        for inst in flip.flattened():  # XOR-combine independent mechanisms
            q = inst.args_copy()[0]
            total = total + q - 2 * total * q
        assert math.isclose(total, 2 / 3 * compose_depolarize1(p, n), rel_tol=1e-9)
    # Semigroup in the real exponent, so non-integer k*dt is consistent.
    a, b = 0.7, 1.9
    lhs = 1 - 4 / 3 * compose_depolarize1(p, a + b)
    rhs = (1 - 4 / 3 * compose_depolarize1(p, a)) * (1 - 4 / 3 * compose_depolarize1(p, b))
    assert math.isclose(lhs, rhs, rel_tol=1e-12)


@pytest.mark.parametrize("dt", [0.5, 1.0, 2.3])  # k*dt = 1.5, 3, 6.9
@pytest.mark.parametrize("basis", ["Z", "X"])
def test_variant_a_idle_rate(basis, dt):
    """Depolarizing: flips in either basis at 2/3 of the composed probability."""
    idle = DepolarizingKOps(k=3.0, p_op=0.01)
    n_q, shots = 100, 4_000
    m = idle_probe_circuit(idle, dt, basis, n_q).compile_sampler(seed=7).sample(shots)
    target = 2 / 3 * compose_depolarize1(0.01, 3.0 * dt)
    assert _within_ci(int(m.sum()), m.size, target), (m.mean(), target)


@pytest.mark.parametrize("f_quad,f_lin,dt", [(0.2, 0.0, 1.0), (0.1, 0.01, 2.0),
                                             (0.0, 0.02, 3.5)])
def test_variant_b_idle_rate(f_quad, f_lin, dt):
    """Dephasing: Z at p_z(dt) (X-basis flips), and no X at all (Z basis)."""
    idle = DephasingH2(f_quad, f_lin)
    n_q, shots = 100, 4_000
    target = math.sin(0.5 * f_quad * dt) ** 2 + f_lin * dt
    assert idle.probability(dt) == pytest.approx(target)
    mx = idle_probe_circuit(idle, dt, "X", n_q).compile_sampler(seed=8).sample(shots)
    assert _within_ci(int(mx.sum()), mx.size, target), (mx.mean(), target)
    mz = idle_probe_circuit(idle, dt, "Z", n_q).compile_sampler(seed=9).sample(shots)
    assert not mz.any()


def test_variant_b_clipping():
    assert dephasing_h2_probability(1.0, 0.0, 0.7) == 0.5
    assert dephasing_h2_probability(1.0, math.pi, 0.1) == 0.5
    assert dephasing_h2_probability(0.0, 0.3, 0.1) == 0.0


# ---------------------------------------------------------------------------
# 4. Basis check: dephasing is invisible to a Z memory
# ---------------------------------------------------------------------------


def test_dephasing_invisible_to_z_memory_visible_to_x_memory():
    d, shots = 3, 5_000
    idle_only = CircuitNoise.dephasing_h2(f_quad=0.0, f_lin=0.05)
    z = build_constant_interval_memory(d, 1.0, 5.0, idle_only, basis="Z")
    x = build_constant_interval_memory(d, 1.0, 5.0, idle_only, basis="X")
    dz, oz = z.sample(shots, seed=21)
    assert not dz.any() and not oz.any()
    assert z.estimate_logical_error(shots, seed=22).failures == 0
    assert x.estimate_logical_error(shots, seed=23).p_L > 0.05


def test_dephasing_adds_nothing_to_z_memory_under_gate_noise():
    """With gate noise on, the Z memory's error model is unchanged by dephasing."""
    gates = dict(p1=1e-3, p2=5e-3, p_init=1e-3, readout=ReadoutNoise.symmetric(2e-3))
    on = CircuitNoise.dephasing_h2(f_quad=0.3, f_lin=0.02, **gates)
    off = CircuitNoise.dephasing_h2(f_quad=0.0, f_lin=0.0, **gates)
    for native in (False, True):
        kw = dict(basis="Z", native_charging=native)
        z_on = build_constant_interval_memory(3, 1.0, 4.0, on, **kw)
        z_off = build_constant_interval_memory(3, 1.0, 4.0, off, **kw)
        assert _error_set(z_on.detector_error_model) == _error_set(z_off.detector_error_model)
    x_on = build_constant_interval_memory(3, 1.0, 4.0, on, basis="X")
    x_off = build_constant_interval_memory(3, 1.0, 4.0, off, basis="X")
    assert _error_set(x_on.detector_error_model) != _error_set(x_off.detector_error_model)


# ---------------------------------------------------------------------------
# 5. Asymmetric readout
# ---------------------------------------------------------------------------


def test_asymmetric_readout_rates():
    p01, p10 = 0.02, 0.08
    n0, n1, shots = 200, 200, 3_000
    m = sample_measurements(readout_probe_circuit(n0, n1), ReadoutNoise(p01, p10),
                            shots, seed=31)
    zeros, ones = m[:, :n0], m[:, n0:]
    assert _within_ci(int(zeros.sum()), zeros.size, p01)          # reported 1 | true 0
    assert _within_ci(int((~ones).sum()), ones.size, p10)         # reported 0 | true 1


def test_asymmetric_readout_in_memory_depends_on_true_outcome():
    """Noiseless active syndromes are 0, so p01 alone fires round detectors and
    p10 alone cannot; the final data readout has 1s, so p10 fires there."""
    d = 3
    n = RotatedSurfaceCode(d).n_checks
    only_10 = build_constant_interval_memory(
        d, 1.0, 4.0, CircuitNoise(BitFlip(0.0), readout=ReadoutNoise(0.0, 0.1)))
    only_01 = build_constant_interval_memory(
        d, 1.0, 4.0, CircuitNoise(BitFlip(0.0), readout=ReadoutNoise(0.1, 0.0)))
    assert only_10.sampling_circuit is not only_10.circuit
    dets, _ = only_10.sample(2_000, seed=41)
    assert not dets[:, :-n].any() and dets[:, -n:].any()
    dets, _ = only_01.sample(2_000, seed=42)
    assert dets[:, :-n].any()
    # Decoder weights use the symmetric average.
    ms = [i for i in only_10.circuit.flattened() if i.name == "M"]
    assert {i.gate_args_copy()[0] for i in ms} == {0.05}


# ---------------------------------------------------------------------------
# Bookkeeping
# ---------------------------------------------------------------------------


def test_round_count_and_metadata_use_shared_schedule():
    m = build_constant_interval_memory(3, 0.1, 1.0, FULL_A)
    assert m.n_rounds == 10  # round(T/dt), not a float accumulator
    meta = m.round_metadata()
    assert len(meta) == 10
    assert meta[3]["t"] == pytest.approx(0.3) and meta[3]["dt"] == pytest.approx(0.1)
    assert meta[0]["p_idle"] == pytest.approx(compose_depolarize1(1e-3, 0.2))
    assert rate_per_unit_time(0.0, 1.0) == 0.0


# ---------------------------------------------------------------------------
# 6. Smoke run (no conclusions)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name,noise", [("depolarizing_kops", FULL_A),
                                        ("dephasing_h2", FULL_B)])
def test_smoke_d3_three_intervals(name, noise, capsys):
    T, shots = 6.0, 4_000
    lines = []
    for dt in (0.5, 1.0, 2.0):
        m = build_constant_interval_memory(3, dt, T, noise)
        est = m.estimate_logical_error(shots, seed=51)
        R = rate_per_unit_time(est.p_L, m.schedule.total_time)
        lines.append(f"{name} basis={m.basis} dt={dt} rounds={m.n_rounds} "
                     f"p_L={est.p_L:.4f} R={R:.5f}/unit time")
        assert 0.0 <= est.p_L <= 1.0
    with capsys.disabled():
        print("\n" + "\n".join(lines))
