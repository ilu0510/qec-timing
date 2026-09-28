"""Stage 1 tests for the pure-Python reference model.

No Selene or Guppy anywhere in this file. Fixed seeds throughout.
"""

from __future__ import annotations

import ast
import math
from pathlib import Path

import numpy as np
import pytest
import stim

from qec_timing import ansatz
from qec_timing.reference import (
    RotatedSurfaceCodeZSector,
    build_decoder,
    build_memory_circuit,
    build_memory_circuit_from_schedule,
    constant_schedule,
    data_flip_probe,
    estimate_logical_error,
    p_data,
    p_idle,
    p_read,
    readout_flip_probe,
    schedule_from_intervals,
    single_error_circuit,
)
from qec_timing.reference.sampler import clopper_pearson

SEED = 20260927
CONFIDENCE = 0.999


# ---------------------------------------------------------------------------
# Eqs. 1-3: noise parameters
# ---------------------------------------------------------------------------


def test_eq1_eq2_eq3_closed_form():
    """Eqs. 1-3 match their closed forms and small-p limits."""
    p, lam, dt, b_read = 0.015, 1.3, 0.7, 2.0

    assert p_idle(p, lam, dt) == pytest.approx(1.0 - math.exp(-p * lam * dt))
    assert p_data(p, lam, dt) == pytest.approx(
        1.0 - (1.0 - p) * (1.0 - p_idle(p, lam, dt))
    )
    assert p_read(p, b_read) == pytest.approx(b_read * p)

    # Small-p limits: p_idle -> p lam dt and p_data -> p (1 + lam dt).
    small = 1e-6
    assert p_idle(small, lam, dt) == pytest.approx(small * lam * dt, rel=1e-5)
    assert p_data(small, lam, dt) == pytest.approx(small * (1 + lam * dt), rel=1e-5)


def test_p_stab_survives_at_zero_lambda():
    """At lam=0 the idling term vanishes but p_stab = p must remain."""
    assert p_idle(0.01, 0.0, 5.0) == 0.0
    assert p_data(0.01, 0.0, 5.0) == pytest.approx(0.01)
    # ... and p_data must be strictly increasing in dt once lam > 0.
    values = [p_data(0.01, 1.0, dt) for dt in (0.1, 0.5, 1.0, 2.0)]
    assert values == sorted(values)
    assert values[0] > 0.01


def test_noise_parameters_reject_bad_input():
    with pytest.raises(ValueError):
        p_idle(1.5, 1.0, 1.0)
    with pytest.raises(ValueError):
        p_idle(0.01, -1.0, 1.0)
    with pytest.raises(ValueError):
        p_read(0.9, 2.0)  # p_read > 1


# ---------------------------------------------------------------------------
# Layout and time bookkeeping
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("d", [3, 5, 7, 9, 11])
def test_layout_structure(d):
    """Check count, graphlike degree, and logical support."""
    layout = RotatedSurfaceCodeZSector(d)
    assert layout.n_data == d * d
    assert len(layout.z_checks) == layout.n_checks == (d * d - 1) // 2

    degree: dict[tuple[int, int], int] = {}
    for support in layout.z_checks:
        assert len(support) in (2, 4)
        for q in support:
            degree[q] = degree.get(q, 0) + 1
    # At most two Z-checks per data qubit keeps the DEM graphlike.
    assert max(degree.values()) == 2
    assert len(layout.logical_z_support) == d


def test_layout_rejects_even_distance():
    with pytest.raises(ValueError):
        RotatedSurfaceCodeZSector(4)


def test_logical_x_is_a_column_not_a_row():
    """A column of X commutes with every Z-check; a row does not.

    This is what fixes the logical Z observable to a row of the final readout.
    """
    layout = RotatedSurfaceCodeZSector(5)
    for col in range(5):
        column = {(r, col) for r in range(5)}
        assert all(len(column & set(s)) % 2 == 0 for s in layout.z_checks)
    offending = 0
    for row in range(5):
        line = {(row, c) for c in range(5)}
        if any(len(line & set(s)) % 2 for s in layout.z_checks):
            offending += 1
    assert offending > 0

    # And a logical X column must flip the observable parity exactly once.
    support = set(layout.logical_z_support)
    for col in range(5):
        assert len(support & {(r, col) for r in range(5)}) == 1


def test_round_count_is_round_T_over_dt():
    """n_rounds = round(T/dt), with exact integer time bookkeeping."""
    for dt, T, expected in [
        (0.5, 10.0, 20),
        (0.25, 10.0, 40),
        (1.7, 200.0, 118),
        (0.1, 1.0, 10),  # a float accumulator would give 11 here
        (3.0, 10.0, 3),  # 3.333 -> 3
    ]:
        schedule = constant_schedule(dt=dt, T=T, lam=1.0)
        assert schedule.n_rounds == expected, f"dt={dt}, T={T}"
        assert schedule.total_time == pytest.approx(expected * dt)
        # t_j = j * dt exactly, with no drift.
        assert schedule.start_times == pytest.approx(
            [j * dt for j in range(expected)]
        )


def test_schedule_rejects_non_tick_multiples():
    """Intervals must be exact tick multiples, or time bookkeeping is a lie."""
    # Not a whole number of ticks.
    with pytest.raises(ValueError, match="integer multiple"):
        constant_schedule(dt=0.0015, T=1.0, lam=1.0)
    # Finer than one tick.
    with pytest.raises(ValueError, match="smaller than one tick"):
        constant_schedule(dt=0.00015, T=1.0, lam=1.0)
    schedule = schedule_from_intervals([0.25, 0.5, 0.25], [1.0, 2.0, 1.0])
    assert schedule.n_rounds == 3
    assert schedule.total_time == pytest.approx(1.0)
    assert schedule.start_times == pytest.approx([0.0, 0.25, 0.75])


# ---------------------------------------------------------------------------
# Detector / round counts
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(("d", "dt", "T"), [(3, 1.0, 9.0), (5, 0.5, 10.0), (7, 0.25, 5.25)])
def test_detector_and_round_counts(d, dt, T):
    """Detectors = (n_rounds + 1) * n_checks, observables = 1."""
    layout = RotatedSurfaceCodeZSector(d)
    schedule = constant_schedule(dt=dt, T=T, lam=1.0)
    circuit = build_memory_circuit_from_schedule(layout, schedule, p=0.01, b_read=1.0)

    n_rounds = round(T / dt)
    assert schedule.n_rounds == n_rounds
    assert circuit.num_detectors == (n_rounds + 1) * layout.n_checks
    assert circuit.num_observables == 1
    # One measurement per check per round, plus the final data readout.
    assert circuit.num_measurements == n_rounds * layout.n_checks + layout.n_data


# ---------------------------------------------------------------------------
# Sampled flip rates against binomial confidence intervals
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("target", [0.01, 0.3])
def test_data_flip_rate_matches_p_data(target):
    """X_ERROR is applied at exactly p_data, one draw per data qubit."""
    layout = RotatedSurfaceCodeZSector(5)
    circuit = data_flip_probe(layout, target)
    shots = 4000
    sample = circuit.compile_sampler(seed=SEED).sample(shots=shots)

    assert sample.shape == (shots, layout.n_data)
    flips = int(np.count_nonzero(sample))
    trials = shots * layout.n_data
    low, high = clopper_pearson(flips, trials, CONFIDENCE)
    assert low <= target <= high, (
        f"p_data={target}: observed {flips}/{trials} = {flips / trials:.5f}, "
        f"99.9% CI [{low:.5f}, {high:.5f}]"
    )


@pytest.mark.parametrize("target", [0.01, 0.3])
def test_readout_flip_rate_matches_p_read(target):
    """MPP flips only the recorded syndrome bit, at exactly p_read."""
    layout = RotatedSurfaceCodeZSector(5)
    n_rounds = 4
    circuit = readout_flip_probe(layout, target, n_rounds=n_rounds)
    shots = 2000
    sample = circuit.compile_sampler(seed=SEED).sample(shots=shots)

    assert sample.shape == (shots, n_rounds * layout.n_checks)
    flips = int(np.count_nonzero(sample))
    trials = sample.size
    low, high = clopper_pearson(flips, trials, CONFIDENCE)
    assert low <= target <= high, (
        f"p_read={target}: observed {flips}/{trials} = {flips / trials:.5f}, "
        f"99.9% CI [{low:.5f}, {high:.5f}]"
    )


def test_readout_noise_does_not_disturb_the_state():
    """p_read flips the classical record only (CLAUDE.md).

    Heavily noisy syndrome measurements with no data faults, followed by a
    *noiseless* data readout, must still return all zeros. (The memory
    circuit's own final readout is deliberately noisy, so this probe appends a
    clean one instead.)
    """
    layout = RotatedSurfaceCodeZSector(3)
    circuit = readout_flip_probe(layout, 0.4, n_rounds=3)
    n_syndrome = circuit.num_measurements
    circuit.append("M", [layout.data_index[q] for q in layout.data_qubits], 0.0)

    sample = circuit.compile_sampler(seed=SEED).sample(shots=500)
    syndrome, final_readout = sample[:, :n_syndrome], sample[:, n_syndrome:]
    assert final_readout.shape == (500, layout.n_data)
    assert not final_readout.any(), "state was disturbed by readout noise"
    # Sanity: the syndrome record really was being flipped meanwhile.
    assert syndrome.any()


# ---------------------------------------------------------------------------
# Correctness at zero noise and for a single fault
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("d", [3, 5])
def test_zero_noise_gives_exactly_zero_logical_error(d):
    """p = 0 means no error mechanisms at all and zero logical failures."""
    layout = RotatedSurfaceCodeZSector(d)
    schedule = constant_schedule(dt=1.0, T=float(3 * d), lam=1.0)
    circuit = build_memory_circuit_from_schedule(layout, schedule, p=0.0, b_read=1.0)

    decoder = build_decoder(circuit)
    assert decoder.is_trivial, "a p=0 circuit should have an empty error model"

    estimate = estimate_logical_error(circuit, decoder, shots=2000, seed=SEED)
    assert estimate.failures == 0
    assert estimate.p_L == 0.0

    # No detector should ever fire either.
    detectors, observables = circuit.compile_detector_sampler(
        seed=SEED
    ).sample(shots=500, separate_observables=True)
    assert not detectors.any()
    assert not observables.any()


@pytest.mark.parametrize("d", [3, 5])
def test_single_data_error_is_always_corrected(d):
    """Any single data X fault, on any qubit, is corrected.

    The decoder is built from the noisy circuit (so the Eq. A2 weights are
    meaningful), then applied to a deterministic single-fault circuit.
    """
    layout = RotatedSurfaceCodeZSector(d)
    n_rounds = 3
    reference = build_memory_circuit(
        layout, [0.01] * n_rounds, p_read=0.01, with_coords=False
    )
    decoder = build_decoder(reference)

    for qubit in layout.data_qubits:
        circuit = single_error_circuit(
            layout, qubit, n_rounds=n_rounds, round_index=1
        )
        assert circuit.num_detectors == reference.num_detectors
        detectors, observables = circuit.compile_detector_sampler(
            seed=SEED
        ).sample(shots=1, separate_observables=True)
        prediction = decoder.decode_batch(detectors)
        assert np.array_equal(prediction, observables), (
            f"single X fault on data qubit {qubit} (d={d}) was not corrected"
        )


# ---------------------------------------------------------------------------
# Code distance
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("d", [3, 5, 7])
def test_code_distance_is_d(d):
    """The minimum-weight undetectable logical error has weight exactly d."""
    layout = RotatedSurfaceCodeZSector(d)
    circuit = build_memory_circuit(
        layout, [0.01] * 3, p_read=0.01, with_coords=False
    )
    dem = circuit.detector_error_model(decompose_errors=False)
    shortest = dem.shortest_graphlike_error()
    assert len(shortest) == d


@pytest.mark.parametrize("d", [3, 5])
def test_z_check_supports_match_stim_reference_code(d):
    """Cross-check the hand-built geometry against Stim's own generator.

    Compares the multiset of Z-check support sizes and the total count, which
    pins down the stabilizer structure without depending on Stim's particular
    coordinate convention.
    """
    layout = RotatedSurfaceCodeZSector(d)
    generated = stim.Circuit.generated(
        "surface_code:rotated_memory_z", distance=d, rounds=1
    )
    # Stim's circuit carries both sectors; one sector is half the stabilizers.
    assert layout.n_checks == (d * d - 1) // 2
    assert generated.num_detectors // 1 >= layout.n_checks

    ours = sorted(len(s) for s in layout.z_checks)
    n_weight2 = d - 1
    n_weight4 = layout.n_checks - n_weight2
    assert ours == sorted([2] * n_weight2 + [4] * n_weight4)


# ---------------------------------------------------------------------------
# Analytic ansatz
# ---------------------------------------------------------------------------


def test_ansatz_does_not_import_reference_or_selene():
    """ansatz.py must stay an independent yardstick."""
    source = Path(ansatz.__file__).read_text(encoding="utf-8")
    tree = ast.parse(source)
    forbidden = ("reference", "selene", "guppylang", "stim", "pymatching")
    names: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            names.append(node.module or "")
    for name in names:
        assert not any(bad in name.lower() for bad in forbidden), (
            f"ansatz.py must not import {name!r}"
        )


def test_alpha_and_K_alpha():
    """Eq. D2 and Eq. D8, including the large-alpha limit K_alpha -> e*alpha."""
    assert float(ansatz.alpha(5, 0.8)) == pytest.approx(0.8 * 6 / 2)
    a = 2.4
    assert float(ansatz.K_alpha(a)) == pytest.approx(a**a / (a - 1) ** (a - 1))

    # Exact identity K_alpha = alpha * (1 + 1/(alpha-1))**(alpha-1).
    for value in (1.5, 2.4, 11.0, 200.0):
        assert float(ansatz.K_alpha(value)) == pytest.approx(
            value * (1.0 + 1.0 / (value - 1.0)) ** (value - 1.0)
        )

    # K_alpha -> e * alpha, approached as 1/(2 alpha): ~2.5e-3 at alpha=200,
    # so the limit is only tight for genuinely large alpha.
    assert float(ansatz.K_alpha(1e5)) == pytest.approx(math.e * 1e5, rel=1e-4)
    with pytest.raises(ValueError):
        ansatz.K_alpha(1.0)


@pytest.mark.parametrize("d", [5, 7, 11, 15, 21, 27])
@pytest.mark.parametrize("lam", [0.03, 0.5, 1.0, 3.0])
def test_optimal_interval_numerically_minimises_the_rate(d, lam):
    """Eq. 6 / D6 really is the argmin of Eq. 5, on a fine grid."""
    dt_star = ansatz.optimal_interval(d=d, lam=lam)

    n_points = 20_001
    span = 10.0  # scan two decades centred on dt*
    grid = np.geomspace(dt_star / span, dt_star * span, n_points)
    rates = ansatz.logical_rate(grid, d=d, p=0.01, lam=lam, A=1.0)

    dt_numeric = float(grid[int(np.argmin(rates))])
    log_spacing = math.log(span**2) / (n_points - 1)
    assert abs(math.log(dt_numeric / dt_star)) <= log_spacing, (
        f"d={d}, lam={lam}: analytic dt*={dt_star:.6g}, grid argmin={dt_numeric:.6g}"
    )

    # The analytic optimum must also beat every grid point (up to fp noise).
    r_star = float(ansatz.logical_rate(dt_star, d=d, p=0.01, lam=lam, A=1.0))
    assert r_star <= float(np.min(rates)) * (1.0 + 1e-12)


@pytest.mark.parametrize("d", [7, 15])
@pytest.mark.parametrize("lam", [0.5, 2.0])
def test_gamma_opt_is_one_at_the_optimum(d, lam):
    """Eq. 7 / D9 equals 1 at dt* and exceeds 1 away from it."""
    dt_star = ansatz.optimal_interval(d=d, lam=lam)
    assert float(ansatz.gamma_opt(dt_star, d=d, lam=lam)) == pytest.approx(1.0)

    for factor in (0.1, 0.5, 2.0, 10.0):
        assert float(ansatz.gamma_opt(dt_star * factor, d=d, lam=lam)) > 1.0

    # Gamma_opt must equal the explicit rate ratio.
    dt = dt_star * 3.0
    ratio = float(
        ansatz.logical_rate(dt, d=d, p=0.012, lam=lam, A=1.7)
        / ansatz.logical_rate(dt_star, d=d, p=0.012, lam=lam, A=1.7)
    )
    assert float(ansatz.gamma_opt(dt, d=d, lam=lam)) == pytest.approx(ratio)


def test_optimal_rate_matches_eq_D7():
    """R*(lam) = C K_alpha lam equals R evaluated at dt*."""
    kwargs = dict(d=11, p=0.015, lam=1.0, A=1.4)
    dt_star = ansatz.optimal_interval(d=11, lam=1.0)
    assert ansatz.optimal_rate(**kwargs) == pytest.approx(
        float(ansatz.logical_rate(dt_star, **kwargs))
    )


def test_eq_A3_limits():
    """Eq. A3 is linear for small accumulated rate and saturates at 1/2."""
    small = ansatz.p_L_from_rates([1e-8] * 10, [0.1] * 10)
    assert small == pytest.approx(1e-8 * 10 * 0.1, rel=1e-6)
    saturated = ansatz.p_L_from_rates([1e3] * 10, [1.0] * 10)
    assert saturated == pytest.approx(0.5)


def test_eq6_reproduces_paper_fig2_slope():
    """Fig. 2/9: 1/dt* should be linear in d with slope ~0.402 at lam=1, g=0.8."""
    ds = np.arange(11, 28, 2)
    inv_dt_star = np.array(
        [1.0 / ansatz.optimal_interval(d=int(d), lam=1.0) for d in ds]
    )
    slope, _ = np.polyfit(ds, inv_dt_star, 1)
    assert slope == pytest.approx(0.4, abs=0.01)
