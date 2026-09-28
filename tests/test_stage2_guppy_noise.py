"""Stage 2 tests: noise injection inside Guppy, static lambda, Stim backend.

Fixed base seeds throughout. The contract under test is docs/noise_spec.md.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from qec_timing.noise import (
    NoiseParams,
    data_probability,
    idle_probability,
    read_probability,
)
from qec_timing.noise.probes import build_data_flip_probe, record_flip_probe
from qec_timing.layout import RotatedSurfaceCodeZSector
from qec_timing.reference.noise import p_data as ref_p_data
from qec_timing.reference.noise import p_idle as ref_p_idle
from qec_timing.reference.noise import p_read as ref_p_read
from qec_timing.reference.sampler import clopper_pearson
from qec_timing.stub_circuit import (
    WORD_BITS,
    build_memory_program,
    unpack_detection_events,
    unpack_packed_words,
)

BASE_SEED = 20260927
CONFIDENCE = 0.999
FORMATS = ("packed", "sparse", "dense")


# ---------------------------------------------------------------------------
# Eqs. 1-3 on the product path, and agreement with the independent reference
# ---------------------------------------------------------------------------


def test_eqs_1_to_3_closed_form():
    p, lam, dt, b_read = 0.015, 1.3, 0.7, 2.0
    assert idle_probability(p, lam, dt) == pytest.approx(
        1.0 - math.exp(-p * lam * dt)
    )
    assert data_probability(p, lam, dt) == pytest.approx(
        1.0 - (1.0 - p) * (1.0 - idle_probability(p, lam, dt))
    )
    assert read_probability(p, b_read) == pytest.approx(b_read * p)


def test_p_stab_survives_at_zero_lambda():
    """p_data(lam=0) == p exactly: the term that creates the dt trade-off."""
    assert idle_probability(0.02, 0.0, 9.0) == 0.0
    assert data_probability(0.02, 0.0, 9.0) == pytest.approx(0.02)


def test_product_and_reference_noise_agree():
    """The two independent implementations of Eqs. 1-3 must not drift apart."""
    for p in (0.0, 0.001, 0.015, 0.1):
        for lam in (0.0, 1.0, 3.0):
            for dt in (0.01, 0.7, 10.0):
                assert idle_probability(p, lam, dt) == pytest.approx(
                    ref_p_idle(p, lam, dt), abs=1e-15
                )
                assert data_probability(p, lam, dt) == pytest.approx(
                    ref_p_data(p, lam, dt), abs=1e-15
                )
        assert read_probability(p, 1.0) == pytest.approx(ref_p_read(p, 1.0))


def test_not_the_linearised_form():
    """Eq. 1 must be the exponential, not p*lam*dt, where they differ."""
    p, lam, dt = 0.015, 1.0, 10.0
    linearised = p * lam * dt
    assert idle_probability(p, lam, dt) == pytest.approx(0.13929, abs=1e-5)
    assert abs(idle_probability(p, lam, dt) - linearised) > 0.01


def test_noise_params_round_count():
    params = NoiseParams(p=0.01, lam=1.0, dt=0.1, b_read=1.0)
    # round(T/dt), not a float accumulator (which would give 11 here).
    assert params.n_rounds(1.0) == 10
    assert params.realised_T(1.0) == pytest.approx(1.0)
    assert NoiseParams(p=0.01, lam=1.0, dt=1.7, b_read=1.0).n_rounds(200.0) == 118


# ---------------------------------------------------------------------------
# Empirical rates inside Guppy
# ---------------------------------------------------------------------------


@pytest.mark.slow
@pytest.mark.parametrize("p_data_value", [0.001, 0.05, 0.3])
def test_data_flip_rate_matches_p_data(p_data_value):
    """apply_data_noise flips each qubit at exactly p_data."""
    n_data = 25
    probe = build_data_flip_probe(n_data)
    n_reps, n_shots = 400, 100

    result = (
        probe.emulator(n_qubits=n_data)
        .stabilizer_sim()
        .with_shots(n_shots)
        .run(p_data=p_data_value, n_reps=n_reps, base_seed=BASE_SEED)
    )
    shots = [s.as_dict() for s in result.results]
    flips = sum(s["flips"] for s in shots)
    trials = sum(s["trials"] for s in shots)
    assert trials == n_shots * n_reps * n_data

    low, high = clopper_pearson(flips, trials, CONFIDENCE)
    assert low <= p_data_value <= high, (
        f"p_data={p_data_value}: {flips}/{trials} = {flips / trials:.6f}, "
        f"99.9% CI [{low:.6f}, {high:.6f}]"
    )


@pytest.mark.slow
@pytest.mark.parametrize("p_read_value", [0.001, 0.05, 0.3])
def test_record_flip_rate_matches_p_read(p_read_value):
    """flip_record flips the classical record at exactly p_read."""
    n_reps, n_shots = 10_000, 100
    result = (
        record_flip_probe.emulator(n_qubits=1)
        .stabilizer_sim()
        .with_shots(n_shots)
        .run(p_read=p_read_value, n_reps=n_reps, base_seed=BASE_SEED)
    )
    shots = [s.as_dict() for s in result.results]
    flips = sum(s["flips"] for s in shots)
    trials = sum(s["trials"] for s in shots)

    low, high = clopper_pearson(flips, trials, CONFIDENCE)
    assert low <= p_read_value <= high, (
        f"p_read={p_read_value}: {flips}/{trials} = {flips / trials:.6f}, "
        f"99.9% CI [{low:.6f}, {high:.6f}]"
    )


@pytest.mark.slow
def test_zero_probability_never_fires():
    """p = 0 must be incapable of firing, not merely unlikely."""
    result = (
        record_flip_probe.emulator(n_qubits=1)
        .stabilizer_sim()
        .with_shots(20)
        .run(p_read=0.0, n_reps=20_000, base_seed=BASE_SEED)
    )
    assert all(s.as_dict()["flips"] == 0 for s in result.results)


# ---------------------------------------------------------------------------
# Memory program: seeding, determinism, wire formats
# ---------------------------------------------------------------------------


def _events(program, *, p_data, p_read, n_rounds, shots, base_seed):
    result = program.run(
        p_data=p_data,
        p_read=p_read,
        n_rounds=n_rounds,
        base_seed=base_seed,
        shots=shots,
    )
    return program.detection_events(result, n_rounds)


@pytest.fixture(scope="module")
def d3_programs():
    layout = RotatedSurfaceCodeZSector(3)
    return {fmt: build_memory_program(layout, emit=fmt) for fmt in FORMATS}


@pytest.mark.slow
def test_all_wire_formats_agree(d3_programs):
    """The three emitters must describe identical detection events."""
    kwargs = dict(p_data=0.15, p_read=0.1, n_rounds=4, shots=8, base_seed=BASE_SEED)
    per_format = {
        fmt: _events(d3_programs[fmt], **kwargs) for fmt in FORMATS
    }
    reference = per_format["packed"]
    assert sum(int(e.sum()) for e in reference) > 0, "test would be vacuous"
    for fmt in ("sparse", "dense"):
        for got, want in zip(per_format[fmt], reference, strict=True):
            assert np.array_equal(got, want), f"{fmt} disagrees with packed"


@pytest.mark.slow
def test_seeding_is_reproducible_and_varies_per_shot(d3_programs):
    """Same base_seed reproduces exactly; different shots differ."""
    program = d3_programs["packed"]
    kwargs = dict(p_data=0.15, p_read=0.1, n_rounds=4, shots=6)

    first = _events(program, base_seed=BASE_SEED, **kwargs)
    again = _events(program, base_seed=BASE_SEED, **kwargs)
    other = _events(program, base_seed=BASE_SEED + 1, **kwargs)

    assert all(np.array_equal(a, b) for a, b in zip(first, again, strict=True))
    # Per-shot seeding: a fixed seed for every shot would make these identical
    # and silently invalidate every statistic.
    assert not np.array_equal(first[0], first[1])
    assert not all(np.array_equal(a, b) for a, b in zip(first, other, strict=True))


@pytest.mark.slow
def test_zero_noise_gives_no_detection_events(d3_programs):
    """p_data = p_read = 0 gives a completely silent detector history."""
    for fmt in FORMATS:
        events = _events(
            d3_programs[fmt],
            p_data=0.0,
            p_read=0.0,
            n_rounds=5,
            shots=4,
            base_seed=BASE_SEED,
        )
        assert all(not e.any() for e in events), f"{fmt} fired at zero noise"


@pytest.mark.slow
def test_ideal_extraction_leaves_logical_state_undisturbed(d3_programs):
    """With no noise, the logical observable must stay 0 in every shot.

    Ideal extraction must not disturb the encoded state: the ancilla CNOTs copy
    only the Z-parity out.
    """
    program = d3_programs["packed"]
    result = program.run(
        p_data=0.0, p_read=0.0, n_rounds=6, base_seed=BASE_SEED, shots=20
    )
    assert all(shot.as_dict()["obs"] == 0 for shot in result.results)


@pytest.mark.slow
def test_packed_round_trip_against_known_pattern():
    """Pack in Guppy, unpack in Python, compare against the injected pattern.

    d = 11 has 60 checks (one word) and d = 15 has 112 (two words), so this
    covers both a single partial word and a full-plus-partial pair.
    """
    for d in (11, 15):
        layout = RotatedSurfaceCodeZSector(d)
        program = build_memory_program(layout, emit="packed")
        n_rounds = 2
        # p_data = 1.0 flips every data qubit every round, which makes the
        # detection-event pattern deterministic and checkable.
        result = program.run(
            p_data=1.0, p_read=0.0, n_rounds=n_rounds, base_seed=BASE_SEED
        )
        tags = result.results[0].collate_tags()
        events = program.detection_events(result, n_rounds)[0]
        assert events.shape == (n_rounds + 1, layout.n_checks)

        # Every data qubit flips each round, so check c sees its full support
        # flip: its syndrome is parity(weight) and round 0 fires iff the weight
        # is odd. All supports here have even weight (2 or 4), so nothing fires.
        weights = np.array([len(s) for s in layout.z_checks])
        assert np.array_equal(events[0], (weights % 2).astype(bool))

        # And the raw words must round-trip through the documented unpacking.
        words = np.asarray(tags["det"], dtype=np.uint64).reshape(
            n_rounds + 1, program.n_words
        )
        assert np.array_equal(unpack_packed_words(words, layout.n_checks), events)
        assert program.n_words == math.ceil(layout.n_checks / WORD_BITS)


@pytest.mark.slow
def test_single_injected_error_fires_expected_detectors():
    """One X on one data qubit must fire exactly the checks containing it.

    Checked against the geometry rather than a golden value, so the test cannot
    silently agree with a wrong layout.
    """
    layout = RotatedSurfaceCodeZSector(5)
    program = build_memory_program(layout, emit="packed")
    n_rounds = 3

    for target in range(layout.n_data):
        qubit_coord = layout.data_qubits[target]
        expected = np.array(
            [qubit_coord in set(support) for support in layout.z_checks]
        )
        result = program.run(
            p_data=0.0,
            p_read=0.0,
            n_rounds=n_rounds,
            base_seed=BASE_SEED,
            inject_round=1,
            inject_qubit=target,
        )
        events = program.detection_events(result, n_rounds)[0]
        assert not events[0].any(), "nothing should fire before the injection"
        assert np.array_equal(events[1], expected), (
            f"data qubit {qubit_coord} fired {np.flatnonzero(events[1])}, "
            f"expected {np.flatnonzero(expected)}"
        )
        # The fault persists, so the next round sees no *change*.
        assert not events[2].any(), "a static fault must not re-fire"
