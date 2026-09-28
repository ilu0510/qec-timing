"""Stage 3 tests: shared scheduling, per-round metadata, end-to-end plumbing.

Fixed seeds throughout. No statistics are claimed here -- the decode check is
plumbing only.
"""

from __future__ import annotations

import numpy as np
import pytest

from qec_timing.layout import RotatedSurfaceCodeZSector
from qec_timing.noise import NoiseParams, schedule_data_probabilities
from qec_timing.reference.circuit import build_memory_circuit, single_error_circuit
from qec_timing.reference.decoder import build_decoder
from qec_timing.reference.noise import data_probabilities
from qec_timing.schedule import (
    DEFAULT_TICK,
    RoundSchedule,
    constant_schedule,
    schedule_from_intervals,
    ticks_for,
)
from qec_timing.stub_circuit import build_memory_program

BASE_SEED = 20260927

# (d, dt, T) combinations for the end-to-end check.
END_TO_END_CASES = [
    (3, 0.5, 3.0),
    (3, 0.25, 2.0),
    (5, 1.0, 5.0),
]


# ---------------------------------------------------------------------------
# 1. Scheduling
# ---------------------------------------------------------------------------


def test_fixed_dt_round_count_is_exact():
    """dt = 0.1, T = 1.0 gives exactly 10 rounds, not the 11 a float sum gives."""
    schedule = constant_schedule(dt=0.1, T=1.0, lam=1.0)
    assert schedule.n_rounds == 10

    # The failure mode this guards against, demonstrated explicitly.
    t, accumulated_rounds = 0.0, 0
    while t < 1.0:
        t += 0.1
        accumulated_rounds += 1
    assert accumulated_rounds == 11
    assert schedule.n_rounds != accumulated_rounds


@pytest.mark.parametrize(
    ("dt", "T", "expected"),
    [(0.5, 10.0, 20), (0.25, 10.0, 40), (1.7, 200.0, 118), (3.0, 10.0, 3)],
)
def test_round_counts_and_times(dt, T, expected):
    schedule = constant_schedule(dt=dt, T=T, lam=1.0)
    assert schedule.n_rounds == expected
    assert schedule.dt == pytest.approx(dt)
    assert schedule.total_time == pytest.approx(expected * dt)
    # t_j = j * dt exactly, with no accumulated drift.
    assert schedule.start_times == pytest.approx([j * dt for j in range(expected)])
    assert schedule.is_constant


def test_integer_tick_round_trip_is_exact():
    """Variable intervals survive the tick clock exactly."""
    intervals = [0.25, 0.5, 0.25, 1.0]
    lams = [1.0, 2.0, 1.0, 0.5]
    schedule = schedule_from_intervals(intervals, lams)

    assert schedule.ticks == (250, 500, 250, 1000)
    assert schedule.intervals == pytest.approx(intervals)
    assert not schedule.is_constant
    assert not schedule.is_static_lambda

    # The clock is integer-exact, and t_j equals the running sum of intervals.
    assert schedule.start_ticks == (0, 250, 750, 1000)
    expected_starts = np.cumsum([0.0] + intervals[:-1])
    assert schedule.start_times == pytest.approx(expected_starts)
    assert schedule.total_ticks == 2000
    assert schedule.total_time == pytest.approx(sum(intervals))

    # A finer tick represents finer intervals exactly too.
    fine = schedule_from_intervals([0.0001, 0.00025], [1.0, 1.0], tick=1e-5)
    assert fine.ticks == (10, 25)
    assert fine.intervals == pytest.approx([0.0001, 0.00025])


def test_ticks_for_round_trip():
    for dt, tick in [(0.5, 1e-3), (1.7, 1e-3), (0.00025, 1e-5), (10.0, 1e-3)]:
        n = ticks_for(dt, tick)
        assert n * tick == pytest.approx(dt)


def test_schedule_guards():
    # dt larger than the whole budget.
    with pytest.raises(ValueError, match="exceeds T"):
        constant_schedule(dt=2.0, T=1.0, lam=1.0)
    # dt finer than one tick.
    with pytest.raises(ValueError, match="smaller than one tick"):
        constant_schedule(dt=0.00015, T=1.0, lam=1.0)
    # dt not a whole number of ticks.
    with pytest.raises(ValueError, match="integer multiple"):
        constant_schedule(dt=0.0015, T=1.0, lam=1.0)
    # non-positive inputs.
    with pytest.raises(ValueError):
        constant_schedule(dt=0.1, T=0.0, lam=1.0)
    with pytest.raises(ValueError):
        constant_schedule(dt=-0.1, T=1.0, lam=1.0)
    # malformed schedules.
    with pytest.raises(ValueError, match="same length"):
        RoundSchedule(ticks=(1, 2), lam=(1.0,))
    with pytest.raises(ValueError, match="at least one round"):
        RoundSchedule(ticks=(), lam=())
    with pytest.raises(ValueError, match="positive number of ticks"):
        RoundSchedule(ticks=(1, 0), lam=(1.0, 1.0))
    with pytest.raises(ValueError, match="non-negative"):
        RoundSchedule(ticks=(1,), lam=(-1.0,))


def test_dt_property_rejects_variable_schedules():
    variable = schedule_from_intervals([0.25, 0.5], [1.0, 1.0])
    with pytest.raises(ValueError, match="varying intervals"):
        _ = variable.dt


def test_reference_and_product_paths_agree_on_the_schedule():
    """Both sides derive identical per-round probabilities from one schedule.

    The schedule is shared; the Eqs. 1-3 implementations behind these two calls
    are deliberately separate.
    """
    for dt, T, lam in [(0.5, 10.0, 1.0), (0.1, 1.0, 0.0), (2.0, 20.0, 3.0)]:
        schedule = constant_schedule(dt=dt, T=T, lam=lam)
        for p in (0.0, 0.001, 0.015):
            reference = data_probabilities(schedule, p)
            product = schedule_data_probabilities(schedule, p)
            assert reference == product
            assert len(reference) == schedule.n_rounds

    # Also for a variable schedule, which Stage 4 does not use but must not break.
    variable = schedule_from_intervals([0.25, 0.5, 1.0], [1.0, 2.0, 0.0])
    assert data_probabilities(variable, 0.015) == schedule_data_probabilities(
        variable, 0.015
    )


def test_default_tick_is_documented_value():
    assert DEFAULT_TICK == 1e-3


# ---------------------------------------------------------------------------
# 2-3. Metadata and end-to-end output shape
# ---------------------------------------------------------------------------


@pytest.mark.slow
@pytest.mark.parametrize(("d", "dt", "T"), END_TO_END_CASES)
def test_end_to_end_shape_and_metadata(d, dt, T):
    """Round count, detector count and output lengths for a fixed-dt run."""
    layout = RotatedSurfaceCodeZSector(d)
    schedule = constant_schedule(dt=dt, T=T, lam=1.0)
    params = NoiseParams(p=0.01, lam=1.0, dt=dt, b_read=1.0)
    program = build_memory_program(layout, emit="packed")

    n_shots = 4
    result = program.run(
        p_data=params.p_data,
        p_read=params.p_read,
        dt=schedule.dt,
        n_rounds=schedule.n_rounds,
        base_seed=BASE_SEED,
        shots=n_shots,
    )

    assert schedule.n_rounds == round(T / dt)
    assert len(result.results) == n_shots

    for shot in result.results:
        tags = shot.collate_tags()

        # Metadata round-trips exactly.
        assert tags["dt"] == [pytest.approx(dt)]
        assert tags["p_data"] == [pytest.approx(params.p_data)]
        assert tags["p_read"] == [pytest.approx(params.p_read)]
        assert tags["n_rounds"] == [schedule.n_rounds]
        assert tags["n_checks"] == [layout.n_checks]
        assert tags["n_words"] == [program.n_words]
        assert len(tags["obs"]) == 1

        # One packed block per round plus the final readout block.
        assert len(tags["det"]) == schedule.n_rounds + 1
        assert all(len(block) == program.n_words for block in tags["det"])

    # Unpacked detector array has the documented shape.
    events = program.detection_events(result, schedule.n_rounds)
    assert len(events) == n_shots
    for shot_events in events:
        assert shot_events.shape == (schedule.n_rounds + 1, layout.n_checks)

    # And it matches the reference circuit's detector count.
    circuit = build_memory_circuit(
        layout,
        list(data_probabilities(schedule, 0.01)),
        params.p_read,
        with_coords=False,
    )
    assert circuit.num_detectors == (schedule.n_rounds + 1) * layout.n_checks
    assert circuit.num_observables == 1


@pytest.mark.slow
def test_t_j_reconstruction_matches_the_schedule():
    """t_j = j * dt, reconstructed from the single emitted dt."""
    layout = RotatedSurfaceCodeZSector(3)
    schedule = constant_schedule(dt=0.25, T=2.0, lam=1.0)
    params = NoiseParams(p=0.01, lam=1.0, dt=0.25, b_read=1.0)
    program = build_memory_program(layout)

    result = program.run(
        p_data=params.p_data,
        p_read=params.p_read,
        dt=schedule.dt,
        n_rounds=schedule.n_rounds,
        base_seed=BASE_SEED,
    )
    tags = result.results[0].collate_tags()
    dt_emitted = tags["dt"][0]
    n_rounds_emitted = tags["n_rounds"][0]

    reconstructed = [j * dt_emitted for j in range(n_rounds_emitted)]
    assert reconstructed == pytest.approx(schedule.start_times)


# ---------------------------------------------------------------------------
# 4. Detector ordering and decode plumbing
# ---------------------------------------------------------------------------


@pytest.mark.slow
@pytest.mark.parametrize("d", [3, 5])
def test_detector_ordering_matches_between_paths(d):
    """One known X, same qubit and round in both paths, fires the same indices.

    This is what lets Guppy detection events be fed straight into a decoder
    built from the reference circuit: flattened index ``j * n_checks + c`` must
    mean the same detector on both sides.
    """
    layout = RotatedSurfaceCodeZSector(d)
    n_rounds, inject_round = 4, 2
    program = build_memory_program(layout, emit="packed")

    for target in (0, layout.n_data // 2, layout.n_data - 1):
        coord = layout.data_qubits[target]

        # Guppy path: all stochastic noise off, one deterministic injection.
        result = program.run(
            p_data=0.0,
            p_read=0.0,
            dt=1.0,
            n_rounds=n_rounds,
            base_seed=BASE_SEED,
            inject_round=inject_round,
            inject_qubit=target,
        )
        events = program.detection_events(result, n_rounds)[0]
        guppy_fired = np.flatnonzero(events.ravel())

        # Reference path: the same fault in the same round.
        circuit = single_error_circuit(
            layout, coord, n_rounds=n_rounds, round_index=inject_round
        )
        detectors = circuit.compile_detector_sampler(seed=BASE_SEED).sample(shots=1)
        reference_fired = np.flatnonzero(detectors[0])

        assert guppy_fired.size > 0, "test would be vacuous"
        assert np.array_equal(guppy_fired, reference_fired), (
            f"d={d}, data qubit {coord}: Guppy fired {guppy_fired.tolist()}, "
            f"reference fired {reference_fired.tolist()}"
        )


@pytest.mark.slow
def test_decode_guppy_events_with_weights_from_emitted_metadata():
    """Plumbing: decode a Guppy run using weights derived from its own metadata.

    The decoder is built from a reference circuit parameterised by the
    ``p_data``/``p_read`` the run itself emitted -- nothing is hardcoded. Stim's
    detector error model carries the Eq. A2 weights ``log((1 - p) / p)``.
    """
    layout = RotatedSurfaceCodeZSector(5)
    schedule = constant_schedule(dt=0.5, T=4.0, lam=1.0)
    params = NoiseParams(p=0.008, lam=1.0, dt=0.5, b_read=1.0)
    program = build_memory_program(layout, emit="packed")

    n_shots = 16
    result = program.run(
        p_data=params.p_data,
        p_read=params.p_read,
        dt=schedule.dt,
        n_rounds=schedule.n_rounds,
        base_seed=BASE_SEED,
        shots=n_shots,
    )

    # Everything the decoder needs comes out of the result stream.
    tags = [shot.collate_tags() for shot in result.results]
    p_data_emitted = tags[0]["p_data"][0]
    p_read_emitted = tags[0]["p_read"][0]
    n_rounds_emitted = tags[0]["n_rounds"][0]
    n_checks_emitted = tags[0]["n_checks"][0]

    circuit = build_memory_circuit(
        layout,
        [p_data_emitted] * n_rounds_emitted,
        p_read_emitted,
        with_coords=False,
    )
    decoder = build_decoder(circuit)
    assert not decoder.is_trivial
    assert decoder.num_detectors == (n_rounds_emitted + 1) * n_checks_emitted
    assert decoder.num_observables == 1

    events = program.detection_events(result, n_rounds_emitted)
    detector_batch = np.stack([e.ravel() for e in events])
    assert detector_batch.shape == (n_shots, decoder.num_detectors)

    predictions = decoder.decode_batch(detector_batch)
    observed = np.array([[bool(t["obs"][0])] for t in tags])
    assert predictions.shape == observed.shape == (n_shots, 1)

    # Plumbing only: at p well below threshold with d=5, decoding should
    # succeed most of the time. No statistical claim is made here.
    failures = int(np.count_nonzero(predictions != observed))
    assert failures <= n_shots // 2, (
        f"{failures}/{n_shots} shots mis-decoded -- suspicious for d=5 at "
        f"p_data={p_data_emitted:.4f}"
    )


@pytest.mark.slow
def test_zero_noise_decodes_cleanly_end_to_end():
    """At zero noise: no events, no predicted flips, observable stays 0."""
    layout = RotatedSurfaceCodeZSector(3)
    schedule = constant_schedule(dt=0.5, T=3.0, lam=1.0)
    program = build_memory_program(layout)

    result = program.run(
        p_data=0.0,
        p_read=0.0,
        dt=schedule.dt,
        n_rounds=schedule.n_rounds,
        base_seed=BASE_SEED,
        shots=4,
    )
    events = program.detection_events(result, schedule.n_rounds)
    assert all(not e.any() for e in events)
    assert all(shot.collate_tags()["obs"][0] == 0 for shot in result.results)
