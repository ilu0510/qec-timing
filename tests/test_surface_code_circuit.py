"""Full-code preparation/extraction, patch isolation and existing decoder seam."""

import itertools

import numpy as np
import pytest

from qec_timing.circuit import build_surface_code_program
from qec_timing.circuit.program import _preparation_corrections
from qec_timing.layout import RotatedSurfaceCode, RotatedSurfaceCodeZSector
from qec_timing.reference.circuit import build_memory_circuit, single_error_circuit
from qec_timing.reference.decoder import build_decoder
from qec_timing.reference.sampler import clopper_pearson
from qec_timing.stub_circuit import build_memory_program
from qec_timing.wire import unpack_packed_words

SEED = 20261001


def check_matrix(layout, checks):
    return np.array([[q in s for q in layout.data_qubits] for s in checks], dtype=int)


@pytest.mark.parametrize("d", [3, 5, 9, 27])
def test_full_css_geometry_and_preparation_inverse(d):
    layout = RotatedSurfaceCode(d)
    assert layout.z_checks == RotatedSurfaceCodeZSector(d).z_checks
    hx, hz = check_matrix(layout, layout.x_checks), check_matrix(layout, layout.z_checks)
    assert hx.shape == hz.shape == ((d*d-1)//2, d*d)
    assert not ((hx @ hz.T) % 2).any()
    lx = np.array([q in layout.logical_x_support for q in layout.data_qubits], dtype=int)
    lz = np.array([q in layout.logical_z_support for q in layout.data_qubits], dtype=int)
    assert not ((hz @ lx) % 2).any()
    assert not ((hx @ lz) % 2).any()
    assert lx.sum() == lz.sum() == d and lx @ lz == 1
    for checks, matrix in [(layout.x_checks, hx), (layout.z_checks, hz)]:
        corrections = np.array(_preparation_corrections(layout, checks), dtype=int).reshape(
            layout.n_checks, layout.n_data)
        assert np.array_equal((matrix @ corrections.T) % 2, np.eye(layout.n_checks))


def test_both_sectors_distance_three_exhaustively():
    layout = RotatedSurfaceCode(3)
    for checks, logical in [(layout.z_checks, layout.logical_z_support),
                            (layout.x_checks, layout.logical_x_support)]:
        matrix = check_matrix(layout, checks)
        observable = np.array([q in logical for q in layout.data_qubits], dtype=int)
        weights = []
        for bits in itertools.product([0, 1], repeat=layout.n_data):
            fault = np.array(bits)
            if not ((matrix @ fault) % 2).any() and fault @ observable % 2:
                weights.append(fault.sum())
        assert min(weights) == 3


@pytest.mark.parametrize("kwargs", [{"d": 2}, {"d": 3.5}, {"d": True},
                                   {"d": 3, "q_d": 0}, {"d": 3, "basis": "Y"}])
def test_invalid_build_arguments(kwargs):
    with pytest.raises(ValueError):
        build_surface_code_program(**kwargs)


@pytest.mark.slow
@pytest.mark.parametrize("basis", ["Z", "X"])
def test_preparation_and_extraction_of_independent_patches(basis):
    program = build_surface_code_program(3, q_d=2, basis=basis)
    assert program.n_qubits == 19
    result = program.run(p_data=0., p_read=0., n_rounds=3, base_seed=SEED, shots=12)
    assert np.asarray(program.detection_events(result, 3)).shape == (12, 2, 4, 4)
    assert not np.asarray(program.detection_events(result, 3)).any()
    assert not program.observables(result).any()
    for shot in result.results:
        tags = shot.collate_tags()
        assert tags["q_d"] == [2]
        assert len(tags["complementary_syn"]) == 3 * 2
        # This verifies normalization of the randomly projected prep syndromes.
        assert not unpack_packed_words(tags["complementary_syn"], 4).any()


@pytest.mark.slow
@pytest.mark.parametrize("basis", ["Z", "X"])
def test_known_fault_detector_order_and_patch_isolation(basis):
    program = build_surface_code_program(5, q_d=2, basis=basis)
    layout = program.layout
    target = layout.data_index[(0, 0)]
    result = program.run(p_data=0., p_read=0., n_rounds=3, base_seed=SEED, shots=4,
                         inject_round=1, inject_qubit=target, inject_patch=1)
    checks = layout.z_checks if basis == "Z" else layout.x_checks
    expected = np.zeros((4, layout.n_checks), dtype=bool)
    expected[1] = [(0, 0) in support for support in checks]
    for events in program.detection_events(result, 3):
        assert not events[0].any()
        assert np.array_equal(events[1], expected)
    assert np.array_equal(program.observables(result), [[False, True]] * 4)
    if basis == "Z":
        reference = single_error_circuit(layout, (0, 0), n_rounds=3, round_index=1)
        assert np.array_equal(reference.compile_detector_sampler().sample(shots=1)[0],
                              expected.ravel())


@pytest.mark.slow
def test_noisy_full_extraction_preserves_validated_stub_stream_and_decoder():
    full = build_surface_code_program(3)
    stub = build_memory_program(full.layout)
    args = dict(p_data=0.09, p_read=0.04, n_rounds=4, base_seed=SEED, shots=512)
    result = full.run(**args)
    legacy = stub.run(**args)
    assert np.array_equal(full.detection_events(result, 4), stub.detection_events(legacy, 4))
    assert np.array_equal(full.observables(result),
                          [s.collate_tags()["obs"] for s in legacy.results])
    for shot in result.results:
        assert not unpack_packed_words(shot.collate_tags()["complementary_syn"], 4).any()
    # Independent reference physics and matching weights are unchanged.
    reference = build_memory_circuit(full.layout, [0.09]*4, 0.04)
    decoder = build_decoder(reference)
    failures = np.count_nonzero(full.decode(result, 4, decoder) != full.observables(result))
    low, high = clopper_pearson(int(failures), 512, confidence=0.999)
    det, obs = reference.compile_detector_sampler(seed=SEED).sample(
        shots=100_000, separate_observables=True)
    rate = np.mean(decoder.decode_batch(det) != obs)
    assert low <= rate <= high


@pytest.mark.slow
@pytest.mark.parametrize("emit", ["sparse", "dense"])
def test_multi_patch_wire_formats(emit):
    program = build_surface_code_program(3, q_d=2, emit=emit)
    result = program.run(p_data=0., p_read=0., n_rounds=2, base_seed=SEED,
                         inject_round=0, inject_qubit=0, inject_patch=1)
    expected = np.zeros((2, 3, 4), dtype=bool)
    expected[1, 0] = [(0, 0) in s for s in program.layout.z_checks]
    assert np.array_equal(program.detection_events(result, 2)[0], expected)
    decoder = build_decoder(build_memory_circuit(program.layout, [0.01]*2, 0.01))
    assert np.array_equal(program.decode(result, 2, decoder), program.observables(result))
