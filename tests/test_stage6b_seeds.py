"""Stage 6b seed allocation: no two shots, within or across points, share a seed."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

SPEC = importlib.util.spec_from_file_location(
    "stage6b_crosscheck", Path(__file__).resolve().parents[1] / "scripts" / "stage6b_crosscheck.py")
stage6b = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = stage6b  # dataclasses need the module registered
SPEC.loader.exec_module(stage6b)


@pytest.mark.parametrize("shots", [1, 34_332, 1_000_000, stage6b.SEED_STRIDE - 11])
def test_injection_seed_ranges_are_disjoint_across_points(shots):
    """program.py seeds shot k with base_seed + k: those ranges must not overlap."""
    seeds = [stage6b.point_seeds(20261009, i, shots) for i in range(len(stage6b.POINTS))]
    ranges = sorted((s["base_seed"], s["base_seed"] + shots) for s in seeds)
    for (_, hi), (lo, _) in zip(ranges, ranges[1:]):
        assert hi <= lo
    # Selene's own seed and the twin's seeds are distinct from every injection seed.
    for s in seeds:
        for other in seeds:
            lo, hi = other["base_seed"], other["base_seed"] + shots
            assert not lo <= s["selene"] < hi
            assert not lo <= s["twin"] < hi and not lo <= s["twin_dup"] < hi
    assert len({s["selene"] for s in seeds}) == len(seeds)


def test_old_stride_overlapped():
    """Regression note: the 1000-stride allocation used before this fix overlapped."""
    base = [20261009 + 1000 * i + 1 for i in range(len(stage6b.POINTS))]
    assert base[1] < base[0] + 34_332


def test_shots_must_fit_the_block():
    with pytest.raises(ValueError):
        stage6b.point_seeds(1, 0, stage6b.SEED_STRIDE)
    with pytest.raises(ValueError):
        stage6b.point_seeds(1, 0, 0)


def test_duplicate_pairs_by_process():
    r = np.zeros((9, 3), dtype=bool)
    r[[0, 1, 4], 0] = True           # shots 0,1,4 identical: procs 0,1,0
    r[[2, 6], 1] = True              # shots 2,6 identical: procs 2,2
    r[3, 2] = True                   # unique
    same, cross = stage6b.duplicate_pairs_by_process(r, 4)
    assert (same, cross) == (2, 2)   # (0,4), (2,6) same; (0,1), (1,4) cross


@pytest.mark.slow
def test_selene_seeds_by_global_shot_index_across_processes():
    """4 processes reproduce 1 process shot for shot: no shot reuses a seed.

    Covers both the program's injection RNG (base_seed + get_current_shot()) and
    the error model's per-shot seeds, which come from compiled Selene code.
    """
    from guppylang import guppy
    from guppylang.std.angles import angle
    from guppylang.std.builtins import comptime, output
    from guppylang.std.qsystem import phased_x
    from guppylang.std.qsystem.random import RNG
    from guppylang.std.qsystem.utils import get_current_shot
    from guppylang.std.quantum import measure, qubit
    from selene_sim import DepolarizingErrorModel

    @guppy
    def probe(base_seed: int) -> None:
        rng = RNG(base_seed + int(get_current_shot()))
        output("shot", int(get_current_shot()))
        output("r", rng.random_float())
        for _ in range(comptime(16)):
            q = qubit()
            for _ in range(comptime(6)):
                phased_x(q, angle(0.0), angle(0.0))
            output("m", measure(q).read())
        rng.discard()

    def run(n_processes):
        em = DepolarizingErrorModel(p_1q=0.3, p_2q=0.0, p_meas=0.0, p_init=0.0)
        inst = (probe.emulator(n_qubits=1).stabilizer_sim().with_error_model(em)
                .with_shots(400).with_seed(424242).with_n_processes(n_processes))
        return [(t["shot"][0], t["r"][0], tuple(t["m"]))
                for t in (s.collate_tags() for s in inst.run(base_seed=99).results)]

    one, four = run(1), run(4)
    assert [o[0] for o in four] == list(range(400))
    assert one == four
    assert len({o[1] for o in four}) == 400
