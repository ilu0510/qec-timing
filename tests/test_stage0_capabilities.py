"""Stage 0 capability checks for the Guppy / Selene / Stim toolchain.

One test per item (a)-(e) of the Stage 0 brief, plus a statistical test for the
in-program RNG. Every test uses fixed seeds and is reproducible.

Timings and API notes are recorded in ``docs/interfaces.md``.
"""

import time

import pytest
from guppylang import guppy
from guppylang.std.builtins import array, nat, output, result
from guppylang.std.qsystem.random import RNG
from guppylang.std.qsystem.utils import get_current_shot
from guppylang.std.quantum import cx, measure, measure_array, qubit, x
from guppylang_internals.error import GuppyError
from scipy.stats import binomtest

# Seeds are fixed for every test in this module.
BASE_SEED = 20260927
SELENE_SEED = 1234

# d=27 needs 729 data + 728 ancilla = 1457 qubits; 1500 gives headroom.
N_SCALE_QUBITS = 1500


# --------------------------------------------------------------------------
# (a) Compile and run a Guppy program on Selene with the Stim backend,
#     for many shots, with a fixed seed.
# --------------------------------------------------------------------------


@guppy
def _bell_pair() -> None:
    """X then CNOT: both qubits deterministically measure 1."""
    q0 = qubit()
    q1 = qubit()
    x(q0)
    cx(q0, q1)
    output("m0", measure(q0).read())
    result("m1", measure(q1).read())  # `result` is a deprecated alias for `output`


@pytest.fixture(scope="module")
def bell_emulator():
    return _bell_pair.emulator(n_qubits=2).stabilizer_sim()


@pytest.mark.slow
def test_a_selene_stim_many_shots_fixed_seed(bell_emulator):
    """(a) Stim backend, many shots, fixed seed, reproducible."""
    n_shots = 25
    res = bell_emulator.with_shots(n_shots).with_seed(SELENE_SEED).run()

    assert len(res.results) == n_shots
    for shot in res.results:
        d = shot.as_dict()
        # X then CNOT puts both qubits in |1>.
        assert d["m0"] == 1
        assert d["m1"] == 1

    # Same seed must reproduce the same stream.
    again = bell_emulator.with_shots(n_shots).with_seed(SELENE_SEED).run()
    assert [s.entries for s in res.results] == [s.entries for s in again.results]


# --------------------------------------------------------------------------
# (b) In-program randomness via the platform RNG under the Stim backend.
# --------------------------------------------------------------------------


@guppy
def _rng_draws() -> None:
    """Per-shot seeding: seed = BASE_SEED + shot index."""
    shot = get_current_shot()
    rng = RNG(BASE_SEED + int(shot))
    output("shot", int(shot))
    output("u0", rng.random_float())
    output("u1", rng.random_float())
    output("ib", rng.random_int_bounded(100))
    rng.discard()  # required: RNG is a non-droppable (linear) type


@guppy
def _bernoulli(p: float, n_draws: int, base_seed: int) -> None:
    """Count successes in n_draws Bernoulli(p) trials.

    p, n_draws and base_seed are runtime entrypoint arguments, so both p values
    in the statistical test reuse a single compiled binary.
    """
    shot = get_current_shot()
    rng = RNG(base_seed + int(shot))
    k = 0
    i = 0
    while i < n_draws:
        if rng.random_float() < p:
            k += 1
        i += 1
    rng.discard()
    output("successes", k)
    output("n_draws", n_draws)


@pytest.fixture(scope="module")
def rng_emulator():
    return _rng_draws.emulator(n_qubits=1).stabilizer_sim()


@pytest.fixture(scope="module")
def bernoulli_emulator():
    return _bernoulli.emulator(n_qubits=1).stabilizer_sim()


@pytest.mark.slow
def test_b_rng_works_under_stim(rng_emulator):
    """(b) The platform RNG runs under the Stim backend and is reproducible."""
    res = rng_emulator.with_shots(4).with_seed(SELENE_SEED).run()
    shots = [s.as_dict() for s in res.results]

    assert [s["shot"] for s in shots] == [0, 1, 2, 3]
    for s in shots:
        assert 0.0 <= s["u0"] < 1.0
        assert 0.0 <= s["u1"] < 1.0
        assert 0 <= s["ib"] < 100

    # Explicit RNG seeding makes draws reproducible run-to-run.
    again = rng_emulator.with_shots(4).with_seed(SELENE_SEED).run()
    assert [s.entries for s in res.results] == [s.entries for s in again.results]


@pytest.mark.slow
def test_b_shot0_and_shot1_differ(rng_emulator):
    """(b) Per-shot seeding gives different draws on different shots."""
    res = rng_emulator.with_shots(2).with_seed(SELENE_SEED).run()
    shot0, shot1 = (s.as_dict() for s in res.results)

    assert shot0["shot"] == 0 and shot1["shot"] == 1
    assert shot0["u0"] != shot1["u0"]
    assert shot0["u1"] != shot1["u1"]


def test_b_discard_is_required():
    """(b) Omitting rng.discard() is a compile-time error, not a silent leak."""
    with pytest.raises(GuppyError) as exc:

        @guppy
        def _leaks_rng() -> None:
            rng = RNG(BASE_SEED)
            output("u", rng.random_float())
            # no rng.discard() -> drop violation

        _leaks_rng.compile()

    message = str(exc.value)
    assert "non-droppable" in message or "leaked" in message


@guppy
def _one_rng() -> None:
    r = RNG(1)
    output("a", r.random_float())
    output("b", r.random_float())
    r.discard()


@guppy
def _two_rngs() -> None:
    """Interleave draws from two separately-seeded RNG instances."""
    r1 = RNG(1)
    r2 = RNG(2)
    output("r1_first", r1.random_float())
    output("r2_first", r2.random_float())
    output("r1_second", r1.random_float())
    r1.discard()
    r2.discard()


@pytest.mark.slow
def test_b_rng_instances_are_not_independent():
    """(b) Two RNG instances alias one global runtime PRNG.

    Enforces the CLAUDE.md rule "use exactly ONE RNG instance per program".
    Constructing a second RNG reseeds the shared state in Selene's runtime
    (`selene_random_seed`), so an existing instance's stream is disturbed --
    including its *first* draw, since Guppy orders operations by dataflow
    rather than by source order.
    """
    single = (
        _one_rng.emulator(n_qubits=1).stabilizer_sim().with_shots(1).run()
    ).results[0].as_dict()
    double = (
        _two_rngs.emulator(n_qubits=1).stabilizer_sim().with_shots(1).run()
    ).results[0].as_dict()

    # With independent generators, RNG(1) would yield the same stream in both
    # programs. It does not.
    assert double["r1_first"] != single["a"]
    assert double["r1_second"] != single["b"]


@pytest.mark.slow
@pytest.mark.parametrize("p", [0.01, 0.3])
def test_b_bernoulli_rate_within_999_ci(bernoulli_emulator, p):
    """(b) Empirical Bernoulli(p) rate agrees with p at 99.9% confidence.

    Uses an exact (Clopper-Pearson) interval rather than a normal
    approximation, which is not trustworthy at p = 0.01.
    """
    n_shots, draws_per_shot = 200, 500
    res = (
        bernoulli_emulator.with_shots(n_shots)
        .with_seed(SELENE_SEED)
        .run(p=p, n_draws=draws_per_shot, base_seed=BASE_SEED)
    )

    shots = [s.as_dict() for s in res.results]
    assert len(shots) == n_shots
    successes = sum(s["successes"] for s in shots)
    n_total = sum(s["n_draws"] for s in shots)
    assert n_total == n_shots * draws_per_shot

    lo, hi = binomtest(successes, n_total).proportion_ci(
        confidence_level=0.999, method="exact"
    )
    assert lo <= p <= hi, (
        f"p={p}: observed {successes}/{n_total} = {successes / n_total:.5f}, "
        f"99.9% CI [{lo:.5f}, {hi:.5f}] excludes p"
    )


# --------------------------------------------------------------------------
# (c) Dynamic while-loop driven by a float accumulator, with dt and T
#     supplied from Python as runtime arguments (no recompilation).
# --------------------------------------------------------------------------


@guppy
def _accumulate(dt: float, T: float) -> None:
    """Count rounds of t += dt while t < T. Both dt and T are runtime args."""
    t = 0.0
    n = 0
    while t < T:
        t += dt
        n += 1
    output("n_rounds", n)
    output("t_final", t)


@pytest.fixture(scope="module")
def accumulate_emulator():
    return _accumulate.emulator(n_qubits=1).stabilizer_sim()


@pytest.mark.slow
def test_c_dynamic_loop_runtime_params(accumulate_emulator):
    """(c) dt and T vary from Python against one compiled binary."""
    cases = [(0.5, 2.0, 4), (0.25, 10.0, 40), (1.7, 200.0, 118)]

    for dt, T, expected_rounds in cases:
        res = accumulate_emulator.with_shots(1).with_seed(SELENE_SEED).run(dt=dt, T=T)
        d = res.results[0].as_dict()
        assert d["n_rounds"] == expected_rounds, f"dt={dt}, T={T}"
        assert d["t_final"] >= T

    # Per-shot parameter sweep against the same binary.
    res = (
        accumulate_emulator.with_shots(2)
        .with_seed(SELENE_SEED)
        .run_per_shot([{"dt": 0.5, "T": 2.0}, {"dt": 0.25, "T": 10.0}])
    )
    assert [s.as_dict()["n_rounds"] for s in res.results] == [4, 40]


@pytest.mark.slow
def test_c_float_accumulator_overshoots(accumulate_emulator):
    """(c) Float accumulation can add a round: 0.1 added 10x is < 1.0 exactly.

    Documents a real hazard for the paper's T = sum(dt_j) protocol, where the
    round count must not be inferred from a float accumulator alone.
    """
    res = accumulate_emulator.with_shots(1).with_seed(SELENE_SEED).run(dt=0.1, T=1.0)
    d = res.results[0].as_dict()
    assert d["n_rounds"] == 11, "expected the off-by-one from binary float error"


# --------------------------------------------------------------------------
# (d) Output bits, ints, floats and arrays per shot, and read them back.
# --------------------------------------------------------------------------


@guppy
def _all_result_types() -> None:
    q = qubit()
    x(q)
    output("bit", measure(q).read())  # bit from a measurement
    output("flag", True)  # bool
    output("an_int", -42)  # int
    output("a_nat", nat(7))  # nat
    output("a_float", 3.5)  # float
    output("arr_bool", array(True, False, True))  # array[bool, 3]
    output("arr_int", array(1, 2, 3, 4))  # array[int, 4]
    output("arr_float", array(0.5, 1.5, 2.5))  # array[float, 3]
    output("rep", 11)  # repeated tag ...
    output("rep", 22)  # ... collated, not overwritten


@pytest.fixture(scope="module")
def types_emulator():
    return _all_result_types.emulator(n_qubits=1).stabilizer_sim()


@pytest.mark.slow
def test_d_result_types_per_shot(types_emulator):
    """(d) Every supported value type round-trips per shot."""
    n_shots = 3
    res = types_emulator.with_shots(n_shots).with_seed(SELENE_SEED).run()
    assert len(res.results) == n_shots

    for shot in res.results:
        d = shot.as_dict()
        assert d["bit"] == 1  # bools arrive as int 0/1, not Python bool
        assert d["flag"] == 1
        assert d["an_int"] == -42
        assert d["a_nat"] == 7
        assert d["a_float"] == pytest.approx(3.5)
        assert d["arr_bool"] == [1, 0, 1]  # arrays arrive as Python lists
        assert d["arr_int"] == [1, 2, 3, 4]
        assert d["arr_float"] == pytest.approx([0.5, 1.5, 2.5])

    # entries preserves order and duplicates; as_dict() keeps only the LAST
    # value of a repeated tag, so collate_tags() is required for per-round data.
    first = res.results[0]
    assert first.as_dict()["rep"] == 22
    assert first.collate_tags()["rep"] == [11, 22]
    assert [tag for tag, _ in first.entries].count("rep") == 2


# --------------------------------------------------------------------------
# (e) Scale: ~1500 qubits on the Stim backend (feasibility of d=27).
# --------------------------------------------------------------------------


@guppy
def _wide_circuit() -> None:
    """Allocate N_SCALE_QUBITS qubits, apply a few X/CNOTs, measure all."""
    qs = array(qubit() for _ in range(N_SCALE_QUBITS))
    for i in range(0, N_SCALE_QUBITS, 250):
        x(qs[i])
    for i in range(0, N_SCALE_QUBITS - 1, 250):
        cx(qs[i], qs[i + 1])
    ms = measure_array(qs)
    ones = 0
    for i in range(N_SCALE_QUBITS):
        if ms[i].read():
            ones += 1
    output("ones", ones)


@pytest.mark.slow
def test_e_scale_1500_qubits(capsys):
    """(e) ~1500 qubits runs on Stim; report per-shot time for d=27 feasibility."""
    t0 = time.perf_counter()
    emu = _wide_circuit.emulator(n_qubits=N_SCALE_QUBITS).stabilizer_sim()
    build_s = time.perf_counter() - t0

    n_shots = 50
    t1 = time.perf_counter()
    res = emu.with_shots(n_shots).with_seed(SELENE_SEED).run()
    run_s = time.perf_counter() - t1

    assert len(res.results) == n_shots
    # 6 X gates at i = 0,250,...,1250, each copied onto i+1 by a CNOT.
    for shot in res.results:
        assert shot.as_dict()["ones"] == 12

    per_shot = run_s / n_shots
    with capsys.disabled():
        print(
            f"\n    [e] {N_SCALE_QUBITS} qubits: build {build_s:.2f} s, "
            f"{n_shots} shots in {run_s:.2f} s = {per_shot * 1000:.1f} ms/shot"
        )
    # Generous ceiling: this is a feasibility check, not a benchmark assertion.
    assert per_shot < 1.0, f"{per_shot:.3f} s/shot is too slow for d=27 sweeps"
