# Interfaces

Toolchain facts established in Stage 0. Everything here was verified against the
installed package source in `.venv/Lib/site-packages/` and by running it, not
from documentation or memory. Tests live in `tests/test_stage0_capabilities.py`.

## Package versions

Python **3.14.0**, win_amd64. Full pins in `requirements.txt`.

| Package | Version | Notes |
| --- | --- | --- |
| guppylang | 1.1.1 | requires Python `>=3.12,<4` |
| guppylang-internals | 1.1.1 | |
| selene-sim | 0.3.2 | |
| selene-core | 0.3.2 | |
| selene-hugr-qis-compiler | 0.5.0 | `cp310-abi3` wheel |
| stim | 1.16.0 | not called directly yet; Stim is reached via the Selene plugin |
| PyMatching | 2.4.0 | |
| numpy | 2.5.3 | |
| scipy | 1.18.1 | |
| pytest | 9.1.1 | |
| pytket | 2.18.4 | transitive, `cp312-abi3` |
| hugr | 0.18.6 | transitive |
| ziglang | 0.16.0 | transitive: Selene's bundled C toolchain/linker |

**Python 3.11 cannot be used** — `guppylang` requires `>=3.12`. Only 3.11 and
3.14 were available on this machine, so 3.14 was chosen; all compiled
dependencies resolve (`abi3` or real `cp314` wheels) and the full suite passes.

**Do not trust `selene_sim.__version__`** — it reports `0.2.10` while the
installed distribution is `0.3.2` (stale in-package version string). Use
`importlib.metadata.version("selene-sim")`.

## (a) Running a Guppy program on Selene with the Stim backend

The high-level path is `guppylang.emulator`, reached via `.emulator()` on a
`@guppy` function. There is no need to touch `selene_sim.build` directly.

```python
from guppylang import guppy
from guppylang.std.builtins import output
from guppylang.std.quantum import qubit, measure, x, cx

@guppy
def prog() -> None:
    q = qubit()
    x(q)
    output("m", measure(q).read())

res = prog.emulator(n_qubits=2).stabilizer_sim().with_shots(1000).with_seed(1234).run()
```

- `.emulator(n_qubits: int, platform="helios")` -> `EmulatorInstance`.
  `n_qubits` is required (it cannot be inferred, since a program may request
  qubits at runtime) unless the function carries
  `@guppylang.decorator.expected_qubits(n)`.
- `.stabilizer_sim()` selects the **Stim** backend. The default is Quest
  (statevector). `.coinflip_sim()` is a no-quantum backend, and
  `.with_simulator(plugin)` accepts any `selene_core.Simulator`.
- Configuration is an immutable chain of `with_*` methods, each returning a new
  instance: `.with_shots(n)` (default **1**), `.with_seed(n)`,
  `.with_shot_offset(n)`, `.with_n_processes(n)`, `.with_error_model(m)`,
  `.with_trace()`, `.with_metrics()`.
- `.run(**entrypoint_args)` -> `EmulatorResult`.

Underneath, `selene_sim.SeleneInstance.run_shots` is what actually runs
(`selene_sim/instance.py:211`):

```python
def run_shots(self, simulator, n_qubits, n_shots=1, error_model=IdealErrorModel(),
              runtime=SimpleRuntime(), event_hook=NoEventHook(), verbose=False,
              timeout=None, results_logfile=None, random_seed=None, shot_offset=0,
              shot_increment=1, n_processes=1, parse_results=True,
              seed_mode="default") -> Iterator[Iterator[TaggedResult] | ...]
```

`seed_mode="default"` (the default) uses `random_seed` to seed an RNG that then
generates each shot's seed; `seed_mode="legacy"` increments the seed by 1 per
shot. The Stim plugin is `selene_stim_plugin.StimPlugin`, exported as
`selene_sim.Stim`, with one parameter `angle_threshold: float = 1e-4`.

## (b) In-program randomness

Works under the Stim backend. Both import paths in the brief are correct:

```python
from guppylang.std.qsystem.random import RNG
from guppylang.std.qsystem.utils import get_current_shot

@guppy
def prog() -> None:
    shot = get_current_shot()           # -> int
    rng = RNG(BASE_SEED + int(shot))    # per-shot seeding
    u = rng.random_float()              # -> float in [0, 1)
    rng.discard()                       # REQUIRED
```

`RNG` methods (`guppylang/std/qsystem/random.py`):

| Method | Signature | Returns |
| --- | --- | --- |
| constructor | `RNG(seed: int)` | `RNG` |
| `random_float` | `(self) -> float` | uniform `[0, 1)` |
| `random_int` | `(self) -> int` | random 32-bit signed int |
| `random_int_bounded` | `(self, bound: int) -> int` | `[0, bound)`, `bound < 2**31` |
| `random_angle` | `(self) -> angle` | `[-pi, pi)` |
| `random_clifford_angle` | `(self) -> angle` | multiple of `pi/2` |
| `random_advance` | `(self, delta: int) -> None` | jump the stream forward/back |
| `shuffle` | `(self, array) -> None` | in-place Fisher-Yates |
| `discard` | `(self @owned) -> None` | consumes the RNG |

Also available: `make_discrete_distribution(weights: array[float, N])` ->
`DiscreteDistribution[N]` with `.sample(rng) -> int`, which is a ready-made
weighted sampler.

**Why it works on Stim.** These are HUGR external ops in
`QSYSTEM_RANDOM_EXTENSION` / `QSYSTEM_UTILS_EXTENSION`, implemented by Selene's
**core runtime** (`selene.dll`: `selene_random_f64`, `selene_random_u32_bounded`,
`selene_random_seed`, `selene_random_advance`, `selene_get_current_shot`), not by
the simulator plugin. Randomness is therefore backend-independent.

Verified behaviour:

- **`discard()` is required**, and enforced at compile time, not at runtime.
  `RNG` is declared `@custom_type(RNGCONTEXT_T, copyable=False, droppable=False)`,
  so omitting it raises `guppylang_internals.error.GuppyError`:
  *"Variable `rng` with non-droppable type `RNG` is leaked"*.
- **Different shots produce different draws** when seeded as
  `base_seed + get_current_shot()`.
- **Draws are independent of `.with_seed()`.** When the RNG is explicitly
  seeded, `.with_seed(1)` and `.with_seed(99999)` give byte-identical draws.
  `.with_seed()` governs the simulator/error-model streams, not an explicitly
  seeded `RNG`. Monte Carlo reproducibility therefore rests on `base_seed`
  alone.
- **`.with_shot_offset(k)` shifts `get_current_shot()`** to `k, k+1, ...`, which
  gives fresh non-overlapping randomness for extending a run.
- **RNG instances are not independent — use exactly one per program.** `RNG(seed)`
  maps onto Selene's *global* runtime PRNG (`selene_random_seed` in
  `selene.dll`), so constructing a second instance reseeds the state the first
  one is drawing from. Measured: `RNG(1)` alone yields
  `0.30447083548642695, 0.896538217086345`, but with an `RNG(2)` also
  constructed, the same `RNG(1)` yields `0.30427282815799117` then
  `0.4489013897255063`. Note that even the **first** draw changes, because Guppy
  orders operations by dataflow rather than source order, so both constructors
  run before either draw. Pass one RNG through all functions. Covered by
  `test_b_rng_instances_are_not_independent`.

The pure-Guppy fallback was not needed. It exists as
`guppylang.std.random.seeded_pcg32(seed: nat) -> PCG32` (PCG32 XSH-RR 64/32),
whose state is a local Guppy value rather than platform-global state, with
`.next_int()` / `.next_int_bounded(bound)`.

## (c) Dynamic loops and runtime parameters

A `while` loop on a float accumulator compiles and runs, and **both `dt` and `T`
can be passed from Python without recompiling**:

```python
@guppy
def accumulate(dt: float, T: float) -> None:
    t = 0.0
    n = 0
    while t < T:
        t += dt
        n += 1
    output("n_rounds", n)

emu = accumulate.emulator(n_qubits=1).stabilizer_sim()
emu.with_shots(1).run(dt=0.25, T=10.0)      # 40 rounds
emu.with_shots(1).run(dt=1.7,  T=200.0)     # 118 rounds, same binary
```

What can and cannot vary without a rebuild:

| Mechanism | Recompile? | Applies to |
| --- | --- | --- |
| **Entrypoint arguments** — `run(**kwargs)` | **No** | `int`, `float`, `bool`, and `array` of those |
| `run_per_shot([{...}, {...}])` | **No** | a different argument mapping per shot; shot count is inferred |
| `get_current_shot()` | No | per-shot values derived in-program |
| Python value capture / `comptime` | **Yes** | anything closed over from Python at decoration time |
| Array lengths, loop bounds of `for ... in range(N)`, `output` tags | **Yes** | these are compile-time (`tag: str @comptime`) |

Entrypoint arguments are emulator-only — the docstring notes this "capability is
not supported when submitting programs to hardware". That is acceptable here
since this project only ever emulates.

**Numerical hazard.** The round count follows binary float accumulation, not
exact arithmetic: `dt=0.1, T=1.0` gives **11 rounds, not 10**, because 0.1 added
ten times is `0.9999999999999999 < 1.0`. The paper's protocol accumulates
`T = sum_j dt_j` (Appendix A), so the number of rounds must be computed
explicitly rather than inferred from a float accumulator, or results will shift
by a whole round at particular `(dt, T)` pairs. Covered by
`test_c_float_accumulator_overshoots`.

## (d) Results

`result` **is a deprecated alias for `output`** — `guppylang/std/platform.py:112`
is literally `result = output`, "deprecated since guppylang v1.0". New code
should use `output`.

```python
def output(tag: str, value) -> None
```

`tag` must be a string **literal** (`tag: str @comptime`). Supported overloads,
as reported by the compiler:

```
output(tag: str @comptime, value: int)              output(tag, value: array[int, n])
output(tag: str @comptime, value: nat)              output(tag, value: array[nat, n])
output(tag: str @comptime, value: bool)             output(tag, value: array[bool, n])
output(tag: str @comptime, value: float)            output(tag, value: array[float, n])
```

A `Measurement` **cannot** be passed directly — call `.read()` on it first,
which blocks until the outcome is available. `guppylang.std.debug.state_output`
handles quantum state output separately.

Reading back per shot: `.run()` returns an `EmulatorResult` whose `.results` is a
list of `QsysShot`, one per shot, in shot order.

| Accessor | Behaviour |
| --- | --- |
| `shot.entries` | list of `(tag, value)` pairs, order and duplicates preserved |
| `shot.as_dict()` | `dict`; **silently keeps only the LAST value of a repeated tag** |
| `shot.collate_tags()` | `dict[str, list]`; gathers repeated tags — use this for per-round data |
| `shot.to_register_bits()` | bitstring view |
| `res.collated_counts()`, `res.register_bitstrings()`, `res.metrics()`, `res.traces()`, `res.circuits()`, `res.partial_states()` | aggregate / analysis views |

Marshalling details that matter:

- `bool` arrives as Python **`int` 0/1**, not `True`/`False`.
- `array[...]` arrives as a Python **`list`**.
- `nat` arrives as `int`.
- `as_dict()` on `output("d", 11); output("d", 22); output("d", 33)` yields
  `{"d": 33}`. Use `collate_tags()` -> `{"d": [11, 22, 33]}` when emitting one
  value per syndrome round.

## (e) Scale — d=27 feasibility

1500 qubits (d=27 needs 729 data + 728 ancilla = **1457**) on the Stim backend,
with a few X/CNOT gates and all qubits measured. Measured on this machine:

| Configuration | Per-shot |
| --- | --- |
| 1500 bits emitted per shot, 1 process | ~60 ms |
| 1 int emitted per shot, 1 process | ~9.5 ms |
| 1 int emitted per shot, `n_processes=4` | ~3.3 ms |

**d=27 is feasible.** Two findings shape how Stage 1 should be written:

1. **Result streaming dominates simulation cost, by roughly 6x.** Emitting 1500
   bits per shot costs ~60 ms/shot; emitting a single reduced integer from the
   identical circuit costs ~9.5 ms. Per-round detector data should be reduced
   in-program, or emitted deliberately, rather than dumped per qubit per round.
2. **`with_n_processes(4)` gives a further ~3x** (9.5 -> 3.3 ms/shot). Beyond 4
   the gain flattens on this machine.

Build (compile + link) cost is separate from run cost and is **not** per shot:

- **First build in a fresh environment: ~69 s** (one-off LLVM/zig toolchain
  warm-up).
- Every subsequent build: **~1.0-3.5 s**.

Since build cost is paid per program and entrypoint arguments avoid rebuilding
for each sweep point, parameter sweeps should reuse one `EmulatorInstance` and
vary `dt`/`T` through `run()` / `run_per_shot()`.

---

# Decoder interface (Stage 1)

Added in Stage 1 with your approval. The production matching decoder is owned by
teammate #3; the reference model decodes against the contract below so that
implementation can be swapped in without touching the sampler. Defined in
`src/qec_timing/reference/decoder.py`.

```python
@runtime_checkable
class Decoder(Protocol):
    @property
    def num_detectors(self) -> int: ...
    @property
    def num_observables(self) -> int: ...
    def decode_batch(
        self, detectors: NDArray[np.bool_]
    ) -> NDArray[np.bool_]: ...
```

- `decode_batch` takes a `(n_shots, num_detectors)` boolean array and returns a
  `(n_shots, num_observables)` boolean array of **predicted observable flips**.
  A shot has failed when the prediction differs from the sampled observable
  (Appendix A: "a shot fails if the correction returned by matching differs from
  the accumulated physical error by a nontrivial logical operator").
- Implementations are constructed from a `stim.DetectorErrorModel`. The DEM
  already carries the Eq. A2 log-likelihood weights
  `w_e = log((1 - p_e) / p_e)` derived from the circuit's error probabilities,
  so **a decoder never needs the noise parameters passed separately.**
- **Believed vs. true noise.** Appendix A specifies that for fixed intervals the
  decoder knows the true `lam(t)`, while for the adaptive controller the decoder
  belief follows the controller state. Both are expressed the same way: build
  the DEM from a circuit carrying the *believed* probabilities, and sample from
  a circuit carrying the *true* ones. The protocol needs no extra surface for
  this; Stages 5-6 will use it.

Reference implementation: `PyMatchingDecoder`, built by
`build_decoder(circuit, decompose_errors=False)`.

- `decompose_errors=False` is correct for this model: every X fault flips at
  most **two** Z-checks (each data qubit lies in at most two Z-checks), so the
  DEM is already graphlike and needs no decomposition. Verified by the layout
  tests.
- **Degenerate case.** At `p = 0` the DEM contains no error mechanisms and
  pymatching cannot be constructed from it. `PyMatchingDecoder` detects this
  (`.is_trivial`) and predicts no flips, which is the correct answer when no
  error can occur. A teammate implementation should handle the same case, or
  document that it requires a non-empty DEM.

## Detector and observable conventions

Fixed by `reference/circuit.py` and relied on by the tests:

| Quantity | Value |
| --- | --- |
| Z-checks per round | `(d**2 - 1) // 2` |
| Detectors | `(n_rounds + 1) * n_checks` |
| Observables | 1 (logical Z) |
| Measurements | `n_rounds * n_checks + d**2` |

- Detector ordering is round-major, then check index within the round; the final
  block of `n_checks` detectors comes from the data readout.
- Round 0 detectors reference a single measurement (the Z-stabilizers are
  deterministic on logical |0>); later rounds compare consecutive rounds
  (Eq. A1, `D_{j,c} = s_{j,c} XOR s_{j-1,c}`).
- The logical Z observable is the parity of the final data readout over **row 0**.
  In this layout the logical X operator is a *column* of X (a column commutes
  with every Z-check; a row does not), and a column meets row 0 exactly once.
- Detectors carry `(row, col, round)` coordinates when
  `with_coords=True` (the default), for debugging.

---

# Detection-event wire format (Stage 2)

The Guppy program computes detection events in-program (Eq. A1) and emits them
**bit-packed into 64-bit words**, one block per round. This is the format the
decoding side should expect. Chosen by measurement; see the benchmark below.

## Layout

Per shot the program emits, under the tag `det`, **`n_rounds + 1` blocks** of
`n_words = ceil(n_checks / 64)` words each:

- blocks `0 .. n_rounds-1` are the syndrome rounds,
- block `n_rounds` is the final data-readout block.

Within a block:

| Rule | Value |
| --- | --- |
| **Bit order** | LSB first. Bit `b` of word `w` is check index `w * 64 + b`. |
| **Word order** | Ascending: word 0 holds checks 0-63, word 1 holds 64-127, ... |
| **Padding** | The final word is **zero-padded in its high bits**. Padding bits are always 0 and carry no meaning. |
| **Word type** | Guppy `nat`, i.e. **unsigned** 64-bit. Values above `2**63` occur and are correct. |
| **Check index** | Position in `RotatedSurfaceCodeZSector.z_checks` |

The program also emits `n_rounds`, `n_checks` and `n_words` once per shot, so a
consumer can reshape the stream without knowing the distance in advance.

## Unpacking (copy this)

```python
import numpy as np

tags = shot.collate_tags()                      # one QsysShot
words = np.asarray(tags["det"], dtype=np.uint64)
words = words.reshape(n_rounds + 1, n_words)

shifts = np.arange(64, dtype=np.uint64)
bits = ((words[:, :, None] >> shifts) & np.uint64(1)).astype(bool)
events = bits.reshape(words.shape[0], -1)[:, :n_checks]
# events[j, c] is True when detector (round j, check c) fired
```

`qec_timing.stub_circuit.unpack_packed_words` is this snippet, and
`MemoryProgram.detection_events(result, n_rounds)` does it for a whole run.
Note the explicit shift-and-mask rather than `np.unpackbits`, which would depend
on host byte order.

**A shot in which nothing ever fires still emits every block** (all-zero words),
so the stream length is fixed. That is not true of the sparse format, where the
`det` tag can be absent entirely.

## Why bit-packed: the benchmark

The choice is a crossover, not a preference. Sparse costs one word per *fired
detector*; packed costs `ceil(n_checks / 64)` words per *round* regardless of
what fires. So sparse only wins below the crossover:

| d | checks | packed words/round | sparse wins below | observed fired/round | verdict |
| --- | --- | --- | --- | --- | --- |
| 11 | 60 | 1 | < 1 per round | 4.8 - 20.8 | never competitive |
| 27 | 364 | 6 | < 6 per round | 28.9 - 128.8 | never competitive |

Measured at the Fig. 2 parameters (`p = 0.015`, `lam = 1`, `dt` from 0.01 to 10,
10 rounds), the observed detection rate sits **5-21x above the crossover** at
both distances, so sparse never wins anywhere in the range the project actually
uses. Emission cost, isolated against a no-emission control:

| Config | baseline (no emit) | packed | sparse | dense |
| --- | --- | --- | --- | --- |
| d=27, dt=0.01 | 244 ms/shot | +9.0 ms (3.6%) | +22.8 ms (8.5%) | +26.2 ms (9.7%) |
| d=27, dt=1.0 | 246 ms/shot | +7.9 ms (3.1%) | +27.3 ms (10.0%) | +21.0 ms (7.9%) |
| d=27, dt=10 | 256 ms/shot | +11.6 ms (4.3%) | **+68.5 ms (21.1%)** | +5.7 ms (2.2%) |

Packed's cost is **flat** (7.9-11.6 ms) across a 4.5x change in detection rate,
while sparse scales with it and reaches 5.9x packed's cost at `dt = 10`. Dense
ties packed on wall time but streams 364 entries per round instead of 6, which
matters for storage and for handing data to the decoder.

**Correction to a Stage 0 extrapolation.** Stage 0 measured ~40 us per emitted
bit and concluded result streaming dominates run time. That held for dumping
1500 loose bits per shot; it does **not** hold here. With bit-packing, emission
is only **3-4% of run time** and the cost is dominated by RNG draws and gate
operations in the extraction circuit (the no-emission control is already
244 ms/shot at d=27). Optimising the wire format further would be wasted effort;
the RNG call rate is where the time goes.

---

# Per-shot output schema (Stage 3, final)

What one shot of the Guppy memory program emits. Read with
`shot.collate_tags()` -- **never `as_dict()`**, which keeps only the last value
of a repeated tag and so would discard all but the final detector block.

| Tag | Type | Count per shot | Meaning |
| --- | --- | --- | --- |
| `det` | list of `n_words` unsigned 64-bit words | `n_rounds + 1` | detection events, one block per round plus the final readout block |
| `obs` | int 0/1 | 1 | logical Z parity of the final data readout |
| `dt` | float | 1 | round duration |
| `p_data` | float | 1 | Eq. 2 data-fault probability actually used |
| `p_read` | float | 1 | Eq. 3 readout-flip probability actually used |
| `n_rounds` | int | 1 | number of syndrome rounds |
| `n_checks` | int | 1 | `(d**2 - 1) // 2` |
| `n_words` | int | 1 | `ceil(n_checks / 64)` |

Bit/word order and padding for `det` are in the *Detection-event wire format*
section above. Booleans arrive as `int` 0/1, not Python `bool`.

## Reconstructing the clock

```python
tags = shot.collate_tags()
dt, n_rounds = tags["dt"][0], tags["n_rounds"][0]
t = [j * dt for j in range(n_rounds)]            # constant dt ONLY
```

> **`t_j = j * dt` is valid only while `dt` is constant.** Stage 3 runs a fixed
> interval, so `dt`, `p_data` and `p_read` are identical in every round and are
> emitted **once per shot** rather than once per round. If a future stage varies
> the interval -- a burst `lam(t)` or an adaptive controller choosing `dt` at
> runtime -- this reconstruction silently becomes wrong, and the program must
> emit `dt`, `t_j`, `p_data` and `p_read` **per round** instead. Anything
> consuming this schema should check `n_rounds` against the number of `det`
> blocks and treat a per-round `dt` as a schema change, not a variant.

`qec_timing.schedule.RoundSchedule` carries the same information on the Python
side and exposes `is_constant`, which is exactly the precondition above;
`start_times` gives `t_j` for the variable case without float drift.

## Decoding a run

Weights come from the emitted metadata, never from hardcoded constants: build a
reference circuit at the emitted `p_data`/`p_read`, and Stim's detector error
model carries the Eq. A2 weights `log((1 - p_e) / p_e)`.

```python
from qec_timing.reference.circuit import build_memory_circuit
from qec_timing.reference.decoder import build_decoder

tags = [shot.collate_tags() for shot in result.results]
circuit = build_memory_circuit(
    layout,
    [tags[0]["p_data"][0]] * tags[0]["n_rounds"][0],
    tags[0]["p_read"][0],
    with_coords=False,
)
decoder = build_decoder(circuit)

events = program.detection_events(result, tags[0]["n_rounds"][0])
predictions = decoder.decode_batch(np.stack([e.ravel() for e in events]))
observed = np.array([[bool(t["obs"][0])] for t in tags])
failures = np.count_nonzero(predictions != observed)
```

**Detector index `j * n_checks + c` means the same detector on both paths** --
the Guppy program and the reference Stim circuit order detectors round-major,
check-index within a round, with the final readout block last. That is asserted
directly: `test_detector_ordering_matches_between_paths` injects one X at a
known qubit and round in both paths and requires the identical detector indices
to fire, at d = 3 and d = 5.

## Shared vs. duplicated modules

| Module | Status | Why |
| --- | --- | --- |
| `qec_timing/layout.py` | **shared** | geometry is combinatorial; both paths must agree exactly |
| `qec_timing/schedule.py` | **shared** | time bookkeeping, no physics |
| `qec_timing/ansatz.py` | **shared** | analytic, imports neither path |
| `noise/params.py` vs `reference/noise.py` | **duplicated on purpose** | Eqs. 1-3 implemented twice so the cross-check is independent; a test asserts they agree |

`schedule.py` holds no physics for this reason: turning a schedule into
probabilities is done by `reference.noise.data_probabilities` or
`noise.params.schedule_data_probabilities`, one per side.
