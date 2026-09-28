# Stage 2 report: noise injection inside Guppy, static lambda

Noise is now injected inside a Guppy program running on Selene's Stim backend.
Everything is driven by one RNG seeded per shot, the detection events are
computed in-program and emitted bit-packed, and the results agree with the
Stage 1 reference model.

**91 tests pass in 95 s** (11 Stage 0, 62 Stage 1, 18 new). The noise contract
teammates need is in `docs/noise_spec.md`; the wire format is in
`docs/interfaces.md`.

## Headline results

| Question | Answer |
| --- | --- |
| Is `exp()` available inside Guppy? | **No** -- `p_data` is computed in Python and passed in as a runtime arg |
| Fastest wire format | **Bit-packed**, adopted; sparse is never competitive |
| Does emission dominate run time? | **No** -- 3-4%, contradicting a Stage 0 extrapolation |
| Guppy vs reference detection rate | **5/5 cases agree** at 99.9% |
| Fig. 6 minimum vs Eq. 6 | **0.7136 +/- 0.0230** and **0.7195 +/- 0.0201** vs `dt* = 0.7143` |

## 0. `exp()` inside Guppy: not available

Checked against the installed guppylang 1.1.1, not from memory:

- `guppylang.std.maths` **does not exist** (`ModuleNotFoundError`).
- `exp` is **not** in `guppylang.std.builtins` (`pow` and `power` are).

So the brief's fallback applies: **`p_data` is computed in Python from Eqs. 1-2
and passed into the program as a runtime `float`.** No linearisation is used
anywhere, and a test asserts that Eq. 1 differs from `p * lam * dt` by more than
1% where they diverge.

This costs nothing, because Stage 0 established that float entrypoint arguments
bind per run without recompiling -- a whole `dt` sweep reuses one binary.

**Worth knowing for Stage 5/6:** `float ** float` *is* available in Guppy, and
`math.e ** (-z)` reproduces `math.exp(-z)` **bit-for-bit** (`0.6065306597126334`
at `z = 0.5`, `0.1353352832366127` at `z = 2.0`). When an adaptive controller
picks `dt` at runtime and `p_data` must be computed in-program, Eq. 1 can be
evaluated exactly as `1.0 - E ** (-p * lam * dt)`. Recorded in `noise_spec.md`.

## 1. Noise primitives (`src/qec_timing/noise/`)

| Item | Notes |
| --- | --- |
| `bernoulli(rng, p)` | `rng.random_float() < p` |
| `apply_data_noise(rng, data, p_data)` | one draw per data qubit, X if true |
| `flip_record(rng, bit, p_read)` | flips the classical bool only |
| `params.py` | Eqs. 1-3 in Python |

Helpers take the RNG as a **borrowed** parameter (`rng: RNG`, not `@owned`), so
one instance threads through the whole program and is discarded exactly once.
That is required rather than stylistic: Stage 0 established that `RNG(seed)` maps
onto Selene's global runtime PRNG, so a second instance would reseed the first.

**`params.py` deliberately re-implements Eqs. 1-3 rather than importing them
from `reference/noise.py`.** CLAUDE.md makes the reference an independent
cross-check; if the product path imported its physics from the reference, an
error in Eqs. 1-2 would appear identically on both sides and the cross-check
could never see it. `test_product_and_reference_noise_agree` asserts the two
agree to 1e-15 across 36 parameter combinations, which catches drift without
creating the coupling.

## 2. Stub extraction circuit (`src/qec_timing/stub_circuit/`)

Marked as a placeholder in the package docstring, the module docstring and this
report. Per Z-check: allocate an ancilla in |0>, CNOT each data qubit of the
support onto it, measure it. All gates and the ancilla measurement are **ideal**,
which is what makes it equivalent to the reference model's ideal MPP.

Geometry is imported from `reference/layout.py` as instructed, so the two paths
cannot drift. **This makes product code depend on the reference package**, which
inverts the direction CLAUDE.md implies; it is deliberate, and the geometry
should move to a shared module when the real circuit lands.

What is compile-time vs runtime:

| Value | Route |
| --- | --- |
| `p_data`, `p_read`, `n_rounds`, `base_seed`, injection controls | **runtime args**, no recompile |
| `d`, check supports, array sizes, `output` tags | **comptime**, one binary per distance |

Getting the geometry into Guppy needed some care: `py(list)` and comprehensions
over `py()` both failed ("Not compile-time evaluatable"). The working route is
`comptime(FLAT)[i]` indexed inside `for i in range(comptime(N))`, with supports
flattened into a fixed-width table padded with `-1`. Also note `measure(data[i])`
is rejected ("Subscript consumed") -- a qubit cannot be moved out of an array,
so `measure_array` is required.

## 3. Output format: bit-packed, decided on the crossover

Detection events are computed in-program (Eq. A1) and emitted as `n_rounds + 1`
blocks of `ceil(n_checks / 64)` unsigned 64-bit words. Bit/word order, padding
and a copy-paste numpy unpacking snippet are specified in `docs/interfaces.md`.

The decision is a crossover, not a single measurement: sparse costs one word per
*fired detector*, packed costs `ceil(n_checks / 64)` words per *round*.

| d | checks | packed words/round | sparse wins below | observed fired/round |
| --- | --- | --- | --- | --- |
| 11 | 60 | 1 | < 1 per round | 4.8 - 20.8 |
| 27 | 364 | 6 | < 6 per round | 28.9 - 128.8 |

Across the whole Fig. 2 range (`p = 0.015`, `lam = 1`, `dt` from 0.01 to 10) the
observed detection rate sits **5-21x above the crossover at both distances**, so
sparse is never competitive anywhere this project operates. Measured emission
cost against a no-emission control:

| Config | baseline | packed | sparse | dense |
| --- | --- | --- | --- | --- |
| d=27, dt=0.01 | 244 ms/shot | +9.0 ms (3.6%) | +22.8 ms (8.5%) | +26.2 ms (9.7%) |
| d=27, dt=1.0 | 246 ms/shot | +7.9 ms (3.1%) | +27.3 ms (10.0%) | +21.0 ms (7.9%) |
| d=27, dt=10 | 256 ms/shot | +11.6 ms (4.3%) | **+68.5 ms (21.1%)** | +5.7 ms (2.2%) |

Packed's cost is **flat** (7.9-11.6 ms) across a 4.5x change in detection rate
while sparse scales with it, reaching **5.9x packed's cost** at `dt = 10`. Dense
ties packed on wall time but streams 364 entries per round instead of 6, which
matters for storage and for handing data to the decoder. **Packed adopted.**

A round-trip test packs in Guppy, unpacks in Python and compares against a known
injected pattern, at d = 11 (60 checks, one partial word) and d = 15 (112
checks, one full word plus a partial). A separate test asserts all three
emitters describe byte-identical detection events for the same `base_seed`.

### Flagged: a Stage 0 extrapolation does not hold

Stage 0 measured ~40 us per emitted bit and I concluded result streaming
dominates run time. **That does not hold for this program.** Against the
no-emission control, bit-packed emission is only **3-4%** of run time; the
remainder is RNG draws and gate operations (the control alone is 244 ms/shot at
d = 27). The Stage 0 measurement was real but described a different workload --
dumping 1500 loose bits per shot. Corrected in `docs/interfaces.md` and
`docs/noise_spec.md`.

The practical consequence: **further optimising the wire format would be wasted
effort.** If Stage 4 needs more throughput, the RNG call rate is the target --
d = 27 draws 729 data Bernoullis plus 364 readout Bernoullis per round.

One measurement caveat: at d = 11 the emission differences are comparable to
run-to-run variance (one control point came out *slower* than the emitting
formats, which is impossible). The d = 27 rows are stable and carry the
conclusion.

## 4. Tests (`tests/test_stage2_guppy_noise.py`, 18 tests)

| Requirement | Result |
| --- | --- |
| Data-flip rate matches `p_data` | Within 99.9% Clopper-Pearson CI at p = 0.001, 0.05, 0.3 (10^6 trials each) |
| Readout-flip rate matches `p_read` | Within 99.9% CI at p = 0.001, 0.05, 0.3 (10^6 trials each) |
| Different shots give different noise | shot 0 != shot 1 asserted explicitly |
| Same `base_seed` reproduces a run | byte-identical detection events |
| `p = 0` gives zero detection events | all three formats, no detector fires |
| Single injected X fires expected detectors | all 25 data qubits at d = 5, checked against the geometry |
| Ideal extraction leaves logical state undisturbed | `obs = 0` in every shot, 20 shots, 6 rounds |

Two details worth calling out. `p = 0` is tested as *incapable* of firing over
400 000 draws, not merely unlikely. And the single-error test derives the
expected detector set from the layout rather than a golden value, so it cannot
silently agree with a wrong geometry; it also asserts the fault does **not**
re-fire in the following round, confirming the record-flip/state distinction.

## 5. Cross-check against Stage 1

Mean detection-event rate per check per round, Guppy vs reference, 400 shots,
round blocks only (the final readout block is excluded on both sides).

| d | p | lam | dt | b_read | Guppy | Reference | Delta | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 3 | 0.010 | 1.0 | 0.5 | 1.0 | 0.05927 | 0.05792 | +2.34% | OK |
| 3 | 0.015 | 1.0 | 2.0 | 1.0 | 0.13948 | 0.14229 | -1.98% | OK |
| 5 | 0.010 | 1.0 | 0.5 | 1.0 | 0.06461 | 0.06542 | -1.23% | OK |
| 5 | 0.015 | 0.0 | 1.0 | 1.0 | 0.07357 | 0.07401 | -0.60% | OK |
| 5 | 0.005 | 3.0 | 1.0 | 2.0 | 0.07893 | 0.07685 | +2.71% | OK |

**5/5 agree**: the exact 99.9% intervals overlap in every case and the
deviations (-2.0% to +2.7%) are consistent with statistical noise at 400 shots.
The cases span `lam = 0` (pure `p_stab`), `lam = 3`, and `b_read = 2`. No
discrepancy to report; the full sweep remains Stage 4's job.

## 6. Carry-over: the Fig. 6 minimum located properly

Stage 1 could only say the minimum lay on one of the two coarse grid points
bracketing `dt*`. Refining inside that bracket (9 points, spacing factor 1.25
instead of 1.78) and fitting a parabola in `log(dt)` -- which uses every point
rather than just the lowest -- gives an actual vertex with an uncertainty.

| p | fitted minimum | Eq. 6 `dt*` | deviation | sigma |
| --- | --- | --- | --- | --- |
| 0.0035 | **0.7136 +/- 0.0230** | 0.7143 | -0.1% | 0.0 |
| 0.001 | **0.7195 +/- 0.0201** | 0.7143 | +0.7% | 0.3 |

**Eq. 6 is confirmed to within 3%** at both physical error rates. 1.9 min,
~11M shots, reference model only.

Note on method: a first pass with 200 000-shot caps put the p = 0.001 vertex at
0.6255 (-12.4%), which looked like a real deviation but was a low-statistics
fluctuation -- the per-point error bars were 8-9% and the curve is very flat
near the optimum. Raising to 1.2M shots per point (3.5% bars) moved it to
0.7195 and halved the vertex uncertainty. Worth remembering for Stage 4: **near
the optimum, `p_L` varies by less than the error bars over a factor ~2 in `dt`**,
so locating minima needs either a parabola fit or a lot of shots.

## Run times

| Task | Time |
| --- | --- |
| Full test suite (91 tests) | 95 s |
| Format benchmark (d = 11 and 27) | ~6 min |
| Cross-check (5 cases, 400 shots) | ~30 s |
| Minimum refinement (18 points, ~11M shots) | 1.9 min |

Guppy-path cost, 10 rounds: **~62 ms/shot at d = 11**, **~250 ms/shot at
d = 27**, dominated by RNG draws and gates.

## Outputs

```
results/stage2/format_benchmark.csv / _summary.json / .log
results/stage2/crosscheck.json
results/stage2/refine_minimum.csv / _summary.json / .png / .log
docs/noise_spec.md          the standalone noise contract
docs/interfaces.md          + detection-event wire format section
```

## Assumptions and open points

1. **Per-round metadata is deferred to Stage 3.** CLAUDE.md asks for per-round
   `dt, t, p_data, p_read`; under a static schedule these are constant, so they
   are emitted once per shot and the round index is implicit in block ordering.
   Since emission turned out to be cheap, this is a redundancy argument, not a
   cost one. Noted in `noise_spec.md` section 8.
2. **`stub_circuit` imports geometry from `reference`**, inverting the intended
   dependency direction. Deliberate, per the brief; revisit when the real
   circuit lands.
3. **One binary per distance.** Geometry must be comptime, so changing `d`
   recompiles. Cheap (build was under a second at d = 27) and it does not affect
   `dt` sweeps, which vary only runtime arguments.
4. The cross-check compares a **single scalar** (mean detection rate). It would
   not catch an error that preserves the mean while changing correlations --
   for example a wrong check-to-check assignment. Stage 4's full comparison,
   including decoding and logical error rates, is what tests that.
5. `n_qubits` is `d**2 + 1`: ancillas are allocated and measured one at a time,
   which Selene reuses. The real circuit will likely want all ancillas live.
