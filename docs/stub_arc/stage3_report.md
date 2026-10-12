# Stage 3 report: time bookkeeping and per-round metadata

Shared geometry and scheduling modules, per-round metadata on the wire, and an
end-to-end Guppy run that decodes through the Stage 1 decoder seam.

**110 tests pass in 128 s** (11 Stage 0, 62 Stage 1, 18 Stage 2, 19 new). No
benchmarks, no sweeps.

## 0. Geometry moved to `src/qec_timing/layout.py`

`RotatedSurfaceCodeZSector` now sits alongside `ansatz.py` as a shared module.
Clean move, no compatibility shim: imports updated in `reference/`,
`stub_circuit/`, both existing test files and four scripts. All 91 pre-existing
tests still pass unchanged.

This removes the wrinkle flagged in Stage 2, where product code imported
geometry from the cross-check package.

## 1. `src/qec_timing/schedule.py`

`RoundSchedule`, `constant_schedule`, `schedule_from_intervals`, `ticks_for`,
promoted out of `reference/noise.py`.

**The shared module contains no physics.** `RoundSchedule` previously had a
`data_probabilities(p)` method that called the reference's `p_data`. Moving the
class into shared code with that method attached would have quietly reunified
Eqs. 1-3, which Stage 2 deliberately implements twice so the cross-check has
something independent to check. So the schedule now carries ticks, `lam` and the
clock only, and each side converts:

- `reference.noise.data_probabilities(schedule, p)`
- `noise.params.schedule_data_probabilities(schedule, p)`

A test asserts the two agree on both constant and variable schedules. The
split is recorded in a table in `docs/interfaces.md`: geometry, scheduling and
the ansatz are shared; the noise equations are not.

| Behaviour | Result |
| --- | --- |
| Fixed `dt` | `n_rounds = round(T / dt)`, integer counter, `t_j = j * dt` |
| Variable `dt` | integer tick counter, default tick `1e-3` |
| `dt > T` | raises ("a single round would overrun the total time budget") |
| `round(T/dt) == 0` | raises (unreachable once `dt <= T`; kept explicit) |
| `dt` below one tick | raises ("smaller than one tick") |
| `dt` not a tick multiple | raises rather than silently snapping |

Refusing to snap matters: a snapped interval would make the emitted `dt`
disagree with the one actually simulated.

The variable-interval path is **built but not used** beyond its round-trip test,
as instructed -- `[0.25, 0.5, 0.25, 1.0]` becomes ticks `(250, 500, 250, 1000)`
with `t_j = (0, 0.25, 0.75, 1.0)` and total 2000 ticks, all exact. Burst `lam(t)`
and the adaptive controller remain out of scope per the roadmap.

The `dt = 0.1, T = 1.0` case is tested against the failure it exists to prevent:
the schedule gives 10 rounds, and the test also runs the float accumulator
in-line and asserts it gives 11.

## 2. Per-round metadata

`dt`, `p_data` and `p_read` are emitted **once per shot**, since a static `lam`
and fixed `dt` make them identical in every round. `t_j = j * dt` is
reconstructed. Read back with `collate_tags()`.

The schema in `docs/interfaces.md` states explicitly that **the `t_j = j * dt`
reconstruction is valid only while `dt` is constant**, and that a stage which
varies the interval must emit `dt`, `t_j`, `p_data` and `p_read` per round
instead -- a schema change, not a variant. `RoundSchedule.is_constant` is
exactly that precondition on the Python side.

## 3. End-to-end at fixed `dt`

Three `(d, dt, T)` combinations: `(3, 0.5, 3.0)`, `(3, 0.25, 2.0)`,
`(5, 1.0, 5.0)`. Each asserts:

- `n_rounds == round(T / dt)`
- `n_rounds + 1` packed blocks of `ceil(n_checks / 64)` words each
- exactly one `obs`, and one each of `dt`, `p_data`, `p_read`, `n_rounds`,
  `n_checks`, `n_words`
- emitted `dt`/`p_data`/`p_read` equal the values passed in
- the unpacked array is `(n_rounds + 1, n_checks)`
- the reference circuit built from the same schedule reports the same detector
  count, `(n_rounds + 1) * n_checks`, and one observable

## 4. Decoding through the Stage 1 seam

A d = 5 run, 16 shots, decoded with **everything the decoder needs taken from
the result stream**: the reference circuit is rebuilt at the emitted `p_data`
and `p_read`, and Stim's detector error model supplies the Eq. A2 weights
`log((1 - p_e) / p_e)`. Nothing is hardcoded. The decoder reports the expected
`num_detectors` and `num_observables`, `decode_batch` returns `(n_shots, 1)`,
and the zero-noise case produces no events, no predicted flips and `obs = 0`.

No statistical claim is made; this is a plumbing check.

### Detector ordering is now asserted directly

The check that makes the above sound: **one known X at a known qubit and round,
injected into both paths, must fire the identical detector indices.** Tested at
d = 3 and d = 5, at three data-qubit positions each (first, middle, last), with
all stochastic noise off. The Guppy detection-event array flattened row-major
and the reference circuit's detector vector agree exactly.

That is what licenses feeding Guppy detection events straight into a decoder
built from the reference circuit: flattened index `j * n_checks + c` denotes the
same detector on both sides, round-major with the final readout block last. A
disagreement here would have produced a decoder that silently mis-corrects
rather than an error, so it is worth pinning.

## CLAUDE.md corrections

Three statements in the "Selene/Guppy facts" section were superseded by Stage 2
measurements and have been corrected, as agreed:

1. *"Result streaming dominates run time (~40 us/bit)"* -> emission is **3-4%**
   of run time once bit-packed, measured against a no-emission control. The
   40 us/bit figure described dumping 1500 loose bits per shot, a different
   workload. RNG draws and gate operations dominate.
2. *"preferably as sparse indices"* -> **bit-packed**; sparse runs 5-21x above
   the crossover at d = 11 and d = 27 and costs 5.9x packed at d = 27, dt = 10.
3. The citation pointed at `docs/stage0_report.md`, which does not exist; it now
   points at `docs/interfaces.md`. No Stage 0 report was backfilled.

## Files

```
src/qec_timing/layout.py         moved from reference/layout.py
src/qec_timing/schedule.py       new, shared, no physics
src/qec_timing/reference/noise.py   Eqs. 1-3 only; re-exports the schedule
src/qec_timing/noise/params.py      + schedule_data_probabilities
src/qec_timing/stub_circuit/program.py  + dt argument, + metadata emission
tests/test_stage3_schedule_metadata.py  19 tests
docs/interfaces.md               + per-shot output schema, + shared/duplicated table
```

## Open points

1. **Per-round metadata is still per-shot.** Correct for a static schedule and
   documented as such, but Stage 4's `dt` sweep varies `dt` *between* runs, not
   within one, so this remains valid there. Only a within-run varying `dt`
   forces the change.
2. `docs/stage1_report.md` and `docs/stage2_report.md` still refer to
   `reference/layout.py`. Left alone deliberately: they are dated records of
   what was done at the time. The living documents -- `noise_spec.md`,
   `interfaces.md`, `CLAUDE.md` -- are current.
3. The decode check uses 16 shots and asserts only that most shots decode. Real
   logical-error-rate comparison against the reference is Stage 4's job.
