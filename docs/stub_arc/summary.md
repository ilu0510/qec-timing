# Noise model: summary for the team

Owner: workstream #2 (error models). Status: **validated**. Details in
`docs/noise_spec.md` (the contract), `docs/interfaces.md` (the wire format),
`docs/performance.md` (costs), and the per-stage reports.

## What the model is

Phenomenological, bit-flip only, Z-check sector, memory in logical |0>. Gates
and stabilizer extraction are **ideal**; every fault is injected explicitly.
Per round of duration `dt` at idling rate `lambda`:

```
p_idle = 1 - exp(-p * lambda * dt)          (Eq. 1)
p_data = 1 - (1 - p) * (1 - p_idle)         (Eq. 2, p_stab = p)
p_read = b_read * p                         (Eq. 3)
```

Three rules that are easy to get wrong and are enforced by tests:

- **One Bernoulli draw per data qubit per round**, at `p_data`. Not two draws
  (one for `p_stab`, one for `p_idle`) -- two flips cancel and give a different
  channel.
- **Read-out noise flips the classical record, never the qubit.** A read-out
  fault must not propagate into the next round's syndrome. Same for the final
  data readout.
- **`p_stab` is never dropped.** At `lambda = 0`, `p_data = p` (to floating-
  point precision; compare with a tolerance, not `==`). It is
  the term that makes frequent measurement costly and creates the `dt`
  trade-off; without it the whole result disappears.

Round counts come from `round(T/dt)` with an integer counter, never a float
accumulator (`dt = 0.1, T = 1.0` gives 10 rounds, not the 11 accumulation
produces).

## How to call it

```python
from qec_timing.layout import RotatedSurfaceCodeZSector
from qec_timing.noise import NoiseParams
from qec_timing.schedule import constant_schedule
from qec_timing.stub_circuit import build_memory_program

layout   = RotatedSurfaceCodeZSector(d=7)
schedule = constant_schedule(dt=0.5, T=10.0, lam=1.0)
params   = NoiseParams(p=0.015, lam=1.0, dt=schedule.dt, b_read=1.0)

program = build_memory_program(layout, emit="packed")   # one binary per d
result  = program.run(
    p_data=params.p_data, p_read=params.p_read, dt=schedule.dt,
    n_rounds=schedule.n_rounds, base_seed=20260927, shots=1000,
)
events = program.detection_events(result, schedule.n_rounds)  # (rounds+1, checks)
```

`p_data`, `p_read`, `dt`, `n_rounds` and `base_seed` are **runtime** arguments,
so a whole sweep reuses one compiled binary. Only `d` is compile-time.

**Two performance rules.** Guppy pays **~9.3 s per `run()` call**, independent of
shot count, so size shots up front and issue **one call per configuration** --
never loop `run()` to accumulate statistics. And per-round cost is unstable to
~3x between runs, so bound work by a time budget rather than a shot count.

## Output schema

Per shot, read with `collate_tags()` (**never `as_dict()`** -- it keeps only the
last value of a repeated tag): `det` (detection events, `n_rounds + 1`
bit-packed blocks), `obs` (logical Z parity), plus `dt`, `p_data`, `p_read`,
`n_rounds`, `n_checks`, `n_words`.

Detection events are **bit-packed into unsigned 64-bit words**, LSB first,
ascending word order, final word zero-padded. **Full schema, bit/word order, and
a copy-paste numpy unpacking snippet are in `docs/interfaces.md`** -- start
there. Note `t_j = j * dt` is only valid while `dt` is constant.

Decoders should depend on the `Decoder` protocol in
`qec_timing.reference.decoder`, not on pymatching directly. Weights come from
the emitted metadata: a detector error model built at the run's own
`p_data`/`p_read` carries the Eq. A2 weights `log((1-p)/p)`.

## Two paths, and why

- **Guppy/Selene is the validation path.** It is the artefact the physics is
  defined by, and it is ~200x more expensive per round because every draw and
  gate goes through the emulator.
- **The pure-Python reference (Stim + pymatching) is the production path** for
  sweeps.

This is a deliberate division of labour. It is sound because the two are
verified to produce the *same* detector stream: identical detector ordering
(asserted by injecting one X at a known qubit and round in both paths) and
agreeing logical error rates. **If the reference changes, re-run the Stage 4a
cross-check before trusting a sweep.**

## Validated results

- **Stage 1** -- Appendix C thresholds reproduced: `p_th` = 0.0353 / 0.0286 /
  0.0230 for `b_read` = 0.5 / 1 / 2, against the paper's 0.0360 / 0.0288 /
  0.0231. The 0.5-2% shortfall is finite size (we stop at d = 19, the paper
  reaches 35); refitting on larger-d subsets moves it monotonically onto the
  paper values.
- **Stage 2** -- noise injection in Guppy: data-flip and read-out-flip rates
  match `p_data` and `p_read` inside 99.9% binomial intervals over 10^6 trials,
  at p from 0.001 to 0.3.
- **Stage 3** -- detector ordering identical between the two paths; a run
  decodes end-to-end with weights taken from its own emitted metadata.
- **Stage 4a** -- **Guppy and the reference agree on `p_L` at 18/18 points**
  across the full `dt` range at d = 5, 7, 9. Each d shows the expected U-shape.
  `1/dt* = (0.376 +/- 0.026) d - (0.821 +/- 0.189)`, against the paper's slope
  0.402 -- agreement at **1.0 sigma**.

**One caveat worth knowing:** the measured optima sit 10-30% *above* Eq. 6
evaluated at the paper's `g = 0.8`. This was tested, not assumed: repeating the
scan at 3x the total time (which removes any finite-rounds artefact, since the
predicted optimum is T-independent) left d = 7 unmoved at 0.1 sigma. The cause
is that the **effective `g` is ~0.70 at d = 5-9**, rising towards 0.8 as d
grows -- the paper fits `g = 0.8` against d = 11-27. If you need `dt*` at small
d, use the local `g`, not 0.8.

## Not covered

- **No gate noise.** Gates and ancilla measurement are ideal by construction.
  Circuit-level noise is a different model, not a parameter change.
- **No burst / time-varying `lambda`** (Eq. 9). `lambda` is static.
- **No adaptive timing** (Appendix F controller). `dt` is fixed within a run.
- **No X-check sector.** Only Z-checks are simulated; the other sector is
  equivalent under X <-> Z.
- **Not validated above d = 9 on the Guppy path**, and not at the paper's
  T = 200. The reference path handles those; Guppy does not, by design.
- `stub_circuit/` is a **placeholder** for workstream #1's extraction circuit,
  chosen to be provably equivalent to an ideal MPP. Replacing it should not
  disturb the noise model, but re-run the Stage 4a cross-check when it lands.
