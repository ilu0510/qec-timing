# Stage 6b report: the Stim twin reproduces Selene's native noise

Scope: show that the Stage 6a circuit-level twin reproduces Selene with its own
native noise model switched on. This plays the same role for the circuit-level
model that the Stage 4a cross-check played for the phenomenological one. No
sweeps beyond the 12 cross-check points.

**Result: 12/12 points pass at 99%, on both the logical error rate and the
per-round detection-event rates.** Pooled over all points, Selene's failure count
is 0.996 +/- 0.029 times the twin's prediction. Selene time was 96 s.

Run: `results/stage6b/20261009T150628Z/` (`points.csv`, `detection_rates.csv`,
`summary.json`, `run.log`), from commit `f39abd4`. Versions: selene-sim 0.3.2,
guppylang 1.1.1, stim 1.16.0, PyMatching 2.4.0, numpy 2.5.3.

## 1. What was verified before running (and what changed)

**Selene's error model.** The noise comes from `selene_sim.DepolarizingErrorModel(p_1q, p_2q, p_meas, p_init)`, the only bundled model with these
parameters. It runs on the Stim backend (`.stabilizer_sim()`) at d = 3 and 5.
Its own counters (`gates_1q`, `errors_1q_x/y/z`, `gates_2q`, `errors_2q_*`)
show:

- every `Rxy` and every `Rz` is charged `p_1q`, with X, Y and Z equally likely: `DEPOLARIZE1`;
- every `Rzz` is charged `p_2q`, spread over the 15 two-qubit Paulis: `DEPOLARIZE2`;
- `p_meas` and `p_init` are charged once per measurement and once per reset.

**Native op counts were measured, not assumed.** I read them from the error
model's input (`CircuitExtractor.get_error_model_input()`) on one noiseless
traced shot of `program.py`:

| source op | compiled | charged single-qubit ops |
| --- | --- | --- |
| `h` | `Rxy`, `Rz` | 2 |
| `cx(c, t)` | `Rxy(t)`; `Rzz(c, t)`; `Rz(c)`; `Rxy(t)`, `Rz(t)` | target 1 before + 2 after; control 0 before + 1 after |
| `x`, `z` | `Rxy` / `Rz` | 1 |
| `phased_x(q, 0, 0)` | `Rxy` | 1 |

**Stage 6a's CNOT assumption was wrong per qubit**, and has been fixed (commit
`3b57dfd`). 6a placed one single-qubit op before and one after the `Rzz` on
both qubits, taking notebook 04c's "2 per data qubit" as the rule. That figure
is the average over the control and target roles. The ancilla gets 1 op as a
control and 3 as a target, not 2. The total per CNOT (4) was right.

**Checked gap by gap.** For every qubit, I compared the number of charged
single-qubit ops between consecutive reset / two-qubit / measure events. This
covered d = 3 and 5, both memory bases and both CNOT orders. Event sequences are
identical, and every gap matches except one: the preparation->round-0 gap on
qubits that receive preparation corrections. That gap is the one place the
twin models an expectation rather than a fixed count (next section).

**Explicit twin option `mirror_guppy_program=True`** (requires
`native_charging`; off by default). It charges the extra ops that `program.py`'s
noise injection and preparation create:

- the per-round `h` wrapper in the X memory (2 ops per `h`);
- the 1-op charge on each injected `x` that fires (exact, as `PAULI_CHANNEL_1`);
- the conditional preparation corrections. Each runs with probability 1/2, so I model it by its per-qubit average, `DEPOLARIZE1(compose(p1/2, c_i))`. This drops an O(p1^2) correlation between qubits in one correction.

**Variant A idle** needed `program.py` to accept idle ops. I added a
compile-time `idle_ops=0` keyword (commit `4e836ff`, approved). Tests show the
default and `idle_ops=0` produce identical shots for the same seed, and that
Selene's `gates_1q` counter grows by exactly `idle_ops x n_data x n_rounds`
(the compiler keeps every op).

## 2. Configuration

Noise (approved): `p_2q = 5e-3`, `p_1q = 1e-4` (ratio 50), `p_meas = 2e-3`,
`p_init = 0`, symmetric readout. T = 12, CNOT order `hook_safe` everywhere.

| | variant A (`depolarizing_kops`) | variant B (`dephasing_h2`) |
| --- | --- | --- |
| memory | Z | X |
| Selene idle | `idle_ops = k dt` charged `phased_x(q, 0, 0)` per data qubit per round, k = 100 | `program.py` injection with `p_data = p_z(dt)`, inside the `h` wrapper |
| twin idle | `DEPOLARIZE1(compose(p_1q, k dt))` | `p_z(dt) = sin^2(0.5 * 0.1 * dt)`, f_lin = 0 (mirrored injection) |
| dt | 1, 4, 12 (100 / 400 / 1200 idle ops per round) | 1, 4, 6 |

The idle noise carries real weight at these settings; without it, the twin's
p_L is much lower. Variant A at d = 3: 0.003 -> 0.092 at dt = 12 and 0.032 -> 0.082 at dt = 1. Variant B, dt = 6: dephasing is the dominant source.

Method:
- **Decoder:** both sides are decoded by the same decoder, built from the twin's detector error model. Only the noise processes differ.
- **Shots:** one Selene `run()` per point, sized for about 100 expected failures, 4 processes. The twin sampled at least 2,000 failures per point.
- **Pass criteria, at 99%:**
  - p_L: two-proportion |z| < 2.576.
  - Detection events: for every round, the mean number of detection events per shot agrees by Welch z, Bonferroni-corrected across the round's n_rounds + 1 rounds.
- **Seeding sanity check:** the number of repeated detection patterns among Selene shots matches that in the same number of twin shots at every point (for example 214 vs 201, 0 vs 0). So the 4 processes do not reuse seeds.

## 3. Results

| d | variant | order | dt | k | rounds | Selene p_L [99% CI] (fail/shots) | twin p_L [99% CI] | z (p_L) | det max abs z / crit |
|---|---|---|---|---|---|---|---|---|---|
| 3 | A depolarizing_kops (Z) | hook_safe | 1 | 100 | 12 | 0.0778 [0.0592, 0.0999] (93/1195) | 0.0837 [0.0802, 0.0873] | -0.72 | 1.62 / 3.36 |
| 3 | A depolarizing_kops (Z) | hook_safe | 4 | 100 | 3 | 0.0824 [0.0640, 0.1039] (108/1311) | 0.0763 [0.0729, 0.0798] | +0.81 | 2.16 / 3.02 |
| 3 | A depolarizing_kops (Z) | hook_safe | 12 | 100 | 1 | 0.1001 [0.0780, 0.1258] (109/1089) | 0.0919 [0.0882, 0.0956] | +0.93 | 1.15 / 2.81 |
| 5 | A depolarizing_kops (Z) | hook_safe | 1 | 100 | 12 | 0.0562 [0.0434, 0.0713] (106/1885) | 0.0531 [0.0502, 0.0560] | +0.60 | 2.85 / 3.36 |
| 5 | A depolarizing_kops (Z) | hook_safe | 4 | 100 | 3 | 0.0569 [0.0434, 0.0730] (97/1704) | 0.0587 [0.0557, 0.0618] | -0.31 | 1.61 / 3.02 |
| 5 | A depolarizing_kops (Z) | hook_safe | 12 | 100 | 1 | 0.0995 [0.0776, 0.1251] (109/1095) | 0.0914 [0.0877, 0.0951] | +0.93 | 0.96 / 2.81 |
| 3 | B dephasing_h2 (X) | hook_safe | 1 | - | 12 | 0.0499 [0.0377, 0.0646] (91/1823) | 0.0549 [0.0520, 0.0579] | -0.91 | 1.75 / 3.36 |
| 3 | B dephasing_h2 (X) | hook_safe | 4 | - | 3 | 0.1162 [0.0889, 0.1481] (94/809) | 0.1237 [0.1178, 0.1299] | -0.64 | 1.26 / 3.02 |
| 3 | B dephasing_h2 (X) | hook_safe | 6 | - | 2 | 0.2238 [0.1767, 0.2767] (107/478) | 0.2092 [0.2019, 0.2168] | +0.78 | 1.01 / 2.94 |
| 5 | B dephasing_h2 (X) | hook_safe | 1 | - | 12 | 0.0264 [0.0199, 0.0343] (91/3447) | 0.0290 [0.0275, 0.0306] | -0.90 | 1.36 / 3.36 |
| 5 | B dephasing_h2 (X) | hook_safe | 4 | - | 3 | 0.1174 [0.0913, 0.1478] (105/894) | 0.1119 [0.1062, 0.1178] | +0.51 | 0.77 / 3.02 |
| 5 | B dephasing_h2 (X) | hook_safe | 6 | - | 2 | 0.2000 [0.1527, 0.2541] (86/430) | 0.2330 [0.2253, 0.2408] | -1.60 | 0.56 / 2.94 |

**Aggregate checks** (not part of the per-point criterion, but they guard
against a consistent bias that each point is too small to resolve):

| check | value | expectation if the models agree |
| --- | --- | --- |
| p_L z-scores, 12 points | sum of z^2 = 8.83 (chi^2, 12 dof, p = 0.72); mean z = -0.04 | chi^2 ~ 12, mean 0 |
| by variant | A: Stouffer z = +0.91 (p = 0.76); B: Stouffer z = -1.13 (p = 0.48) | |
| pooled failures | Selene 1196 vs twin-expected 1200.6: ratio **0.996 +/- 0.029** | 1 |
| detection rates, 78 rounds | sum of z^2 = 69.0 (78 dof, p = 0.76); Selene/twin mean ratio 1.006 (sd 0.050); largest abs z = 2.85, below its 3.36 threshold | |

So the z-scores look like N(0, 1) draws, with no sign pattern by variant,
distance or dt.

**Resolution.** At about 100 failures per point, each point resolves the twin
to roughly +/-20-30% (99% interval half-width). The pooled ratio bounds any
common bias to about +/-7.5% at 99%. That is enough to rule out the kind of
error 6b was designed to catch: the old CNOT split, a missing `h` wrapper, or
uncharged idle ops each shift the per-round charge by tens of percent. It is not
yet a precision match. The 1000-failure rerun (section 5) would cut the pooled
bound to about +/-2.5%.

## 4. Measured Selene cost

The smoke-era cost model (8 us per op, single process) was far too high. I
refit it to the 12 timed `run()` calls. With 4 processes, every point lands
within about 25% of:

```
run() wall time = shots * (0.181 us * ops + 0.90 ms)
```

In single-process terms (about 3x speed-up at 4 processes, per
`interfaces.md`) that is **0.54 us per charged native op plus 2.7 ms per shot**,
plus **3.4 s per point** outside `run()` (build and decode). The Stage 4a
figure of 9.3 s per `run()` call did not appear here. The per-shot term is RNG
draws, emission and shot setup, and it dominates the small variant B circuits.
That is why their "per op" figure in `points.csv` (`measured_us_per_op`, which
divides everything by ops) runs 1-5.6 us, against 0.57-0.75 us for variant A.
Variant A's cost is almost entirely the idle ops: 11-36k charged ops per shot,
of which k T n_data = 10.8k (d = 3) or 30k (d = 5) are `phased_x`.

The script now uses these values (`US_PER_OP_DEFAULT = 0.54`,
`US_PER_SHOT_DEFAULT = 2700`, `SECONDS_PER_RUN = 3.4`). The ~3x run-to-run
instability noted in CLAUDE.md is handled at run time. Each run's
measured/predicted ratio raises a slowdown factor, which then resizes the
remaining points' shots to fit the budget.

## 5. Not covered, and next

- **Not exercised:**
  - asymmetric readout (Selene's `p_meas` is symmetric);
  - `p_init`, which was 0;
  - the `row_major` order, which was dropped by decision;
  - non-integer `k dt`;
  - d > 5.
  These are tested on the twin alone in 6a.
- **Rerun at 1000 failures per point** (same 12 points, 30 min budget; planned at about 10 min of Selene time):
  `.venv\Scripts\python.exe scripts\stage6b_crosscheck.py --target-failures 1000`.
  The twin's sample now defaults to 10x the Selene target (at least 2,000), so
  its error stays negligible.
- Stage 6c (the dt\* scaling study under circuit-level noise) has not been started.
