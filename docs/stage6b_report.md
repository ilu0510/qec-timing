# Stage 6b report: the Stim twin reproduces Selene's native noise

Scope: show that the Stage 6a circuit-level twin reproduces Selene with its own
native noise model switched on. This plays the same role for the circuit-level
model that the Stage 4a cross-check played for the phenomenological one. No
sweeps beyond the 12 cross-check points.

**Final result (section 9): the twin reproduces Selene at all 12 points.**
The aggregates use only the 13 independently seeded measurements, all at about
1000 failures or more:

- p_L: chi^2 = 13.81 on 13 degrees of freedom (p = 0.39).
- Pooled Selene/twin failure ratio: **0.985 +/- 0.008** (variant A 0.996 +/- 0.013, variant B 0.978 +/- 0.011).
- Detection-event rates: chi^2 = 74.2 on 81 (p = 0.69).

One point (B, d = 5, dt = 6) failed once, at z = -2.95. Its pre-registered
3000-failure repeat in a fresh seed block gave z = -0.82, so the -2.95 is
recorded as a fluctuation (section 8). Variant B shows a mild, non-significant
tendency for Selene to sit below the twin. It is carried forward as a bound,
not a correction (section 9).

**Stage 6b is closed.**

Runs, in order:

| run | what |
| --- | --- |
| `20261009T150628Z` | 100 failures/point, 12 points |
| `20261009T154022Z` | 1000 failures/point, 12 points |
| `20261009T160001Z` | B d=5 dt=1, fixed seeds |
| `20261009T161404Z` | the other 5 B points, fixed seeds |
| `20261009T161923Z` | B d=5 dt=6 at 3000 failures, fresh seed block |

The first run (`20261009T150628Z`, from commit `f39abd4`) is described in
sections 1-5. Each run directory holds `points.csv`, `detection_rates.csv`,
`summary.json` and `run.log`. Versions: selene-sim 0.3.2,
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
- The 1000-failure rerun is in section 6, and the seed investigation it prompted is in section 7.
- Stage 6c (the dt\* scaling study under circuit-level noise) has not been started.

## 6. Rerun at 1000 failures per point

Run: `results/stage6b/20261009T154022Z/`, from commit `68643c0`
(`--target-failures 1000`; the twin sampled at least 10,000 failures per point).
**12/12 pass.** Selene time was 689 s (627 s inside `run()`, the rest build and
decode), against a planned 607 s.

| d | variant | order | dt | k | rounds | Selene p_L [99% CI] (fail/shots) | twin p_L [99% CI] | z (p_L) | det max abs z / crit | Selene s |
|---|---|---|---|---|---|---|---|---|---|---|
| 3 | A depolarizing_kops (Z) | hook_safe | 1 | 100 | 12 | 0.0838 [0.0775, 0.0905] (1025/12231) | 0.0818 [0.0799, 0.0837] | +0.79 | 2.65 / 3.36 | 41 |
| 3 | A depolarizing_kops (Z) | hook_safe | 4 | 100 | 3 | 0.0738 [0.0681, 0.0798] (983/13324) | 0.0751 [0.0733, 0.0769] | -0.54 | 1.41 / 3.02 | 34 |
| 3 | A depolarizing_kops (Z) | hook_safe | 12 | 100 | 1 | 0.0909 [0.0840, 0.0982] (998/10980) | 0.0911 [0.0889, 0.0932] | -0.06 | 2.15 / 2.81 | 26 |
| 5 | A depolarizing_kops (Z) | hook_safe | 1 | 100 | 12 | 0.0526 [0.0485, 0.0569] (999/19010) | 0.0526 [0.0513, 0.0539] | -0.03 | 2.40 / 3.36 | 160 |
| 5 | A depolarizing_kops (Z) | hook_safe | 4 | 100 | 3 | 0.0589 [0.0543, 0.0638] (988/16768) | 0.0596 [0.0582, 0.0611] | -0.38 | 0.94 / 3.02 | 109 |
| 5 | A depolarizing_kops (Z) | hook_safe | 12 | 100 | 1 | 0.0890 [0.0822, 0.0962] (985/11068) | 0.0904 [0.0882, 0.0925] | -0.48 | 1.34 / 2.81 | 68 |
| 3 | B dephasing_h2 (X) | hook_safe | 1 | - | 12 | 0.0505 [0.0464, 0.0547] (955/18929) | 0.0528 [0.0515, 0.0541] | -1.40 | 2.23 / 3.36 | 26 |
| 3 | B dephasing_h2 (X) | hook_safe | 4 | - | 3 | 0.1236 [0.1144, 0.1333] (1010/8171) | 0.1224 [0.1197, 0.1251] | +0.32 | 2.46 / 3.02 | 4 |
| 3 | B dephasing_h2 (X) | hook_safe | 6 | - | 2 | 0.2118 [0.1967, 0.2274] (1012/4779) | 0.2093 [0.2050, 0.2136] | +0.40 | 1.00 / 2.94 | 2 |
| 5 | B dephasing_h2 (X) | hook_safe | 1 | - | 12 | 0.0293 [0.0270, 0.0318] (1007/34332) | 0.0291 [0.0284, 0.0299] | +0.21 | 1.89 / 3.36 | 140 |
| 5 | B dephasing_h2 (X) | hook_safe | 4 | - | 3 | 0.1153 [0.1068, 0.1243] (1030/8931) | 0.1120 [0.1094, 0.1146] | +0.96 | 1.04 / 3.02 | 11 |
| 5 | B dephasing_h2 (X) | hook_safe | 6 | - | 2 | 0.2332 [0.2168, 0.2503] (1002/4296) | 0.2328 [0.2284, 0.2373] | +0.06 | 1.01 / 2.94 | 4 |
| 5 | B, rerun with fixed seeds (section 7) | hook_safe | 1 | - | 12 | 0.0279 [0.0257, 0.0303] (954/34150) | 0.0293 [0.0286, 0.0300] | -1.41 | 1.66 / 3.36 | 416 |

| aggregate | as recorded | with B d=5 dt=1 replaced by its rerun |
| --- | --- | --- |
| p_L z, sum of z^2 (12 dof) | 4.49 (p = 0.97) | 6.44 (p = 0.89) |
| mean z | -0.01 | -0.15 |
| Stouffer z, A / B | -0.28 / +0.23 | -0.28 / -0.43 |
| pooled Selene / twin failures | **0.9994 +/- 0.0091** (11,994 vs 12,001) | 0.9950 +/- 0.0091 |
| detection rates, 78 rounds | sum of z^2 = 85.7 (p = 0.26); largest abs z 2.65 (threshold 3.36); mean Selene/twin ratio 0.9993 | |

At 1000 failures each point resolves the twin to roughly +/-8% (99% interval
half-width). The pooled ratio now bounds any common bias to **about +/-2.3% at
99%**. The low chi^2 as recorded (P(chi^2 <= 4.49) = 0.03) is not significant
after the rerun substitution (0.11). Section 7 explains why the B points of
this run are also not fully independent of each other.

**Cost.** 689 s against the planned 607 s (+13%) confirms the refit cost
model. The isolated rerun of B d=5 dt=1 then took 416 s against 140 s for
the identical configuration in this run. That is the ~3x run-to-run
instability noted in CLAUDE.md, seen directly. The script's slowdown factor
exists for this reason.

## 7. Duplicate-pattern anomaly at B d=5 dt=1, and the seed allocation

**Observation.** In the 1000-failure run, B d=5 dt=1 had 680 repeated
detection patterns among Selene's shots against 287 in the twin. Every other
point matched. It was also the largest run (34,332 shots).

**Within a point, no shot reuses a seed.** This was checked three ways:

1. *Shot partitioning* (installed `selene_sim/instance.py`, `run_shots`):
   process j of 4 runs the global shots k with k = j (mod 4) (offset j,
   increment 4), so the four shot sets are disjoint. `program.py` seeds its
   injection RNG with `base_seed + k`. For this point that is base_seed
   20270010 (old run) or 110261010 (rerun), giving seeds {base + j, base + j + 4, ...}
   for process j: four disjoint residue classes.
2. *Seeding is by global shot index* (new slow test
   `test_selene_seeds_by_global_shot_index_across_processes`). A 4-process run
   reproduces the 1-process run exactly, shot for shot. That covers the
   injection RNG draws and the error model's noise, whose per-shot seeds come
   from compiled Selene code. Reported shot indices are exactly 0..N-1.
3. *Full records split by process* (rerun, all rounds + observable, shots with
   at least one detection event): identical-record pairs split **376 same-process
   / 1,177 cross-process**. Random coincidences over 4 interleaved processes
   should split 1 : 3, i.e. 388 : 1,165. Shared seeds between processes would
   show up as an excess of cross-process pairs; there is none.

**Where the extra duplicates came from: the comparison, not Selene.** The
twin's duplicate count was taken from its first sampling batch, which is capped
at 20,000 shots. So 34,332 Selene shots were compared with 20,000 twin shots.
Repeated patterns grow roughly with the square of the shot count (a birthday
count), so the scaled expectation is about 287 x (34,332 / 20,000)^2, roughly
850. This was the only point with more than 20,000 shots, which is why only it
looked anomalous. The script now draws an independent twin sample of exactly
the Selene shot count. On the rerun: **Selene 735, twin 719** at 34,150 shots
each.

**A real seed problem found on the way: overlap ACROSS points.** Seeds were
allocated as `seed + 1000 i` for point i, with `base_seed = seed + 1000 i + 1`.
Every point had more than 1,000 shots, so the injection-seed ranges of
different points overlapped. All 66 pairs overlapped in the 1000-failure run;
for example, B d=3 dt=1 and B d=5 dt=1 shared 15,929 seeds. Consequences:

- Variant A is unaffected. It runs with `p_data = p_read = 0`, so its RNG draws never inject anything.
- The error model's seeds were always distinct per point (`with_seed(seed + 1000 i)`).
- Each point's own pass/fail verdict is valid, since its shots are all distinct.
- The 15 B-B pairs shared injection streams. Their injected-dephasing faults were therefore correlated between points: for the same d, a dt=1 fault (u < p_z(1)) implies a dt=4 fault (u < p_z(4)) on paired shots. The cross-point aggregates in sections 3 and 6 assumed independent B points, so their B parts are slightly optimistic. The verdict does not depend on them.
- It also contradicts CLAUDE.md's rule of distinct base_seeds per batch.

**Fix** (`scripts/stage6b_crosscheck.py`): point i now owns the seed block
`[seed + i * 10^7, seed + (i + 1) * 10^7)`. The Selene seed sits at the block
start, the injection seeds run from start + 1 + shot, and the twin seeds sit at
the block end. A point whose shots would not fit is refused. `--only LABEL`
reruns a single point within its own block. Tests are in
`tests/test_stage6b_seeds.py` (8 passing): disjoint ranges for all 12 points
up to the full block, refusal of oversized runs, the pair-counting helper, and
the Selene global-shot seeding test above.

**Rerun of B d=5 dt=1 only, with the fixed allocation**
(`results/stage6b/20261009T160001Z/`; Selene seed 110261009, injection seeds
[110261010, 110295160), now disjoint from every other point): **pass**, with
p_L 0.0279 [0.0257, 0.0303] against twin 0.0293, z = -1.41; detection max
abs z 1.66 against 3.36. The other 5 B points still share injection seeds with
each other in run `20261009T154022Z`. If strictly independent aggregates are
wanted, rerun them with
`--target-failures 1000 --only B_d3_hook_safe_dt1 B_d3_hook_safe_dt4 B_d3_hook_safe_dt6 B_d5_hook_safe_dt4 B_d5_hook_safe_dt6`
(about 1-3 min of Selene time, allowing for the instability). That rerun is
`20261009T161404Z` (section 8).

**The earlier 100-failure run is not independent of the 1000-failure run.** Both used
the same per-point seeds (base 20261009, stride 1000). Selene's output depends
only on the seeds and the global shot index (the test above), so each point's
100-failure shots are exactly the first shots of its 1000-failure run. That
follows from the construction; the raw records were not kept to compare
directly. `20261009T150628Z` is therefore a subset of `20261009T154022Z`, not a
second sample, and it is excluded from the final aggregates.

## 8. Independent rerun of the B points, and the pre-registered repeat

`20261009T161404Z` reran the 5 remaining B points at 1000 failures, each in its
own 10^7 seed block (Selene seeds 80261009-130261009). **4 pass; B d=5 dt=6
fails at z = -2.95** (Selene 0.2163 [0.2003, 0.2331], 916/4234, against twin
0.2362). The same point had passed at z = +0.06 in `20261009T154022Z`. Its
detection rates agree (max abs z 0.75 against 2.94), and its Selene/twin
duplicate counts match (136 / 139).

**Decision rule, fixed before the repeat:** rerun that point once at 3000
target failures in a fresh seed block. If |z| < 3, record the failure as a
fluctuation and close 6b. If |z| >= 3, diagnose the dephasing path at large dt
(f_quad dt = 0.6) before any 6c work.

**Repeat** (`20261009T161923Z`; `--seed 1000000000`, seed block [1.11e9, 1.12e9),
disjoint from every earlier seed, all of which lie below 1.41e8): Selene
0.2328 [0.2232, 0.2426] (2959/12712) against twin 0.2360 [0.2331, 0.2390],
**z = -0.82**; detection max abs z 2.01 against 2.94; duplicates 855 / 788, with
same/cross-process pairs 467 / 1433 (1 : 3 expected). **|z| < 3, so z = -2.95 is
recorded as a fluctuation.**

**Look-elsewhere effect.** A single point reaches |z| >= 2.95 with probability
0.0032. Across the 12 independently seeded points that existed when it was
observed, at least one such point appears with probability 0.037, and at least
one 99% failure with probability 0.114. One failure at this level, which did
not reproduce, is unremarkable.

## 9. Final aggregates, independently seeded measurements only

Included: the 6 A points of `20261009T154022Z` (their injection draws never
fire, and their error-model seeds are distinct per point), plus every B
measurement made in a disjoint seed block. That is 13 measurements of 12
points. Both B d=5 dt=6 measurements are kept: dropping the failing one would
bias the aggregate.

| run | d | variant | dt | Selene p_L [99% CI] (fail/shots) | twin p_L [99% CI] | z | det max abs z / crit | Selene seed |
|---|---|---|---|---|---|---|---|---|
| 154022Z | 3 | A | 1 | 0.0838 [0.0775, 0.0905] (1025/12231) | 0.0818 [0.0799, 0.0837] | +0.79 | 2.65 / 3.36 | 20261009 |
| 154022Z | 3 | A | 4 | 0.0738 [0.0681, 0.0798] (983/13324) | 0.0751 [0.0733, 0.0769] | -0.54 | 1.41 / 3.02 | 20262009 |
| 154022Z | 3 | A | 12 | 0.0909 [0.0840, 0.0982] (998/10980) | 0.0911 [0.0889, 0.0932] | -0.06 | 2.15 / 2.81 | 20263009 |
| 154022Z | 5 | A | 1 | 0.0526 [0.0485, 0.0569] (999/19010) | 0.0526 [0.0513, 0.0539] | -0.03 | 2.40 / 3.36 | 20264009 |
| 154022Z | 5 | A | 4 | 0.0589 [0.0543, 0.0638] (988/16768) | 0.0596 [0.0582, 0.0611] | -0.38 | 0.94 / 3.02 | 20265009 |
| 154022Z | 5 | A | 12 | 0.0890 [0.0822, 0.0962] (985/11068) | 0.0904 [0.0882, 0.0925] | -0.48 | 1.34 / 2.81 | 20266009 |
| 160001Z | 5 | B | 1 | 0.0279 [0.0257, 0.0303] (954/34150) | 0.0293 [0.0286, 0.0300] | -1.41 | 1.66 / 3.36 | 110261009 |
| 161404Z | 3 | B | 1 | 0.0510 [0.0470, 0.0552] (1004/19689) | 0.0508 [0.0495, 0.0521] | +0.12 | 1.65 / 3.36 | 80261009 |
| 161404Z | 3 | B | 4 | 0.1228 [0.1136, 0.1324] (1012/8241) | 0.1213 [0.1187, 0.1240] | +0.39 | 0.58 / 3.02 | 90261009 |
| 161404Z | 3 | B | 6 | 0.2040 [0.1892, 0.2194] (975/4780) | 0.2092 [0.2050, 0.2135] | -0.86 | 0.74 / 2.94 | 100261009 |
| 161404Z | 5 | B | 4 | 0.1115 [0.1030, 0.1204] (984/8825) | 0.1133 [0.1108, 0.1159] | -0.52 | 1.29 / 3.02 | 120261009 |
| 161404Z | 5 | B | 6 | 0.2163 [0.2003, 0.2331] (916/4234) | 0.2362 [0.2317, 0.2407] | -2.95 | 0.75 / 2.94 | 130261009 |
| 161923Z | 5 | B | 6 | 0.2328 [0.2232, 0.2426] (2959/12712) | 0.2360 [0.2331, 0.2390] | -0.82 | 2.01 / 2.94 | 1110000000 |

The pooled ratios below include the twin's own sampling error.

| aggregate | value |
| --- | --- |
| p_L, all 13: sum of z^2 | **13.81 on 13 dof (p = 0.39)**; mean z -0.52; Stouffer -1.87 |
| pooled Selene/twin failures, all 13 | **0.9854 +/- 0.0084** (14,782 vs 15,001; -1.7 sigma, p = 0.08) |
| variant A (6) | ratio 0.996 +/- 0.013 (p = 0.78); sum of z^2 = 1.28 |
| variant B (7) | ratio **0.978 +/- 0.011** (-2.0 sigma, p = 0.043); Stouffer -2.29 |
| variant B without the z = -2.95 measurement (6) | ratio 0.986 +/- 0.012 (p = 0.22); Stouffer -1.27 |
| detection rates, 81 rounds | sum of z^2 = 74.2 (p = 0.69); largest abs z 2.65 |

**Reading.**

- **Variant A** agrees with no visible bias, to +/-3.4% at 99%.
- **Variant B** leans slightly low: Selene's p_L is about 2% below the twin's. This is nominally p = 0.04, but it depends on the one measurement that did not reproduce. Without it, p = 0.22, and it is not significant once the look-elsewhere effect over the variant and aggregate cuts is counted.
- It is **carried forward as a bound, not a correction**: at 99% the twin over-predicts variant B's p_L by at most about 5%.
- There is one known variant-B-specific approximation in the twin: the preparation corrections are modelled by their per-qubit average, dropping an O(p_1q^2) correlation. Its expected size is far below 2%, so it is not offered as the explanation.
- For Stage 6c, a common 2-5% scale error in R would not move dt\*, which depends on the shape of R(dt). It would only matter if it varied with dt, so I tested that.
  - Per-dt Selene/twin ratios, A: 1.012 / 0.986 / 0.992 (+/-0.023) at dt = 1 / 4 / 12.
  - Per-dt ratios, B: 0.979 / 0.998 (+/-0.023) / 0.970 (+/-0.014) at dt = 1 / 4 / 6.
  - A weighted linear fit of B's ratio against dt gives a slope of -0.0025 +/- 0.0054 per unit dt (z = -0.47); without the z = -2.95 measurement, +0.0008 +/- 0.0056.
  - So no dt dependence is resolved. At 99%, a B trend of up to about 1.4% per unit dt over dt = 1-6 cannot be excluded. This limit is what to carry into 6c.

**Stage 6b is closed.** The twin, with `order="hook_safe"`, `native_charging=True` and
(for Selene comparisons) `mirror_guppy_program=True`, is the validated model
for Stage 6c.
