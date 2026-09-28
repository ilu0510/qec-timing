# Stage 4a report: the noise model is validated

Scope: a **validated noise model**, not a reproduction of Fig. 2. All three
success criteria are met.

Total compute: **22 min** for the main run, plus ~10 min for three follow-up
checks (seed sensitivity, RNG timing split, and the T = 30 test that identifies
why dt* sits above Eq. 6). Everything below is at p = 0.015,
lambda = 1, b_read = 1, **T = 10**, d = 5, 7, 9.

## Results against the criteria

| # | Criterion | Result |
| --- | --- | --- |
| 1 | Guppy and reference agree within error bars across the dt range | **18/18 points overlap** at 99.9%; offset consistent with zero across seeds |
| 2 | U-shaped curve per d with a minimum near Eq. 6's dt* | **Yes** at all three d; minima sit 17-31% above Eq. 6, traced to an effective g ~ 0.70 rather than 0.8 |
| 3 | 1/dt* roughly linear in d, slope with real error bars | **1/dt\* = (0.376 +/- 0.026) d - (0.821 +/- 0.189)**; paper's 0.402 is 1.0 sigma away |

## 1. Agreement (the main result)

Six dt points per d, spanning dt\*/5 to 4 dt\*, compared point by point against
the reference sampler. **Every one of the 18 comparisons overlaps** at 99.9%
confidence, with per-point Guppy statistics of 19-98 logical failures.

### The residual pattern, and why it is not a bias

The signs are not random: at d = 5 Guppy is high on all six points (+8 to +19%),
at d = 7 low on all six (-4 to -18%), at d = 9 mixed. Treating the points as
independent, that aggregates to about **+3.0 sigma at d = 5 and -3.3 sigma at
d = 7** -- which would look alarming if taken at face value.

It is an artefact of shared randomness, not a model error. All dt points at a
given d were run from one `base_seed`, so their noise realisations are
correlated and the per-d aggregate is not a valid independent test. Re-running
d = 5 under three different base seeds:

| base_seed | mean relative offset | combined z |
| --- | --- | --- |
| 20260927 | +13.5% | +2.98 |
| 20261927 | **-4.0%** | -0.99 |
| 20262927 | +3.6% | +0.85 |

The offset changes sign and magnitude with the seed. Averaged over seeds it is
**+4.4% with a spread of about +/-9%**, entirely consistent with zero given
~10% per-point error bars. A genuine defect in the noise model could not flip
sign between seeds, nor between d = 5 and d = 7.

**Lesson for future comparisons:** vary `base_seed` per dt point, or treat the
per-d aggregate as a single correlated measurement. Do not multiply per-point
sigmas together across a shared-seed grid.

## 2. U-shape and the location of the minimum

Every d shows a clear minimum. Located by a **parabola fit in log(dt)** on the
refined reference grid -- never by taking the lowest sampled point -- and fitted
to **log R** rather than log p_L, with Eq. A3 inverted first
(`R = -ln(1 - 2 p_L) / (2 T_actual)`), so the ceiling cannot distort curvature
and the fitted quantity is the one Eq. 6 actually minimises.

| d | fitted dt* | Eq. 6 dt* | deviation |
| --- | --- | --- | --- |
| 5 | 0.9343 +/- 0.0676 | 0.7143 | +30.8% |
| 7 | 0.5546 +/- 0.0114 | 0.4545 | +22.0% |
| 9 | 0.3888 +/- 0.0109 | 0.3333 | +16.7% |

The minima are systematically **high**, by an amount that **shrinks as d grows**.
That is not noise -- the d = 7 and d = 9 deviations are many sigma.

It has a clean interpretation. Eq. 6 is `dt* = 1/(lambda(alpha - 1))` with
`alpha = g(d+1)/2`, evaluated at the paper's fitted **g = 0.8**. Inverting each
measured dt* for the effective g:

| d | alpha_eff | implied g_eff |
| --- | --- | --- |
| 5 | 2.070 | 0.690 +/- 0.026 |
| 7 | 2.803 | 0.701 +/- 0.009 |
| 9 | 3.572 | 0.714 +/- 0.014 |

So at these small distances the effective **g is about 0.69-0.71 and rising
towards 0.8 with d**. The paper fits g = 0.8 against d = 11-27; we are well
below that range. The discrepancy is a small-d effect in the *ansatz parameter*,
not a disagreement about where the optimum is -- and the trend points the right
way.

### Which explanation: tested, and it is mostly (a)

Two candidates for the offset: **(a)** the effective g really is below 0.8 at
small d, or **(b)** a finite-rounds boundary effect -- at T = 10 the large-dt
points run very few rounds (d = 9 at dt = 1.33 uses only 8, fewer than d), so
the final-readout detector block is a large share of the history and the
paper's T/dt = O(d) condition is violated.

**The discriminator is T.** In the ansatz dt* does not depend on T at all: R(dt)
is T-independent and `p_L = 1/2(1 - exp(-2 R T))` is monotone in `R T`, so
`argmin_dt p_L = argmin_dt R` for any fixed T. Tripling T triples the round
count at every dt while leaving the predicted optimum untouched. If the minima
move toward Eq. 6 it was (b); if they stay put it is (a). Re-ran the scan at
**T = 30** on the reference path (3.2 min):

| d | dt* at T=10 | dt* at T=30 | Eq. 6 | deviation T=10 -> T=30 | shift |
| --- | --- | --- | --- | --- | --- |
| 7 | 0.5546 | 0.5565 +/- 0.0102 | 0.4545 | +22.0% -> +22.4% | **0.1 sigma** |
| 9 | 0.3888 | 0.3673 +/- 0.0079 | 0.3333 | +16.6% -> +10.2% | 1.6 sigma |

At T = 30 every sampled point has at least 22 rounds (>= 2.4d), so
round-starvation is eliminated. **d = 7 did not move at all** (0.1 sigma), and
d = 9 moved only 1.6 sigma -- not significant, though in the direction (b)
predicts, which is unsurprising since d = 9 was the case with rounds < d at
T = 10.

**Conclusion: the offset is predominantly (a), an effective g below 0.8.** A
boundary effect cannot account for it -- removing the boundary entirely leaves
both minima well above Eq. 6 (+22.4% and +10.2%). Implied g:

| d | g_eff at T=10 | g_eff at T=30 |
| --- | --- | --- |
| 7 | 0.701 | 0.699 |
| 9 | 0.714 | 0.745 |

Stable at d = 7 and still clearly below the paper's 0.8, with the same rise
towards 0.8 as d increases. A modest (b) contribution at d = 9 is plausible but
not established.

Practical consequence: **Eq. 6 with g = 0.8 under-predicts dt* by 10-30% at
d = 5-9.** Using the locally fitted g ~ 0.70 instead reproduces the measured
optima. Nothing here affects the validation result; it is a statement about the
ansatz parameter outside the range it was fitted in.

## 3. 1/dt* versus d

Weighted linear fit over the three distances:

> **1/dt\* = (0.376 +/- 0.026) d - (0.821 +/- 0.189)**

- The paper's slope is **0.402** (Fig. 9). Ours is **1.0 sigma** away. Agreement.
- The theory fixes the slope to `g/2`, so our slope implies **g = 0.752 +/-
  0.053**, consistent with the paper's 0.8 within 0.9 sigma.
- The theory also fixes the intercept to `z = eta - 1 = -0.624`. Measured
  -0.821 +/- 0.189, **1.04 sigma** away. The fit is internally consistent with
  the ansatz's structure, not just its slope.

Note the tension worth recording: the slope route gives g = 0.752 +/- 0.053
while the per-d inversions give 0.69-0.71. Both are below 0.8 and both trend
towards it; with only three closely spaced distances the slope has limited
leverage. Nothing here needs resolving for a validated noise model, but it would
be the first thing to look at if the model were ever pushed to larger d.

## Choices and why

**T = 10 for every d.** Chosen so p_L stays unsaturated -- measured 0.048 to
0.27 across every sampled point, against Eq. A3's ceiling of 0.5. Saturation
would compress the U and bias the vertex fit. The same value also gives ~3d
rounds at the optimum (14, 22, 30 for d = 5, 7, 9), matching Appendix C's
convention and the paper's requirement that T/dt = O(d). One flat T satisfying
both across d = 5-9 is a convenience, not a coincidence to rely on at larger d.

**d = 5, 7, 9; not 11.** d = 11 was costed at ~74 min of Guppy time on its own,
against ~18 min for d = 5-9 together. Dropped as agreed.

**~100 failures per point, capped by time.** Guppy per-round cost is unstable to
~3x between runs (docs/performance.md), so shots were capped by a **per-point
time budget** of 150 s rather than a shot count -- bounding total runtime by
construction. The cap bound on 5 of 18 points, all at the expensive small-dt end
of d = 7 and d = 9, where Guppy statistics fell to 19-36 failures (~17-23% error
bars). Those points are the weakest part of the agreement test and are marked in
`stage4a_crosscheck.csv` via `shots_time_capped`.

**Guppy ran the coarse grid; the reference ran the refinement.** The two-path
division of labour doing real work: the 9-point refined grid per d costs seconds
on the reference and would have cost tens of minutes on Guppy, and criterion 1
has already established the two agree.

## Follow-ups recorded elsewhere

- **`run()` overhead is per call, not per shot** -- 9.3 s fixed, independent of
  shot count. Shots are sized up front and each point is a single call. Full
  measurement in `docs/performance.md`.
- **Geometric-gap sampling rejected, with the timing split now measured.** RNG
  is only **6-9%** of the per-round cost (a bare Bernoulli draw costs 0.48 us);
  the data-noise draws alone are **4-6%**, so the best case is a **1.04-1.06x**
  speed-up. This is *weaker* than the ~2x I first estimated from operation
  counts -- timing shows the quantum operations dominate, not the draws.
  Against a ~10,000x gap it changes nothing.

## Files

```
results/stage4a/stage4a.json              full results
results/stage4a/stage4a_crosscheck.csv    the 18 comparisons
results/stage4a/stage4a_points.csv        every sampled point, both paths
results/stage4a/stage4a.png               p_L vs dt, and 1/dt* vs d
results/stage4a/seed_sensitivity.json     the three-seed check
results/stage4a/rng_split.json            measured RNG vs gates
docs/performance.md                       cost model and the two-path rationale
```

## Status

The noise model is validated: the Guppy implementation and the independent
reference agree everywhere tested, the physics shows the expected U-shape, and
the optimal-interval scaling reproduces the paper's slope within 1 sigma. The
remaining deviations are understood (small-d g, shared-seed correlation) rather
than open.
