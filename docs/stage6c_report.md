# Stage 6c report: how dt\* scales with d on the circuit-level twin

Run: `results/stage6c/20261010T010451Z/` (`points.csv`, `fits.json`,
`R_vs_dt.png`, `dtstar_vs_d.png`, `run.log`, `summary.json`), script
`scripts/stage6c_sweep.py` at commit `bcd164f`, in the twin configuration that
6b validated against Selene (`order="hook_safe"`, `native_charging=True`).
Twin only: no Selene. Versions: stim 1.16.0, PyMatching 2.4.0, numpy 2.5.3.

**Answer.**

- **Variant A (linear idle):** dt\* falls as 1/(λ(α-1)). 1/dt\* is linear in d (decisively preferred over 1/dt\*²). The log-log slope over d = 3-11 is -1.29 ± 0.02: steeper than the asymptotic -1, as the α-1 offset predicts at small d.
- **Variant B (quadratic idle):** dt\* falls with slope -0.45 ± 0.015 over d = 3-9, or -0.49 ± 0.03 over d = 5-9, consistent with -1/2. Over this d range, though, the two linearised forms fit B equally well, so B's 1/√d law rests on the slope and the per-d shape fits, not on the linear-form comparison.

## 1. Budget overrun: what was skipped

The sweep had a 45 min (2700 s) cap. It ran 2879 s in total (2853 s of work after setup). **123 of 192 planned points ran; 69 were skipped by the budget:**

| skipped | count | consequence |
| --- | --- | --- |
| all 6-point refinements around each minimum (A and B, d = 3-11) | 60 | dt\* comes from the 12-point grid alone (spacing ×1.52). Since the shape fit uses every point, this costs precision only. |
| both 2T spot checks (A and B at d = 5) | 2 | **T-independence of R is not verified in this run.** |
| variant B, d = 11: grid points 6-11 (dt ≥ 0.678) and the baseline dt | 7 | **No points to the right of the minimum, so B's d = 11 dt\* is an extrapolation and is excluded below.** |

Why:

- **Slow run:** the run's own plan predicted 3831 s against the 2700 s budget, while my `--plan-only` before handing over had predicted 1806 s. The run's pilots took 329 s against 141 s in that plan, so everything ran about 2× slower than planned.
- **Skip logic:** the script skips a point whose predicted cost exceeds what is left (as designed), but a point already running cannot be stopped. The last stage-1 point overran the cap by about 150 s.
- **Ordering flaw (mine):** refinements and spot checks ran after every stage-1 grid point, so they were the first to be cut. A rerun should order them before the costliest stage-1 points, or give them a reserved share of the budget.

Every point that did run has at least 1027 failures (the target was ≥ 1000).

## 2. Configuration (fixed before the sweep)

- **Gate noise:** p_2q = 4e-3, p_1q = 8e-5 (ratio 50, as in 6b), p_meas = 1.6e-3, p_init = 0. With idle off this gives a suppression factor Λ = 2.6-3.2 per d+2 in both memories. At p_2q = 5e-3, Λ falls to 1.9 at d = 3->5.
- **Idle:**
  - A: k = 100 charged DEPOLARIZE1(p_1q) ops per unit time, in a Z memory.
  - B: f_quad = 0.1, f_lin = 0 (p_z = sin²(0.05 dt)), in an X memory. At the minima f_quad·dt ≤ 0.11, so this is the quadratic regime.
- **Grid:** 12 log-spaced points over [0.1, 10] × the coarse dt\* (from a pre-sweep scan) per (variant, d), plus a d-independent baseline dt (2.0 for A, 1.1 for B).
- **Per point:**
  - rounds chosen so p_L ≈ 0.2 (at least 2d, at most 60);
  - sampled until ≥ 1000 failures;
  - points with p_L > 0.4 flagged and excluded;
  - R from Eq. A3 inverted, as in Stage 4a;
  - each point in its own disjoint seed block (recorded in `points.csv`; a run that would overflow its block is refused).
- **dt\* estimator:** a per-d fit of the ansatz shape, log R = a - log dt + α_d log(1 + c_d u), with u = dt (A) or dt² (B), and a 400-sample binomial bootstrap for 95% intervals. A symmetric parabola was rejected as the main estimator: on synthetic ansatz data it is biased 5-13% low. It is kept as a cross-check.

## 3. dt\* per d

| variant | d | dt\* [95%] | α_d | idle coeff. at this d | points / right of minimum | parabola cross-check |
| --- | --- | --- | --- | --- | --- | --- |
| A | 3 | 2.168 [2.069, 2.283] | 1.47 | λ_d = 0.98 | 12 / 5 | 2.025 |
| A | 5 | 1.172 [1.136, 1.212] | 1.82 | 1.03 | 13 / 6 | 1.017 |
| A | 7 | 0.753 [0.728, 0.781] | 2.21 | 1.10 | 13 / 7 | n/a (window) |
| A | 9 | 0.544 [0.529, 0.560] | 2.81 | 1.02 | 13 / 7 | n/a |
| A | 11 | 0.411 [0.400, 0.422] | 3.32 | 1.05 | 13 / 7 | n/a |
| B | 3 | 1.089 [1.060, 1.120] | 1.57 | κ_d = 0.40 | 11 / 5 | 1.036 |
| B | 5 | 0.885 [0.859, 0.907] | 1.84 | 0.48 | 11 / 5 | 0.977 |
| B | 7 | 0.761 [0.744, 0.777] | 2.41 | 0.45 | 11 / 5 | 0.732 |
| B | 9 | 0.664 [0.651, 0.677] | 2.88 | 0.48 | 11 / 5 | 0.694 |
| B | 11 | *0.598 [0.545, 1.333]: not bracketed, excluded* | (fit at bound) | | 6 / 0 | |

"n/a" means fewer than four points fell inside the parabola's factor-2 window, so the parabola cross-check was not computed for that d.

Two consistency checks the ansatz implies, done per d without assuming the ansatz across d:

- **The idle coefficient is d-independent:** λ_d = 0.98-1.10 for A, κ_d = 0.40-0.48 for B (d = 3-9).
- **α_d grows with d, but not as g(d+1)/2 with one g:** the implied g_d = 2α_d/(d+1) is 0.74 at d = 3 and 0.55-0.61 for d ≥ 5 (A), and 0.78 then 0.58-0.61 (B). The effective g is about 0.55-0.6 at circuit level, against the paper's phenomenological 0.8, and d = 3 is an outlier.

## 4. Scaling of dt\* with d

Fit uncertainties come from the bootstrap spread of log dt\* at each d. Weighted fits; χ² is quoted with its degrees of freedom.

| variant | d range | log-log slope | 1/dt\* linear in d: χ² (p) | 1/dt\*² linear in d: χ² (p) |
| --- | --- | --- | --- | --- |
| A | 3-11 | **-1.292 ± 0.020** (χ² 3.3/3) | 31.5/3 (p < 0.001) | 372.9/3 (p < 0.001) |
| A | 5-11 | -1.324 ± 0.027 | **4.9/2 (p = 0.085)** | 62.8/2 (p < 0.001) |
| B | 3-9 | **-0.451 ± 0.015** (χ² 3.7/2) | 0.73/2 (p = 0.69) | 2.19/2 (p = 0.34) |
| B | 5-9 | -0.492 ± 0.028 | 0.04/1 (p = 0.85) | 0.85/1 (p = 0.36) |

**Which form fits better, in numbers:**

- **A:** 1/dt\* linear in d fits far better than 1/dt\*² (χ² 31.5 vs 373 on 3 degrees of freedom; 4.9 vs 62.8 without d = 3), which is the expected form. The remaining misfit over d = 3-11 comes from d = 3, the same outlier as in g_d.
- **B:** both forms fit (p = 0.69 and 0.34), so the data do not separate them. dt\* only changes by 1.6× over d = 3-9, so a straight line in d fits both 1/dt\* and 1/dt\*². The quadratic-law evidence for B is the slope (-0.45 ± 0.015, against A's -1.29 in the same pipeline) and the per-d fits (quadratic idle term with a d-independent κ).

**About the slopes:** neither -1 nor -1/2 is the right yardstick at d = 3-11. With the finite-d offsets (α-1 and 2α-1), even data generated exactly from the ansatz give steeper slopes; synthetic ansatz data gave -1.32 and -0.53 at d = 3-9. The measured -1.29 and -0.45 are in line with that. The linearised forms test the laws themselves.

## 5. Global ansatz fit (Eq. 5, all d at once)

| | A (linear idle) | B (quadratic idle) |
| --- | --- | --- |
| A (prefactor) | 0.0238 [95%: e^-3.80, e^-3.69] | 0.0151 |
| β | 0.12 [0.02, 0.20] | -0.31 [-0.41, -0.20] |
| p/p_th (fitted, single number) | 0.342 | 0.310 |
| g | **0.643** [0.629, 0.655] | **0.686** [0.662, 0.711] |
| idle coefficient | λ = 0.826 | κ = 0.374 |
| χ²/dof | 660.6/59 = 11.2 | 363.5/45 = 8.1 |

The global fits are poor (χ²/dof ≈ 8-11). A single g cannot capture the d = 3 outlier from section 3, and B's β comes out unphysical (negative). Measured dt\* against the global-fit prediction:

| d | A measured | A predicted | B measured | B predicted |
| --- | --- | --- | --- | --- |
| 3 | 2.17 | 4.25 | 1.09 | 1.24 |
| 5 | 1.17 | 1.31 | 0.885 | 0.927 |
| 7 | 0.753 | 0.771 | 0.761 | 0.772 |
| 9 | 0.544 | 0.547 | 0.664 | 0.676 |
| 11 | 0.411 | 0.424 | (excluded) | 0.608 |

For d ≥ 7 the global ansatz predicts dt\* within 2-3%. At d = 3 it fails for A by 2× (α(3) = 1.29 sits close to 1, where dt\* = 1/(λ(α-1)) blows up). The per-d fits, not the global ansatz, are the measurement.

## 6. Advantage of the optimal interval (Eq. D9 analogue)

Γ = R(dt_base) / R(dt\*), against a fixed, d-independent baseline: the coarse dt\* at d = 3 (A: dt_base = 2.0, B: 1.1). Bootstrap 95% intervals; the prediction uses the global ansatz.

| d | A Γ measured | A predicted | B Γ measured | B predicted |
| --- | --- | --- | --- | --- |
| 3 | 1.02 [0.96, 1.07] | 1.07 | 1.00 [0.95, 1.06] | 1.01 |
| 5 | 1.11 [1.04, 1.17] | 1.04 | 1.01 [0.97, 1.06] | 1.02 |
| 7 | 1.25 [1.17, 1.32] | 1.33 | 1.12 [1.06, 1.18] | 1.13 |
| 9 | 1.86 [1.75, 1.98] | 1.90 | 1.28 [1.21, 1.36] | 1.29 |
| 11 | **2.68 [2.53, 2.86]** | 2.86 | (no baseline point) | 1.51 |

The gain from tuning dt grows with d, and faster for linear idle than for quadratic idle (2.7× vs 1.3× by d = 9-11). That matches dt\* falling faster for A.

## 7. Systematics and caveats

- **Twin vs Selene (from 6b):** for variant B the twin may over-predict p_L by up to about 5% (99%), with any dt-trend bounded by about 1.4% per unit dt. A common scale error does not move dt\*. A trend of that size would shift ln dt\* by -ε·dt\*/curvature, using the curvature of ln R at the minimum from the per-d fits: about 1.1% at d = 3 and 0.6% at d = 9. That is comparable to the bootstrap intervals, so it does not change the scaling conclusions.
- **T-independence of R is unverified here** (both spot checks were skipped). Rounds were chosen so p_L ≈ 0.2 with at least 2d rounds, which keeps the boundary contribution small, but it is not measured in this run.
- **B, d = 11, is excluded** (section 1); B's scaling uses d = 3-9.
- **Coarse grid only** (spacing ×1.52, no refinement). The bootstrap intervals include this.
- **Effective g ≈ 0.55-0.6** (d ≥ 5) at circuit level, not the paper's 0.8, and d = 3 does not follow the large-d form.

## 8. To complete the run (not done)

1. Rerun the 7 skipped B d = 11 grid points, which brackets B's d = 11 minimum.
2. Run the two 2T spot checks (cheap; about 15 s on the planning run).
3. Optionally, the refinements. The script needs a `--only`-style option or a resume from `points.csv` to do these without repeating stage 1, and its plan should reserve budget for refinements and spot checks ahead of the costliest grid points.
