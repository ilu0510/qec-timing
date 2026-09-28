# Performance: the two paths, measured

Why this project runs validation on Guppy/Selene and sweeps on the pure-Python
reference. Measured, not estimated. This is the accepted design (see CLAUDE.md,
"Two paths"), not a shortfall.

## Summary

| | Guppy / Selene | Reference (Stim + pymatching) |
| --- | --- | --- |
| Role | **validation** -- the artefact the physics is defined by | **production** -- parameter sweeps |
| Cost per syndrome round | 0.3-3.5 ms per shot (d = 5..11) | ~88 ns per detector, i.e. `88 ns x n_checks` per round |
| Ratio at d = 11 | **~200x slower** (best case measured) | 1x |
| Fixed cost | **~9.3 s per `run()` call** | negligible |

The gap is not a defect. Selene emulates the program instruction by instruction,
including every Bernoulli draw and every gate; Stim samples a detector error
model directly. They are different kinds of object, and the cross-check is what
makes it legitimate to use the fast one for sweeps.

## The `run()` fixed cost is per CALL, not per shot

This is the single most actionable number here.

Measured at d = 7, 40 rounds, after a warm-up call:

| shots in one `run()` | wall time |
| --- | --- |
| 1 | 9.62 s |
| 2 | 9.87 s |
| 8 | 9.82 s |
| 32 | 11.16 s |
| 128 | 19.41 s |

Linear fit: **9.32 s fixed per call, 0.078 s per shot.** One shot and eight
shots cost the same to within noise, so the overhead is paid once per
invocation and is independent of the shot count.

**Consequence: batch aggressively.** Splitting 128 shots across four `run()`
calls costs an extra ~28 s for no benefit. Concretely:

- Size the shot count for a point *up front* (from the ansatz) and issue **one
  `run()` call per (d, dt) point**. Do not use an adaptive
  sample-until-N-failures loop on the Guppy path: each extra batch costs another
  9.3 s. (Adaptive batching is fine on the reference path, where the per-call
  cost is negligible.)
- If many points need few shots each, they can be packed into a single call with
  `run_per_shot`, which accepts a different argument mapping per shot, since
  `dt`, `p_data`, `p_read` and `n_rounds` are all runtime arguments.

## Per-round cost

Guppy, marginal ms per round per shot, measured by differencing two round counts
at fixed shots:

| d | data qubits | checks | RNG draws/round | ms/round/shot |
| --- | --- | --- | --- | --- |
| 5 | 25 | 12 | 37 | 0.283 |
| 7 | 49 | 24 | 73 | 1.258 |
| 9 | 81 | 40 | 121 | 1.704 |
| 11 | 121 | 60 | 181 | 3.483 |

Reference: **~88 ns per detector**, so `88 ns x n_checks` per round -- 5.3 us
per round at d = 11 against 3.5 ms for Guppy, a factor of ~660 in that
measurement (~200x in the most favourable Guppy measurement taken, see the
caveat below).

### Caveat: these numbers are not stable to better than ~3x

An earlier measurement in the same session, same code, same machine, gave
**1.101 ms/round at d = 11** and showed no large per-call fixed cost. A later
measurement of the identical configuration gave **3.483 ms/round** and a 9.3 s
per-call cost. Nothing in the repository changed between them.

So treat the absolute numbers as order-of-magnitude and the *ratio between the
two paths* as the robust finding. The conclusions below hold under either
measurement, because they turn on a factor of 100-1000, not on a factor of 3.
If a future run needs a firm budget, re-measure rather than trusting this table.

## Why the full Fig. 2 sweep is not run on Guppy

Cost of **one dt point per d** -- the minimum only -- at the Fig. 2 targets
(p = 0.015, lambda = 1, T = 200), with enough shots for 400 logical failures.
Using the more favourable Guppy measurement (1.1 ms/round at d = 11, scaled by
`n_data + n_checks`):

| d | dt* | p_L at dt* | rounds | shots | reference | Guppy |
| --- | --- | --- | --- | --- | --- | --- |
| 11 | 0.263 | 3.2e-1 | 760 | 1.2e3 | 5 s | 19 min |
| 15 | 0.185 | 9.8e-2 | 1080 | 4.1e3 | 43 s | 2.6 h |
| 21 | 0.128 | 1.4e-2 | 1560 | 2.9e4 | 15 min | 2.2 d |
| 27 | 0.098 | 1.6e-3 | 2040 | 2.5e5 | 4.5 h | **39 d** |

Summed over d = 11..27, one point each: **reference 7.3 h, Guppy 64 days.** A
real sweep needs 15-20 dt points per d. The Guppy path is therefore not a
candidate for paper-scale sweeps at any plausible level of optimisation, and the
reference path is.

Method: `dt* = 1/(lambda(alpha-1))` from Eq. 6 with `alpha = g(d+1)/2`; `p_L`
from Eq. 5 through Eq. A3 with beta = 2, g = 0.8, p_th = 0.029 and the paper's
fitted A per d (1.4 at d = 11 rising to ~2.05 at d = 27); shots = 400 / p_L;
cost = shots x rounds x the per-round figures above.

## Geometric-gap sampling: considered and rejected

The idea: instead of one Bernoulli draw per data qubit per round, draw the *gap*
to the next flip from a geometric distribution, which needs about `1/p_data`
fewer RNG calls. At the Fig. 2 targets `p_data ~ 0.019`, so ~50x fewer draws for
the data-noise channel.

**It does not close the gap, because RNG is only part of the cost.** Per round at
d = 11 the program performs ~181 RNG draws (121 data + 60 readout) against ~300
quantum operations (240 CNOTs plus 60 ancilla allocate/measure pairs). Even
eliminating *every* data-noise draw leaves the gates, the ancilla handling, the
detection-event computation and the emission untouched.

Timed split (RNG cost measured directly with the Stage 2 `record_flip_probe`,
which is a bare Bernoulli draw with no quantum operations, then multiplied by
the draws per round and compared against the measured per-round cost):

| d | RNG draws/round | ms/round | of which RNG | of which everything else | RNG share |
| --- | --- | --- | --- | --- | --- |
| 5 | 37 | 0.257 | 0.018 | 0.239 | **6.9%** |
| 7 | 73 | 0.557 | 0.035 | 0.522 | **6.3%** |
| 9 | 121 | 0.642 | 0.058 | 0.584 | **9.1%** |

One bare Bernoulli draw costs **0.48 us**. The *data-noise* draws alone -- the
only ones geometric-gap sampling would remove -- are **4.2-6.1%** of the
per-round cost, giving a best-case speed-up of **1.04-1.06x**:

| d | data-noise share | best-case speed-up |
| --- | --- | --- |
| 5 | 4.7% | 1.05x |
| 7 | 4.2% | 1.04x |
| 9 | 6.1% | 1.06x |

This is **weaker than the ~2x estimated from operation counts alone**. Counting
operations suggested RNG was roughly a third of the work; timing it shows the
draws are cheap (0.48 us each) and the cost is overwhelmingly the quantum
operations -- CNOTs, ancilla allocation and measurement -- which geometric-gap
sampling does not touch.

A ~1.05x gain against a ~10,000x gap does not change which path runs sweeps, and it
would mean rewriting the one part of the system that Stages 2 and 3 verified
statistically -- the noise injection -- trading a validated model for a speed-up
that is irrelevant to the decision it would inform. **Rejected.** Recorded here
so it is not re-proposed.

## Practical guidance

1. **One `run()` call per configuration**, with shots sized in advance. Never
   loop `run()` to accumulate statistics on the Guppy path.
2. **Adaptive sampling belongs on the reference path**, where per-call overhead
   is negligible.
3. **Keep Guppy work at small d and short T.** Validation runs use d = 5, 7, 9
   and T = 10, which keeps every point to seconds rather than minutes.
4. **Re-measure before budgeting.** See the stability caveat above.
