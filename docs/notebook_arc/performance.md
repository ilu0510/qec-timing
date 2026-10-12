# Performance: what the notebook arc costs to run

The counterpart of [`../stub_arc/performance.md`](../stub_arc/performance.md). Less was measured
here: the notebooks record results, not timings. What follows comes from two places, and each
number says which:

- **Selene logs**: `results/circuit_level/logs/d{5,7}_20261001_*/worker*.log`, which time every
  block of the 06b/06c Selene sweeps.
- **Cache timestamps**: the `saved` field of every `results/circuit_level/cache/*/*/meta.json`.
  The span from first to last save of an experiment is an **upper bound** on its run time, and is
  only meaningful when the run was continuous.

Machine: one workstation, 16 worker processes (`N_WORKERS = min(16, cpu_count)`).

## Summary

| Path | Used in | Cost |
| --- | --- | --- |
| **Selene** (`surface_memory.Geometry`) | 03–06c | 0.50–0.65 ms per round per shot at d = 5, 0.98–1.19 at d = 7 |
| **Stim** (all `*Memory` classes) | 07–09 | a whole notebook's sweep (~1,000–1,600 cached jobs) in ~15–20 min |
| Native gate extraction (`native_memory_ops`) | 07b, 09 | one Selene run per `(d, rounds)`, cached; ~70 min for all 58 the first time |

As in the stub arc, Selene is the expensive path and Stim the cheap one. The difference is what
each is for: here Selene is not the validated artefact, only the source of native gate sequences
and of the 03–06c results. From 07 on, everything is sampled in Stim.

## Selene: per-round cost (logs)

06b/06c sweep: rotated code, `T = 96` idle units, 10,000 shots per block, 46 blocks per `τ`, both
decoders run on each block. Per-block wall time over the round count:

| τ | R | d = 5: s/block | ms/round/shot | d = 7: s/block | ms/round/shot |
| --- | --- | --- | --- | --- | --- |
| 1 | 96 | 480 | 0.500 | 936 | 0.975 |
| 2 | 48 | 252 | 0.524 | 492 | 1.025 |
| 3 | 32 | 178 | 0.557 | 338 | 1.058 |
| 4 | 24 | 146 | 0.608 | 265 | 1.105 |
| 6 | 16 | 103 | 0.646 | 191 | 1.193 |
| 8 | 12 | 76 | 0.629 | 137 | 1.145 |
| 12 | 8 | 51 | 0.640 | 90 | 1.128 |
| 16 | 6 | 38 | 0.633 | — | — |

- The per-round cost rises ~25% from τ = 1 to τ ≥ 6: each Δt unit adds one `phased_x` per data
  qubit to the round.
- Whole runs: d = 5 took 68 min, d = 7 took 127 min (`run.info`, 46 blocks × 16 workers each).
- This program measures **both** check sectors with a fresh ancilla array per round, so it does
  roughly twice the stub's work per round (stub: Z-checks only, one reused ancilla; 0.28 ms/round
  at d = 5 in the stub arc's measurement). The two were measured on different runs and code; do
  not read the ratio as more than order of magnitude.
- Each block is one Selene `run_shots` call, so the stub arc's ~9.3 s fixed cost per call is
  under 2% of a block here.

## Stim: per-notebook cost (cache timestamps)

| Experiment | Cached entries | First → last save | Continuous? |
| --- | --- | --- | --- |
| `07b_native_dephasing` | 1,116 | 13 min | yes |
| `08_uniform_noise_premise` | 1,605 | 20 min | yes |
| `09_native_uniform_noise` | 1,525 | 18 min | yes |
| `07_dephasing_interval` | 1,036 | 6 days | **no**: built up over several sessions |
| `stim_memory` | 183 | 2.9 days | no |

Each entry is one job: a `(d, idle setting, Δt, seed)` point at 10⁵ shots, an edge measurement at
4×10⁵ shots, or a phenomenological point. 07b and 09 also need the native gate sequences, cached
separately under `native_memory_ops` (58 entries, extracted from Selene over 70 min on first use).

**Cost model used to schedule jobs** (`run_all` in 08/09): `d² × R × shots`, largest first. It
orders work; it is not a calibrated time estimate.

## Re-running from scratch

`results/circuit_level/cache/` is gitignored (204 MB here). On a fresh clone every notebook
recomputes. From the numbers above, on a 16-worker machine:

- **07b, 08, 09:** ~15–20 min each, plus ~70 min once for native gate extraction.
- **06b/06c Selene sweeps:** ~1 h at d = 5, ~2 h at d = 7.
- **07:** not known; its cache was built over several sessions.
- **05:** four parameter sets at 3 × 10⁵ shots per point, 850–1,060 s each (the notebook's own
  printed timings).

Later notebooks read earlier ones' cache entries (09 §5 reads 07, 07b and 08), so they must be
run in order.

## Practical guidance

1. **Keep the cache.** It is the only record of a run's environment (git commit, package
   versions), and re-running costs hours.
2. **Extract native gate sequences before forking.** Selene is touched only in the parent
   process; workers read the cache (09 `run_all`).
3. **Re-measure before budgeting.** None of these numbers was taken as a benchmark; the stub arc
   found Guppy timings unstable to ~3× between runs on the same machine.
