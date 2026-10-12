# Notebook arc: summary for the team

Scope: `notebooks/03`–`09` and `src/qec_timing/circuit_level/`. Status: **premise holds in
simulation; no device parameter has been measured.** Details in
[`noise_spec.md`](noise_spec.md) (the noise models), [`interfaces.md`](interfaces.md) (the
code), [`performance.md`](performance.md) (what a run costs), [`notation.md`](notation.md)
(symbols and time units) and [`hardware_tests.md`](hardware_tests.md) (what a device can test).
The notebooks themselves are the per-step reports.

## The question

The stub arc (`docs/stub_arc/`) imposes the paper's phenomenological model and checks that the
implementation reproduces the paper. This arc asks the opposite: **start from a compiled
surface-code circuit with realistic noise, and see whether the paper's model comes out of it.**
The paper's model has two knobs per round:

```
p_data = p_stab ⊕ p_idle(Δt)        a flip on each data qubit, once per round
p_read = b_read · p                 a flip of each recorded syndrome bit
```

with `λ` (in `p_idle`) and `b_read` assumed. Here neither is assumed; both are measured from the
circuit's detection events.

## The ladder

Each step adds one ingredient and keeps everything else fixed.

| Step | Circuit | Idle noise | Memory | Finding |
| --- | --- | --- | --- | --- |
| **03** | single qubits in Selene | — | — | Selene's `p_meas` is the paper's `p_read`. Selene has no pure bit-flip and no duration term |
| **04, 04c** | Selene, compiled rotated/unrotated code | explicit `phased_x(0,0)` ops, charged `p_1q` | Z | effective `λ` ≈ 0.10–0.15, `b_read` ≈ 1.35–1.50: set by circuit structure, not hardware |
| **04d** | Selene | as 04 | Z | readout error composes as two independent flips `a ⊕ b` |
| **05** | paper's model | effective vs paper parameters | Z | effective parameters lengthen Δt\* 2.6–3.0×, almost all through `λ` |
| **06a** | Stim copy of Selene's circuit | as 04 | Z | gate-for-gate copy, statistically indistinguishable from Selene |
| **06b** | Selene shots, two decoders | as 04 | Z | decoder changes `p_L` up to 7–8× and decides whether an optimum appears |
| **06c** | Selene, full-DEM decoding | as 04 | Z | interior optimum: Δt\* = 12.0, 5.1, 3.4 at d = 3, 5, 7 |
| **07** | Stim textbook circuit | **Z only, linear `γ` + quadratic `c`, on data once per interval** | **X** | interior optimum on all 21 curves; at H2's estimated ratio, back to back is best except a small gain at d = 3 |
| **07b** | Selene's native gates, in Stim | as 07 | X | Δt\* 0–7% longer than textbook; 07's conclusions stand |
| **08** | Stim textbook circuit | **every qubit, through the round (7 layers) and the wait** | X | premise holds (P1–P3, below) |
| **09** | native gates | every qubit, through the round (18 layers) and the wait | X | premise holds; numbers barely move from 08 |

Why the basis changes at 07: the idle noise trapped ions have is Z (dephasing), and a Z-basis
memory decoded on Z-checks cannot see Z errors. 07 §9 finding 1: turning dephasing on moves the
Z-basis `p_L` by at most 2.1σ, the X-basis `p_L` by up to 304σ. The Z-basis results of 04–06c
are therefore about Selene's depolarizing idle ops, not about dephasing.

## The premise (08, 09)

| | Claim | 08 (textbook gates) | 09 (native gates) |
| --- | --- | --- | --- |
| **P1** | an interior optimum Δt\* exists | all 9 curves | all 9 curves |
| **P2** | `p_read` does not depend on Δt | 0.996–1.010 × its Δt = 1 value | 0.995–1.003 × |
| **P2** | `p_data(Δt) = p_data(1) ⊕ q_wait(Δt − 1)` | 0.996–1.004 × | 0.996–1.002 × |
| — | measured `b_read`, gates only | 0.98–1.02 | 0.95–0.99 |
| — | measured `b_read`, with linear dephasing | 1.08–1.15 | 1.05–1.09 |
| **P3** | Δt\* predicted from one calibration at Δt = 1 | +16–18% (d = 3), +3–6% (d = 5), −5…−2% (d = 7) | +17–21%, +4–9%, −2…+1% |
| — | `p_L` at the optimum, prediction / circuit | 0.95–1.11, 0.74, 0.51 | 0.94–1.09, 0.73–0.76, 0.51 |

Settings: H2-1E gate ratios scaled ×4.76 (`p_2q = 5×10⁻³`), `T = 192`, d = 3, 5, 7, idle noise
linear (`γ = 0.3·p_2q`), quadratic (`c = 0.01·p_2q`) or both.

**For the pitch: predict Δt\*, not `p_L`.** The phenomenological model's Δt\* error shrinks with
distance while its `p_L` error grows (08 finding 5).

## What it means for H2 (07 §7b, 09 §6)

At H2's **estimated** idle-to-gate ratio (idle error over one round `q(1) ≈ 1.1·p_2q`, from an
assumed 6 layers per round):

- **Linear idle noise:** a short wait helps only at d = 3 (Δt ≈ 2.3, +11–16%). Nothing at d ≥ 5
  with H2's real gate rates.
- **Quadratic idle noise:** 07 says back to back is best; 09 says a short wait (Δt ≈ 1.3) helps
  13–18%. The difference is whether the quadratic term restarts at every layer, which is
  unresolved. **Do not quote the quadratic H2 rows** (09 finding 4).

## Two things that are easy to get wrong

- **The decoder.** A hand-built phenomenological matching graph misses mid-round and propagated
  faults, and at d = 3 it erased the optimum entirely (06b). Decode with the circuit's full
  detector error model.
- **Results that hold by construction.** In these Stim models the ancilla is reset every round
  and the wait is a separate memoryless channel, so P2 is largely structural: the 08/09 checks
  show the parameters are *measurable* from detection events, not that a device keeps them
  *stable*. That is what [`hardware_tests.md`](hardware_tests.md) T9–T10 would test.

## Not covered

- **Every idle parameter is assumed.** `γ`, `c`, how `c` splits over time, `t_R` and the layers
  per round have never been measured on a device (`hardware_tests.md` §2).
- **Equal layer durations.** Measurement, transport and gates each take `t_R / L`. On ions
  measurement and transport are slower.
- **Not in any model:** leakage, crosstalk, asymmetric readout (the average `p_meas` is used),
  coherent errors (the quadratic term is Pauli-twirled), `Rz` noise (off in 07b/09).
- **Placement of `p_2q`** in the Selene mirror is by analogy, not read off Selene's log
  (`selene_mirror`, `allow_unverified_p2q`).
- **Distances 3–7 only**, rotated code from 07 on.
- **No user-controlled delay exists** in Guppy, tket or pytket (`hardware_tests.md` §3), so
  none of the waiting can yet be written for a device.
