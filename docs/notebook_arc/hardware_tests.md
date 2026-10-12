# Notebook arc: what can be tested on hardware

Status: **plan, nothing run on hardware yet.** Scope is the circuit-level arc,
`notebooks/03`–`09` (with `noise_analysis`, `qec_timestep` and
`capability_audit` as groundwork). That arc is simulated in Selene and Stim.
This document lists which of its claims a device can test, how, and in what
order. The goal is to calibrate the Stim model on the device, then check that
the model's notebook-by-notebook evolution conforms to it.

Labels: **notebook** = a finding in that notebook; **checked** = verified on
2026-10-12 against installed or latest packages and Quantinuum's docs;
**interpretation** = my reading, not tested.

## 1. What the arc predicts

- **Premise holds as an effective model** (08, 09). A circuit with:
  - trapped-ion idle noise (Z only, linear γ and time² c);
  - native gates;
  - the same noise on every qubit;
  - spreading two-qubit errors

  still looks like the paper's model through its detection events. p_read is
  flat in Δt (±0.5–1%). p_data(Δt) = p_data(1) ⊕ p_idle(wait) (±0.4%). b_read
  ≈ 1 at H2's gate ratios.
- **The optimum can be predicted from one calibration** (09). Within +17–21%
  (d=3), +4–9% (d=5) and ±2% (d=7). The model's p_L, however, is
  underestimated more as d grows (predicted / simulated p_L = 0.51 at d=7).
- **At H2's estimated idle-to-gate ratio, back to back is best** (07, 07b, 09).
  - Under linear idle noise, a wait helps only at d=3 (Δt ≈ 2.3, 11–16%).
  - Under time² idle noise the answer depends on an unresolved assumption, and
    09 says not to quote it.
- **These rest on parameters never measured on a device:**
  - γ and c: estimated from q(1) ≈ 1.1·p2, assuming 6 layers per round;
  - how the time² term splits over time;
  - the round duration t_R.

## 2. Parameters to calibrate

| Parameter | Used in | Current source |
|---|---|---|
| p1, p2 | 07–09 gate noise | H2-1E emulator sheet (2025-07-16) |
| p_meas (0→1, 1→0), p_init | 07–09 | H2-1E sheet; readout averaged |
| γ (`dephase_linear`) | 07–09 idle | **assumed** |
| c (`dephase_quad`) and its time split | 07–09 idle | **assumed / unresolved** |
| idle X/Y part (`idle_depol`) | 07 assumes 0 | **assumed** (from the emulator docs' Z-only form) |
| t_R, layers per round | 08, 09 time grid | **unknown** |

## 3. Blocker: there is no user-controlled delay (checked)

| Where | Finding |
|---|---|
| guppylang 1.0.3 (installed) and **1.1.3** (latest) std library | No sleep, delay or wait; only `barrier` |
| tket-exts 0.14.3 Helios ops | `PhasedX, Rz, ZZPhase, Reset, LazyMeasure*, QFree, TryQAlloc, RuntimeBarrier, FutureToMeasurement`. `RuntimeBarrier` takes no duration |
| Selene, by experiment | `h; barrier; barrier; h` executes with the barriers removed |
| Selene runtime header | Has `sleep_time` on barriers and `selene_simulate_delay` "for time-based noise modelling"; Guppy cannot reach either |
| pytket-quantinuum 0.59.3 | No sleep or delay |
| [Helios performance docs](https://docs.quantinuum.com/systems/trainings/helios/getting_started/helios_performance_maximization.html) | "Runtime behavior is not specifiable by the user." The advice there is to *reduce* idle time |
| H2 (QASM / pytket) | No documented `sleep` found; not confirmed either way |

**Resolve by:**
1. **Ask Quantinuum.** Is there a delay on Helios or H2? Is there an approved
   way to idle that the compiler keeps? Does Nexus report executed duration?
2. **Fallback: spectator idling.** The target qubit idles while k gate layers
   run on *other* qubits, ordered with barriers. The emulator docs say
   dephasing accrues "while idling while other qubits are being gated".
   - Waits then come in gate-layer units, not chosen durations.
   - t_R comes from the duration report, if Nexus provides one.
3. **Validate the fallback on the hosted emulator first.** It models idling by
   duration; local Selene does not.

## 4. Testable claims

**Cost** is my estimate: Low ≈ a few qubits and ≲ 10⁴ shots in total;
Medium ≈ 10⁴–10⁵; High ≈ 10⁵ or more per point.
**Wait?** = whether the test needs controllable idle time.

| # | Claim | From | Hardware test | Observable | Cost | Wait? |
|---|---|---|---|---|---|---|
| T1 | Readout flips only the record, not the qubit | 03, 04d | Measure a qubit mid-circuit, then again; Bell-pair disagreement as in 03 | Second readout agrees with the first; disagreement = 2p(1−p) | Low | No |
| T2 | Readout is asymmetric (0→1 ≠ 1→0) | 07 uses the average | Prepare \|0⟩ and \|1⟩, then measure | p_meas0, p_meas1 | Low | No |
| T3 | Idle error is Z only (no X/Y) | 07 | Idle \|0⟩, measure Z; idle \|+⟩, measure X | Z-basis flip ≈ 0; X-basis flip grows with wait | Low | Yes |
| T4 | Idle Z error = linear γ + time² c | 07 | Idle \|+⟩ for a log grid of waits (~0.1–30 t_R) | Fit P(t) = γt + sin²(√c·t) | Low–Med | Yes |
| T5 | How the time² term splits: does it restart at gates? | 09 (open) | Same total wait, as 1 block versus k blocks with a gate between | Error ∝ c·t² (accumulates) or c·t²/k (restarts) | Low–Med | Yes |
| T6 | H2's idle-to-gate ratio is q(1) ≈ 1.1·p2 | 07 §7 | Follows from T4 plus t_R | γ·t_R + c·t_R² against p2 | from T4 | Yes |
| T7 | Ancillas and data share one noise model; idling continues during measurement of neighbours | 08, 09 | Spectator in \|+⟩ while neighbours measure and reset, against a quiet wait of equal length | Extra spectator error; measurement crosstalk (not modelled) | Low–Med | Yes |
| T8 | A Z-basis memory does not see idle error; an X-basis one does | 07 §2 | d=3 memory in both bases, short versus long wait | Detection rate rises with wait only in X basis | Med | Yes |
| T9 | p_read is flat in Δt | 08, 09 (P2) | d=3 X-basis memory at Δt = 1 plus 1–2 longer waits; `timed_memory.edge_probabilities` | Time-like edge probability flat within error | Med (~10⁴ shots, estimate) | Yes |
| T10 | p_data(Δt) = p_data(1) ⊕ p_idle(wait) | 08, 09 (P2) | Same runs as T9 | Space-like (+ diagonal) edges against T4's p_idle | Med | Yes |
| T11 | b_read ≈ 1 from the gates alone at H2 ratios | 08, 09 | Same runs as T9, at Δt = 1 | p_read / p_data(1) | Med | No |
| T12 | The full error-model decoder beats the hand-built graph | 06b | Decode the T9 shots both ways | Logical failures, both decoders | free (offline) | — |
| T13 | The device compiles a CX as Selene does (1 Rzz + 4 single-qubit gates) | 06a, 07b | Inspect Nexus's compiled program, if it is exposed | Native gate list per round | Low | No |
| T14 | Δt\* follows from one calibration; back to back is best at H2 | 09 (P3), 07 §7b | p_L at Δt = 1 against Δt ≈ 2.3, d=3 X basis | p_L ratio; predicted effect 11–16% | **High** (10⁵+ per point) | Yes |

**Not testable on hardware, and why:**
- **04 / 04c / 05 / 06c:** λ from counting operations, and Selene's Δt\*.
  These describe Selene's per-operation noise model, not a device.
- **07b's native-versus-textbook comparison:** a device has only one
  decomposition.
- **09's p_L underestimate at d=7:** needs 73+ qubits and very large shot
  counts.

## 5. Order of work

0. **Settle the delay question** (section 3) and get t_R. This blocks T3–T10
   and T14.
1. **No-wait tests:** T1, T2, T11, T13. Cheap, and runnable now.
2. **Single-qubit idle:** T3 → T4 → T5, giving T6. This is the step that
   decides the H2 conclusions.
3. **Spectator idle:** T7.
4. **Small-code detection events:** T8 → T9, T10. T12 runs offline on the same
   shots.
5. **Conformance ladder.** Feed sections 2 / T1–T7 into each model and compare
   predicted edge probabilities with T9–T10:
   - **07:** textbook gates, idle on data only;
   - **07b:** native gates;
   - **08:** uniform noise on every qubit;
   - **09:** native gates and uniform noise together.

   The first model that matches is sufficient. Where a step fails, the missing
   ingredient is identified.
6. **Only then** compute Δt\* and the depth of the U at calibrated rates in
   Stim. Decide whether T14 (expensive) is worth running.

## 6. Expected outcomes (interpretation)

- **T4/T5 move the conclusions most.**
  - If γ ≪ 1.1·p2 / t_R, "back to back is best" gets stronger.
  - If the idle error is larger, or mostly time², 09's open quadratic case
    decides the answer.
- **If T9–T10 hold,** the two-number paper model calibrated on the device is
  sufficient. The p_L sweep can then stay in Stim, and T14 becomes a spot check
  rather than the experiment.
- **T7** may reveal measurement crosstalk, which no notebook models. If it is
  large, it is the next ingredient to add to 09.
