# Stub to device: running the injected-noise memory on Quantinuum hardware

Status: **plan, nothing run on hardware yet.** Scope is the stub architecture
(`src/qec_timing/stub_circuit/`) and a hardware version of the same experiment.
Numbers in sections 5-6 are from the reference path, which is verified
shot-for-shot against the stub (Stage 4a), at the paper's injected parameters.

## 1. What this experiment can and cannot show

The stub reproduces the paper's model: perfect gates, and all noise injected
in the program. The hardware version keeps the injection and runs on a real
device.

- **It can show:** the *injected* Δt trade-off reproduced on a device. The
  device's own errors add a background that does not depend on Δt.
- **It cannot show:** the device's *own* idle trade-off. Δt changes only the
  injected probability. The program never waits, so every round lasts the same
  physical time whatever Δt is.

## 2. The stub today

| Item | Implementation |
|---|---|
| Model | X (bit-flip) only, memory in \|0_L>, Z-checks only |
| p_idle, p_data, p_read | Computed in Python (`noise/params.py`); passed to `run()` as `float` arguments |
| Data noise | `apply_data_noise`: one `rng.random_float() < p_data` draw per data qubit per round, then `x()`. Applied before extraction |
| Readout noise | `flip_record` on the classical bit from every ancilla measurement **and** every final data measurement, at one `p_read` |
| Randomness | One `RNG(base_seed + get_current_shot())`, `rng.discard()` at the end |
| Circuit | d² data qubits + **1 reused ancilla**; Z-checks measured one after another |
| Gates | Written as `std.quantum` `cx` / `x` / `measure`. Compiled for `platform="helios"` to native `Rzz`, `Rxy`, `Rz`, `Reset`, `MeasureRequest` (one CX = 1 `Rzz` + 4 single-qubit gates). Simulated with no error, so they have no effect on results |

## 3. What changes on a device

**The native gates become a real noise source.** They contribute to the same
model parameters, but not exactly in the model's form:

1. **Not uniform.** A data qubit sits in 1 or 2 Z-checks, so it gets 1 or 2 CX
   per round.
2. **Split between parameters.** Gate errors on data qubits act as extra p_stab.
   Errors on ancillas (up to 4 CX each) act as extra p_read.
3. **Only part of each error counts.** Only the X and Y parts flip a data qubit
   in a Z-basis memory.
4. **Correlated.** One `Rzz` fault can hit a data qubit and its ancilla together.

**Rough size of the background.** Estimated from the H2-1E emulator rates, not
measured: p_1q = 1.9e-5, p_2q = 1.05e-3, readout ≈ 1e-3.

| Model term | Device background | Relative to injected p = 0.015 |
|---|---|---|
| p_stab | ≤ 2 × (8/15) × p_2q ≈ 1e-3 | ~7% |
| p_read | readout + ≤ 4 CX on the ancilla ≈ 3e-3 | ~20% |
| idle | mostly phase (Z), which does not show in a Z-basis memory | small |
| hook errors | none: data is the CX control, ancilla the target | none |

## 4. Steps

In order. Step 1 decides whether this route works at all.

1. **Confirm randomness on hardware.** Check that `qsystem.random.RNG` and
   `get_current_shot()` are supported, and that the shot number changes from
   shot to shot. Use a small test job, e.g. a few qubits, emitting
   `random_float()` and the shot number. If the shot number is constant, every
   shot gets identical noise and all statistics are invalid.
2. **Fix parameters at compile time.** Hardware does not accept `run()`
   arguments (guppylang `emulator/__init__.py` docstring). Pass `p_data`,
   `p_read`, `n_rounds` and `base_seed` as `comptime` values, which means one
   compiled program per sweep point.
3. **Parallel Z-ancillas.** Allocate (d²−1)/2 ancillas per round instead of
   reusing one. Schedule the CXs in **4 layers** so no data qubit is in two CXs
   at once. Keep data as control and ancilla as target. Measure all ancillas,
   then compute detection events as now. The detector order must stay identical
   to the reference (re-run the single-fault injection test from Stage 3).
4. **Submission and results.** Add Nexus/`qnexus` submission, and an adapter
   from Nexus results to the `collate_tags()` form that `detection_events()`
   expects. Repeated tags must stay in order.
5. **Leakage.** Flag leaked qubits with `measure_leaked` (in the Helios gate
   set) on the final readout at least. Discard or count flagged shots.
6. **Pre-hardware check.** Run the hardware version locally under Selene's
   `DepolarizingErrorModel` at H2-like rates, or on the hosted emulator. Confirm
   it builds, the detector order is correct, and the curve shifts as section 3
   predicts.
7. **Background calibration on the device.** Run with zero injection
   (`p_data = p_read = 0`) at the production d and a few round counts. Fit the
   background p_stab_bg and p_read_bg from detection-event rates.
8. **Set the injected amounts** so the totals hit the target:
   `p_inj = 1 − (1 − p_target) / (1 − p_bg)` for each of p_stab and p_read.
   p_idle is injected as before. Re-run the section 6 scan on the reference
   path with the totals, and adjust the grid if the minimum has moved.
9. **Decoder weights.** Use injected + measured background, not the injected
   values alone.
10. **Production run** on the section 7 grid. Decode on the reference
    machinery, and compare with the reference at the same n_rounds and totals.

## 5. Resource counts (Z-only, parallel ancillas)

Rotated code: d² data qubits, (d²−1)/2 Z-checks: (d−1)²/2 of weight 4 and
d−1 of weight 2.

| d | Qubits (3d²−1)/2 | Fits H2 (56)? | CX = Rzz per round, 2d(d−1) | Measurements per round | RNG draws per round |
|---|---|---|---|---|---|
| 3 | 13 | yes | 12 | 4 | 13 |
| 5 | 37 | yes | 40 | 12 | 37 |
| 7 | 73 | no | 84 | 24 | 73 |

Every shot also has d² final data measurements. The measurement count per
round is the same whether ancillas are reused or parallel. Parallel only
shortens the round. How much of the parallelism the device really achieves
depends on its gate and measurement zones and on the compiler. This is not
verified.

## 6. Expected behaviour (reference path, injected parameters only)

Settings: p = 0.015, λ = 1, b_read = 1. 200k shots per point.
R = −ln(1 − 2 p_L) / (2 T_actual) is the logical error rate per unit time.

**d = 5, T = 10.** Clear U-shape, minimum near Δt ≈ 1:

| Δt | 0.25 | 0.35 | 0.5 | 0.7 | 1.0 | 1.4 | 2.0 | 2.8 | 4.0 | 5.0 |
|---|---|---|---|---|---|---|---|---|---|---|
| rounds | 40 | 29 | 20 | 14 | 10 | 7 | 5 | 4 | 2 | 2 |
| p_L | 0.187 | 0.165 | 0.144 | 0.132 | 0.132 | 0.134 | 0.146 | 0.178 | 0.149 | 0.199 |
| R × 100 | 2.34 | 1.97 | 1.70 | 1.56 | 1.54 | 1.59 | 1.73 | 1.97 | 2.20 | 2.54 |

**d = 3, T = 10.** Minimum near Δt ≈ 2, but the curve to its right is almost
flat:

| Δt | 0.25 | 0.5 | 1.0 | 1.4 | 2.0 | 2.8 | 4.0 | 5.0 |
|---|---|---|---|---|---|---|---|---|
| rounds | 40 | 20 | 10 | 7 | 5 | 4 | 2 | 2 |
| R × 100 | 4.91 | 3.20 | 2.44 | 2.26 | 2.17 | 2.22 | 2.32 | 2.45 |

R rises only 7–13% from Δt = 2 to 5, so resolving the right side needs
percent-level error bars (~10⁵ shots per point). **d = 3 is a smoke test, not
the experiment.**

**Choice of T.**
- **T = 10.** At d=5, R at T = 10 and T = 20 agrees within ~2% at every Δt, as
  the ansatz predicts (Δt* does not depend on T), so T = 20 only doubles the
  depth.
- **T = 5 is too short.** It leaves 1–2 rounds at large Δt, and rounding causes
  jitter (Δt = 2.8 rounds to 2 rounds, so T_actual = 5.6).
- **Grid points:** choose Δt = T/n with a whole number of rounds n, so that
  T_actual = T exactly.

## 7. Shot allocation: d = 5, T = 10

Reference values at 400k shots per point. Shots are sized for ±5% on R.

| n rounds | Δt | p_L | R / R_min | Shots | Rzz per shot |
|---|---|---|---|---|---|
| 30 | 0.333 | 0.164 | 1.29 | 3,100 | 1,200 |
| 14 | 0.714 | 0.134 | 1.01 | 3,600 | 560 |
| 10 | 1.000 | 0.133 | 1.00 | 3,600 | 400 |
| 7 | 1.429 | 0.135 | 1.02 | 3,500 | 280 |
| 4 | 2.500 | 0.155 | 1.20 | 3,200 | 160 |
| 3 | 3.333 | 0.169 | 1.34 | 3,000 | 120 |

- **Total: ~20,000 shots, ~9 × 10⁶ `Rzz`.** The n = 30 point alone is ~40% of
  the `Rzz`.
- **The ends carry the result.** They are 20–34% above the minimum, about
  4–7σ at ±5%. This shows that a U-shape exists.
- **The bottom is flat.** The three middle points are within 2% of each other.
  Pinning Δt* down needs roughly 4× more shots there (estimate, not fitted).
- **Calibration (step 7) and the d = 3 smoke test** are extra. They are not
  sized here.
- **Cost per shot in hardware credits** is not estimated.

Shot formula: relative error on R ≈ k · √((1 − p_L) / (N p_L)), with
k = 2 p_L / ((1 − 2 p_L) · (−ln(1 − 2 p_L))).

## 8. How much gate error the effect tolerates

**Question.** At what ratio of gate error to injected noise does the device
background wash out the U-shape?

**Method.** A Stim model of the planned hardware circuit, using the same
parameters as section 7:
- d = 5, parallel Z-ancillas, CXs in 4 layers;
- `DEPOLARIZE2(p2)` after every CX, readout error p_meas = 0.95 · p2 (H2's
  ratio);
- the stub's injection on top: p = 0.015, λ = 1, T = 10, n = 30, 20, 14, 10,
  7, 5, 4, 3;
- decoded with pymatching on the circuit's exact error model, 150k shots per
  point;
- Δt* from a parabola fit of log R against log Δt around the lowest point.

Sanity check: at p2 = 0 it reproduces the reference (R_min = 0.01538, against
0.01541 in section 6).

Contrast below is R at each end of the grid divided by R_min.

**Uncompensated: inject p = 0.015, background on top.**

| p2 / p_inj | p_L at min | Δt* | Contrast, n=30 end | Contrast, n=3 end |
|---|---|---|---|---|
| 0 | 0.132 | 0.98 | 1.31 | 1.35 |
| 0.07 (H2) | 0.154 | 1.24 | 1.40 | 1.23 |
| 0.14 | 0.173 | 1.34 | 1.45 | 1.15 |
| 0.28 | 0.203 | 1.75 | 1.63 | 1.07 |
| 0.56 | 0.247 | 2.45 | 1.90 | **1.01: right half gone** |

**Compensated (step 8): injected amounts lowered so the average totals stay at
p = 0.015.**

| p2 / p_inj | p_L at min | Δt* | Contrast, n=30 end | Contrast, n=3 end |
|---|---|---|---|---|
| 0 | 0.132 | 0.98 | 1.31 | 1.35 |
| 0.07 (H2) | 0.136 | 1.07 | 1.33 | 1.29 |
| 0.14 | 0.139 | 1.05 | 1.34 | 1.27 |
| 0.28 | 0.141 | 1.13 | 1.29 | 1.22 |
| 0.56 | — | — | — | — |

At 0.56, compensation is **infeasible**: the background readout error already
exceeds the target p_read.

**What "washing out" means.** The background adds a cost every round that
does not depend on Δt. So the optimum moves to longer Δt. Uncompensated, it
leaves the right edge of the grid, and the right half of the U goes flat.
Compensation holds the totals at the model's values, and the U keeps its shape.

**Feasibility region** (with p_meas ≈ p2):

- **Compensated: p2 / p_inj ≲ 0.36.** Readout sets the bound. An ancilla's
  background readout error is about p_meas + (average check weight) × (8/15) ·
  p2 ≈ 2.7 · p2 at d = 5 (average weight 3.3). It must stay below
  p_read = b_read · p_inj. Inside the bound, Δt* moves by ≲ 15%. That residual
  comes from the four non-model effects in section 3.
- **Uncompensated: p2 / p_inj ≲ 0.15** keeps ≥ 15% contrast on the right side
  within this grid. Beyond that the grid must move to longer Δt.
- **H2 at p = 0.015 sits at 0.07,** inside both regions. Uncompensated, Δt*
  shifts about 27%. Compensated, it shifts about 9%.

**Consequence for the plan.** Compensation (steps 7–8) is what makes the run
faithful to the model. Without it, the run still shows a U at H2's ratio, but
at a shifted Δt*. Lowering the injected p to save depth or shots moves the
ratio up: at H2's p2 the compensated bound is reached near p_inj ≈ 3 × 10⁻³.

**Limits of this check.**
- The Δt* values have no error bars; at 150k shots, the 5–15% compensated
  shifts may be partly noise.
- The decoder is given the exact noise model; on hardware the weights come from
  calibration.
- Not included: single-qubit gate, reset and leakage errors.
- d = 5 only.

## 9. Open questions and caveats

- **Randomness on hardware** (step 1) is unverified, and blocks everything else.
- **Section 6–7 numbers use injected noise only.** The device background shifts
  them. Redo the scan after step 7.
- **The large-Δt points have 3–4 rounds, fewer than d.** Comparing the device
  against the reference at the same n is still valid, because both see the same
  edge effects. But R there is not the steady per-round rate.
- **The 4 background effects in section 3** (non-uniform, split, partial,
  correlated) are not in the model. They limit how close device and reference
  can agree.
- **The hardware compiler's gate sequence and scheduling** may differ from
  Selene's.
- **Machine:** which machine is available (H2 or Helios) decides whether d = 7
  is possible (73 qubits).

## Provenance

- **Reference scans:** run 2026-10-12 with `qec_timing.reference`
  (`build_memory_circuit`, `build_decoder`, `estimate_logical_error`) under
  stim 1.16.0 and PyMatching 2.4.0.
- **Native-gate lowering:** checked the same day by logging the operations of a
  `qubit` / `x` / `cx` / `measure` program through `.emulator().stabilizer_sim()`
  under guppylang 1.0.3 / selene-sim 0.3.0 (the stub itself did not build in
  that environment; it needs guppylang 1.1.1).
- **Background tolerance (section 8):** a standalone Stim + PyMatching model of
  the section 4 step-3 circuit, run on the same day with the same versions. It
  is not yet a repo script.
- **H2-1E rates:** Quantinuum H2 emulator parameters, 2025-07-16.
