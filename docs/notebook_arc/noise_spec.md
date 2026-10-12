# Noise specification: the circuit-level models

The noise models of the notebook arc, in one place. The counterpart of
[`../stub_arc/noise_spec.md`](../stub_arc/noise_spec.md), with one difference in kind: there
the paper's model is the **contract**; here it is the **hypothesis**, and the circuit's noise is
the input. Symbols follow [`notation.md`](notation.md); paper equations are in
[`../paper_equations.md`](../paper_equations.md).

There are two families. **Selene models** (03–06c) run Selene's own `DepolarizingErrorModel`.
**Stim models** (07–09) build the circuit in Stim and add the idle noise Selene cannot express.

## 1. Scope

- **Circuit-level.** Every gate, measurement and reset carries its own fault; nothing is
  injected as a phenomenological flip. The paper's `p_data` and `p_read` are *measured*
  (section 6), never set.
- **Pauli channels only.** Stim samples Pauli channels; a coherent rotation is Pauli-twirled
  (section 3). No leakage, crosstalk, transport or asymmetric readout.
- **Memory basis follows the idle noise.** Selene models: Z basis (`|0_L⟩`, Z-checks). Stim
  models with dephasing: **X basis** (`|+_L⟩`, X-checks), because a Z-basis memory cannot see Z
  errors (07 §9 finding 1).

## 2. Selene models (notebooks 03–06c)

Selene's `DepolarizingErrorModel` charges noise **per native operation** and has **no duration
term** (03, `capability_audit`). What it offers:

| Parameter | Acts as | Placement (read off Selene's log, `selene_mirror`) |
| --- | --- | --- |
| `p_1q` | uniformly random Pauli | **before** every native single-qubit gate (`Rxy`, `Rz`) |
| `p_2q` | two-qubit depolarizing | on each `Rzz`; placement **by analogy, not verified** |
| `p_meas` | flip of the reported bit only (the paper's `p_read` form) | each measurement |
| `p_init` | — | not used (0): resets noiseless |

**Idling is written in as operations.** One Δt unit is `k` explicit `phased_x(q, 0, 0)` ops on
each data qubit, each charged `p_1q`. So idle noise is **depolarizing** (X, Y and Z), and `λ` is
a ratio of op counts, not a hardware constant:

```
p_2q = 0 (04):   λ = k / n_1                        p_1q cancels
p_2q > 0 (04c):  λ(ρ, k) = k / (n_1 + 0.8·ρ·n_2)    ρ = p_2q / p_1q
```

`n_1` and `n_2` are the charged single-qubit ops and `Rzz` events reaching a data qubit per
round, measured in 04/04c and stored in `config/surface_extraction_measured.yaml` and
`config/surface_extraction_p2q.yaml` (`extraction_constants.yaml` is the retired Steane version).

**Readout composition (04d).** Ancilla faults and `p_meas` combine as two independent flips,
`p_read = p_anc ⊕ p_meas`, with `a ⊕ b = a + b − 2ab`. 05 and 06c used `1 − (1−a)(1−b)`, which
is wrong by `p_anc·p_meas` (6×10⁻⁶ at the operating point; no conclusion changes). The **final
data readout** carries `p_meas` only, not `p_read`: the data are measured directly.

Operating point of 06c: `p_1q = 5×10⁻⁴`, `p_meas = 2×10⁻³`, total idle 96, rotated code.

## 3. The idle channel (Stim models, 07–09)

Form from Quantinuum's emulator noise model: memory error during idling is **Pauli Z only**,
two terms applied together. Implemented once, in `stim_memory.idle_channel`:

```
linear       Z with probability  γ·t                                  (capped at 1/2)
quadratic    RZ(κ·t), Pauli-twirled: Z with probability  sin²(κ·t/2)
combined     p_Z = p_lin ⊕ p_quad
```

The notebooks pass `κ = 2√c`, so the quadratic term is `sin²(√c·t) ≈ c·t²`, the form in
[`notation.md`](notation.md). `idle_channel` also accepts an `idle_depol` (X, Y, Z) term; 07
onward set it to 0, assuming the emulator docs' Z-only form.

**Time unit.** `t_R`, the duration of one extraction round. `τ` (written Δt in the notebooks) is
the **full round period**, so `τ ≥ 1` and the wait between rounds is `τ − 1`.

**Settings used throughout 08–09.** 07 and 07b sweep seven: `γ/p_2q` ∈ {0.15, 0.3, 0.6}, `c/p_2q` ∈ {0.0025, 0.01, 0.04}, and the mixed point, which 07 holds out of calibration.

| key | `γ` | `c` |
| --- | --- | --- |
| `lin 0.3` | 0.3·`p_2q` | 0 |
| `quad 0.01` | 0 | 0.01·`p_2q` |
| `mixed` | 0.3·`p_2q` | 0.01·`p_2q` |

### Where the idle channel is applied

| Model | Class | Inside the round | Over the wait | Ancillas |
| --- | --- | --- | --- | --- |
| 07 | `StimMemory` | — | `idle_channel(τ)` once per round, data only | none |
| 07b | `NativeMemory` | — | `idle_channel(τ)` once per round, data only | none |
| 08 | `TimedMemory` | `idle_channel(1/7)` after each of 7 `TICK` layers | `idle_channel(τ − 1)`, every qubit | same as data |
| 09 | `TimedNativeMemory` | `idle_channel(1/L)` after each of `L` ASAP layers (`L` = 18 at d = 3, 5, 7) | `idle_channel(τ − 1)`, every qubit | same as data |

`TimedMemory(ancilla_idle=False)` exempts the ancillas, as a diagnostic only (08 §5). On a
freshly reset ancilla a Z does nothing; that is a consequence, not a rule.

### Open: does the quadratic term restart?

`sin²(√c·t)` is not additive in time. Applied per layer (08, 09), one round costs about `c/L`;
applied once over an unbroken round (07) it costs about `c`. Which is right depends on whether
the coherent phase behind it is reset between layers on hardware. This decides the quadratic
H2 conclusion (09 finding 4) and is test T5 in [`hardware_tests.md`](hardware_tests.md). The
linear term adds up in time and has no such issue.

## 4. Gate, measurement and reset noise (Stim models)

| Parameter | Textbook circuit (`StimMemory`, `TimedMemory`) | Native circuit (`NativeMemory`, `TimedNativeMemory`) |
| --- | --- | --- |
| `p1` | `DEPOLARIZE1` after each `H` | `DEPOLARIZE1` on each `Rxy`; on `Rz` only if `rz_noise` (off in 07b, 09) |
| `p2` | `DEPOLARIZE2` after each `CX` | `DEPOLARIZE2` on each `Rzz` |
| `p_meas` | flip before each measurement | the record flipped |
| `p_init` | flip after each reset | `X_ERROR` after each reset |

The native gate sequence is Selene's compiled program (each CNOT = one `Rzz` + four single-qubit
gates, 06a), read off one noiseless shot and rebuilt gate for gate by `selene_mirror`. In the
native circuit `p2` is a modelling choice (H2's two-qubit error on each `Rzz`), not a copy of
Selene's `p_2q`.

**Rates.** H2-1E emulator sheet (2025-07-16), with readout averaged over its two directions:

| | `p1` | `p2` | `p_meas` | `p_init` |
| --- | --- | --- | --- | --- |
| H2-1E | 1.89×10⁻⁵ | 1.05×10⁻³ | 9.95×10⁻⁴ (avg of 6×10⁻⁴, 1.39×10⁻³) | 4×10⁻⁵ |
| **scaled ×4.76** (used in 07–09) | 9.0×10⁻⁵ | **5×10⁻³** | 4.74×10⁻³ | 1.9×10⁻⁴ |

The ×4.76 scale keeps H2's ratios and is there for measurability (07).
Real H2 rates are used only in the H2-ratio checks (07 §7b, 09 §6).

**H2's idle-to-gate ratio** is an estimate: idle error over one round `q(1) ≈ 1.1·p_2q`, from
the H2 data sheet's memory error per depth-1 circuit and an assumed 6 layers per round (07 §7).

## 5. Time bookkeeping

`T` is fixed (192 in 07–09; 96 in 05–06c) and the round count varies with Δt. Grids are built
from whole round counts, so the realised time equals `T` exactly:

```
dt_grid(d):  R = round(T / x) for x log-spaced in [1, T/d], keep R ≥ d;  Δt = T / R
```

Logical error per unit time follows Eq. A3 inverted: `R = −ln(1 − 2 p_L) / (2T)`. Points with
`p_L ≥ 0.45` are excluded from the fit. Δt\* is the vertex of a weighted parabola in
`log Δt` through the points within 30% of the minimum, with a bootstrap 16–84% band.

## 6. Reading the circuit as the paper's model

`timed_memory.edge_probabilities` measures the two phenomenological parameters from detection
events with the `p_ij` estimator (Google, Nature 595, 383 (2021), supplement), exact when every
error flips at most two detectors:

```
p_ij = 1/2 − 1/2 · sqrt( 1 − 4(⟨x_i x_j⟩ − ⟨x_i⟩⟨x_j⟩) / (1 − 2⟨x_i⟩ − 2⟨x_j⟩ + 4⟨x_i x_j⟩) )
```

| Measured quantity | Detector pairs, bulk rounds only |
| --- | --- |
| `p_read` | mean time-like edge: check `c` in rounds `r` and `r+1` |
| `p_data` | mean over check pairs sharing a data qubit of space-like ⊕ both diagonal edges |
| `p_stab` | `p_data` of the same circuit with idle noise off |
| `b_read` | `p_read / p_stab`, the paper's definition with both sides measured |

Diagonal edges are faults between two checks' CNOTs; they are 37–40% of the per-round data error
at Δt = 1 and fold into the per-round term (08 finding 2).

**Breakdown.** The estimator assumes rare, graph-like errors. When the data error per round
nears ~0.2 it returns negative edges; such intervals are dropped (one case: d = 3, mixed,
Δt = 64).

**The paper's model as a comparison circuit:** `timed_memory.phenomenological_circuit` — ideal
gates, one Z flip per data qubit per round at `p_data`, every recorded bit (syndrome and final
readout) flipped at `p_read`.

## 7. Decoding

Always on the circuit's own detector error model: `circuit.detector_error_model(
decompose_errors=True)` into PyMatching. A hand-built phenomenological graph misses mid-round and
propagated faults and misstates `p_L` by up to 7–8× (06b). Selene's own shots are decoded against
the error model of their Stim mirror (`selene_mirror.decode_mirror`).

## 8. Randomness and seeds

- Stim runs: seeds 1 and 2, 10⁵ shots each, per point (07–09). Edge measurements: seed 1,
  4×10⁵ shots, 20 rounds.
- **Selene seeds overlap.** Selene seeds shot `k` with `random_seed + shot_offset + k`, so runs
  with overlapping seed ranges share their noise (04d §1b). 06c's `p_meas` sweep ran every
  `p_meas` at seed 1.
