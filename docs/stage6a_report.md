# Stage 6a report: circuit-level Stim twin, built and tested

Scope: **build and test** a circuit-level Stim model of the rotated surface-code
memory with a switchable idle channel. No sweeps, no conclusions about dt\*.
This stage goes beyond the CLAUDE.md roadmap (team-approved) and deliberately
departs from the Stage-1 physics, which has no gate noise. Stage-1 code
(`reference/`, `noise/`, `stub_circuit/`) and `docs/interfaces.md` are untouched.

Versions: **stim 1.16.0**, **PyMatching 2.4.0** (Python 3.14.0, win_amd64).
Every Stim API used was checked in the installed package first.

## What was built

| File | Contents |
| --- | --- |
| `src/qec_timing/circuit_level/noise.py` | noise settings, idle channels, depolarizing composition, native op counts |
| `src/qec_timing/circuit_level/memory.py` | circuit builder, sampler, decoder, single-channel probe circuits, rate per unit time |
| `src/qec_timing/circuit/ordering.py` | selectable CNOT order, shared by the Guppy program and the twin |
| `src/qec_timing/circuit/program.py` | one keyword-only argument `order="row_major"`; nothing else changed |
| `tests/test_stage6a_circuit_level.py`, `tests/test_circuit_ordering.py` | 75 tests |

Noise model:

- **Gates:** `DEPOLARIZE1(p1)` after each single-qubit gate, `DEPOLARIZE2(p2)` after each CNOT.
- **Reset:** `X_ERROR(p_init)` after every reset.
- **Readout:** a classical flip of every record (ancilla and final data), possibly asymmetric.
- **Idle:** applied to every data qubit once per round, for the round's `dt`, before extraction.

| Idle variant | Channel | Default memory basis |
| --- | --- | --- |
| `depolarizing_kops` (A) | one `DEPOLARIZE1(p_n)`, `p_n = 3/4 (1 - (1 - 4p/3)^(k dt))` | Z |
| `dephasing_h2` (B) | `Z_ERROR(p_z)`, `p_z = sin^2(0.5 f_quad dt) + f_lin dt`, clipped to [0, 0.5] | X |
| `bitflip` (bridge only) | `X_ERROR(p)` per round | Z |

Other details:

- **Rounds:** counted with `qec_timing.schedule` (`n_rounds = round(T/dt)`).
- **Decoder:** built from `detector_error_model(decompose_errors=True)` through the existing `Decoder` / `PyMatchingDecoder`.
- **Figure of merit:** the Stage-4a rate per unit time, `R = -ln(1 - 2 p_L) / (2 T_actual)`.

## The extraction order, and why it is now selectable

The twin follows `circuit/program.py` qubit for qubit. Checks run one at a time
on one reused ancilla: every active check `ci`, then the complementary check
`ci`, with both measured (noisily) every round. The program's original order
is NW, NE, SW, SE for both check types. This is not Stim's generated ordering.

**That order costs the X memory distance under gate noise.** If a Z-check
ancilla fails after two of its four CNOTs, it leaves a horizontal pair of Z
errors on the data. That pair runs along the logical Z (a row). Measured graphlike
circuit distance (`shortest_graphlike_error`, gate noise on, both variants,
with and without native charging):

| order | Z memory d=3/5/7 | X memory d=3/5/7 |
| --- | --- | --- |
| `row_major` (default, unchanged) | 3 / 5 / 7 | **2 / 3 / 4** |
| `hook_safe` | 3 / 5 / 7 | **3 / 5 / 7** |

Variant B needs the X memory, so with the default order its d-scaling would be
distorted. `hook_safe` uses the standard shapes: N (NW, SW, NE, SE) for Z-checks,
so their Z hooks are vertical, and Z (NW, NE, SW, SE) for X-checks, so their X
hooks are horizontal. The shapes are derived from the layout's logical supports.
Weight-2 boundary checks keep the order of the two time slots they occupy, which
matches the layout's order in both modes. **The default is unchanged.** The
support tables are byte-identical to the pre-change code for d = 3, 5, 7, 9, 27
(sha256 snapshot taken before the edit), and the existing circuit, Stage 5 and
Stage 1 suites pass unmodified (87 tests). Recommendation for 6b: use
`order="hook_safe"` for any X-memory run.

The twin and the Guppy builder go through the same `order_supports` helper. A
test also checks the twin against Selene's actual execution: one noiseless
traced shot per (basis, order). Each CNOT compiles to one `ZZPhase` on
(ancilla, data), and the executed sequence matches the twin's CNOT sequence
qubit for qubit in all four cases.

Two deliberate differences from `program.py`:

- In the X memory, `program.py` wraps its per-round X-noise in H gates. That is a noise-injection trick, so those H gates are not in the twin.
- The preparation corrections are conditional gates. Stim cannot charge noise on a gate that only sometimes runs, so they are noiseless here. They only act on the checks that are not decoded.

## Native charging (`native_charging=True`, default off)

The counts were **taken from committed notebook outputs, not measured here**.
No Selene runs were allowed for this, and `config/*.yaml` is not in the repo.

| Op | Charged as | Source |
| --- | --- | --- |
| `h` | 2 single-qubit ops -> one `DEPOLARIZE1(compose(p1, 2))` | notebook 03, section 2 (k = 2 fits within 1.6 sigma; k = 1, 3 miss by 33, 21 sigma) |
| `cx` | one `Rzz` (`p2`) + 2 single-qubit ops per qubit -> `DEPOLARIZE1(p1)`, CX, `DEPOLARIZE2(p2)`, `DEPOLARIZE1(p1)` | notebook 04c, section 2 (`CircuitExtractor` read-back) |
| idle op | 1 single-qubit op (`phased_x(q, 0, 0)`) | notebook 03, section 1 |

Merging and moving these channels is exact, because single- and two-qubit
depolarizing noise is unchanged by local Cliffords. **Two assumptions need
checking in 6b:**

- Notebook 04c only read back data qubits, so 2 ops are assumed at the ancilla end too.
- The 2 single-qubit ops are split as one before and one after the `Rzz`.

The Selene trace above already exposes the compiled sequence (`PhasedX`, `Rz`,
`ZZPhase`), so both can be checked directly in 6b.

## Approximations

- **Asymmetric readout.** Stim's `M(p)` flips a record with the same probability whatever the outcome, and records combine only by XOR. When `p01 != p10`, the twin samples with perfect readout, flips each recorded bit in numpy with `p01` or `p10` according to its value (**exact**), then rebuilds detectors with `compile_m2d_converter`. Only the decoder's weights are approximate: they use `(p01 + p10)/2`. Side effect: the preparation corrections see the record before the flip, which only touches the checks that are not decoded.
- **Non-integer k dt.** The merged probability above is exact for whole numbers of ops. In between, it is the continuous-time depolarizing channel that agrees at every integer. Selene can only charge whole ops, so that rounding is a 6b question.

## Tests (all fixed seeds; 99.9% Clopper-Pearson intervals)

1. **Noiseless:** zero detections and zero logical errors, both bases, both variants, with native charging on and off. The complementary checks also read +1 every round, which confirms the preparation step.
2. **Stage 4a bridge** (`test_reduces_to_phenomenological_reference_stage4a_bridge`): with gate and init noise off and the idle channel set to `X_ERROR(p_data(p, lam, dt))`, the detector error model has exactly the same error mechanisms as `qec_timing.reference`. Sampled p_L agrees within 4 sigma at 4 (d, p, dt) points, one of them with readout noise on.
3. **Idle rates:** variant A flips in both bases at 2/3 p_n, including k dt = 1.5 and 6.9. Variant B applies Z at p_z(dt) and never X. The merging formula exactly matches repeated `DEPOLARIZE1` in Stim's error model.
4. **Basis check:** with only dephasing noise, the Z memory shows zero detections and zero failures, while the X memory has p_L > 0.05. With gate noise on, the Z memory's error model is identical with and without dephasing.
5. **Asymmetric readout:** the flip rates for true 0 and true 1 each fall within their intervals. Inside a memory run, a flip that only hits true 1s shows up only in the final data readout.
6. **Order:** byte-identical default tables; hook_safe shapes; full hook_safe distance (d = 3, 5, 7 x Z/X x A/B x native on/off); the twin matches the traced Guppy circuit. `test_program_py_order_reduces_x_memory_circuit_distance` is kept as is.
7. **Smoke run**, d = 3, T = 6, default order, 4,000 shots. **No conclusions are drawn from these numbers:**

| variant | basis | dt = 0.5 | dt = 1.0 | dt = 2.0 |
| --- | --- | --- | --- | --- |
| depolarizing_kops (k=2, p1=1e-3, p2=5e-3) | Z | p_L 0.0377, R 0.0065 | 0.0192, 0.0033 | 0.0143, 0.0024 |
| dephasing_h2 (f_quad=0.1, f_lin=1e-3) | X | p_L 0.256, R 0.060 | 0.202, 0.043 | 0.180, 0.037 |

## Open for 6b

- Use `order="hook_safe"` for any X-memory comparison.
- Check the two native-charging assumptions against the Selene trace.
- Decide how Selene rounds a non-integer k dt.
