# Stage 7a report: H2-Emulator feasibility (no full build)

Branch `stage7-h2`. Target: the H2-Emulator on Nexus via Guppy 0.x -> HUGR ->
QIR. Nothing was submitted to Nexus: no upload, no costing job, no emulator
run. The main `.venv`, `circuit/program.py` and all validated code are
untouched. New code lives in `src/qec_timing/h2/` (the main venv never imports
the package) and in the `scripts/stage7a_*.py` files.

## 0. Environment, and one correction to the premise

`.venv-h2` (Python 3.11.9; full pins in `requirements-h2.txt`):

| package | version |
| --- | --- |
| guppylang | 0.21.11 |
| guppylang-internals | 0.32.0 |
| hugr | 0.15.5 |
| hugr-qir | 0.1.2 |
| pyqir | 0.12.6 |
| quantinuum-qircheck | 0.7.0 |
| qnexus | 0.51.1 |
| quantinuum-schemas | 7.8.2 |
| pytket | 2.18.5 |

**"hugr-qir is not compatible with Guppy v1" is out of date.** The PyPI release
history shows the dependency switched in late August:

- hugr-qir 0.1.2 (2026-07-02) is the **last** release for `guppylang ~= 0.21.11`.
- hugr-qir **0.2.0 (2026-08-27) and later require `guppylang >= 1.0.1, < 2`**. The current release, 0.3.0, also needs Python >= 3.12.

I followed the instruction and pinned the valid 0.x pair (0.21.11 + 0.1.2).
A Guppy-1.x route to QIR now exists, though. Whether to move to it is your
call. It would bring the H2 program onto the same Guppy major version as the
validated code, at the price of re-checking the conversion limits on 0.3.0.

Every API used was checked in the installed source. Where it differs from the
brief:

- **Cost:** `qnx.circuits.cost` only takes pytket `CircuitRef`s. For QIR the call is `qnx.qir.cost`.
- **RNG:** `RandomFloat` is not lowered to QIR (section 3).

## 1. H2 idle model and dt*

### Unit convention

`src/qec_timing/h2/noise_params.py` uses the H2-1E defaults (H2-1E table dated 2025-07-16):

- gates and SPAM: p1 = 1.89e-5, p2 = 1.05e-3, p_meas = 6e-4 (0 outcome) / 1.39e-3 (1 outcome), p_init = 4e-5;
- dephasing: quadratic 0.043, linear 0.0028, coherent_to_incoherent_factor 2.

The noise-model page states, verbatim:

- **Linear:** "probability of linear_dephasing_rate multiplied by duration ... where linear_dephasing_rate is per second (s^-1), and duration is in units of seconds." So **p_lin = r·d**, with d in seconds.
- **Incoherent quadratic** (the stabilizer-simulator default): "probability, sin(fd/2)^2, where frequency is equal to quadratic_dephasing_rate multiplied by the coherent_to_incoherent_factor (all in units of 2 pi radians per second)." So the factor scales the **frequency**, not the probability.

The one ambiguity is "units of 2π radians per second":

- **"cycles" (the default):** read literally, 0.043 is in cycles per second. The angle is then 2π·0.043·2·d.
- **"radians":** read as a loose label, the printed formula takes f = 0.086 rad/s directly.

The quadratic terms differ by (2π)² ≈ 39. I default to "cycles" because it follows the units statement literally, and it is also the conservative (noisier) choice. Both readings are implemented and both are reported below.

**Timing convention:** dt is the syndrome **period**, so dt ≥ t_round. Data qubits are charged memory noise for the whole period, and the emulator applies it during both idling and transport. In twin units (time unit = t_round = 2 ms by default, a parameter):

- f_quad = ω·t_round: 1.081e-3 rad per round ("cycles") or 1.72e-4 ("radians");
- f_lin = r·t_round = 5.6e-6 per round;
- p_z for one 2 ms round: 5.9e-6 ("cycles") or 5.6e-6 ("radians"). Both are about 200× below p2.

The twin's `DephasingH2(f_quad, f_lin)` adds the two terms. My H2 helper composes them as independent Z channels. The twin is therefore higher by exactly 2·q·l: about 3e-7 at the d=3 dt*, and 0.5% of p_z at the far end of the grid. Tests pin this in `tests/test_stage7a_h2_params.py`.

### dt* with the Stage 6c twin

Configuration: variant B, X memory, `hook_safe`, native charging, H2-1E p1 / p2 / asymmetric p_meas / p_init, and 16 log-spaced periods from 1 round upward. Each point targets 400 failures and has its own seed block. dt* comes from the 6c per-d shape fit with a 300-sample binomial bootstrap (95%). Run: `results/stage7a/20261010T025251Z-dtstar/dtstar.json`.

| convention | d | dt* (rounds) [95%] | dt* (t_round = 2 ms) | R(every round) / R(dt*) |
| --- | --- | --- | --- | --- |
| cycles (default) | 3 | **44.1** [41.9, 46.2] | **88 ms** | 23 |
| cycles (default) | 5 | **35.9** [34.7, 37.2] | **72 ms** | 20 |
| radians | 3 | 181 [173, 190] | 362 ms | 85 |
| radians | 5 | 138 [132, 144] | 275 ms | 70 |

All four minima are interior. The points flagged as saturated sit well to the right of each minimum and are excluded from the fit.

Readings:

- **Idle on H2 is so much weaker than a round's gate error that the optimum is to extract rarely**: tens of rounds' worth of time between extractions, with a 20× (or larger) gain over extracting back to back.
- **The 2π ambiguity moves dt\* by 4×.** That is less than the (2π)¹ expected from the quadratic law alone, because the linear term also contributes. It has to be settled, either by asking Quantinuum or by measuring it on the emulator, before any H2 dt* claim.
- dt* falls from d=3 to d=5, as 6c predicts for quadratic idle.

These are twin estimates with these approximations:

- **Rz:** native charging charges p1 on every Rz, while H2's Rz is virtual (error-free). The effect is small, since p1 = 1.89e-5.
- **Not modelled:** emission, crosstalk, leakage, and ancilla dephasing during extraction.
- **t_round:** dt* in rounds was computed at t_round = 2 ms. dt* in seconds changes with t_round, because the gate error per round is fixed.

## 2. Minimal d=3 test program

`src/qec_timing/h2/program.py`, `build_program(d=3, n_rounds=2, p_inject=0.05, base_seed=7_000_000, checks="Z")`:

- **Shape:** the entrypoint is `@guppy.comptime` (traced), so every Python loop unrolls into straight-line code. Qubits live in a Python list, not a Guppy array. There are no array-typed functions, and the entrypoint takes no arguments, so every parameter is compile-time.
- **Qubits:** 9 data qubits and 1 ancilla, allocated once. The ancilla is reused via `project_z` + `reset`.
- **CNOT order:** hook_safe, from `src/qec_timing/h2/geometry.py`. It re-implements `ordering.py`'s logic rather than importing it. A main-venv test (loading the file by path) shows it matches `layout` + `order_supports(..., "hook_safe")` for both check types at d = 3-11.
- **Seeding:** one `RNG(base_seed + get_current_shot())`, discarded at the end.
- **Injection:** each round, before extraction, Z on each data qubit with compile-time probability, inside a small non-traced `@guppy` helper, because traced code cannot branch on a runtime value.
- **Outputs:** `result("s{round}_{check}", bit)` for syndromes and `result("m{i}", bit)` for the final data readout.

**The spec as written cannot see its own injected noise.** It extracts Z-checks while injecting Z errors, and Z-checks are blind to Z: the same silent-blindness case as Stage 6a test 4. `checks="Z"` (the default) implements the spec exactly; its syndromes will carry only the emulator's own noise. `checks="X"` (data in |+⟩, X-check extraction) does detect the injected Z. Use `checks="X"` for any run whose syndromes should mean something.

## 3. QIR conversion (`hugr_to_qir`, validation on)

| what | outcome |
| --- | --- |
| Guppy -> HUGR | works |
| HUGR -> QIR, first attempt | **failed**: `Unknown op: RandomFloat`. Strings in the hugr-qir 0.1.2 binary show it lowers `RandomInt`, `RandomIntBounded`, `NewRNGContext` and `GetCurrentShot`, but not `RandomFloat`. |
| workaround | integer Bernoulli: Z if `random_int_bounded(2^20) < round(p·2^20)` (CLAUDE.md's named fallback). The effective p is recorded (0.05 -> 0.0500002). |
| d=3, 2 rounds, Z and X | **converts and validates**. Z: 8.1 KB bitcode; X: 9.0 KB. 0.3 s. |
| larger programs | **failed in validation**: `quantinuum_qircheck` recurses once per basic block, and every conditional Z adds blocks. It hits Python's recursion limit beyond about 600 injections (d=5 with 32 rounds). |
| workaround | validation kept **on**: recursion limit raised to 10^6, conversion run in a worker thread with a 64 MiB stack. (Windows rejects stack sizes of 256 MiB and above.) |

The d=3 QIR uses only native gates (`rz`, `phasedx`, `rzz`) plus `mz`, `reset`, the RNG and shot intrinsics, and `bool_record_output`. The counts are as expected: 24 `rzz`, 18 RNG draws (9 data × 2 rounds) and 17 outputs (8 syndromes + 9 data).
Files: `results/stage7a/d3_r2_{Z,X}.{bc,ll,json}`.

## 4. Commands (nothing below has been run)

Log in once, then build. Cost runs by default; the emulator runs only with `--submit`:

```
.venv-h2\Scripts\qnx.exe login
.venv-h2\Scripts\python.exe scripts\stage7a_build_qir.py --checks X
.venv-h2\Scripts\python.exe scripts\stage7a_submit.py results\stage7a\d3_r2_X.bc
.venv-h2\Scripts\python.exe scripts\stage7a_submit.py results\stage7a\d3_r2_X.bc --submit --shots 100 --wait
```

- **Cost:** the first submit command uploads the QIR and runs `qnx.qir.cost(..., system_name="H2-1")`, which is a costing job on the **H2-1SC syntax checker**. That job shows up in the Nexus portal. The command then prints the HQC estimate.
- **Run:** `--submit` runs on `H2-1E` with `QuantinuumConfig(simulator="stabilizer", noisy_simulation=True, error_params=UserErrorParams(idle_dephasing=False))`. The backend config was validated offline.
- **Transport dephasing:** it stays at the emulator default unless you pass `--transport-dephasing-off`.
- **Shot cap:** shot counts above 1000 are refused unless `--max-shots` is raised.

## 5. What will not scale (measured; `results/stage7a/20261010T030607Z-qir-scaling/`)

| d | rounds | qubits | HUGR | QIR bitcode | QIR calls | HUGR + QIR time |
| --- | --- | --- | --- | --- | --- | --- |
| 3 | 2 | 10 | 14 KB | 9 KB | 302 | 0.3 s |
| 3 | 8 | 10 | 27 KB | 23 KB | 1,010 | 0.9 s |
| 3 | 32 | 10 | 83 KB | 80 KB | 3,842 | 8.7 s |
| 3 | 64 | 10 | 151 KB | 168 KB | 7,618 | 19.9 s |
| 5 | 2 | 26 | 26 KB | 21 KB | 918 | 2.2 s |
| 5 | 8 | 26 | 67 KB | 64 KB | 3,138 | 6.7 s |
| 5 | 32 | 26 | 224 KB | 264 KB | 12,018 | 33.4 s |
| 7 | 2 | 50 | 46 KB | 39 KB | 1,862 | 4.3 s |
| 7 | 8 | 50 | 125 KB | 137 KB | 6,410 | 14.7 s |
| 7 | 16 | 50 | 228 KB | 273 KB | 12,474 | 34.2 s |

1. **Qubits cap d at 7.** H2-1E's stabilizer simulator holds 56 qubits. With one reused ancilla, d=7 needs 50 qubits and d=9 needs 82.
2. **Program size and compile time grow linearly in d² × rounds**, because everything is unrolled (no loops allowed). The rate is about 0.25 MB and 35 s per 10^4 QIR calls. A dt\* sweep at fixed total time needs many rounds at short periods; at d=5 and 100 rounds, extrapolation gives about 0.8 MB of bitcode and about 2 min per program, one program per (d, dt). Building is feasible; whether Nexus accepts or prices programs this large is unknown until the first cost run.
3. **Validator recursion:** worked around (section 3). It will grow with program size again if the stack or recursion settings are ever dropped.
4. **The physics of the experiment, not the tooling, is the open question for 7b.**
   - With `idle_dephasing` off, the program's RNG-injected Z stands in for idle noise during the period.
   - The emulator still adds its own transport dephasing, which the twin does not model. It can be switched off with `--transport-dephasing-off` for a like-for-like check.
   - Making the emulator itself idle for a chosen dt would need a delay or wait primitive in the QIR path. I have not checked whether one exists; that is a 7b item.
5. **Cost (HQC):** unknown until `stage7a_submit.py` is run with your Nexus login.

## Open items

- Decide between the pinned Guppy 0.x path and hugr-qir ≥ 0.2 on Guppy 1.x.
- Settle the 2π convention (it changes dt\* by 4×).
- Run the cost command.
- Use `checks="X"` for any emulator run meant to see the injected dephasing.
- Stage 6c: the full sweep finished (`results/stage6c/20261010T010451Z/`, untracked, on the 6c working tree) and its report is still to be written.
