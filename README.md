# qec-timing — idle + readout noise models for syndrome timing

The **error-model** workstream for our syndrome-timing project. Reproduces the setup of
Haug, Bharti & Aolita (2026), *Exponential logical-error reduction in quantum memories via
optimal syndrome-measurement timing* — see [`docs/paper/SOURCE.md`](docs/paper/SOURCE.md).
Equations are in [`docs/paper_equations.md`](docs/paper_equations.md), symbols in
[`docs/notebook_arc/notation.md`](docs/notebook_arc/notation.md).

The repo holds **two arcs** that ask opposite questions about the paper's model:

| | **Stub arc** (phenomenological, Stages 0–4a) | **Notebook arc** (circuit-level, `notebooks/03`–`09`) |
| --- | --- | --- |
| Question | *If* the noise is exactly the paper's model, does the implementation reproduce the paper? | *Does* a realistic circuit's own noise reduce to the paper's model? |
| Direction | model → circuit: the model is imposed | circuit → model: the model is tested |
| Noise | X flips injected by the program; gates ideal | gate noise + trapped-ion idle dephasing; nothing injected |
| Δt | a number that sets `p_idle` | real waiting time, in units of the round duration `t_R` |
| `p`, `λ`, `b_read` | inputs (Eqs. 1–3) | outputs, measured from detection events |
| Memory | Z basis, Z-checks | X basis (Z idle noise is invisible to a Z-basis memory) |
| Validated by | two independent implementations agreeing (Guppy vs Stim reference) | the paper's model, fed measured parameters, predicting the circuit's Δt\* |
| Code | `noise/`, `stub_circuit/`, `reference/` | `circuit_level/` |
| Status | **validated** | premise holds in simulation; device parameters unmeasured |

The two share only `ansatz.py` (Eqs. 5–7, analytic). Neither arc's code imports the other's.

## Scope — please read this first

**I own the noise model only.** The syndrome-extraction circuit and the matching decoder are
other people's workstreams. Three directories here look like they overlap with that work and
do not:

- **`src/qec_timing/stub_circuit/`** — a placeholder extraction circuit, written to be provably
  equivalent to an ideal MPP so the noise model has something to run against. It is **test
  scaffolding, not a competing circuit**, and is marked as such in the code. Replace it with
  the real one when it lands.
- **`src/qec_timing/reference/`** — an independent pure-Python sampler (Stim + pymatching) used
  to cross-check the Guppy implementation, plus a `Decoder` protocol defining the seam the real
  decoder should implement. It is **a cross-check, not a decoder project**.
- **`src/qec_timing/circuit_level/`** — builds its own surface-code extraction circuit and
  decodes with pymatching directly, **for the notebook arc's analysis only**. It does not use
  `layout.py` or the `Decoder` protocol, and it is not the validated Stage 1–4 model.

## Setup

Requires **Python 3.12+** (`guppylang` will not install on 3.11; 3.14.0 is what the stub arc
used).

```bash
python -m venv .venv
.venv/Scripts/python -m pip install -r requirements.txt   # Linux/macOS: .venv/bin/python
```

Verify — **110 tests** for the stub arc, a few minutes (they build and run real Selene
programs), plus `tests/test_circuit_level.py` for the notebook arc:

```bash
.venv/Scripts/python -m pytest
```

**Environment caveat.** `requirements.txt` pins guppylang 1.1.1 / selene-sim 0.3.2. The
notebook arc's cached runs were recorded under **guppylang 1.0.3 / selene-sim 0.3.0**
(Python 3.12), and the stub does not build there. Whether `circuit_level/` runs unchanged
under the pinned versions has not been checked.

---

## Stub arc: the phenomenological model (Stages 0–4a)

A phenomenological idle + readout noise model for a rotated surface-code memory, implemented
in Guppy and run on Quantinuum's Selene emulator (Stim backend).

The noise is **bit-flip (X) only**, so the memory is prepared in logical |0⟩ and **only the
Z-check sector is exercised**. The X sector is equivalent under X ↔ Z and is not simulated.
There is **no gate noise** — gates and ancilla measurement are ideal by construction, and all
faults are injected explicitly.

### Where to start

1. **[`docs/stub_arc/summary.md`](docs/stub_arc/summary.md)** — one page: what the model is, how to call it, what
   is validated, what is not covered. Start here.
2. **[`docs/stub_arc/interfaces.md`](docs/stub_arc/interfaces.md)** — the per-shot **output schema** and the
   bit-packed detection-event wire format, with a copy-paste numpy unpacking snippet. This is
   the contract if you are consuming the output.
3. [`docs/stub_arc/noise_spec.md`](docs/stub_arc/noise_spec.md) — the full noise contract (Eqs. 1–3,
   one-Bernoulli-draw rule, record-flip readout, time bookkeeping).
4. [`docs/stub_arc/performance.md`](docs/stub_arc/performance.md) — why sweeps run on the reference path and
   validation runs on Guppy, with measured costs.
5. Stage reports: [`stage1`](docs/stub_arc/stage1_report.md), [`stage2`](docs/stub_arc/stage2_report.md),
   [`stage3`](docs/stub_arc/stage3_report.md), [`stage4a`](docs/stub_arc/stage4a_report.md).

### Results

**Already computed — nobody needs to re-run anything.** [`results/stage4a/`](results/stage4a/)
holds the validation outputs: the Guppy-vs-reference comparison (18/18 points agree), the
`1/Δt*` vs `d` fit, the plot, and the three follow-up checks (`T_test.json`,
`seed_sensitivity.json`, `rng_split.json`). Read alongside
[`docs/stub_arc/stage4a_report.md`](docs/stub_arc/stage4a_report.md).

Earlier stages' outputs are in `results/stage1/` and `results/stage2/`; those are not tracked
by git — regenerate with the `stage1_*` / `stage2_*` scripts in `scripts/` if needed.

---

## Notebook arc: does the premise survive a real circuit?

The paper's model has two knobs per round, a data flip `p_data = p_stab ⊕ p_idle(Δt)` and a
record flip `p_read = b_read·p`, with `λ` and `b_read` assumed. The notebooks start from a
compiled surface-code circuit and add realism one ingredient at a time, asking at each step
whether an interior optimum Δt\* survives and whether the paper's model still describes the
circuit.

### The narrative

| Notebook | Question | Finding |
| --- | --- | --- |
| **03** noise channel audit | Can Selene express the paper's channels? | `p_meas` is the paper's `p_read` (a record flip). Selene has no pure bit-flip and no duration term |
| **04 / 04c** extraction cost | What does one extraction round cost the data? | Idling written as explicit ops gives an effective `λ` ≈ 0.10–0.15 and `b_read` ≈ 1.35–1.50, set by circuit structure (04) or the gate-error ratio `p_2q/p_1q` (04c) |
| **04d** readout composition | Does readout error add the way 05–07 assume? | Yes, as two independent flips `a ⊕ b`; two bookkeeping errors in 05/06c, neither changes a conclusion |
| **05** paper vs effective parameters | What do the effective parameters do to Δt\*? | Lengthen it 2.6–3.0×, almost all through `λ`, as Eq. 6 says |
| **06a** Stim mirror | Can Selene's compiled circuit be copied gate for gate into Stim? | Yes, statistically indistinguishable; gives Selene's full detector error model |
| **06b** decoder comparison | Does the decoder matter? | Up to 7–8× in `p_L`, and it decides whether an optimum appears. 06c onward use the full error model |
| **06c** Selene interval | Is there an optimum under Selene's own noise? | Yes: Δt\* = 12.0, 5.1, 3.4 at d = 3, 5, 7, shortening with d as predicted |
| **07** dephasing | Idle noise of trapped-ion form (Z only, linear `γ` + quadratic `c`) | Needs an X-basis memory. Interior optimum on all 21 curves; at H2's estimated idle-to-gate ratio, back to back is best except a small gain at d = 3 |
| **07b** native gates | Does Selene's native CNOT decomposition matter? | Δt\* 0–7% longer; every conclusion of 07 stands |
| **08** uniform noise | Same noise on every qubit, ancillas included; `b_read` measured, not assumed | **Premise holds**: (P1) interior optimum; (P2) `p_read` flat in Δt, `p_data` = per-round part ⊕ wait; (P3) Δt\* predicted from one calibration |
| **09** full model | Native gates and uniform noise together | Same as 08. Δt\* from one calibration: +17–21% (d = 3), +4–9% (d = 5), −2…+1% (d = 7); `p_L` underestimated more with d. The quadratic H2 case depends on an unresolved convention |

Groundwork, not part of the numbered arc: `noise_analysis` (does Selene separate idling from
measurement?), `qec_timestep` (one Steane round, by hand) and `capability_audit` (what Selene
can be asked at all). Superseded notebooks are kept in `notebooks/retired/`.

### Where to start

1. **[`docs/notebook_arc/summary.md`](docs/notebook_arc/summary.md)** — one page: the ladder,
   the premise (P1–P3), what it means for H2, what is not covered. Start here.
2. [`docs/notebook_arc/noise_spec.md`](docs/notebook_arc/noise_spec.md) — the noise models:
   Selene's, the idle channel, where it is applied in each notebook, the rates, and how `p_data`
   and `p_read` are measured from detection events.
3. [`docs/notebook_arc/interfaces.md`](docs/notebook_arc/interfaces.md) — the `circuit_level`
   classes, detector conventions and run cache.
4. [`docs/notebook_arc/performance.md`](docs/notebook_arc/performance.md) — what a run costs.
5. The notebooks are the per-step reports. Read
   **[`09_native_uniform_noise.ipynb`](notebooks/09_native_uniform_noise.ipynb)** first (the full
   model and the 07 → 07b → 08 → 09 ladder), then
   [`08`](notebooks/08_uniform_noise_premise.ipynb) (the premise and the edge estimator) and
   [`07`](notebooks/07_dephasing_interval.ipynb) (why X basis; the H2-ratio conclusion).
6. [`docs/notebook_arc/notation.md`](docs/notebook_arc/notation.md) — symbols and time units
   across 04–09.

### Results

Each notebook is saved with its outputs, and its figures are in
[`notebooks/figures/`](notebooks/figures/). Simulation results are cached by
`circuit_level/runcache` under `results/circuit_level/cache/`, which is **gitignored**: on a
fresh clone, re-running a notebook recomputes everything: ~15–20 min each for 07b–09 on 16
workers, 1–2 h per distance for the Selene sweeps of 06b/06c (see
[`performance.md`](docs/notebook_arc/performance.md)). Run them in order: later notebooks read
earlier ones' cache.

---

## Hardware plans (nothing run on hardware yet)

| Plan | Arc | First blocker |
| --- | --- | --- |
| [`docs/stub_arc/stub_to_device.md`](docs/stub_arc/stub_to_device.md) | stub | Does the platform RNG / `get_current_shot()` work on the device? |
| [`docs/notebook_arc/hardware_tests.md`](docs/notebook_arc/hardware_tests.md) | notebook | There is no user-controlled delay instruction |

Slide versions of both plans and of the notebook arc are in [`ppts/`](ppts/).

## Layout

```
src/qec_timing/
  layout.py        rotated surface-code geometry (stub arc, shared with reference)
  schedule.py      round scheduling, integer-tick clock (stub arc, no physics)
  ansatz.py        analytic Eqs. 5, 6, 7/D9, A3 (imports nothing; used by both arcs)
  noise/           Eqs. 1-3 + the Guppy noise primitives   <- the stub-arc deliverable
  stub_circuit/    PLACEHOLDER extraction circuit + wire format
  reference/       independent Stim/pymatching cross-check + Decoder protocol
  circuit_level/   notebook arc: Selene circuit, Stim mirror, timed/native memories, run cache
config/            notebook arc: measured extraction constants (yaml)
notebooks/         notebook arc 03-09, groundwork, figures/, retired/
docs/
  paper_equations.md, paper/   shared: the paper's equations and source
  stub_arc/        summary, interfaces, noise spec, performance, stage reports, stub_to_device
  notebook_arc/    summary, interfaces, noise spec, performance, notation, hardware_tests
scripts/           stage1_*/stage2_*/stage4*: stub arc;  circuit_level_*: notebook arc
tests/             test_stage0-3: stub arc (110);  test_circuit_level: notebook arc
results/           stage4a/ (committed);  circuit_level/ (gitignored cache)
ppts/              slide decks
```
