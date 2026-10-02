# qec-timing — idle + readout noise model

The **error-model** workstream for our syndrome-timing project: a
phenomenological idle + readout noise model for a rotated surface-code memory,
implemented in Guppy and run on Quantinuum's Selene emulator (Stim backend).

Reproduces the setup of Haug, Bharti & Aolita (2026), *Exponential logical-error
reduction in quantum memories via optimal syndrome-measurement timing* — see
[`docs/paper/SOURCE.md`](docs/paper/SOURCE.md).

## Scope — please read this first

**I own the noise model only.** The syndrome-extraction circuit and the matching
decoder are other people's workstreams. Two directories here look like they
overlap with that work and do not:

- **`src/qec_timing/stub_circuit/`** — a placeholder extraction circuit, written
  to be provably equivalent to an ideal MPP so the noise model has something to
  run against. It is **test scaffolding, not a competing circuit**, and is
  marked as such in the code. Replace it with the real one when it lands.
- **`src/qec_timing/reference/`** — an independent pure-Python sampler (Stim +
  pymatching) used to cross-check the Guppy implementation, plus a
  `Decoder` protocol defining the seam the real decoder should implement. It is
  **a cross-check, not a decoder project**.

The paper's noise is **bit-flip (X) only**, so the memory is prepared in logical
|0⟩ and the decoded detector stream uses the Z-check sector. The new
**`src/qec_timing/circuit/`** implementation prepares full code states and
measures both stabilizer types on `q_d` independent patches. Start with
`build_surface_code_program(d=5, q_d=2)`; see
[the circuit guide](docs/surface_code.md) for usage and decoding. Optional X
memory tests the complementary Pauli sector separately. Historical validation
scripts still use the retained stub. There is **no gate noise** — gates and ancilla
measurement are ideal by construction, and all faults are injected explicitly.

## Setup

Requires **Python 3.12+** (`guppylang` will not install on 3.11; 3.14.0 is what
was used here).

```bash
python -m venv .venv
.venv/Scripts/python -m pip install -r requirements.txt   # Linux/macOS: .venv/bin/python
```

Run the test suite (some tests build and run real Selene programs):

```bash
.venv/Scripts/python -m pytest
```

## Where to start

1. **[`docs/summary.md`](docs/summary.md)** — one page: what the model is, how to
   call it, what is validated, what is not covered. Start here.
2. **[`docs/interfaces.md`](docs/interfaces.md)** — the per-shot **output
   schema** and the bit-packed detection-event wire format, with a copy-paste
   numpy unpacking snippet. This is the contract if you are consuming the
   output.
3. [`docs/noise_spec.md`](docs/noise_spec.md) — the full noise contract
   (Eqs. 1–3, one-Bernoulli-draw rule, record-flip readout, time bookkeeping).
4. [`docs/performance.md`](docs/performance.md) — why sweeps run on the
   reference path and validation runs on Guppy, with measured costs.

## Results

**Already computed — nobody needs to re-run anything.**
[`results/stage4a/`](results/stage4a/) holds the validation outputs: the
Guppy-vs-reference comparison (18/18 points agree), the `1/Δt*` vs `d` fit, the
plot, and the three follow-up checks (`T_test.json`, `seed_sensitivity.json`,
`rng_split.json`). Read alongside
[`docs/stage4a_report.md`](docs/stage4a_report.md).

Earlier stages' outputs are in `results/stage1/` and `results/stage2/`; those
are not tracked by git — regenerate with the scripts in `scripts/` if needed.

[`results/stage5/`](results/stage5/) holds full-circuit timing validation with
both X and Z stabilizer extraction on two independent patches, including a
low-statistics smoke test and a d=5,7 scan at fixed storage time T=10. The
scientific scan selected intervals Δt=1 for d=5 and Δt=0.5 for d=7 against a
Δt=2 baseline. Fresh-shot confirmation was **inconclusive on all four patches**
at the adjusted confidence level; timing improvement was not demonstrated.
All full-circuit/reference 99.9% intervals overlapped, a consistency diagnostic
rather than proof of equivalence. Gates remain ideal and noise is X-only.

With the project environment active, run:

```bash
python scripts/stage5_validation.py --quick --q-d 2
python scripts/stage5_validation.py --q-d 2
```

Each invocation saves counts, confidence intervals, confirmation summaries,
a plot, and a log in a new UTC-timestamped directory. See
[the Stage 5 guide](docs/surface_code.md#timing-premise-validation-stage-5)
for configuration and interpretation.

## Layout

```
src/qec_timing/
  layout.py        rotated surface-code geometry (shared)
  schedule.py      round scheduling, integer-tick clock (shared, no physics)
  ansatz.py        analytic Eqs. 5, 6, 7/D9, A3 (imports nothing else)
  noise/           Eqs. 1-3 + the Guppy noise primitives   <- the deliverable
  circuit/         full X/Z extraction on independent logical patches
  stub_circuit/    PLACEHOLDER extraction circuit + wire format
  reference/       independent Stim/pymatching cross-check + Decoder protocol
docs/              summary, interfaces, noise spec, performance, stage reports
scripts/           the runs that produced results/
tests/             unit, statistical, and Selene integration tests
```
