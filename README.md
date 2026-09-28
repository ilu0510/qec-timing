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

The noise is **bit-flip (X) only**, so the memory is prepared in logical |0⟩ and
**only the Z-check sector is exercised**. The X sector is equivalent under
X ↔ Z and is not simulated. There is **no gate noise** — gates and ancilla
measurement are ideal by construction, and all faults are injected explicitly.

## Setup

Requires **Python 3.12+** (`guppylang` will not install on 3.11; 3.14.0 is what
was used here).

```bash
python -m venv .venv
.venv/Scripts/python -m pip install -r requirements.txt   # Linux/macOS: .venv/bin/python
```

Verify — **110 tests**, a few minutes (they build and run real Selene programs):

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

## Layout

```
src/qec_timing/
  layout.py        rotated surface-code geometry (shared)
  schedule.py      round scheduling, integer-tick clock (shared, no physics)
  ansatz.py        analytic Eqs. 5, 6, 7/D9, A3 (imports nothing else)
  noise/           Eqs. 1-3 + the Guppy noise primitives   <- the deliverable
  stub_circuit/    PLACEHOLDER extraction circuit + wire format
  reference/       independent Stim/pymatching cross-check + Decoder protocol
docs/              summary, interfaces, noise spec, performance, stage reports
scripts/           the runs that produced results/
tests/             110 tests
```
