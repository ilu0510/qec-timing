# Interfaces: the `circuit_level` package

What the notebook arc's code exposes, and the conventions it relies on. The counterpart of
[`../stub_arc/interfaces.md`](../stub_arc/interfaces.md). Unlike that file, nothing here is a
contract with teammates #1 (circuits) or #3 (decoding): `circuit_level` builds its own extraction
circuit and decodes with PyMatching directly, for the notebooks' analysis only. It imports nothing
from `noise/`, `reference/`, `stub_circuit/`, `layout.py` or `schedule.py`.

## Package versions

The cached runs (`results/circuit_level/cache/*/meta.json`, ~7,700 of them) record:

| Package | Version |
| --- | --- |
| Python | 3.12.11 |
| guppylang | 1.0.3 |
| selene-sim | 0.3.0 |
| selene-hugr-qis-compiler | 0.4.2 |
| stim | 1.16.0 |
| PyMatching | 2.4.0 |
| numpy | 2.5.3 |
| PyYAML | (needed by `surface_memory`) |

These are **not** the `requirements.txt` pins (guppylang 1.1.1, selene-sim 0.3.2), and the
package has not been run under those. Each cache entry records the git commit and whether the
tree was dirty, so a result can be traced to the code that made it.

## Modules

| Module | What it is | Used in |
| --- | --- | --- |
| `surface_memory` | `Geometry`: the Selene surface-code memory — Guppy program generation, build, run, detectors, the hand-built matching graph | 03, 04c, 04d, 06a–06c, 07b; 09 through `native_memory` |
| `effective_noise` | per-round charged-op counts (`n_1`, `n_2`) measured in Selene; 04/04c packaged as `measure(code, d, with_p2q)` | `scripts/circuit_level_interval.py` |
| `memory_shots` | raw Selene shots (`det`, `obs`), bit-packed and stored once, decoded any number of times | 04d, 06a–06c |
| `selene_mirror` | Selene's native gate sequence rebuilt gate for gate in Stim; decoders on the full DEM | 04d, 06a–06c, 07b; 09 through `native_memory` |
| `stim_memory` | `idle_channel` and `StimMemory`: Stim's generated circuit with emulator-style idle noise on the data | 04, 04c, 07–09 |
| `native_memory` | `NativeMemory` (native gates + `idle_channel`) and `TimedNativeMemory` (time through the round, every qubit) | 07b, 09 |
| `timed_memory` | `TimedMemory` (textbook gates, time through the round, every qubit), the paper's model as a circuit, and the `p_ij` edge estimator | 08, 09 |
| `runcache` | content-addressed result cache | all |

## Memory classes

All five take the same rate arguments and are built per `(tau, rounds)`. `tau` is the full round
period in units of `t_R` (`>= 1` for the timed classes); see [`noise_spec.md`](noise_spec.md) §3.

```python
from qec_timing.circuit_level.timed_memory import TimedMemory
from qec_timing.circuit_level.native_memory import TimedNativeMemory

M = TimedNativeMemory("rotated", 5, basis="X", p1=9e-5, p2=5e-3, p_meas=4.7e-3,
                      p_init=1.9e-4, dephase_linear=1.5e-3, dephase_quad=0.0, rz_noise=False)
c = M.circuit(tau=4.0, rounds=48)                       # stim.Circuit
pL, sigma, fails = M.logical_error(4.0, 48, shots=100_000, seed=1)
```

| Class | Gates | Idle noise | Bases | Extra arguments |
| --- | --- | --- | --- | --- |
| `StimMemory` | textbook (Stim generated) | data, once per interval | Z, X | `idle_depol` |
| `TimedMemory` | textbook | every qubit, 7 layers + wait | Z, X (default X) | `ancilla_idle` (diagnostic) |
| `NativeMemory` | Selene native | data, once per interval | Z, X (default X) | `idle_depol`, `rz_noise` |
| `TimedNativeMemory` | Selene native, ASAP layers | every qubit, `L` layers + wait | Z, X (default X) | `rz_noise`; `.schedule(rounds)` |
| `surface_memory.Geometry` | Selene itself | `tau` explicit `phased_x(0,0)` ops (integer) | Z (hand decoder), X (mirror only) | `constants` |

`logical_error` returns `(p_L, sigma, failures)` for the timed and native classes and
`StimMemory`; decoding is always PyMatching on `circuit.detector_error_model(decompose_errors=True)`.

**Dephasing argument.** `dephase_quad` is `κ` in `sin²(κt/2)`. The notebooks pass
`κ = 2√c` to get `sin²(√c·t)`.

**Native gate sequences** are extracted from Selene once per `(code, d, basis, rounds)` and
cached (`native_memory_ops`). With multiprocessing, extract before forking (09 `run_all`): the
workers then only read the cache.

## The paper's model, and measuring it

```python
from qec_timing.circuit_level.timed_memory import (phenomenological_circuit,
    phenomenological_error, edge_probabilities, layout_edges, xor)

e = edge_probabilities(M.circuit(1.0, 20), shots=400_000, seed=1,
                       edges=layout_edges(M.G.sup, 20))      # native circuits: no coordinates
e["p_data"], e["p_read"], e["p_diag"], e["p_space"], e["p_data_se"], e["p_read_se"]

pL, sigma, fails = phenomenological_error("rotated", 5, rounds=48, p_data=e["p_data"],
                                          p_read=e["p_read"], shots=100_000)
```

- `edge_probabilities` defaults to `check_edges(c)`, which reads detector coordinates from a Stim
  **generated** circuit. Circuits from `selene_mirror` (`NativeMemory`, `TimedNativeMemory`) carry
  none: pass `edges=layout_edges(M.G.sup, rounds)`.
- Only bulk rounds are used (not the first round, not the final readout layer). Moments are
  accumulated in chunks of 50,000 shots, so memory does not grow with `shots`.
- `xor(*ps)` is the probability of an odd number of independent flips, used to compose
  `p_data(1) ⊕ q_wait`.

## Detector and observable conventions

There are **two** detector orderings, and neither is the stub arc's
(`RotatedSurfaceCodeZSector`, `docs/stub_arc/interfaces.md`):

| Circuits | Ordering | Coordinates |
| --- | --- | --- |
| `StimMemory`, `TimedMemory`, `phenomenological_circuit` | Stim's generated `surface_code:*_memory_*` circuits | `(x, y, round)` |
| `Geometry`, `selene_mirror`, `NativeMemory`, `TimedNativeMemory` | round-major over the decoded checks in `Geometry.cols` order: detector `j * n_checks + i` compares check `i` in rounds `j` and `j − 1`; final layer from the data readout | none |

Common to both: `D[j, c] = s[j, c] ⊕ s[j−1, c]` over the checks of the memory basis (Z-checks for
Z, X-checks for X), round 0 against the deterministic initial value, one observable (parity of the
logical support in the final readout).

`Geometry`'s Selene program allocates a fresh ancilla array **every round** and measures
**both** check sectors, decoding only the memory basis's. The X-basis program, CNOT order and
logical support come from Stim's `*_memory_x` circuit.

**Shots from Selene:** `Geometry.run` returns `(syndromes [shots, rounds, n_anc],
data [shots, n_data])`; `memory_shots.get` returns `(det [shots, n_det], obs [shots])` as `uint8`,
from the store or by running Selene once.

## The hand-built decoder (retired from use)

`Geometry.priors` + `Geometry.matching` build a phenomenological graph (data edges in a round,
measurement edges between rounds) from the measured op counts. Z basis only. It misses mid-round
and propagated faults and misstates `p_L` by up to 7–8× (06b); it is kept for that comparison
(`selene_mirror.decode_hand`). Use `selene_mirror.decode_mirror` / `matcher(circuit)`.

## `selene_mirror.to_stim`

```python
to_stim(ops, layout, p_1q=0.0, p_2q=0.0, p_meas=0.0, allow_unverified_p2q=False, *,
        idle=None, p_init=0.0, rz_noise=True) -> stim.Circuit
```

- Defaults reproduce Selene's `DepolarizingErrorModel`: a uniform Pauli **before** each native
  single-qubit gate, `p_meas` as a record flip.
- `p_2q > 0` is refused unless `allow_unverified_p2q=True`: its placement (`DEPOLARIZE2` before
  each `Rzz`) is by analogy and has not been read off Selene's log.
- Each native gate is converted through its unitary, and the round trip is checked: Stim's
  `Tableau.from_unitary_matrix` silently returns a wrong Clifford for a non-Clifford gate, so a
  mismatch raises.
- `ops` may include `{"op": "Idle", "qubits": [...], "channel": (pX, pY, pZ)}` entries, which is
  how `TimedNativeMemory` places idle noise inside the round.

## Run cache

```python
from qec_timing.circuit_level import runcache

hit = runcache.load(EXPERIMENT, params)          # dict, or None on a miss
if hit is None:
    hit = compute(params)
    runcache.save(EXPERIMENT, params, hit)
runcache.index("09_native_uniform_noise")        # every cached run, newest first
runcache.summary()
```

- Layout: `results/circuit_level/cache/<experiment>/<key>/{meta.json, result.json[, arrays.npz]}`,
  gitignored. `key` is a 16-hex SHA-1 of the parameters, with floats formatted `%.12g`, so `0.002`
  and `2e-3` share a key.
- `meta.json` holds the parameters, the save time and the environment (git commit, dirty flag,
  Python, host, package versions). `result.json` is written last, so an interrupted save reads as
  a miss.
- Array results (`save_arrays` / `load_arrays`) go in `arrays.npz`.
- The notebooks key each job as `(kind, sorted params)` with `kind` in `full`, `edges`, `pheno`,
  and later notebooks read earlier ones' entries directly (09 §5 reads 07, 07b and 08).

## Paths the package assumes

- `surface_memory.CONFIG` and `runcache.DATA` are both `parents[3]` of the module file, i.e. the
  repo's `config/` and `results/circuit_level/cache/`. **Moving `circuit_level/` to a different
  depth breaks both.** Assign `runcache.DATA` to keep runs elsewhere.
- Notebooks import the package with `sys.path.insert(0, str(pathlib.Path.cwd().parent / 'src'))`,
  so they must be run from `notebooks/`.
