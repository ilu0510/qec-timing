# Noise specification

The phenomenological noise contract for this project, in one place. Written for
teammates #1 (circuits) and #3 (decoding); it is meant to stand alone, so it
repeats what it needs from the paper rather than pointing at it.

Equation numbers refer to `docs/paper_equations.md`. Implemented twice, on
purpose: `src/qec_timing/noise/params.py` (product path, used by the Guppy
program) and `src/qec_timing/reference/noise.py` (independent pure-Python
reference). A test asserts the two agree; keeping them separate is what makes
the reference an independent check rather than a mirror.

## 1. Scope

- **Phenomenological model, no gate noise.** Gates and stabilizer extraction are
  ideal. Every fault is injected explicitly.
- **Bit-flip (X) noise only.** The memory is prepared in logical |0>, so only
  the **Z-check sector** matters. Z-checks detect X-type data faults. The X
  sector is equivalent under X <-> Z and is not simulated.
- Ancilla measurement in the extraction circuit is **ideal**; readout noise
  enters only as a classical flip of the recorded bit (section 4).

## 2. Noise parameters (Eqs. 1-3)

For a round of duration `dt` at idling rate `lam`, with base physical error
rate `p` and read-out coefficient `b_read`:

```
p_idle = 1 - exp(-p * lam * dt)                       (Eq. 1)
p_stab = p                                            (Eq. 2)
p_data = 1 - (1 - p_stab) * (1 - p_idle)              (Eq. 2)
p_read = b_read * p                                   (Eq. 3)
```

Non-negotiable points:

- **`p_stab` must never be dropped.** At `lam = 0`, `p_idle = 0` and
  `p_data = p` (algebraically exact; in floating point it differs by ~1e-17, so
  compare with a tolerance rather than `==`). `p_stab` is the measurement-induced term that makes
  frequent measurement costly; without it there is no trade-off in `dt` and the
  whole result disappears.
- **Never substitute the linearised `p * lam * dt` for Eq. 1.** The linearisation
  is only valid for small `p * lam * dt`, and the `dt` sweeps run to
  `lam * dt = 10` and beyond, where it is wrong by tens of percent.
- `p_read` is independent of `dt`: read-out noise does not accumulate with
  waiting time.

### Where `p_data` is computed

**`exp()` is not available inside Guppy.** There is no `guppylang.std.maths`
module and no `exp` in `guppylang.std.builtins` (checked against the installed
guppylang 1.1.1). Therefore:

> **Stage 2 route: `p_data` and `p_read` are computed in Python from Eqs. 1-3
> and passed into the Guppy program as runtime `float` arguments.**

This costs nothing: Stage 0 established that `float` entrypoint arguments are
bound per run without recompiling, so a whole `dt` sweep reuses one compiled
binary.

For later stages, where a controller may choose `dt` at runtime and `p_data`
must then be computed *inside* the program: `float ** float` **is** available in
Guppy, and `math.e ** (-z)` was verified to reproduce `math.exp(-z)`
bit-for-bit (`0.6065306597126334` at `z=0.5`, `0.1353352832366127` at `z=2.0`).
So Eq. 1 can be evaluated in-program as `1.0 - E ** (-p * lam * dt)` with `E`
captured as a comptime constant. Use that rather than a linearisation.

## 3. Data faults: one Bernoulli draw per qubit per round

Before the extraction circuit, each data qubit independently suffers X with
probability `p_data`:

```python
for i in range(n_data):
    if rng.random_float() < p_data:      # ONE draw
        x(data[i])
```

**One draw, not two.** `p_data` already combines the idling and
measurement-induced mechanisms through Eq. 2. Drawing separately for `p_stab`
and `p_idle` and applying X twice is a different channel: two flips cancel, so
the effective rate becomes `p_stab + p_idle - 2 p_stab p_idle` instead of the
specified `p_stab + p_idle - p_stab p_idle`. The difference is second order but
systematic, and it grows exactly where the `dt` trade-off lives.

Faults are applied **before** extraction, so a fault in round `j` is visible in
round `j`'s syndrome.

## 4. Read-out faults: flip the record, not the state

```python
bit = measure(ancilla).read()            # ideal measurement
if rng.random_float() < p_read:
    bit = not bit                        # classical flip only
```

The qubit is **not** touched. A read-out fault must not propagate into the next
round's syndrome; it produces a pair of detection events separated in time (one
at round `j`, one at `j+1`), which is what makes it distinguishable from a data
fault by the decoder.

The same rule applies to the **final data readout**: each measured data bit is
flipped with probability `p_read` before being used in the final detectors and
the logical observable.

## 5. Time bookkeeping

**Never loop on a float accumulator.** Accumulating `t += dt` until `t < T`
silently changes the round count: with `dt = 0.1, T = 1.0` the accumulated sum
after ten steps is `0.9999999999999999 < 1.0`, giving **11 rounds instead of
10**.

- **Fixed `dt`:** `n_rounds = round(T / dt)`, loop over an integer counter,
  `t_j = j * dt`.
- **Variable `dt` (adaptive stages):** every interval is an integer multiple of
  a tick (default `1e-3`); carry time as an integer tick counter and convert
  only for output.
- `T` is fixed; the number of rounds varies with `dt`. Because `n_rounds` is
  rounded, the realised total time is `n_rounds * dt`, which may differ slightly
  from nominal `T`. **Use the realised value** whenever a rate per unit time is
  computed.

## 6. Detector and observable conventions

Shared by the reference model and the Guppy path, so the two can be compared
shot for shot.

| Quantity | Value |
| --- | --- |
| Data qubits | `d**2` |
| Z-checks per round | `(d**2 - 1) // 2` |
| Detector rounds emitted | `n_rounds + 1` |
| Observables | 1 (logical Z) |

- **Detection events** follow Eq. A1: `D_{j,c} = s_{j,c} XOR s_{j-1,c}`.
- **Round 0** compares against the deterministic initial value: the Z-stabilizers
  are all +1 on logical |0>, so `D_{0,c} = s_{0,c}`.
- **The final block** (`j = n_rounds`) comes from the data readout: for each
  check, the parity of the final measured bits over its support, XORed with that
  check's last syndrome.
- **Logical Z observable** is the parity of the final data readout over **row 0**
  of the `d x d` grid. In this layout the logical X operator is a *column* of X
  (a column commutes with every Z-check; a row does not), and a column meets
  row 0 in exactly one qubit, so a logical X fault flips this parity.
- Geometry is defined once, in `src/qec_timing/layout.py`
  (`RotatedSurfaceCodeZSector`), and imported by **both** the reference model
  and the Guppy path so the two cannot drift. Geometry is shared deliberately;
  the noise equations in section 2 are deliberately **not** (see the header).

The wire format for detection events (bit-packed words, bit/word order,
padding, and a numpy unpacking snippet) is specified in `docs/interfaces.md`.

## 7. Randomness

- Use the platform RNG, `guppylang.std.qsystem.random.RNG`. `Bernoulli(p)` is
  `rng.random_float() < p`.
- **Exactly one RNG instance per program**, threaded through every helper.
  Instances are not independent: `RNG(seed)` maps onto Selene's *global* runtime
  PRNG, so constructing a second one reseeds the first. Helpers take the RNG as
  a borrowed parameter (`rng: RNG`), which does not consume it.
- **Seed per shot:** `base_seed + get_current_shot()`. A seed that does not vary
  with the shot gives byte-identical noise in every shot and silently invalidates
  every statistic computed from the run.
- The RNG is a linear resource: call `rng.discard()` exactly once, at the end.
  This is enforced at compile time, not at runtime -- omitting it is a
  "Drop violation" error, not a leak.
- Draws are **independent of Selene's `.with_seed()`** when the RNG is seeded
  explicitly. Reproducibility rests on `base_seed` alone.

## 8. What is constant in Stage 2

Stage 2 is **static `lam`** and fixed `dt`, so `p_data` and `p_read` are the same
in every round. They are emitted **once per shot** rather than per round,
together with `n_rounds`, `n_checks` and `n_words`.

CLAUDE.md asks for per-round `dt, t, p_data, p_read`. With a static schedule
those four values are constant, so emitting them once is information-equivalent
and the round index is implicit in the block ordering of the detection-event
stream. **Stage 3 introduces per-round metadata** once `dt` varies, at which
point the per-round form becomes necessary and this section should be revised.

For the record, the cost argument for deferring it turned out to be weak: the
Stage 2 benchmark found that emission accounts for only **3-4% of run time**
with bit-packing, not the dominant share a Stage 0 extrapolation suggested
(measured against a no-emission control; see `docs/interfaces.md`). Per-round
metadata would therefore have been affordable. It is deferred because it is
redundant under a static schedule, not because it is expensive.
