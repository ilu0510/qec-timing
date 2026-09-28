# Project: optimal syndrome-measurement timing on Quantinuum Selene

Reproducing results from Haug, Bharti & Aolita (2026), "Exponential logical-error reduction in
quantum memories via optimal syndrome-measurement timing" (docs/paper/). Clean equations are in
docs/paper_equations.md; prefer that file over the PDF text for any formula.

## Team split (I own #2)
1. Circuits: rotated surface code syndrome extraction (teammate)
2. Error models: idle + readout noise (me, this repo's main focus)
3. Decoding: matching decoder (teammate)
Interfaces between these live in docs/interfaces.md. Do not change them without asking me.

## Two paths: validation vs production (accepted design, Stage 4a)
The Guppy/Selene path is the VALIDATION path. The pure-Python reference is the PRODUCTION path for
parameter sweeps. This is a deliberate division of labour, not a shortfall.
- Selene/Guppy simulates the noise model as it would run on hardware, and is the artefact the physics
  is defined by. It costs ~200x more per syndrome round than the reference (measured; see
  docs/performance.md), because every Bernoulli draw and gate goes through the emulator runtime.
- The reference (Stim + pymatching) reproduces the same detector stream ~200x faster, and is verified
  against the Guppy path shot-for-shot: identical detector ordering, agreeing logical error rates.
- Therefore: sweeps that need many shots per point run on the reference. Guppy validates the model at
  d = 5, 7, 9 with T = 10, where the two can be compared directly within error bars in minutes.
- What makes this sound is the cross-check, so it must not be skipped or weakened. If the reference is
  ever changed, re-run the Stage 4a cross-check before trusting a sweep.
- Do NOT optimise the Guppy path to close the gap. Geometric-gap sampling for the data noise was
  considered and rejected in Stage 4a: it buys ~1.05x (measured) against a ~10,000x gap and puts a validated noise
  model at risk. Measured RNG-vs-gate split in docs/performance.md.
- Guppy pays ~9.3 s PER run() CALL, independent of shot count (1 shot and 8 shots cost the same).
  So: size shots up front and issue ONE run() call per configuration. Never loop run() to accumulate
  statistics -- each extra batch costs another ~9 s. Adaptive sampling belongs on the reference path.
  Guppy per-round cost is unstable to ~3x between runs, so cap work by a time budget, not a shot count.

## Stage-1 physics (phenomenological model, NO gate noise)
- Gates and stabilizer extraction are ideal. All noise is injected explicitly in the Guppy program.
- Noise is bit-flip (X) only. Memory is prepared in logical |0>; only the Z-check sector matters.
- Per round of duration dt, at dimensionless time t:
  - p_idle = 1 - exp(-p * lam(t) * dt)
  - p_data = 1 - (1 - p) * (1 - p_idle)  -> apply X to each data qubit w.p. p_data,
    as ONE Bernoulli draw (not two independent flips), BEFORE the extraction circuit
  - p_read = b_read * p -> XOR each recorded syndrome bit w.p. p_read
    (flip only the classical record, not the qubit)
- The final data-qubit readout also has classical flips w.p. p_read.
- Round count: never loop on a float accumulator (float drift adds rounds). Fixed dt: n_rounds = round(T/dt), loop over an integer counter, t_j = j*dt. Variable dt (adaptive stages): all intervals are integer multiples of a tick (default 1e-3); keep time as an integer tick counter. Total time T is fixed; the number of rounds varies with dt.
- p_stab (the p inside p_data) must NEVER be dropped: it is what creates the dt trade-off.
- Figure of merit: logical error per unit time at fixed T, not per round.
- Each run must also output, per round: dt, t, p_data and p_read (the decoder needs them for weights).

## Randomness inside Guppy
- Use the platform RNG: guppylang.std.qsystem.random.RNG. Bernoulli(p) = rng.random_float() < p.
- Use exactly ONE RNG instance per program, passed through all functions (instances are not independent).
- Seed per shot: seed = base_seed + get_current_shot() (guppylang.std.qsystem.utils). A fixed
  seed gives identical noise in every shot and invalidates all statistics.
- The RNG is a linear resource: call rng.discard() at the end.
- Fallback if the platform RNG is unsupported: guppylang.std.random.seeded_pcg32 (integers only;
  Bernoulli via next_int_bounded compared against round(p * bound)).

## Selene/Guppy facts (verified in Stage 0, see docs/interfaces.md)
- Run: prog.emulator(n_qubits=N).stabilizer_sim().with_shots(K).with_seed(S).run(**args)
- Emit results with output(tag, value); result() is a deprecated alias. Read with collate_tags(), never as_dict()
  (as_dict drops repeated tags).
- Reduce in-program: emit detection events (syndrome XOR previous round), BIT-PACKED into 64-bit words
  (wire format in docs/interfaces.md); never dump raw per-qubit-per-round data at large d. Sparse indices
  lose: the observed detection rate runs 5-21x above the sparse/packed crossover at both d=11 and d=27,
  and sparse costs 5.9x packed at d=27, dt=10 (Stage 2 benchmark).
- Emission is NOT the bottleneck once packed: ~3-4% of run time, measured against a no-emission control.
  RNG draws and gate operations dominate. (The earlier ~40 us/bit figure was real but described dumping
  1500 loose bits per shot, not this workload.) If more throughput is needed, target the RNG call rate.
- Noise randomness depends only on base_seed (with_seed does not affect RNG draws). Use distinct base_seeds per
  batch, or with_shot_offset() to extend a run.
- Only ONE RNG instance may exist: a second RNG(seed) reseeds the global runtime PRNG.
- Get package versions with importlib.metadata.version, not __version__.


## Engineering rules
- Guppy/Selene APIs change quickly. Before using any guppylang or selene-sim function, check it
  exists in the INSTALLED package source (pip show / inspect the site-packages source). Do not rely
  on memory. Record the versions used in docs/interfaces.md.
- Use Selene's stabilizer (Stim) simulator backend. Everything is Clifford; d up to 27 is needed.
- Every stochastic component gets a statistical test: empirical rate within a binomial confidence
  interval of the target, with fixed seeds.
- src/qec_timing/reference/ is an independent pure-Python (numpy + pymatching) phenomenological
  sampler used ONLY to cross-check the Selene implementation. It is not the product.
- src/qec_timing/stub_circuit/ is a placeholder for the teammate's circuit and must be marked as such.
- Work in small steps. Stop and report at the end of each stage; don't start the next one unasked.

## Roadmap (scope: Stages 0-4; 5-6 deferred to a team decision)
0. Environment + Selene capability checks [done]
1. Pure-Python reference model; Appendix C thresholds [done]
2. Noise injection in Guppy (static lam) + statistical tests [done]
3. Time bookkeeping: dynamic rounds, per-round metadata output
4. dt sweep reproducing Fig. 2a; Selene vs reference cross-check. FINAL STAGE.
Deferred, not in scope: burst lam(t) (Eq. 9) and the adaptive controller (Appendix F). The scheduling
and noise code should not block these later, but do not build them.