"""Circuit-level extension: the extraction circuit compiled and noised in Selene, and in Stim.

NOT the validated Stage 1-4 model. That model is phenomenological (ideal gates, explicit
bit-flip injection) and lives in ``qec_timing.noise`` and ``qec_timing.reference``. This package
runs the real syndrome-extraction circuit under Selene's ``DepolarizingErrorModel`` (charged per
native operation) or under circuit-level Stim noise, and measures what one round costs.

It builds its own extraction circuit (``surface_memory``), which overlaps teammate #1's
workstream; it is used only by the analysis in ``notebooks/``.

Copied from ``Quanti/noise_qec/source`` (commit 696c0f2) and repointed at this repo:
configs in ``config/``, the run cache in ``results/circuit_level/cache/`` (gitignored).

    runcache         content-addressed cache of simulation results
    surface_memory   surface-code memory in Selene: program, run, decode
    effective_noise  per-round charged-op counts (notebooks 04 and 04c, packaged)
    stim_memory      circuit-level Stim memory with emulator-style idle channels
"""
