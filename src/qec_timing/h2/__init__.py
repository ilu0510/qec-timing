"""Stage 7: H2-Emulator (Nexus) path via Guppy 0.x -> HUGR -> QIR.

ISOLATED from the main code base:

* ``program.py`` needs Guppy 0.21.x (``.venv-h2``); the main ``.venv`` (Guppy
  1.x) must never import it.
* ``noise_params.py`` and ``geometry.py`` are standard-library only. They do
  not import anything from ``qec_timing`` either: the hook_safe order is
  re-implemented from ``circuit/ordering.py``'s logic, not imported.

This ``__init__`` deliberately imports nothing, so loading one submodule never
pulls in Guppy.
"""
