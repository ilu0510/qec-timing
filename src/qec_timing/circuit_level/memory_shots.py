"""Raw Selene memory shots, stored once and decoded any number of times.

One block = one Selene run of ``surface_memory.Geometry`` at fixed (tau, rounds, noise, seed,
shot_offset). What is kept per shot is everything a decoder needs:

    det   detection events, [shots, (rounds + 1) * n_z], Geometry.detectors' definition and order
    obs   the logical observable (parity of the logical support in the final data readout)

Both are bit-packed (``np.packbits``, rows along axis 1 for ``det``) in
``results/circuit_level/cache/selene_memory_shots/<key>/arrays.npz``, with the usual
``meta.json`` (parameters, environment) and a ``result.json`` summary. The key does not depend
on which notebook asked, so 06 and 07 share blocks with identical parameters.

Same ``seed`` and ``shot_offset`` reproduce Selene's shots exactly, which is what makes a stored
block equivalent to the cached failure counts it was decoded into.
"""

from __future__ import annotations

import numpy as np

from . import runcache

EXPERIMENT = "selene_memory_shots"


def params(G, tau: int, rounds: int, shots: int, offset: int, seed: int,
           p_1q: float, p_meas: float, p_2q: float = 0.0) -> dict:
    return dict(code=G.code, distance=G.d, tau=tau, rounds=rounds, shots=shots,
                shot_offset=offset, seed=seed, p_1q=p_1q, p_2q=p_2q, p_meas=p_meas,
                circuit="surface_memory.Geometry: tau idle ops, then both check sectors")


def get(G, tau: int, rounds: int, shots: int, offset: int = 0, seed: int = 1,
        p_1q: float = 0.0, p_meas: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    """(det uint8 [shots, n_det], obs uint8 [shots]) -- from the store, or run Selene once."""
    prm = params(G, tau, rounds, shots, offset, seed, p_1q, p_meas)
    n_det = (rounds + 1) * len(G.zcols)
    arr = runcache.load_arrays(EXPERIMENT, prm)
    if arr is not None:
        det = np.unpackbits(arr["det"], axis=1, count=n_det)
        obs = np.unpackbits(arr["obs"], count=shots)
        assert det.shape == (shots, n_det) and obs.shape == (shots,)
        return det, obs
    syn, dat = G.run(tau, rounds, shots, offset=offset, seed=seed, p_1q=p_1q, p_meas=p_meas)
    det = G.detectors(syn, dat).astype(np.uint8)
    obs = G.observable(dat).astype(np.uint8)
    runcache.save_arrays(EXPERIMENT, prm, det=np.packbits(det, axis=1), obs=np.packbits(obs))
    runcache.save(EXPERIMENT, prm,
                  {"shots": int(shots), "n_detectors": int(n_det),
                   "mean_detector_rate": float(det.mean()),
                   "raw_observable_flips": int(obs.sum())},
                  extra={"format": "arrays.npz: det = packbits(axis=1) of [shots, n_det]; "
                                   "obs = packbits of [shots]"})
    return det, obs
