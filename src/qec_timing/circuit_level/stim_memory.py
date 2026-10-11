"""Circuit-level surface-code memory in Stim, with emulator-style noise channels.

The Selene arm is limited by what Selene's `DepolarizingErrorModel` can express (see
notebooks/03_noise_channel_audit.ipynb). This arm relaxes those limits where Stim can:

    channel                    Selene arm                      this arm
    idling                     tau explicit ops at p_1q        a Pauli channel over the idle
                               (symmetric, op-count)           DURATION: depolarizing + linear
                                                               and quadratic dephasing
    memory basis               |0_L>, Z-checks only            Z or X (X is needed to see Z-type
                                                               dephasing at all)
    decoder weights            per-qubit marginals, no         Stim's exact detector error model:
                               correlations                    hook errors, correlated faults
    reset error                off                             p_init

It keeps Stim's own notion of a gate: a CX is one native two-qubit gate charged `p2`, and the
H on each X-ancilla is charged `p1`. There is no native-gate decomposition, crosstalk,
leakage, coherent over-rotation or asymmetric readout -- those need PECOS or the hosted
emulator.

    M = StimMemory("rotated", 5, basis="X", p1=1e-3, p2=5e-3, p_meas=2e-3,
                   idle_depol=1e-3, dephase_linear=2e-3)
    pL, sigma = M.logical_error(tau=8, rounds=8, shots=100_000)
"""

from __future__ import annotations

import numpy as np
import pymatching
import stim

SPEC = {("rotated", "Z"): "surface_code:rotated_memory_z",
        ("rotated", "X"): "surface_code:rotated_memory_x",
        ("unrotated", "Z"): "surface_code:unrotated_memory_z",
        ("unrotated", "X"): "surface_code:unrotated_memory_x"}

# Placeholder probabilities that tag each noise location in Stim's generated circuit, so it
# can be rewritten with the channels below. Distinct values, never used as real rates.
_IDLE, _GATE, _MEAS, _RESET = 0.0123457, 0.0234568, 0.0345679, 0.0456781


def _is(x, tag):
    return abs(x - tag) < 1e-9


def idle_channel(tau: float, idle_depol: float, dephase_linear: float,
                 dephase_quad: float) -> tuple[float, float, float]:
    """(pX, pY, pZ) for one idle period of duration tau.

    Per unit of time: depolarizing at `idle_depol` (X, Y, Z at idle_depol/3 each), composed
    over tau units. Over the whole period, the emulator's two dephasing terms:
    linear       Z with probability dephase_linear * tau
    quadratic    RZ(dephase_quad * tau), Pauli-twirled: Z with probability sin^2(f*tau/2)
    Channels compose exactly through their Pauli eigenvalues (lambda_X, lambda_Y, lambda_Z).
    """
    e = (1 - 4 * idle_depol / 3) ** tau                # depolarizing, per axis
    pz_lin = min(dephase_linear * tau, 0.5)
    pz_quad = np.sin(dephase_quad * tau / 2) ** 2
    pz = pz_lin + pz_quad - 2 * pz_lin * pz_quad       # two independent Z flips, mod 2
    lx = ly = e * (1 - 2 * pz)                         # Z dephases the X and Y axes
    lz = e
    px = (1 + lx - ly - lz) / 4
    py = (1 - lx + ly - lz) / 4
    pz_tot = (1 - lx - ly + lz) / 4
    return tuple(float(max(v, 0.0)) for v in (px, py, pz_tot))


class StimMemory:
    """One code, one basis, one set of rates; circuits and decoders are built per (tau, R)."""

    def __init__(self, code: str, d: int, basis: str = "Z", p1: float = 0.0,
                 p2: float = 0.0, p_meas: float = 0.0, p_init: float = 0.0,
                 idle_depol: float = 0.0, dephase_linear: float = 0.0,
                 dephase_quad: float = 0.0):
        if (code, basis) not in SPEC:
            raise ValueError(f"unknown code/basis {code}/{basis}")
        self.code, self.d, self.basis = code, d, basis
        self.p1, self.p2, self.p_meas, self.p_init = p1, p2, p_meas, p_init
        self.idle = dict(idle_depol=idle_depol, dephase_linear=dephase_linear,
                         dephase_quad=dephase_quad)

    def circuit(self, tau: float, rounds: int) -> stim.Circuit:
        """Stim's generated memory circuit, with each tagged noise location rewritten."""
        tagged = stim.Circuit.generated(
            SPEC[(self.code, self.basis)], distance=self.d, rounds=rounds,
            before_round_data_depolarization=_IDLE, after_clifford_depolarization=_GATE,
            before_measure_flip_probability=_MEAS, after_reset_flip_probability=_RESET)
        px, py, pz = idle_channel(tau, **self.idle)

        def rewrite(block: stim.Circuit) -> stim.Circuit:
            out = stim.Circuit()
            for inst in block:
                if isinstance(inst, stim.CircuitRepeatBlock):
                    out.append(stim.CircuitRepeatBlock(inst.repeat_count,
                                                       rewrite(inst.body_copy())))
                    continue
                name, args, tg = inst.name, inst.gate_args_copy(), inst.targets_copy()
                if name == "DEPOLARIZE1" and _is(args[0], _IDLE):
                    out.append("PAULI_CHANNEL_1", tg, [px, py, pz])
                elif name == "DEPOLARIZE1" and _is(args[0], _GATE):
                    out.append("DEPOLARIZE1", tg, [self.p1])
                elif name == "DEPOLARIZE2" and _is(args[0], _GATE):
                    out.append("DEPOLARIZE2", tg, [self.p2])
                elif name in ("X_ERROR", "Z_ERROR") and _is(args[0], _MEAS):
                    out.append(name, tg, [self.p_meas])
                elif name in ("X_ERROR", "Z_ERROR") and _is(args[0], _RESET):
                    out.append(name, tg, [self.p_init])
                else:
                    out.append(inst)
            return out

        c = rewrite(tagged)
        left = [v for v in (_IDLE, _GATE, _MEAS, _RESET) if f"{v}" in str(c)]
        assert not left, f"untranslated noise tags remain: {left}"
        return c

    def logical_error(self, tau: float, rounds: int, shots: int, seed: int = 1):
        """(pL, sigma), decoded with PyMatching on the circuit's detector error model."""
        c = self.circuit(tau, rounds)
        dem = c.detector_error_model(decompose_errors=True)
        matcher = pymatching.Matching.from_detector_error_model(dem)
        det, obs = c.compile_detector_sampler(seed=seed).sample(
            shots, separate_observables=True)
        fails = int(np.sum(matcher.decode_batch(det)[:, 0] != obs[:, 0]))
        e = fails / shots
        return e, float(np.sqrt(max(e * (1 - e), 1e-12) / shots)), fails
