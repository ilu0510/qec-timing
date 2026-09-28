"""Stim circuits for the phenomenological rotated-surface-code memory.

Structure per Appendix A, Z-check sector only, memory prepared in logical |0>:

1. For each round ``j``: ``X_ERROR(p_data_j)`` on every data qubit, then an
   *ideal* Z-stabilizer measurement whose recorded outcome is flipped with
   probability ``p_read``, then detectors ``D_{j,c} = s_{j,c} XOR s_{j-1,c}``
   (Eq. A1). Round 0 compares against the deterministic initial value, so its
   detectors reference a single measurement.
2. A final data readout with classical flip probability ``p_read``, final
   detectors comparing each check's support parity against its last syndrome,
   and the logical Z observable.

The stabilizer measurement is ideal in the sense required by CLAUDE.md: MPP
flips only the classical record, never the state. Verified in Stage 0 --
``MPP(p)`` reports a flipped outcome with probability ``p`` while a following
noiseless measurement is undisturbed.

Note on the final readout: no ``X_ERROR`` is applied before it, so the total
elapsed time is exactly ``sum_j dt_j``. Idling noise is charged once per round.
"""

from __future__ import annotations

import stim

from ..layout import RotatedSurfaceCodeZSector
from ..schedule import RoundSchedule
from .noise import data_probabilities, p_read as p_read_from

__all__ = [
    "build_memory_circuit",
    "build_memory_circuit_from_schedule",
    "data_flip_probe",
    "readout_flip_probe",
    "single_error_circuit",
]


def _mpp_targets(layout: RotatedSurfaceCodeZSector, support) -> list:
    """Z-basis product-measurement targets for one check."""
    targets: list = []
    for k, q in enumerate(support):
        if k:
            targets.append(stim.target_combiner())
        targets.append(stim.target_z(layout.data_index[q]))
    return targets


def build_memory_circuit(
    layout: RotatedSurfaceCodeZSector,
    data_probabilities: tuple[float, ...] | list[float],
    p_read: float,
    *,
    with_coords: bool = True,
) -> stim.Circuit:
    """Build the memory-experiment circuit from explicit per-round probabilities.

    Args:
        layout: the Z-sector layout.
        data_probabilities: ``p_data_j`` for each round; its length sets the
            number of rounds.
        p_read: syndrome and final-readout flip probability.
        with_coords: annotate detectors with ``(row, col, round)`` coordinates.

    Taking probabilities directly -- rather than ``(p, lam, b_read)`` -- lets the
    data-fault and readout-fault channels be exercised independently, which is
    what the rate tests need, since ``p_data`` and ``p_read`` both scale with p.
    """
    if not data_probabilities:
        raise ValueError("need at least one round")
    for j, pd in enumerate(data_probabilities):
        if not 0.0 <= pd <= 1.0:
            raise ValueError(f"data_probabilities[{j}] = {pd} is not a probability")
    if not 0.0 <= p_read <= 1.0:
        raise ValueError(f"p_read = {p_read} is not a probability")

    checks = layout.z_checks
    n_checks = len(checks)
    n_data = layout.n_data
    data_targets = [layout.data_index[q] for q in layout.data_qubits]

    circuit = stim.Circuit()
    if with_coords:
        for ci in range(n_checks):
            row, col = layout.check_coords(ci)
            circuit.append("QUBIT_COORDS", [n_data + ci], [row, col])

    for j, pd in enumerate(data_probabilities):
        circuit.append("X_ERROR", data_targets, pd)
        for support in checks:
            circuit.append("MPP", _mpp_targets(layout, support), p_read)
        for ci in range(n_checks):
            recs = [stim.target_rec(-n_checks + ci)]
            if j > 0:
                # Eq. A1: compare with the same check in the previous round.
                recs.append(stim.target_rec(-2 * n_checks + ci))
            args = [*layout.check_coords(ci), j] if with_coords else []
            circuit.append("DETECTOR", recs, args)
        circuit.append("TICK")

    # Final noisy data readout.
    circuit.append("M", data_targets, p_read)
    n_rounds = len(data_probabilities)
    for ci, support in enumerate(checks):
        recs = [stim.target_rec(-n_data + layout.data_index[q]) for q in support]
        recs.append(stim.target_rec(-n_data - n_checks + ci))
        args = [*layout.check_coords(ci), n_rounds] if with_coords else []
        circuit.append("DETECTOR", recs, args)

    circuit.append(
        "OBSERVABLE_INCLUDE",
        [stim.target_rec(-n_data + layout.data_index[q]) for q in layout.logical_z_support],
        0,
    )
    return circuit


def build_memory_circuit_from_schedule(
    layout: RotatedSurfaceCodeZSector,
    schedule: RoundSchedule,
    p: float,
    b_read: float,
    *,
    with_coords: bool = True,
) -> stim.Circuit:
    """Build the circuit from a schedule of ``dt_j``/``lam_j`` and ``(p, b_read)``.

    Derives ``p_data_j`` via Eqs. 1-2 and ``p_read`` via Eq. 3.
    """
    return build_memory_circuit(
        layout,
        data_probabilities(schedule, p),
        p_read_from(p, b_read),
        with_coords=with_coords,
    )


# ---------------------------------------------------------------------------
# Probe circuits: isolate one noise channel so its rate can be measured.
# ---------------------------------------------------------------------------


def data_flip_probe(
    layout: RotatedSurfaceCodeZSector, p_data_value: float
) -> stim.Circuit:
    """One ``X_ERROR(p_data)`` on all data qubits, then a noiseless readout.

    Each measurement is 1 exactly when that data qubit was flipped, so the mean
    over qubits and shots estimates ``p_data`` directly.
    """
    circuit = stim.Circuit()
    targets = [layout.data_index[q] for q in layout.data_qubits]
    circuit.append("X_ERROR", targets, p_data_value)
    circuit.append("M", targets, 0.0)
    return circuit


def readout_flip_probe(
    layout: RotatedSurfaceCodeZSector, p_read_value: float, n_rounds: int = 1
) -> stim.Circuit:
    """Noisy Z-check measurements with **no** data faults.

    Every Z-stabilizer is deterministically +1 on the all-zero state, so a
    reported 1 is exactly a readout flip and the mean estimates ``p_read``.
    """
    circuit = stim.Circuit()
    for _ in range(n_rounds):
        for support in layout.z_checks:
            circuit.append("MPP", _mpp_targets(layout, support), p_read_value)
    return circuit


def single_error_circuit(
    layout: RotatedSurfaceCodeZSector,
    qubit: tuple[int, int],
    *,
    n_rounds: int = 3,
    round_index: int = 1,
) -> stim.Circuit:
    """Noiseless circuit with one deterministic X fault on ``qubit``.

    Used to check that any single data fault is corrected. All stochastic
    channels are zero, so the only error is the injected ``X_ERROR(1.0)``.
    """
    if not 0 <= round_index < n_rounds:
        raise ValueError(f"round_index must lie in [0, {n_rounds})")

    checks = layout.z_checks
    n_checks = len(checks)
    n_data = layout.n_data
    data_targets = [layout.data_index[q] for q in layout.data_qubits]

    circuit = stim.Circuit()
    for j in range(n_rounds):
        if j == round_index:
            circuit.append("X_ERROR", [layout.data_index[qubit]], 1.0)
        for support in checks:
            circuit.append("MPP", _mpp_targets(layout, support), 0.0)
        for ci in range(n_checks):
            recs = [stim.target_rec(-n_checks + ci)]
            if j > 0:
                recs.append(stim.target_rec(-2 * n_checks + ci))
            circuit.append("DETECTOR", recs)
    circuit.append("M", data_targets, 0.0)
    for ci, support in enumerate(checks):
        recs = [stim.target_rec(-n_data + layout.data_index[q]) for q in support]
        recs.append(stim.target_rec(-n_data - n_checks + ci))
        circuit.append("DETECTOR", recs)
    circuit.append(
        "OBSERVABLE_INCLUDE",
        [stim.target_rec(-n_data + layout.data_index[q]) for q in layout.logical_z_support],
        0,
    )
    return circuit
