"""Stage 7a: upload the QIR to Nexus, estimate its HQC cost, and (only with
--submit) run it on the H2-Emulator. Run in .venv-h2 ONLY. Requires a Nexus
login first:  .venv-h2\\Scripts\\qnx.exe login

What it does (qnexus 0.51.1 APIs, checked in the installed source):
1. ``qnx.projects.get_or_create(name=--project)``
2. ``qnx.qir.upload(qir=<bitcode>, name, project)``
3. ``qnx.qir.cost(programs, n_shots, project, system_name=--cost-system)``
   -> prints the HQC estimate. NOTE: ``qnx.circuits.cost`` only takes pytket
   CircuitRefs; for QIR the cost call is ``qnx.qir.cost``. It runs a costing
   job on the "<system>SC" syntax checker (e.g. H2-1SC), which is visible in
   the Nexus portal. qnexus marks qir.cost as deprecated for Helios-or-newer
   systems only; H2 is older.
4. Only with ``--submit``: ``qnx.start_execute_job`` on ``--device`` (H2-1E)
   with ``QuantinuumConfig(simulator="stabilizer", noisy_simulation=True,
   error_params=UserErrorParams(idle_dephasing=False))``. Idle dephasing is
   off because the program injects its own dephasing (RNG-conditioned Z).
   Transport dephasing stays at the emulator default unless
   ``--transport-dephasing-off``. Shot counts above --max-shots are refused.

Usage (from the repo root):
    .venv-h2\\Scripts\\python.exe scripts\\stage7a_submit.py results\\stage7a\\d3_r2_X.bc
    .venv-h2\\Scripts\\python.exe scripts\\stage7a_submit.py results\\stage7a\\d3_r2_X.bc --submit --shots 100 --wait
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import sys
from pathlib import Path

import qnexus as qnx
from qnexus.models import QuantinuumConfig, UserErrorParams

ROOT = Path(__file__).resolve().parents[1]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("qir", type=Path, help="QIR bitcode (.bc) from stage7a_build_qir.py")
    ap.add_argument("--project", default="qec-timing-stage7a")
    ap.add_argument("--name", default=None, help="Nexus program name (default: file stem)")
    ap.add_argument("--shots", type=int, default=100)
    ap.add_argument("--max-shots", type=int, default=1000)
    ap.add_argument("--cost-system", default="H2-1", help="costed on <system>SC")
    ap.add_argument("--device", default="H2-1E")
    ap.add_argument("--transport-dephasing-off", action="store_true")
    ap.add_argument("--submit", action="store_true", help="actually run on the emulator")
    ap.add_argument("--wait", action="store_true", help="with --submit: wait and save results")
    args = ap.parse_args(argv)

    if args.shots < 1 or args.shots > args.max_shots:
        raise SystemExit(f"--shots must be in [1, {args.max_shots}] (raise --max-shots explicitly)")
    qir = args.qir.read_bytes()
    meta_path = args.qir.with_suffix(".json")
    meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
    name = args.name or args.qir.stem
    stamp = _dt.datetime.now(_dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")

    try:
        project = qnx.projects.get_or_create(name=args.project)
    except Exception as e:  # noqa: BLE001
        raise SystemExit(f"Nexus call failed ({e!r}). Log in first: .venv-h2\\Scripts\\qnx.exe login")
    qir_ref = qnx.qir.upload(qir=qir, name=f"{name}-{stamp}", project=project,
                             description="Stage 7a surface-code feasibility test",
                             properties={"d": meta.get("d", -1), "rounds": meta.get("rounds", -1)})
    print(f"uploaded {args.qir} ({len(qir)} B) as {qir_ref.id} in project {args.project!r}")

    hqc = qnx.qir.cost(programs=qir_ref, n_shots=args.shots, project=project,
                       system_name=args.cost_system)
    print(f"HQC estimate on {args.cost_system}SC for {args.shots} shots: {hqc}")

    record = dict(qir=str(args.qir), qir_id=str(qir_ref.id), project=args.project,
                  shots=args.shots, cost_system=args.cost_system, hqc=hqc, submitted=False,
                  meta=meta)
    if not args.submit:
        print("not submitted (pass --submit to run on the emulator)")
    else:
        params = UserErrorParams(idle_dephasing=False,
                                 **({"transport_dephasing": False}
                                    if args.transport_dephasing_off else {}))
        config = QuantinuumConfig(device_name=args.device, simulator="stabilizer",
                                  noisy_simulation=True, error_params=params)
        job = qnx.start_execute_job(programs=[qir_ref], n_shots=[args.shots],
                                    backend_config=config, name=f"{name}-{stamp}",
                                    project=project)
        print(f"submitted job {job.id} on {args.device} ({args.shots} shots)")
        record.update(submitted=True, job_id=str(job.id), device=args.device,
                      error_params=params.model_dump(exclude_none=True))
        if args.wait:
            qnx.jobs.wait_for(job)
            result = qnx.jobs.results(job)[0].download_result()
            out = ROOT / "results" / "stage7a" / f"{name}-{stamp}-result.json"
            payload = result.to_dict() if hasattr(result, "to_dict") else str(result)
            out.write_text(json.dumps(payload, indent=2, default=str))
            print(f"results -> {out}")
    log = ROOT / "results" / "stage7a" / f"{name}-{stamp}-nexus.json"
    log.write_text(json.dumps(record, indent=2, default=str))
    print(f"record -> {log}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
