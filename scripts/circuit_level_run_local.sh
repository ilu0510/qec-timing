#!/usr/bin/env bash
# Run circuit_level_generate.py for one distance on this machine, detached-safe.
#
#   setsid nohup scripts/circuit_level_run_local.sh <distance> <blocks> <workers> [python] &
#
# Writes results/circuit_level/logs/<run>/: verify.log, worker<i>.log, STATUS
# STATUS is RUNNING, VERIFY_FAILED, DONE, or FAILED (some worker exited non-zero).
set -u
D=$1; BLOCKS=$2; WORKERS=$3; PY=${4:-python}
REPO=$(cd "$(dirname "$0")/.." && pwd)
cd "$REPO"
LOG=results/circuit_level/logs/d${D}_$(date +%Y%m%d_%H%M%S)
mkdir -p "$LOG"
echo RUNNING > "$LOG/STATUS"
echo "d=$D blocks=$BLOCKS workers=$WORKERS started $(date -Is) pid $$" > "$LOG/run.info"

if ! "$PY" scripts/circuit_level_generate.py --distance "$D" --verify-only > "$LOG/verify.log" 2>&1; then
    echo VERIFY_FAILED > "$LOG/STATUS"; exit 1
fi

pids=()
for ((w = 0; w < WORKERS; w++)); do
    "$PY" scripts/circuit_level_generate.py --distance "$D" --blocks "$BLOCKS" \
        --worker "$w" --workers "$WORKERS" > "$LOG/worker$w.log" 2>&1 &
    pids+=($!)
done
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=1; done
if [ "$fail" = 0 ]; then echo DONE > "$LOG/STATUS"; else echo FAILED > "$LOG/STATUS"; fi
echo "finished $(date -Is)" >> "$LOG/run.info"
