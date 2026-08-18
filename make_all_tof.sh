#!/usr/bin/env bash
cd "$(dirname "$0")"

BIN=./build/bin/MakeTOF
EVTDIR=evt
LOG=logs/tof_batch.log
mkdir -p logs

for run in $(seq 1000 1197); do
  padded=$(printf "%04d" "$run")
  rundir="${EVTDIR}/run${run}"

  if [ -f "tof/tof${padded}.tof" ]; then
    echo "run ${run}: already have a .tof, skipping" >> "$LOG"
    continue
  fi

  if [ ! -d "$rundir" ]; then
    echo "run ${run}: no evt directory found (${rundir}), skipping" >> "$LOG"
    continue
  fi

  echo "run ${run}: processing" >> "$LOG"
  "$BIN" "${rundir}"/run-"${run}"-*.evt >> "$LOG" 2>&1
done

echo "done" >> "$LOG"
