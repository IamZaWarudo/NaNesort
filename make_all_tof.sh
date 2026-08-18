#!/usr/bin/env bash
cd "$(dirname "$0")"
shopt -s nullglob

BIN=./build/bin/MakeTOF
EVTDIR=../data
LOG=logs/tof_batch.log
mkdir -p logs

for run in $(seq 1000 1197); do
  padded=$(printf "%04d" "$run")
  files=("${EVTDIR}"/run-"${run}"-*.evt)

  if [ -f "tof/tof${padded}.tof" ]; then
    echo "run ${run}: already have a .tof, skipping" >> "$LOG"
    continue
  fi

  if [ ${#files[@]} -eq 0 ]; then
    echo "run ${run}: no evt files found, skipping" >> "$LOG"
    continue
  fi

  echo "run ${run}: processing ${#files[@]} file(s)" >> "$LOG"
  "$BIN" "${files[@]}" >> "$LOG" 2>&1
done

echo "done" >> "$LOG"
