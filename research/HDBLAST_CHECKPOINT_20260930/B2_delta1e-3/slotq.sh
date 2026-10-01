#!/bin/bash
# usage: slotq.sh <jobfile>  launches each job (via run_chain.sh) when fewer than P (default 2) evolve processes of THIS folder run
cd "$(dirname "$0")"; here=$(pwd)
mine() { n=0; for p in $(pgrep -f "evolve_a1.py --out"); do [ "$(readlink /proc/$p/cwd)" = "$here" ] && n=$((n+1)); done; echo $n; }
grep -v '^#' "$1" | grep -v '^$' | while read -r line; do
  while [ "$(mine)" -ge "${P:-2}" ]; do sleep 10; done
  echo "$(date) launch: $line"
  nohup ./run_chain.sh $line > /dev/null 2>&1 &
  sleep 20
done
