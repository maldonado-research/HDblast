#!/bin/bash
# usage: run_chain.sh <outdir> <tag> <evolve_a1.py args...>
# Runs one job in processes of at most --wall seconds (checkpointed); continues with --resume while the stop reason is the
# wall-clock limit and the accumulated runtime is below MAXCPU seconds (default 3600).  One core.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
out=$1; tag=$2; shift 2
mkdir -p "$out"
python3 -W ignore evolve_a1.py --out "$out" --tag "$tag" "$@" >> "$out/$tag.log" 2>&1
while true; do
  s="$out/${tag}_summary.json"
  [ -f "$s" ] || break
  cont=$(python3 -c "import json,sys; s=json.load(open('$s')); print(1 if s['stop_reason'].startswith('wall-clock') and s['runtime_s'] < ${MAXCPU:-3600} else 0)")
  [ "$cont" = "1" ] || break
  echo "--- resume $(date)" >> "$out/$tag.log"
  python3 -W ignore evolve_a1.py --out "$out" --tag "$tag" "$@" --resume >> "$out/$tag.log" 2>&1
done
