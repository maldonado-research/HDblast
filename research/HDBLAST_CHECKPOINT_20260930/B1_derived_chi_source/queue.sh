#!/bin/bash
# usage: queue.sh <jobfile>   each line: <outdir> <tag> <evolve_b1.py args...>; runs ONE job at a time (1 core assigned)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
grep -v '^#' "$1" | grep -v '^$' | while read -r out tag args; do
  if [ -f "$out/${tag}_summary.json" ]; then echo "skip $tag (done)"; continue; fi
  mkdir -p "$out"
  echo "$(date -u +%H:%M:%S) start $tag"
  python3 -W ignore evolve_b1.py --out "$out" --tag "$tag" $args > "$out/$tag.log" 2>&1
  echo "$(date -u +%H:%M:%S) end $tag"
done
