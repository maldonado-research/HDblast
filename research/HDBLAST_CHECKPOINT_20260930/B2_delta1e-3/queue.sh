#!/bin/bash
# usage: queue.sh <jobfile>   each line: <outdir> <tag> <evolve_a1.py args...> ; P (default 2) jobs at a time, one core each,
# each job through run_chain.sh (checkpointed processes of <= --wall seconds, continued with --resume)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
grep -v '^#' "$1" | grep -v '^$' | xargs -P ${P:-2} -L 1 ./run_chain.sh
