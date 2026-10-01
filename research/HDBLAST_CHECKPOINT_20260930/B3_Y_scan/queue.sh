#!/bin/bash
# usage: queue.sh <jobfile>   each line: <outdir> <tag> <evolve_a1.py args...> ; ONE job at a time (1 core assigned)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
grep -v '^#' "$1" | grep -v '^$' | xargs -P ${P:-1} -L 1 sh -c 'out=$0; tag=$1; shift 1; mkdir -p $out; python3 -W ignore evolve_a1.py --out $out --tag $tag "$@" > $out/$tag.log 2>&1'
