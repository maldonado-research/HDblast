#!/bin/bash
# usage: launch.sh <outdir> <tag> <evolve_a1.py args...>   (background, one core, log next to the outputs)
cd "$(dirname "$0")"
out=$1; tag=$2; shift 2
mkdir -p "$out"
OMP_NUM_THREADS=1 nohup python3 -W ignore evolve_a1.py --out "$out" --tag "$tag" "$@" >> "$out/$tag.log" 2>&1 &
echo "launched $tag pid $!"
