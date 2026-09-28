#!/bin/bash
# Re-runs a cheap decisive subset of the workstream's runs into verify/rerun/ (one core).
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
RC=../run_case.py
r(){ out=$1; shift; [ -f "$out.json" ] || timeout 1200 python3 -B $RC "$@" --out "$out" > "$out.log" 2>&1; echo "$out exit $?" >> rerun/exit_codes.txt; }
D="--kappa 10 --damp-off 0.02 0.05"
r rerun/base_h4e-4 --hmin 4e-4
r rerun/proj_h4e-4 --hmin 4e-4 --project
r rerun/k10_h4e-4  --hmin 4e-4 $D
r rerun/k10p_h4e-4 --hmin 4e-4 --project $D
r rerun/o6p_h4e-4  --hmin 4e-4 --project --order 6
r rerun/km5_h4e-4  --hmin 4e-4 --kappa -5
r rerun/proj_h2e-4 --hmin 2e-4 --project
r rerun/short_base_h5e-5 --hmin 5e-5 --tf 0.1
r rerun/short_proj_h5e-5 --hmin 5e-5 --tf 0.1 --project
echo DONE >> rerun/exit_codes.txt
