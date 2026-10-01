#!/bin/bash
# A1 tuned-model runs (MODEL CHANGE: sigma = 2W + delta(1 + c phi + d* phi^2/2)), delta = 0.1.  Two runs in parallel (xargs -P 2).
# Stage 1: coarse old-chart pre-runs set the chart parameter xc (REGISTRATION.md section 4).
# Stage 2: main runs in the bounded proper-clock chart, two shell spacings, two seeds, Y in {0.3, 1, 3} (+ Y = 0 baseline).
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
P=${P:-2}
stage=${1:-all}
YS="0 0.3 1 3"; DCS="1e-2 1e-4"
if [ "$stage" = pre ] || [ "$stage" = all ]; then
  for Y in $YS; do for dc in $DCS; do
    echo "--dstar --Y $Y --dc $dc --dzf 2e-3 --dzc 8e-3 --L 12 --Tf 20 --prerun --out runs/pre --tag pre_dstar_Y${Y}_dc${dc}"
  done; done | xargs -P $P -I{} sh -c 'python3 -W ignore evolve_a1.py {} > runs/pre/$(echo "{}" | sed "s/.*--tag //").log 2>&1'
fi
if [ "$stage" = main ] || [ "$stage" = all ]; then
  for sp in "1e-3 4e-3" "5e-4 2e-3"; do set -- $sp
    for Y in $YS; do for dc in $DCS; do
      echo "--dstar --Y $Y --dc $dc --dzf $1 --dzc $2 --L 16 --Tf 14 --kappa 10 --project 15.9 --xc_from runs/pre/pre_dstar_Y${Y}_dc${dc}_summary.json --out runs/main --tag main_dstar_Y${Y}_dc${dc}_dzf$1"
    done; done
  done | xargs -P $P -I{} sh -c 'python3 -W ignore evolve_a1.py {} > runs/main/$(echo "{}" | sed "s/.*--tag //").log 2>&1'
fi
