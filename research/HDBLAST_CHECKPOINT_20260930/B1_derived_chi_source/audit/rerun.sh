#!/bin/bash
# audit reruns (1 core): exact rerun of one b_cut grid run, plus a third spacing
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
A="--dstar --source chi --dc 1e-2 --phistar 0.5 --G 100 --y 1 --b 6.49e-6 --L 16 --Tf 14 --kappa 10 --project 15.9 --xc_from runs/pre/pre_ps0.5_G100_y1_bcut_summary.json"
python3 -W ignore evolve_b1.py --out audit/runs --tag rerun_ps0.5_G100_y1_bcut_dzf1e-3 $A --dzf 1e-3 --dzc 4e-3 > audit/runs/rerun1.log 2>&1
python3 -W ignore evolve_b1.py --out audit/runs --tag rerun_ps0.5_G100_y1_bcut_dzf7.5e-4 $A --dzf 7.5e-4 --dzc 3e-3 > audit/runs/rerun2.log 2>&1
python3 -W ignore evolve_b1.py --out audit/runs --tag rerun_ps0.5_G100_y1_bcut_wrongsign_dzf1e-3 $A --dzf 1e-3 --dzc 4e-3 --signj -1 > audit/runs/rerun3.log 2>&1
echo done > audit/runs/DONE
