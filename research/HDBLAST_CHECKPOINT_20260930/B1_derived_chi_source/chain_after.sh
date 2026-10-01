#!/bin/bash
# runs after the grid queue: exploratory cell, then the blow-up diagnostics (one process at a time; 1 core)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
while pgrep -f "queue.sh jobs_grid.txt" > /dev/null; do sleep 10; done
./queue.sh jobs_explore.txt > queue_explore.out 2>&1
mkdir -p diag
python3 -W ignore diag_blowup.py diag_ps0.5_G100_y1_bneed_dzf1e-3 0.05 -- --dstar --source chi --dc 1e-2 --phistar 0.5 --G 100 --y 1 --b 0.010 --dzf 1e-3 --dzc 4e-3 --L 16 --Tf 14 --kappa 10 --project 15.9 --xc_from runs/pre/pre_ps0.5_G100_y1_bneed_summary.json > diag/diag_bneed.log 2>&1
python3 -W ignore diag_blowup.py diag_Y0_A1equiv_dzf1e-3 0.05 -- --dstar --source friction --Y 0 --dc 1e-2 --dzf 1e-3 --dzc 4e-3 --L 16 --Tf 14 --kappa 10 --project 15.9 --xc_from /home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260928/A1_tuned_vacuum/runs/pre/pre_dstar_Y0_dc1e-2_summary.json > diag/diag_Y0.log 2>&1
echo done > chain_after.done
