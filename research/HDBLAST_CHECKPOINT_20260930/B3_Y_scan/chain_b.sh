#!/bin/bash
# runs after the main queue: restarts + controls, dc = 1e-4 repeats (Y = 2, 3), then Y = 5 (cfl 0.25)
cd "$(dirname "$0")"
while pgrep -f "queue.sh jobs_main.txt" > /dev/null; do sleep 20; done
./queue.sh jobs_chain2.txt > queue_chain2.out 2>&1
./queue.sh jobs_dc4_pre.txt > queue_dc4_pre.out 2>&1
./queue.sh jobs_dc4_main.txt > queue_dc4_main.out 2>&1
M="--L 10 --Tf 20 --kappa 10 --project 9.9 --stop_recollapse --cfl 0.25"
: > jobs_y5_main.txt
echo "runs/main main_Y5_dc1e-2_dzf1e-3_cfl0.25 --dstar --Y 5 --dc 1e-2 --dzf 1e-3 --dzc 4e-3 $M --save_T 0.25 --xc_from runs/pre/pre_dstar_Y5_dc1e-2_cfl0.25_summary.json" >> jobs_y5_main.txt
echo "runs/main main_Y5_dc1e-2_dzf5e-4_cfl0.25 --dstar --Y 5 --dc 1e-2 --dzf 5e-4 --dzc 2e-3 $M --save_T 0.25 --xc_from runs/pre/pre_dstar_Y5_dc1e-2_cfl0.25_summary.json" >> jobs_y5_main.txt
./queue.sh jobs_y5_main.txt > queue_y5_main.out 2>&1
touch chain_b.done
