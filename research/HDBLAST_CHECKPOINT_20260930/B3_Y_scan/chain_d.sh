#!/bin/bash
cd "$(dirname "$0")"
while [ ! -f chain_c.done ]; do sleep 20; done
echo "runs/main main_Y5_dc1e-2_dzf5e-4_cfl0.25_ext --dstar --Y 5 --dc 0.01 --dzf 5e-4 --dzc 2e-3 --L 10 --Tf 20 --kappa 10 --xc 8.7 --stop_recollapse --cfl 0.25 --save_T 0.25 --restart runs/main/states/main_Y5_dc1e-2_dzf5e-4_cfl0.25_state_T8.75.npz" > jobs_chain_d.txt
./queue.sh jobs_chain_d.txt > queue_chain_d.out 2>&1
touch chain_d.done
