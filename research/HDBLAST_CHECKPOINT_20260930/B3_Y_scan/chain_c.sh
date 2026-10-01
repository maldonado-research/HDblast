#!/bin/bash
# after chain_b: Y = 2 dc = 1e-4 secondary restart (if triggered), Y = 3 tertiary hybrid, Y = 3 dc = 1e-4 pair
cd "$(dirname "$0")"
while [ ! -f chain_b.done ]; do sleep 20; done
export OMP_NUM_THREADS=1
: > jobs_chain_c.txt
python3 secondary_trigger.py > /dev/null 2>&1
if python3 -c "import json,sys; c=json.load(open('SECONDARY_TRIGGER.json'))['cases'].get('main_Y2_dc1e-4',{}); sys.exit(0 if c.get('trigger') else 1)"; then
  python3 make_restart_job.py main_Y2_dc1e-4_dzf5e-4 runs/fine fine_Y2_dc1e-4_dzf2.5e-4 2.5e-4 1e-3 >> jobs_chain_c.txt 2>>chain_c.err
fi
echo "runs/fine fine_Y3_dc1e-2_dzf1.25e-4 --dstar --Y 3 --dc 0.01 --dzf 1.25e-4 --dzc 1e-3 --L 10 --Tf 20 --kappa 10 --xc 6.1 --stop_recollapse --cfl 0.5 --wall 1500 --restart runs/main/states/main_Y3_dc1e-2_dzf5e-4_state_T6.25.npz" >> jobs_chain_c.txt
cat jobs_dc4_main_Y3.txt >> jobs_chain_c.txt
./queue.sh jobs_chain_c.txt > queue_chain_c.out 2>&1
touch chain_c.done
