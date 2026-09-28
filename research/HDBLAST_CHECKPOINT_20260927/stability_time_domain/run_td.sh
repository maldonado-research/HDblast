#!/bin/bash
# Time-domain batch A: calibration on the original unstable shell + wrong-formula/perturbed controls.
# Single thread, sequential.  Run together with run_td2.sh (batch B) on two cores.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
run(){ timeout 1150 python3 td_evolve.py "$@" || echo "FAILED $*"; }
D=results/td
for h in 2e-4 1e-4 5e-5; do run --bg original --hmin $h --tf 5 --linear --out $D/original_lin_h$h; done
run --bg original --hmin 1e-4 --tf 5 --eps 1e-8 --out $D/original_nl_eps1e-8_h1e-4
run --bg original --hmin 5e-5 --tf 5 --eps 1e-8 --out $D/original_nl_eps1e-8_h5e-5
# Wrong-formula control: sign of sigma''(phi_b) flipped in the scalar junction (+1 background).
for h in 2e-4 1e-4; do run --bg plus --hmin $h --tf 0.2 --linear --s20-scale -1 --max-f 1e300 --out $D/control_plus_s20flip_lin_h$h; done
# Control: sigma'' term removed (Neumann-like scalar junction).
run --bg plus --hmin 1e-4 --tf 6 --linear --s20-scale 0 --out $D/control_plus_s20zero_lin_h1e-4
# Perturbed parameter: +1 branch at detuning 0.003.
run --bg plus --hmin 1e-4 --tf 6 --linear --tdet 3e-3 --out $D/plus_tdet3e-3_lin_h1e-4
