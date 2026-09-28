#!/bin/bash
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
O=verify/rerun
python3 td_evolve.py --bg original --hmin 2e-4 --tf 5 --linear --out $O/original_lin_h2e-4
python3 td_evolve.py --bg original --hmin 1e-4 --tf 5 --linear --out $O/original_lin_h1e-4
python3 td_evolve.py --bg original --hmin 5e-5 --tf 5 --linear --out $O/original_lin_h5e-5
python3 td_evolve.py --bg plus --hmin 1e-4 --tf 0.2 --linear --s20-scale -1 --max-f 1e300 --out $O/control_plus_s20flip_lin_h1e-4
python3 td_evolve.py --bg plus --hmin 1e-4 --tf 6 --linear --z1 -0.01 --w1 0.008 --z2 -0.3 --w2 0.1 --out $O/plus_lin_bumpB_h1e-4
