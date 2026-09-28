#!/bin/bash
# Time-domain batch B: the +1 static branch.  Single thread, sequential.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
run(){ timeout 1150 python3 td_evolve.py "$@" || echo "FAILED $*"; }
D=results/td
for h in 4e-4 2e-4 1e-4 5e-5; do run --bg plus --hmin $h --tf 6 --linear --out $D/plus_lin_h$h; done
# Shell-localised initial disturbance (different initial profile).
for h in 2e-4 1e-4 5e-5; do run --bg plus --hmin $h --tf 6 --linear --z1 -0.01 --w1 0.008 --z2 -0.3 --w2 0.1 --out $D/plus_lin_bumpB_h$h; done
# Longer domain.
run --bg plus --hmin 1e-4 --L 10 --tf 6 --linear --out $D/plus_lin_L10_h1e-4
# Unchanged nonlinear RHS, two amplitudes, two resolutions.
run --bg plus --hmin 1e-4 --tf 6 --eps 1e-6 --out $D/plus_nl_eps1e-6_h1e-4
run --bg plus --hmin 1e-4 --tf 6 --eps 2e-6 --out $D/plus_nl_eps2e-6_h1e-4
run --bg plus --hmin 2e-4 --tf 6 --eps 1e-6 --out $D/plus_nl_eps1e-6_h2e-4
