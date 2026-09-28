#!/bin/bash
# Shift-invert eigenvalue searches on finer grids (single thread; sequential).
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
run(){ timeout 1150 python3 spectra.py "$@" || echo "FAILED $*"; }
S=1.6572,1.2,0.5,0.1,-0.5,-1.0,-1.3,-1.45,-1.55,-2.0,-2.9,-4.0
for h in 2e-4 1e-4 5e-5; do
  run --bg plus --hmin $h --mode shift --k 6 --shifts $S --out results/spectra/shift_plus_h$h
  run --bg original --hmin $h --mode shift --k 6 --shifts 1.6572,1.0,0.1 --out results/spectra/shift_original_h$h
done
run --bg plus --hmin 1e-4 --L 8 --mode shift --k 6 --shifts $S --out results/spectra/shift_plus_h1e-4_L8
run --bg plus --hmin 1e-4 --tdet 3e-3 --mode shift --k 6 --shifts $S --out results/spectra/shift_plus_tdet3e-3_h1e-4
run --bg plus --hmin 1e-4 --tdet 1e-2 --mode shift --k 6 --shifts $S --out results/spectra/shift_plus_tdet1e-2_h1e-4
run --bg plus --hmin 1e-4 --mode shift --k 2 --s20-scale -1 --shifts 167.5 --out results/spectra/control_shift_plus_s20flip_h1e-4
run --bg plus --hmin 5e-5 --mode shift --k 2 --s20-scale -1 --shifts 167.5 --out results/spectra/control_shift_plus_s20flip_h5e-5
