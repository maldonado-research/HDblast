#!/bin/bash
# Dense discrete-operator spectra (single thread each; run sequentially).
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
run(){ timeout 1150 python3 spectra.py "$@" || echo "FAILED $*"; }
run --bg plus --hmin 4e-4 --L 6 --mode dense --out results/spectra/plus_h4_L6
run --bg original --hmin 4e-4 --L 6 --mode dense --out results/spectra/original_h4_L6
run --bg plus --hmin 3e-4 --L 6 --mode dense --out results/spectra/plus_h3_L6
run --bg original --hmin 3e-4 --L 6 --mode dense --out results/spectra/original_h3_L6
run --bg plus --hmin 4e-4 --L 6 --mode dense --phi-poly --out results/spectra/plus_h4_L6_phipoly
run --bg plus --hmin 4e-4 --L 8 --mode dense --out results/spectra/plus_h4_L8
run --bg plus --hmin 4e-4 --L 10 --mode dense --out results/spectra/plus_h4_L10
run --bg plus --hmin 4e-4 --L 6 --stretch 0.1 --mode dense --out results/spectra/plus_h4_L6_st10
run --bg plus --hmin 4e-4 --L 6 --mode dense --s20-scale -1 --out results/spectra/control_plus_s20flip_h4
run --bg plus --hmin 4e-4 --L 6 --mode dense --s20-scale 0 --out results/spectra/control_plus_s20zero_h4
run --bg plus --hmin 2e-4 --L 6 --mode dense --out results/spectra/plus_h2_L6
run --bg original --hmin 2e-4 --L 6 --mode dense --out results/spectra/original_h2_L6
