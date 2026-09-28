#!/bin/bash
# Audit re-runs of the audited scripts (unchanged), outputs redirected into verify/rerun/.
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
O=verify/rerun
S=1.6572,1.2,0.5,0.1,-0.5,-1.0,-1.3,-1.45,-1.55,-2.0,-2.9,-4.0
python3 spectra.py --bg plus --hmin 4e-4 --L 6 --mode dense --out $O/plus_h4_L6
python3 spectra.py --bg plus --hmin 5e-5 --mode shift --k 6 --shifts $S --out $O/shift_plus_h5e-5
python3 spectra.py --bg original --hmin 5e-5 --mode shift --k 6 --shifts 1.6572,1.0,0.1 --out $O/shift_original_h5e-5
python3 spectra.py --bg plus --hmin 5e-5 --mode shift --k 2 --s20-scale -1 --shifts 167.5 --out $O/control_shift_plus_s20flip_h5e-5
python3 spectra.py --bg plus --hmin 1e-4 --tdet 1e-2 --mode shift --k 6 --shifts $S --out $O/shift_plus_tdet1e-2_h1e-4
python3 td_evolve.py --bg plus --hmin 4e-4 --tf 6 --linear --out $O/plus_lin_h4e-4
python3 td_evolve.py --bg plus --hmin 2e-4 --tf 6 --linear --out $O/plus_lin_h2e-4
python3 td_evolve.py --bg plus --hmin 1e-4 --tf 6 --linear --out $O/plus_lin_h1e-4
python3 td_evolve.py --bg plus --hmin 5e-5 --tf 6 --linear --out $O/plus_lin_h5e-5
