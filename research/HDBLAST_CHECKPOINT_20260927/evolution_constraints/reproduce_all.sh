#!/bin/bash
# Regenerates every run used in README.md, sequentially on one core.
# Approximate single-core wall times (4-CPU Linux box, numpy 2.3.5 / scipy 1.16.3):
# h=4e-4: 20-30 s, 2e-4: 1.5 min, 1e-4: 5-6 min, 5e-5: 17-25 min per run to t=1.
# Total roughly 4-5 hours.  Run from this folder.  Each run writes <out>.json/.npz/.log.
set -e
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
r(){ out=$1; shift; mkdir -p "$(dirname "$out")"; [ -f "$out.json" ] || python3 -B run_case.py "$@" --out "$out" > "$out.log" 2>&1; }

# 0. the copied checkpoint wrapper itself (reproduction of the archived coarse wide run)
mkdir -p runs/baseline_original
[ -f runs/baseline_original/balanced_epsp01_wide_h2.json ] || python3 -B src/evolve_balanced.py --epsilon .01 --hmin .0002 \
   --stretch 2 --L 3 --tf 1 --output runs/baseline_original/balanced_epsp01_wide_h2 > runs/baseline_original/balanced_epsp01_wide_h2.log 2>&1

# 1. main matrix: epsilon=0.01, L=3, stretch=2 (the checkpoint's "wide" map), t in [0,1]
D="--kappa 10 --damp-off 0.02 0.05"
for h in 4e-4 2e-4 1e-4 5e-5; do
  r runs/final/base_h$h --hmin $h                       # baseline (order 4, no damping)
  r runs/final/proj_h$h --hmin $h --project             # remedy A: discrete-constraint projection
  r runs/final/k10_h$h  --hmin $h $D                    # remedy B: outgoing C+ damping, off near the shell
  r runs/final/k10p_h$h --hmin $h --project $D          # A + B
  r runs/final/o6p_h$h  --hmin $h --project --order 6   # remedy C: order-6 interior + projection
done

# 2. short (t<=0.1) diagnostics of the early near-shell generation at h=5e-5
for a in "base_h1e-4 --hmin 1e-4" "base_h5e-5 --hmin 5e-5" "xt_h5e-5 --hmin 5e-5 --xtol 1e-300" \
         "proj_h5e-5 --hmin 5e-5 --project" "ko02_h5e-5 --hmin 5e-5 --ko 0.2" "ko02p_h5e-5 --hmin 5e-5 --ko 0.2 --project" \
         "o6_h5e-5 --hmin 5e-5 --order 6" "o6p_h5e-5 --hmin 5e-5 --order 6 --project" "cfl02_h5e-5 --hmin 5e-5 --cfl 0.2"; do
  set -- $a; n=$1; shift; r runs/short/$n --tf 0.1 "$@"
done

# 3. damping controls at h=4e-4
r runs/damp/o4k10_h4e-4 --hmin 4e-4 --kappa 10                        # uniform kappa (shell instability)
r runs/ctrl/o4km5_h4e-4 --hmin 4e-4 --kappa -5                        # anti-damping control
r runs/ctrl/o4honly10_h4e-4 --hmin 4e-4 --kappa 10 --damp-form honly  # naive H-only damping (unstable)
r runs/damp/k30off_h4e-4 --hmin 4e-4 --kappa 30 --damp-off 0.02 0.05

# 4. calibrations
for h in 4e-4 2e-4 1e-4; do   # uncompensated single pulse (corner-incompatible control)
  r runs/calib/gauss_h$h --hmin $h --epsilon 0 --gauss 0.01 -0.3 0.05 --project
done
for h in 4e-4 2e-4 1e-4 5e-5; do  # smooth corner-free bulk pulse, width 0.02 (order calibration)
  [ -f runs/calib3/g_h$h.json ] || python3 -B calib.py --w 0.02 --hmin $h --out runs/calib3/g_h$h > runs/calib3/g_h$h.log 2>&1
done
for h in 4e-4 2e-4 1e-4; do       # same with width 0.05 (reaches the floor ~1e-10)
  [ -f runs/calib2/g_h$h.json ] || python3 -B calib.py --w 0.05 --hmin $h --out runs/calib2/g_h$h > runs/calib2/g_h$h.log 2>&1
done

# 5. controls and analysis
python3 -B controls.py > CONTROLS.log 2>&1
python3 -B analyze.py > ANALYSIS.log 2>&1
