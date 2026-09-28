#!/bin/zsh
cd "/Users/ricardomaldonado/Documents/D-Blast 3/untitled folder 150/HDBLAST_CHAT13_NONLINEAR_ROLLOFF_20260921/referee_numerics"
P=python3
# group 1 (parallel): seed tests + KO variations + far boundary
$P rolloff5d_ref.py --tdet 0.1 --dc=-5e-5 --dz 2e-3 --L 12 --tf 9.5 --phimin -40 --tag runs/ref_minus_half > runs/ref_minus_half.log 2>&1 &
$P rolloff5d_ref.py --tdet 0.1 --dc 1e-4 --dz 2e-3 --L 12 --tf 12 --ko 0.02 --tag runs/ref_plus_ko002 > runs/ref_plus_ko002.log 2>&1 &
$P rolloff5d_ref.py --tdet 0.1 --dc 1e-4 --dz 2e-3 --L 12 --tf 12 --ko 0.1 --tag runs/ref_plus_ko01 > runs/ref_plus_ko01.log 2>&1 &
$P rolloff5d_ref.py --tdet 0.1 --dc 1e-4 --dz 2e-3 --L 16 --tf 16 --tag runs/ref_plus_L16 > runs/ref_plus_L16.log 2>&1 &
wait
# group 2
$P rolloff5d_ref.py --tdet 0.1 --dc 1e-4 --dz 2e-3 --L 12 --tf 12 --pslope 1 --tag runs/ref_plus_pslope > runs/ref_plus_pslope.log 2>&1 &
$P rolloff5d_ref.py --tdet 0.1 --dc=-1e-4 --dz 2e-3 --L 12 --tf 9 --phimin -40 --ko 0.1 --tag runs/ref_minus_ko01 > runs/ref_minus_ko01.log 2>&1 &
$P rolloff5d_ref.py --tdet 0.1 --dc=-1e-4 --dz 2e-3 --L 12 --tf 9 --phimin -40 --pslope 1 --tag runs/ref_minus_pslope > runs/ref_minus_pslope.log 2>&1 &
# t = 1e-3 resolution scan (small domain, linear stage only)
$P rolloff5d_ref.py --tdet 1e-3 --dc 1e-4 --dz 1e-3 --L 3 --tf 3 --tag runs/ref_t1e3_dz1e-3 > runs/ref_t1e3_dz1e-3.log 2>&1 &
wait
$P rolloff5d_ref.py --tdet 1e-3 --dc 1e-4 --dz 5e-4 --L 3 --tf 3 --tag runs/ref_t1e3_dz5e-4 > runs/ref_t1e3_dz5e-4.log 2>&1 &
$P rolloff5d_ref.py --tdet 1e-3 --dc 1e-4 --dz 2.5e-4 --L 3 --tf 3 --tag runs/ref_t1e3_dz2.5e-4 > runs/ref_t1e3_dz2.5e-4.log 2>&1 &
$P rolloff5d_ref.py --tdet 1e-3 --dc 1e-4 --dz 5e-4 --L 3 --tf 3 --pslope 1 --tag runs/ref_t1e3_dz5e-4_pslope > runs/ref_t1e3_dz5e-4_pslope.log 2>&1 &
wait
echo BATCH_DONE > runs/BATCH_DONE
