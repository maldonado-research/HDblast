import sys, time, math, cProfile, pstats
sys.path.insert(0, '.'); sys.path.insert(0, '/home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927/mechanisms/pilot_5d')
import numpy as np, evolve_a1 as E
run = E.Run(delta=1e-3, d=-3.1942416958680835, dc=1e-2, Y=1.0, dzf=2e-4, dzc=4e-3, zfine=0.06, L=10, kappa=10, log=print)
Yv = np.zeros((6, run.n)); R = 0.0
t0 = time.time()
for k in range(40): run.rhs(0.1*k, Yv, R)
print('rhs per call %.2f ms, n=%d' % ((time.time()-t0)/40*1e3, run.n))
t0 = time.time()
for k in range(40): run.ref(0.1*k)
print('ref per call %.2f ms' % ((time.time()-t0)/40*1e3))
cProfile.run('for k in range(30): run.rhs(0.1*k, Yv, R)', 'diag/prof.out')
pstats.Stats('diag/prof.out').sort_stats('tottime').print_stats(12)
