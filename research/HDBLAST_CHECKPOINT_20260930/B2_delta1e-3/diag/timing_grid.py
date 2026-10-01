"""B2 dev: grid sizes and rhs cost for the rho-adapted grid at delta = 1e-3 (tuned d*), identity and bounded charts."""
import sys, time, math, json
sys.path.insert(0, '.'); sys.path.insert(0, '/home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927/mechanisms/pilot_5d')
import numpy as np, evolve_a1 as E
t0 = time.time()
r0 = E.Run(delta=1e-3, d=-3.1942416958680835, dc=1e-2, Y=1.0, dzf=2e-4, dzc=2e-2, L=20, kappa=10, log=print, grid_kind='rho')
print('shell setup %.0fs' % (time.time() - t0))
out = []
for dzf, dzc, L, xc in [(4e-4, 4e-2, 24, math.inf), (2.5e-4, 2e-2, 20, math.inf), (1.7e-4, 2e-2, 20, math.inf), (1.25e-4, 1.5e-2, 20, math.inf), (2.5e-4, 2e-2, 20, 20.0)]:
    r = E.Run(delta=1e-3, d=-3.1942416958680835, dc=1e-2, Y=1.0, dzf=dzf, dzc=dzc, L=L, kappa=10, log=lambda s: None, grid_kind='rho',
              shells=(r0.bg, r0.seed), xc=xc)
    Yv = np.zeros((6, r.n)); R = 0.0
    r.rhs(0.0, Yv, R)
    t1 = time.time(); k = 0
    while time.time() - t1 < 3: r.rhs(1e-3*k, Yv, R); k += 1
    ms = (time.time() - t1)/k*1e3
    z = r.z; zp = r.zp
    row = dict(dzf=dzf, dzc=dzc, L=L, xc=xc, n=r.n, rhs_ms=ms, dt=0.5*zp.min(), pts_per_wall=(1/r.rb)/zp[-1],
               dz_at={str(zz): float(np.interp(zz, z, zp)) for zz in [-0.05, -0.1, -0.2, -0.5, -1, -2, -5]},
               est_min_per_T=4*ms/1e3/(0.5*zp.min())/60)
    print(json.dumps(row)); out.append(row)
json.dump(out, open('diag/TIMING_GRID.json', 'w'), indent=1)
