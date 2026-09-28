import json, math, numpy as np, time
from gn_lib import *
t0=time.time()
Ip = I_plus(); c = 2/Ip - 4/3; t = 1e-3
sol, bg = polish_shell(t, c, 8.785514694954969e-07, 8.228201550199884, n=4000)
print(sol); C = bg_coeffs(bg)
mu2 = np.concatenate([np.linspace(-400, -30, 371)[:-1], np.linspace(-30, 2.2499, 3226)])
X, rec = integrate_pert(mu2, C, sol['phi_h'], record_every=500)
M, chis, eps = mismatch(X, C, sol)
Es, psis, chis2, eps2 = metric_mismatch(X, C, sol)
np.save('scan.npy', np.array([mu2, M, chis, eps, Es, psis, chis2]))
for name, f in [('M (scalar-junction mismatch, eps by J1)', M), ('chi_s', chis), ('E_s prime (eps by J2)', Es), ('psi_s (eps by J2)', psis)]:
    sg = np.sign(f); idx = np.where(sg[1:]*sg[:-1] < 0)[0]
    print(name, 'sign changes at mu2 ~', [(round(mu2[i],4), 'pole?' if abs(f[i])+abs(f[i+1]) > 1e3*np.median(abs(f)) else '') for i in idx])
# constraint monitors
for (i, Xi) in rec:
    r = residuals(i, Xi, mu2, C)
    print('y=%.3f  max|HC|=%.2e  max|EVpsi|=%.2e  max|YY|=%.2e' % (C['y'][i], *[np.max(np.abs(v)) for v in r]))
print('time', time.time()-t0)
