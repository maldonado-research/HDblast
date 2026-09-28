#!/usr/bin/env python3
"""Slow-roll survey of natural detuning families on the HJ shell-modulus field space (derivatives taken in y_b)."""
import numpy as np, json
d = np.load("eft_landscape.npz"); ys, phi, f, ZEy = d["y"], d["phi"], d["f"], d["ZEy"]
m = (ys <= 9) & (ys >= -16); ys, phi, f, ZEy = ys[m], phi[m], f[m], ZEy[m]
ys, phi, f, ZEy = ys[::20], phi[::20], f[::20], ZEy[::20]          # dy = 2e-3
Ty = -np.sqrt(ZEy)                                                 # dTheta/dy (Theta increases toward phi=+1, i.e. decreasing y)
W  = lambda p: 1 - p + p**3/3
c_star = 2/1.0357712571566782 - 4/3
fam = {
 "registered: 1 + c* phi":        lambda p: 1 + c_star*p,
 "(1-phi)":                       lambda p: 1 - p,
 "(1-phi)^2":                     lambda p: (1 - p)**2,
 "W - W(+1) = (1-phi)^2(2+phi)/3": lambda p: W(p) - 1/3,
 "W_phi^2 = (1-phi^2)^2":          lambda p: (1 - p*p)**2,
 "(1-phi^2)":                     lambda p: 1 - p*p,
 "(1+phi)":                       lambda p: 1 + p,
 "(1+phi)^2":                     lambda p: (1 + p)**2,
 "W(-1) - W = (1+phi)^2(2-phi)/3": lambda p: 5/3 - W(p),
}
out = {}
for name, Vf in fam.items():
    V = Vf(phi)/f**2
    Vy = np.gradient(V, ys); Vth = Vy/Ty; Vthth = np.gradient(Vth, ys)/Ty
    with np.errstate(divide='ignore', invalid='ignore'):
        eps = 0.5*(Vth/V)**2; eta = Vthth/V
    good = (V > 1e-12) & np.isfinite(eps) & np.isfinite(eta); good[:5] = False; good[-5:] = False
    print("\n== V ∝ %s" % name); rec = dict(stationary=[], windows=[])
    sg = np.sign(Vy); idx = np.where((sg[:-1]*sg[1:] < 0) & good[:-1])[0]
    for i in idx:
        print("   stationary: y_b=%+.3f phi_b=%+.5f  V_E/t=%.5f  eta_V=%+.4f  mu2=%+.4f" % (ys[i], phi[i], V[i], eta[i], 3*eta[i]))
        rec["stationary"].append(dict(y_b=float(ys[i]), phi_b=float(phi[i]), mu2=float(3*eta[i])))
    sr = good & (eps < 0.1) & (np.abs(eta) < 0.05)
    if sr.any():
        ii = np.where(sr)[0]; br = np.where(np.diff(ii) > 1)[0]; st = np.r_[ii[0], ii[br+1]]; en = np.r_[ii[br], ii[-1]]
        for a, b in zip(st, en):
            if b - a < 20: continue
            seg = slice(a, b+1)
            N = abs(np.trapz(np.abs(Ty[seg])/np.sqrt(2*eps[seg] + 1e-300), ys[seg]))
            print("   slow-roll window: y_b in [%+.2f,%+.2f], phi_b in [%+.5f,%+.5f], e-folds available ~ %.1f, eta in [%+.4f,%+.4f], eps max %.2e" % (
                ys[b], ys[a], phi[b], phi[a], N, eta[seg].min(), eta[seg].max(), eps[seg].max()))
            rec["windows"].append(dict(y_lo=float(ys[b]), y_hi=float(ys[a]), N=float(N), eta_min=float(eta[seg].min()), eta_max=float(eta[seg].max())))
    else:
        print("   no window with eps<0.1 and |eta|<0.05")
    out[name] = rec
json.dump(out, open("EFT_FAMILIES_SURVEY.json", "w"), indent=1)
