#!/usr/bin/env python3
"""HDBLAST Chat 9 - replay of the orchestrator's core calculations (python3 + numpy only; about 3-6 minutes).
Usage:  python3 REPLAY.py          Writes REPLAY_RESULT.json and exits non-zero on any failed gate.
Stages: exact identities (fractions) -> background vs certified M462 root -> 5D tachyon eigenvalue -> HJ closed form
        -> instanton entropy identity.  Floating-point stages are diagnostics with explicit tolerances, not interval proofs."""
import json, math, os, subprocess, sys
here = os.path.dirname(os.path.abspath(__file__)); out = {}; ok = True
def gate(name, cond, info):
    global ok
    out[name] = dict(passed=bool(cond), **info); ok = ok and bool(cond)
    print("[%s] %s  %s" % ("PASS" if cond else "FAIL", name, json.dumps(info)), flush=True)

# 1. exact identities
r = subprocess.run([sys.executable, "verify_exact_identities.py"], cwd=os.path.join(here, "exact_checks"), capture_output=True, text=True)
gate("exact_identities", r.returncode == 0 and "ALL PASS" in r.stdout, dict(tail=r.stdout.strip().splitlines()[-1] if r.stdout else r.stderr[-200:]))

# 2. background versus the certified M462 root
sys.path.insert(0, os.path.join(here, "background")); import hdblast_background as bg
Ip = bg.I_plus(); c = 2/Ip - 4/3; t = 1e-3
sol = bg.solve_shell(t, c, guess=(math.log10(8.7855e-7), 8.2282), dv=4e-4)
delta = sol["alpha"]/t**1.8; eta = (sol["h"]/(t/(6*Ip)) - 1)/t
gate("background_reproduces_M462_root", abs(delta - 0.2206821504276063) < 5e-8 and abs(eta - 0.12006349132995513) < 5e-6,
     dict(delta=delta, eta=eta, certified_delta=0.2206821504276063, certified_eta=0.12006349132995513, h=sol["h"], rho_b=sol["rho_b"]))

# 3. 5D scalar-sector eigenvalue (vectorised shooting)
sys.path.insert(0, os.path.join(here, "stability", "orchestrator_derivation"))
import numpy as np
from t_scan_orchestrator import shoot_vec, roots_from_grid
f = lambda m: shoot_vec(m, sol["phi_h"], sol["v_b"], t, c, n=6000)
grid = np.linspace(-14, 2.2, 82); roots = roots_from_grid(grid, f(grid), f)
gate("exactly_one_bound_state_and_it_is_tachyonic", len(roots) == 1 and abs(roots[0] + 7.7178716) < 2e-6, dict(roots=[float(x) for x in roots], expected=-7.7178716))

# 4. Hamilton-Jacobi closed form and O(t) consistency
mu0 = -4*(3*c*c - 4*c + 8)/(c*(3*c + 4))
gate("HJ_closed_form_matches_5D_to_order_t", abs((roots[0] - mu0)/t - 1.9243) < 2e-3 if roots else False, dict(mu0_closed_form=mu0, slope_estimate=(roots[0] - mu0)/t if roots else None, expected_slope=1.9243))

# 5. instanton entropy identity (exact identity; the direct 5D evaluation cancels two terms of order 1e11,
#    so the floating ratio differs from 1 by ~1e-5 at this step size: tolerance 1e-4)
yend, rec = bg.integrate(sol["phi_h"], sol["v_b"], dv=4e-4, record=True)
y, rho, phi = rec[:, 0], rec[:, 1], rec[:, 2]
trap = lambda a, x: float(np.sum(0.5*(a[1:] + a[:-1])*np.diff(x)))
Om4 = 8*math.pi**2/3; sig = 2*bg.W(sol["phi_b"]) + t*(1 + c*sol["phi_b"])
SE = -(4/3)*Om4*(trap(rho**4*bg.U(phi), y) + rec[0, 0]**5/5*bg.U(sol["phi_h"])) - Om4/3*sol["rho_b"]**4*sig
SdS = 16*math.pi**2*(trap(rho**2, y) + rec[0, 0]**3/3)
gate("instanton_action_equals_minus_dS_entropy", abs(SE/(-SdS) - 1) < 1e-4, dict(S_E_5D=SE, minus_16pi2_int_rho2=-SdS, ratio=SE/(-SdS)))

json.dump(dict(all_pass=ok, stages=out), open(os.path.join(here, "REPLAY_RESULT.json"), "w"), indent=1)
print("REPLAY:", "ALL PASS" if ok else "FAILED"); sys.exit(0 if ok else 1)
