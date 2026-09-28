#!/usr/bin/env python3
"""Negative controls for verify_linearised_equations.py: each deliberately wrong coefficient must make checks FAIL.
Run:  ../../.hdblast_venv/bin/python negative_controls.py"""
import json, os, subprocess, sys, tempfile
src = open("verify_linearised_equations.py").read()
mutations = {
 "M1_master_(2+mu2)->(3+mu2)":        [("(2 + mu2)/rho**2", "(3 + mu2)/rho**2")],
 "M2_master_4/3->5/3":                [("-sp.Rational(4, 3)*pp**2", "-sp.Rational(5, 3)*pp**2"), ("-sp.Rational(4, 3)*sp.diff(ph0, y)**2", "-sp.Rational(5, 3)*sp.diff(ph0, y)**2")],
 "M3_Codazzi_chi_-3->-2":             [("chi_expr = -3*(", "chi_expr = -2*(")],
 "M4_xi=-2psi->xi=-psi":              [("e = e.subs(xi, -2*psi)", "e = e.subs(xi, -psi)")],
 "M5_master_4H->3H_in_zeroth_term":   [("- 4*H*g_ +", "- 3*H*g_ +"), ("- 4*sp.diff(rho, y)/rho*(U1 - 4*sp.diff(rho, y)/rho*sp.diff(ph0, y))/sp.diff(ph0, y) +", "- 3*sp.diff(rho, y)/rho*(U1 - 4*sp.diff(rho, y)/rho*sp.diff(ph0, y))/sp.diff(ph0, y) +")],
 "M6_shell_lam=mu2+4->mu2+3":         [("3*(mu2 + 4)*psi/(rho**2*pp)", "3*(mu2 + 3)*psi/(rho**2*pp)")],
 "M7_dS_curvature_wrong_(flat_slices_e^{2tau}->e^{3tau})": [("gam = [-1, sp.exp(2*tau), sp.exp(2*tau), sp.exp(2*tau)]", "gam = [-1, sp.exp(3*tau), sp.exp(3*tau), sp.exp(3*tau)]")],
}
out = {}
for name, reps in mutations.items():
    s = src
    for a, b in reps:
        assert a in s, (name, a); s = s.replace(a, b)
    s = s.replace('open("SYMBOLIC_VERIFICATION_RESULT.json", "w")', 'open(os.devnull, "w")').replace("import json, sys, time", "import json, os, sys, time")
    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as f: f.write(s); path = f.name
    r = subprocess.run([sys.executable, path], capture_output=True, text=True); os.unlink(path)
    fails = [l for l in r.stdout.splitlines() if l.startswith("[FAIL]")]
    out[name] = dict(exit_code=r.returncode, n_fail=len(fails), first_fail=(fails[0][:90] if fails else None))
    print("%-58s exit=%d  failed checks=%d" % (name, r.returncode, len(fails)), flush=True)
ok = all(v["exit_code"] != 0 and v["n_fail"] > 0 for v in out.values())
json.dump(dict(all_mutations_rejected=ok, mutations=out), open("NEGATIVE_CONTROLS_RESULT.json", "w"), indent=1)
print("ALL MUTATIONS REJECTED" if ok else "SOME MUTATION WAS NOT DETECTED"); sys.exit(0 if ok else 1)
