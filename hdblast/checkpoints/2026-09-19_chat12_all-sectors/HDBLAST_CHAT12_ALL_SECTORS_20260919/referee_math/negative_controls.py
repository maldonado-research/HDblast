#!/usr/bin/env python3
"""Referee (math) - negative controls for the three Chat 12 scripts by TEXT MUTATION of copies.
Each mutant must (a) apply (the search string occurs), (b) make the named check(s) FAIL, exit code != 0.
Originals are never touched; mutants run in ./negctl/<name>/ .   Run with the venv python."""
import json, os, shutil, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); SRC = os.path.dirname(HERE)
PY = sys.executable
M = [  # (tag, script, old, new, expected failing check prefix)
 ("T_4H_to_3H", "verify_tensor_sector.py", "4*H*sp.diff(h, y)", "3*H*sp.diff(h, y)", "T1"),
 ("T_m2_to_2m2", "verify_tensor_sector.py", "+ m2*h/rho**2)", "+ 2*m2*h/rho**2)", "T1"),
 ("T_9over16_to_1over2", "verify_tensor_sector.py", "9*pp**2/16", "pp**2/2", "T3"),
 ("T_Uover8_to_Uover6", "verify_tensor_sector.py", "pp**2/2)/8))", "pp**2/2)/6))", "T3"),
 ("T_Q_coeff_3half_to_1", "verify_tensor_sector.py", "Q = lambda f: dz(f) - sp.Rational(3, 2)*Hc*f", "Q = lambda f: dz(f) - sp.Rational(1, 1)*Hc*f", "T2"),
 ("T_u_weight_3half_to_2", "verify_tensor_sector.py", "hh = u/rho**sp.Rational(3, 2)", "hh = u/rho**2", "T2"),
 ("T_Urange_3over20_to_1over5", "verify_tensor_sector.py", "PU.count_roots(-1, sp.Rational(3, 20))", "PU.count_roots(-1, sp.Rational(1, 5))", "T4a"),
 ("T_Ubound_range_[-1,0]_to_[-1,1/10]", "verify_tensor_sector.py", "inner.count_roots(-1, 0) == 0 and inner.eval(-sp.Rational(1, 2)) > 0", "inner.count_roots(-1, sp.Rational(1,10)) == 0 and (U + sp.Rational(1, 6)).subs(p, sp.Rational(1, 10)) <= 0", "T4b"),
 ("V_2H_to_3H", "verify_vector_sector.py", "C1 = sp.diff(Bv, y) + 2*H*Bv", "C1 = sp.diff(Bv, y) + 3*H*Bv", "V2a"),
 ("V_Box+3_to_Box+2", "verify_vector_sector.py", "Bv/(2*rho**2)*op[2]) == 0", "Bv/(2*rho**2)*(op[2] - Vlow[2])) == 0", "V2b"),
 ("V_Box+3_to_Box+4", "verify_vector_sector.py", "Bv/(2*rho**2)*op[2]) == 0", "Bv/(2*rho**2)*(op[2] + Vlow[2])) == 0", "V2b"),
 ("V_gauge_rho2_to_rho", "verify_vector_sector.py", "Lie[0, m+1] - rho**2*sp.diff(zeta, y)*Vlow[m]", "Lie[0, m+1] - rho*sp.diff(zeta, y)*Vlow[m]", "V1"),
 ("V_nontransverse_harmonic", "verify_vector_sector.py", "Vlow = [0, 0, T*E, 0]", "Vlow = [0, T*E, 0, 0]", "V0"),
 ("S_Hessian_sign", "verify_special_harmonics.py", "sp.simplify(Hs[m, n_] + gam[m, n_]*Yf) == 0", "sp.simplify(Hs[m, n_] - gam[m, n_]*Yf) == 0", "S1"),
 ("S_wrong_harmonic_e^-tau", "verify_special_harmonics.py", '"e^tau": sp.exp(tau),', '"e^tau": sp.exp(-tau),', "S1"),
 ("S_mu2_-4_to_-3", "verify_special_harmonics.py", "\nmu2 = -4\n", "\nmu2 = -3\n", "S2c"),
 ("S_Codazzi_3_to_2", "verify_special_harmonics.py", "red(chi_g + 3*(sp.diff(psi_g, y)", "red(chi_g + 2*(sp.diff(psi_g, y)", "S2d"),
 ("S_gauge_ODE_4_to_3", "verify_special_harmonics.py", "b2 - (-(4*rho*sp.diff(rho, y)*sp.diff(b, y) + 2*b)/rho**2)", "b2 - (-(3*rho*sp.diff(rho, y)*sp.diff(b, y) + 2*b)/rho**2)", "S2b"),
 ("S_B_sigma2_half_to_full", "verify_special_harmonics.py", "Bq = sp.diff(ph0, y, 2)/sp.diff(ph0, y) + s2/2", "Bq = sp.diff(ph0, y, 2)/sp.diff(ph0, y) + s2", "S3b"),
 ("S_B_plus_one", "verify_special_harmonics.py", "Bq = sp.diff(ph0, y, 2)/sp.diff(ph0, y) + s2/2", "Bq = sp.diff(ph0, y, 2)/sp.diff(ph0, y) + s2/2 + 1", "S3b"),
 ("S_bg_junction_-2phi'_to_-phi'", "verify_special_harmonics.py", "e = e.subs(s0, 6*sp.diff(rho, y)/rho).subs(s1, -2*sp.diff(ph0, y))", "e = e.subs(s0, 6*sp.diff(rho, y)/rho).subs(s1, -1*sp.diff(ph0, y))", "S3a"),
 ("S_bg_junction_6H_to_4H", "verify_special_harmonics.py", "e = e.subs(s0, 6*sp.diff(rho, y)/rho)", "e = e.subs(s0, 4*sp.diff(rho, y)/rho)", "S3a"),
]
out = {}; allok = True
root = os.path.join(HERE, "negctl"); shutil.rmtree(root, ignore_errors=True); os.makedirs(root)
for tag, script, old, new, expect in M:
    d = os.path.join(root, tag.replace("/", "_").replace("'", "")); os.makedirs(d)
    shutil.copy(os.path.join(SRC, "lin_gr.py"), d)
    src = open(os.path.join(SRC, script)).read()
    applied = src.count(old) >= 1
    open(os.path.join(d, script), "w").write(src.replace(old, new, 1))
    p = subprocess.run([PY, script], cwd=d, capture_output=True, text=True, timeout=900)
    fails = [l.split()[1] for l in p.stdout.splitlines() if l.startswith("[FAIL]")]
    rejected = applied and p.returncode != 0 and any(f.startswith(expect) for f in fails)
    crashed = ("Traceback" in p.stderr)
    out[tag] = dict(applied=applied, returncode=p.returncode, failed_checks=[f[:40] for f in fails], expected=expect, rejected=rejected, crashed=crashed)
    allok &= rejected and not crashed
    print("%-34s applied=%s rc=%d rejected=%s fails=%s%s" % (tag, applied, p.returncode, rejected, [f[:12] for f in fails], "  CRASH" if crashed else ""), flush=True)
json.dump(dict(all_mutants_rejected=allok, mutants=out), open(os.path.join(HERE, "NEGATIVE_CONTROLS_RESULT.json"), "w"), indent=1)
print("ALL MUTANTS REJECTED" if allok else "SOME MUTANT SURVIVED")
