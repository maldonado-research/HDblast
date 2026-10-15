#!/usr/bin/env python3
"""Fresh rational first-exit/majorant feasibility calculation at ell=1e-6.

This program does not execute or import a historical producer.  It uses only
the retained exact shell-box JSON and the equations stated in the pinned note.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
INPUTS = {
    "inputs/previous/archive/ANALYTIC_SHELL_BOX_DELTA_1E3.json": "8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8",
    "inputs/previous/archive/ACTUAL_BACKGROUND_AND_RIGOROUS_COLLAR_BOUNDS.md": "0832a16c4bdf18cc188861a2c2860e7e6ff736dec39426d6ca7c589dc0434d5c",
    "sol/UNIFORM_WAVEFORM_CERTIFICATE.md": "a9b0354c2ecb8a68d04e9441259de7873c3c93fb7cb4369257478b31bfffcd56",
}


class I:
    def __init__(self, lo, hi=None):
        self.lo = F(lo)
        self.hi = self.lo if hi is None else F(hi)
        if self.lo > self.hi:
            raise ValueError("reversed rational interval")

    def __add__(self, other):
        other = interval(other)
        return I(self.lo+other.lo, self.hi+other.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self+-interval(other)

    def __rsub__(self, other):
        return interval(other)+-self

    def __mul__(self, other):
        other = interval(other)
        ends = [a*b for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
        return I(min(ends), max(ends))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = interval(other)
        if other.lo <= 0 <= other.hi:
            raise ValueError("division by interval containing zero")
        return self*I(1/other.hi, 1/other.lo)

    def __rtruediv__(self, other):
        return interval(other)/self

    def square(self):
        return I(0 if self.lo <= 0 <= self.hi else min(self.lo**2, self.hi**2),
                 max(self.lo**2, self.hi**2))

    def absolute_upper(self):
        return max(abs(self.lo), abs(self.hi))

    def record(self):
        return [str(self.lo), str(self.hi)]


def interval(value):
    return value if isinstance(value, I) else I(value)


def polynomial(coefficients, x):
    result = I(0)
    for coefficient in reversed(coefficients):
        result = result*x+coefficient
    return result


def product(a, b):
    out = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def derivative(a):
    return [i*a[i] for i in range(1, len(a))]


def majorant(ell, m, M, B, Gamma):
    T = 3*ell
    a, r, s = 1+T, M*T*T/2, B*T*T/(2*m)
    determinant = (1-r)*(1-s)-T*r/m
    if not (m > 0 and 1-r > 0 and 1-s > 0 and determinant > 0):
        raise ValueError("majorant matrix positivity failed")
    A = a*(1-s+T/m)/determinant
    V = a/(m*determinant)
    if not ((1-r)*A-T*V == a and (1-s)*V-r*A/m == a/m):
        raise ValueError("coupled majorant equations failed")
    eta = ell*(M*T*A+(1+Gamma*T)*V)
    Vtail = ell/m**2+ell/m+r*A/m+s*V
    tail = ell*(M*T*A+Vtail+Gamma*T*V)
    return {"ell": ell, "T": T, "m": m, "M": M, "B": B, "Gamma": Gamma,
            "r": r, "s": s, "determinant": determinant, "A_upper": A,
            "V_upper": V, "eta_upper": eta, "V_tail_upper": Vtail,
            "eta_tail_upper": tail, "L2_error_squared_upper": eta**2+tail**2}


def run(source_root=ROOT):
    checks = []

    def check(name, predicate):
        if not predicate:
            raise RuntimeError(name)
        checks.append(name)

    pins = {}
    for path, expected in INPUTS.items():
        data = (source_root/path).read_bytes()
        got = sha256(data).hexdigest()
        check("authenticated mathematical input "+path, got == expected)
        pins[path] = {"sha256": got, "bytes": len(data)}
    source = json.loads((source_root/next(iter(INPUTS))).read_text())
    check("same exact detuning", source["parameters_exact"]["delta"] == "1/1000")
    check("same exact branch parameter", source["parameters_exact"]["c"] == "96714321425692/161839258930731")
    outer = {"rho": I(129, 131), "eta": I(F(-9, 10**5), F(-8, 10**5)),
             "p": I(F(-14, 10**5), F(-12, 10**5)), "a": I(F(111, 1000), F(112, 1000))}
    initial = {name: I(*source["boxes"][name+"_b"]["exact"])
               for name in ("rho", "eta", "p", "a")}
    rho, eta, p, a = (outer[name] for name in ("rho", "eta", "p", "a"))
    # U=(W_eta)^2/2 - 2 W^2/3, W=1/3+eta^2+eta^3/3.
    superpotential = [F(1, 3), F(0), F(1), F(1, 3)]
    superprime = derivative(superpotential)
    halfprime = [v/2 for v in product(superprime, superprime)]
    wsquared = [2*v/3 for v in product(superpotential, superpotential)]
    U = [(halfprime[i] if i < len(halfprime) else 0)-wsquared[i]
         for i in range(len(wsquared))]
    Up, Upp = derivative(U), derivative(derivative(U))
    check("independent exact potential polynomial", U == [F(-2,27),0,F(14,9),F(50,27),F(-1,6),F(-4,9),F(-2,27)])
    check("independent potential derivative", Up == [0,F(28,9),F(50,9),F(-2,3),F(-20,9),F(-4,9)])
    up, upp = polynomial(Up, eta), polynomial(Upp, eta)
    rhs = {"rho": -rho.square()*a, "eta": -rho*p,
           "p": -rho*(up-4*a*p), "a": 1/rho+rho*p.square()/3}
    ell = F(1, 10**6)
    L = 8*ell
    first_exit = {}
    depth_thresholds = []
    for name in outer:
        shell, box, flow = initial[name], outer[name], rhs[name]
        check("initial strictly inside outer "+name, box.lo < shell.lo <= shell.hi < box.hi)
        image = shell+I(0,L)*flow
        lower, upper = image.lo-box.lo, box.hi-image.hi
        check("strict full-prefix first-exit lower "+name, lower > 0)
        check("strict full-prefix first-exit upper "+name, upper > 0)
        limits = {}
        if flow.lo < 0:
            limits["lower"] = (shell.lo-box.lo)/(-flow.lo)
        if flow.hi > 0:
            limits["upper"] = (box.hi-shell.hi)/flow.hi
        depth_thresholds.extend((value, name, side) for side, value in limits.items())
        first_exit[name] = {"initial": shell.record(), "outer": box.record(), "rhs": flow.record(),
                            "integrated_full_prefix_image": image.record(), "lower_margin": str(lower),
                            "upper_margin": str(upper), "depth_thresholds": {k:str(v) for k,v in limits.items()}}
    limiting_depth, limiting_name, limiting_side = min(depth_thresholds)
    check("target depth below limiting rectangle threshold", L < limiting_depth)
    # This wrong control asks too much of this same rectangle; failure to
    # certify depth 5e-5 does not assert nonexistence of the actual background.
    too_deep = initial["eta"]+I(0,F(1,20000))*rhs["eta"]
    check("overextended same rectangle rejected", too_deep.hi >= outer["eta"].hi)
    check("smooth continuation denominators bounded away from zero", rho.lo > 0 and p.hi < 0)

    h = up/p
    C = rho*(2*h-5*a)
    potential = I(F(-35,4), 0)+rho.square()*(2*h.square()-10*a*h+F(35,4)*a.square()-upp-p.square()/6)
    mu = I(*source["boxes"]["mu"]["exact"])
    gamma = I(*source["boxes"]["C_b"]["exact"])/2
    beta = mu*I(F(-25,4),F(5,2))+gamma
    check("retained coarse mu on widened collar", mu.lo > F(1,1000))
    check("retained coarse W on widened collar", potential.absolute_upper() < 10**6)
    check("retained coarse beta on widened collar", beta.absolute_upper() < 1000)
    check("retained coarse gamma on widened collar", gamma.absolute_upper() < 250)
    check("stronger rounded mu", mu.lo > F(1,500))
    check("stronger rounded W", potential.absolute_upper() < 102000)
    check("stronger rounded beta", beta.absolute_upper() < 250)
    check("stronger rounded gamma", gamma.absolute_upper() < 250)
    coarse = majorant(ell,F(1,1000),F(10**6),F(1000),F(250))
    tighter = majorant(ell,F(1,500),F(102000),F(250),F(250))
    for label, values, eta_bound, tail_bound, norm_bound in (
        ("retained_coarse",coarse,F(100378,10**8),F(477,10**8),F(10038,10**7)),
        ("stronger_rounded",tighter,F(500684,10**9),F(9325,10**10),F(500685,10**9))):
        check(label+" strict global error", values["eta_upper"] < eta_bound)
        check(label+" strict tail error", values["eta_tail_upper"] < tail_bound)
        check(label+" strict support-aware L2 error", values["L2_error_squared_upper"] < norm_bound**2)
        values.update(global_display_upper=eta_bound,tail_display_upper=tail_bound,L2_display_upper=norm_bound)
        # Check that the inverse genuinely retains the feedback coupling.
        wrong_det = (1-values["r"])*(1-values["s"])
        wrong_V = (1+values["T"])/(values["m"]*wrong_det)
        wrong_A = (1+values["T"])*(1-values["s"]+values["T"]/values["m"])/wrong_det
        check(label+" omitted determinant coupling rejected", (1-values["r"])*wrong_A-values["T"]*wrong_V < 1+values["T"])

    def encoded(value):
        if isinstance(value,F):
            return str(value)
        if isinstance(value,dict):
            return {k:encoded(v) for k,v in value.items()}
        if isinstance(value,(tuple,list)):
            return [encoded(v) for v in value]
        return value

    return encoded({
        "status": "PASS_BOUNDED_WIDER_STATE_FEASIBILITY",
        "arithmetic": "Exact integer/Fraction interval arithmetic; no sampled trajectory or imported producer",
        "python_optimization": sys.flags.optimize,
        "input_pins": pins,
        "input_resolution": "Three hash-pinned paths relative to --source-root; no producer import or execution",
        "target": {"ell":ell,"observation_T":3*ell,"collar_depth":L,"depth_enlargement_factor":10000},
        "first_exit": first_exit,
        "same_rectangle_strict_depth_threshold": {"depth":limiting_depth,"limiting_coordinate":limiting_name,
                                                    "limiting_side":limiting_side,"ell_if_L_equals_8ell":limiting_depth/8,
                                                    "interpretation":"Sufficient strict first-exit threshold of this fixed rectangle, not maximal actual existence depth"},
        "independently_recomputed_coefficients": {"U_coefficients":U,"U_prime_coefficients":Up,"U_second_coefficients":Upp,
                    "U_prime_outer":up.record(),"U_second_outer":upp.record(),"h_outer":h.record(),"C_outer":C.record(),
                    "W_all_real_t_outer":potential.record(),"mu_shell":mu.record(),"gamma_shell":gamma.record(),"beta_all_real_t_shell":beta.record()},
        "uniform_majorants": {"retained_coarse":coarse,"stronger_rounded":tighter},
        "checks":checks,"check_count":len(checks),
        "scope":"Bounded feasibility deduction for a next wider reconstruction, conditional on retained analytic actual-background and coupled linear evolution hypotheses. No wider waveform export, response C interval, effective nonlinear amplitude, observational scale or publication claim.",
        "reference_only_producer_executed":False,
        "own_checker_sha256":sha256(Path(__file__).read_bytes()).hexdigest(),
    })


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output",default="WIDER_FEASIBILITY_NORMAL.json")
    parser.add_argument("--source-root",type=Path,default=ROOT,
                        help="Extracted unchanged waveform capsule containing the three pinned input paths")
    args = parser.parse_args()
    if Path(args.output).name != args.output:
        raise ValueError("output must be a filename inside this owned planning directory")
    result = run(args.source_root.resolve())
    (HERE/args.output).write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("PASS:",result["check_count"],"exact wider-state feasibility checks")
