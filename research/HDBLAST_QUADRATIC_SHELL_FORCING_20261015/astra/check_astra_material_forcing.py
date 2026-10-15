#!/usr/bin/env python3
"""Fresh symbolic checks of material constitutive shell sources.

No prior scientific implementation is imported or executed. Explicit failures
remain active under python -O. --inject changes a candidate mathematical
formula before comparison with separately differentiated exact expressions.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import sympy as s


INJECTIONS = (
    "none", "reverse_one_normal", "double_tension", "half_hessian",
    "wrong_projection", "omit_unknown_response", "omit_clock_transport",
    "drop_metric_stress",
)


class FormulaFailure(Exception):
    pass


def run(inject):
    checks = []

    def equal(name, lhs, rhs):
        residue = s.factor(s.expand(lhs - rhs))
        if residue != 0:
            raise FormulaFailure(f"{name}: nonzero symbolic residue {residue}")
        checks.append(name)

    e, z = s.symbols("epsilon eta")
    eta, delta, c = s.symbols("eta_b delta c")
    psi1, psi2 = s.symbols("psi1 psi2")
    pp1, pm1, pp2, pm2 = s.symbols("p_plus1 p_minus1 p_plus2 p_minus2")
    sigma = 2 * (s.Rational(1, 3) + z**2 + z**3 / 3) + delta*(1+c+c*z)
    sig = [s.diff(sigma, z, k).subs(z, eta) for k in range(4)]
    equal("registered sigma second derivative", sig[2], 4+4*eta)
    equal("registered sigma third derivative", sig[3], 4)
    eta_e = eta + e*psi1 + e**2*psi2/2
    sigma_e = sigma.subs(z, eta_e)
    pp_e = -sig[1]/2 + e*pp1 + e**2*pp2/2
    pm_e = -sig[1]/2 + e*pm1 + e**2*pm2/2
    residual = pp_e + pm_e + s.diff(sigma, z).subs(z, eta_e)
    equal("outward sum background shell law", residual.subs(e, 0), 0)
    equal("first scalar sum row", s.diff(residual, e).subs(e, 0), pp1+pm1+sig[2]*psi1)
    scalar2 = pp2+pm2+sig[2]*psi2+sig[3]*psi1**2
    if inject == "reverse_one_normal":
        scalar2 = pp2-pm2+sig[2]*psi2+sig[3]*psi1**2
    elif inject == "double_tension":
        scalar2 = pp2+pm2+2*sig[2]*psi2+2*sig[3]*psi1**2
    elif inject == "half_hessian":
        scalar2 = pp2+pm2+sig[2]*psi2+sig[3]*psi1**2/2
    elif inject == "omit_unknown_response":
        scalar2 = pp2+pm2+sig[3]*psi1**2
    equal("second scalar sum derivative", s.diff(residual, e, 2).subs(e, 0), scalar2)
    equal("factorial coefficient is half the derivative", s.expand(residual).coeff(e, 2), scalar2/2)
    p2 = s.symbols("p2")
    equal("identical copy division", scalar2.subs({pp2:p2, pm2:p2})/2,
          p2+sig[2]*psi2/2+2*psi1**2)
    equal("tension second derivative", s.diff(sigma_e, e, 2).subs(e, 0),
          sig[1]*psi2+sig[2]*psi1**2)

    h0, h1, h2 = s.symbols("h0 h1 h2", nonzero=True)
    he = h0+e*h1+e**2*h2/2
    stress = -sigma_e*he
    stress0 = -sig[0]*h0
    stress1 = -sig[0]*h1-sig[1]*psi1*h0
    stress2 = -sig[0]*h2-2*sig[1]*psi1*h1-(sig[1]*psi2+sig[2]*psi1**2)*h0
    if inject == "drop_metric_stress":
        stress2 += 2*sig[1]*psi1*h1
    equal("lower-index stress product", s.diff(stress,e,2).subs(e,0), stress2)
    inverse0, inverse1, inverse2 = 1/h0, -h1/h0**2, 2*h1**2/h0**3-h2/h0**2
    equal("inverse-metric cancellation in mixed stress",
          inverse0*stress2+2*inverse1*stress1+inverse2*stress0,
          -sig[1]*psi2-sig[2]*psi1**2)
    equal("four-dimensional mixed stress trace", s.diff(-4*sigma_e,e,2).subs(e,0),
          -4*sig[1]*psi2-4*sig[2]*psi1**2)
    israel2 = s.diff(sigma_e*he/6,e,2).subs(e,0)
    equal("material Israel product", israel2,
          sig[0]*h2/6+sig[1]*psi2*h0/6+sig[1]*psi1*h1/3+sig[2]*psi1**2*h0/6)

    vel, dpsi1, dpplus1 = s.symbols("clock_velocity dpsi1 dpplus1")
    dpminus1 = -sig[2]*dpsi1-dpplus1
    psi2hat = psi2+2*vel*dpsi1
    if inject == "omit_clock_transport":
        psi2hat = psi2
    transported = pp2+2*vel*dpplus1+pm2+2*vel*dpminus1+sig[2]*psi2hat
    equal("clock covariance using differentiated first junction", transported,
          pp2+pm2+sig[2]*psi2)
    equal("tension response clock transport",
          sig[1]*psi2hat+sig[2]*psi1**2-(sig[1]*psi2+sig[2]*psi1**2),
          2*vel*sig[1]*dpsi1)

    x = s.symbols("X1:5")
    r2 = sum(t*t for t in x)
    Y = x[0]**2-x[1]**2
    J = x[0]**2+x[1]**2
    F2 = J-r2/2
    F4 = Y**2-r2*J/2+r2**2/12
    laplace = lambda f: sum(s.diff(f, t, 2) for t in x)
    equal("degree two homogeneous harmonic", laplace(F2), 0)
    equal("degree four homogeneous harmonic", laplace(F4), 0)
    projection2 = s.Rational(1, 2)
    if inject == "wrong_projection":
        projection2 = s.Rational(1, 3)
    equal("homogeneous square decomposition", Y**2, r2**2/6+projection2*r2*F2+F4)

    def mean(poly):
        out = 0
        for powers, coeff in s.Poly(s.expand(poly), *x).terms():
            if any(k % 2 for k in powers):
                continue
            aa = [k//2 for k in powers]
            out += coeff*s.prod(s.rf(s.Rational(1,2), k) for k in aa)/s.rf(2,sum(aa))
        return s.factor(out)

    for name, polynomial, expected in [
        ("sphere mean Y squared",Y**2,s.Rational(1,6)),
        ("sphere mean Y fourth",Y**4,s.Rational(3,40)),
        ("degree two norm mean",F2**2,s.Rational(1,12)),
        ("degree four norm mean",F4**2,s.Rational(19,720)),
        ("degree two zero mean",F2,0),
        ("degree four zero mean",F4,0),
        ("cross degree orthogonality",F2*F4,0),
        ("Y orthogonal to its square",Y**3,0),
        ("degree two projection coefficient",mean(Y**2*F2)/mean(F2**2),s.Rational(1,2)),
        ("degree four projection coefficient",mean(Y**2*F4)/mean(F4**2),1),
    ]:
        equal(name, mean(polynomial), expected)
    equal("projection squared norms sum", s.Rational(1,18)+s.Rational(1,24)+s.Rational(19,360),s.Rational(3,20))

    rho, p, ell, ch, H = s.symbols("rho_b p_b ell cosh_t H", nonzero=True)
    gain = 3/(rho**3*p)
    Ssq = gain**2*H**2/(ell**2*ch**3)
    measure_factor = rho**4*ch**3*ell
    equal("proper-volume cosh and width cancellation", s.pi**2*Ssq*measure_factor/3,
          3*s.pi**2*H**2/(rho**2*p**2*ell))
    m1, m2 = s.symbols("measure1 measure2")
    measure = 1+e*m1+e**2*m2/2
    equal("integrated tension retains unknown measure and response",
          s.diff(sigma_e*measure,e,2).subs(e,0),
          sig[1]*psi2+sig[2]*psi1**2+2*sig[1]*psi1*m1+sig[0]*m2)
    equal("squared response removes second jet and measure at order two",
          s.diff((eta_e-eta)**2*measure,e,2).subs(e,0), 2*psi1**2)

    Cmin, Cmax = s.Integer(606000000000), s.Integer(744000000000)
    sigmin, sigmax = 4-4*s.Rational(9,100000), 4-4*s.Rational(8,100000)
    bounds = {
        "sigma_second_closed_bounds":[str(sigmin),str(sigmax)],
        "scalar_sum_forcing_open_bounds":[str(-4*Cmax),str(-4*Cmin)],
        "per_copy_scalar_forcing_open_bounds":[str(-2*Cmax),str(-2*Cmin)],
        "tension_constitutive_open_bounds":[str(sigmin*Cmin),str(sigmax*Cmax)],
        "stress_trace_constitutive_open_bounds":[str(-4*sigmax*Cmax),str(-4*sigmin*Cmin)],
    }
    equal("tension lower rational product",sigmin*Cmin,s.Integer(2423781840000))
    equal("tension upper rational product",sigmax*Cmax,s.Integer(2975761920000))
    return checks,bounds


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inject",choices=INJECTIONS,default="none")
    args = parser.parse_args()
    result = {
        "schema":"hdblast-astra-material-forcing-symbolic-v1",
        "python_optimized":bool(sys.flags.optimize),
        "sympy_version":s.__version__,
        "injection":args.inject,
        "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "prior_code_imported_or_executed":False,
        "scope":"Specified symbolic identities and exact rational corollaries; no full second-order solve.",
    }
    try:
        names,bounds = run(args.inject)
        result.update(status="PASS",checks=len(names),check_names=names,derived_bounds=bounds)
        print(json.dumps(result,indent=2))
        return 0
    except FormulaFailure as exc:
        result.update(status="FAIL",failure=str(exc))
        print(json.dumps(result,indent=2))
        return 1


if __name__ == "__main__":
    sys.exit(main())
