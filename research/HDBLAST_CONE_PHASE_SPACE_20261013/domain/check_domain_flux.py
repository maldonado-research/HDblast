#!/usr/bin/env python3
"""Exact controls for domain/null-limit audit; no full-background existence claim."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    checks = []

    def zero(name: str, expr: sp.Expr) -> None:
        reduced = sp.simplify(expr)
        if reduced != 0:
            raise RuntimeError(f"{name}: nonzero residual {reduced}")
        checks.append({"name": name, "result": "pass"})

    def reject(name: str, expr: sp.Expr) -> None:
        reduced = sp.simplify(expr)
        if reduced == 0:
            raise RuntimeError(f"{name}: deliberately wrong formula was accepted")
        checks.append({"name": name, "result": "rejected as intended", "residual": str(reduced)})

    y, x, t, v, u = sp.symbols("y x t v u", positive=True)
    kappa = sp.symbols("kappa", positive=True)
    r2, p2 = sp.symbols("r2 p2", real=True)
    rho = y * (1 + r2 * y**2)
    p = kappa * y * (1 + p2 * y**2)
    s = rho * (sp.Rational(3, 2) * sp.diff(rho, y) / rho + sp.diff(p, y) / p)
    VL = -rho * sp.diff(s, y) + s**2 + sp.Rational(2, 3) * rho**2 * p**2
    zero("s expansion through y^2", sp.series(s, y, 0, 4).removeO() - sp.Rational(5, 2) - (sp.Rational(11, 2) * r2 + 2 * p2) * y**2)
    zero("VL expansion through y^2", sp.series(VL, y, 0, 4).removeO() - sp.Rational(25, 4) - (sp.Rational(33, 2) * r2 + 6 * p2) * y**2)
    zero("h expansion through relative y^2", sp.series((1 + r2*y**2)**sp.Rational(-3, 2)/(1+p2*y**2), y, 0, 4).removeO() - (1-(sp.Rational(3, 2)*r2+p2)*y**2))
    zero("x_y expansion", sp.series(-1/rho, y, 0, 3).removeO() - sp.diff(-sp.log(y) + r2*y**2/2, y))
    S = sp.Rational(3, 2) * sp.diff(rho, y)
    VT = S**2 + rho * sp.diff(S, y)
    zero("tensor potential expansion", sp.series(VT, y, 0, 4).removeO() - sp.Rational(9, 4) - sp.Rational(45, 2) * r2 * y**2)
    reject("wrong scalar Liouville threshold 9/4", sp.Rational(25, 4)-sp.Rational(9, 4))
    reject("wrong first s correction", (sp.Rational(11, 2)*r2+2*p2)-(3*r2+2*p2))

    # Exact coefficient normalization behind the orientation of the Green form.
    R, K = sp.symbols("R K", positive=True)
    alpha = 1/(R**2*K**2)
    h = 1/(R**sp.Rational(3, 2)*K)
    zero("alpha/(rho*h^2)=1", alpha/(R*h**2)-1)
    f1, f2, fx1, fx2, sx = sp.symbols("f1 f2 fx1 fx2 sx")
    P1, P2 = f1/h, f2/h
    Py1, Py2 = -(fx1-sx*f1)/(R*h), -(fx2-sx*f2)/(R*h)
    zero("Green orientation", alpha*(P1*Py2-Py1*P2)+(f1*fx2-fx1*f2))
    reject("wrong Green orientation", alpha*(P1*Py2-Py1*P2)-(f1*fx2-fx1*f2))

    # Removing global dS friction, for any radial potential/operator term.
    Phi = sp.Function("Phi")(t)
    temporal_scale = sp.cosh(t)**sp.Rational(-3, 2)
    psi = temporal_scale * Phi
    temporal = (sp.diff(psi, t, 2)+3*sp.tanh(t)*sp.diff(psi, t))/temporal_scale
    target = sp.diff(Phi, t, 2)+(-sp.Rational(9, 4)+sp.Rational(3, 4)/sp.cosh(t)**2)*Phi
    zero("global dS friction removal", sp.trigsimp(temporal-target))
    reject("wrong residual threshold cancellation", sp.Rational(25, 4)-4-sp.Rational(3, 2))

    A1, A2 = sp.Function("A1")(v), sp.Function("A2")(v)
    def D(A):
        return 4*A+v*sp.diff(A, v)
    def L(A):
        return 4*A-4*v*sp.diff(A, v)-v**2*sp.diff(A, v, 2)
    zero("characteristic D/L identity", L(A1)-D(A1)+v*sp.diff(D(A1), v))
    zero("characteristic physical-current identity", L(A1)*D(A2)-L(A2)*D(A1)-v*(D(A1)*sp.diff(D(A2), v)-D(A2)*sp.diff(D(A1), v)))
    X1, X2 = -3*D(A1)/kappa, -3*D(A2)/kappa
    actual_flux = 9*v**2/(8*kappa**2)*(L(A1)*D(A2)-L(A2)*D(A1))
    scalar_null_flux = v**3/8*(X1*sp.diff(X2, v)-X2*sp.diff(X1, v))
    zero("leading grav-scalar current equals scalar null current", actual_flux-scalar_null_flux)
    reject("wrong equal-Lambda packet replacement", actual_flux-9*v**2/(8*kappa**2)*(A1*D(A2)-A2*D(A1)))
    Aconst, Alin = -kappa/12, -kappa*v/30
    zero("leading X1 comparator", -3*D(Aconst)/kappa-1)
    zero("leading X2 comparator", -3*D(Alin)/kappa-v/2)
    zero("leading comparator flux", 9*v**2/(8*kappa**2)*(L(Aconst)*D(Alin)-L(Alin)*D(Aconst))-v**3/16)

    # Check radial indicial expression directly in the flat leading model.
    sigma, Lambda = sp.symbols("sigma Lambda")
    zero("scalar indicial polynomial", sp.expand(-sigma*(sigma-1)+4*sigma-Lambda)-(-sigma**2+5*sigma-Lambda))
    P_of_y = y**4 * sp.Function("A")(y*sp.exp(t))
    leading_tau = -y**2*sp.diff(P_of_y,y,2)+4*y*sp.diff(P_of_y,y)
    # Replacing t derivatives by the same symbolic composite is simpler and exact.
    zero("flat leading radial/tangential equation", leading_tau-(4*P_of_y-3*sp.diff(P_of_y,t)-sp.diff(P_of_y,t,2)))

    eps = sp.symbols("eps", positive=True)
    t_of_v = sp.log(v/eps)
    ch = (v/eps+eps/v)/2
    sh = (v/eps-eps/v)/2
    timelike_flux = eps**4*ch**3*sh/v
    zero("flat exact timelike-to-null density", timelike_flux-(v+eps**2/v)**3*(v-eps**2/v)/(16*v))
    zero("flat exact characteristic flux limit", sp.limit(timelike_flux,eps,0,dir="+")-v**3/16)
    primitive = v**4/64+eps**2*v**2/16+eps**6/(16*v**2)+eps**8/(64*v**4)
    zero("exact flux primitive", sp.diff(primitive,v)-timelike_flux)
    zero("fixed-time flux primitive", sp.diff(sp.cosh(t)**4/4,t)-sp.cosh(t)**3*sp.sinh(t))
    zero("null scalar flux convention", ((v/2)**3)*(sp.Rational(1,2))-v**3/16)
    reject("wrong interchange of whole-time and finite-time limits", (sp.Rational(2)**4-1)/64)

    # Actual differentiated densities, not copied thresholds: near zero an
    # integral of y**b is finite precisely when b > -1.
    decay = sp.symbols("decay", positive=True)
    profile = y**(-decay)
    radial_norm_density = y**2*profile**2
    radial_gradient_density = y**4*sp.diff(profile,y)**2
    inertial_gradient_density = y**3*sp.diff(profile,y)**2
    def power(density):
        return sp.simplify(y*sp.diff(density,y)/density)
    zero("derived radial norm exponent",power(radial_norm_density)-(2-2*decay))
    zero("derived lapse-weighted radial gradient exponent",power(radial_gradient_density)-(2-2*decay))
    zero("derived inertial gradient exponent",power(inertial_gradient_density)-(1-2*decay))
    radial_threshold = sp.solve(power(radial_gradient_density)+1,decay)[0]
    inertial_threshold = sp.solve(power(inertial_gradient_density)+1,decay)[0]
    reject("wrong identical radial and inertial energy thresholds",inertial_threshold-radial_threshold)

    # Concrete values demonstrate the two limits. These are controls, not
    # a numerical approximation of the registered background.
    eps_values = [sp.Rational(1,2),sp.Rational(1,4),sp.Rational(1,10),sp.Rational(1,100),sp.Rational(1,1000)]
    moving = primitive.subs(v,2)-primitive.subs(v,1)
    fixed = eps**4*(sp.cosh(1)**4-1)/4
    rows = [{"epsilon": float(e), "fixed_t_0_to_1":float(fixed.subs(eps,e)), "moving_v_1_to_2":float(moving.subs(eps,e)), "moving_limit":float(sp.Rational(15,64))} for e in eps_values]

    source_root = Path(__file__).resolve().parents[1]/"basis"/"prior-checkpoint"
    source_names = ["manuscript/JACOBI_RESPONSE_TO_SPECTRAL_POSITIVITY.md","domain/SCALAR_SELF_ADJOINT_DOMAIN_AND_GAP.md","operator/PHYSICAL_OPERATOR_AUDIT_20261013.md","tensor/TENSOR_PARTNER_GAP.md","NEXT_CALCULATIONS.md"]
    pins = [{"path":name,"sha256":hashlib.sha256((source_root/name).read_bytes()).hexdigest()} for name in source_names]
    note = Path(__file__).with_name("DOMAIN_AND_CHARACTERISTIC_FLUX_AUDIT.md")
    report = {"scope":"Exact local identities and noncommuting cone-limit controls; no registered hyperbolic existence/completeness theorem.","sympy_version":sp.__version__,"checks":checks,"check_count":len(checks),"all_passed":True,"flux_controls":rows,"source_pins":pins,"checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"note_sha256":hashlib.sha256(note.read_bytes()).hexdigest()}
    if args.output:
        args.output.write_text(json.dumps(report,indent=2)+"\n")
    print(f"{len(checks)} exact identity/wrong-formula controls passed (SymPy {sp.__version__}).")
    for row in rows:
        print(f"epsilon={row['epsilon']:.3g}: fixed-time={row['fixed_t_0_to_1']:.10g}, moving-time={row['moving_v_1_to_2']:.10g}, null-limit={row['moving_limit']:.10g}")


if __name__ == "__main__":
    main()
