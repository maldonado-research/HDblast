#!/usr/bin/env python3
"""Exact local current checks; no physical-completeness claim.

Run with the existing research environment. All checks raise on failure under
ordinary and optimized Python. Receipts identify this script and frozen inputs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import sympy as s


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("CURRENT_CHECKS.json"))
    args = parser.parse_args()
    checks = []

    def equal(name, expression):
        residual = s.factor(s.simplify(expression))
        if residual != 0:
            raise RuntimeError(f"FAILED {name}: {residual}")
        checks.append({"name": name, "kind": "identity", "residual": "0"})

    def rejects(name, expression):
        residual = s.factor(s.simplify(expression))
        if residual == 0:
            raise RuntimeError(f"WRONG FORMULA SURVIVED {name}")
        checks.append({"name": name, "kind": "wrong-formula control", "nonzero_residual": str(residual)})

    rho, p = s.symbols("rho p", nonzero=True, real=True)
    a, pp = s.symbols("a pprime", real=True)
    e1, e2 = s.symbols("e1 e2", real=True)
    ps1, ps2, n1, n2, x1, x2 = s.symbols("psi1 psi2 N1 N2 X1 X2", real=True)
    ps1p, ps2p, x1p, x2p = s.symbols("psi1prime psi2prime X1prime X2prime", real=True)
    dn1, dn2, dp1, dp2, dx1, dx2 = s.symbols("DN1 DN2 Dpsi1 Dpsi2 DX1 DX2", real=True)
    psi, lapse = e1*ps1+e2*ps2, e1*n1+e2*n2
    k = 4*(a+e1*ps1p+e2*ps2p)*s.exp(-lapse)

    def radial_potential(ni, psi_i, psi_ip, xi):
        return rho**4*s.exp(4*psi)*(k*(ni-psi_i)-4*psi_ip*s.exp(-lapse)
            -(p+e1*x1p+e2*x2p)*xi*s.exp(-lapse))

    omega_y = (s.diff(radial_potential(n2,ps2,ps2p,x2), e1)
               -s.diff(radial_potential(n1,ps1,ps1p,x1), e2)).subs({e1:0,e2:0})
    target = rho**4*(12*((ps1p-a*n1)*ps2-(ps2p-a*n2)*ps1)
        +x1*x2p-x2*x1p+p*(n1*x2-n2*x1)-4*p*(ps1*x2-ps2*x1))
    equal("nonlinear Lee-Wald radial potential field-space curl", omega_y-target)

    dN, dpsi, dpsip, nn, KK = s.symbols("deltaN deltapsi deltapsiprime n K", nonzero=True)
    equal("GHY potential cancellation",
        KK*(dN-dpsi)-4*dpsip/nn + (4*dpsi*KK+4*dpsip/nn-KK*dN)-3*KK*dpsi)

    def tangent_potential(ni, psi_i, xi, dni, dpi):
        return rho**2*s.exp(lapse+2*psi)*(-dni-3*dpi
            +(psi_i-ni)*(e1*dn1+e2*dn2)-(e1*dx1+e2*dx2)*xi)

    omega_mu = (s.diff(tangent_potential(n2,ps2,x2,dn2,dp2),e1)
              -s.diff(tangent_potential(n1,ps1,x1,dn1,dp1),e2)).subs({e1:0,e2:0})
    constrained_mu = omega_mu.subs({n1:-2*ps1,n2:-2*ps2,dn1:-2*dp1,dn2:-2*dp2})
    equal("direct tangential Lee-Wald current on generic slice",
        constrained_mu-rho**2*(x1*dx2-x2*dx1+6*(ps1*dp2-ps2*dp1)))

    constrained_y = omega_y.subs({n1:-2*ps1,n2:-2*ps2,
        ps1p:-2*a*ps1-p*x1/3,ps2p:-2*a*ps2-p*x2/3})
    physical_y = rho**4*(x1*x2p-x2*x1p-2*p*(ps1*x2-ps2*x1))
    equal("gravity-scalar cross cancellation",constrained_y-physical_y)
    rejects("omitting gravitational radial current",
        rho**4*(x1*x2p-x2*x1p-6*p*(ps1*x2-ps2*x1))-physical_y)

    l1,l2 = s.symbols("Lambda1 Lambda2",real=True)
    P1,P2,V1,V2 = s.symbols("P1 P2 P1prime P2prime",real=True)
    alpha=1/(rho**2*p**2)
    q=2/(3*rho**2)
    w=1/(rho**4*p**2)
    spectral_y=physical_y.subs({x1p:pp*x1/p+(3*l1/(rho**2*p)-2*p)*ps1,
        x2p:pp*x2/p+(3*l2/(rho**2*p)-2*p)*ps2})
    spectral_y=spectral_y.subs({ps1:P1/rho**2,ps2:P2/rho**2,
        x1:-3*V1/(rho**2*p),x2:-3*V2/(rho**2*p)})
    green=9*alpha*(l1*P1*V2-l2*V1*P2)
    equal("unequal-spectral physical radial current",spectral_y-green)
    equal("common-spectral Green normalization",green.subs(l2,l1)-9*l1*alpha*(P1*V2-P2*V1))
    rejects("common-Lambda shortcut for unequal profiles",green-9*l1*alpha*(P1*V2-P2*V1))
    equal("physical tangential energy polarization",
        rho**2*(9*V1*V2/(rho**4*p**2)+6*P1*P2/rho**4)-9*(alpha*V1*V2+q*P1*P2))

    # Differentiate alpha-weighted products, using only the two radial ODEs.
    derivative=9*((l1-l2)*alpha*V1*V2+l1*P1*(q-l2*w)*P2-l2*(q-l1*w)*P1*P2)
    equal("unequal-spectral current conservation",derivative+9*(l2-l1)*(alpha*V1*V2+q*P1*P2))
    beta,co=s.symbols("beta codd",positive=True)
    equal("even shell flux for unequal spectral values",green.subs({V1:l1*P1/(beta*alpha),V2:l2*P2/(beta*alpha)}))
    odd_flux=2*green.subs({V1:-co*P1/alpha,V2:-co*P2/alpha})
    shell_div=18*co*(l2-l1)*P1*P2
    equal("odd bulk-shell current conservation",odd_flux-shell_div)
    rejects("missing moving-shell current",odd_flux)
    zeta=s.symbols("zeta",real=True)
    edge_action=2*(3*a*rho**2*zeta**2+3*rho**2*zeta*(-a*zeta))-(6*a)*rho**2*zeta**2/2
    equal("doubled moving-domain GHY and single-shell kinetic accounting",edge_action+3*a*rho**2*zeta**2)

    F1,F2,F1x,F2x,G1,G2,sl=s.symbols("F1 F2 F1x F2x G1 G2 sL",real=True)
    # h_L^-1=rho^(3/2)|p| and (h_L^-1)'=(sL/rho)h_L^-1.
    factor=rho**s.Rational(3,2)*s.Abs(p)
    P1y=factor*(-F1x+sl*F1)/rho
    P2y=factor*(-F2x+sl*F2)/rho
    equal("Liouville Green current orientation",alpha*(factor*F1*P2y-factor*F2*P1y)+(F1*F2x-F2*F1x))
    equal("general packet Liouville current",alpha*(factor*G1*P2y-factor*G2*P1y)
        -(-G1*F2x+G2*F1x+sl*(G1*F2-G2*F1)))

    u,v=s.symbols("u v",positive=True)
    yy=s.sqrt(u*v)
    tt=s.log(v/u)/2
    jac=s.det(s.Matrix([[s.diff(yy,u),s.diff(yy,v)],[s.diff(tt,u),s.diff(tt,v)]]))
    equal("null coordinate density Jacobian",jac-1/(2*yy))
    Z=s.symbols("scalar_null_Wronskian",real=True)
    jy=v**4*Z/8
    jt=-v**4*Z/(8*yy)
    ju=(u/yy*jy-u*jt)*jac
    equal("true null current density",ju-v**3*Z/8)
    equal("fixed-radius pullback equals true cone pullback",jy/v-ju)
    rejects("omitting density transformation Jacobian",u/yy*jy-u*jt-v**3*Z/8)

    A=s.Function("A")(v)
    # P=y^4 A(v): partial_y at fixed t equals (u partial_u+v partial_v)/y.
    P=yy**4*A
    py=(u*s.diff(P,u)+v*s.diff(P,v))/yy
    equal("characteristic reconstruction constraint",py/yy**3-(4*A+v*s.diff(A,v)))
    equal("smooth leading metric uu coefficient",-6*(u*v*A)*(v/(2*yy))**2+s.Rational(3,2)*v**2*A)

    B=s.Function("B")(v)
    D=lambda f:4*f+v*s.diff(f,v)
    L=lambda f:4*f-4*v*s.diff(f,v)-v**2*s.diff(f,v,2)
    equal("characteristic mixed-spectral operator cancellation",L(A)*D(B)-L(B)*D(A)
        -v*(D(A)*s.diff(D(B),v)-D(B)*s.diff(D(A),v)))
    f=s.Function("f")(v)
    equal("nonzero smooth candidate-data Wronskian",f*s.diff(v*f,v)-v*f*s.diff(f,v)-f**2)

    frozen=Path(__file__).resolve().parents[1]/"basis"/"prior-checkpoint"
    source_names=["operator/PHYSICAL_OPERATOR_AUDIT_20261013.md",
        "operator/ODD_KINETIC_COMPLETION_20261013.md",
        "domain/SCALAR_SELF_ADJOINT_DOMAIN_AND_GAP.md",
        "manuscript/JACOBI_RESPONSE_TO_SPECTRAL_POSITIVITY.md"]
    receipt={"scope":"Exact local current identities, not a global phase-space theorem",
        "sympy_version":s.__version__,"script_sha256":sha256(Path(__file__)),
        "frozen_source_sha256":{name:sha256(frozen/name) for name in source_names},
        "passed":len(checks),"identities":sum(x["kind"]=="identity" for x in checks),
        "wrong_formula_controls":sum(x["kind"]=="wrong-formula control" for x in checks),
        "checks":checks}
    args.output.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(f"PASS {receipt['passed']} checks ({receipt['identities']} identities, "
          f"{receipt['wrong_formula_controls']} wrong-formula controls).")
    print(f"Receipt: {args.output}")


if __name__ == "__main__":
    main()
