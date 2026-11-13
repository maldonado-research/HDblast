#!/usr/bin/env python3
"""Finite exact algebra witnesses, not an analytic existence/bounds proof.

Production formulas use parameter derivatives (the e^2 coefficient is half the
second derivative).  Oracles differentiate actual polynomial/rational families.
The metric oracle uses the coordinate Christoffel/Ricci definition, while the
production calculation uses background-covariant connection differences.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import sys

import mpmath
import sympy as s

SCHEMA = "hdblast-proper-second-geometry-formal-claims-v1"
FACTOR_NAMES = (
    "scalar_cross", "defect_transport", "defect_velocity", "quotient_cross",
    "inverse_quadratic", "connection_mixed", "ricci_quadratic",
)
F = {name: s.Symbol("factor_" + name) for name in FACTOR_NAMES}
e = s.Symbol("e")


class VerificationError(ValueError):
    pass


def require(ok, message):
    if not ok:
        raise VerificationError(message)


def canonical(value):
    return (json.dumps(value, sort_keys=True, indent=2) + "\n").encode()


def digest(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def clean(expr):
    return s.cancel(s.expand(expr))


def d2(expr):
    return clean(s.diff(expr, e, 2).subs(e, 0))


def trunc(expr):
    p = s.Poly(s.expand(expr), e)
    return s.Add(*(coefficient * e**power[0]
                   for power, coefficient in p.terms() if power[0] <= 2))


def scalar_checks():
    x, y, a, b, c, d = s.symbols("x y a b c d")
    coeff = s.symbols("w0:12")
    w1 = coeff[0] + coeff[1]*x + coeff[2]*y + coeff[3]*x*x + coeff[4]*x*y + coeff[5]*y*y
    w2 = coeff[6] + coeff[7]*x + coeff[8]*y + coeff[9]*x*x + coeff[10]*x*y + coeff[11]*y*y
    path = {x: x + e*a + e*e*b/2, y: y + e*c + e*e*d/2}
    scalar = 17 + e*w1 + e*e*w2/2
    oracle = d2(scalar.subs(path, simultaneous=True))
    production = w2 + F["scalar_cross"]*(a*s.diff(w1, x) + c*s.diff(w1, y))
    # A nonconstant background must fail the constant-shell formula.
    nonconstant = d2((x*x + y + e*w1 + e*e*w2/2).subs(path, simultaneous=True))
    background_missing = clean(nonconstant - oracle)
    require(background_missing != 0, "scalar background premise witness is degenerate")
    return [("constant_shell_scalar_pullback", [production-oracle])], {
        "nonconstant_background_omitted_terms": str(background_missing),
        "scope": "Generic quadratic scalar jets and generic second-order chart path; shell background is constant.",
    }


def defect_checks():
    x, y = s.symbols("x y")
    a, b, c, d, u, v, r, z, u2, v2 = s.symbols("a b c d u v r z u2 v2")
    E1 = s.Matrix([[1+x*x+y, x*y+2*x], [x*y+2*x, 3+x-y*y]])
    E2 = s.Matrix([[x-y, x*x+y], [x*x+y, 1+x*y]])
    xi = s.Matrix([a, c]); v0 = s.Matrix([u, v]); vc1 = s.Matrix([r, z])
    path = {x: x+e*a+e*e*b/2, y: y+e*c+e*e*d/2}
    vel = v0+e*vc1+e*e*s.Matrix([u2, v2])/2
    tensor = (e*E1+e*e*E2/2).subs(path, simultaneous=True)
    oracle = d2((vel.T*tensor*vel)[0])
    # A fixed nonzero torsion-free background connection.  This algebraic
    # conversion holds for every such connection, including the actual LC one.
    gamma = [[[s.Rational(1+i+2*j+3*k, 7) for k in range(2)]
              for j in range(2)] for i in range(2)]
    for i in range(2):
        gamma[i][1][0] = gamma[i][0][1]
    correction = s.Matrix([sum(gamma[i][j][k]*xi[j]*v0[k]
                               for j in range(2) for k in range(2)) for i in range(2)])
    eta = vc1+correction
    covE1 = s.zeros(2)
    for i in range(2):
        for j in range(2):
            covE1[i,j] = sum(xi[k]*(s.diff(E1[i,j], (x,y)[k])-
                sum(gamma[l][k][i]*E1[l,j]+gamma[l][k][j]*E1[i,l] for l in range(2))) for k in range(2))
    production = ((v0.T*E2*v0)[0] + F["defect_transport"]*(v0.T*covE1*v0)[0]
                  + F["defect_velocity"]*(eta.T*E1*v0)[0])
    coordinate_mixed = ((v0.T*E2*v0)[0]+2*(v0.T*covE1*v0)[0]+4*(vc1.T*E1*v0)[0])
    mismatch = clean(coordinate_mixed-oracle)
    require(mismatch != 0, "covariant eta conversion witness is degenerate")
    return [("background_defect_second_contraction", [production-oracle])], {
        "coordinate_velocity_in_covariant_formula_rejected": True,
        "scope": "E0 is identically zero; tensor and path jets are formal polynomials. Nonzero background connection tests eta conversion.",
    }


def quotient_checks():
    S, D0, D1, D2, B2, J = s.symbols("S D0 D1 D2 B2 J", nonzero=True)
    # Finite weighted integration is encoded by its independent moments.  A0=0
    # and A1=S*Y imply N0=0, N1=S D0, N2=B2+2 S D1 exactly.
    # The direct rational quotient is differentiated, not its derivative rule.
    numerator = e*S*D0+e*e*(B2+2*S*D1)/2
    denominator = D0+e*D1+e*e*D2/2
    oracle = d2(numerator/denominator)
    production = (B2+2*S*D1)/D0-F["quotient_cross"]*(S*D0)*D1/D0**2
    # For a general first mode, J=int Y A1 j1 need not equal S D1.
    general = d2((e*S*D0+e*e*(B2+2*J)/2)/denominator)
    missing = clean(general-B2/D0)
    require(missing != 0, "generic first-mode area witness is degenerate")
    # A second oracle starts before integration: multiply the actual scalar
    # and area polynomials at four weighted nodes, sum, then differentiate the
    # resulting rational quotient.  In particular its cross factor is not
    # supplied through a preconstructed N2 coefficient.
    nodeY = [-2,-1,1,3]
    nodeW = [1,2,3,5]
    j1 = s.symbols("j1_0:4")
    j2 = s.symbols("j2_0:4")
    a2 = s.symbols("a2_0:4")
    areas = [1+e*j1[i]+e*e*j2[i]/2 for i in range(4)]
    fields = [e*S*nodeY[i]+e*e*a2[i]/2 for i in range(4)]
    nodeN = sum(nodeW[i]*nodeY[i]*fields[i]*areas[i] for i in range(4))
    nodeD = sum(nodeW[i]*nodeY[i]**2*areas[i] for i in range(4))
    nodeD0 = nodeD.subs(e,0)
    nodeD1 = s.diff(nodeD,e).subs(e,0)
    nodeB2 = sum(nodeW[i]*nodeY[i]*a2[i] for i in range(4))
    node_oracle = d2(nodeN/nodeD)
    node_production = ((nodeB2+2*S*nodeD1)/nodeD0-
                       F["quotient_cross"]*S*nodeD1/nodeD0)
    return [("selected_normalized_moment_area_cross_cancellation", [production-oracle, oracle-B2/D0]),
            ("weighted_node_scalar_area_quotient_direct_expansion", [node_production-node_oracle,node_oracle-nodeB2/nodeD0])], {
        "generic_first_mode_retained_area_term": str(missing),
        "scope": "Exact rational moment algebra. No nonzero-amplitude denominator margin is asserted.",
    }


def geometry_checks():
    coords = s.symbols("t x y z")
    t,x,y,z = coords
    n = 4
    g0 = s.diag(-1, 1+t+t*t, 1, 1)
    p1 = s.Matrix([[t+x+x*x, 1+t*y, 1+t+z, x*y],
                   [1+t*y, x+y*y, y+t*z, t*z],
                   [1+t+z, y+t*z, 1+x*z, y],
                   [x*y, t*z, y, t+x*x]])
    p2 = s.Matrix([[1+t*y, t+x, y*z, z],
                   [t+x, 1+z*z, t+y, x],
                   [y*z, t+y, x+y*z, 1+t],
                   [z, x, 1+t, 1+y]])
    at = {q: s.Integer(0) for q in coords}
    A = g0.inv()
    Gamma0 = [[[clean(sum(A[k,l]*(s.diff(g0[l,j],coords[i])+s.diff(g0[l,i],coords[j])-
                                  s.diff(g0[i,j],coords[l])) for l in range(n))/2)
                for j in range(n)] for i in range(n)] for k in range(n)]

    def S(p,d,i,j):
        def cov(k,a,b):
            return (s.diff(p[a,b],coords[k]) - sum(Gamma0[l][k][a]*p[l,b]+
                    Gamma0[l][k][b]*p[a,l] for l in range(n)))
        return cov(i,d,j)+cov(j,d,i)-cov(d,i,j)

    C1 = [[[clean(sum(A[k,l]*S(p1,l,i,j) for l in range(n))/2)
            for j in range(n)] for i in range(n)] for k in range(n)]
    C2 = [[[clean(sum(A[k,l]*S(p2,l,i,j) for l in range(n))/2 -
                    F["connection_mixed"]*sum(A[k,l]*p1[l,m]*C1[m][i][j]
                    for l in range(n) for m in range(n)))
            for j in range(n)] for i in range(n)] for k in range(n)]

    def covC(C,a,k,i,j):
        return s.diff(C[k][i][j],coords[a])+sum(Gamma0[k][a][l]*C[l][i][j]-
            Gamma0[l][a][i]*C[k][l][j]-Gamma0[l][a][j]*C[k][i][l] for l in range(n))

    Ric2 = s.Matrix(n,n,lambda i,j: clean(sum(covC(C2,a,a,i,j)-covC(C2,j,a,i,a) for a in range(n))+
        F["ricci_quadratic"]*sum(C1[a][a][b]*C1[b][i][j]-C1[a][j][b]*C1[b][i][a]
                                for a in range(n) for b in range(n))))
    B2 = F["inverse_quadratic"]*A*p1*A*p1*A-A*p2*A

    # Independent coordinate oracle.  Invert the actual e-dependent 4x4
    # matrix first by SymPy's exact rational inverse.  Truncation only removes
    # orders >=3.  No production inverse/C1/C2/Ric2 formula enters this path.
    ge = g0+e*p1+e*e*p2/2
    G = ge.subs(at)
    direct_inverse = G.inv()
    oracle_inv2 = direct_inverse.applyfunc(d2)
    inv = direct_inverse.applyfunc(lambda f: clean(f.subs(e,0))+clean(s.diff(f,e).subs(e,0))*e+d2(f)*e*e/2)
    dg = [ge.diff(q).subs(at) for q in coords]
    ddg = [[ge.diff(q).diff(r).subs(at) for r in coords] for q in coords]
    dinv = [(-inv*dg[a]*inv).applyfunc(trunc) for a in range(n)]

    def W(l,i,j):
        return dg[i][l,j]+dg[j][l,i]-dg[l][i,j]

    gamma = [[[trunc(sum(inv[k,l]*W(l,i,j) for l in range(n))/2)
               for j in range(n)] for i in range(n)] for k in range(n)]

    def dgamma(a,k,i,j):
        return trunc(sum(dinv[a][k,l]*W(l,i,j)+inv[k,l]*(ddg[a][i][l,j]+
                     ddg[a][j][l,i]-ddg[a][l][i,j]) for l in range(n))/2)

    ric = s.Matrix(n,n,lambda i,j: trunc(sum(dgamma(a,a,i,j)-dgamma(j,a,i,a) for a in range(n))+
        sum(gamma[a][a][b]*gamma[b][i][j]-gamma[a][j][b]*gamma[b][i][a]
            for a in range(n) for b in range(n))))
    # Background curvature is explicitly nonzero, so this does not accidentally
    # reduce all background-covariant connection bookkeeping to a flat chart.
    background_ric = ric.subs(e,0)
    require(any(v != 0 for v in background_ric), "metric witness background is flat")
    inv_residual = list(B2.subs(at)-oracle_inv2)
    connection_residual = [C2[k][i][j].subs(at)-d2(gamma[k][i][j])
                           for k in range(n) for i in range(n) for j in range(n)]
    ricci_residual = list(Ric2.subs(at)-ric.applyfunc(d2))
    # Check the Einstein-defect subtraction separately with a generic constant
    # lambda=3.  The witness background itself need not satisfy E0=0; the
    # E0=0 contraction algebra is checked independently above.
    E2_production = Ric2.subs(at)-3*p2.subs(at)
    E2_oracle = (ric-3*ge.subs(at)).applyfunc(d2)
    return [("four_dimensional_inverse_second_jet", inv_residual),
            ("four_dimensional_connection_second_jet", connection_residual),
            ("four_dimensional_Ricci_second_jet", ricci_residual),
            ("four_dimensional_Einstein_defect_second_jet", list(E2_production-E2_oracle))], {
        "dimension": n,
        "background_signature_at_germ": [-1,1,1,1],
        "background_Ricci_at_germ": [[str(v) for v in row] for row in background_ric.tolist()],
        "offblock_metric_perturbations": True,
        "coordinate_germ": ["t=0","x=0","y=0","z=0"],
        "scope": "Exact formal e jets of a four-dimensional polynomial metric germ with a curved Lorentzian background; finite witness, not a general analytic proof.",
    }


def conformal_check():
    c=s.Symbol("c",nonzero=True)
    exact_R=-3*s.exp(-2*c*e)
    oracle=d2(exact_R)
    # For g_e=exp(2ce)g0, Ric(g0)=3g0, V_e=exp(-ce)V0:
    # E2(V0,V0)=12c^2, nabla E1=0, E1(eta1,V0)=-6c^2.
    production=12*c*c-F["defect_velocity"]*6*c*c
    return [("independent_constant_conformal_scaling_R2",[production-oracle])], {
        "exact_R_family":"-3*exp(-2*c*e)","exact_R2":"-12*c**2",
        "scope":"Auxiliary algebraic scaling family only; not the registered initial-anchor family or its quiet collar.",
    }


def read_claims(path):
    claims = json.loads(Path(path).read_text())
    validate_claims(claims)
    return claims


def validate_claims(claims):
    require(set(claims) == {"schema","factors"}, "invalid claims keys")
    require(claims["schema"] == SCHEMA, "invalid claims schema")
    require(set(claims["factors"]) == set(FACTOR_NAMES), "invalid factor keys")
    require(all(type(v) is int for v in claims["factors"].values()), "factors must be integers")


def evaluate_residuals(checks, claims):
    replacements = {F[name]: value for name,value in claims["factors"].items()}
    results = []
    failures = []
    for name, residuals in checks:
        bad = [clean(v.subs(replacements)) for v in residuals]
        bad = [v for v in bad if v != 0]
        results.append({"name":name,"components":len(residuals),"status":"pass" if not bad else "fail"})
        if bad:
            # Small exact witness: first substitute deterministic rational data,
            # retaining the full symbolic nonzero decision independently.
            expr = bad[0]
            symbols = sorted(expr.free_symbols,key=str)
            value = clean(expr.subs({q:s.Rational(i+2,i+3) for i,q in enumerate(symbols)}))
            if value == 0:
                value = expr
            failures.append({"check":name,"nonzero_residual_witness":str(value)})
    return results, failures


def build_checks():
    checks=[]; premises={}
    for name, function in (("scalar",scalar_checks),("defect",defect_checks),
                           ("quotient",quotient_checks),("metric",geometry_checks),
                           ("conformal",conformal_check)):
        new, metadata = function()
        checks.extend(new); premises[name]=metadata
    return checks,premises


def replay(claims, checks=None, premises=None):
    validate_claims(claims)
    require(s.__version__ == "1.14.0", "requires pinned sympy==1.14.0")
    require(mpmath.__version__ == "1.3.0", "requires pinned mpmath==1.3.0")
    if checks is None:
        checks,premises=build_checks()
    results,failures=evaluate_residuals(checks,claims)
    return {"schema":"hdblast-proper-second-geometry-formal-replay-v1",
            "status":"pass" if not failures else "fail", "engine":"sympy==1.14.0",
            "arithmetic":"exact rational and polynomial; no floating-point exponentiation",
            "claims_sha256":digest(claims),"checks":results,"failures":failures,
            "premises_and_scope":premises,
            "analytic_scope":"Only finite algebraic identities and fault witnesses are replayed. This does not prove PDE bounds, existence, source authentication, family selection, domain continuation, or a finite-amplitude radius."}


def controls(claims, checks, premises):
    require(replay(claims,checks,premises)["status"] == "pass", "controls require a passing baseline")
    records=[]
    for name in FACTOR_NAMES:
        mutated=copy.deepcopy(claims)
        mutated["factors"][name] += 1
        # Rehashing authenticates the altered bytes only.  The exact independent
        # mathematical oracle remains unchanged and must reject the bad factor.
        self_hash=digest(mutated)
        report=replay(mutated,checks,premises)
        require(report["status"] == "fail", "mathematical factor fault survived: "+name)
        records.append({"fault":name+"_plus_one", "self_rehashed_claims_sha256":self_hash,
                        "self_hash_valid":self_hash == digest(mutated),
                        "rejected_by":"independent exact mathematical residual",
                        "failures":report["failures"]})
    # This separate mutation explicitly recomputes the whole oracle after the
    # changed claim has been serialized and rehashed.  It demonstrates that a
    # refreshed content hash cannot make a mathematical factor error valid.
    mutated=copy.deepcopy(claims)
    mutated["factors"]["connection_mixed"] = 1
    serialized=canonical(mutated)
    reloaded=json.loads(serialized)
    self_hash=hashlib.sha256(serialized).hexdigest()
    fresh_report=replay(reloaded)
    require(fresh_report["status"] == "fail", "fresh rehashed connection factor fault survived")
    records.append({"fault":"connection_mixed_one_fresh_oracle", "self_rehashed_claims_sha256":self_hash,
                    "self_hash_valid":self_hash == digest(reloaded),
                    "oracle_recomputed_after_rehash":True,
                    "rejected_by":"independent exact mathematical residual", "failures":fresh_report["failures"]})
    return {"schema":"hdblast-proper-second-geometry-formal-controls-v1",
            "status":"pass","controls":records,
            "scope":"All seven independently rehashed mathematical coefficient mutations are rejected after exact oracle comparison; these are not hash-failure controls."}


def evaluate(symbolic_claims, check):
    """Closed-protocol adapter; check(name, boolean) must reject false checks.

    The caller supplies the complete claims dictionary.  This function always
    builds the independent oracle afresh and raises ValueError for bad math,
    independently of the caller's check callback behavior.
    """
    report=replay(symbolic_claims)
    for result in report["checks"]:
        check("formal_geometry/"+result["name"],result["status"] == "pass")
    if report["status"] != "pass":
        raise ValueError("formal geometry oracle rejects claims: "+json.dumps(report["failures"],sort_keys=True))
    return report


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--claims",type=Path,default=Path(__file__).with_name("claims.json"))
    parser.add_argument("--report",type=Path)
    parser.add_argument("--controls",type=Path)
    args=parser.parse_args(argv)
    try:
        claims=read_claims(args.claims)
        checks,premises=build_checks()
        report=replay(claims,checks,premises)
        payload=canonical(report)
        if args.report:
            args.report.write_bytes(payload)
        if args.controls:
            args.controls.write_bytes(canonical(controls(claims,checks,premises)))
        sys.stdout.buffer.write(payload)
        return 0 if report["status"] == "pass" else 1
    except (VerificationError,ValueError,KeyError,OSError) as exc:
        sys.stderr.write("verification failed: "+str(exc)+"\n")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
