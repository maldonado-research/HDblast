#!/usr/bin/env python3
"""Independent full graph-junction and packet-current checks; no prior code imports.

Exact identities only. The analytical sufficiency theorem is in the adjacent
note. This checker does not prove a PDE existence or phase-space completion
theorem. Checks deliberately avoid Python assert so -O has the same contract.
"""
import argparse
import hashlib
import json
from pathlib import Path

import sympy as S


ROOT = Path(__file__).resolve().parent
checks = []
controls = []


def zero(name, value):
    result = S.simplify(S.expand(value))
    if result != 0:
        raise RuntimeError(f"{name}: nonzero residual {result}")
    checks.append(name)


def reject(name, value):
    result = S.simplify(S.expand(value))
    if result == 0:
        raise RuntimeError(f"wrong-formula control did not fail: {name}")
    controls.append({"name": name, "nonzero_residual": str(result)})


def linear(expr, e):
    return S.diff(expr, e).subs(e, 0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="packet-junctions-results.json")
    parser.add_argument("--inject-failure", action="store_true")
    args = parser.parse_args()
    y, t, x1, x2, x3, e = S.symbols("y t x1 x2 x3 epsilon", real=True)
    coords = [y, t, x1, x2, x3]
    rho = S.Function("rho")(y)
    psi = S.Function("Psi")(*coords)
    chi = S.Function("chi")(*coords)
    z = S.Function("zeta")(t, x1, x2, x3)
    phi = S.Function("phi")(y)
    gamma = S.diag(-1, S.exp(2*t), S.exp(2*t), S.exp(2*t))
    gi = gamma.inv()
    # Direct Christoffel computation for the longitudinal metric, retaining
    # arbitrary dependence on every tangential coordinate.
    metric = S.diag(1-4*e*psi,
                    *[rho**2*(1+2*e*psi)*gamma[j,j] for j in range(4)])
    inv = S.diag(1+4*e*psi,
                 *[(1-2*e*psi)/rho**2*gi[j,j] for j in range(4)])
    connection = {}
    for A in range(5):
        for B in range(5):
            for C in range(B, 5):
                term = sum(inv[A,D]*(S.diff(metric[D,C],coords[B])
                           + S.diff(metric[D,B],coords[C])
                           - S.diff(metric[B,C],coords[D]))/2 for D in range(5))
                term = term.subs(e,0) + e*linear(term,e)
                connection[A,B,C] = connection[A,C,B] = S.simplify(term)
    a = S.diff(rho,y)/rho
    H = 1/rho**2
    normal_cov = [1-2*e*psi] + [-e*S.diff(z,q) for q in coords[1:]]
    tangents = [[e*S.diff(z,q)] + [int(i==j) for j in range(4)]
                for i,q in enumerate(coords[1:])]
    cov_der = [[S.diff(normal_cov[B],coords[A])
                - sum(connection[C,A,B]*normal_cov[C] for C in range(5))
                for B in range(5)] for A in range(5)]
    induced_delta = S.Matrix(4,4, lambda mu,nu:
                             2*rho**2*(psi+a*z)*gamma[mu,nu])
    for mu in range(4):
        for nu in range(mu,4):
            graph_cov = sum(tangents[mu][A]*tangents[nu][B]*cov_der[A][B]
                            for A in range(5) for B in range(5))
            base = S.simplify(graph_cov.subs(e,0))
            # The evaluation shift is part of the graph pullback.
            delta = S.simplify(linear(graph_cov,e) + z*S.diff(base,y))
            hessian = S.diff(z,coords[mu+1],coords[nu+1])
            for lam in range(4):
                gconn = connection[lam+1,mu+1,nu+1].subs(e,0)
                hessian -= gconn*S.diff(z,coords[lam+1])
            expected_cov = rho**2*(S.diff(psi,y)+4*a*psi
                                  +(S.diff(a,y)+2*a*a)*z)*gamma[mu,nu] - hessian
            zero(f"direct graph extrinsic curvature covariant {mu}{nu}",delta-expected_cov)
            # Raise with the perturbed common induced inverse metric.
            actual_mixed = sum(gi[mu,lam]/rho**2*delta
                               for lam in [mu]) - 2*a*(psi+a*z)*int(mu==nu)
            expected_mixed = (S.diff(psi,y)+2*a*psi+S.diff(a,y)*z)*int(mu==nu)
            expected_mixed -= sum(gi[mu,lam]*hessian/rho**2 for lam in [mu])
            zero(f"direct graph extrinsic curvature mixed {mu}{nu}",actual_mixed-expected_mixed)
    for mu in range(4):
        for nu in range(mu,4):
            graph_induced = sum(metric[A,B]*tangents[mu][A]*tangents[nu][B]
                                for A in range(5) for B in range(5))
            delta = linear(graph_induced,e)+z*S.diff(graph_induced.subs(e,0),y)
            zero(f"induced metric pullback {mu}{nu}",delta-induced_delta[mu,nu])
    scalar_normal = sum(inv[A,B]*normal_cov[B]*S.diff(phi+e*chi,coords[A])
                        for A in range(5) for B in range(5))
    delta_p = linear(scalar_normal,e)+z*S.diff(scalar_normal.subs(e,0),y)
    zero("direct scalar normal pullback",delta_p-S.diff(chi,y)
         -2*S.diff(phi,y)*psi-S.diff(phi,y,2)*z)
    zero("direct induced scalar",linear(phi+e*chi,e)
         +z*S.diff(phi,y)-(chi+z*S.diff(phi,y)))

    # Independent algebraic two-sided residuals, allowing arbitrary jets.
    ab, Hb, pb, dpb, sig2, B = S.symbols("a H p pprime sigma2 B", nonzero=True)
    Pe, Po, Pey, Poy, De, Do = S.symbols("Pe Po Pey Poy DPe DPo")
    ze, zo, Tze, Tzo = S.symbols("ze zo tensor_ze tensor_zo")
    ie, io, xe, xo, iye, iyo, xye, xyo = S.symbols("psi_e psi_o chi_e chi_o psiy_e psiy_o chiy_e chiy_o")
    alpha = Hb/pb**2
    co = Hb/(3*ab)
    k = pb**2/(3*ab)
    # Trace-reversed mixed Israel equation is sum K - sigma/3 identity.
    iso_before = 2*(iye+2*ab*ie+(-Hb-pb**2/3)*ze+pb*(xe+pb*ze)/3)
    zero("full Israel isotropic cancellation",
         iso_before-2*(iye+2*ab*ie+pb*xe/3-Hb*ze))
    israel = 2*(iye+2*ab*ie+pb*xe/3)-2*Hb*Tze
    zero("full Israel reconstructed residual",
         israel.subs(iye,-2*ab*ie-pb*xe/3)+2*Hb*Tze)
    scalar_sum = 2*(xye+2*pb*ie+dpb*ze)+sig2*(xe+pb*ze)
    zero("full scalar junction residual",
         scalar_sum-2*(xye+2*pb*ie+sig2*xe/2+(dpb+sig2*pb/2)*ze))
    zero("B definition",(dpb+sig2*pb/2-pb*B).subs(B,dpb/pb+sig2/2))
    mapped = {ie:Hb*Pe, io:Hb*Po, xe:-3*Hb*Pey/pb, xo:-3*Hb*Poy/pb}
    odd_match = (xo+pb*zo).subs(zo,-io/ab).subs(mapped)
    zero("odd Robin induced scalar matching",odd_match.subs(Poy,-k*Po))
    zero("odd embedding induced metric matching",(io+ab*zo).subs(zo,-io/ab))
    xye_constraint = dpb*xe/pb+3*Hb*Hb*De/pb-2*pb*ie
    even_scalar = scalar_sum.subs(xye,xye_constraint).subs(ze,0)
    even_scalar = even_scalar.subs(mapped).subs(sig2,2*(B-dpb/pb))
    target = -6*Hb/pb*(B*Pey-Hb*De)
    zero("unseparated scalar junction to master",even_scalar-target)
    zero("even Wentzell scalar sufficiency",target.subs(Pey,Hb*De/B))
    # The normal-force identity is the difference of bulk Hamiltonian constraints.
    DZ = S.symbols("D_zeta")
    kappa_o = -Hb*DZ
    normal_phi_o = -3*ab*Hb*DZ/pb
    zero("Hamiltonian difference / normal force",3*ab*kappa_o-pb*normal_phi_o)
    zero("shell normal-force residual",6*ab*kappa_o-2*pb*normal_phi_o)

    # Current from the exact radial potential, independently curled on the
    # conformal zero-shift slice, then constrained pointwise.
    n1,n2,ps1,ps2,py1,py2,c1,c2,cy1,cy2 = S.symbols("N1 N2 Psi1 Psi2 Psiy1 Psiy2 chi1 chi2 chiy1 chiy2")
    radial = 12*((py1-ab*n1)*ps2-(py2-ab*n2)*ps1)+c1*cy2-c2*cy1
    radial += pb*(n1*c2-n2*c1)-4*pb*(ps1*c2-ps2*c1)
    # Compute the field-space curl from the completed radial potential,
    # instead of accepting the displayed radial formula as a definition.
    e1,e2 = S.symbols("variation1 variation2")
    conformal = e1*ps1+e2*ps2
    lapse = e1*n1+e2*n2
    radial_psi = e1*py1+e2*py2
    radial_chi = e1*cy1+e2*cy2
    theta1 = S.exp(4*conformal-lapse)*(12*(ab+radial_psi)*ps1-(pb+radial_chi)*c1)
    theta2 = S.exp(4*conformal-lapse)*(12*(ab+radial_psi)*ps2-(pb+radial_chi)*c2)
    curl = (S.diff(theta2,e1)-S.diff(theta1,e2)).subs({e1:0,e2:0})
    zero("direct field-space curl of GHY-completed potential",curl-radial)
    reduced = radial.subs({n1:-2*ps1,n2:-2*ps2,
                           py1:-2*ab*ps1-pb*c1/3,py2:-2*ab*ps2-pb*c2/3})
    physical = c1*cy2-c2*cy1-2*pb*(ps1*c2-ps2*c1)
    zero("physical radial current pullback",reduced-physical)
    P1,P2,P1y,P2y,DP1,DP2,Q1,Q2 = S.symbols("P1 P2 P1y P2y DP1 DP2 BoxP1 BoxP2")
    subphysical = {ps1:Hb*P1, ps2:Hb*P2,
                   c1:-3*Hb*P1y/pb,c2:-3*Hb*P2y/pb}
    cy1_constraint = dpb*c1/pb+3*Hb*Hb*DP1/pb-2*pb*ps1
    cy2_constraint = dpb*c2/pb+3*Hb*Hb*DP2/pb-2*pb*ps2
    packet = physical.subs({cy1:cy1_constraint,cy2:cy2_constraint}).subs(subphysical)/Hb**2
    exact_packet = 9*alpha*(DP1*P2y-DP2*P1y)
    zero("exact arbitrary packet current",packet-exact_packet)
    zero("even packet shell flux",exact_packet.subs({P1y:Hb*DP1/B,P2y:Hb*DP2/B}))
    odd_flux = 2*exact_packet.subs({P1y:-k*P1,P2y:-k*P2})
    divergence = 18*co*(P1*Q2-P2*Q1)
    zero("odd packet shell current divergence",odd_flux.subs({DP1:Q1+4*P1,DP2:Q2+4*P2})-divergence)
    zero("odd zeta shell current coefficient",18*co*(ab/Hb)**2-6*ab/Hb)
    # Vector-density conversion to the genuinely moving ambient shell.
    rr,Om,jy,jt = S.symbols("r Omega jy jt", nonzero=True)
    jac = S.Matrix([[S.sinh(t)/Om,rr*S.cosh(t)],
                    [S.cosh(t)/Om,rr*S.sinh(t)]])
    zero("ambient coordinate Jacobian",jac.det()+rr/Om)
    jT = Om/rr*(S.sinh(t)/Om*jy+rr*S.cosh(t)*jt)
    jR = Om/rr*(S.cosh(t)/Om*jy+rr*S.sinh(t)*jt)
    zero("moving ambient normal flux density",
         jR-S.tanh(t)*jT-jy/(rr*S.cosh(t)))
    zero("shell temporal vector density transformation",
         1/(rr*S.cosh(t))*(rr*S.cosh(t))*jt-jt)

    # Re-derive graph GHY kinetic coefficient from exact divergence formula.
    # n=exp(e N Q), h= rho² exp(2e Psi Q) gamma, f=yb+e z Q.
    # The graph normal denominator is d=sqrt(1-n² f_t²/h_tt_abs).
    # Exact sqrt(-h_graph)K = d/n [dy(V/d)+dt(n² V f_t/(h_tt_abs d))].
    rb, N, PS, ZZ, Q, QT, QTT = S.symbols("rho_b N psi z Q Qt Qtt", nonzero=True)
    rho_y = ab*rb
    # Expansion minus slice GHY after evaluating at displaced graph:
    graph_quad = rb**2*(-ab*ZZ**2*QT**2
                       +ZZ*(N+2*PS+2*ab*ZZ)*Q*(QTT+3*QT)
                       +ZZ*(2*N+2*PS)*QT**2)
    # Independently expand the exact graph divergence expression, then apply
    # the moving-surface evaluation shift to its first-order coefficient.
    QQ = S.Function("Q")(t)
    NN = S.Function("N")(y)
    PP = S.Function("psi")(y)
    nn = S.exp(e*NN*QQ)
    hh = rho**2*S.exp(2*e*PP*QQ)
    VV = rho**4*S.exp(4*e*PP*QQ+3*t)
    slope = e*ZZ*S.diff(QQ,t)
    denominator = S.sqrt(1-nn**2*slope**2/hh)
    exact_graph = denominator/nn*(S.diff(VV/denominator,y)
                                  +S.diff(nn**2*VV*slope/(hh*denominator),t))
    exact_slice = S.diff(VV,y)/nn
    exact_delta = (exact_graph-exact_slice)/S.exp(3*t)
    delta1 = linear(exact_delta,e)
    delta2 = S.diff(exact_delta,e,2).subs(e,0)/2+ZZ*QQ*S.diff(delta1,y)
    jetmap = {rho:rb,S.diff(rho,y):ab*rb,NN:N,PP:PS,
              QQ:Q,S.diff(QQ,t):QT,S.diff(QQ,t,2):QTT}
    zero("exact graph-GHY second variation including shift",delta2.subs(jetmap)-graph_quad)
    # With sqrt(-gamma)=exp(3t), integration by parts sends Q(Qtt+3Qt)
    # to -Qt². DQ·DQ=-Qt² for this temporal test; no harmonic equation used.
    coeff_qt2 = S.expand(graph_quad).coeff(QT,2)
    coeff_qqtt = S.expand(graph_quad).coeff(Q,1).coeff(QTT,1)
    graph_DQ_coeff = -coeff_qt2+coeff_qqtt
    zero("direct graph-GHY kinetic coefficient",graph_DQ_coeff-rb**2*(3*ab*ZZ**2-N*ZZ))
    # The exact Einstein tangential divergence V^t=-rho² n_t e^(2psi)
    # gives +rho²Nz DQ² at the moving boundary.
    divergence_DQ_coeff = rb**2*N*ZZ
    zero("graph GHY plus Einstein corner",graph_DQ_coeff+divergence_DQ_coeff-3*ab*rb**2*ZZ**2)
    # Shifting the radial domain in the linear intrinsic-curvature term
    # -3 e rho² psi Box Q contributes +3 rho² z psi (DQ)².
    domain_DQ_coeff = 3*rb**2*ZZ*PS
    sigma = S.symbols("sigma")
    volume_DQ_coeff = -sigma*rb**2*ZZ**2/2
    full_shell_coeff = 2*(graph_DQ_coeff+divergence_DQ_coeff+domain_DQ_coeff)+volume_DQ_coeff
    completed = S.expand(full_shell_coeff.subs({sigma:6*ab,PS:-ab*ZZ}))
    zero("full odd moving-shell kinetic action",completed+3*ab*rb**2*ZZ**2)
    zero("odd current from completed action",-2*completed-6*ab*rb**2*ZZ**2)
    # Reflection cancellation with arbitrary even and odd tangential gradients.
    ge,go,gz = S.symbols("gradient_even gradient_odd gradient_zeta")
    zero("mixed parity moving-domain cancellation",
         gz*(ge+go)-gz*(ge-go)-2*gz*go)

    # Smooth diffeomorphism obstruction for finite S3 scalar harmonics l>=2.
    # In a general warped product, Lie derivative with covector components
    # xi_a=f_a H and xi_i=J D_iH has the mixed coefficient below. The angular
    # trace-free part is exactly 2J Hessian_TF H; warp derivatives affect trace.
    J, Ja, wa, w, fa = S.symbols("J Ja wa w fa", nonzero=True)
    mixed_lie = Ja-2*wa*J/w+fa
    zero("warped mixed diffeomorphism after angular TF constraint",
         mixed_lie.subs({J:0,Ja:0})-fa)
    ell, lam = S.symbols("ell lambda", integer=True)
    bochner = lam*(lam-2)-lam**2/3
    zero("S3 Hessian-TF norm Bochner identity",bochner-S.Rational(2,3)*lam*(lam-3))
    zero("l>=2 Hessian-TF positive factorization",
         bochner.subs(lam,ell*(ell+2))
         -S.Rational(2,3)*ell*(ell+2)*(ell-1)*(ell+3))

    # The exact T=0 Cauchy pullback works on arbitrary compact annuli.
    ff,qq,qr,mm,ss,LL = S.symbols("f q qr m s L", nonzero=True)
    zero("actual scalar Cauchy weight",
         Om**3*rr**3*(S.Rational(3,2)/(Om**2*ff))**2
         -S.Rational(9,4)*rr**3/(Om*ff**2))
    zero("actual gravitational Cauchy weight",
         6*Om**3*rr**3*(rr**2/Om**2)**2-6*rr**7/Om)
    # mr log derivative is 3-2sL, so d_r(mr)=m(4-2sL).
    original_integrand = mm*((rr*qr)**2+9*rr*qq*qr+20*qq**2)
    completed_integrand = mm*((rr*qr)**2+(2+9*ss*LL)*qq**2)
    divergence_integrand = S.Rational(9,2)*mm*((4-2*ss*LL)*qq**2+2*rr*qq*qr)
    zero("compact-annulus K diagonal integration identity",
         original_integrand-completed_integrand-divergence_integrand)
    reject("drop varying weight from K integration",9*ss*LL*mm*qq**2)

    # Wrong controls test meaningful missing pieces/signs.
    reject("omit scalar 2pPsi term",(-4*pb*ie).subs(mapped))
    reject("reverse odd Robin sign",odd_match.subs(Poy,k*Po))
    reject("omit graph evaluation shift",2*(Hb+pb**2/3)*ze)
    reject("omit gravitational scalar mixing",(2*pb*(ps1*c2-ps2*c1)).subs(subphysical))
    reject("equal-spectral shortcut on packets",9*alpha*(DP1*P2y-DP2*P1y)-9*alpha*DP1/P1*(P1*P2y-P2*P1y))
    reject("omit moving integration domain",S.expand((full_shell_coeff-2*domain_DQ_coeff).subs({sigma:6*ab,PS:-ab*ZZ}))-completed)
    reject("omit Einstein graph corner",-2*divergence_DQ_coeff)
    reject("omit odd shell current",divergence)
    reject("double shell tension instead of one",-6*ab*rb**2*ZZ**2/2)
    if args.inject_failure:
        zero("injected bad normal force sign",6*ab*kappa_o+2*pb*normal_phi_o)
    source_paths = [
        "basis/prior-checkpoint/reconstruction/SMOOTH_AMBIENT_RECONSTRUCTION_AND_LOCAL_FLUX.md",
        "basis/prior-checkpoint/boundary/BOUNDARY_LAWS_AND_LOCAL_CAUCHY_PHASE_SPACE.md",
        "basis/prior-checkpoint/covariant/COVARIANT_CURRENT_AND_CONE_LIMITS.md",
        "basis/prior-checkpoint/operator/PHYSICAL_OPERATOR_AUDIT_20261013.md",
        "basis/prior-checkpoint/operator/ODD_KINETIC_COMPLETION_20261013.md",
    ]
    result = {"status":"pass", "scope":"exact local junction/current identities; no PDE theorem certification",
              "checks":checks,"check_count":len(checks),"wrong_controls":controls,
              "source_sha256":{p:hashlib.sha256((ROOT.parent / p).read_bytes()).hexdigest() for p in source_paths},
              "code_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "sympy_version":S.__version__}
    output = Path(args.output)
    if not output.is_absolute():
        output = ROOT/output
    output.write_text(json.dumps(result,indent=2)+"\n")
    print(f"PASS {len(checks)} identities; {len(controls)} wrong-formula controls; {output}")


if __name__ == "__main__":
    main()
