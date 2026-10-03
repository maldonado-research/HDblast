#!/usr/bin/env python3
"""Independent variational/algebra checks. No physical field or shell evolution.

Uses a lapse-dependent FRW action before gauge fixing, and independently checks
the doubled-action junction, Codazzi, shell flux and Weyl ledgers. Explicit
exceptions survive python -O. All positive and negative residuals are retained.
"""
import argparse
import hashlib
import json
import platform
from pathlib import Path

import sympy as s


def verify():
    rows = []

    def check(name, expression, negative=False, note=""):
        residual = s.simplify(s.expand(expression))
        passed = (residual != 0) if negative else (residual == 0)
        rows.append(dict(name=name, negative_control=negative,
                         passed=bool(passed), residual=str(residual), note=note))
        if not passed:
            raise RuntimeError(f"{name}: residual {residual}")

    t = s.symbols("t", real=True)
    a, N, x = (s.Function(q)(t) for q in ("a", "N", "x"))
    v0, v1, v2, v3, f0, f1, f2, alpha = s.symbols("v0 v1 v2 v3 f0 f1 f2 alpha")
    V = v0 + v1*x + v2*x**2/2 + v3*x**3/6
    F = f0 + f1*x + f2*x**2/2
    R = 6*(s.diff(a,t,2)/(a*N**2) + s.diff(a,t)**2/(a**2*N**2)
           - s.diff(a,t)*s.diff(N,t)/(a*N**3))
    L = N*a**3*(-V + F*R + alpha*R**2)

    def EL(variable, order):
        return sum((-1)**j*s.diff(s.diff(L, s.diff(variable,t,j)),t,j)
                   for j in range(order+1))

    gauge = {N: 1, **{s.diff(N,t,j): 0 for j in range(1,5)}}
    H, d, e, f, z, zz, xx, aa = s.symbols("H d e f z zz xx aa")
    jets = {a:aa, s.diff(a,t):aa*H, s.diff(a,t,2):aa*(H**2+d),
            s.diff(a,t,3):aa*(H**3+3*H*d+e),
            s.diff(a,t,4):aa*(H**4+6*H**2*d+3*d**2+4*H*e+f),
            x:xx, s.diff(x,t):z, s.diff(x,t,2):zz}

    def cosmic(expr):
        return s.simplify(expr.subs(gauge).subs(jets))

    rho_var = cosmic(-EL(N,1)/a**3)
    p_var = cosmic(EL(a,2)/(3*a**2))
    S_var = cosmic(-EL(x,0)/a**3)
    Vj, Fj = V.subs(x,xx), F.subs(x,xx)
    Vx, Fx = s.diff(Vj,xx), s.diff(Fj,xx)
    Fdot, Fddot = Fx*z, s.diff(Fx,xx)*z**2+Fx*zz
    Rj = 6*(d+2*H**2)
    rho_expected = Vj-6*Fj*H**2-6*H*Fdot+alpha*(36*d**2-216*H**2*d-72*H*e)
    p_expected = (-Vj+2*Fj*(2*d+3*H**2)+2*Fddot+4*H*Fdot
                  +alpha*(108*d**2+216*H**2*d+144*H*e+24*f))
    S_expected = Vx-Fx*Rj
    check("lapse_variation_density", rho_var-rho_expected,
          note="N varied before fixing N=1; no Ward equation used to infer pressure")
    check("scale_factor_variation_pressure", p_var-p_expected)
    check("mass_variation_current", S_var-S_expected)
    D = lambda expr: (s.diff(expr,H)*d+s.diff(expr,d)*e+s.diff(expr,e)*f
                      +s.diff(expr,xx)*z+s.diff(expr,z)*zz)
    check("nonlinear_V_F_and_R2_Ward", D(rho_expected)+3*H*(rho_expected+p_expected)-z*S_expected)
    check("consistently_doubled_local_action_also_conserves", D(2*rho_expected)+3*H*(2*rho_expected+2*p_expected)-z*2*S_expected,
          note="Intentional positive identity: conservation cannot select the action coefficient")
    check("nonlinear_F_pressure_requires_Fxx", -6*H*f2*z**2, True,
          note="Residual if the 2 F_xx xdot^2 pressure term is omitted")
    check("R2_pressure_requires_H_triple_dot", -72*alpha*H*f, True)
    check("local_current_requires_curvature_piece", -z*Fx*Rj, True)
    check("regular_at_H_zero", p_var.subs(H,0)-p_expected.subs(H,0))

    phi, phistar, G, m0, Q = s.symbols("phi phi_star G m0 Q")
    coupling = m0**2+G**2*(phi-phistar)**2
    Jphi = s.diff(coupling,phi)*Q/2
    check("mass_coupling_chain_rule", Jphi-G**2*(phi-phistar)*Q)
    check("half_factor_omission_detected", s.diff(coupling,phi)*Q-Jphi, True)

    kap2, sig, sig1, rho, p, J, vv = s.symbols("kappa5sq sigma sigma_phi rho p J v", nonzero=True)
    lam = sig/kap2
    shell_stress = s.diag(-lam-rho, -lam+p, -lam+p, -lam+p)
    traceS = s.trace(shell_stress)
    K = kap2/2*(shell_stress-s.eye(4)*traceS/3)
    ks = (sig+kap2*rho)/6
    k0 = (sig-kap2*(2*rho+3*p))/6
    w = -(sig1+kap2*J)/2
    check("trace_reversed_Israel_spatial", K[1,1]-ks)
    check("trace_reversed_Israel_time", K[0,0]-k0)
    check("doubled_scalar_boundary_variation", -2*w/kap2-sig1/kap2-J)
    check("single_copy_scalar_variation_fails", -w/kap2-sig1/kap2-J, True)
    check("double_shell_count_fails", -2*w/kap2-2*sig1/kap2-2*J, True)
    check("tension_only_time_junction_fails", sig/6-k0, True)

    HH, vdot, Jdot, sig2, Ud, Weyl = s.symbols("H vdot Jdot sigma_phiphi U_phi W")
    rhodot = -3*HH*(rho+p)+J*vv
    ksdot = (sig1*vv+kap2*rhodot)/6
    wdot = -(sig2*vv+kap2*Jdot)/2
    check("Codazzi_all_sources", ksdot+HH*(ks-k0)+w*vv/3)
    shell_gain = sig1*vv/kap2+rhodot+3*HH*(rho+p)
    check("paired_tension_matter_bulk_flux", shell_gain+2*w*vv/kap2)
    check("wrong_flux_orientation_detected", shell_gain-2*w*vv/kap2, True)
    check("omitted_scalar_source_detected", shell_gain+2*(-sig1/2)*vv/kap2, True)
    check("wrong_scalar_source_sign_detected", shell_gain+2*(-(sig1-kap2*J)/2)*vv/kap2, True)
    rdot, Bdot = s.symbols("rhodot Bdot")
    ks_dot_unpaired = (sig1*vv+kap2*rdot)/6
    momentum_shell = (-3*(Bdot*ks+ks_dot_unpaired)-3*HH*ks
                      +3*HH*k0+3*ks*Bdot-vv*w)
    check("shell_momentum_constraint_equals_Ward_residual", momentum_shell+kap2/2*(rdot+3*HH*(rho+p)-J*vv))

    # Einstein-scalar orthonormal stress: n is spacelike, u future timelike.
    U = s.symbols("U")
    metric = s.diag(-1,1,1,1,1)
    grad = s.Matrix([vv,0,0,0,w])
    Tbulk = (grad*grad.T-metric*((-vv**2+w**2)/2+U))/kap2
    check("bulk_scalar_flux_normalization", Tbulk[4,0]-w*vv/kap2)
    traceT = s.trace(metric*Tbulk)
    F00 = s.Rational(2,3)*kap2*(Tbulk[0,0]-(Tbulk[4,4]-traceT/4))
    Fii = s.Rational(2,3)*kap2*(Tbulk[1,1]+Tbulk[4,4]-traceT/4)
    check("Gauss_scalar_projection_time", F00-(vv**2/4-w**2/4+U/2))
    check("Gauss_scalar_projection_space", Fii-(5*vv**2/12+w**2/4-U/2))
    Hsq = ks**2+vv**2/12-w**2/12+U/6+Weyl
    Hdot = ks*(k0-ks)-vv**2/3-2*Weyl
    check("Gauss_spatial_projection_Raychaudhuri", -(2*Hdot+3*Hsq)-(Fii-ks**2-2*ks*k0+Weyl))
    Wdot_from_gauss = 2*HH*Hdot-2*ks*ksdot-vv*vdot/6+w*wdot/6-Ud*vv/6
    Wbalance = (4*ks*w*vv-4*HH*vv**2-vv*vdot+w*wdot-Ud*vv)/6
    check("paired_Weyl_balance", Wdot_from_gauss+4*HH*Weyl-Wbalance)
    check("wrong_Weyl_redshift_detected", Wdot_from_gauss+3*HH*Weyl-Wbalance, True)

    Vlocal, Flocal, E00, Eii, D00, Dii = s.symbols("Vloc Floc E00 Eii D00 Dii")
    rho_l, p_l = Vlocal-2*D00-2*alpha*E00, -Vlocal-2*Dii-2*alpha*Eii
    # Tensor equality: moving F and R2 to the left and V into sigma
    # must remove the same pieces from the metric right-hand source.
    k_alg = (sig+kap2*(rho+rho_l))/6
    check("potential_absorption_spatial", k_alg-(sig+kap2*Vlocal+kap2*(rho-2*D00-2*alpha*E00))/6)
    check("double_count_local_density_detected", (sig+kap2*Vlocal+kap2*(rho+rho_l))/6-k_alg, True)

    # Coordinate T is not shell conformal time: d eta/dT = exp(B-A).
    AT, BT, ATT, kk, BB, AA, massx = s.symbols("A_T B_T A_TT k B A x")
    Ac, Bc, uc = (s.Function(q)(t) for q in ("Achart", "Bchart", "u"))
    eta_derivative = lambda expr: s.exp(Ac-Bc)*s.diff(expr,t)
    clock_jets = {Ac:AA, Bc:BB, s.diff(Ac,t):AT,
                  s.diff(Bc,t):BT, s.diff(Ac,t,2):ATT}
    transformed = eta_derivative(eta_derivative(uc))/s.exp(2*(Ac-Bc))
    check("mode_clock_first_derivative", transformed-s.diff(uc,t,2)-(s.diff(Ac,t)-s.diff(Bc,t))*s.diff(uc,t))
    curvature = (eta_derivative(eta_derivative(s.exp(Ac)))/s.exp(Ac)).subs(clock_jets)
    mode_coefficient = s.exp(2*(BB-AA))*(kk**2+s.exp(2*AA)*massx-curvature)
    check("mode_clock_frequency", mode_coefficient-(s.exp(2*(BB-AA))*kk**2+s.exp(2*BB)*massx-ATT-2*AT**2+AT*BT))
    check("using_T_as_eta_is_wrong", (AT-BT), True,
          note="Generic inherited chart has nonconstant d eta/dT")
    iprod = s.symbols("jprod_T")
    check("inherited_scalar_velocity_helper_missing_piece", -s.exp(BB)*iprod/2, True,
          note="B1 omission in gFt helper only; current constraints consume gAt, not gFt")

    return dict(status="PASS", scope="Exact algebra and lapse variation only; no coupled experiment",
                positive_checks=sum(not q['negative_control'] for q in rows),
                negative_controls=sum(q['negative_control'] for q in rows),
                python=platform.python_version(), sympy=s.__version__, checks=rows,
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = verify()
    blob = json.dumps(result, indent=2, sort_keys=True)+"\n"
    if args.output:
        args.output.write_text(blob)
    print(json.dumps({q:result[q] for q in ("status","scope","positive_checks","negative_controls","source_sha256")},sort_keys=True))
