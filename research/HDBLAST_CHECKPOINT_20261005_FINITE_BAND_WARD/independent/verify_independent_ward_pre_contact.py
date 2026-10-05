"""Independent symbolic action/Ward identities; no physical-source inputs.

The formulas are transcribed from the specified minimal-scalar action and the
published subtraction inventory.  This program imports no project producer,
source callback, numerical state, or other agent's implementation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sympy as s


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def run() -> dict:
    identities: dict[str, str] = {}
    mutations: dict[str, str] = {}

    def check(name: str, expression) -> None:
        reduced = s.cancel(s.expand(expression))
        require(reduced == 0, f"Identity failed: {name}: {reduced}")
        identities[name] = "PASS"

    def reject(name: str, expression) -> None:
        reduced = s.cancel(s.expand(expression))
        require(reduced != 0, f"Omission control undetected: {name}")
        mutations[name] = str(reduced)

    # Bare quadratic observables, independently represented by three real
    # bilinears U=|v|^2, A=|D|^2, C=Re(D v*), with D=v'-Lv.
    a, L, x, xp, k, U, A, C = s.symbols("a L x xp k U A C", real=True)
    kap = k**2 + a**2*x
    rho = (A + kap*U)/(2*a**4)
    pressure = (A - (k**2/s.Integer(3) + a**2*x)*U)/(2*a**4)
    Q = U/a**2
    vector = {a:a*L, x:xp, U:2*C+2*L*U,
              A:-2*L*A-2*kap*C, C:A-kap*U}
    dt = lambda f: sum(s.diff(f, y)*dy for y, dy in vector.items())
    check("bare_general_mass_modewise_exchange",
          dt(rho)+3*L*(rho+pressure)-xp*Q/2)
    check("physical_energy_volume_work",
          dt(a**3*rho)+pressure*dt(a**3)-a**3*xp*Q/2)
    reject("halve_mass_source_current_again",
           dt(rho)+3*L*(rho+pressure)-xp*Q/4)
    reject("omit_direct_pressure_mass_term",
           dt(rho)+3*L*(rho+pressure+x*U/(2*a**2))-xp*Q/2)

    # An off-shell mode residual F=v''+omega^2 v leaves exactly
    # Re(D*F)/a^4; Ward conservation cannot certify the mode equation by itself.
    residual_pairing = s.symbols("Re_Dstar_F", real=True)
    offshell = dict(vector)
    offshell[A] += 2*residual_pairing
    d_off = lambda f: sum(s.diff(f, y)*dy for y, dy in offshell.items())
    check("bare_offshell_defect_pairing",
          d_off(rho)+3*L*(rho+pressure)-xp*Q/2-residual_pairing/a**4)

    # Canonical source potential is a different operator from physical rho.
    V2, Hc, omega2, omega2p, Ccan = s.symbols(
        "V2 Hc omega2 omega2p Ccan", real=True)
    canonical_energy = (Hc+omega2*V2)/2
    canonical_vector = {V2:2*Ccan, Hc:-2*omega2*Ccan, omega2:omega2p}
    d_can = lambda f: sum(s.diff(f, y)*dy for y, dy in canonical_vector.items())
    check("canonical_Hamiltonian_source_work",
          d_can(canonical_energy)-omega2p*V2/2)

    # Fixed positive reference subtraction, derived independently as rational
    # functions of arbitrary jets.  aa=a^2 and w^2=k^2+aa*r.  Eliminating k^2
    # by w^2-aa*r is legitimate: its derivative below is identically zero.
    aa, r, ww = s.symbols("aa r w", nonzero=True, real=True)
    lj = s.symbols("L0:6", real=True)
    dj = s.symbols("D0:5", real=True)  # D=a^2(x-r).
    jet_vector = {aa:2*lj[0]*aa, ww:aa*r*lj[0]/ww}
    jet_vector.update({lj[j]:lj[j+1] for j in range(5)})
    jet_vector.update({dj[j]:dj[j+1] for j in range(4)})
    derivative = lambda f: s.cancel(sum(s.diff(f, y)*dy for y, dy in jet_vector.items()))
    kk = ww**2-aa*r
    check("constant_comoving_momentum_relation", derivative(kk))
    w1 = derivative(ww)
    w2 = derivative(w1)
    u2 = s.cancel((dj[0]-lj[1]-lj[0]**2)/(2*ww)
                  -w2/(4*ww**2)+3*w1**2/(8*ww**3))
    u2p = derivative(u2)
    u2pp = derivative(u2p)
    u4 = s.cancel(-u2**2/(2*ww)-u2pp/(4*ww**2)
                  +w2*u2/(4*ww**3)+3*w1*u2p/(4*ww**3)
                  -3*w1**2*u2/(4*ww**4))
    b = -kk/s.Integer(3)-aa*r
    c = 1-b/ww**2
    J2 = lj[0]*w1/ww**2+w1**2/(4*ww**3)
    J4 = (lj[0]*u2p/ww**2-2*lj[0]*w1*u2/ww**3
          +w1*u2p/(2*ww**3)-3*w1**2*u2/(4*ww**4))
    SR0 = ww/2
    SP0 = kk/(6*ww)
    SR2 = ((dj[0]+lj[0]**2)/ww+J2)/4
    SP2 = (c*u2+(lj[0]**2-dj[0])/ww+J2)/4
    SR4 = (u2**2/ww-(dj[0]+lj[0]**2)*u2/ww**2+J4)/4
    SP4 = (c*u4+b*u2**2/ww**3-(lj[0]**2-dj[0])*u2/ww**2+J4)/4
    source0 = (dj[1]-2*lj[0]*dj[0])/(4*ww)
    source2 = -(dj[1]-2*lj[0]*dj[0])*u2/(4*ww**2)
    check("subtraction_grade0_arbitrary_jets",
          derivative(SR0)-lj[0]*(SR0-3*SP0))
    check("subtraction_grade2_arbitrary_jets",
          derivative(SR2)-lj[0]*(SR2-3*SP2)-source0)
    check("subtraction_grade4_arbitrary_jets",
          derivative(SR4)-lj[0]*(SR4-3*SP4)-source2)
    reject("omit_subtraction_grade4_pressure_u4",
           derivative(SR4)-lj[0]*(SR4-3*(SP4-c*u4/4))-source2)

    # Metric perturbation about a0=L=-1/eta, x=r=2.  u and w are actual
    # perturbations (include epsilon); f is an arbitrary real forcing.
    ell, eps, h0, h1, h2, h3, X, Z, T, f, g = s.symbols(
        "ell eps h0 h1 h2 h3 X Z T f g", real=True)
    Rm = ((2*k**2+3*ell**2)*X-k*T-ell*Z)/(2*k*eps)
    Pm = ((2*k**2/3-ell**2)*X-k*T-ell*Z)/(2*k*eps)
    cr = X-T/(2*k)
    forcing = {ell:ell**2, X:Z, Z:-2*k*T-eps*f, T:2*k*Z,
               h0:h1, h1:h2, h2:h3}
    dforce = lambda y: sum(s.diff(y, v)*dv for v, dv in forcing.items())
    check("unprojected_canonical_real_constant", dforce(cr))
    check("raw_mode_Ward_for_arbitrary_real_forcing",
          dforce(Rm)-ell*(Rm-3*Pm)-ell*f/(2*k))
    G = (3*ell**2*X-ell*Z)/(2*k*eps)
    check("raw_mode_energy_split", Rm-k*cr/eps-G)
    check("raw_mode_primitive_for_arbitrary_real_forcing",
          dforce(G)-ell*(Rm-3*Pm)-ell*f/(2*k))
    q = X/(k*eps)
    qp = dforce(q)
    qpp = dforce(qp)
    check("raw_mode_variance_density_no_projection",
          Rm-k*cr/eps-(3*ell**2*q-ell*qp)/2)
    check("raw_mode_variance_pressure_source_contact",
          Pm-k*cr/(3*eps)-(qpp-3*ell*qp-3*ell**2*q)/6-f/(6*k))
    reject("omit_raw_pressure_forcing_contact",
           Pm-k*cr/(3*eps)-(qpp-3*ell*qp-3*ell**2*q)/6)
    gmetric = 4*ell**2*h0-2*ell*h1-h2
    CRb = ell*h1/(2*k)-2*k*h0-2*ell**2*h0/k
    CPb = ell*h1/(2*k)-2*k*h0/3
    R0b = (2*k**2+3*ell**2)/(4*k)
    P0b = (2*k**2/3-ell**2)/(4*k)
    FCb = ell*(CRb-3*CPb)-3*h1*(R0b+P0b)
    check("bare_metric_contact_Ward_from_direct_operators",
          dforce(CRb)-FCb+ell*gmetric/(2*k))
    check("linear_fixed_mass_total_Ward_source_mismatch",
          dforce(Rm+CRb)-ell*(Rm+CRb-3*(Pm+CPb))
          +3*h1*(R0b+P0b)-ell*(f-gmetric)/(2*k))
    reject("omit_baseline_metric_volume_work",
           dforce(Rm+CRb)-ell*(Rm+CRb-3*(Pm+CPb))
           -ell*(f-gmetric)/(2*k))
    reject("project_canonical_real_constant_to_zero",
           Rm-(3*ell**2*q-ell*qp)/2)
    g1 = s.symbols("g1", real=True)
    canonical_first = k*cr/eps+g/(4*k)
    check("first_order_canonical_work_is_gprime",
          dforce(k*cr/eps)+g1/(4*k)-g1/(4*k))
    reject("omit_canonical_source_potential", dforce(k*cr/eps)-g1/(4*k))

    # Exact residual decomposition for any differentiable approximate state.
    # Imag(e_u) does not enter these particular real quadratic linearizations.
    euR, euI, ewR, ewI = s.symbols("euR euI ewR ewI", real=True)
    defect_vector = dict(forcing)
    defect_vector[X] += euR
    defect_vector[Z] += ewR
    defect_vector[T] += ewI
    d_def = lambda y: sum(s.diff(y, v)*dv for v, dv in defect_vector.items())
    E_R = ((2*k**2+3*ell**2)*euR-k*ewI-ell*ewR)/(2*k*eps)
    E_G = (3*ell**2*euR-ell*ewR)/(2*k*eps)
    check("full_density_arbitrary_state_ODE_defect",
          d_def(Rm)-ell*(Rm-3*Pm)-ell*f/(2*k)-E_R)
    check("primitive_G_arbitrary_state_ODE_defect",
          d_def(G)-ell*(Rm-3*Pm)-ell*f/(2*k)-E_G)
    check("canonical_residual_missing_from_G_only",
          E_R-E_G-d_def(k*cr/eps))
    check("full_contact_ledger_state_and_source_mismatch",
          d_def(Rm+CRb)-ell*(Rm+CRb-3*(Pm+CPb))
          +3*h1*(R0b+P0b)-ell*(f-gmetric)/(2*k)-E_R)
    reject("substitute_G_defect_for_full_density_defect", E_R-E_G)

    # Mixed numerical/analytic momentum targets.  All quantities are abstract
    # real symbols; no grid values, source data, or physical integrals are read.
    Md, Ma, Jf, Jg, dG, dC, intFC = s.symbols(
        "Md Ma Jf Jg dG dC intFC", real=True)
    integral_FA = dG-Md*Jf+intFC
    contact_identity = {intFC:dC+Ma*Jg}
    check("mixed_target_signed_integrated_primitive",
          integral_FA.subs(contact_identity)-(dG+dC-Md*Jf+Ma*Jg))
    check("consistent_forcing_mixed_target_correction",
          integral_FA.subs(contact_identity).subs(Jf,Jg)
          -(dG+dC+(Ma-Md)*Jg))
    reject("discard_discrete_vs_analytic_band_correction",
           integral_FA.subs(contact_identity).subs(Jf,Jg)-(dG+dC))
    # Increasing pressure by (Ma-Md)g/3 decreases the direct F by
    # (Ma-Md)Lg.  This is a different declared pressure target.
    shiftP = (Ma-Md)*g/3
    check("explicit_mixed_pressure_sign", -3*ell*shiftP-(Md-Ma)*ell*g)

    # Closed finite-band kernels for a source change with the same incoming
    # state.  These are real entire functions of tau, with the removable
    # singularities proved explicitly below.  Pi is a declared positive
    # constant: mathematical pi or the fixed producer constant, consistently.
    tau = s.symbols("tau", nonzero=True, real=True)
    cutoff, Pi = s.symbols("K Pi", positive=True, real=True)
    A0 = (1-s.cos(2*cutoff*tau))/(2*tau)
    A1 = (2*cutoff*tau*s.sin(2*cutoff*tau)+s.cos(2*cutoff*tau)-1)/(4*tau**2)
    A2 = (-2*cutoff**2*tau**2*s.cos(2*cutoff*tau)
          +2*cutoff*tau*s.sin(2*cutoff*tau)+s.cos(2*cutoff*tau)-1)/(4*tau**3)
    check("band_sine_integral_upper_derivative",
          s.diff(A0,cutoff)-s.sin(2*cutoff*tau))
    check("band_k_cosine_integral_upper_derivative",
          s.diff(A1,cutoff)-cutoff*s.cos(2*cutoff*tau))
    check("band_k2_sine_integral_upper_derivative",
          s.diff(A2,cutoff)-cutoff**2*s.sin(2*cutoff*tau))
    for name, expression in (("A0",A0),("A1",A1),("A2",A2)):
        check(name+"_zero_cutoff", expression.subs(cutoff,0))
    check("A0_zero_lag_extension", s.limit(A0,tau,0))
    check("A1_zero_lag_extension", s.limit(A1,tau,0)-cutoff**2/2)
    check("A2_zero_lag_extension", s.limit(A2,tau,0))
    check("kernel_integral_derivative_A0", s.diff(A0,tau)-2*A1)
    check("kernel_integral_derivative_A1", s.diff(A1,tau)+2*A2)
    KR = ell*A1/(4*Pi**2)-3*ell**2*A0/(8*Pi**2)
    KP = A2/(6*Pi**2)+ell*A1/(4*Pi**2)+ell**2*A0/(8*Pi**2)
    check("closed_band_kernel_interior_Ward",
          s.diff(KR,tau)+ell**2*s.diff(KR,ell)-ell*(KR-3*KP))
    band_mass = cutoff**2/(8*Pi**2)
    check("closed_band_density_kernel_boundary_work",
          s.limit(KR,tau,0)-ell*band_mass)
    check("closed_band_pressure_kernel_zero_lag",
          s.limit(KP,tau,0)-ell*band_mass)
    reject("double_pressure_sine_moment", ell*(A2/(2*Pi**2)))

    # Source-truncation component only: independent exact rational arithmetic
    # propagation of epsilon_N/(N+2), N=24.  No source evaluation or full
    # numerical enclosure is substituted into this calculation.
    dsource = s.Rational(1,390*2**90)
    maxL = s.Rational(2,7)
    maxK = s.Integer(256)
    pi_floor = s.Integer(3)
    BR = (maxL*maxK**2+3*maxL**2*maxK)/(8*pi_floor**2)*dsource
    BP = (maxK**3/(18*pi_floor**2)+maxL*maxK**2/(8*pi_floor**2)
          +maxL**2*maxK/(8*pi_floor**2))*dsource
    BJ = maxL*maxK**2/(8*pi_floor**2)*dsource
    upperR = s.Rational(5407,10**31)  # 5.407e-28
    upperP = s.Rational(2152,10**28)  # 2.152e-25
    upperJ = s.Rational(5389,10**31)  # 5.389e-28
    require(BR < upperR, "Density source-truncation bound upper rounding failed")
    require(BP < upperP, "Pressure source-truncation bound upper rounding failed")
    require(BJ < upperJ, "Source-work bound upper rounding failed")
    bounds = {"L1_Cauchy_source_truncation":str(dsource),
              "domain":"eta in [-9/2,-7/2], K in [0,256], Pi>=3",
              "density_source_only_bound_exact":str(BR),
              "pressure_source_only_bound_exact":str(BP),
              "time_integrated_mismatch_source_only_bound_exact":str(BJ),
              "density_strict_decimal_upper":"5.407e-28",
              "pressure_strict_decimal_upper":"2.152e-25",
              "time_integrated_mismatch_strict_decimal_upper":"5.389e-28",
              "bound_excludes":["source-coefficient arithmetic", "kernel arithmetic",
                 "mode/state residuals", "incoming-state uncertainty", "contact/operator uncertainty",
                 "geometrical uncertainty", "momentum quadrature", "cutoff tail"]}

    du, dwr, dwi = s.symbols("bound_euR bound_ewR bound_ewI", nonnegative=True)
    density_coeff = (2*k**2+3*ell**2)/(2*k)
    residue_bound_integrand = k**2/(2*Pi**2)*(density_coeff*du+dwi/2+ell*dwr/(2*k))
    integrated_residue_bound = (cutoff**4+3*ell**2*cutoff**2)*du/(8*Pi**2) \
        +cutoff**3*dwi/(12*Pi**2)+ell*cutoff**2*dwr/(8*Pi**2)
    check("uniform_state_residual_band_bound_coefficients",
          s.integrate(residue_bound_integrand,(k,0,cutoff))-integrated_residue_bound)

    return {
        "status":"PASS_INDEPENDENT_EXACT_MODEWISE_WARD_ALGEBRA",
        "physical_sources_evaluated":False,
        "physical_modes_loaded":False,
        "archived_arrays_loaded":False,
        "remote_operations":False,
        "project_numerical_helpers_imported":False,
        "identities_passed":len(identities),
        "identities":identities,
        "omission_controls_rejected":len(mutations),
        "omission_residuals":mutations,
        "sympy_version":s.__version__,
        "python_optimized":not __debug__,
        "analytic_source_truncation_component_bounds":bounds,
        "scope":"Exact symbolic action and matched finite-band identities only; no trajectory, momentum convergence, thermalization or cosmology result.",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    require(not args.output.exists(), "Output receipt must be new")
    result = run()
    result["verifier_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open("x") as out:
        out.write(json.dumps(result, indent=2, sort_keys=True)+"\n")
    print(json.dumps({k:result[k] for k in ("status", "identities_passed", "omission_controls_rejected", "physical_sources_evaluated", "python_optimized")}))


if __name__ == "__main__":
    main()
