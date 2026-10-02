#!/usr/bin/env python3
"""Pure exact algebra. No physical source, mode, root, or spectrum evaluation."""
import argparse
import hashlib
import json
import platform
from pathlib import Path
import sympy as s


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    eta, eps, k, x, r, H = s.symbols('eta eps k x r H', real=True)
    a, h, n, dx, T, v, Om, Q = [s.Function(z)(eta) for z in ('a','h','n','dx','T','v','Om','Q')]
    L = s.diff(a, eta)/a
    c = (3*h-n)/2
    N = a*(1+eps*n)
    A = a*(1+eps*h)
    F = a*(1+eps*h)**s.Rational(3,2)*(1+eps*n)**s.Rational(-1,2)
    full = N**2*(k**2/A**2+x+eps*dx)-s.diff(F,eta,2)/F
    derived = s.diff(full, eps).subs(eps, 0)
    expected = 2*k**2*(n-h)+2*a**2*x*n+a**2*dx-s.diff(c,eta,2)-2*L*s.diff(c,eta)
    checks = []
    def eq(name, lhs, rhs):
        residual = s.factor(s.simplify(lhs-rhs))
        if residual != 0:
            raise RuntimeError(f'{name}: nonzero residual {residual}')
        checks.append({'name':name, 'status':'PASS', 'exact_residual':'0'})
    eq('general_lapse_scale_canonical_frequency', derived, expected)
    conformal = expected.subs(n,h).doit()
    eq('conformal_gauge_frequency', conformal, a**2*dx+2*a**2*x*h-s.diff(h,eta,2)-2*L*s.diff(h,eta))
    de_sitter = conformal.subs(a,-1/(H*eta)).doit()
    eq('arbitrary_root_conformal_frequency', de_sitter, dx/(H**2*eta**2)+2*x*h/(H**2*eta**2)-s.diff(h,eta,2)+2*s.diff(h,eta)/eta)
    eq('exact_plane_wave_calibration_frequency', de_sitter.subs(x,2*H**2), dx/(H**2*eta**2)+4*h/eta**2-s.diff(h,eta,2)+2*s.diff(h,eta)/eta)
    Omega = k**2+a**2*x-s.diff(a,eta,2)/a
    gauge = expected.subs({h:-L*T,n:-s.diff(T,eta)-L*T,dx:0}).doit()
    gauge_expected = -T*s.diff(Omega,eta)-2*s.diff(T,eta)*Omega-s.diff(T,eta,3)/2
    eq('arbitrary_background_time_diffeomorphism_frequency', gauge, gauge_expected)
    vg = -T*s.diff(v,eta)+s.diff(T,eta)*v/2
    mode_residual = s.diff(vg,eta,2)+Om*vg
    mode_residual = mode_residual.subs(s.diff(v,eta,3),-s.diff(Om,eta)*v-Om*s.diff(v,eta))
    mode_residual = mode_residual.subs(s.diff(v,eta,2),-Om*v)
    eq('time_diffeomorphism_mode_solution',mode_residual,(T*s.diff(Om,eta)+2*s.diff(T,eta)*Om+s.diff(T,eta,3)/2)*v)
    cg=s.diff(T,eta)/2-L*T
    qg=-2*cg*Q-T*s.diff(a**2*Q,eta)/a**2+s.diff(T,eta)*Q
    eq('bare_or_finite_band_variance_covariance',qg,-T*s.diff(Q,eta))
    # General canonical Wronskian normalization gives the Green jump, not a
    # special Hankel/plane-wave state assumption.
    vv, vp, vc, vcp=s.symbols('vv vp vc vcp')
    eq('retarded_green_unit_jump',s.I*(vp*vc-vcp*vv),-s.I*(vv*vcp-vp*vc))
    eq('retarded_green_jump_from_unit_wronskian',(-s.I*(vv*vcp-vp*vc)).subs(vv*vcp-vp*vc,s.I),1)
    # Expand the general minimal FRW trace remainder at a constant de Sitter
    # root. The perturbation dH is the proper-Hubble invariant.
    dH=s.Function('dH')(eta)
    HH=H+eps*dH
    xx=x+eps*dx
    RR=6*(s.diff(HH,eta)+2*HH**2)
    boxR=-s.diff(RR,eta,2)-3*HH*s.diff(RR,eta)
    boxx=-s.diff(xx,eta,2)-3*HH*s.diff(xx,eta)
    anomaly=(xx-r-RR/6)**2/2-HH**2*(s.diff(HH,eta)+HH**2)/15+boxR/30-boxx/6
    dR=6*s.diff(dH,eta)+24*H*dH
    M=x-r-2*H**2
    anomaly_expected=M*dx-(M/6+H**2/90)*dR-(s.diff(dR,eta,2)+3*H*s.diff(dR,eta))/30+(s.diff(dx,eta,2)+3*H*s.diff(dx,eta))/6
    eq('proper_time_metric_and_mass_trace_remainder',s.diff(anomaly,eps).subs(eps,0),anomaly_expected)
    # Fixed-r variance subtraction: derive its complete variation from the
    # same order-two expression, without dropping metric derivatives of w.
    w, dw, Delta, dDelta=[s.Function(z)(eta) for z in ('w','dw','Delta','dDelta')]
    W=w+eps*dw
    DD=Delta+eps*dDelta
    U=DD/(2*W)-s.diff(W,eta,2)/(4*W**2)+3*s.diff(W,eta)**2/(8*W**3)
    U0=U.subs(eps,0)
    dU=dDelta/(2*w)-Delta*dw/(2*w**2)-s.diff(dw,eta,2)/(4*w**2)+s.diff(w,eta,2)*dw/(2*w**3)+3*s.diff(w,eta)*s.diff(dw,eta)/(4*w**3)-9*s.diff(w,eta)**2*dw/(8*w**4)
    eq('complete_metric_variation_order_two_subtraction',s.diff(U,eps).subs(eps,0),dU)
    SQ=(1/(2*W)-U/(2*W**2))/(a*(1+eps*h))**2
    SQ0=SQ.subs(eps,0)
    dSQ=-2*h*SQ0+(-dw/(2*w**2)-dU/(2*w**2)+U0*dw/w**3)/a**2
    eq('complete_metric_variation_variance_subtraction',s.diff(SQ,eps).subs(eps,0),dSQ)
    eq('fixed_geometry_subtraction_reduction',dSQ.subs({h:0,dw:0,dDelta:a**2*dx}).doit(),-dx/(4*w**3))
    # Linear Ward identity: N cancels between both proper-time derivatives.
    rho0,Q0=s.symbols('rho0 Q0')
    drho,dp=[s.Function(z)(eta) for z in ('drho','dp')]
    rho=rho0+eps*drho
    pressure=-rho0+eps*dp
    ward=s.diff(rho,eta)+3*s.diff(A,eta)/A*(rho+pressure)-s.diff(xx,eta)*(Q0+eps*Q)/2
    eq('linear_static_root_ward_has_no_metric_background_term',s.diff(ward,eps).subs(eps,0),s.diff(drho,eta)+3*L*(drho+dp)-Q0*s.diff(dx,eta)/2)
    # Radial brane bending cancels the coordinate variation of shell values.
    bulk_phi,wp,zeta,Y=s.symbols('bulk_phi wp zeta Y')
    eq('moving_shell_scalar_radial_gauge_invariant',(bulk_phi-wp*Y)+wp*(zeta+Y),bulk_phi+wp*zeta)
    bulk_h,qp=s.symbols('bulk_h qp')
    eq('moving_shell_induced_metric_radial_gauge_invariant',(bulk_h-qp*Y)+qp*(zeta+Y),bulk_h+qp*zeta)
    # Wrong-formula witnesses are exact rational-polynomial substitutes, not
    # additional quantum states or physical experiments.
    mutations=[]
    def reject(name, lhs, rhs, replacement):
        residual=s.simplify((lhs-rhs).subs(replacement).doit())
        if residual == 0 or residual.has(s.Derivative):
            raise RuntimeError(f'{name}: ineffective or unresolved mutation {residual}')
        mutations.append({'name':name,'status':'DETECTED','exact_nonzero_residual':str(residual)})
    simple={a:s.exp(eta),h:eta**2,n:eta**3,dx:eta, k:s.Integer(2),x:s.Integer(3)}
    reject('omit_lapse_gradient_channel',expected-2*k**2*(n-h),expected,simple)
    reject('omit_canonical_metric_second_derivative',expected+s.diff(c,eta,2),expected,simple)
    reject('replace_general_canonical_rescaling_by_h',2*k**2*(n-h)+2*a**2*x*n+a**2*dx-s.diff(h,eta,2)-2*L*s.diff(h,eta),expected,simple)
    reject('wrong_time_pullback_mode_contact',-T*s.diff(v,eta)-s.diff(T,eta)*v/2,vg,{T:eta**2,v:1+eta})
    reject('drop_variance_metric_prefactor',qg+2*cg*Q,-T*s.diff(Q,eta),{a:s.exp(eta),T:eta**2,Q:1+eta})
    reject('move_reference_with_metric_subtraction',dU.subs(dw,0).doit(),dU,{w:2+eta,dw:eta,Delta:eta,dDelta:1+eta})
    reject('drop_box_curvature_trace_term',anomaly_expected+(s.diff(dR,eta,2)+3*H*s.diff(dR,eta))/30,anomaly_expected,{dH:eta**3,dx:eta,H:1,x:2,r:2})
    reject('omit_brane_bending_scalar',bulk_phi,bulk_phi+wp*zeta,{bulk_phi:0,wp:2,zeta:3})
    result={'status':'PASS','scope':'Exact symbolic prerequisite identities and synthetic formula mutations only. No physical evaluation or registered response experiment.','identity_count':len(checks),'mutation_count':len(mutations),'identities':checks,'mutations':mutations,'python':platform.python_version(),'sympy':s.__version__,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:result[k] for k in ('status','identity_count','mutation_count','source_sha256')}))


if __name__=='__main__':
    main()
