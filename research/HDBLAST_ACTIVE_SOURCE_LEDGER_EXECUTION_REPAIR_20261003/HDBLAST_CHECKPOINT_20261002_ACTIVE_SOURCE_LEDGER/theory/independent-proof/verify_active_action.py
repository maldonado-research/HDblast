#!/usr/bin/env python3
"""First-principles formal action/stress audit; no physical data or integrations.

The only inputs are symbolic scalar variables and the minimal scalar mode
action. Neither conservation nor this audit defines the pressure.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import sympy as s


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def audit():
    identities, mutations = [], []

    def zero(name, expression):
        residue = s.cancel(s.expand(expression))
        require(residue == 0, name + ': ' + str(residue))
        identities.append(name)

    def reject(name, expression):
        residue = s.cancel(s.expand(expression))
        require(residue != 0, 'Mutation invisible: ' + name)
        mutations.append({'name': name, 'nonzero_symbolic_residue': str(s.factor(residue))})

    # Nonlinear unrenormalized mode stress from the minimally coupled action.
    # M=m^2*a^2, C=a''/a, ell=a'/a. No state normalization is assumed.
    k, M = s.symbols('k M', positive=True)
    ell, C = s.symbols('ell C', real=True)
    ar, ai, br, bi = s.symbols('v_R v_I vp_R vp_I', real=True)
    kinetic = (br-ell*ar)**2+(bi-ell*ai)**2
    variance = ar**2+ai**2
    Rb = (kinetic+(k**2+M)*variance)/2
    Pb = (kinetic-(k**2/s.Integer(3)+M)*variance)/2
    def generic_D(q):
        return s.expand(br*s.diff(q, ar)+bi*s.diff(q, ai)
            -(k**2+M-C)*(ar*s.diff(q, br)+ai*s.diff(q, bi))
            +(C-ell**2)*s.diff(q, ell)+2*ell*M*s.diff(q, M))
    zero('minimal_action_nonlinear_bare_stress_Ward_arbitrary_state',
         generic_D(Rb)-ell*Rb+3*ell*Pb)
    zero('minimal_action_arbitrary_state_Wronskian_is_constant',
         generic_D(ar*bi-ai*br))
    reject('hold_canonical_M_constant_despite_fixed_physical_mass',
           generic_D(Rb)-2*ell*M*s.diff(Rb,M)-ell*Rb+3*ell*Pb)

    # a=L*(1+tau*h), a0=L=-1/eta, L'=L^2, physical m^2=2.
    # tau differentiates a response scaled by epsilon, while u,w contain
    # the actual mode variation. Four real phase-space coordinates remain
    # independent; no numerical Wronskian condition is imposed.
    L, epsilon = s.symbols('L epsilon', positive=True)
    X,Y,Z,T = s.symbols('u_R u_I w_R w_I', real=True)
    h = s.symbols('h0:6', real=True)
    tau = s.symbols('tau', real=True)
    u, w = X+s.I*Y, Z+s.I*T
    uc, wc = X-s.I*Y, Z-s.I*T
    g = 4*L**2*h[0]-2*L*h[1]-h[2]
    Mnew = 2*L**2*(1+2*tau*h[0])
    Cnew = 2*L**2+tau*(h[2]+2*L*h[1])
    zero('metric_directional_canonical_potential_g',
         s.diff(Mnew-Cnew,tau)-g)
    # v0'/v0=-ik and |v0|^2=1/(2k), the incoming BD branch at r=2.
    vfactor = 1+tau*u/epsilon
    vconjfactor = 1+tau*uc/epsilon
    Dfactor = -L-s.I*k+tau*((w-(L+s.I*k)*u)/epsilon-h[1])
    Dconjfactor = -L+s.I*k+tau*((wc-(L-s.I*k)*uc)/epsilon-h[1])
    abs_v = vfactor*vconjfactor/(2*k)
    abs_D = Dfactor*Dconjfactor/(2*k)
    # Direct rho and p variations, including a^-4. Pressure is the action
    # operator above, not a value recovered from a density derivative.
    rho_scaled = (1-4*tau*h[0])*(abs_D+(k**2+Mnew)*abs_v)/2
    p_scaled = (1-4*tau*h[0])*(abs_D-(k**2/s.Integer(3)+Mnew)*abs_v)/2
    R0 = s.expand(rho_scaled.subs(tau,0))
    P0 = s.expand(p_scaled.subs(tau,0))
    R = s.expand(s.diff(rho_scaled,tau).subs(tau,0))
    P = s.expand(s.diff(p_scaled,tau).subs(tau,0))
    Rmode = ((2*k**2+3*L**2)*X/k-T-L*Z/k)/(2*epsilon)
    Pmode = ((2*k**2/s.Integer(3)-L**2)*X/k-T-L*Z/k)/(2*epsilon)
    Roperator = L*h[1]/(2*k)
    Poperator = Roperator
    Rmass = L**2*h[0]/k
    Pmass = -Rmass
    Rscale = -4*h[0]*R0
    Pscale = -4*h[0]*P0
    Rc = Roperator+Rmass+Rscale
    Pc = Poperator+Pmass+Pscale
    zero('BD_bare_background_density',R0-(2*k**2+3*L**2)/(4*k))
    zero('BD_bare_background_pressure',P0-(2*k**2/s.Integer(3)-L**2)/(4*k))
    zero('direct_quadratic_Taylor_density_mode_operator_mass_scale',R-Rmode-Rc)
    zero('direct_quadratic_Taylor_pressure_mode_operator_mass_scale',P-Pmode-Pc)
    zero('density_contact_simplification',Rc-(L*h[1]/(2*k)-2*k*h[0]-2*L**2*h[0]/k))
    zero('pressure_contact_simplification',Pc-(L*h[1]/(2*k)-2*k*h[0]/3))

    # v=v0*(1+u): delta mode equation is u''-2ik u'=-epsilon*g.
    # Consequently u'=w, w'=2ik*w-epsilon*g.
    def D(q, forcing_sign=-1, Lprime=None):
        return s.expand((L**2 if Lprime is None else Lprime)*s.diff(q,L)
            +Z*s.diff(q,X)+T*s.diff(q,Y)
            +(-2*k*T+forcing_sign*epsilon*g)*s.diff(q,Z)
            +2*k*Z*s.diff(q,T)
            +sum(h[j+1]*s.diff(q,h[j]) for j in range(5)))
    zero('real_imag_forced_ODE',D(w)-2*s.I*k*w+epsilon*g)
    zero('amplitude_Wronskian_unprojected_constant',D(2*X-T/k))
    zero('finite_K_bare_background_Ward',D(R0)-L*R0+3*L*P0)
    zero('active_mode_only_Ward_force_term',D(Rmode)-L*Rmode+3*L*Pmode-L*g/(2*k))
    zero('active_local_contacts_cancel_force_term',D(Rc)-L*Rc+3*L*Pc
         +3*h[1]*(R0+P0)+L*g/(2*k))
    zero('full_active_bare_Ward_arbitrary_unprojected_amplitudes',
         D(R)-L*R+3*L*P+3*h[1]*(R0+P0))
    zero('forcing_first_derivative',D(g)-(8*L**3*h[0]+2*L**2*h[1]-2*L*h[2]-h[3]))
    zero('forcing_second_derivative',D(D(g))-(24*L**4*h[0]+12*L**3*h[1]-2*L*h[3]-h[4]))
    Fmode = L*(Rmode-3*Pmode)
    Fc = L*(Rc-3*Pc)-3*h[1]*(R0+P0)
    zero('direct_active_bare_Ward_integrand_contacts_expanded',
         Fc-(-2*k*h[1]-s.Rational(5,2)*L**2*h[1]/k-2*L**3*h[0]/k))
    zero('active_bare_integrand_direct_action_match',Fmode+Fc-D(R))
    euR,euI,ewR,ewI=s.symbols('defect_u_R defect_u_I defect_w_R defect_w_I',real=True)
    density_defect=((2*k**2+3*L**2)*euR/k-ewI-L*ewR/k)/(2*epsilon)
    offshell_D_R=(D(R)+euR*s.diff(R,X)+euI*s.diff(R,Y)
                  +ewR*s.diff(R,Z)+ewI*s.diff(R,T))
    zero('full_bare_Ward_projects_continuous_ODE_defects_exactly',
         offshell_D_R-L*R+3*L*P+3*h[1]*(R0+P0)-density_defect)

    # Integrating factor coordinates follow from the ODE, not a stress Ward
    # relation. They expose the forcing needed for an independent ledger.
    d = w/(2*s.I*k)
    c = u-d
    zero('variation_of_constants_d_ODE',D(d)-2*s.I*k*d-epsilon*s.I*g/(2*k))
    zero('variation_of_constants_c_ODE',D(c)+epsilon*s.I*g/(2*k))
    creal = X-T/(2*k)
    zero('real_constant_unprojected_mode_coefficient',D(creal))
    co, si = s.symbols('cos_phase sin_phase',real=True)
    def Dphase(q):
        return s.expand(D(q)-2*k*si*s.diff(q,co)+2*k*co*s.diff(q,si))
    rotated_w = w*(co-s.I*si)
    zero('rotating_w_ODE',Dphase(rotated_w)+epsilon*g*(co-s.I*si))

    # Renormalized Ward is the difference of two separately action-derived
    # operator identities. This is a conditional formal subtraction check,
    # not a calculation of actual WKB coefficients or integrals.
    NR,NP,nr,np,NRp,NPp,nrp,npp=s.symbols('NR NP delta_NR delta_NP NRprime NPprime delta_NRprime delta_NPprime')
    SR,SP=nr-4*h[0]*NR,np-4*h[0]*NP
    SRp=nrp-4*h[1]*NR-4*h[0]*NRp
    unscaled_sub_ward=nrp-L*nr+3*L*np-h[1]*NR+3*h[1]*NP
    baseline_sub_ward=NRp-L*NR+3*L*NP
    zero('physical_subtraction_metric_Ward_from_separate_operator_premises',
         SRp-L*SR+3*L*SP+3*h[1]*(NR+NP)
         -unscaled_sub_ward+4*h[0]*baseline_sub_ward)
    zero('renormalized_Ward_preserves_finite_K_baseline_difference',
         (D(R)-L*R+3*L*P+3*h[1]*(R0+P0))
         -(SRp-L*SR+3*L*SP+3*h[1]*(NR+NP))
         -(D(R)-SRp-L*(R-SR)+3*L*(P-SP)
           +3*h[1]*(R0-NR+P0-NP)))

    def ward(r,p,baseline=True,forcing_sign=-1,Lprime=None):
        return D(r,forcing_sign,Lprime)-L*r+3*L*p+(3*h[1]*(R0+P0) if baseline else 0)
    reject('omit_finite_K_baseline_contact',ward(R,P,baseline=False))
    reject('reverse_forcing_sign',ward(R,P,forcing_sign=1))
    reject('omit_operator_contact_from_both_stresses',ward(R-Roperator,P-Poperator))
    reject('omit_fixed_mass_contact_from_both_stresses',ward(R-Rmass,P-Pmass))
    reject('omit_a_minus_four_scale_contact_from_both_stresses',ward(R-Rscale,P-Pscale))
    reject('omit_metric_contact_only_from_pressure',ward(R,P-Pc))
    reject('incorrectly_treat_L_as_constant',ward(R,P,Lprime=s.Integer(0)))
    reject('force_Wronskian_zero_in_direct_density',R.subs(X,T/(2*k))-R)
    reject('use_source_free_rotating_w_ODE_inside_source',Dphase(rotated_w))
    reject('omit_mode_only_force_term',D(Rmode)-L*Rmode+3*L*Pmode)
    reject('pretend_arbitrary_ODE_defects_imply_Ward_accuracy',density_defect)

    formulas={
      'convention':'Per mode; multiply by k^2/(2*pi^2) and fixed momentum weights. R=a0^4*delta_rho/epsilon, P=a0^4*delta_p/epsilon; u,w retain epsilon.',
      'action_canonical_mode_ODE':"v''+(k^2+a^2*m^2-a''/a)*v=0",
      'g':str(g),'R0':str(R0),'P0':str(P0),
      'Rmode':str(Rmode),'Pmode':str(Pmode),'Rcontact':str(s.factor(Rc)),
      'Pcontact':str(s.factor(Pc)),'Fcontact':str(s.factor(Fc)),
      'Fmode':str(s.factor(Fmode)),
      'full_bare_Ward':"R'=L*(R-3*P)-3*h1*(R0+P0)",
      'mode_Ward':"Rmode'=L*(Rmode-3*Pmode)+L*g/(2*k)",
      'contact_Ward':"Rcontact'=L*(Rcontact-3*Pcontact)-3*h1*(R0+P0)-L*g/(2*k)",
      'd_coordinate':"d=w/(2*i*k); d'=2*i*k*d+i*epsilon*g/(2*k)",
      'c_coordinate':"c=u-w/(2*i*k); c'=-i*epsilon*g/(2*k); Re(c)'=0",
      'rotating_w':"[w*exp(-2*i*k*(t-a))]'=-epsilon*g*exp(-2*i*k*(t-a))",
      'ODE_defect_density_projection':str(density_defect),
      'ODE_defect_definition':"defect_u=u'-w; defect_w=w'-2*i*k*w+epsilon*g",
    }
    result={'status':'PASS_PURE_ACTION_IDENTITIES','identity_count':len(identities),
            'identities':identities,'mutation_count':len(mutations),'mutations_rejected':mutations,
            'physical_evaluations':0,'physical_arrays_loaded':0,
            'physical_source_integrals_evaluated':0,'numeric_profile_values_evaluated':0,
            'scope':'Exact symbolic action/ODE identities. No physical archive, trajectory, source value, or quadrature is read or evaluated. Renormalized subtraction transfer is conditional on its separately proved operator premises.'}
    return result,formulas


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--formulas',type=Path)
    args=parser.parse_args()
    require(not args.output.exists(),'Refusing to overwrite a proof receipt')
    require(not args.formulas or not args.formulas.exists(),'Refusing to overwrite formulas')
    result,formulas=audit()
    result.update(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  sympy_version=s.__version__,python_optimization=sys.flags.optimize)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    if args.formulas:
        args.formulas.parent.mkdir(parents=True,exist_ok=True)
        args.formulas.write_text(json.dumps(formulas,indent=2,sort_keys=True)+'\n')
    print(json.dumps({key:result[key] for key in ('status','identity_count','mutation_count','physical_evaluations','physical_arrays_loaded','python_optimization')}))


if __name__=='__main__':
    main()
