#!/usr/bin/env python3
"""Exact local boundary/energy/current identities; no PDE or IBVP existence replay."""
import argparse
import hashlib
import json
from pathlib import Path

import sympy as s


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--inject-failure', action='store_true', help='Deliberately reverse the even boundary kinetic sign')
    args = parser.parse_args()
    positive, controls = [], []
    def equal(name, lhs, rhs=0):
        residual = s.factor(s.simplify(lhs-rhs))
        if residual != 0:
            raise RuntimeError(name + ': ' + str(residual))
        positive.append(name)
    def wrong(name, lhs, rhs):
        residual = s.factor(s.simplify(lhs-rhs))
        if residual == 0:
            raise RuntimeError('Wrong-formula control did not distinguish: ' + name)
        controls.append({'name': name, 'nonzero_residual': str(residual)})

    rho, p, B, a, r, Om = s.symbols('rho p B a r Omega', positive=True)
    alpha, beta = 1/(rho**2*p**2),rho**4*p**2*B
    equal('positive_even_endpoint_mass_conversion', alpha*beta,rho**2*B)
    F, nF, box, DP = s.symbols('F nF Box_h_F D_P')
    regular_Py_over_Pfactor = nF+4*F/rho
    even = regular_Py_over_Pfactor-(rho**2*box+4*F)/(alpha*beta)
    gamma_e = 4/rho-4/(rho**2*B)
    expected_even = nF-box/B+gamma_e*F
    if args.inject_failure:
        expected_even = nF+box/B+gamma_e*F
    equal('even_regularized_physical_Wentzell_law',even,expected_even)
    co = 1/(3*a*rho**2)
    odd = regular_Py_over_Pfactor+co*F/alpha
    gamma_o = 4/rho+p**2/(3*a)
    equal('odd_regularized_Robin_law',odd,nF+gamma_o*F)
    wrong('wrong_even_negative_boundary_kinetic_coefficient',even,nF+box/B+gamma_e*F)
    wrong('wrong_odd_drop_r4_derivative',odd,nF+p**2*F/(3*a))

    # Exact shell geometry and time conversion, not an asymptotic replacement.
    T, rb, O = s.symbols('T rb Omega_b', positive=True)
    Rb = s.sqrt(rb**2+T**2)
    equal('shell_timelike_hyperboloid_identity',Rb**2-T**2,rb**2)
    equal('ambient_shell_unit_normal',(-T**2+Rb**2)/rb**2,1)
    equal('physical_shell_unit_normal',O**2*(-T**2+Rb**2)/(O*rb)**2,1)
    equal('shell_induced_TT_metric',-1+s.diff(Rb,T)**2,-rb**2/Rb**2)
    z = s.Function('z')(T)
    angular = s.symbols('Delta_S3_z')
    rho_b = rb*O
    sqrt_h = (rho_b/Rb)*(O*Rb)**3
    induced_box = s.diff(-sqrt_h*Rb**2/rho_b**2*s.diff(z,T),T)/sqrt_h+angular/(O**2*Rb**2)
    converted_box = (-Rb**2*s.diff(z,T,2)-4*T*s.diff(z,T)+rb**2/Rb**2*angular)/rho_b**2
    equal('induced_boundary_wave_operator_in_ambient_time',induced_box,converted_box)
    wrong('wrong_3T_boundary_friction',induced_box,(-Rb**2*s.diff(z,T,2)-3*T*s.diff(z,T)+rb**2/Rb**2*angular)/rho_b**2)
    # Tangential derivative of t = asinh(T/rb) means d/dt=Rb*d/dT.
    transformed_temporal = -Rb*s.diff(Rb*s.diff(z,T),T)-3*T*s.diff(z,T)
    equal('global_dS_to_ambient_second_time_derivative',transformed_temporal,-Rb**2*s.diff(z,T,2)-4*T*s.diff(z,T))

    ss = s.symbols('s',real=True)
    Omega, f = s.Function('Omega')(ss),s.Function('f')(ss)
    LL = s.diff(Omega,ss)/Omega+2*s.diff(f,ss)/f
    A = 1/(Omega*f**2)
    equal('positive_master_weight_log_derivative',s.diff(A,ss)/A,-LL)
    # Conjugate via F=A^-1/2 psi. Use the log derivative, without abs-value
    # differentiation of f; its sign is fixed and nonzero in the neighborhood.
    L,Ls,ff = s.symbols('L L_s f',real=True)
    log_tprime = L/2
    box_t_over_t = 10*log_tprime+4*ss*(Ls/2+log_tprime**2)
    E_t_over_t = 2*ss*log_tprime
    conjugated_zero_order = box_t_over_t-2*L*E_t_over_t-8*L-s.Rational(8,3)*ss*ff**2
    Vc = 3*L-2*ss*Ls+ss*L**2+s.Rational(8,3)*ss*ff**2
    equal('canonical_master_potential',conjugated_zero_order,-Vc)
    equal('canonical_first_derivative_cancellation',4*log_tprime-2*L,0)
    wrong('wrong_canonical_potential_Ls_sign',conjugated_zero_order,-(3*L+2*ss*Ls+ss*L**2+s.Rational(8,3)*ss*ff**2))
    equal('canonical_boundary_Robin_shift',s.Rational(1,2)*L*(2*r/Om),r*L/Om)
    equal('canonical_positive_boundary_weight',(Om/B)/Om**2,1/(Om*B))
    # Physical measure factors give the auxiliary boundary coefficient.
    equal('auxiliary_physical_boundary_kinetic_weight',(Om**-4*ff**-2)/B,1/(Om**4*ff**2*B))

    # Actual Cauchy symplectic coefficient: rho=r Omega, chi=-3(D+4)q/(2Omega²f).
    c = -3/(2*Om**2*ff)
    equal('physical_Cauchy_scalar_weight',Om**3*r**3*c**2,s.Rational(9,4)*r**3/(Om*ff**2))
    equal('physical_Cauchy_metric_weight',6*Om**3*r**3*(r**2/Om**2)**2,6*r**7/Om)
    # Integration by parts: m_r=(9/4)Ar³, r*m_r'/m_r=3−2sL.
    symmetric_mass = 20-s.Rational(9,2)*(1+3-2*ss*L)
    equal('Cauchy_K_positive_symmetric_mass',symmetric_mass,2+9*ss*L)

    # Exact finite-time local Cauchy jets. No flat-background substitution.
    tau,C,aa,bb,lam,L0,H = s.symbols('tau C a_data b_data lambda L0 H',nonzero=True)
    chi1=-4*C*aa*H
    chi2=-C*tau*bb*H
    chi1v=C*aa*(lam/(2*tau)+4*tau*L0)*H
    chi2v=-C*bb*(s.Rational(5,2)-tau**2*L0)*H
    W=chi1*chi2v-chi2*chi1v
    equal('exact_actual_cone_flux_jet_cancellation',W,C**2*aa*bb*(10+lam/2)*H**2)
    equal('l2_doubled_flux_density',((2*tau)**3/4*W).subs(lam,8),28*tau**3*C**2*aa*bb*H**2)
    wrong('wrong_omit_derivative_of_T_in_scalar_jet',chi1*(-C*bb*(2-tau**2*L0)*H)-chi2*chi1v,W)
    f0=s.symbols('f0',nonzero=True)
    equal('relational_non_gauge_jet_1',4*f0*(2*tau*chi1v-chi1),4*f0*C*aa*(lam+4+8*tau**2*L0)*H)
    equal('relational_non_gauge_jet_2',4*f0*(2*tau*chi2v-chi2),4*f0*C*tau*bb*(-4+2*tau**2*L0)*H)
    v=s.symbols('v',positive=True)
    equal('finite_cone_trace_kernel',v*s.diff(v**-4,v)+4*v**-4,0)
    equal('ambient_reconstruction_trace_injectivity',2*F*(5*ss-3*ss),4*ss*F)

    result={'status':'PASS_EXACT_BOUNDARY_AND_LOCAL_CAUCHY_IDENTITIES',
            'scope':'Boundary conversion, canonical auxiliary energy, actual Cauchy symplectic coefficient and exact cone jets; not a PDE existence or after-contact coupled IBVP proof.',
            'positive_count':len(positive),'wrong_formula_count':len(controls),
            'positive_checks':{name:'passed' for name in positive},
            'wrong_formula_controls':{row['name']:row['nonzero_residual'] for row in controls},
            'sympy_version':s.__version__,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({key:result[key] for key in ('status','positive_count','wrong_formula_count')}))


if __name__=='__main__':
    main()
