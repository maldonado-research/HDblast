#!/usr/bin/env python3
"""Exact manufactured controls, independent of producer code and data.

The explicit polynomial frequency-window packet is C^2 after zero extension,
not C-infinity. Its finite-order decay is a separate manufactured control of
the mechanism; it is not substituted for the theorem's all-N smoothness proof.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sympy as sp

checks=[]
def require(condition,detail):
    if not condition:raise RuntimeError(detail)
def exact(name,value):
    reduced=sp.simplify(value)
    require(reduced==0,(name,str(reduced)))
    checks.append({'name':name,'kind':'EXACT_SYMBOLIC','residual':'0'})
def negative(name,value):
    reduced=sp.simplify(value)
    require(reduced!=0,(name,'negative control lost sensitivity'))
    checks.append({'name':name,'kind':'NEGATIVE_CONTROL','nonzero_witness':str(reduced)})

p,c,s,omega=sp.symbols('p c s omega',positive=True)
g=sp.symbols('g',real=True,nonzero=True)
D=sp.symbols('D',real=True)
matrix=sp.Matrix([[sp.I*p-c,g],[-g,D]])
forcing=sp.Matrix([sp.I*p+c,g])
solution=matrix.inv()*forcing
R,Q=map(sp.simplify,solution)
den=sp.simplify(matrix.det())
exact('uneliminated_matrix_solution_boundary',((matrix*sp.Matrix([R,Q])-forcing)[0]))
exact('uneliminated_matrix_solution_brane',((matrix*sp.Matrix([R,Q])-forcing)[1]))
exact('reflection_negative_conjugate_denominator',R+sp.conjugate(den)/den)
exact('strict_denominator_real_imaginary_decomposition',sp.re(den)**2+sp.im(den)**2-((g*g-c*D)**2+p*p*D*D))
exact('g_sign_reflection_invariant',R.subs(g,-g)-R)
exact('g_sign_brane_amplitude_reverses',Q.subs(g,-g)+Q)

# Derive the flux from the actual quadratic current for arbitrary reflection,
# rather than defining net flux by its expected incoming/outgoing difference.
u,v=sp.symbols('u v',real=True)
r=u+sp.I*v
field=1+r;normal_derivative=sp.I*p*(r-1)
direct_flux=-sp.re((-sp.I*omega*field)*sp.conjugate(normal_derivative))/2
exact('direct_Noether_current_interference',direct_flux-omega*p*(u*u+v*v-1)/2)
exact('full_scattering_flux_zero',direct_flux.subs({u:sp.re(R),v:sp.im(R)}))
phase_wrong=-R
negative('unit_modulus_alone_does_not_validate_boundary',
         (sp.I*p*(phase_wrong-1)-c*(1+phase_wrong)+g*Q).subs({p:2,c:3,g:1,D:1}))
exact('wrong_phase_still_has_unit_modulus',phase_wrong*sp.conjugate(phase_wrong)-1)

exact('D_zero_matrix_has_nonzero_determinant',den.subs(D,0)-g*g)
exact('D_zero_reflection',R.subs(D,0)+1)
exact('D_zero_brane_response',Q.subs(D,0)-2*sp.I*p/g)
zero_coupling_degeneracy=matrix.subs({g:0,D:0})
require(zero_coupling_degeneracy.rank()==1,'exact g=D=0 must leave one free oscillator amplitude')
checks.append({'name':'exact_g_zero_D_zero_is_underdetermined','kind':'EXACT_RANK','rank':1})
G=1/(c-sp.I*p)
bareR=(sp.I*p+c)/(sp.I*p-c)
chi=1/(D-g*g*G)
exact('Robin_input_susceptibility_equivalence',Q-chi*g*(1+bareR))
exact('outgoing_amplitude_superposition',R-bareR-g*G*Q)

# Orthogonality is an extended bulk+brane Hilbert-space statement. Its
# bulk-only term is nonzero, so omitting q cannot establish zero occupation.
bound_relation=g*g/(c+s)-s*s-p*p
bulk_pair=g/(c+s)*(1/(s+sp.I*p)+R/(s-sp.I*p))
exact('full_bound_continuum_orthogonality',(Q+bulk_pair).subs(D,bound_relation))
negative('bulk_only_pairing_not_orthogonal',bulk_pair.subs(D,bound_relation).subs({p:2,c:3,g:1,s:1}))
exact('velocity_projection_also_zero',(-sp.I*omega*(Q+bulk_pair)).subs(D,bound_relation))

# A smooth boundary-compatible datum with q=qdot=0 can still contain the
# coupled bound oscillator. M=2,c=3,g=1,m^2=13/4 has s=1,Z=32/33.
y=sp.symbols('y',nonnegative=True)
Z=sp.Rational(32,33)
bulk_initial=(1+4*y)*sp.exp(-y)
bound_bulk=sp.sqrt(Z)*sp.exp(-y)/4
exact('zero_initial_q_profile_satisfies_Robin_boundary',sp.diff(bulk_initial,y).subs(y,0)-3*bulk_initial.subs(y,0))
overlap=sp.integrate(bound_bulk*bulk_initial,(y,0,sp.oo))
exact('nonzero_bound_overlap_of_zero_initial_q',overlap-3*sp.sqrt(Z)/8)
negative('zero_q_and_velocity_do_not_imply_zero_bound_projection',overlap)
exact('persistent_q_component_from_preexisting_projection',sp.sqrt(Z)*overlap-sp.Rational(4,11))

# Explicit compact frequency-window packet, derived independently using
# repeated exact endpoint integration. This demonstrates finite-order decay
# and a nonzero transient, with a declared C^2 rather than C-infinity window.
w=sp.symbols('w',real=True)
t=sp.symbols('t',real=True,nonzero=True)
B=(w-3)**3*(4-w)**3
for order in range(3):
    exact('packet_left_endpoint_derivative_'+str(order),sp.diff(B,w,order).subs(w,3))
    exact('packet_right_endpoint_derivative_'+str(order),sp.diff(B,w,order).subs(w,4))
negative('packet_finite_smoothness_boundary',sp.diff(B,w,3).subs(w,3))
def exact_window_transform(poly):
    return sp.expand(sum(-(sp.diff(poly,w,j).subs(w,4)*sp.exp(-4*sp.I*t)
                           -sp.diff(poly,w,j).subs(w,3)*sp.exp(-3*sp.I*t))
                         /(sp.I*t)**(j+1)
                         for j in range(sp.degree(poly,w)+1)))
qt=exact_window_transform(B)
qdt=exact_window_transform(-sp.I*w*B)
exact('closed_packet_transform_time_derivative',sp.diff(qt,t)-qdt)
exact('nonzero_packet_transient_at_t_zero',sp.integrate(B,(w,3,4))-sp.Rational(1,140))
exact('packet_zero_time_frequency_moment',sp.integrate(w*B,(w,3,4))-sp.Rational(1,40))
exact('packet_closed_form_removable_zero_limit',sp.limit(qt,t,0)-sp.Rational(1,140))
Cq=sum(abs(sp.diff(B,w,j).subs(w,a)) for j in range(3,7) for a in (3,4))
Cqd=sum(abs(sp.diff(w*B,w,j).subs(w,a)) for j in range(3,8) for a in (3,4))
energy_bound=sp.simplify((Cqd*Cqd+sp.Rational(13,4)*Cq*Cq)/2)
require(Cq>0 and Cqd>0 and energy_bound>0,'exact decay majorant must be positive')
checks.append({'name':'explicit_finite_order_packet_decay','kind':'EXACT_FINITE_BOUNDARY_EXPANSION',
               'window_support_omega':['3','4'],'window_regular_after_zero_extension':'C2',
               'q_bound_for_abs_t_at_least_one':str(Cq)+'*abs(t)^-4',
               'qdot_bound_for_abs_t_at_least_one':str(Cqd)+'*abs(t)^-4',
               'brane_energy_bound_for_abs_t_at_least_one':str(energy_bound)+'*abs(t)^-8',
               'scope':'Manufactured finite-order packet; not a claim of the all-N C-infinity rate.'})

parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
report={'status':'PASS_EXACT_INDEPENDENT_MANUFACTURED_CONTROLS','checks':checks,
        'check_count':len(checks),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'sympy_version':sp.__version__,'producer_code_imported':False,'protected_input_reads':0,
        'physical_source_callbacks':0,'physical_target_evaluations':0,'network_calls':0,
        'numerical_quadrature_or_interval_claim':False}
args.output.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'status':report['status'],'check_count':len(checks),'output':str(args.output)}))
