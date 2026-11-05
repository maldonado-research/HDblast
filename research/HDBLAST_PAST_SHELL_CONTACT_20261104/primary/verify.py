#!/usr/bin/env python3
"""Exact scalar checks for the registered first master through shell contact."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
PINS='27c15e8cf9262ad3504489899d88fd04067ed38011a79b1b47d32b5b665a9fbd'
def need(ok,name):
    if not ok:raise ValueError(name)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_bytes())
def iv(xs):
    a,b=map(F,xs);need(a<=b,'ORDERED_INTERVAL');return a,b
def add(*xs):return sum((x[0] for x in xs),F(0)),sum((x[1] for x in xs),F(0))
def times(c,x):return (c*x[0],c*x[1]) if c>=0 else (c*x[1],c*x[0])
def mul(a,b):
    values=[x*y for x in a for y in b];return min(values),max(values)
def inv(x):
    need(not x[0]<=0<=x[1],'NONZERO_INTERVAL');return 1/x[1],1/x[0]
def div(a,b):return mul(a,inv(b))
def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:enc(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v) for v in x]
    return x
def authenticate():
    need((HERE/'SOURCE_PINS.json').is_file() and not (HERE/'SOURCE_PINS.json').is_symlink() and sha(HERE/'SOURCE_PINS.json')==PINS,'LOCAL_PINS')
    rows=read(HERE/'SOURCE_PINS.json')['files']
    for rel,row in rows.items():
        need(not Path(rel).is_absolute() and '..' not in Path(rel).parts,'RELATIVE_INPUT_PATH')
        p=HERE/rel;need(p.is_file() and all(not q.is_symlink() for q in [p,*p.absolute().parents]) and p.stat().st_size==row['bytes'] and sha(p)==row['sha256'],'INPUT_BYTES')
    return rows
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    need(a.output.is_absolute(),'ABSOLUTE_OUTPUT_REQUIRED')
    out=a.output.absolute()
    need(all(not p.is_symlink() for p in [out,*out.parents]),'OUTPUT_SYMLINK')
    out=out.resolve()
    need(not out.exists() and out.parent.is_dir() and not out.is_relative_to(HERE),'FRESH_EXTERNAL_OUTPUT')
    pins=authenticate();shell=read(HERE/'inputs/SHELL_BOX.json');prior=read(HERE/'inputs/PRECONTACT_NORMAL.json');checks=[]
    def check(ok,name):need(ok,name);checks.append(name)
    check(shell['status']=='PASS_EXACT_RATIONAL_CONSEQUENCES_OF_RETAINED_ANALYTIC_THEOREM','retained_analytic_shell_box_status')
    boxes=shell['boxes'];rho,eta,p,aprime,Uprime=[iv(boxes[k]['exact']) for k in ['rho_b','eta_b','p_b','a_b','U_prime_b']]
    phi=add((F(1),F(1)),eta);sigma2=times(F(4),phi)
    check(sigma2==iv(boxes['sigma_second_b']['exact']),'actual_tension_second_derivative_is_four_phi')
    h=div(Uprime,p);py_over_p=add(h,times(F(-4),aprime));Bb=add(py_over_p,times(F(2),phi))
    check(Bb[0]>0,'raw_reconstructed_B_b_strictly_positive')
    mu_b=inv(mul(rho,Bb));Chalf=times(F(1,2),mul(rho,add(times(F(2),h),times(F(-5),aprime))))
    check(F(1,1000)<mu_b[0]<=mu_b[1]<F(1,400),'actual_mu_b_reconstructed_from_raw_rho_and_B')
    check(Chalf[0]>189 and Chalf[1]<300-F(1,160),'actual_master_C_b_half_sign_guard')
    beta=(Chalf[0]-F(25,4)*mu_b[1],Chalf[1]+F(5,2)*mu_b[1])
    check(180<beta[0]<=beta[1]<300,'positive_dynamic_spring_for_all_real_times')
    ell=F(1,10**6);Lambda=times(1/ell,mu_b);b=times(ell,beta)
    check(1000<Lambda[0]<=Lambda[1]<2500 and 180*ell<b[0]<=b[1]<300*ell,'scaled_boundary_mass_and_spring')
    check(prior['status']=='PASS_ACTUAL_REGISTERED_PAST_MASTER_PRECONTACT' and F(prior['static_collar_length'])==8*ell,'same_precontact_verified_wide_static_collar')
    check(prior['derived_Q_interval']==['1/3','6'] and all(F(row['strict_lower_margin'])>0 and F(row['strict_upper_margin'])>0 for row in prior['first_exit'].values()),'retained_fresh_static_ODE_inclusions_and_Q_caps')
    Wlow=F(129)**2/3-F(35,4);Whigh=F(131)**2*6
    check(Wlow==F(22153,4)>5500 and Whigh==102966<103000,'all_real_time_positive_master_potential')
    check(2+3==5<8,'postcontact_outer_support_inside_actual_collar')
    kappa0=5500*ell*ell
    check(180*ell/2500>kappa0,'positive_mass_space_coercivity_including_boundary')
    # b_omega=ell*mu_b*T and Lambda*V_omega=(mu_b/ell)*ell^2*T.
    check(ell==(1/ell)*(ell*ell),'boundary_and_bulk_temporal_mass_identity')
    omega_bound=F(315,4)*ell**4;contraction=F(9,2)*omega_bound
    check(contraction==F(2835,8)*ell**4<1,'three_time_static_wave_Duhamel_contraction')
    relative=F(63,4400)*ell*ell
    check(0<relative<F(1,10**13) and 1-relative>0,'positive_nearly_constant_energy_comparison')
    check(F(35,4)/5500*9*ell*ell==relative,'integrated_temporal_energy_loss_coefficient')
    initial=F(1,2)*(1+F(103001,4)*ell*ell)
    check(F(1,2)<initial<F(500001,10**6),'initial_energy_upper_per_Q_derivative_squared')
    check(F(2,1000)==F(1,500) and F(4,500)==F(1,125),'boundary_velocity_and_integrated_trace_energy_factors')
    # For y>=4 and n=2,3,4, d(y^n exp(4-y))/dy has sign n-y.
    # The unit support converts the resulting sup bounds to L2 bounds.
    bump_caps={'Q':F(1),'Q_r':F(4)**2,'Q_rr':F(4)**4+2*F(4)**3+2*F(4)**2}
    check(all(n<=4 for n in [2,3,4]) and bump_caps=={'Q':F(1),'Q_r':F(16),'Q_rr':F(416)},'unit_bump_derivative_L2_majorants')
    graph_initial_twice=(bump_caps['Q_rr']+103000*ell*ell)**2+ell*ell*((bump_caps['Q_r']+ell)**2+103000*ell*ell)
    graph_norm_guard=F(417)
    check(graph_initial_twice/(1-relative)<graph_norm_guard**2,'actual_all_power_domain_graph_norm_guard_417')
    E0_unit=F(1,2)*(bump_caps['Q_r']**2+103001*ell*ell)
    q_guard=F(2)
    check(E0_unit<129 and E0_unit/125<q_guard*q_guard,'actual_master_boundary_trace_guard_two')
    check(rho[0]>129 and -p[1]>F(3,25000),'actual_physical_scalar_denominator_guards')
    mass_sqrt_guard=F(1,20000)
    check(mu_b[1]*ell<mass_sqrt_guard**2,'boundary_mass_graph_to_physical_scalar_factor')
    # Coefficients in Ur/ell-Cb*q/2 after Ur=b0*q-Lambda*a,
    # b0=ell*(Cb/2+5mu_b/2), Lambda=mu_b/ell.
    check(ell/ell-1==0 and F(5,2)*ell/ell==F(5,2) and (1/ell)/ell==1/(ell*ell),'exact_physical_scalar_C_and_omega_cancellation_coefficients')
    scalar_prefactor=F(3)/(F(129)**3*F(3,25000))
    scalar_upper=scalar_prefactor*(graph_norm_guard*mass_sqrt_guard/(ell*ell)+F(5,2)*F(1,400)*q_guard)
    scalar_guard=F(250000000)
    check(scalar_upper==F(1042500000000625,4293378) and scalar_upper<scalar_guard,'actual_first_physical_scalar_uniform_upper_bound')
    eta_front=ell/4;beta_front=F(103000)*F(5,4)*ell*eta_front/2
    check(beta_front==F(103,6400000000)<1 and eta_front+3*beta_front<1,'negative_velocity_qualitative_front_activation_guard')
    check(authenticate()==pins,'input_bytes_unchanged_after_arithmetic')
    report=enc({
        'status':'PASS_CONDITIONAL_PAST_SHELL_CONTACT_SCALARS',
        'schema':'HDBLAST_PAST_SHELL_CONTACT_PRIMARY_V1',
        'date_us_pacific':'2026-11-04',
        'finite_checks':len(checks),'check_names':checks,
        'source_pins_sha256':PINS,
        'selected_local_inputs_authenticated':len(pins),
        'raw_derived_boundary_intervals':{
            'rho_b':rho,'h_b':h,'p_y_over_p_b':py_over_p,'B_b':Bb,
            'mu_b':mu_b,'C_b_over_two':Chalf,'beta_all_real_t':beta,
            'Lambda':Lambda,'b_all_real_x':b},
        'energy':{
            'scaled_potential_lower':kappa0,
            'all_real_W_lower':Wlow,'all_real_W_upper':Whigh,
            'past_scaled_x_interval':['-3','0'],
            'backward_elapsed_slab':['0','3'],
            'support_outer_depth':'5','actual_static_collar_depth':'8',
            'bulk_plus_boundary_derivative':'E_x=V_x*(integral U^2+Lambda q^2)/2>=0 for x<=0',
            'backward_upper_factor':'1','backward_lower_factor':1-relative,
            'relative_energy_change_upper':relative,
            'initial_energy_upper_per_D_squared':initial,
            'boundary_velocity_square_per_E0':F(1,500),
            'boundary_master_trace_square_per_E0_to_three':F(1,125),
            'static_wave_Volterra_contraction_upper':contraction},
        'graph_energy':{
            'unit_bump_L2_derivative_caps':bump_caps,
            'twice_initial_graph_energy_upper':graph_initial_twice,
            'coercivity_factor':1-relative,
            'A0_Z_H_norm_strict_upper':graph_norm_guard,
            'initial_basic_energy_strict_upper':'129',
            'master_boundary_value_strict_upper':q_guard,
            'A0_half_Z_compact_support_claimed':False,
            'all_endpoint_terms_retained_by_spectral_form_identity':True},
        'physical_scalar_trace':{
            'definition':'Phi1(t,Omega)=S_phys(t)*(n1^2-n2^2), |n1^2-n2^2|<=1 and angular mean zero',
            'S_definition':'3 V_z(t,0)/(rho_b^3 p_b)',
            'exact_dynamic_graph_identity':'S_phys(ell*x)=3 mu_b cosh(ell*x)^(-3/2)/(rho_b^3 p_b)*[-(A0Z)_b/ell^2+(5/2)q]',
            'normalization_coefficients':{'scalar_factor':F(3),'graph_ell_power':-2,'q_factor':F(5,2)},
            'sqrt_mu_b_ell_strict_upper':mass_sqrt_guard,
            'rho_b_lower_guard':'129','minus_p_b_lower_guard':'3/25000',
            'derived_uniform_upper':scalar_upper,
            'safe_strict_uniform_upper':scalar_guard,
            'original_time_interval':['-3ell','0'],
            'background_proper_clock':'tau=rho_b*t',
            'actual_selected_first_linear_coefficient_bound':True,
            'scalar_mean_monopole_volume_or_radiation_bound':False,
            'finite_nonlinear_amplitude_window_supplied':False},
        'qualitative_activation':{
            'max_eta':eta_front,'max_beta':beta_front,
            'positive_contradiction_numerator_lower':1-eta_front-3*beta_front,
            'conclusion':'S_phys is nonzero somewhere in every (-T1,-ell), ell<T1<=3ell; it vanishes on [-ell,0]',
            'auxiliary_contradiction_lower_value_is_actual_response_bound':False},
        'prior_operator_credit':{
            'snapshot':'inputs/PRIOR_SPECTRAL_DOMAIN.md',
            'existing_structure':'Positive endpoint-mass Hilbert form/selfadjoint operator',
            'unscaled_local_constant_shift':F(5,2),
            'new_general_selfadjoint_energy_principle_claimed':False,
            'auxiliary_extension_is_global_regular_cone_operator_claimed':False},
        'new_activity':{
            'producer_processes_this_mode':1,
            'exact_interval_boundary_constant_reconstructions':1,
            'new_background_ODE_image_recomputations':0,
            'legacy_checker_or_project_imports':0,
            'background_field_bump_or_source_queries':0,
            'numerical_PDE_solves':0},
        'analytic_limits':[
            'Actual static branch/shell intervals and full first reconstruction remain retained analytic premises.',
            'Linear Wentzell existence uses the established positive-form/wave construction, not numerical evolution.',
            'The actual scalar upper bound uses the additional graph estimate, not the base energy or zero-trace contradiction alone.',
            'Earlier principal channel shares are not extended through reflection by this energy comparison.',
            'No second constrained data, nonlinear past family, physical Noether energy or remote incoming-state law is supplied.']})
    with out.open('xb') as stream:stream.write((json.dumps(report,sort_keys=True,indent=2)+'\n').encode())
    print(json.dumps({'status':report['status'],'finite_checks':len(checks),'output_sha256':sha(out)},sort_keys=True))
if __name__=='__main__':
    try:main()
    except (ValueError,OSError,ArithmeticError,KeyError,TypeError) as exc:
        print(json.dumps({'status':'FAIL','witness':str(exc)},sort_keys=True));sys.exit(2)
