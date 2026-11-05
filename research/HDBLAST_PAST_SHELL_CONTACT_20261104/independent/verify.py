"""Exact independent raw-coefficient, past-contact and physical-trace margins."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
PINS_SHA='25a496354e4bb278faefc7713dc1973a32c1673077773f14c7860f6544f029f6'

def need(ok,message):
    if not ok:raise ValueError(message)
def sha(raw):return hashlib.sha256(raw).hexdigest()
def regular(p):return p.is_file() and all(not q.is_symlink() for q in [p,*p.parents])
def add(a,b):return a[0]+b[0],a[1]+b[1]
def scale(c,a):return (c*a[0],c*a[1]) if c>=0 else (c*a[1],c*a[0])
def mul(a,b):
    v=[x*y for x in a for y in b];return min(v),max(v)
def poly(c,x):return sum((v*x**n for n,v in c.items()),F(0))
def encode(v):
    if isinstance(v,F):return str(v)
    if isinstance(v,dict):return {k:encode(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [encode(x) for x in v]
    return v

def calculate():
    pins_path=HERE/'SOURCE_PINS.json';need(regular(pins_path),'SOURCE_PINS_REGULAR')
    pins_raw=pins_path.read_bytes();need(sha(pins_raw)==PINS_SHA,'SOURCE_PINS_BYTES');pins=json.loads(pins_raw)['files'];inputs={}
    def authenticate():
        for name,row in pins.items():
            p=HERE/name;need(regular(p),'INPUT_REGULAR:'+name);raw=p.read_bytes()
            need(len(raw)==row['bytes'] and sha(raw)==row['sha256'],'INPUT_BYTES:'+name)
            inputs[name]=raw
    authenticate();shell=json.loads(inputs['inputs/SHELL_BOX.json']);wide=json.loads(inputs['inputs/WIDE_BACKGROUND.json'])
    checks=[]
    def gate(name,condition):
        need(name not in checks,'DUPLICATE_GATE');need(condition,name);checks.append(name)
    ell=F(1,10**6);length=8*ell
    gate('explicit_wide_static_collar',F(wide['candidate_state_exact']['Z'])==length and F(wide['candidate_state_exact']['ell'])==ell)
    outer={'rho':(F(129),F(131)),'eta':(-F(9,100000),-F(8,100000)),
           'p':(-F(7,50000),-F(3,25000)),'a':(F(111,1000),F(14,125))}
    rho,eta,p,a=(outer[k] for k in ['rho','eta','p','a'])
    U1={1:F(28,9),2:F(50,9),3:-F(2,3),4:-F(20,9),5:-F(4,9)}
    U2={0:F(28,9),1:F(100,9),2:-F(2),3:-F(80,9),4:-F(20,9)}
    gate('polynomial_endpoint_monotonicity',F(100,9)-F(80,3)*eta[0]**2>0)
    U2box=tuple(poly(U2,t) for t in eta);U1box=tuple(poly(U1,t) for t in eta)
    gate('wide_U2_and_h_enclosures',F(31,10)<U2box[0]<U2box[1]<F(25,8)
         and F(7,4)<U1box[1]/p[0]<U1box[0]/p[1]<F(7,3))
    rhs={'rho':scale(-1,mul((rho[0]**2,rho[1]**2),a)),
         'eta':scale(-1,mul(rho,p)),
         'p':scale(-1,mul(rho,add(U1box,scale(-4,mul(a,p))))),
         'a':add((1/rho[1],1/rho[0]),scale(F(1,3),mul(rho,(p[1]**2,p[0]**2))))}
    first_exit={}
    for name in outer:
        initial=tuple(F(v) for v in shell['boxes'][name+'_b']['exact'])
        tube=add(initial,mul((F(0),length),rhs[name]))
        gate('fresh_first_exit_'+name,outer[name][0]<tube[0]<=tube[1]<outer[name][1])
        first_exit[name]={'tube':tube,'lower_margin':tube[0]-outer[name][0],'upper_margin':outer[name][1]-tube[1]}
    gate('physical_collar_endpoint_margin',rho[1]*length<9)
    gate('joint_Q_monotonicity',4*F(7,4)-10*a[1]>0 and -10*F(7,4)+F(35,2)*a[1]<0)
    qlo=2*F(7,4)**2-10*a[1]*F(7,4)+F(35,4)*a[1]**2-F(25,8)-p[0]**2/6
    qhi=2*F(7,3)**2-10*a[0]*F(7,3)+F(35,4)*a[0]**2-F(31,10)
    gate('fresh_joint_background_Q_bounds',qlo>1 and qhi<F(16,3))
    gate('all_real_positive_W_and_auxiliary_extension',rho[0]**2-F(35,4)>5500
         and rho[1]**2*F(16,3)<103000 and F(5547)<6000<102966
         and F(5547)-F(35,4)>5500)
    raw=lambda name:tuple(F(x) for x in shell['boxes'][name]['exact'])
    rb,etab,pb,ab=(raw(k+'_b') for k in ['rho','eta','p','a'])
    ub=tuple(poly(U1,t) for t in etab)
    hb=(ub[1]/pb[0],ub[0]/pb[1])
    sigma2=scale(4,add((F(1),F(1)),etab))
    gate('tension_second_derivative_four_phi',sigma2==raw('sigma_second_b'))
    B=add(add(hb,scale(-4,ab)),scale(F(1,2),sigma2))
    C_half=mul(rb,add(hb,scale(-F(5,2),ab)))
    gate('raw_endpoint_B_positive',B[0]>0)
    mu=(1/(rb[1]*B[1]),1/(rb[0]*B[0]))
    gate('raw_reconstructed_boundary_mass',F(1,1000)<mu[0]<mu[1]<F(1,400))
    gate('longitudinal_C_half_bounds',189<C_half[0]<C_half[1]<300-F(1,160))
    beta=(C_half[0]-F(25,4)*mu[1],C_half[1]+F(5,2)*mu[1])
    gate('positive_boundary_spring_all_real_times',180<beta[0]<beta[1]<300)
    Lambda=scale(1/ell,mu);b=scale(ell,beta)
    gate('normal_and_time_scaling_guards',1000<Lambda[0]<Lambda[1]<2500 and 180*ell<b[0]<b[1]<300*ell)
    gap=5500*ell**2;boundary_gap=b[0]/Lambda[1]
    gate('bulk_and_boundary_H_coercivity',boundary_gap>gap and F(180)*ell/2500>gap)
    omega_max=F(315,4)*ell**4;contraction=F(9,2)*omega_max
    gate('three_time_static_wave_contraction',contraction==F(2835,8)*ell**4 and contraction<1)
    epsilon=F(63,4400)*ell**2
    gate('direct_integrated_relative_energy_factor',omega_max/gap==epsilon and 0<epsilon<1)
    E0_per_D2=(1+103001*ell**2/4)/2
    gate('initial_energy_per_D_squared',E0_per_D2>F(1,2) and E0_per_D2<F(501,1000))
    gate('trace_velocity_and_history_factors',F(2,1000)==F(1,500) and F(4,500)==F(1,125))
    Q1=F(16);Q2=F(4)**4+2*F(4)**3+2*F(4)**2
    gate('explicit_bump_jet_envelope_constants',Q1==16 and Q2==416)
    V0_hi=103000*ell**2
    graph_zero_twice=(Q2+V0_hi)**2+ell**2*(Q1+ell)**2+103000*ell**4
    graph_norm_squared=graph_zero_twice/(1-epsilon)
    gate('higher_graph_norm_below_417',graph_norm_squared<417**2)
    E0_abs=(Q1**2+103001*ell**2)/2;q_squared=E0_abs/125
    gate('actual_master_trace_below_two',E0_abs<129 and q_squared<4)
    root_weight=F(1,20000)
    gate('boundary_coordinate_weight',root_weight**2==F(1,400)*ell)
    gate('physical_scalar_temporal_cancellation_constant',F(35,4)-F(25,4)==F(5,2))
    graph_term=root_weight*417/ell**2;trace_term=F(5,2)*F(1,400)*2
    S_upper=3*(graph_term+trace_term)/(F(129)**3*F(3,25000))
    gate('actual_first_physical_scalar_coefficient_bound',graph_term==20850000000 and trace_term==F(1,80)
         and S_upper==F(1042500000000625,4293378)<250000000)
    activation_eta=ell/4;activation_beta=F(5,32)*103000*ell**2
    activation_numerator=1-activation_eta-3*activation_beta
    gate('qualitative_backward_front_activation',activation_beta<1 and activation_numerator>0)
    authenticate();need(sha(pins_path.read_bytes())==PINS_SHA,'SOURCE_PINS_BYTES_AFTER')
    return encode({
      'schema':'HDBLAST_PAST_CONTACT_INDEPENDENT_EXACT_V1','status':'PASS_CONDITIONAL_PAST_CONTACT_AND_ACTUAL_FIRST_SCALAR_BOUND',
      'date_label':'4 November 2026, America/Los_Angeles','source_pins_sha256':PINS_SHA,
      'first_exit':first_exit,'source_route':'Raw shell eta,p,rho,a plus exact polynomial endpoint monotonicity; stored derived Uprime/h/B/mu/C boxes and old precontact reports are not numerical premises.',
      'raw_boundary':{'U1_from_eta':ub,'h':hb,'B_b':B,'C_b_over_two':C_half,'mu_b':mu,'beta_all_real':beta,'Lambda':Lambda,'b_all_real':b},
      'coefficient_guards':{'mu_b':['1/1000','1/400'],'beta':['180','300'],'Lambda':['1000','2500'],
        'actual_collar_W':['5500','103000'],'static_collar_scaled':'8','postcontact_outer_support_scaled':'5',
        'joint_Q_lower':qlo,'joint_Q_upper':qhi,'gap':gap,'boundary_gap_lower':boundary_gap},
      'energy':{'scaled_time_slab':['-3','0'],'static_wave_Volterra_norm_upper':contraction,
        'relative_energy_loss_upper':epsilon,'backward_energy_factor_lower':1-epsilon,'backward_energy_factor_upper':'1',
        'monotonicity_identity':'E_x=omega_x*(integral U^2+Lambda q^2)/2>=0 on x<=0',
        'initial_energy_upper_per_D_squared':E0_per_D2,'boundary_velocity_squared_per_E0':'1/500','boundary_trace_squared_per_E0':'1/125'},
      'higher_graph':{'Q_L2_upper':'1','Qr_L2_upper':Q1,'Qrr_L2_upper':Q2,'twice_initial_energy_upper':graph_zero_twice,
        'A0Z_norm_squared_upper':graph_norm_squared,'A0Z_norm_strict_upper':'417','base_initial_energy_upper':E0_abs,
        'q_squared_upper':q_squared,'q_absolute_strict_upper':'2','mu_over_sqrt_Lambda_upper':root_weight,
        'physical_scalar_coefficient_absolute_upper':S_upper,'safe_strict_physical_scalar_coefficient_upper':'250000000',
        'actual_pointwise_first_variation':'Phi1=S(t)Y, |Y|<=1','nonlocal_half_power_support_assumed':False,
        'identity':'S=3mu_b cosh(t)^(-3/2)/(rho_b^3 p_b)*[-(A0Z)_b/ell^2+(5/2)q]'},
      'activation':{'backward_observation_time':'s_obs=-t=ell+eta','max_eta':activation_eta,'max_beta':activation_beta,
        'contradiction_numerator_lower':activation_numerator,'positive_auxiliary_value_is_actual_response_lower_bound':False},
      'verification':{'finite_scalar_gates':len(checks),'check_names':checks},
      'activity':{'fresh_interval_ODE_images':4,'raw_boundary_reconstructions':1,'legacy_checker_or_project_imports':0,
        'background_field_bump_or_source_queries':0,'numerical_PDE_solves':0,'producer_processes_this_mode':1},
      'limits':['The closed form, smooth branch, physical master/reconstruction and gauge identification remain analytic premises.',
        'The graph identity uses functional calculus, not finite propagation for the nonlocal half-power state.',
        'The bound is for the actual first coefficient in fixed model units, not finite nonlinear amplitude or SI calibration.',
        'No physical Noether energy, second constrained preparation, remote incoming law or nonlinear parameter window is supplied.',
        'Precontact principal-channel shares are not propagated through reflection by this energy result.'],
    })

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args()
    need(a.output.is_absolute(),'ABSOLUTE_OUTPUT');original=a.output.absolute()
    need(all(not q.is_symlink() for q in [original,*original.parents]),'OUTPUT_SYMLINK')
    out=original.resolve();need(not out.exists() and out.parent.is_dir(),'FRESH_REGULAR_OUTPUT')
    need(out!=HERE and HERE not in out.parents,'OUTPUT_OUTSIDE_COMPONENT')
    result=calculate();raw=(json.dumps(result,sort_keys=True,indent=2)+'\n').encode()
    with out.open('xb') as stream:stream.write(raw)
    print(json.dumps({'status':result['status'],'finite_scalar_gates':result['verification']['finite_scalar_gates'],'output_sha256':sha(raw)},sort_keys=True))

if __name__=='__main__':
    try:main()
    except (ValueError,OSError,KeyError,TypeError,ArithmeticError) as e:
        print(json.dumps({'status':'FAIL','witness':str(e)},sort_keys=True));raise SystemExit(2)
