"""Quantum state counts and kinetic-rate benchmarks in the units of the charged wall."""
import argparse
import json
import math
from pathlib import Path

PI=math.pi
DEFAULT_DYNAMICS=Path(__file__).resolve().parents[1]/'coupled'/'DECAY.json'
point={}
def read_primary(filename):
    model=json.loads(filename.read_bytes())
    pars=model['parameters']['primary']
    run=model['runs']['primary']
    trace=[dict(zip(model['trace_columns'],x,strict=True)) for x in run['trace']]
    point.update(a=run['maximum_intrinsic_radius'],theta=run['peak_daughter_theta'],
        D=run['peak_daughter_rest_density'],D_parent_initial=pars['D'],initial_a=2.,
        Ri=pars['R_low'],Rf=pars['R_high'],Q=pars['Q'],K=trace[-1]['K'],
        kappa=pars['kappa'],nu=pars['nu'],flavors=3,spins_per_charge=2,
        peak_time=run['maximum_intrinsic_tau'],peak_daughter_velocity=run['peak_daughter_velocity'],
        peak_converted_fraction=run['peak_daughter_fraction'])
    reconstructed=pars['D']/point['Ri']*(2/point['a'])**3*point['peak_converted_fraction']
    if not math.isclose(reconstructed,point['D'],rel_tol=2e-10):raise ValueError('Daughter number balance')
    return pars,trace
def fermion_gap(mu,a):
    return math.hypot(1.,1.5/(mu*a))

def spinor_capacity(x,mu,a):
    """Free positive-energy Dirac states of one charge/flavor on a round S^3."""
    top=math.floor(x*mu*a-1.5)
    return 0 if top<0 else 2*(top+1)*(top+2)*(top+3)//3

def sample(mu,alpha=None,S0=None):
    a,T,K,Q=point['a'],point['theta'],point['K'],point['Q']
    D=point['D']
    kappa=point['kappa']
    if (S0 is None)==(alpha is None):raise ValueError('Choose a wall action or a charge normalization')
    pref=3*Q**2/(4*PI*K)
    if S0 is None:S0=pref*mu**2/alpha
    else:alpha=pref*mu**2/S0
    W=S0/mu**4
    x_b=math.sqrt(point['Rf']**2-1)
    g=point['spins_per_charge']
    flavors=point['flavors']
    total_species=2*flavors
    xF=(6*PI**2*D*W/(total_species*g))**(1/3)
    total_particles=2*PI**2*a**3*D*S0/mu
    capacities=total_species*spinor_capacity(x_b,mu,a)
    k=3*T
    v=k/math.hypot(k,1.)
    omega_D=math.sqrt(D*Q**2/(K*T))
    # The target density for one independent Abelian force is one opposite
    # sign population (n_total / 6), with Moller velocity 2v.
    phase_cap=(D*W*mu/total_species)*2*v*4*PI/(k*k)
    b=(omega_D/(2*mu*k))**2
    coulomb_log=math.log1p(1/b)-1/(1+b)
    vector_born=phase_cap*.5*alpha**2*coulomb_log
    # Standard spin-averaged s-channel scalar contribution to two other
    # degenerate flavor final states. The total width sums all three flavors.
    y2=4*PI*kappa/(mu*flavors)
    sprime=4*(1+k*k)
    width_over_m=kappa/mu*point['Rf']*(1-point['Rf']**-2)**1.5
    M=2*point['Rf']
    mass_factor=(sprime-4)**2/((sprime-M*M)**2+(M*width_over_m)**2)
    sigma_y=(flavors-1)*y2*y2*mass_factor/(16*PI*sprime)
    scalar_rate=(D*W*mu/total_species)*2*v*sigma_y
    thermal_swave_drag=(D*W*mu/total_species)*math.sqrt(1-point['peak_daughter_velocity']**2)*PI/(T*T)
    thermal_vector_bound=thermal_swave_drag*.5*alpha**2*math.log1p(20*(mu*T/omega_D)**2)
    thermal_swave_collisions=2*thermal_swave_drag
    mfp=1/point['nu']*v
    # If even the typical de Broglie momentum fails to resolve the sphere,
    # the dilute local transport closures should not be evaluated.
    local_continuum=(mu*T*a>10.) and (omega_D/ (mu*T)<.1)
    yukawa_perturbative=kappa/(mu*flavors)<1.
    # Compare number to massless, zero-chemical-potential FD gas at T.
    fd_n_ratio=(D*W/total_species)/(g*3*1.2020569031595943*T**3/(4*PI**2))
    return dict(mu=mu,S0=S0,q=math.sqrt(3)*Q*mu,late_alpha=alpha,
        hot_local_resolution_and_scale_separation=local_continuum,
        per_vertex_yukawa_alpha_below_one=yukawa_perturbative,
        early_alpha_if_same_charged_species=alpha*K/.001,
        total_yukawa_alpha=kappa/mu,per_flavor_yukawa_squared=y2,used_daughter_rest_density=D,
        total_carrier_count=total_particles,cold_fermi_momentum_over_m=xF,
        all_cold_fraction_upper_continuum=min(1.,(x_b/xF)**3),
        free_S3_threshold_capacity=capacities,
        free_S3_capacity_over_required=capacities/total_particles,
        free_S3_gap_over_m=fermion_gap(mu,a),one_flavor_hot_N_over_mu0_N=fd_n_ratio,
        ideal_Debye_inverse_length_times_L=omega_D,
        ideal_Debye_particles_in_sphere=(4*PI/3)*(D*S0/mu/flavors)/(omega_D**3),
        k_cm_over_m=k,coulomb_log_bracket=coulomb_log,
        S_wave_rate_ceiling_times_L=phase_cap if local_continuum else None,
        required_s_wave_fraction=point['nu']/phase_cap if local_continuum else None,
        thermal_swave_drag_ceiling_times_L=thermal_swave_drag if local_continuum else None,
        thermal_swave_collision_frequency_ceiling_times_L=thermal_swave_collisions if local_continuum else None,
        required_thermal_swave_fraction=point['nu']/thermal_swave_drag if local_continuum else None,
        thermal_screened_vector_drag_upper_times_L=thermal_vector_bound if local_continuum else None,
        screened_vector_transport_rate_times_L=vector_born if local_continuum else None,
        scalar_other_flavors_rate_times_L=scalar_rate if local_continuum and yukawa_perturbative else None,
        scalar_rate_to_requested=scalar_rate/point['nu'] if local_continuum and yukawa_perturbative else None,
        mfp_over_wall_radius_at_requested_rate=mfp/a)
def gauss(n,lower,upper):
    pts=[]
    for j in range(1,(n+1)//2+1):
        z=math.cos(PI*(j-.25)/(n+.5))
        for _ in range(30):
            p,q=1.,z
            for k in range(2,n+1):
                p,q=q,((2*k-1)*z*q-(k-1)*p)/k
            der=n*(p-z*q)/(1-z*z)
            dz=q/der
            z-=dz
            if abs(dz)<3e-16:break
        w=(upper-lower)/((1-z*z)*der*der)
        pts.extend([((upper+lower-z*(upper-lower))/2,w),
                    ((upper+lower+z*(upper-lower))/2,w)])
    return sorted(pts)


def thermal_weights(n):
    # Momentum p/T has density x^2 e^-x / 2. The pair's direction cosine
    # is uniform on [-1,1]. All omitted Boltzmann tails are above p/T=50.
    nodes=[(x,w*x*x*math.exp(-x)/2) for x,w in gauss(n,0.,50.)]
    angle=gauss(n,-1.,1.)
    weights=[]
    for i,(x,wx) in enumerate(nodes):
        for j in range(i+1):
            y,wy=nodes[j]
            multiplicity=1 if i==j else 2
            for c,wc in angle:
                W=x/y+y/x-2*c
                weights.append((x*y*(1-c)/2,wx*wy*wc*.5*multiplicity*W/3))
    normalization=math.fsum(w for _,w in weights)
    assert abs(normalization-1)<2e-8
    return weights,normalization

def transport_factor(weights,k0_over_temperature):
    b=k0_over_temperature*k0_over_temperature
    return math.fsum(w*k/(k+b) for k,w in weights)

def elastic_kernel(n,row):
    weights,normalization=thermal_weights(n)
    want=point['nu']/row['thermal_swave_drag_ceiling_times_L']
    lo,hi=0.,100.
    for _ in range(45):
        mid=(lo+hi)/2
        if transport_factor(weights,mid)>want:lo=mid
        else:hi=mid
    val=(hi+lo)/2
    return dict(k0_over_T=val,k0_times_L=val*row['mu']*point['theta'],
        k0_over_m=val*point['theta'],fitted_nu_times_L=transport_factor(weights,val)*row['thermal_swave_drag_ceiling_times_L'],
        phase_space_quadrature_normalization=normalization)


def fixed_kernel_history(trace,corner,fit):
    weights,_=thermal_weights(40)
    out=[]
    for t in (11.1,11.4,11.5,11.7,12.,12.4):
        z=min(trace,key=lambda r:abs(r['t']-t))
        gamma=math.hypot(1,z['u'])
        v=z['u']/gamma
        if abs(v)>.1 or z['theta']<100:raise ValueError('Outside near-isotropic hot benchmark')
        # One opposite-sign flavor density, transformed from lab to the stream rest frame.
        nt_over_m3=z['rest']*corner['S0']/(6*gamma*corner['mu']**4)
        cap=PI*nt_over_m3*corner['mu']/z['theta']**2
        physical_scale_ratio=fit['k0_over_m']/z['theta']
        nu=cap*transport_factor(weights,physical_scale_ratio)
        out.append(dict(t=z['t'],theta=z['theta'],stream_velocity=v,converted_fraction=z['f'],
             target_proper_density_over_m3=nt_over_m3,k0_over_local_T=physical_scale_ratio,
             thermal_drag_times_L=nu,fixed_fluid_drag_times_L=point['nu']))
    return out


def phase_supply(trace,pars,reference):
    """Collisionless-shell crossing capacity for instantaneous isotropic births."""
    rows=[]
    for t in (10.36,10.4,10.5,10.7,11.0,12.0):
        z=min(trace,key=lambda r:abs(r['t']-t))
        R=z['R']
        if R<=1:continue
        x=(z['psi']-pars['psi_star'])/pars['width_phi']+.5
        dR=0.
        if 0<x<1:
            step=(R-pars['R_low'])/(pars['R_high']-pars['R_low'])
            dR=((pars['R_high']-pars['R_low'])*step*(1-step)*
                (1/(1-x)**2+1/x**2)*z['pi']/pars['width_phi'])
        birth=math.sqrt(R*R-1)
        speed=abs(R*dR/birth+z['H']*birth)
        force=abs(point['Q']*z['X'])
        mean=speed if speed>=force else (speed*speed+force*force)/(2*force)
        supply=z['rate']*(z['parent']/R)*reference['S0']/reference['mu']**4
        capacity=12/(2*PI**2)*birth**2*mean
        rows.append(dict(t=z['t'],birth_momentum_over_m=birth,drift_speed_in_p_over_m=force,
                    fraction_converted=z['f'],requested_over_collisionless_capacity=supply/capacity))
    return rows


def calculate(path=DEFAULT_DYNAMICS):
    parameters,trace=read_primary(path)
    a=point['a'];Q=point['Q'];K=point['K'];T=point['theta']
    D=point['D']
    x_b=math.sqrt(point['Rf']**2-1); g=point['spins_per_charge']
    k=3*T;v=k/math.hypot(k,1.)
    c0=g*x_b**3/(PI**2*D) # cold allows S0 <= c0 * mu^4
    c1=3*point['nu']*k**2/(4*PI*D*v) # hard rate needs S0 >= c1 * mu^3
    min_mu=c1/c0
    thermal_c1=6*point['nu']*T*T/(PI*D)
    mu_Q1=1/(math.sqrt(3)*Q)
    rows={'small_phase_charge':sample(mu_Q1,S0=1.),
          'visible_coupling_mu100':sample(100.,alpha=1/137.036),
          'visible_coupling_mu1e5':sample(1.e5,alpha=1/137.036),
          'visible_coupling_mu1000':sample(1000.,alpha=1/137.036),
          'weak_yukawa_and_perturbative_vector':sample(1000.,alpha=.01),
          'conditional_large_mass_phase_shift':sample(1.e11,S0=1.e47)}
    # Increase mu and S0/mu^3 if necessary to actually realize an s-wave
    # elastic rate in the last row while passing the cold packing test.
    corner=rows['conditional_large_mass_phase_shift']
    if not (corner['required_s_wave_fraction']<1 and corner['all_cold_fraction_upper_continuum']==1):
        raise ValueError('The stated illustrative normalization is inconsistent')
    fit=elastic_kernel(40,corner)
    fit_coarse=elastic_kernel(24,corner)
    fit['change_24_to_40']=abs(fit['k0_over_T']-fit_coarse['k0_over_T'])
    assert fit['change_24_to_40']<.001 and abs(fit['fitted_nu_times_L']-point['nu'])<1e-10
    assert abs(rows['visible_coupling_mu100']['late_alpha']-1/137.036)<1.e-15
    assert all(row['S_wave_rate_ceiling_times_L'] is None or math.isfinite(row['S_wave_rate_ceiling_times_L']) for row in rows.values())
    return dict(evaluation_point=point,physical_conventions='hbar=c=1, mu=mL, S0=sigma L^4',
        closed_free_spinor_lowest_momentum_times_aL=1.5,
        y_threshold_width_times_L=point['kappa']*point['Rf']*(1-point['Rf']**-2)**1.5,
        daughter_birth_momentum_over_m_at_Rf=x_b,
        min_mu_for_free_pair_at_Rf_on_fieldfree_sphere=1.5/(a*x_b),
        c0_cold_S0_upper_per_mu4=c0,c1_hard_swave_S0_lower_per_mu3=c1,
        min_mu_for_simultaneously_cold_packed_and_swave_rate=min_mu,
        c1_thermal_swave_drag_S0_lower_per_mu3=thermal_c1,
        min_mu_for_cold_packed_and_thermal_swave_drag=thermal_c1/c0,
        calibrated_additional_elastic_kernel=fit,
        fixed_elastic_kernel_along_hot_trace=fixed_kernel_history(trace,corner,fit),
        sampled_birth_supply_visible_mu100=phase_supply(trace,parameters,rows['visible_coupling_mu100']),rows=rows)
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--dynamics',type=Path,default=DEFAULT_DYNAMICS)
    args=parser.parse_args()
    args.output.write_text(json.dumps(calculate(args.dynamics),indent=2,sort_keys=True,allow_nan=False)+'\n')
