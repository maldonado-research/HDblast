"""Compare the finite closed wall's turnaround energy budget with present photon entropy."""
import argparse
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent
HBARC=1.973269804e-16       # GeV m
HBAR=6.582119569e-25        # GeV s
KB=8.617333262145e-14       # GeV/K
C_LIGHT=299792458.           # m/s
MPC=3.085677581491367e22    # m
MPL=2.435e18                 # reduced Planck mass, GeV
TEMPERATURE=2.72548
HUBBLE=67.4
GSM=106.75


def hubble_photons(temperature,hubble):
    radius=C_LIGHT/(hubble*1000./MPC)
    entropy_density=4*math.pi**2/45*(KB*temperature/HBARC)**3
    return 4*math.pi/3*radius**3*entropy_density


def entropy_on_wall(a,e,action,g):
    # Equilibrium ideal gas, zero chemical potentials, complete round 3-sphere.
    temp_length=(30*e*action/(math.pi**2*g))**.25
    return 2*math.pi**2*a**3 * 2*math.pi**2*g/45 * temp_length**3


def entropy_at_temperature(aH, target, temp_GeV=.005, g_energy=10.75, g_entropy=10.75):
    # R^2 rho_rad <= 3 M_Pl^2 [1+(aH)^2] after matching conventional G.
    coefficient=4*math.pi*90**1.5/45 * g_entropy/g_energy**1.5
    at_rest=coefficient*(MPL/temp_GeV)**3
    factor=(target/at_rest)**(1/3)
    cap=at_rest*(1+aH*aH)**1.5
    return dict(test_temperature_GeV=temp_GeV,g_rho=g_energy,g_s=g_entropy,
                measured_aH=aH,entropy_upper_in_kB=cap,photon_entropy_over_upper=target/cap,
                minimum_aH_at_this_temperature=math.sqrt(max(0.,factor**2-1.)))


def calculate(row,mu,action,target):
    a=row['maximum_radius']
    t=row['maximum_tau']
    if not row['stop_turnaround'] or not (a>0 and t>0):raise ValueError("Missing turnaround")
    e_bound=1/(2*a*a)
    upper=entropy_on_wall(a,e_bound,action,GSM)
    # S = C_g (a^4 e S0)^(3/4). Both algebraic forms must agree.
    constant=entropy_on_wall(1,1,1,GSM)
    if not math.isclose(upper,constant*(a*a*action/2)**.75,rel_tol=2e-14):
        raise ValueError("Energy-to-entropy conversion disagrees")
    length=math.sqrt(action/6.)/MPL
    a_d,r_d,theta=row['doubling_radius'],row['doubling_rest'],row['doubling_theta']
    intrinsic=r_d/a_d**3*(1.5*theta+math.hypot(1.,1.5*theta)-1.)
    e_d=intrinsic/row['doubling_hot_fraction']
    temp_eq=(30*action*e_d/(math.pi**2*GSM))**.25/length
    required=2*(target/constant)**(4/3)
    r=dict(action=action,daughter_mass_L=mu,turn_dimensionless_radius=a,
           tight_energy_ceiling=1/(a*a*(math.sqrt(1+1/a**2)+1)),
           radiation_entropy_upper_in_kB=upper,photon_entropy_over_upper=target/upper,
           action_needed_at_same_turn_radius=required/a**2,
           curvature_fraction_magnitude_at_hot_doubling=1/row['doubling_aH']**2,
           thermal_benchmark_entropy_at_doubling_in_kB=entropy_on_wall(a_d,e_d,action,GSM),
           with_observed_G=dict(length_m=length*HBARC,daughter_mass_GeV=mu/length,
                     initial_to_turn_seconds=t*length*HBAR,turn_radius_m=a*length*HBARC,
                     doubling_T_kinetic_GeV=theta*mu/length,doubling_T_equilibrium_GeV=temp_eq,
                     minimum_closed_turn_radius_for_target_Mpc=math.sqrt(required/6)*HBARC/MPL/MPC))
    r['with_observed_G']['five_MeV_at_doubling']=entropy_at_temperature(row['doubling_aH'],target)
    if 'first_hot_aH' in row:
        r['with_observed_G']['five_MeV_at_first_hot']=entropy_at_temperature(row['first_hot_aH'],target)
    if 'max_stored_a4e' in row:
        gmax=row['max_curvature_budget']
        emax=row['max_stored_a4e']
        if row['min_wall_mass_to_turn']<=0 or 2*emax>gmax*(1+2e-5):
            raise ValueError("Energy and curvature inventories disagree")
        r.update(largest_sampled_equilibrium_entropy_before_turn_in_kB=constant*(action*emax)**.75,
                 sampled_maximum_time=row['max_stored_time'],
                 curvature_based_pre_turn_entropy_ceiling_sampled_in_kB=constant*(action*gmax/2.)**.75)
    if 'turn_wall_mass' in row:
        m=row['turn_wall_mass'];e=row['turn_wall_e']
        residual=abs(e*(2+e)+m/a**4-1/a**2)
        if not (m>=0 and 0<e<=e_bound and residual<1e-7):raise ValueError("Turnaround energy identity failed")
        r.update(measured_turn_bulk_mass=m, measured_turn_wall_density=e,
                 turnaround_constraint_residual=residual,
                 thermal_benchmark_entropy_at_turn_in_kB=entropy_on_wall(a,e,action,GSM))
    elif row['min_bulk_mass']<=0:raise ValueError("Background bulk-mass assumption failed")
    return r


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base',type=Path,default=HERE.parent/'prior')
    parser.add_argument('--feedback',type=Path,default=HERE.parent/'feedback/FEEDBACK.json')
    parser.add_argument('--evolution',type=Path,default=HERE.parent/'boost/EVOLUTION.json')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    env=json.loads((args.base/'envelope/ENVELOPE.json').read_text())
    feed=json.loads(args.feedback.read_text())
    sweep=json.loads((args.base/'kinetics/SWEEP.json').read_text())['runs']
    target=hubble_photons(TEMPERATURE,HUBBLE)
    result={'constants':dict(Tgamma_K=TEMPERATURE,H0_km_s_Mpc=HUBBLE,radiation_dof=GSM,
                reduced_Planck_mass_GeV=MPL,hbarc_GeV_m=HBARC,hbar_GeV_s=HBAR,
                Boltzmann_GeV_K=KB,c_m_s=C_LIGHT,Mpc_m=MPC),
            'photon_entropy_Hubble_sphere_in_kB':target,
            'endpoint_sensitivity':dict(lower=hubble_photons(TEMPERATURE-.00057,HUBBLE+.5),
                                        upper=hubble_photons(TEMPERATURE+.00057,HUBBLE-.5)),
            'runs':{}}
    for case in ('mass_4e10','mass_4e12'):
        scale=sweep[case]['action']
        mu=env['parameters'][case]['daughter_mass']
        ref=env['runs'][case]
        rows={'unlimited_archived':calculate(ref,mu,scale,target)}
        for name,model in feed[case].items():
            for res,row in model.items():
                rows[name+'_'+res]=calculate(row,mu,scale,target)
        # Use the most permissive angle-integrated sharp-shell upper action on this
        # exact archived trajectory, even where a more restrictive angular test fails.
        ceiling=sweep[case]['runs']['half_step']['continuum_pointwise_S0_upper']
        rows['fixed_archived_background_at_central_sweep_action_ceiling']=calculate(ref,mu,ceiling,target)
        result['runs'][case]=rows
    pole=json.loads(args.evolution.read_bytes())
    result['explicit_pole_runs']={
        name:calculate(row,pole['daughter_mass_L'],pole['action'],target)
        for name,row in pole['runs'].items()
        if name.startswith('pole_') and 'doubling_radius' in row and row.get('hot_radius_factor',0.)>=2.
    }
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')

if __name__=='__main__':
    main()
