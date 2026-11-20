"""Cold charged tracers on an imported closed-wall Maxwell trajectory (stdlib only)."""
import argparse
import bisect
import csv
import gzip
import json
import math
import random
from pathlib import Path
from typing import NamedTuple

START = 4.6
MASS = 0.1
CHARGE = 1.0
PROJECTED = CHARGE / math.sqrt(3)
SEED_FRACTION = 1e-4
RATES = (0.0, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0)


class Knot(NamedTuple):
    t: float
    a: float
    H: float
    potential: float
    adot_potential: float
    density: float
    sheet: float
    K: float


def interpolate(left, right, t):
    f = (t - left.t) / (right.t - left.t)
    return Knot(t, *(x + f * (y - x) for x, y in zip(left[1:], right[1:])))


def trajectory(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as source:
        rows = []
        previous = None
        for item in csv.DictReader(source):
            t, a, h, b, x, e, sheet, k = (float(item[n]) for n in
                                          ('t','a','H','Bmag','X','e','sheet','K'))
            current = Knot(t, a, h, a*a*b/2, a*x, e, sheet, k)
            if previous and previous.t < START <= current.t:
                rows.append(interpolate(previous, current, START))
            if t > START:
                if rows and current.H <= 0 < rows[-1].H:
                    dt = current.t-rows[-1].t
                    turn = rows[-1].t + dt * rows[-1].H / (rows[-1].H-current.H)
                    rows.append(interpolate(rows[-1],current,turn))
                    break
                rows.append(current)
            previous = current
    if len(rows)<3 or abs(rows[0].t-START)>1e-10 or abs(rows[-1].H)>1e-8:
        raise ValueError('The input must bracket the start and the subsequent expansion maximum.')
    return rows


def deriv(y, b, nu, a, h):
    u, s, d, c, _, _, _, _ = y
    par2 = (s+2*d)/3
    return (b-nu*u, 2*b*u, 2*b*u-nu*d, b*(s+2*par2)-nu*c, 4*b*c,
            2*nu*u*u, b*u/(MASS*a*a), h*s/(MASS*a*a))


def step(y, b, nu, a0, a1, dt):
    h = math.log(a1/a0)/dt
    amid = math.sqrt(a0*a1)
    k1 = deriv(y,b,nu,a0,h)
    k2 = deriv([v+0.5*dt*k for v,k in zip(y,k1)],b,nu,amid,h)
    k3 = deriv([v+0.5*dt*k for v,k in zip(y,k2)],b,nu,amid,h)
    k4 = deriv([v+dt*k for v,k in zip(y,k3)],b,nu,a1,h)
    return [v+dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6 for i,v in enumerate(y)]


def run(rows, nu):
    y = [0.0]*8
    start, end = rows[0], rows[-1]
    n0 = SEED_FRACTION*start.density/MASS      # both signs combined
    peak_seed, peak_pickup, peak_u = 0.0, 0.0, 0.0
    total_variation = 0.0
    for x,z in zip(rows[:-1],rows[1:]):
        dt=z.t-x.t
        b=PROJECTED*(z.potential-x.potential)/dt
        # Keep time-integration checks valid even on a coarse, irregular input.
        count=max(1, math.ceil(nu*dt/.15))
        for j in range(count):
            a0=x.a*(z.a/x.a)**(j/count)
            a1=x.a*(z.a/x.a)**((j+1)/count)
            y=step(y,b,nu,a0,a1,dt/count)
        total_variation += abs(z.potential-x.potential)
        u,s=y[:2]
        n=n0*(start.a/z.a)**3
        peak_seed=max(peak_seed,n*(MASS+s/(2*MASS*z.a*z.a))/z.density)
        peak_pickup=max(peak_pickup,n*s/(2*MASS*z.a*z.a*z.sheet))
        peak_u=max(peak_u,PROJECTED*total_variation/(MASS*z.a))
    u,s,d,c,fourth,var_integral,work,redshift=y
    par2=(s+2*d)/3
    var=max(0.0,s-u*u)
    kurtosis=fourth-4*u*c+2*u*u*s+4*u*u*par2-3*u**4
    n=n0*(start.a/end.a)**3
    scale=2*MASS*end.a**2
    return dict(rate=nu, mean_comoving_momentum=u, mean_squared_comoving_momentum=s,
                comoving_momentum_variance=var, random_share=var/s if s and nu else 0.,
                mean_directed_energy_per_carrier=u*u/scale,
                random_energy_per_carrier=var/scale if nu else 0.,
                random_energy_over_total_EM=n*var/(scale*end.sheet) if nu else 0.,
                total_kinetic_over_total_EM=n*s/(scale*end.sheet),
                central_pressure_difference_over_trace=(d-u*u)/var if nu and var else None,
                central_fourth_moment_over_Maxwell_value=3*kurtosis/(5*var*var) if nu and var else None,
                probability_no_deflection=math.exp(-nu*(end.t-start.t)),
                work_per_carrier=work, expansion_pressure_work_per_carrier=redshift,
                energy_balance_absolute=abs(work-redshift-s/scale),
                variance_balance_absolute=abs(var-var_integral), max_tracer_over_wall=peak_seed,
                max_kinetic_over_total_EM=peak_pickup, max_momentum_over_mass_bound=peak_u)


def add_checks(rows, results):
    total_variation=sum(abs(y.potential-x.potential) for x,y in zip(rows[:-1],rows[1:]))
    integral=sum((y.t-x.t)*(abs(x.adot_potential)+abs(y.adot_potential))/2
                 for x,y in zip(rows[:-1],rows[1:]))
    start,end=rows[0],rows[-1]
    raw=results[0]
    expected=PROJECTED*(end.potential-start.potential)
    if not math.isclose(raw['mean_comoving_momentum'],expected,abs_tol=3e-12):
        raise ValueError('Collisionless impulse failed.')
    if max(r['variance_balance_absolute'] for r in results)>2e-10:
        raise ValueError('Variance identity failed.')
    if max(r['energy_balance_absolute'] for r in results)>2e-9:
        raise ValueError('Expansion/work closure failed.')
    # Exact monoenergetic pulse followed by field-free isotropic deflections.
    rate=.7
    y=[1.,1.,1.,1.,1.,0.,0.,0.]
    for _ in range(200):
        y=step(y,0.,rate,1.,1.,.01)
    target=[math.exp(-2*rate),1.,math.exp(-2*rate),math.exp(-2*rate),1.,1-math.exp(-4*rate),0.,0.]
    delta=max(abs(a-b) for a,b in zip(target,y))
    if delta>1e-9:
        raise ValueError('Monoenergetic pulse scattering check failed.')
    # For dilute passive cold particles, |P| <= |q| TV(A)/sqrt(3) even relativistically.
    end_boost=total_variation/(math.sqrt(3)*end.a)
    rest_over_em=SEED_FRACTION*start.density*(start.a/end.a)**3/end.sheet
    def scale_for(share):
        x=share/rest_over_em
        return math.sqrt(x*(2+x))/end_boost
    xmax=CHARGE/MASS*end_boost
    maxkin=rest_over_em*xmax*xmax/(math.hypot(1,xmax)+1)
    return dict(start_time=start.t,turn_time=end.t,start_radius=start.a,turn_radius=end.a,
                start_potential=start.potential,turn_potential=end.potential,
                seed_rest_fraction_of_start_surface_matter=SEED_FRACTION,mass=MASS,charge_norm=CHARGE,
                projected_charge=PROJECTED,start_surface_density=start.density,
                turn_EM_density=end.sheet,rest_over_EM_at_turn=rest_over_em,
                potential_variation=total_variation,integral_abs_aX=integral,
                monoenergetic_scatter_control_max_error=delta,
                turn_max_p_over_m=xmax,all_rate_relativistic_kinetic_over_EM_bound=maxkin,
                charge_to_mass_needed_for_1pct=scale_for(.01),
                charge_to_mass_needed_for_50pct=scale_for(.5),
                charge_to_mass_needed_for_90pct=scale_for(.9))


def trajectory_controls(rows):
    # The integrated test current is bounded without choosing a scattering rate.
    start=rows[0]
    n0=SEED_FRACTION*start.density/MASS
    impulse=charge_current=old_current=peak_displacement=0.0
    for x,z in zip(rows[:-1],rows[1:]):
        impulse+=abs(PROJECTED*(z.potential-x.potential))
        n=n0*(start.a/z.a)**3
        vmax=impulse/math.hypot(MASS*z.a,impulse)
        current=z.a*z.a*CHARGE*n*vmax
        charge_current+=(z.t-x.t)*(old_current+current)/2
        old_current=current
        peak_displacement=max(peak_displacement,abs(z.a*z.K*z.adot_potential)/math.sqrt(3))
    return dict(wall_current_comoving_integral_bound=charge_current,
                peak_projected_background_electric_displacement=peak_displacement,
                ratio_of_impulse_bound_to_peak_displacement=charge_current/peak_displacement)


def particle_control(rows, reference, seed=20261120, samples=65536):
    # Draw the independent random flights directly. A magnetic field rotates the
    # transverse plane and therefore leaves all the following observables unchanged.
    rng=random.Random(seed)
    times=[row.t for row in rows]
    start,end=rows[0],rows[-1]
    def potential(t):
        i=min(len(times)-2,max(0,bisect.bisect_right(times,t)-1))
        x,z=rows[i],rows[i+1]
        r=(t-x.t)/(z.t-x.t)
        return x.potential+r*(z.potential-x.potential)
    summed=[0.0]*6
    untouched=0
    rate=reference['rate']
    for _ in range(samples):
        t=start.t
        prev=start.potential
        par=perp2=0.0
        count=0
        while True:
            next_t=t-math.log1p(-rng.random())/rate
            if next_t>=end.t:
                par+=PROJECTED*(end.potential-prev)
                break
            current=potential(next_t)
            par+=PROJECTED*(current-prev)
            mag2=par*par+perp2
            par=math.sqrt(mag2)*(2*rng.random()-1)
            perp2=max(0.,mag2-par*par)
            prev=current
            t=next_t
            count+=1
        untouched+=int(count==0)
        mag2=par*par+perp2
        p2=mag2/end.a**2
        exact_energy=p2/(math.sqrt(MASS*MASS+p2)+MASS)
        for i,v in enumerate((par,mag2,par*par,mag2*par,mag2*mag2,exact_energy)):
            summed[i]+=v
    u,s,par2,c,fourth,exact=[x/samples for x in summed]
    var=s-u*u
    kurtosis=3*(fourth-4*u*c+2*u*u*s+4*u*u*par2-3*u**4)/(5*var*var)
    random_energy=var/(2*MASS*end.a**2)
    differences=dict(random_energy=abs(random_energy/reference['random_energy_per_carrier']-1),
                     kurtosis=abs(kurtosis/reference['central_fourth_moment_over_Maxwell_value']-1),
                     unscattered_absolute=abs(untouched/samples-reference['probability_no_deflection']))
    if max(differences.values())>.025:
        raise ValueError('Direct random-flight control disagreed with the exact moment evolution.')
    return dict(particles=samples,seed=seed,rate=rate,random_energy_per_carrier=random_energy,
                raw_fourth_central_over_Maxwell_value=kurtosis,unscattered_fraction=untouched/samples,
                exact_relativistic_kinetic_per_carrier=exact,
                relative_relativistic_reduction=1-exact/(s/(2*MASS*end.a**2)),
                comparison_to_moments=differences)


def calculate(path):
    rows=trajectory(path)
    results=[run(rows,nu) for nu in RATES]
    primary=add_checks(rows,results)
    # Skip eight of nine output times, retaining both interval endpoints.
    coarse=rows[::9]
    if coarse[-1]!=rows[-1]:
        coarse.append(rows[-1])
    coarse_row=run(coarse, .3)
    hi=next(r for r in results if r['rate']==.3)
    ratios={k:abs(coarse_row[k]/hi[k]-1) for k in
            ('random_energy_per_carrier','mean_directed_energy_per_carrier',
             'central_fourth_moment_over_Maxwell_value')}
    peak=max(results,key=lambda r:r['random_energy_per_carrier'])
    return dict(setup_and_bounds=primary,rate_scan=results,
                direct_particle_control=particle_control(rows,hi),
                test_current_control=trajectory_controls(rows),
                highest_random_among_tested_rates=peak['rate'],
                trace_points_used=len(rows),coarsened_trace_points=len(coarse),
                relative_change_under_trace_coarsening=ratios)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=Path(__file__).resolve().parents[1]/'background/trace.csv.gz')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.write_text(json.dumps(calculate(args.input),indent=2,allow_nan=False)+'\n')

