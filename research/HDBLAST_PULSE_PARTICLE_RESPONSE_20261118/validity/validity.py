"""Audit the B1 mass scale, finite-band work envelope and local crossing durations."""
import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import zipfile

ARCHIVE_SHA256 = '9a7069e44c379e68555511e82fe841e0f01a72b23e77fa5ec69caa1373beaf70'
SUMMARY_SHA256 = 'bf162614a84fd207f2f2736b8b77cd1cb173c2db012469cf1d4fe1ccc50c8d0f'
PRIMARY_CUT = 0.15
BAND = (5.0, 18.0)
NEAR_CROSSING = 0.0012


def authenticate(path, expected):
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        raise ValueError('Unexpected numerical input: '+path.name)


def trapezoid(t, y):
    return math.fsum((b-a)*(u+v)/2 for a,b,u,v in zip(t,t[1:],y,y[1:]))


def tension(phi, parameters):
    return (2*(1-phi+phi**3/3) + parameters['delta']*(
        1+parameters['c']*phi+parameters['d']*phi**2/2))


def interpolate(a,b,f):
    return {k:a[k]+f*(b[k]-a[k]) for k in a}


FIELDS = [('A_wall', 'A_T_wall'), ('B_wall', 'B_T_wall'), ('phi_wall', 'phi_T_wall')]


def controls(left, right, field, slope):
    dt = right['T'] - left['T']
    return (left[field], left[field] + dt * left[slope] / 3,
            right[field] - dt * right[slope] / 3, right[field])


def hermite(left, right, fraction):
    dt = right['T'] - left['T']
    row = {'T': left['T'] + fraction * dt}
    for field, slope in FIELDS:
        secant = (right[field] - left[field]) / dt
        c2 = 3 * secant - 2 * left[slope] - right[slope]
        c3 = left[slope] + right[slope] - 2 * secant
        row[field] = left[field] + dt * fraction * (left[slope] + fraction * (c2 + fraction*c3))
        row[slope] = left[slope] + fraction * (2*c2 + 3*fraction*c3)
    return row


def elapsed(left, right, fraction, rho):
    # Two-point Gauss quadrature on the same cubic Hermite lapse used for the modes.
    acc = 0.
    for node in (0.5 - math.sqrt(3)/6, 0.5 + math.sqrt(3)/6):
        acc += math.exp(hermite(left, right, fraction*node)['B_wall']) / rho
    return fraction * (right['T']-left['T']) * acc/2


def series(archive, name, rho):
    with archive.open(name) as member:
        rows = [{k:float(r[k]) for k in ('T','phi_wall','A_wall','B_wall',
                    'phi_T_wall','A_T_wall','B_T_wall')} for r in csv.DictReader(io.TextIOWrapper(member))]
    if rows[0]['T'] != 0 or any(b['T']<=a['T'] for a,b in zip(rows,rows[1:])):
        raise ValueError('Invalid time grid')
    first = rows[0]['A_wall']
    out, clock = [rows[0]], 0.
    out[0]['s'] = 0.
    for left, right in zip(rows, rows[1:]):
        ds = elapsed(left, right, 1., rho)
        if clock + ds >= PRIMARY_CUT:
            lo, hi = 0., 1.
            for _ in range(48):
                mid = (lo+hi)/2
                if clock + elapsed(left, right, mid, rho) < PRIMARY_CUT:
                    lo = mid
                else:
                    hi = mid
            end = hermite(left, right, (lo+hi)/2)
            end['s'] = PRIMARY_CUT
            out.append(end)
            break
        clock += ds
        right['s'] = clock
        out.append(right)
    if out[-1]['s'] != PRIMARY_CUT:
        raise ValueError('Requested cut not covered')
    for r in out:
        lapse = math.exp(r['B_wall'])/rho
        r.update(lapse=lapse, loga=r['A_wall']-first,
                 h=r['A_T_wall']/lapse, v=r['phi_T_wall']/lapse)
    return out


def certify_intervals(rows, G, gap, scale):
    # The cubic and its derivative both lie inside their Bernstein control hulls.
    initial_a = rows[0]['A_wall']
    min_slope, min_ceiling = math.inf, math.inf
    for left, right in zip(rows, rows[1:]):
        dt = right['T'] - left['T']
        ca = controls(left, right, 'A_wall', 'A_T_wall')
        cp = controls(left, right, 'phi_wall', 'phi_T_wall')
        min_expansion = min(3*(v-u)/dt for u,v in zip(ca,ca[1:]))
        max_speed = max(abs(3*(v-u)/dt) for u,v in zip(cp,cp[1:]))
        if max(cp) >= gap:
            raise ValueError('Cubic reaches the mass zero')
        min_slope = min(min_slope, min_expansion - max_speed/(gap-max(cp)))
        mass_upper = G*(gap-min(cp))
        margin_sq = scale**2-mass_upper**2
        if margin_sq <= 0:
            raise ValueError('Cubic reaches the EFT mass scale')
        # If the ceiling a*sqrt(scale^2-m^2) increases, use its left endpoint.
        # Otherwise use the conservative product of two interval lower bounds.
        if min_expansion - G*mass_upper*max_speed/margin_sq >= 0:
            floor = math.exp(left['A_wall']-initial_a)*math.sqrt(scale**2-(G*(gap-left['phi_wall']))**2)
        else:
            floor = math.exp(min(ca)-initial_a)*math.sqrt(margin_sq)
        min_ceiling = min(min_ceiling, floor)
    if min_slope <= 0:
        raise ValueError('Monotonicity of a*m was not certified')
    return min_slope, min_ceiling


def crossing_tests(rows,G,threshold):
    # Log the history at a shifted crossing as a kinematic diagnostic only.
    crossings=[]
    for a,b in zip(rows,rows[1:]):
        if (a['phi_wall']-threshold)*(b['phi_wall']-threshold)<0:
            event=interpolate(a,b,(threshold-a['phi_wall'])/(b['phi_wall']-a['phi_wall']))
            q=G*abs(event['v'])
            crossings.append(dict(s=event['s'],T=event['T'],q=q,h=event['h'],
                nominal_event_duration=1/math.sqrt(q),
                nonadiabatic_field_width=math.sqrt(q)/G,
                field_ramp_width=3*math.sqrt(q)/G,
                scale_factor_rate_times_duration=event['h']/math.sqrt(q)))
    outcome={'threshold':threshold,'crossings':crossings}
    if len(crossings)==2:
        a,b=crossings
        excursion=max(r['phi_wall']-threshold for r in rows if a['s']<=r['s']<=b['s'])
        outcome.update(crossing_separation=b['s']-a['s'],positive_side_excursion=excursion,
            first_to_next_in_event_units=(b['s']-a['s'])/a['nominal_event_duration'],
            maximum_first_outward_excursion_over_first_event_width=excursion/a['nonadiabatic_field_width'])
    return outcome


def analyze(rows, model):
    pars=model['params']
    rho=model['rho_b']
    G, gap, b, y = (pars['matter'][k] for k in ('G','phistar','b','y'))
    scale=b**(-1/3)
    masses=[G*abs(r['phi_wall']-gap) for r in rows]
    if min(masses)<=0 or max(masses)>=scale:
        raise ValueError('This noncrossing audit presupposes a regular subcutoff mass')
    # Upper allowed comoving k is a*sqrt((M5/H0)^2 - m^2/H0^2).
    scale_factors=[math.exp(r['loga']) for r in rows]
    am=[a*m for a,m in zip(scale_factors,masses)]
    slope_floor, lower_k_cap = certify_intervals(rows, G, gap, scale)
    if not all(b>a for a,b in zip(am,am[1:])):
        raise ValueError('Nodal monotonicity of a*m failed')
    results={}
    mass_rates=[math.copysign(1,r['phi_wall']-gap)*G*r['v'] for r in rows]
    for k in (0.,5.,18.,20.,40.):
        frequencies=[math.hypot(k/a,m) for m,a in zip(masses,scale_factors)]
        eps=[abs(m*(r['h']*m+dm))/omega**3
             for m,dm,r,omega in zip(masses,mass_rates,rows,frequencies)]
        initial_wkb_occupation=eps[0]**2/16
        results[str(k)] = dict(max_frequency=max(frequencies),max_frequency_over_M5=max(frequencies)/scale,
            inside_sampled_physical_frequency_cut=max(frequencies)<=scale,
            min_abs_frequency_derivative_indicator=min(eps),max_abs_frequency_derivative_indicator=max(eps),
            instantaneous_n_of_initial_WKB_tangent=initial_wkb_occupation,
            conditional_decay_optical_depth=trapezoid([r['s'] for r in rows],
                [y*y*m*m/(8*math.pi*omega) for m,omega in zip(masses,frequencies)]))
    a_final=scale_factors[-1]
    bounds={}
    for k0,k1,label in [(BAND[0],BAND[1],'primary'),(0.,lower_k_cap,'all_common_initially_admitted')]:
        omega0=math.hypot(k0,am[0])
        omegaf=math.hypot(k0,am[-1])
        max_squeeze=.5*math.log(omegaf/omega0)
        eps_initial=abs(masses[0]*(rows[0]['h']*masses[0]+mass_rates[0]))/math.hypot(k0,masses[0])**3
        initial_wkb_squeeze=math.asinh(eps_initial/4)
        choices={}
        for basis,initial in [('instantaneous_vacuum',0.),('WKB_tangent',initial_wkb_squeeze)]:
            max_n=math.sinh(initial+max_squeeze)**2
            # The more conservative constant M5 bounds energy per quantum over this band.
            density=max_n*scale*(k1**3-k0**3)/(6*math.pi**2*a_final**3)
            norm=b*density
            rad=rho*tension(rows[-1]['phi_wall'],pars)*norm/18+norm**2/36
            choices[basis]=dict(max_occupation=max_n,rho_excitation_over_H0_fourth=density,
                                hypothetical_radiation_share=rad/rows[-1]['h']**2)
        bounds[label]=dict(k_min=k0,k_max=k1,upper=choices)
    return dict(T_at_s_cut=rows[-1]['T'],s_cut=PRIMARY_CUT,
        log_scale_factor_at_cut=rows[-1]['loga'],h_final=rows[-1]['h'],
        M5_over_H0=scale,mass_min=min(masses),mass_max=max(masses),max_mass_over_M5=max(masses)/scale,
        certified_comoving_k_ceiling=lower_k_cap,
        hermite_lower_bound_d_log_am_dT=slope_floor,min_log_am_proper_slope=min(
            r['h']+dm/m for r,dm,m in zip(rows,mass_rates,masses)),
        max_mass_variation_over_mass_squared=max(abs(dm)/m**2 for dm,m in zip(mass_rates,masses)),
        modes=results,excitation_envelopes=bounds,
        shifted_threshold_kinematics=crossing_tests(rows,G,NEAR_CROSSING))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive',type=Path)
    parser.add_argument('summary',type=Path)
    args=parser.parse_args()
    authenticate(args.archive,ARCHIVE_SHA256)
    authenticate(args.summary,SUMMARY_SHA256)
    model=json.loads(args.summary.read_text())
    with zipfile.ZipFile(args.archive) as z:
        report={grid:{name:analyze(series(z,f'runs/{grid}/{name}.csv',model['rho_b']),model)
                for name in ('quiet','delayed','delayed_negative')} for grid in ('g4','g8')}
    print(json.dumps(report,indent=2,allow_nan=False))


if __name__=='__main__':
    main()
