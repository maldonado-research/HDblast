"""Integrate two-sided interface work on the published pulse-minus-quiet trajectories."""
import argparse
from bisect import bisect_right
import csv
from fractions import Fraction
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
WINDOWS = (0.10, 0.15, 0.18)
# Uniform one-sided derivative fits at the wall using strictly interior samples.
STENCILS = {
    2: (Fraction(5, 2), -4, Fraction(3, 2)),
    4: (Fraction(77, 12), Fraction(-107, 6), Fraction(39, 2),
        Fraction(-61, 6), Fraction(25, 12)),
}


def trap(t, y):
    return math.fsum((b-a)*(v+w)/2 for a,b,v,w in zip(t,t[1:],y,y[1:]))


def interpolate(t, y, cut):
    if not t[0] <= cut <= t[-1]:
        raise ValueError('Cut lies outside the saved history')
    if cut == t[-1]:
        return y[-1]
    j = bisect_right(t,cut)-1
    return y[j] + (y[j+1]-y[j])*(cut-t[j])/(t[j+1]-t[j])


def integrate(t, y, cut):
    stop = bisect_right(t,cut)
    if t[stop-1] == cut:
        return trap(t[:stop],y[:stop])
    return trap(t[:stop]+[cut], y[:stop]+[interpolate(t,y,cut)])


def tension_difference(phi, initial, model):
    delta,c,d = (model[name] for name in ('delta','c','d'))
    return (phi-initial)*(-2 + delta*c + delta*d*(phi+initial)/2
                          + 2*(phi**2+phi*initial+initial**2)/3)


def history(path, h, params):
    with path.open(newline='') as f:
        rows = [{k:float(v) for k,v in row.items()} for row in csv.DictReader(f)]
    t = [r['T'] for r in rows]
    if t[0] != 0 or not all(x < y for x,y in zip(t,t[1:])):
        raise ValueError('Invalid coordinate clock')
    a0, phi0 = rows[0]['A_wall'], rows[0]['phi_wall']
    lapse, clock = [], [0.]
    rates = {k:[] for k in ('quadratic','quartic','absolute','junction','volume_work','weighted_delta')}
    for r in rows:
        volume = math.exp(3*(r['A_wall']-a0))
        normal_factor = math.exp(-r['B_wall'])
        phi, v = r['phi_wall'],r['phi_T_wall']
        for degree,name in ((2,'quadratic'),(4,'quartic')):
            # Subtract the value at the first interior point to reduce cancellation.
            first = r['phi_minus_1']
            grad = math.fsum(w*(r[f'phi_minus_{j}']-first) for j,w in
                             enumerate(STENCILS[degree],1))/h
            rates[name].append(-2*volume*normal_factor*v*grad)
        rates['absolute'].append(abs(rates['quartic'][-1]))
        derivative = 2*(phi**2-1)+params['delta']*(params['c']+params['d']*phi)
        rates['junction'].append(volume*derivative*v)
        weighted_change = volume*tension_difference(phi,phi0,params)
        rates['weighted_delta'].append(weighted_change)
        rates['volume_work'].append(3*r['A_T_wall']*weighted_change)
        lapse.append(math.exp(r['B_wall'])/params['rho_b'])
    for i in range(1,len(t)):
        clock.append(clock[-1]+(t[i]-t[i-1])*(lapse[i]+lapse[i-1])/2)
    return t,clock,rates


def budget(t,clock,rates,cut):
    flux = integrate(t,rates['quartic'],cut)
    endpoint = interpolate(t,rates['weighted_delta'],cut) - integrate(t,rates['volume_work'],cut)
    absolute = integrate(t,rates['absolute'],cut)
    return dict(T=cut,s=interpolate(t,clock,cut),flux_into_shell=flux,
                absolute_flux=absolute,second_order_flux=integrate(t,rates['quadratic'],cut),
                junction=integrate(t,rates['junction'],cut),pressure_corrected_tension=endpoint,
                flux_minus_tension=flux-endpoint)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',nargs='?',type=Path,default=ROOT/'evolution/runs')
    parser.add_argument('--half',type=Path,default=ROOT/'evolution/half_runs',
                        help='additional grids with the half-amplitude matched pair')
    parser.add_argument('--mass',type=Path,default=ROOT/'evolution/mass_runs',
                        help='additional grids for the mass-matched family')
    args=parser.parse_args()
    model=json.loads((ROOT/'B1_summary.json').read_text())
    params=dict(model['params'],rho_b=model['rho_b'])
    out={}
    for folder in sorted(args.directory.glob('g*')):
        h=json.loads((folder/'metadata.json').read_text())['spacing']
        t,s,rates=history(folder/'delayed.csv',h,params)
        qt,qs,quiet=history(folder/'quiet.csv',h,params)
        reverse=history(folder/'delayed_negative.csv',h,params)
        out[folder.name]={'windows':[]}
        for cut in WINDOWS:
            p=budget(t,s,rates,cut)
            q=budget(qt,qs,quiet,cut)
            quiet_cut=interpolate(qs,qt,p['s']) if p['s']<=qs[-1] else None
            paired_clock_flux=None
            if quiet_cut is not None:
                paired_clock_flux=p['flux_into_shell']-integrate(qt,quiet['quartic'],quiet_cut)
            out[folder.name]['windows'].append(dict(pulse=p,quiet=q,
                excess_flux_at_common_T=p['flux_into_shell']-q['flux_into_shell'],
                excess_flux_at_common_proper_clock=paired_clock_flux))
        # Compare the two pulse signs at the very same shell proper time as the quiet cut.
        target_s=interpolate(qt,qs,0.15)
        at_quiet_cut=integrate(qt,quiet['quartic'],0.15)
        def common_excess(obs):
            time,proper,measure=obs
            return integrate(time,measure['quartic'],interpolate(proper,time,target_s))-at_quiet_cut
        positive,negative=(common_excess(obs) for obs in ((t,s,rates),reverse))
        def parity(a,b):
            return dict(positive=a,negative=b,even=(a+b)/2,odd=(a-b)/2)
        out[folder.name]['parity_0.001']=parity(positive,negative)
        half=args.half/folder.name
        if half.is_dir():
            out[folder.name]['parity_0.0005']=parity(
                common_excess(history(half/'delayed_positive.csv',h,params)),
                common_excess(history(half/'delayed_negative.csv',h,params)))
        mass=args.mass/folder.name
        if mass.is_dir():
            qm=history(mass/'quiet.csv',h,params)
            pm=history(mass/'mass_matched.csv',h,params)
            out[folder.name]['mass_family_at_T_0.15']=dict(
                pulse=budget(*pm,.15),quiet=budget(*qm,.15))
    if len(out)!=4:
        raise ValueError('Expected four paired grids')
    print(json.dumps(out,indent=2,allow_nan=False))


if __name__=='__main__':
    main()
