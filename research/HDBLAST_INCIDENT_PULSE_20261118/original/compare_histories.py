"""Compare the published RK4 implementation with the shrinking-domain discretization."""
import argparse
from bisect import bisect_right
import csv
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent


def read(path):
    with path.open(newline='') as f:
        return [{k:float(v) for k,v in row.items()} for row in csv.DictReader(f)]


def interp(t,y,x):
    if x<t[0]-1e-12 or x>t[-1]+1e-12:
        raise ValueError('Uncovered comparison interval')
    if x>=t[-1]:
        return y[-1]
    j=max(0,bisect_right(t,x)-1)
    return y[j]+(y[j+1]-y[j])*(x-t[j])/(t[j+1]-t[j])


def comparison(runs):
    result={}
    finest=runs/'g8'
    for grid in ('g1','g2','g4'):
        folder=ROOT/'runs'/grid
        p,q=(read(folder/(name+'.csv')) for name in ('delayed','quiet'))
        times=[row['T'] for row in p]
        if times!=[row['T'] for row in q] or abs(times[-1]-.18)>1e-11:
            raise ValueError('Unmatched or incomplete RK4 histories')
        y=[a['phi_b']-b['phi_b'] for a,b in zip(p,q)]
        if any(row[key] != 0 for row in p+q for key in ('R','X')):
            raise ValueError('Nonzero material density')
        pulse_peak=max(y)
        stats={'rows':len(times),'peak_delta':pulse_peak,'peak_T':times[y.index(pulse_peak)],
               'max_phi':max(row['phi_b'] for row in p)}
        for label,other in (('same_grid',runs/grid),('finest_local_grid',finest)):
            pa,qa=(read(other/(name+'.csv')) for name in ('delayed','quiet'))
            t=[a['T'] for a in pa]
            local=[a['phi_wall']-b['phi_wall'] for a,b in zip(pa,qa)]
            gaps=[abs(v-interp(t,local,x)) for x,v in zip(times,y)]
            stats[label]={'max_abs_gap':max(gaps),'relative_to_peak':max(gaps)/abs(pulse_peak)}
        result[grid]=stats
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('runs',type=Path)
    print(json.dumps(comparison(p.parse_args().runs),indent=2))
