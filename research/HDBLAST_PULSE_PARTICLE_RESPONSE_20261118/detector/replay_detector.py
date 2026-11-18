"""Independent standard-library RK4 and Simpson replay of the screened detector."""
import argparse
from bisect import bisect_right
import cmath
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import zipfile

CUT = .15
GAP = 20.
HISTORY_SHA = "9a7069e44c379e68555511e82fe841e0f01a72b23e77fa5ec69caa1373beaf70"
SUMMARY_SHA = "bf162614a84fd207f2f2736b8b77cd1cb173c2db012469cf1d4fe1ccc50c8d0f"


def drive(raw, rho, steps):
    rows = list(csv.DictReader(io.StringIO(raw.decode())))
    t = [float(r['T']) for r in rows]
    a0 = float(rows[0]['A_wall'])
    loga = [float(r['A_wall'])-a0 for r in rows]
    m = [100*(.5-float(r['phi_wall'])) for r in rows]
    n = [math.exp(float(r['B_wall']))/rho for r in rows]
    clock = [0.]
    for j in range(len(t)-1):
        clock.append(clock[-1]+(t[j+1]-t[j])*(n[j+1]+n[j])/2)
    if clock[-1] < CUT or min(m)<=0:
        raise ValueError('Incomplete or crossing history')
    grid = []
    for i in range(2*steps+1):
        s = CUT*i/(2*steps)
        j = bisect_right(clock,s)-1
        x = (s-clock[j])/(clock[j+1]-clock[j])
        a = math.exp(loga[j]+x*(loga[j+1]-loga[j]))
        mass = m[j]+x*(m[j+1]-m[j])
        if .10<s<CUT:
            switch = math.sin(math.pi*(s-.10)/.05)**4
        else:
            switch = 0.
        # d/ds(v,pi,I) = (pi/a, -(k^2/a + a*m^2)*v, coeff*v).
        grid.append((1/a,a*mass*mass,switch*cmath.exp(-1j*GAP*s)/a))
    return grid,m[0]


def solution(grid, mass, k):
    w = math.hypot(k,mass)
    v = 1/math.sqrt(2*w)
    p = -1j*w*v
    integral = 0j
    steps = (len(grid)-1)//2
    h = CUT/steps
    def slopes(coeff, q, momentum):
        a_inv, mass_term, response = coeff
        return a_inv*momentum,-(k*k*a_inv+mass_term)*q,response*q
    for j in range(steps):
        l,mid,r = grid[2*j:2*j+3]
        v1,p1,i1 = slopes(l,v,p)
        v2,p2,i2 = slopes(mid,v+h*v1/2,p+h*p1/2)
        v3,p3,i3 = slopes(mid,v+h*v2/2,p+h*p2/2)
        v4,p4,i4 = slopes(r,v+h*v3,p+h*p3)
        v += h*(v1+2*v2+2*v3+v4)/6
        p += h*(p1+2*p2+2*p3+p4)/6
        integral += h*(i1+2*i2+2*i3+i4)/6
    return integral,v,p


def integrate_band(raw,rho,steps,count):
    if count%2 or count<=0:
        raise ValueError('Simpson rule requires a positive even interval count')
    grid,mass=drive(raw,rho,steps)
    dk=13/count
    terms=[]
    for j in range(count+1):
        k=5+j*dk
        response=solution(grid,mass,k)[0]
        weight=1 if j==0 or j==count else (4 if j%2 else 2)
        terms.append(weight*k*k*abs(response)**2)
    _,v,p=solution(grid,mass,10.)
    return dk*math.fsum(terms)/(6*math.pi**2),v,p


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive',type=Path)
    parser.add_argument('summary',type=Path)
    args=parser.parse_args()
    data=args.summary.read_bytes()
    if hashlib.sha256(data).hexdigest()!=SUMMARY_SHA or hashlib.sha256(args.archive.read_bytes()).hexdigest()!=HISTORY_SHA:
        raise ValueError('Wrong source revision')
    rho=json.loads(data)['rho_b']
    cases={'quiet':'runs/g8/quiet','plus':'runs/g8/delayed','minus':'runs/g8/delayed_negative',
           'half_plus':'half_runs/g8/delayed_positive','half_minus':'half_runs/g8/delayed_negative'}
    out={}
    with zipfile.ZipFile(args.archive) as archive:
        for steps,n in ((600,32),(1200,64)):
            records={name:integrate_band(archive.read(path+'.csv'),rho,steps,n) for name,path in cases.items()}
            quiet,vq,pq=records['quiet']
            result={}
            for name,(prob,v,p) in records.items():
                result[name]={'probability_per_coupling_squared':prob,'difference':prob-quiet,
                              'relative_difference':(prob-quiet)/quiet,'relative_squeeze_k10':abs(vq*p-pq*v)**2}
            if not 5.30e-7<result['plus']['difference']<5.34e-7 or not -.0031<result['minus']['relative_difference']<-.0030:
                raise ValueError('Detector response outside published numerical range')
            out[f's{steps}_k{n}']=result
    print(json.dumps(out,indent=2,allow_nan=False))


if __name__=='__main__':
    main()
