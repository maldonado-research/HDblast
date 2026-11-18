"""In/out evolution of the externally returned mass profile in the stationary test."""
import cmath
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent


def spline(x):
    x=abs(x)
    if x>=2:
        return 0.
    if x>=1:
        return (2-x)**3/6
    return (4-6*x*x+3*x**3)/6


def solve(omega,profile,dt):
    q,v=1.+0j,-1j*omega
    energy_work=0.
    last=0.
    for potential in profile:
        energy_work+=(potential-last)*abs(q)**2/(4*omega)
        last=potential
        w=math.sqrt(omega*omega+potential)
        cs,sn=math.cos(w*dt),math.sin(w*dt)
        q,v=cs*q+sn*v/w,-w*sn*q+cs*v
    energy_work-=last*abs(q)**2/(4*omega)
    alpha=(q+1j*v/omega)/2
    beta=(q-1j*v/omega)/2
    n=abs(beta)**2
    return n,abs(abs(alpha)**2-n-1),abs(energy_work-omega*n)


def calculate(c,j,h,m,steps):
    dt=.05/steps
    vals=[math.fsum(a*spline((.1+(i+.5)*dt)/h-b) for a,b in zip(c,j)) for i in range(steps)]
    spectrum=[]
    for i in range(53):
        k=5+i*.25
        omega=math.hypot(k,m)
        n,defect,work=solve(omega,vals,dt)
        spectrum.append([k,n,defect,work])
    def integral(energy):
        ys=[k*k*n*(math.hypot(k,m) if energy else 1)/(2*math.pi**2) for k,n,*_ in spectrum]
        return .25/3*(ys[0]+ys[-1]+4*math.fsum(ys[1:-1:2])+2*math.fsum(ys[2:-1:2]))
    return dict(spectrum=[[k,n] for k,n,*_ in spectrum],band_number=integral(False),band_energy=integral(True),
                max_wronskian_error=max(r[2] for r in spectrum),max_energy_work_error=max(r[3] for r in spectrum))


def main():
    d=json.loads((HERE/'PULSES.json').read_text())
    bounds=json.loads((HERE/'BOUNDS.json').read_text())
    m=100*(.5-float(d['phi0']))
    # A square mass step has an analytic occupation even when it is large.
    w,potential,width=51.,-7.,.025
    exact=(potential*math.sin(math.sqrt(w*w+potential)*width)/(2*w*math.sqrt(w*w+potential)))**2
    tested=solve(w,[potential],width)[0]
    if abs(tested-exact)>1e-16:
        raise ValueError('Square-step benchmark failed')
    out={}
    for name,coeff in d['protocols'].items():
        c=list(map(float,coeff))
        coarse,fine=(calculate(c,d['centers_in_h'],float(d['h']),m,n) for n in (400,800))
        fine['max_two_grid_occupation_gap']=max(abs(a[1]-b[1]) for a,b in zip(coarse['spectrum'],fine['spectrum']))
        b=bounds[name]
        if not all(float(b['occupation_lower'])<=n<=float(b['occupation_upper']) for k,n in fine['spectrum']):
            raise ValueError('Incompatible with certified entire-band bound')
        out[name]=fine
    print(json.dumps(dict(square_step=dict(theory=exact,numerical=tested),protocols=out),indent=2))


if __name__=='__main__':
    main()
