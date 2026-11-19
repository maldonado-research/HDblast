"""Independent tagged-parent number and energy in a switched conformal Dirac channel."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import zipfile

import numpy as np
from scipy.integrate import cumulative_simpson, simpson
from scipy.interpolate import CubicHermiteSpline, PchipInterpolator

ARCHIVE = Path('research/HDBLAST_SOURCE_WORK_20261118/coupled/WALLS.zip')
SUMMARY = Path('research/HDBLAST_PULSE_PARTICLE_RESPONSE_20261118/classical/runs/fine/quiet/rolling_summary.json')
SHA_ARCHIVE = '7f3289631a80291f6180e51f5b648a135c6e33c6366d4f97c32a4b90630cab15'
SHA_SUMMARY = 'abec0887ee4de00ba0757e7d96987610f36b7157a85bc4f00696fee8064a6810'
WINDOWS = {'wide': (1.15,2.05), 'narrow':(1.25,1.95)}
CUT = 2.10
M5 = 6.49e-6**(-1/3)


def pin(path, digest):
    data=path.read_bytes()
    if hashlib.sha256(data).hexdigest()!=digest:
        raise ValueError(path)
    return data


def wall(z, arm, origin):
    rows=[]
    for name in ('pre',arm):
        with z.open('fine/y1_b1/'+name+'.csv') as stream:
            rows.append(np.genfromtxt(stream,delimiter=',',names=True))
    rec=np.concatenate(rows)
    s=rec['H0tau']
    if not np.all(np.diff(s)>0):
        raise ValueError('Clock')
    return (s[0], CubicHermiteSpline(s,rec['A_b']-origin,rec['H_over_H0']),
            CubicHermiteSpline(s,rec['phi_b'],rec['v_over_H0']))


def transfer(curves, steps):
    first, logscale, field=curves
    # Invert the integral d eta = ds/a separately from solving the field equation.
    ss=np.linspace(first,CUT,32769)
    clock=cumulative_simpson(np.exp(-logscale(ss)),x=ss,initial=0.)
    inverse=PchipInterpolator(clock,ss)
    eta=np.linspace(0.,clock[-1],steps+1)
    s=inverse(eta)
    am=np.exp(logscale(inverse((eta[1:]+eta[:-1])/2)))
    mm=100*(field(inverse((eta[1:]+eta[:-1])/2))-.5)
    dt=eta[1]-eta[0]
    phase=am*mm*dt
    cs=np.cos(phase)
    time=dt*np.sinc(phase/np.pi)
    upper=time/am**2
    lower=-time*am**4*mm**2
    a=np.exp(logscale(s))
    mu0=100*abs(field(first)-.5)
    q=np.empty(steps+1,dtype=complex)
    p=np.empty_like(q)
    q[0]=1/np.sqrt(2*a[0]**3*mu0)
    p[0]=-1j*a[0]**3*mu0*q[0]
    for i in range(steps):
        q[i+1]=cs[i]*q[i]+upper[i]*p[i]
        p[i+1]=lower[i]*q[i]+cs[i]*p[i]
    error=float(np.max(abs(1j*(q.conj()*p-p.conj()*q)-1)))
    if error>2e-9:
        raise ValueError('Mode normalization')
    return eta,s,a,logscale(s,1),q,p,error


def calculate(curves, steps, padding):
    eta,s,a,h,chi,p,wr=transfer(curves,steps)
    out={}
    dt=eta[1]-eta[0]
    for name,(lo,hi) in WINDOWS.items():
        z=np.pi*(s-lo)/(hi-lo)
        active=(s>=lo)&(s<=hi)
        f=np.where(active,np.sin(z)**4,0.)
        fd=np.where(active,(4*np.pi/(hi-lo))*np.sin(z)**3*np.cos(z),0.)
        v=a*chi
        dv=a**2*h*chi+p/a
        g=f*v
        dfv=fd*a*v
        fdv=f*dv
        norm=lambda y:float(simpson(abs(y)**2,x=eta)/(4*np.pi))
        num=norm(dfv+fdv)
        n_switch,n_parent=norm(dfv),norm(fdv)
        out[name]=dict(local_number=num,switch_number=n_switch,parent_number=n_parent,
                       interference=num-n_switch-n_parent,spectrum={})
        for factor in padding:
            n_fft=factor*steps
            dw=(2*np.pi)/(n_fft*dt)
            omega=2*np.pi*np.fft.fftfreq(n_fft,d=dt)
            mag=abs(dt*np.fft.fft(g,n=n_fft))**2
            ndens=omega**2*mag*dw/(8*np.pi**2)
            edens=ndens*abs(omega)/a[-1]
            nf=float(ndens.sum())
            energy=float(edens.sum())
            if abs(nf/num-1)>2e-5:
                raise ValueError(('Parseval disagreement',nf,num))
            # Bound refers to individual massless daughters; the total pair energy is 2k/a.
            limit=2*np.exp(curves[1](lo))*M5
            high=abs(omega)>limit
            out[name]['spectrum'][str(factor)]=dict(number=nf,energy=energy,mean_pair_energy=energy/nf,
                   delta_omega=float(dw),
                   negative_frequency_number=float(ndens[omega<0].sum()),
                   positive_frequency_number=float(ndens[omega>0].sum()),
                   negative_frequency_energy=float(edens[omega<0].sum()),
                   positive_frequency_energy=float(edens[omega>0].sum()),
                   over_screen_number_fraction=float(ndens[high].sum()/nf),
                   over_screen_energy_fraction=float(edens[high].sum()/energy))
    return dict(steps=steps,a_cut=float(a[-1]),eta_span=float(eta[-1]),wronskian=wr,windows=out)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('repository',type=Path)
    args=parser.parse_args()
    root=args.repository
    origin=json.loads(pin(root/SUMMARY,SHA_SUMMARY))['matter']['events'][0]['A_star']
    with zipfile.ZipFile(io.BytesIO(pin(root/ARCHIVE,SHA_ARCHIVE))) as z:
        result={}
        for steps in (4096,8192,16384,32768):
            arms={a:calculate(wall(z,a,origin),steps,(4,8,16,32)) for a in ('quiet','plus')}
            contrasts={}
            for window in WINDOWS:
                pw=arms['plus']['windows'][window]
                qw=arms['quiet']['windows'][window]
                contrasts[window]=dict(number_percent=100*(pw['local_number']/qw['local_number']-1),
                  energy_percent=100*(pw['spectrum']['32']['energy']/qw['spectrum']['32']['energy']-1))
            result[str(steps)]=dict(arms=arms,contrasts=contrasts)
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__=='__main__':
    main()
