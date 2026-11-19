"""Independent frozen Gaussian channels and exact step-integrated radiation for fixed shell histories."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import zipfile

import numpy as np
from scipy.interpolate import CubicHermiteSpline

WALL = Path('research/HDBLAST_SOURCE_WORK_20261118/coupled/WALLS.zip')
INPUTS = Path('research/HDBLAST_PULSE_PARTICLE_RESPONSE_20261118/classical/runs')
WALL_SHA = '7f3289631a80291f6180e51f5b648a135c6e33c6366d4f97c32a4b90630cab15'
SUMMARY_SHA = 'abec0887ee4de00ba0757e7d96987610f36b7157a85bc4f00696fee8064a6810'
END = 2.10


def checked(path, digest):
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != digest:
        raise ValueError(f'Unexpected source file: {path}')
    return raw


def history(z, arm, origin):
    prefix = 'fine/y1_b1/'
    rows = []
    for name in ('pre', arm):
        with z.open(prefix + name + '.csv') as stream:
            rows.append(np.genfromtxt(stream, delimiter=',', names=True))
    data = np.concatenate(rows)
    s = data['H0tau']
    if not np.all(np.diff(s)>0) or not s[0] < END <= s[-1]:
        raise ValueError('Proper clock or final time')
    ln_a = CubicHermiteSpline(s, data['A_b'] - origin, data['H_over_H0'])
    phi = CubicHermiteSpline(s, data['phi_b'], data['v_over_H0'])
    roots = phi.solve(.5,extrapolate=False)
    roots = roots[(roots > s[0]) & (roots < END)]
    if len(roots)!=1:
        raise ValueError('Expected one mass crossing')
    return ln_a,phi,float(s[0]),float(roots[0])


def propagate(curves, bath, nodes, subdivisions, preparation):
    la, phi, start, crossing = curves
    # Concentrate time resolution near the narrow low-momentum crossing.
    edges = (start, crossing-.035, crossing+.035, END)
    segments = [np.linspace(l,r,1+int(np.ceil((r-l)*subdivisions*n)))
                for l,r,n in zip(edges,edges[1:],(1,8,1))]
    s = np.concatenate([segments[0],segments[1][1:],segments[2][1:]])
    mid = (s[1:]+s[:-1])/2
    dt = np.diff(s)
    a = np.exp(la(mid))
    h = la(mid,1)
    mu = 100*(phi(mid)-.5)
    dmu = 100*phi(mid,1)
    nodes_std, weights = np.polynomial.legendre.leggauss(nodes)
    k = 4*(nodes_std+1)
    nu = 4*weights*k*k/(2*np.pi**2)
    a0 = np.exp(la(start))
    e0 = np.hypot(k/a0,100*(phi(start)-.5))
    M0 = a0**3*e0
    mu0 = 100*(phi(start)-.5)
    h0 = la(start,1)
    dmu0 = 100*phi(start,1)
    j0 = (h0+(h0*mu0**2+mu0*dmu0)/(2*e0**2))/e0 if preparation=='B' else np.zeros_like(k)
    qq = .5/M0
    pp = M0*(1+j0*j0)/2
    r = -j0/2
    radiation = np.zeros(2)
    least_redshifted = 0.
    time_least = start
    least_density = 0.
    time_least_density = start
    negative_loss = 0.
    min_count_pulse = float("inf")
    B_integral = 0.
    for i, step in enumerate(dt):
        E = np.hypot(k/a[i],mu[i])
        M = a[i]**3*E
        n = .5*(M*qq+pp/M-1)
        z = .5*(M*qq-pp/M)+1j*r
        gamma = mu[i]**2/(8*np.pi*E)
        D = h[i] + (h[i]*mu[i]**2+mu[i]*dmu[i])/(2*E**2)
        j = D/E if bath == 'B' else np.zeros_like(k)
        nt = j*j/4
        zt = -nt - .5j*j
        alpha = gamma + 2j*E
        settled = gamma*zt/alpha
        loss = -np.expm1(-gamma*step)
        rotate = -np.expm1(-alpha*step)
        # Integrals of gamma*n_H and gamma*z over a frozen step.
        gn = gamma*step*nt + (n-nt)*loss
        gz = gamma*(step*settled + (z-settled)*rotate/alpha)
        signed = (n-nt)*loss
        counted = gn + j*gz.imag + .5*j*j*(gn+.5*gamma*step+gz.real)
        incr = a[i]*np.array([nu@(E*signed),nu@(E*counted)])
        radiation += incr
        min_count_pulse = min(min_count_pulse,float(incr[1]))
        negative_loss += min(0.,float(incr[0]))
        density = radiation[0]/np.exp(4*la(s[i+1]))
        if density < least_density:
            least_density = float(density)
            time_least_density = float(s[i+1])
        B_integral += (a[i]*nu) @ (E*nt*gamma)*step
        if radiation[0] < least_redshifted:
            least_redshifted = float(radiation[0])
            time_least = float(s[i+1])
        n += (nt-n)*loss
        z += (settled-z)*rotate
        qq = (n+.5+z.real)/M
        pp = M*(n+.5-z.real)
        r = z.imag
    af = float(np.exp(la(END)))
    hf, mf, mf_dot = float(la(END,1)),float(100*(phi(END)-.5)),float(100*phi(END,1))
    ef = np.hypot(k/af,mf)
    M = af**3*ef
    n = .5*(M*qq+pp/M-1)
    c = .5*(M*qq-pp/M)
    jf = (hf+(hf*mf*mf+mf*mf_dot)/(2*ef**2))/ef if bath=='B' else 0.
    gamma = mf**2/(8*np.pi*ef)
    nb = n + jf*r + .5*jf*jf*(n+.5+c)
    if min_count_pulse < -1e-12:
        raise ValueError('Counted deposition decreased')
    if min(n)<-5e-9 or min(nb)<-5e-9 or min(qq*pp-r*r)<.25-5e-9:
        raise ValueError('Negative occupation or uncertainty violation')
    return dict(steps=len(dt),nodes=nodes,cut_s=END,scale_end=af,start_s=start,
                prep=preparation, bath=bath,scale_start=float(a0),min_R_H=least_density,
                s_min_R_H=time_least_density,integrated_negative_a4_Q_H=negative_loss,
                crossing_s=crossing, R_H=float(radiation[0]/af**4),R_count=float(radiation[1]/af**4),
                comoving_number_H=float(nu@n),comoving_number_count=float(nu@nb),
                Q_H_final=float(nu@(gamma*ef*(n-.25*jf*jf))/af**3),
                Q_count_final=float(nu@(gamma*ef*nb)/af**3),
                min_a4_R_H=least_redshifted,s_min_a4_R_H=time_least,
                target_H_ground_energy_offset_radiation=float(B_integral/af**4))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('repo',type=Path)
    args = p.parse_args()
    origin = json.loads(checked(args.repo/INPUTS/'fine/quiet/rolling_summary.json',SUMMARY_SHA))['matter']['events'][0]['A_star']
    with zipfile.ZipFile(io.BytesIO(checked(args.repo/WALL,WALL_SHA))) as z:
        result={}
        for case, nodes, count in [('step',24,2500),('time_refined',24,5000),('joint_refined',48,10000)]:
            result[case]={}
            for arm in ('quiet','plus'):
                curves=history(z,arm,origin)
                for bath, preparation, suffix in [('H','H','H'),('B','H','B'),('B','B','B_prepared')]:
                    result[case][f'{arm}_{suffix}']=propagate(curves,bath,nodes,count,preparation)
        print(json.dumps(result,indent=2,allow_nan=False))


if __name__=='__main__':
    main()
