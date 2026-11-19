"""Independent moving-parent tests with canonical transfer and direct Fourier integration."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import zipfile

import numpy as np
from scipy.integrate import cumulative_simpson, simpson
from scipy.interpolate import CubicHermiteSpline, PchipInterpolator


WALL = Path('research/HDBLAST_SOURCE_WORK_20261118/coupled/WALLS.zip')
REF = Path('research/HDBLAST_PULSE_PARTICLE_RESPONSE_20261118/classical/runs/fine/quiet/rolling_summary.json')
WALL_SHA = '7f3289631a80291f6180e51f5b648a135c6e33c6366d4f97c32a4b90630cab15'
REF_SHA = 'abec0887ee4de00ba0757e7d96987610f36b7157a85bc4f00696fee8064a6810'
WINDOWS = {'wide': (1.15, 2.05), 'narrow': (1.25, 1.95)}
CUT = 2.10
NAMES = ('annihilation_number','stimulated_number','annihilation_energy','stimulated_energy')


def pinned(path, digest):
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != digest:
        raise ValueError(f'Historical input changed: {path}')
    return data


def historical(archive, arm, ref):
    samples = []
    for part in ('pre', arm):
        with archive.open(f'fine/y1_b1/{part}.csv') as stream:
            samples.append(np.genfromtxt(stream, delimiter=',', names=True))
    data = np.concatenate(samples)
    s = data['H0tau']
    if np.any(np.diff(s) <= 0) or not s[0] < CUT <= s[-1]:
        raise ValueError('Proper clock')
    loga = CubicHermiteSpline(s, data['A_b'] - ref, data['H_over_H0'])
    mass = CubicHermiteSpline(s, 100*(data['phi_b']-.5), 100*data['v_over_H0'])
    grid = np.linspace(s[0], CUT, 32769)
    times = cumulative_simpson(np.exp(-loga(grid)), x=grid, initial=0.)
    return PchipInterpolator(times, grid), loga, mass, times[-1]


def modes(background, q, steps):
    clock, loga, mass, elapsed = background
    eta = np.linspace(0, elapsed, steps+1)
    s = clock(eta)
    amid = np.exp(loga(clock((eta[:-1]+eta[1:])/2)))
    mmid = mass(clock((eta[:-1]+eta[1:])/2))
    a = np.exp(loga(s))
    omega = np.sqrt((amid*mmid)[:,None]**2 + q**2)
    step = eta[1]-eta[0]
    cosine = np.cos(omega*step)
    # Exact harmonic transfer with (a,m) fixed at each conformal midpoint.
    sine = step*np.sinc(omega*step/np.pi)
    upper = sine / amid[:,None]**2
    lower = -sine*(omega*amid[:,None])**2
    chi = np.empty((len(eta),len(q)), dtype=complex)
    pi = np.empty_like(chi)
    inertia = a[0]**3*np.hypot(q/a[0], mass(s[0]))
    chi[0] = 1/np.sqrt(2*inertia)
    pi[0] = -1j*inertia*chi[0]
    for i in range(steps):
        chi[i+1] = cosine[i]*chi[i] + upper[i]*pi[i]
        pi[i+1] = lower[i]*chi[i] + cosine[i]*pi[i]
    error = float(np.max(abs(1j*(chi.conj()*pi-pi.conj()*chi)-1)))
    if error > 2e-9:
        raise ValueError('Canonical normalization')
    return eta, s, a, a[:,None]*chi, (a*a*loga(s,1))[:,None]*chi+pi/a[:,None], error


def nodes(lo, hi, order):
    x, w = np.polynomial.legendre.leggauss(order)
    return lo+(x+1)*(hi-lo)/2, w*((hi-lo)/2)


def integrate_band(eta, g, q, af, w, weights):
    # Direct positive and negative Fourier integrals. Our support excludes
    # the conformal-grid endpoints, so their trapezoid weights vanish.
    phase = np.exp(1j*np.outer(w, eta))
    dt = eta[1]-eta[0]
    if q.size == 1 and g.shape[1] == 1:
        # Avoid multithreaded matrix-vector calls on the low-frequency panels.
        plus = np.sum(phase*g[:,0],axis=1)[:,None]
        minus = np.sum(phase.conj()*g[:,0],axis=1)[:,None]
    else:
        plus, minus = phase@g, phase.conj()@g
    measure = weights[:,None]*(w[:,None]**2-q**2) * (dt**2/(8*np.pi**2))
    def moment(F, energy):
        integrand = measure*abs(F)**2
        if energy:
            integrand *= (w/af)[:,None]
        return np.sum(integrand,axis=0)
    return np.vstack((moment(plus,False),moment(minus,False),moment(plus,True),moment(minus,True)))


def coefficients(spectrum):
    values = {k:float(v) for k,v in zip(NAMES,spectrum,strict=True)}
    values['number'] = values['annihilation_number'] + values['stimulated_number']
    values['energy'] = values['annihilation_energy'] + values['stimulated_energy']
    values['net_loss'] = values['annihilation_number'] - values['stimulated_number']
    return values


def window(curves, q, limits, spectrum_steps, frequency_order):
    eta, s, a, v, vp, _ = curves
    lo, hi = limits
    arg = np.pi*(s-lo)/(hi-lo)
    f = np.where((s>=lo)&(s<=hi), np.sin(arg)**4, 0.)
    df = np.where((s>=lo)&(s<=hi), 4*np.pi/(hi-lo)*np.sin(arg)**3*np.cos(arg), 0.)
    g = f[:,None]*v
    gp = (df*a)[:,None]*v+f[:,None]*vp
    signed_time = simpson(abs(gp)**2-(q**2)*abs(g)**2,x=eta,axis=0)/(4*np.pi)
    stride = (len(s)-1)//spectrum_steps
    if spectrum_steps*stride!=len(s)-1:
        raise ValueError('Spectral sampling does not divide the clock')
    # Restrict time to the nonzero switch. Retain one sample outside either end.
    active = np.flatnonzero(f[::stride]>0)
    sl = slice(max(0,active[0]-1),min(spectrum_steps+1,active[-1]+2))
    eta = eta[::stride][sl]
    g = g[::stride][sl]
    thresholds = np.zeros((4,len(q)))
    space = np.zeros(len(q))
    for i, k in enumerate(q):
        if k < 8:
            w,c = nodes(k,8,frequency_order)
            thresholds[:,i] = integrate_band(eta,g[:,i:i+1],q[i:i+1],a[-1],w,c)[:,0]
        if k > 0:
            w,c = nodes(0,k,frequency_order)
            below = integrate_band(eta,g[:,i:i+1],q[i:i+1],a[-1],w,c)
            space[i] = -np.sum(below[:2])
    def high(lo,hi):
        # Common absolute-frequency nodes amortize a matrix multiplication across q.
        left = np.arange(lo,hi,16.)
        panels = [nodes(x,min(x+16,hi),frequency_order) for x in left]
        w,c = (np.concatenate(cols) for cols in zip(*panels))
        return integrate_band(eta,g,q,a[-1],w,c)
    screen = thresholds+high(8.,60.)
    mid = screen+high(60.,192.)
    first = mid+high(192.,384.)
    full = first+high(384.,768.)
    time_res = float(np.max(abs(full[0]+full[1]-signed_time-space)))
    tail = np.sum(full[2:]-first[2:],axis=0)/np.sum(full[2:],axis=0)
    if time_res>4e-6 or max(tail)>2e-8 or np.any(full < -1e-12):
        raise ValueError(('Fourier and time comparison',time_res,max(tail)))
    out = {}
    for name,arr in (('screen',screen),('full',full)):
        averages = {
            'uniform': 3*simpson(arr*q**2,x=q,axis=1)/8**3,
            'tapered': (105/(8*8**3))*simpson(arr*q**2*(1-(q/8)**2)**2,x=q,axis=1),
        }
        out[name] = {key:coefficients(value) for key,value in averages.items()}
        for k in (0,2,4,6,8):
            index = np.flatnonzero(q==k).item()
            out[name][f'q{k}'] = coefficients(arr[:,index])
    return dict(data=out, max_number_parseval_error=time_res,
                max_fraction_energy_384_to_768=float(max(tail)))


def run(repo,configs):
    ref = json.loads(pinned(repo/REF,REF_SHA))['matter']['events'][0]['A_star']
    result = {}
    with zipfile.ZipFile(io.BytesIO(pinned(repo/WALL,WALL_SHA))) as archive:
        for steps, parent_intervals, freq_order, spectrum_steps in configs:
            q = np.linspace(0,8,parent_intervals+1)
            arms = {}
            for arm in ('quiet','plus'):
                curves = modes(historical(archive,arm,ref),q,steps)
                arms[arm] = dict(wronskian=curves[-1], a_cut=float(curves[2][-1]), windows={
                    name: window(curves,q,limits,spectrum_steps,freq_order)
                    for name,limits in WINDOWS.items()})
            contrasts = {}
            for name in WINDOWS:
                quiet = arms['quiet']['windows'][name]['data']
                plus = arms['plus']['windows'][name]['data']
                contrasts[name] = {band:{state:
                    {k:100*(plus[band][state][k]/v-1) for k,v in row.items()}
                    for state,row in states.items()} for band,states in quiet.items()}
            density = {}
            for arm, results in arms.items():
                density[arm] = {}
                for name, w in results['windows'].items():
                    density[arm][name] = {}
                    for band, states in w['data'].items():
                        density[arm][name][band] = {}
                        for state, count_integral in (('uniform',8**3/3),('tapered',8**3*8/105)):
                            factor = count_integral/(2*np.pi**2*results['a_cut']**3)
                            density[arm][name][band][state] = {
                                k:factor*states[state][k] for k in ('number','energy')}
            result[f'{steps}_{parent_intervals}'] = dict(steps=steps,parent_intervals=parent_intervals,
                frequency_order=freq_order,spectrum_steps=spectrum_steps,arms=arms,
                contrasts_percent=contrasts, physical_densities=density)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('repository',type=Path)
    parser.add_argument('--quick',action='store_true')
    opt=parser.parse_args()
    settings=[(4096,16,10,1024),(8192,32,12,2048),(16384,64,14,2048),(32768,128,16,2048)]
    print(json.dumps(run(opt.repository,settings[:1] if opt.quick else settings),indent=2,allow_nan=False))
