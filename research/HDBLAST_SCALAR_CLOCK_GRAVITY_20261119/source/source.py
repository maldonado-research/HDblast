"""Fixed-trajectory, order-y^2 scalar-clock Yukawa sources for native live RK.

All exported rho,p,j and their proper-time derivatives are hatted physical
density coefficients per lambda=y^2. Multiply by lambda and b/r_b once at the
native matter interface. These are forcing tables on archived trajectories,
not a live evolution of the causal moments on the perturbed geometry.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import time
from types import SimpleNamespace
import zipfile

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import CubicSpline, CubicHermiteSpline

HERE = Path(__file__).resolve().parent
WINDOW = (1.35, 2.05)
ARMS = ('quiet', 'plus', 'minus', 'far')
COMPONENTS = ('parent', 'daughter', 'interaction', 'counterterm')


def read_native(native, grid, arm='quiet'):
    """Native archive reader, from the preceding capsule's budget.py."""
    if arm not in ARMS:
        raise ValueError('Unknown native arm')
    native = Path(native)
    fork = json.loads((native/'FORK.json').read_bytes())[grid]
    with zipfile.ZipFile(native/'WALLS.zip') as archive:
        with archive.open(f'{grid}/retain/{arm}/wall.csv') as stream:
            wall = np.genfromtxt(stream, delimiter=',', names=True)
    s = wall['H0tau']
    loga = CubicHermiteSpline(s, wall['A_b']-fork['reference_A'],
                              wall['H_over_H0'], extrapolate=False)
    mass = CubicHermiteSpline(s, 100*(wall['phi_b']-.5),
                              100*wall['v_over_H0'], extrapolate=False)
    q, weights = map(np.array, (fork['q'],fork['weight']))
    node, dw = leggauss(len(q))
    np.testing.assert_allclose(q,4*(node+1),rtol=1e-13)
    np.testing.assert_allclose(weights,4*dw*q*q/(2*np.pi**2),rtol=2e-12)
    np.testing.assert_allclose([s[0],np.exp(loga(s[0]))],[fork['H0tau'],fork['a']],rtol=2e-13)
    states = {name:tuple(np.asarray(fork['states'][name][key]) for key in ('X','C','P'))
              for name in ('retain','erase')}
    inertia = fork['a']**3*np.hypot(q/fork['a'],mass(s[0]))
    states['fork_vacuum'] = (1/(2*inertia),np.zeros_like(q),inertia/2)
    return q,weights,states,(float(s[0]),loga,mass)


def budget_module():
    """Compatibility helper used by the local independent verification."""
    return SimpleNamespace(read_native=read_native)


class ScalarClock:
    """One C2 inverse, calibrated on fine quiet and reused on every grid.

    A cubic spline interpolates s(phi) at 8193 equally spaced fine-quiet
    proper times. F(phi)=sin^4(pi*(s(phi)-lo)/(hi-lo)) inside its field support,
    and is exactly zero outside. F, F_phi and F_phiphi vanish at endpoints.
    """
    def __init__(self, quiet_mass, knots=8193):
        at = np.linspace(*WINDOW, knots)
        phi = quiet_mass(at)/100 + .5
        if np.min(np.diff(phi)) <= 0:
            raise ValueError('Quiet scalar is not a monotone clock')
        self.phi_knots = phi
        self.s_knots = at
        self.inverse = CubicSpline(phi, at, bc_type=((1, 100/quiet_mass(at[0], 1)),
                                                    (1, 100/quiet_mass(at[-1], 1))),
                                   extrapolate=False)
        self.support = (float(phi[0]), float(phi[-1]))

    def __call__(self, phi):
        phi = np.asarray(phi)
        active = (phi > self.support[0]) & (phi < self.support[1])
        clipped = np.clip(phi, *self.support)
        clock = self.inverse(clipped)
        u1, u2 = self.inverse(clipped, 1), self.inverse(clipped, 2)
        angle = np.pi*(clock-WINDOW[0])/(WINDOW[1]-WINDOW[0])
        k = np.pi/(WINDOW[1]-WINDOW[0])
        sn, cs = np.sin(angle), np.cos(angle)
        gate = sn**4
        g1 = 4*k*sn**3*cs
        g2 = 4*k*k*(3*sn*sn*cs*cs-sn**4)
        return tuple(np.where(active, x, 0.) for x in
                     (gate, g1*u1, g2*u1*u1+g1*u2))


def integrate(native, grid, arm, cutoff, steps=4096, order=64, output=HERE/'tables'):
    """Evolve free, bare-response, static-response and causal cross moments."""
    started = time.monotonic()
    budget = budget_module()
    q, weights, states, background = budget.read_native(native, grid, arm)
    start, loga, mass = background
    quiet = budget.read_native(native, 'fine', 'quiet')[-1][-1]
    clock = ScalarClock(quiet)
    end = float(mass.x[-1])
    # Include source support events and the published readout exactly. Each
    # interval is subdivided to no more than the requested uniform step size.
    from scipy.optimize import brentq
    events = [brentq(lambda s: mass(s)/100+.5-field, start, end, xtol=1e-14)
              for field in clock.support]
    # The archived Hermite background is C1: its second derivatives may jump
    # at native wall knots. Include every such knot so source interpolation
    # never crosses an unresolved background derivative event.
    bounds = sorted(set([float(start), *events, 2.1, end, *mass.x, *loga.x]))
    max_step = (end-start)/steps
    at = np.concatenate([np.linspace(l, r, int(np.ceil((r-l)/max_step))+1)[:-1]
                         for l, r in zip(bounds[:-1], bounds[1:])] + [np.array([end])])
    time_grid = np.empty(2*len(at)-1)
    time_grid[::2] = at
    time_grid[1::2] = (at[:-1]+at[1:])/2
    a = np.exp(loga(time_grid))
    h = loga(time_grid, 1)
    m, dm = mass(time_grid), mass(time_grid, 1)
    phi, velocity = m/100+.5, dm/100
    F, Fphi, Fphiphi = clock(phi)
    fprime = Fphi*velocity
    A = a**-3
    B = a[:, None]*q**2 + (a**3*m**2)[:, None]

    node, quadrature = leggauss(order)
    omega = q[:, None]+(cutoff-q[:, None])*(node+1)/2
    domega = (cutoff-q[:, None])*quadrature/2
    omega2 = omega**2
    g2 = omega*(omega2-q[:, None]**2)*domega/(4*np.pi**2)
    ct = np.sum(g2/omega2, axis=-1)
    z = np.zeros((10, len(q)))
    z[:3] = states['retain']
    mix = np.zeros((4, len(q), order))

    def weighted(v):
        return np.einsum('qj,qj->q', v, g2, optimize=False)

    def rhs(i, state, mixed):
        X, C, P, X2, C2, P2, Xc, Cc, Pc, Dconf = state
        Qx, Qp, Px, Pp = mixed
        alpha, beta, f, ai = A[i], B[i], F[i], a[i]
        potential = ai*f*f*ct
        mx, mp, px = map(weighted, (Qx, Qp, Px))
        dz = np.empty_like(state)
        dz[:3] = (2*alpha*C, alpha*P-beta*X, -2*beta*C)
        dz[3:6] = (2*alpha*C2, alpha*P2-beta*X2-f*mx,
                   -2*(beta*C2+f*px))
        dz[6:9] = (2*alpha*Cc, alpha*Pc-beta*Xc-potential*X,
                   -2*(beta*Cc+potential*C))
        dz[9] = -f*mp
        dt = np.empty_like(mixed)
        dt[0] = alpha*Px+Qp/ai
        dt[1] = alpha*Pp-omega2*Qx/ai-f*X[:, None]
        dt[2] = -beta[:, None]*Qx+Pp/ai-f/(2*omega)
        dt[3] = -beta[:, None]*Qp-omega2*Px/ai-f*C[:, None]
        return dz, dt

    fields = ('rho', 'p', 'j', 'drho', 'dp', 'dj')
    tables = {scheme: {key: np.empty(len(at)) for key in fields}
              for scheme in ('bare', 'static')}
    for table in tables.values():
        for component in COMPONENTS:
            for field in fields:
                table[field+'_'+component] = np.empty(len(at))
        for key in ('rad', 'drad', 'energy_comoving', 'daughter_comoving'):
            table[key] = np.empty(len(at))
    max_ward = {scheme: 0. for scheme in tables}
    max_canonical = 0.

    def save(k):
        nonlocal max_canonical
        i = 2*k
        dz, dt = rhs(i, z, mix)
        X, C, P = z[:3]
        alpha, ai, hi, mi, dmi = A[i], a[i], h[i], m[i], dm[i]
        f, fp, fpp, fdot, phidot = F[i], Fphi[i], Fphiphi[i], fprime[i], velocity[i]
        k2 = (q/ai)**2
        T, dT = weighted(mix[0]), weighted(dt[0])
        D, dD = z[9], dz[9]
        zero = np.zeros_like(q)
        r_d = D/ai**4
        dr_d = (dD-4*hi*D)/ai**4
        r_i = f*T/ai**3
        dr_i = (fdot*T+f*dT)/ai**3-3*hi*r_i
        j_i = fp*T/ai**3
        dj_i = (fpp*phidot*T+fp*dT)/ai**3-3*hi*j_i
        for scheme, use in (('bare', 0), ('static', 1)):
            X2, P2 = z[3]+use*z[6], z[5]+use*z[8]
            dX2, dP2 = dz[3]+use*dz[6], dz[5]+use*dz[8]
            r_p = .5*(P2/ai**6+(mi*mi+k2)*X2)
            p_p = .5*(P2/ai**6-(mi*mi+k2/3)*X2)
            j_p = 100*mi*X2
            dr_p = .5*((dP2-6*hi*P2)/ai**6+(2*mi*dmi-2*hi*k2)*X2+(mi*mi+k2)*dX2)
            dp_p = .5*((dP2-6*hi*P2)/ai**6-(2*mi*dmi-2*hi*k2/3)*X2-(mi*mi+k2/3)*dX2)
            dj_p = 100*(dmi*X2+mi*dX2)
            r_c = use*.5*f*f*ct*X/ai**2
            dr_c = use*ct/ai**2*(f*fdot*X+.5*f*f*dz[0]-hi*f*f*X)
            j_c = use*f*fp*ct*X/ai**2
            dj_c = use*ct/ai**2*(phidot*(fp*fp+f*fpp)*X+f*fp*dz[0]-2*hi*f*fp*X)
            components = {
                'parent': (r_p, p_p, j_p, dr_p, dp_p, dj_p),
                'daughter': (r_d, r_d/3, zero, dr_d, dr_d/3, zero),
                'interaction': (r_i, zero, j_i, dr_i, zero, dj_i),
                'counterterm': (r_c, -r_c/3, j_c, dr_c, -dr_c/3, dj_c),
            }
            table = tables[scheme]
            for number, field in enumerate(fields):
                total = 0.
                for component, values in components.items():
                    value = float(values[number] @ weights)
                    table[field+'_'+component][k] = value
                    total += value
                table[field][k] = total
            table['rad'][k] = table['rho_daughter'][k]
            table['drad'][k] = table['drho_daughter'][k]
            table['energy_comoving'][k] = table['rho'][k]*ai**3
            table['daughter_comoving'][k] = table['rad'][k]*ai**3
            ward = table['drho'][k]+3*hi*(table['rho'][k]+table['p'][k])-table['j'][k]*phidot
            max_ward[scheme] = max(max_ward[scheme], abs(ward))
        canonical = np.max(abs(X*P-C*C-(states['retain'][0]*states['retain'][2]-states['retain'][1]**2)))
        max_canonical = max(max_canonical, float(canonical))

    save(0)
    for k, delta in enumerate(np.diff(at)):
        i = k*2
        s1, t1 = rhs(i, z, mix)
        s2, t2 = rhs(i+1, z+delta/2*s1, mix+delta/2*t1)
        s3, t3 = rhs(i+1, z+delta/2*s2, mix+delta/2*t2)
        s4, t4 = rhs(i+2, z+delta*s3, mix+delta*t3)
        z += delta/6*(s1+2*s2+2*s3+s4)
        mix += delta/6*(t1+2*t2+2*t3+t4)
        save(k+1)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    result = dict(grid=grid, arm=arm, cutoff=cutoff, steps=len(at)-1, order=order,
                  start=float(start), end=end, gate_on=events[0], gate_off=events[1],
                  max_step=float(np.max(np.diff(at))), max_abs_pointwise_ward=max_ward,
                  max_free_covariance_determinant_drift=max_canonical, files={})
    for scheme, table in tables.items():
        path = output/f'{grid}_{arm}_{cutoff}_{scheme}.npz'
        np.savez_compressed(path, s=at, a=a[::2], h=h[::2], phi=phi[::2], v=velocity[::2],
                            F=F[::2], Fphi=Fphi[::2], Fphiphi=Fphiphi[::2], **table)
        readout = int(np.argmin(abs(at-2.1)))
        result['files'][scheme] = dict(path=str(Path(output.name)/path.name), sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            readout={key:float(table[key][readout]) for key in ('rho','p','j','rad','energy_comoving','daughter_comoving')})
    result['seconds'] = time.monotonic()-started
    return result


class SourceTable:
    """Cubic Hermite evaluation with supplied analytic derivatives, zero before fork.

    No extrapolation after the archived endpoint is allowed. Evaluate returns
    (rho,p,j,rad) and derivative returns their proper-time derivatives.
    """
    def __init__(self, path):
        with np.load(path) as table:
            self.s = table['s'].copy()
            values = np.array([table[key] for key in ('rho','p','j','rad')]).T
            slopes = np.array([table[key] for key in ('drho','dp','dj','drad')]).T
        self.spline = CubicHermiteSpline(self.s, values, slopes, extrapolate=False)

    def evaluate(self, s, derivative=0):
        s = np.asarray(s)
        if np.any(s > self.s[-1]+1e-12):
            raise ValueError('Source evaluation beyond archived endpoint')
        value = self.spline(np.clip(s, self.s[0], self.s[-1]), derivative)
        return np.where((s < self.s[0])[..., None], 0., value)


def _job(args):
    return integrate(*args)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native', type=Path, required=True,
                        help='Prior capsule native/ containing FORK.json and WALLS.zip')
    parser.add_argument('--output', type=Path, default=HERE/'tables')
    parser.add_argument('--grid', choices=('fine','coarse','all'), default='all')
    parser.add_argument('--arm', choices=(*ARMS,'all'), default='all')
    parser.add_argument('--cutoff', type=int, choices=(60,120), nargs='+', default=[60,120])
    parser.add_argument('--steps', type=int, default=4096)
    parser.add_argument('--order', type=int, default=64)
    parser.add_argument('--jobs', type=int, default=4)
    options = parser.parse_args()
    grids = ('fine','coarse') if options.grid == 'all' else (options.grid,)
    arms = ARMS if options.arm == 'all' else (options.arm,)
    jobs = [(options.native,g,a,c,options.steps,options.order,options.output)
            for g in grids for a in arms for c in options.cutoff]
    report = dict(schema=1, approximation='fixed-trajectory explicit order-y^2 forcing',
                  units='hatted physical density coefficient per lambda=y^2',
                  derivative='d/ds, s=H0*tau', scalar_clock='fine quiet C2 cubic inverse, 8193 knots',
                  gate_window=list(WINDOW), native_inputs={name:hashlib.sha256((options.native/name).read_bytes()).hexdigest()
                                                          for name in ('FORK.json','WALLS.zip')}, runs={})
    if options.jobs == 1:
        results = map(_job, jobs)
        for result in results:
            key = f"{result['grid']}_{result['arm']}_{result['cutoff']}"
            elapsed = result.pop('seconds')
            report['runs'][key] = result
            print(key, elapsed, flush=True)
    else:
        with ProcessPoolExecutor(max_workers=options.jobs) as pool:
            for result in pool.map(_job, jobs):
                key = f"{result['grid']}_{result['arm']}_{result['cutoff']}"
                elapsed = result.pop('seconds')
                report['runs'][key] = result
                print(key, elapsed, flush=True)
    options.output.mkdir(parents=True, exist_ok=True)
    (options.output.parent/'SOURCES.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
