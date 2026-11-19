"""Short-time source test using a fundamental matrix for each spectral dyad."""
import argparse
import io
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import brentq

ON = 1.35
WIDTH = 0.70
R = 4


def native_start(native, grid='fine', activation_mass=None):
    fork = json.loads((native / 'FORK.json').read_bytes())[grid]
    with ZipFile(native / 'WALLS.zip') as zip_file:
        raw = zip_file.read(f'{grid}/retain/quiet/wall.csv')
    wall = np.genfromtxt(io.BytesIO(raw), delimiter=',', names=True)
    s = wall['H0tau']
    loga = CubicHermiteSpline(s, wall['A_b']-fork['reference_A'], wall['H_over_H0'])
    mass = CubicHermiteSpline(s, 100*(wall['phi_b']-.5), 100*wall['v_over_H0'])
    q, weights = np.asarray(fork['q']), np.asarray(fork['weight'])
    initial = np.array([fork['states']['retain'][k] for k in ['X', 'C', 'P']])
    def rhs(t, vector):
        x, c, p = vector.reshape(3,-1)
        a = np.exp(loga(t))
        invm, stiffness = a**-3, a*q*q+a**3*mass(t)**2
        return np.r_[2*invm*c, invm*p-stiffness*x, -2*stiffness*c]
    activation = ON if activation_mass is None else brentq(
        lambda t: mass(t)-activation_mass, s[0], s[-1], xtol=1e-14)
    result = solve_ivp(rhs, (float(s[0]), activation), initial.ravel(), method='DOP853',
                       rtol=2e-12, atol=1e-14, max_step=.002)
    if not result.success:
        raise RuntimeError(result.message)
    return q, weights, result.y[:,-1].reshape(3,-1), loga, mass, activation


def matrices(omega, q, a, mass, f_scaled):
    # Unit coupling to one oscillator. The quadrature weights g_j² enter the
    # trace after evolving its first and second response coefficients.
    m = np.zeros((len(q),4,4))
    m[:,0,1] = a**-3
    m[:,1,0] = -(a*q*q+a**3*mass**2)
    m[:,2,3] = 1/a
    m[:,3,2] = -omega**2/a
    kick = np.zeros_like(m)
    kick[:,1,2] = kick[:,3,0] = -f_scaled
    ct = np.zeros_like(m)
    ct[:,1,0] = -a*f_scaled**2 / omega**2
    return m, kick, ct


def dyad_run(native_data, cutoff, interval, order, steps):
    q, qw, moments, loga, mass, _ = native_data
    z, wg = np.polynomial.legendre.leggauss(order)
    om = q[:,None]+(cutoff-q[:,None])*(1+z)/2
    dw = (cutoff-q[:,None])*wg/2
    strength = (qw[:,None]*dw*om*(om**2-q[:,None]**2)/(4*np.pi**2)).ravel()
    q = np.repeat(q,order)
    om = om.ravel()
    x0,c0,p0 = np.repeat(moments,order,axis=1)
    v = np.zeros((len(q),4,4))
    v[:,0,0],v[:,0,1],v[:,1,0],v[:,1,1] = x0,c0,c0,p0
    v[:,2,2],v[:,3,3] = .5/om,.5*om
    # Each row evolves a 4x4 fundamental matrix, initially (I,0,0,0).
    # The rows are the free solution and the linear, bare quadratic and
    # explicit counterterm responses divided by u^(r+1),u^(2r+2),u^(2r+1).
    state = np.zeros((4,len(q),4,4))
    state[0] = np.eye(4)
    def rhs(t, state):
        s = ON+interval*t
        a = np.exp(loga(s))
        f = (np.sin(np.pi*interval*t/WIDTH)/interval)**R
        m, kick, ct = matrices(om,q,a,mass(s),f)
        dm = (interval*m) @ state
        dm[1] += kick @ state[0]
        dm[2] += kick @ state[1]
        dm[3] += ct @ state[0]
        return dm
    h=1/steps
    for i in range(steps):
        t=i*h
        k1=rhs(t,state)
        k2=rhs(t+h/2,state+h*k1/2)
        k3=rhs(t+h/2,state+h*k2/2)
        k4=rhs(t+h,state+h*k3)
        state += h*(k1+2*k2+2*k3+k4)/6
    free, one, two, ct = state
    trans = lambda z:z.swapaxes(-1,-2)
    second = (two @ v) @ trans(free)
    second = second+trans(second)+(one @ v) @ trans(one)
    mixed = (one @ v) @ trans(free)
    mixed = (mixed+trans(mixed))[:,0,2]
    seagull = (ct @ v) @ trans(free)
    seagull += trans(seagull).copy()
    free_x = ((free @ v) @ trans(free))[:,0,0]
    s=ON+interval
    a=np.exp(loga(s))
    m=mass(s)
    mprime=mass(s,1)
    f=(np.sin(np.pi*interval/WIDTH)/interval)**R
    df=R*np.pi/WIDTH*(np.sin(np.pi*interval/WIDTH)/interval)**(R-1)*np.cos(np.pi*interval/WIDTH)
    def parent_energy(v):
        return .5*(v[:,1,1]/a**3+(a*q**2+a**3*m**2)*v[:,0,0])
    def parent_pressure_comoving(v):
        return .5*(v[:,1,1]/a**3-(a*q**2/3+a**3*m**2)*v[:,0,0])
    parent = parent_energy(second)
    daughter=.5*(second[:,3,3]+om**2*second[:,2,2])/a
    interaction=f*mixed/interval
    bare=parent+daughter+interaction
    ct_energy=.5*a*f*f*free_x/om**2
    stat=interval**2*bare+interval*parent_energy(seagull)+ct_energy
    p_bare=parent_pressure_comoving(second)+daughter/3
    p_static=interval**2*p_bare+interval*parent_pressure_comoving(seagull)-ct_energy/3
    j_bare=interval*a**3*m*mprime*second[:,0,0]+df*mixed/interval
    j_static=interval**2*(j_bare+a**3*m*mprime*seagull[:,0,0])+a*f*df*free_x/om**2
    observed = dict(parent=parent,daughter=daughter,bare_total=bare,
                    bare_pressure=p_bare,bare_scalar_power=j_bare,
                    static_total=stat,static_pressure=p_static,static_scalar_power=j_static)
    a0=np.exp(loga(ON))
    alpha=(np.pi/WIDTH)**R
    gp=1/(2*om*a0**3)
    gd=x0/a0
    up=alpha**2*gp/(2*(R+1)**2)
    ud=alpha**2*gd/(2*(R+1)**2)
    leading_ct=.5*a0*alpha**2*x0/om**2
    expected=dict(parent=up,daughter=ud,bare_total=-R*(up+ud)/(R+2),
                  bare_pressure=up+ud/3,bare_scalar_power=-(R/(R+2))*(2*R+2)*(up+ud),
                  static_total=leading_ct,static_pressure=-leading_ct/3,
                  static_scalar_power=2*R*leading_ct)
    weighted=lambda v:float(np.dot(strength,v))
    ratios={k:weighted(observed[k])/weighted(val) for k,val in expected.items()}
    return {'interval':interval, 'bare_total_over_free':weighted(bare)/weighted(parent+daughter),
            'bare_interaction_over_free':weighted(interaction)/weighted(parent+daughter),
            'static_pressure_over_density':weighted(p_static)/weighted(stat),
            'ratio_to_leading':ratios,
            'leading_comoving_coefficients':{k:weighted(x) for k,x in expected.items()}}


def band_threshold(native_data):
    q, w, moments, loga, _, activation = native_data
    a = float(np.exp(loga(activation)))
    x = moments[0]
    # Continuous frequency integrals with the same finite 24-node parent band.
    Polynomial = np.polynomial.Polynomial
    parent = Polynomial([w @ (2*q**3/3), -(w @ q**2), 0., sum(w)/3]) / (8*np.pi**2*a**3)
    daughter = Polynomial([w @ (x*q**4), 0., -2*(w @ (x*q*q)), 0., w @ x]) / (16*np.pi**2*a)
    difference = parent-daughter
    shifted = difference(Polynomial([8,1])).coef
    if not (min(shifted[:3]) > 0 > max(shifted[3:])):
        raise ValueError("The polynomial does not certify one positive root in the shifted interval")
    zero = brentq(difference,8.,60.,xtol=1e-13)
    # Independent finite Gauss check on the unexpanded spectral integrals.
    z, wg = np.polynomial.legendre.leggauss(32)
    om = q[:,None]+(zero-q[:,None])*(1+z)/2
    dw = (zero-q[:,None])*wg/2
    g2 = w[:,None]*dw*om*(om*om-q[:,None]**2)/(4*np.pi**2)
    up = np.sum(g2/(2*om*a**3))
    ud = np.sum(g2*x[:,None]/a)
    rel = float((up-ud)/(up+ud))
    if abs(rel) > 1e-11:
        raise ValueError("Independent spectral integration at the root")
    return {"activation_time": float(activation), "threshold": float(zero),
            "shifted_polynomial_ascending": shifted.tolist(), "relative_integral_difference_at_root": rel,
            "enthalpy_over_free_gain": {str(c):float(difference(c)/(3*(parent(c)+daughter(c))))
                                         for c in (8,12,18,24,40,60,120)}}


def compare_tables(native_data, tables):
    # Compare scaled Taylor coefficients of two independent integrations on
    # the same output rows. This avoids subtracting extensive vacuum energies.
    by_scheme = {}
    for scheme in ("bare", "static"):
        with np.load(tables / f"fine_quiet_60_{scheme}.npz") as table:
            by_scheme[scheme] = {k: table[k] for k in
                                 ("s", "a", "v", "rho_parent", "rho_daughter", "rho", "p", "j")}
    bare, static = by_scheme["bare"], by_scheme["static"]
    if not np.array_equal(bare["s"], static["s"]):
        raise ValueError("The two tables must share a clock")
    origin = int(np.argmin(abs(bare["s"] - ON)))
    if abs(bare["s"][origin] - ON) > 1e-12:
        raise ValueError("Expected the gate onset to be recorded")
    report = {}
    for cells in (2, 4, 8, 16, 32):
        i = origin + cells
        delay = bare["s"][i] - ON
        volume = bare["a"][i]**3
        velocity = bare["v"][i]
        direct = dyad_run(native_data, 60., delay, 32, 72)
        values = {
            "parent": bare["rho_parent"][i]*volume/delay**10,
            "daughter": bare["rho_daughter"][i]*volume/delay**10,
            "bare_total": bare["rho"][i]*volume/delay**10,
            "bare_pressure": bare["p"][i]*volume/delay**10,
            "bare_scalar_power": bare["j"][i]*volume*velocity/delay**9,
            "static_total": static["rho"][i]*volume/delay**8,
            "static_pressure": static["p"][i]*volume/delay**8,
            "static_scalar_power": static["j"][i]*volume*velocity/delay**7,
        }
        expected = direct["leading_comoving_coefficients"]
        ratios = direct["ratio_to_leading"]
        errors = {k: float(val/(ratios[k]*expected[k]) - 1) for k, val in values.items()}
        report[str(cells)] = {"delay": float(delay), "relative_difference": errors,
                              "max_relative_difference": max(map(abs, errors.values()))}
    if report["32"]["max_relative_difference"] > 2e-6:
        raise ValueError("Dyad integration disagrees with injected matter")
    return report


def experiment(native):
    initial=native_start(native)
    thresholds = {"fine": band_threshold(initial),
                  "coarse": band_threshold(native_start(native, 'coarse', float(initial[4](ON))))}
    rows=[dyad_run(initial,60.,u,16,48) for u in (.01,.003,.001,.0003,.0001,.00003)]
    refined=dyad_run(initial,60.,.00003,24,96)
    tol=.003
    if (max(abs(x-1) for x in rows[-1]['ratio_to_leading'].values())>tol or
        abs(rows[-1]['bare_total_over_free']+2/3)>tol or
        abs(rows[-1]['static_pressure_over_density']+1/3)>tol):
        raise ValueError('Short-time limit or sign failed')
    if max(abs(refined['ratio_to_leading'][k]-x) for k,x in rows[-1]['ratio_to_leading'].items())>2e-5:
        raise ValueError('Refinement failed')
    return {'wall':'fine/retain/quiet','native_parent_band':[0,8],'cutoff':60,
            'switch_power':R,'clock_on':ON,'base_steps':48,'base_frequency_nodes':16,'enthalpy_threshold':thresholds,
            'refined_steps':96,'refined_frequency_nodes':24,'delays':rows,'refined':refined}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('native',type=Path)
    parser.add_argument("--tables", type=Path, help="Compare against extracted published source tables")
    args = parser.parse_args()
    result = (compare_tables(native_start(args.native), args.tables)
              if args.tables else experiment(args.native))
    print(json.dumps(result, indent=2, allow_nan=False))
