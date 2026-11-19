"""Independent spectral energy, Ward, component, and table interpolation checks."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import cumulative_simpson, solve_ivp
from scipy.interpolate import CubicSpline
from scipy.special import roots_legendre

sys.path.insert(0, str(Path(__file__).resolve().parent))
from source import HERE, ScalarClock, SourceTable, budget_module


def free_parent(q, state, background):
    """Independent complex free modes, as in the preceding clock capsule."""
    start, loga, mass = background
    a0 = np.exp(loga(start))
    inertia = a0**3*np.hypot(q/a0,mass(start))
    initial = 1/np.sqrt(2*inertia)
    X,C,P = state
    n = (inertia*X+P/inertia-1)/2
    c = (inertia*X-P/inertia)/2+1j*C
    size = len(q)
    def rhs(t,y):
        a = np.exp(loga(t))
        return np.r_[y[size:2*size]/a**3,-a**3*((q/a)**2+mass(t)**2)*y[:size],1/a]
    sol = solve_ivp(rhs,(start,2.1),np.r_[initial,-1j*inertia*initial,0.],
                    method='DOP853',dense_output=True,rtol=1e-11,atol=5e-14,max_step=.004)
    if not sol.success:
        raise RuntimeError(sol.message)
    u,p = sol.y[:size],sol.y[size:2*size]
    error = float(np.max(abs(u*p.conj()-u.conj()*p-1j)))
    assert error < 3e-9
    return sol.sol,n,c,error


def pair_energy(q,weights,n,c,modes,eta,scale,gate_weight,cutoff,order):
    """Independent inclusive spectral-amplitude daughter energy."""
    nodes,wg = roots_legendre(order)
    total = 0.
    for k,wq,pop,anomalous,wave in zip(q,weights,n,c,modes,strict=True):
        omega = k+(cutoff-k)*(nodes+1)/2
        dweight = wg*(cutoff-k)/2*(omega**2-k**2)/(8*np.pi**2)
        phase = np.exp(1j*omega[:,None]*eta)
        amp_a = np.einsum('ij,j->i',phase,wave*gate_weight)
        amp_b = np.einsum('ij,j->i',phase.conj(),wave*gate_weight)
        probability = pop*abs(amp_a)**2+(pop+1)*abs(amp_b)**2+2*(anomalous*amp_a*amp_b).real
        assert np.min(probability) > -1e-13
        total += wq*np.dot(dweight*omega/scale,probability)
    return float(total)


def verify(native, tables=HERE/'tables', output=HERE/'VERIFICATION.json'):
    budget = budget_module()
    report = {'tables': {}, 'independent_spectral_energy': {}, 'quiet_calibration': {}}
    mquiet = budget.read_native(native, 'fine', 'quiet')[-1][-1]
    clock = ScalarClock(mquiet)
    test = np.linspace(1.35, 2.05, 200001)
    report['quiet_calibration']['max_abs_gate_error'] = float(np.max(abs(clock(mquiet(test)/100+.5)[0]-np.sin(np.pi*(test-1.35)/.7)**4)))
    assert report['quiet_calibration']['max_abs_gate_error'] < 2e-12
    paths = sorted(Path(tables).glob('*.npz'))
    assert len(paths) == 32
    for path in paths:
        grid, arm, cutoff, scheme = path.stem.split('_')
        _, _, _, (_, loga, mass) = budget.read_native(native, grid, arm)
        with np.load(path) as z:
            s = z['s']
            h, v = z['h'], z['v']
            ward = z['drho']+3*h*(z['rho']+z['p'])-z['j']*v
            component_errors = {field:float(np.max(abs(z[field]-sum(z[field+'_'+component] for component in ('parent','daughter','interaction','counterterm')))))
                                for field in ('rho','p','j','drho','dp','dj')}
            pre = z['F'] == 0
            pre &= s < 1.35
            assert all(np.all(z[field][pre] == 0) for field in ('rho','p','j','drho','dp','dj','rad','drad'))
            work = cumulative_simpson(z['a']**3*(z['j']*v-3*h*z['p']), x=s, initial=0)
            work_error = float(np.max(abs(z['energy_comoving']-work)))
            slopes = np.array([z[k] for k in ('drho','dp','dj')]).T
            value_spline = CubicSpline(s, np.array([z[k] for k in ('rho','p','j')]).T)
            slope_difference = value_spline(s[2:-2],1)-slopes[2:-2]
            slope_scale = np.maximum(np.max(abs(slopes),axis=0), 1.)
            slope_error = np.max(abs(slope_difference),axis=0)/slope_scale
        table = SourceTable(path)
        middle = (s[:-1]+s[1:])/2
        val, deriv = table.evaluate(middle), table.evaluate(middle,1)
        hh, vv = loga(middle,1), mass(middle,1)/100
        interpolation_ward = deriv[:,0]+3*hh*(val[:,0]+val[:,1])-val[:,2]*vv
        report['tables'][path.stem] = dict(max_abs_pointwise_ward=float(np.max(abs(ward))),
            max_abs_midpoint_interpolated_ward=float(np.max(abs(interpolation_ward))),
            max_abs_integrated_comoving_work_error=work_error,
            max_component_sum_error=component_errors,
            independently_reconstructed_derivative_relative_max=dict(zip(('rho','p','j'),map(float,slope_error))))
        assert np.max(abs(ward)) < 1e-8
        assert max(component_errors.values()) == 0
        assert np.max(abs(interpolation_ward)) < 1e-5
        assert work_error < 1e-6

    # This is the prior independent oscillator-mode spectral-amplitude
    # calculation, rather than the causal covariance ODE used for the tables.
    q, weights, states, background = budget.read_native(native, 'fine', 'quiet')
    free, n, c, commutator = free_parent(q, states['retain'], background)
    at = np.linspace(1.35, 2.05, 2049)
    y = free(at)
    w = np.ones(len(at)); w[1:-1:2] = 4; w[2:-1:2] = 2
    w *= (at[-1]-at[0])/(len(at)-1)/3
    gate_weight = np.sin(np.pi*(at-1.35)/.7)**4*w
    for cutoff in (60,120):
        energy = pair_energy(q,weights,n,c,y[:len(q)],y[-1].real,
                    np.exp(background[1](2.1)),gate_weight,cutoff,96)
        with np.load(Path(tables)/f'fine_quiet_{cutoff}_bare.npz') as z:
            index = np.argmin(abs(z['s']-2.1))
            actual = float(z['daughter_comoving'][index])
        report['independent_spectral_energy'][str(cutoff)] = dict(spectral=energy,causal_table=actual,
            difference=actual-energy,free_canonical_commutator_error=commutator)
        np.testing.assert_allclose(actual,energy,rtol=3e-9,atol=1e-9)
    Path(output).write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({key:value for key,value in report.items() if key!='tables'},indent=2))
    print('All32 source tables passed component, exact-zero, Ward and interpolation checks.')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native',type=Path,required=True)
    parser.add_argument('--tables',type=Path,default=HERE/'tables')
    parser.add_argument('--output',type=Path,default=HERE/'VERIFICATION.json')
    args = parser.parse_args()
    verify(args.native,args.tables,args.output)
