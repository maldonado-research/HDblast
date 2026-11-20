"""Linear transverse Maxwell transfer on the specified axion wall history."""
import argparse
from bisect import bisect_right
from decimal import Decimal, localcontext
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
HARMONICS = (2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 12288, 16384, 24576, 32768, 49152, 65536)


def smooth(x):
    if x <= 0:
        return 0.
    if x >= 1:
        return 1.
    z = 1/(1-x)-1/x
    exp = math.exp(-abs(z))
    return 1/(1+exp) if z >= 0 else exp/(1+exp)


def trajectory(path, scalar_start, width, smoothed_stride=1):
    record = json.loads(path.read_text())
    rows = [dict(zip(record['columns'], row, strict=True)) for row in record['rows']]
    data = [tuple(row[k] for k in ('t','a','psi','pi','K')) for row in rows]
    brackets = [i for i, row in enumerate(data) if 1e-10 < (row[2]-scalar_start)/width < 1.-1e-10]
    if not brackets or min(brackets) == 0 or max(brackets) == len(data)-1:
        raise ValueError('Provide samples spanning the complete scalar step and both free endpoints.')
    data = data[min(brackets)-1:max(brackets)+2]
    if any(r1[0] <= r0[0] or r1[2] <= r0[2] for r0, r1 in zip(data, data[1:])):
        raise ValueError('Requires a single monotone passage through the selected scalar step.')
    km, kp = data[0][4], data[-1][4]
    k_error = max(abs((km+(kp-km)*smooth((psi-scalar_start)/width))-meas)/meas
                  for _, _, psi, _, meas in data)
    if k_error > 0.001:
        raise ValueError(f'Trace and constitutive law disagree (relative {k_error}).')
    data = data[::smoothed_stride] + ([data[-1]] if (len(data)-1)%smoothed_stride else [])
    times = [v[0] for v in data]

    def state(t):
        j = min(len(data)-2,max(0,bisect_right(times,t)-1))
        t0,a0,psi0,pi0,_ = data[j]
        t1,a1,psi1,pi1,_ = data[j+1]
        dt = t1-t0
        x = (t-t0)/dt
        # Cubic Hermite interpolation anchored by the proper scalar velocity.
        psi = psi0 + x*(dt*pi0+x*(3*(psi1-psi0)-dt*(2*pi0+pi1)
                                      +x*(dt*(pi0+pi1)-2*(psi1-psi0))))
        return a0+x*(a1-a0), smooth((psi-scalar_start)/width)
    return state, dict(start=times[0], stop=times[-1], K_in=km, K_out=kp,
                       samples=len(data), trace_K_law_max_relative_error=k_error)


def layers(state, meta, theta_jump, n):
    t0 = meta['start']
    dt = (meta['stop']-t0)/n
    k0, dk = meta['K_in'], meta['K_out']-meta['K_in']
    prev = state(t0)[1]
    result = []
    for i in range(n):
        a, mid = state(t0+(i+.5)*dt)
        current = state(t0+(i+1)*dt)[1]
        result.append((dt/a,k0+dk*mid,theta_jump*(current-prev)))
        prev = current
    return result


def transfer(grid, k0, kf, k, helicity, precise=False, seed=None):
    """Cayley steps on (sqrt(k K_in) A, Pi/sqrt(k K_in)), Pi=K A'."""
    cast = Decimal.from_float if precise else float
    a,b,c,d = cast(1.),cast(0.),cast(0.),cast(1.)
    ref = k0
    for deta, K, dtheta in grid:
        u = cast(k*deta*ref/K)
        v = cast((k*deta*K + helicity*dtheta)/ref)
        z = u*v/4
        inv = 1/(1+z)
        diag = (1-z)*inv
        up = u*inv
        dn = -v*inv
        a,b,c,d = diag*a+up*c,diag*b+up*d,dn*a+diag*c,dn*b+diag*d
    ratio = cast(kf/k0)
    total = float(ratio*(a*a+b*b)+(c*c+d*d)/ratio)
    det = float(a*d-b*c-1)
    if total < 2-1e-7:
        raise ValueError("Unphysical energy gain for a real symplectic transfer.")
    total = max(2., total)
    result = dict(max_energy_gain=(total+math.sqrt((total-2)*(total+2)))/2,
                  phase_averaged_energy_gain=total/2, determinant_minus_one=det)
    if seed is not None:
        b0, e0 = map(cast,seed)
        magnetic = ratio*(a*b0+b*e0)**2
        electric = (c*b0+d*e0)**2/ratio
        result['seed_gain'] = float((magnetic+electric)/(b0*b0+e0*e0))
        result['seed_final_electric_share'] = float(electric/(magnetic+electric))
    return result


def gain_errors(reference, test):
    return max(abs(result['max_energy_gain']/test[k][h]['max_energy_gain']-1)
               for k, helicities in reference.items() for h, result in helicities.items())


def independent_checks():
    # Frozen oscillator. The Cayley propagator is a rotation for arbitrary step size.
    grid = [(.013,7.,0.)]*19
    static = transfer(grid,7.,7.,6,1,True)
    if abs(static['max_energy_gain']-1) > 2e-7:
        raise ValueError('Static energy-normalized wave must preserve energy.')
    # Impulsive limit at fixed A, Pi_out=Pi_in-h k DeltaTheta A for either helicity.
    ratio, delta, kin = 35.,9.,.4
    instant = transfer([(0.,kin,delta)],kin,ratio*kin,4,-1,True)
    expected = (ratio + 1/ratio + (delta/kin)**2/ratio)/2
    if abs(instant['phase_averaged_energy_gain']/expected-1)>1e-11:
        raise ValueError('Sudden-jump symplectic shear disagrees with analytic transfer.')
    return dict(frozen_gain=static['max_energy_gain'], sudden_jump_avg_gain=expected)


def calculate(path, scalar_start, width, theta_jump, n):
    with localcontext() as context:
        context.prec = 70
        checks = independent_checks()
        def run(steps, skip, precise=False):
            state, meta = trajectory(path, scalar_start, width, skip)
            grid = layers(state, meta, theta_jump, steps)
            out = {str(k): {str(h):transfer(grid, meta['K_in'], meta['K_out'], k,h,precise)
                            for h in (-1,0,1)} for k in HARMONICS}
            return out, meta, grid
        full, meta, grid = run(n,1,True)
        coarse = run(n//2,1)[0]
        sparser = run(n,2)[0]
        fields = json.loads(path.read_text())
        initial,final = [dict(zip(fields['columns'],r,strict=True)) for r in
                         (fields['rows'][0],fields['rows'][-1])]
        if initial['t']!=meta['start'] or final['t']!=meta['stop']:
            raise ValueError('The k=2 benchmark requires the same initial and final field epochs.')
        seed = transfer(grid,meta['K_in'],meta['K_out'],2,1,True,(initial['B'],initial['X']))
        measured = ((final['a']/initial['a'])**4 * final['K']/initial['K']
                    * (final['B']**2+final['X']**2)/(initial['B']**2+initial['X']**2))
        homogeneous = dict(spectator_energy_gain=seed['seed_gain'],actual_wall_energy_gain=measured,
                           spectator_final_electric_share=seed['seed_final_electric_share'],
                           actual_wall_final_electric_share=final['X']**2/(final['B']**2+final['X']**2))
        diagnostics = dict(
            cayley_half_steps_max_relative_gain_change=gain_errors(full,coarse),
            half_source_samples_max_relative_gain_change=gain_errors(full,sparser),
            max_abs_wronskian_error=max(abs(r['determinant_minus_one']) for hs in full.values() for r in hs.values()),
            max_sampled_abs_theta_prime_over_K=max(abs(dt)/(eta*K) for eta,K,dt in grid),
            conformal_duration=math.fsum(row[0] for row in grid), steps=n, **checks)
        if diagnostics['max_abs_wronskian_error']>1e-30:
            raise ValueError('High-precision Wronskian failed.')
    return dict(background=meta, theta_jump=theta_jump, scalar_start=scalar_start, width=width,
                numerical_checks=diagnostics, homogeneous_seed_comparison=homogeneous, modes=full)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--trace', type=Path, default=HERE.parent/'axion/HANDOFF.json')
    p.add_argument('--scalar-start', type=float, default=.369)
    p.add_argument('--width', type=float, default=.0006)
    p.add_argument('--theta', type=float, default=28000.)
    p.add_argument('--steps', type=int, default=16384)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    result = calculate(args.trace,args.scalar_start,args.width,args.theta,args.steps)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
