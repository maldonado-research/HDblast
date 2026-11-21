"""Necessary transverse incoming-characteristic exposure at the wall, with arbitrary Maxwell modes."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
from pathlib import Path
import runpy
import shlex
import subprocess
import tempfile

import numpy as np

BOUND = runpy.run_path(str(Path(__file__).resolve().parents[1]/'global/clumping.py'))
COLS = ('t','a','H','K','E','B','current','incoming','phi','D')


def capture_solver(prior, folder):
    fb = runpy.run_path(str(prior/'feedback/feedback.py'))
    parent = prior/'prior'
    source = runpy.run_path(str(parent/'envelope/build.py'))['transform'](
        (parent/'coupled/joined_forward.cpp').read_text(), centered=True)
    source = runpy.run_path(str(prior/'boost/evolve.py'))['make_source'](source,fb,(prior/'feedback/angular.hpp').read_text())
    change = fb['change']
    source = '#include <fstream>\n' + source
    source = change(source, 'if (argc != 34)', 'if (argc != 36)')
    source = change(source, 'int collision_mode = 0;', 'double branching=1.;\nint collision_mode = 0;')
    source = change(source, 'collision_mode=std::stoi(argv[33]);', '''collision_mode=std::stoi(argv[33]);
    branching=std::stod(argv[34]); std::ofstream trace(argv[35],std::ios::binary);
    trace.exceptions(std::ios::badbit|std::ios::failbit);''')
    source = change(source, 'return proper_density*boosted_rate(', 'return branching*proper_density*boosted_rate(')
    anchor = '        const double margin=gas.intrinsic_random-.9*wnew.e;'
    source = change(source,anchor,'''        {
            const double rec[]={t+dt,a,h,wnew.k,wnew.x,wnew.b,gas.current,y[s.index(0,I)],y[PHI],gas.comoving_rest};
            trace.write(reinterpret_cast<const char*>(rec),sizeof(rec));
        }
'''+anchor)
    file = folder/'capture.cpp'; file.write_text(source)
    exe = folder/'capture'
    subprocess.run(['g++','-std=c++17','-O3','-I',str(prior/'boost'),str(file),'-o',str(exe)],check=True)
    return exe


def integrals(weights, suffix_time, exponent):
    """Vector of sums over cells u of weights[u] * S[max(u,cut)]**exponent."""
    left = np.r_[0.,np.cumsum(weights)[:-1]]
    right = BOUND['suffix'](weights * suffix_time**exponent)
    return left * suffix_time**exponent + right


def lookback(summary, meta, spec, marker, rows):
    end = spec['markers'][marker]['time']
    j = int(np.searchsorted(summary['t'],end))
    d_j = {k:v[:j+1] for k,v in summary.items()}
    t, a, d = (d_j[k] for k in ('t','a','D'))
    dt = np.diff(t)
    K = BOUND['K']
    q = meta['charge_to_mass']
    c = int(np.searchsorted(rows['t'],t[0]))
    ii = np.interp(t,rows['t'],rows['incoming'])
    if np.max(np.abs(np.interp(t,rows['t'],rows['K'])-K))>1.e-10:
        raise ValueError("Gauge coefficient varies after first births")
    maxI = np.maximum(np.abs(ii[:-1]),np.abs(ii[1:]))
    if np.any(np.diff(a)<0.) or d[0]!=0.:
        raise ValueError("Invalid expanding startup")
    scale_d = spec['channel_fraction']*meta['action']*BOUND['TOOLS']['pair_kernel_upper']()/(6.*meta['daughter_mass_L']**3)
    hazard = BOUND['suffix'](scale_d*dt*d[1:]/a[:-1]**3)
    birth = np.diff(d)/d[-1];birth[-1]=0.
    floor = BOUND['suffix'](birth*np.exp(-hazard))
    thermal = BOUND['ThermalCost'](min(d_j['theta'][-2:]),24576)
    weights = dt * np.sqrt(d[1:]/d[-2]) / a[:-1]
    time_kernel = BOUND['suffix'](weights)
    current_field = q / K * np.cumsum(dt*d[1:]/a[:-1])
    self_impulse = math.sqrt(2.)*q/a[-2]*BOUND['suffix'](weights*current_field)
    f = 4.*math.sqrt(2.)*q/(K*a[-2])
    ga = a[1:]**2
    unit = f * integrals(dt*ga,time_kernel,1.)
    l2 = f * np.sqrt(integrals(dt*ga**2,time_kernel,2.))
    prof = f * integrals(dt*ga*maxI,time_kernel,1.)
    pre_t = rows['t'][:c+1]
    pre_i = rows['incoming'][:c+1]
    pre_dt = np.diff(pre_t)
    pre_fluence = float(np.dot(pre_dt,(pre_i[:-1]**2+pre_i[:-1]*pre_i[1:]+pre_i[1:]**2)/3.))
    window_fluence = float(np.dot(dt,(ii[:-1]**2+ii[:-1]*ii[1:]+ii[1:]**2)/3.))
    current = np.interp(t,rows['t'],rows['current'])
    other = 2.*K*a*d_j['B'] + a*a*(2.*d_j['X'] + current)
    incident = 4.*a*a*ii
    lhs = K*(a[-1]**2*d_j['X'][-1] - a[0]**2*d_j['X'][0])
    inferred = lhs + .5*np.dot(dt,other[:-1]+other[1:])
    observed = .5*np.dot(dt,incident[:-1]+incident[1:])
    magnitude = .5*np.dot(dt,abs(incident[:-1])+abs(incident[1:]))
    flux_error = abs(inferred-observed)/magnitude
    if flux_error>5.e-4:
        raise ValueError(f"Direct and reconstructed incoming signal differ: {flux_error}")
    result = dict(time=end, birth_start=t[0], duration_upper=t[-1]-t[0],
        inferred_signed_incoming_canonical_integral=inferred, direct_signed_incoming_canonical_integral=observed,
        independent_birth_window_incoming_relative_error=flux_error,
        incoming_peak_birth_window=float(maxI.max()), incoming_peak_prior_to_birth=float(np.max(abs(pre_i))),
        incoming_rms_birth_window_enlarged=math.sqrt(window_fluence/dt.sum()),
        prebirth_incoming_squared_integral=pre_fluence, startup_incoming_squared_integral=window_fluence,
        bounds={})
    for eps in (0.,.01,.05):
        required=thermal(np.maximum(0.,floor-eps))
        deficit=np.maximum(0.,required-self_impulse)
        peak_ratio=deficit/unit
        rms_ratio=deficit/(l2*math.sqrt(t[-1]-t[0]))
        shape_ratio=deficit/prof
        h=int(np.argmax(shape_ratio))
        output=dict(minimum_peak_spatial_rms_incoming=float(np.max(peak_ratio)),
                minimum_spacetime_rms_incoming=float(np.max(rms_ratio)),
                minimum_incoming_profile_multiple=float(shape_ratio[h]),
                profile_witness_age_L=float(end-t[h]), required_impulse_lower=float(required[h]),
                self_current_impulse_upper=float(self_impulse[h]), unit_profile_impulse_upper=float(prof[h]))
        result['bounds'][str(eps)]=output
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prior',type=Path,required=True);p.add_argument('--history',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--jobs',type=int,default=2)
    opt=p.parse_args()
    prior=opt.prior.resolve()
    meta=json.loads((opt.history/'markers.json').read_text())
    config=json.loads((prior/'prior/envelope/ENVELOPE.json').read_text())['parameters']['mass_4e12']
    run=runpy.run_path(str(prior/'prior/coupled/joined.py'))['run']
    # Independent direct O(N²) control of the optimized integrated response to each source cell.
    support=np.array([6.,3.5,1.5,.4]); w=np.array([.125,.21,.64,1.1])
    for exponent in (1.,2.):
        direct=np.array([sum(w[k]*support[max(i,k)]**exponent for k in range(len(w))) for i in range(len(w))])
        if not np.allclose(integrals(w,support,exponent),direct,rtol=1.e-13,atol=1.e-14):
            raise ValueError("Integrated response-kernel check failed")
    with tempfile.TemporaryDirectory(prefix='hdblast-incoming-') as wd:
        tmp=Path(wd);exe=capture_solver(prior,tmp)
        def case(name):
            spec=meta['runs'][name]
            path=tmp/(name+'.bin');wrapper=tmp/name
            args=['2',repr(meta['action']),str(prior/'boost/BOOST_TABLE.dat'),'1',repr(spec['channel_fraction']),str(path)]
            wrapper.write_text('#!/bin/sh\nexec '+shlex.quote(str(exe))+' "$@" '+shlex.join(args)+'\n');wrapper.chmod(0o755)
            with np.load(opt.history/(name+'.npz')) as tab: data={k:tab[k] for k in tab.files}
            try: run(wrapper,dict(config,cfl=config['cfl']*(.5 if name.endswith('half') else 1.),sample=0,end=data['t'][-1]+.0002))
            except subprocess.CalledProcessError as e: raise RuntimeError(e.stderr) from e
            raw=np.fromfile(path,dtype='=f8').reshape(-1,len(COLS))
            # The initial collar and wall have vanishing EM fields and incoming characteristic.
            initial = np.array([[0.,2.,0.,raw[0,3],0.,0.,0.,0.,0.,0.]])
            raw=np.vstack((initial,raw))
            z={k:raw[:,i] for i,k in enumerate(COLS)}
            guide_at_nodes=np.interp(data['t'],z['t'],z['B'])
            if not np.allclose(data['B'],guide_at_nodes,rtol=2e-6,atol=2e-8):
                raise ValueError("Guide reconstruction differs")
            j=z['a']**2*(4.*z['incoming']-2.*z['E']-z['current'])-2.*z['K']*z['a']*z['B']
            dt=np.diff(z['t'])
            canonical=z['a']**2*z['K']*z['E']
            diff=canonical[1:]-np.cumsum(.5*dt*(j[1:]+j[:-1]))
            control={'integrated_maxwell_residual_over_peak_canonical':float(np.max(abs(diff))/np.max(abs(canonical))),
                'max_abs_guide_difference':float(np.max(abs(data['B']-guide_at_nodes)))}
            if control['integrated_maxwell_residual_over_peak_canonical']>1e-4:
                raise ValueError(control)
            return name,dict(control=control,markers={k:lookback(data,meta,spec,k,z) for k in spec['markers']})
        with ThreadPoolExecutor(opt.jobs) as pool: rows=dict(pool.map(case,meta['runs']))
    opt.output.write_text(json.dumps(rows,indent=2,sort_keys=True,allow_nan=False)+'\n')


if __name__=='__main__':main()
