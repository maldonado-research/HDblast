"""Transfer of a prescribed tangential incident characteristic to an initially empty wall mode."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
from pathlib import Path
import runpy
import subprocess
import tempfile

import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
CURL = np.array([2., 5., 10., 20., 50.])
PHYSICAL = np.array([.5, 1., 2., 5., 10., 25., 50., 100.])
DELAY = (.10, .50)
KINDS = ('whole', 'before_birth', 'after_birth', 'delayed_0.10', 'delayed_0.50')


def run_solver(exe, folder, trajectory, modes, birth, keep, tag):
    # The raw input and C++ output stay native-endian on the machine running the calculation.
    inputs = folder / (tag + '.input')
    spec = folder / (tag + '.spec')
    nodes = folder / (tag + '.nodes')
    out = folder / (tag + '.out')
    np.asarray(trajectory, dtype='=f8').tofile(inputs)
    np.savetxt(spec, modes, fmt=('%.17g', '%d'))
    np.savetxt(nodes, keep, fmt='%d')
    subprocess.run([str(exe), str(inputs), str(spec), repr(float(birth)), str(out), str(nodes)], check=True)
    actual = np.fromfile(out, dtype='=f8').reshape(len(keep), len(modes), 7)
    if not np.isfinite(actual).all():
        raise ValueError('Nonfinite response')
    return actual


def checks(exe, folder):
    # Closed form with constant radius, K and drive. The augmented matrix exponential
    # supplies a phase-sensitive solution without any time integration.
    a, K, k, inc, T = 2.7, .83, 7., .12, 1.7
    matrix = np.array([[-2./K, -k*K/a, 4.*a*a*inc], [k/(K*a), 0., 0.], [0., 0., 0.]])
    ref = (expm(T*matrix) @ np.array([0., 0., 1.]))[:2] / np.array([K*a*a, a*a])
    errors = []
    for n in (40, 80, 160):
        t = np.linspace(0., T, n+1)
        trajectory = np.column_stack((t,np.full_like(t,a),np.full_like(t,K),np.full_like(t,inc),
                                      np.zeros_like(t), np.zeros_like(t)))
        out = run_solver(exe, folder, trajectory, [(k,0)], T, [n], 'check')
        errors.append(float(np.max(abs(out[0,0,:2]-ref))))
    if min(np.array(errors[:-1])/errors[1:])<3.9:
        raise ValueError('Midpoint refinement failed')
    return dict(constant_medium_exact_EB=ref.tolist(), constant_medium_errors=errors)


def fields(values, time_index, modes, scale, limit):
    # Mode energy is the physical spatial mean multiplied by a^4. Only compare
    # time-local field norms with a physical threshold in the same normalization.
    result = {}
    for j, (k, drive) in enumerate(modes):
        rec = values[time_index, j]
        e, b, energy, work, leak, scalar, peak = map(float, rec)
        rms = math.hypot(e,b)/math.sqrt(2.)
        name = (format(k,'.0f') if k in CURL else 'physical_'+format(k/scale,'.6g'))
        event = dict(k_physical=k/scale, E_cos=e, B_sin=b, spatial_rms=rms,
                     canonical_energy=energy, incident_work=work, leakage=leak, scalar_work=scalar,
                     peak_canonical_energy=peak)
        if limit is not None and rms>0.:
            event['fraction_for_initial_rms_floor'] = limit/rms
        result.setdefault(name,{})[KINDS[drive]]=event
    return result


def solve_case(prior, rawdir, tmp, outdir, exe, label, spec, floor, dense):
    data = np.memmap(rawdir / (label + '.bin'), dtype='=f8', mode='r').reshape(-1,10)
    initial = np.array([[0., 2., 0., data[0,3], 0., 0., 0., 0., 0., 0.]])
    raw = np.vstack((initial, data))
    t, a, K, inc = (raw[:,j] for j in (0,1,3,7))
    t_birth = spec['first_birth_left_knot']
    hist = prior / 'history' / (label + '.npz')
    with np.load(hist) as z:
        window = z['t']
    after = np.searchsorted(t,window)
    if after[0] < 1 or after[-1] >= len(t):
        raise ValueError('Source does not cover the particle history')
    before = after-1
    ids = np.where(abs(t[before]-window)<abs(t[after]-window),before,after)
    if (np.max(abs(t[ids]-window)) > 1.e-8 or not np.all(np.diff(ids)==1)
            or abs(t[ids[0]]-t_birth)>1.e-8 or raw[ids[0],9]!=0.):
        raise ValueError('Reference birth knots differ')
    # Floating-point contraction can shift a proper-time stamp by a few ulps.
    # Use the captured knot for actual time integration and verify it is the same
    # last empty population state as the input particle history.
    birth_index = ids[0]
    t_birth = float(t[birth_index])
    a_birth = float(a[birth_index])
    delays = [np.interp(t-delay,t,inc,left=0.) for delay in DELAY]
    traj = np.column_stack((t,a,K,inc,*delays))
    ks = np.r_[CURL, a_birth*PHYSICAL]
    modes = [(k,kind) for k in ks for kind in range(len(KINDS))]
    # Monitoring all prebirth fields at these nodes catches sign, normalization
    # and variable-K errors, including the early sharp maximum of the incident wave.
    monitor = np.unique(np.r_[np.arange(0,birth_index,64),np.argmax(abs(inc[:birth_index]))])
    keep = np.unique(np.r_[monitor, ids])
    at_birth = np.searchsorted(keep,birth_index)
    i_trace = np.searchsorted(keep,ids)
    result = run_solver(exe,tmp,traj,modes,t_birth,keep,label)
    # Prebirth-only and postbirth-only fields must add back to the solution with both drives.
    powers = result.reshape(len(keep),len(ks),len(KINDS),7)
    whole = powers[:,:,0]
    split_error = float(np.max(abs(whole[:,:,:2]-powers[:,:,1,:2]-powers[:,:,2,:2])))
    normalized_balance = abs(result[:,:,2]-result[:,:,3]+result[:,:,4]-result[:,:,5])/np.maximum(result[:,:,6],1.e-10)
    balance_error = float(np.max(normalized_balance))
    control_nodes = keep[:at_birth+1]
    axial = raw[control_nodes,4:6]
    error = np.max(abs(result[:at_birth+1,0,:2]-axial),axis=0)
    relative = error/np.max(abs(axial),axis=0)
    if max(relative)>1.e-4 or split_error>1.e-7 or balance_error>1.e-7:
        raise ValueError((label,'canonical control',relative,split_error,balance_error))
    trace_dir = outdir/'traces'
    trace_dir.mkdir(exist_ok=True)
    trace = dict(t=t[ids], k=ks/a_birth, k_comoving=ks, incoming=inc[ids],
                 E_cos=whole[i_trace,:,0].T, B_sin=whole[i_trace,:,1].T,
                 E_pre_cos=powers[i_trace,:,1,0].T, B_pre_sin=powers[i_trace,:,1,1].T)
    np.savez_compressed(trace_dir/(label+'.npz'),**trace)
    # The 0.1% benchmark is the more demanding initial RMS floor, under uniform density,
    # no new incoming drive after birth and arbitrary allowed self currents.
    birth_floor = floor[label]['0.001']['general_fields_eps_0.01']['initial_rms_field_lower_by_density_cap']['1']
    summaries = {'birth':fields(result,at_birth,modes,a_birth,birth_floor)}
    for mark, item in spec['markers'].items():
        # Linear interpolation on native knots. This output reports fields and separate work
        # integrals; the orbit calculation consumes the complete time series.
        split = np.searchsorted(t[keep],item['time'])
        left = split-1
        alpha = (item['time']-t[keep[left]])/(t[keep[split]]-t[keep[left]])
        endpoints=np.array([(1-alpha)*result[left]+alpha*result[split]])
        summaries[mark]=fields(endpoints,0,modes,a_birth,None)
    # Gauge-variation counterfactual at fixed boundary forcing and geometry. Choose the
    # same low modes used in the radial probe, with no change in incident waveform.
    response_controls = {}
    compare_modes = [(k,0) for k in ks]
    if label.endswith('half'):
        raw_before = traj[:birth_index+1]
        for which, constant in [('initial', K[0]), ('final', K[birth_index])]:
            fixed = raw_before.copy()
            fixed[:,2] = constant
            constant_vals = run_solver(exe,tmp,fixed,compare_modes,t_birth,[birth_index],label+'-K'+which)[0]
            response_controls['fixed_K_'+which] = {
                'K': float(constant),
                'opening_rms': (np.hypot(constant_vals[:,0],constant_vals[:,1])/math.sqrt(2.)).tolist()}
        # When was the seed deposited? Split the actual signed waveform by its cumulative
        # squared amplitude at the wall. A weak late component can offset an earlier kick.
        ti=t[:birth_index+1]
        profile=inc[:birth_index+1]
        signal=np.cumsum(np.diff(ti)*(profile[1:]**2+profile[1:]*profile[:-1]+profile[:-1]**2)/3.)
        pulse_splits={}
        split_modes=[(k,kind) for k in CURL for kind in (1,2)]
        for fraction in (.50,.90,.99,.9999):
            index=1+int(np.searchsorted(signal,fraction*signal[-1]))
            part=run_solver(exe,tmp,raw_before,split_modes,ti[index],[birth_index],label+'-seed')[0]
            vectors=part[:,:2].reshape(len(CURL),2,2)
            if np.max(abs(vectors.sum(axis=1)-whole[at_birth,:len(CURL),:2]))>1.e-7:
                raise ValueError('Seed epoch split failed superposition')
            pulse_splits[str(fraction)] = {
                'cut_time': float(ti[index]), 'enclosed_fraction':float(signal[index-1]/signal[-1]),
                'comoving_modes':CURL.tolist(), 'early_opening_EB':vectors[:,0].tolist(),
                'later_opening_EB':vectors[:,1].tolist()}
        response_controls['waveform_partition']=pulse_splits
        # Coarsen the exact same source history by two; retain the birth knot if odd.
        sparse_ids=np.unique(np.r_[np.arange(0,birth_index+1,2),birth_index])
        coarser=run_solver(exe,tmp,raw_before[sparse_ids],compare_modes,t_birth,[len(sparse_ids)-1],label+'-coarse')[0]
        field_ref=whole[at_birth,:,:2]
        absolute=np.linalg.norm(field_ref-coarser[:,:2],axis=1)
        relative=absolute/np.linalg.norm(field_ref,axis=1)
        response_controls['same_input_stride2']={'opening_EB_vector_absolute_error':absolute.tolist(),
                                                 'opening_EB_vector_relative_error':relative.tolist()}
        if np.max(relative[:5])>1.e-4 or np.max(absolute)>1.e-4:
            raise ValueError('Same-source resolution check failed')
    # Separate mode scan on log-spaced k and all eight named physical scales.
    scan = {}
    if dense:
        grid = np.unique(np.r_[np.geomspace(.1,100.,161)*a_birth,ks])
        scan_modes = [(k,0) for k in grid]
        scan_vals=run_solver(exe,tmp,traj,scan_modes,t_birth,[birth_index],label+'-sweep')[0]
        norm = np.hypot(scan_vals[:,0],scan_vals[:,1])/math.sqrt(2.)
        peak = int(np.argmax(norm))
        scan=dict(k_physical=(grid/a_birth).tolist(), opening_rms=norm.tolist(),
                  sampled_best_k_physical=float(grid[peak]/a_birth),sampled_best_rms=float(norm[peak]),
                  sampled_minimum_fraction_for_early_floor=float(birth_floor/norm[peak]))
    return dict(a_birth=a_birth, time_birth=t_birth, kappa2_EB_max_error=error.tolist(),
                kappa2_EB_max_relative_error=relative.tolist(),
                linear_split_max_absolute_error=split_error, maximum_energy_balance_relative_error=balance_error,
                initial_rms_floor_at_0_1_percent=birth_floor, named_comoving_modes=ks.tolist(),
                response_controls=response_controls, events=summaries, scan=scan)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prior',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--raw-dir',type=Path,help='reuse the reference captures written by capture.py')
    parser.add_argument('--jobs',type=int,default=2)
    args=parser.parse_args()
    prior=args.prior.resolve()
    out=args.output_dir.resolve()
    out.mkdir(parents=True,exist_ok=True)
    meta=json.loads((prior/'history/markers.json').read_text())
    limits=json.loads((prior/'global/CLUMPING.json').read_text())
    with tempfile.TemporaryDirectory() as folder:
        tmp=Path(folder)
        raw=args.raw_dir
        if raw is None:
            raw=tmp/'raw'
            runpy.run_path(str(HERE/'capture.py'))['capture'](prior,raw,args.jobs)
        exe=tmp/'evolve'
        subprocess.run(['g++','-std=c++17','-O3',str(HERE/'evolve.cpp'),'-o',str(exe)],check=True)
        validation=checks(exe,tmp)
        def one(label):
            return label,solve_case(prior,raw,tmp,out,exe,label,meta['runs'][label],limits,label.endswith('half'))
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            rows=dict(pool.map(one,meta['runs']))
    (out/'transfer.json').write_text(json.dumps(dict(controls=validation,runs=rows), sort_keys=True,
                                                indent=2,allow_nan=False)+'\n')

if __name__=='__main__':
    main()
