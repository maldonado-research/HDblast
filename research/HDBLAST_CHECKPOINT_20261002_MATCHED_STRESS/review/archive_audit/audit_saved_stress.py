#!/usr/bin/env python3
"""Post-run second implementation audit of saved modes and direct stress histories."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
import traceback
import numpy as np

LD = np.longdouble
NAMES = ('q','q_prime','q_second','rho','p','Q0','anomaly','current')
FREEZE = '4a5dad8dda6d0a57a9cf88cc8b4a9b8feeaea45d'
REGISTRATION = '4a1533fa7bf9df41cc635d76fa61a22b3aa119d8fbf215349fe09cf3ee54c893'
ARITHMETIC = LD('2e-11')


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def close(actual, expected, label, allowance=ARITHMETIC):
    actual, expected = np.asarray(actual), np.asarray(expected)
    require(actual.shape == expected.shape, label+': shape')
    require(np.all(np.isfinite(actual)) and np.all(np.isfinite(expected)), label+': nonfinite')
    error = np.max(np.abs(actual-expected), initial=LD(0))
    scale = max(LD(1), np.max(np.abs(expected), initial=LD(0)))
    require(error <= allowance+LD(1024)*np.finfo(LD).eps*scale, label+': arithmetic mismatch')
    return float(error)


def audit(checkpoint, modes_path):
    require(sha(checkpoint/'FULL_REGISTRATION.json') == REGISTRATION, 'Public registration anchor mismatch')
    registration = load(checkpoint/'FULL_REGISTRATION.json')
    for name, digest in registration['files'].items():
        path = (checkpoint/name).resolve()
        require(path.is_relative_to(checkpoint) and sha(path) == digest, 'Frozen file pin mismatch: '+name)
    experiment = load(checkpoint/'EXPERIMENT.json')
    gates = experiment['gates']
    modes = load(modes_path)
    require(modes['freeze_commit'] == FREEZE and modes['status'] == 'passed_internal_gates', 'Mode freeze/status')
    manifest_path = checkpoint/'independent/MANIFEST.json'
    require(modes['provenance']['manifest_sha256'] == sha(manifest_path), 'Mode manifest pin')
    require(modes['provenance']['files'] == load(manifest_path)['files'], 'Mode pin list')
    require(sys.version_info[:2] == (3,12) and np.__version__ == '2.2.6', 'Pinned audit runtime')
    require(np.finfo(LD).nmant >= 63, 'Extended arithmetic required')
    eps = LD(str(experiment['source']['epsilon']))
    summaries = {(r['source'],r['eta']):r for r in modes['rows']}
    expected = {(s,t) for s in experiment['source']['ids'] for t in experiment['observation_eta']}
    require(set(summaries) == expected and len(summaries) == len(modes['rows']), 'Summary grid')
    require(len(modes['runs']) == 4 and {(r['source'],r['setting']) for r in modes['runs']} ==
            {(s,t) for s in experiment['source']['ids'] for t in ('coarse','fine')}, 'Run coverage')
    report_runs = []
    ledger_lookup = {}
    maximums = {'phase_free_bare_operator':0.,'archived_integrand_reconstruction':0.,
               'momentum_sum_reconstruction':0.,'ledger_integrand_reconstruction':0.,
               'ledger_endpoint_reconstruction':0.,'fine_registered_ward_residual':0.,
               'fine_all_even_history_ward_residual':0.,'mode_amplitude_wronskian':0.,
               'panel_polynomial_moment':0.,'canonical_energy_integral':0.}
    n_snapshots = n_sums = n_ward = n_moments = 0
    for run in modes['runs']:
        source, setting = run['source'], run['setting']
        archive_path = (modes_path.parent/run['archive']['path']).resolve()
        require(archive_path.is_relative_to(modes_path.parent.resolve()) and
                sha(archive_path) == run['archive']['sha256'], 'Raw archive pin')
        fine = setting == 'fine'
        dt = LD(1)/(256 if fine else 128)
        width = LD(1)/(4 if fine else 2)
        steps, n_nodes = (1152,16384) if fine else (576,8192)
        item = {'source':source,'setting':setting,'archive_sha256':sha(archive_path),
                'snapshots':[],'registered_ward_endpoints':[]}
        with np.load(archive_path, allow_pickle=False) as arc:
            k, weights = arc['k'], arc['momentum_weights']
            require(k.dtype == weights.dtype == np.dtype(LD) and k.shape == weights.shape == (n_nodes,), 'Momentum dtype/grid')
            require(np.all(k > 0) and np.all(np.diff(k)>0) and np.all(weights>0), 'Momentum order/positivity')
            left = np.arange(n_nodes//16,dtype=LD)[:,None]*width
            z = 2*(k.reshape(-1,16)-left)/width-1
            local_weights = weights.reshape(-1,16)*2/width
            power = np.ones_like(z)
            for degree in range(32):
                target = LD(2)/(degree+1) if degree%2==0 else LD(0)
                moment_error = float(np.max(np.abs(np.sum(local_weights*power,axis=1,dtype=LD)-target)))
                require(moment_error < 1e-11, 'GL16 polynomial moment')
                maximums['panel_polynomial_moment'] = max(maximums['panel_polynomial_moment'],moment_error)
                n_moments += len(left)
                power *= z
            times, values, jets = arc['history_eta'],arc['history_values'],arc['history_source_jet']
            require(tuple(arc['history_quantity_names']) == NAMES, 'History quantity order')
            require(values.shape == (steps+1,3,8) and jets.shape == (steps+1,6), 'History shape')
            close(times, LD(-6)+np.arange(steps+1,dtype=LD)*dt, 'History time grid', LD(0))
            require(np.all(np.isfinite(values)) and np.all(np.isfinite(jets)), 'History nonfinite')
            a = -1/times
            d_prime = eps*(jets[:,1]-2*a*jets[:,0])/(a*a)
            close(arc['history_d_prime'],d_prime,'Mass derivative record')
            source_work = a[:,None]**4*values[:,:,5]*d_prime[:,None]/(2*eps)
            trace = -values[:,:,3]+3*values[:,:,4]
            flux = source_work-a[:,None]*trace
            residual = close(flux,arc['history_ledger_integrand'],'Saved Ward integrand')
            maximums['ledger_integrand_reconstruction'] = max(maximums['ledger_integrand_reconstruction'],residual)
            # Integrate disjoint two-step panels, then prefix-sum them. This
            # differs from the producer's global odd/even endpoint grouping.
            panels = dt/3*(flux[:-2:2]+4*flux[1:-1:2]+flux[2::2])
            cumulative = np.vstack((np.zeros((1,3),dtype=LD),np.cumsum(panels,axis=0,dtype=LD)))
            direct_even = values[::2,:,3]-values[0,:,3]
            all_time_error = float(np.max(np.abs(cumulative-direct_even)))
            item['all_even_history_ward_residual'] = all_time_error
            if fine:
                maximums['fine_all_even_history_ward_residual'] = max(maximums['fine_all_even_history_ward_residual'],all_time_error)
            source_free = times <= -5
            require(np.count_nonzero(jets[source_free]) == 0 and
                    np.count_nonzero(values[source_free][:,:,[0,1,2,3,4,6,7]]) == 0,
                    'Source-free initial history')
            post = times >= -3
            require(np.count_nonzero(jets[post]) == 0 and np.count_nonzero(values[post,:,6]) == 0,
                    'Postpulse local source/anomaly must vanish')
            close(values[post,:,7],values[post,:,0],'Postpulse current contact')
            for i,row in enumerate(run['rows']):
                eta = LD(str(row['eta'])); L = -1/eta
                index = int((eta-times[0])/dt)
                u,w = arc[f'u_{i}'],arc[f'w_{i}']
                jet = arc[f'source_jet_{i}']
                require(u.dtype == w.dtype == np.dtype(np.clongdouble) and u.shape == w.shape == k.shape, 'Mode dtype/shape')
                close(jet,jets[index],'Observation/history jet',LD(0))
                s = eps*jet[0]
                A, A1 = u.real/k,w.real/k
                # Exact phase-free forms of the TWO minimal complex operators.
                # Im(w) is retained as evolved, with no Wronskian replacement.
                variation_D = (k*k+L*L)*A-w.imag-L*A1
                bare_rho = (variation_D+(k*k+2*L*L)*A+s/(2*k))/(2*eps)
                bare_p = (variation_D-(k*k/3+2*L*L)*A-s/(2*k))/(2*eps)
                for key,term in (('rho',bare_rho),('p',bare_p)):
                    error = close(term,arc[f'bare_{key}_{i}'],'Phase-free bare '+key)
                    maximums['phase_free_bare_operator'] = max(maximums['phase_free_bare_operator'],error)
                measure = k*k/(2*np.arccos(LD(-1))**2)
                rebuilt = {'rho':measure*(bare_rho-arc[f'sub_rho_{i}']),
                           'p':measure*(bare_p-arc[f'sub_p_{i}'])}
                for key,array in rebuilt.items():
                    error = close(array,arc[f'integrand_{key}_{i}'],'Combined '+key)
                    maximums['archived_integrand_reconstruction'] = max(maximums['archived_integrand_reconstruction'],error)
                W = float(np.max(np.abs(2*u.real-w.imag/k))/eps)
                require(W <= gates['wronskian_over_epsilon'], 'Observation Wronskian')
                maximums['mode_amplitude_wronskian'] = max(maximums['mode_amplitude_wronskian'],W)
                canonical = measure*(2*k*u.real-w.imag)/(2*eps)
                if eta >= -3:
                    require(np.count_nonzero(jet) == 0 and
                            all(np.count_nonzero(arc[f'{name}_{i}']) == 0
                                for name in ('sub_rho','sub_p','delta_W2','delta_W4')),
                            'Postpulse stress contains a local subtraction/contact')
                snapshot = {'eta':float(eta),'wronskian_scaled':W,'postpulse':bool(eta>=-3),'finite_k':[]}
                for j,point in enumerate(row['finite_k']):
                    cutoff = point['K']; require(cutoff == (64,128,256)[j], 'Cutoff grid')
                    mask = k < cutoff
                    terms = {}
                    for name in NAMES:
                        array = rebuilt[name] if name in rebuilt else arc[f'integrand_{name}_{i}']
                        # Pairwise panel integration, followed by a reverse panel
                        # sum, independently groups the archived weighted nodes.
                        panel = np.sum((array[mask]*weights[mask]).reshape(-1,16),axis=1,dtype=LD)
                        value = np.sum(panel[::-1],dtype=LD)
                        error = close(np.array([value,value]),np.array([LD(point['values'][name]),values[index,j,NAMES.index(name)]]), 'Cutoff reconstruction '+name)
                        maximums['momentum_sum_reconstruction'] = max(maximums['momentum_sum_reconstruction'],error)
                        terms[name] = float(value)
                        n_sums += 1
                    canonical_integral = float(np.sum(canonical[mask]*weights[mask],dtype=LD))
                    maximums['canonical_energy_integral'] = max(maximums['canonical_energy_integral'],abs(canonical_integral))
                    ledger = cumulative[index//2,j]
                    ledger_lookup[(source,setting,float(eta),cutoff)] = float(ledger)
                    ledger_error = close(np.array([ledger]),np.array([LD(point['ward']['ledger'])]),'Ward endpoint grouping')
                    maximums['ledger_endpoint_reconstruction'] = max(maximums['ledger_endpoint_reconstruction'],ledger_error)
                    ward_residual = float(abs(ledger-values[index,j,3]))
                    if fine:
                        require(ward_residual <= gates['ward_endpoint'], 'Fine registered Ward gate')
                        maximums['fine_registered_ward_residual'] = max(maximums['fine_registered_ward_residual'],ward_residual)
                    item['registered_ward_endpoints'].append({'eta':float(eta),'K':cutoff,'residual':ward_residual})
                    snapshot['finite_k'].append({'K':cutoff,'values':terms,'canonical_energy_variation_integral':canonical_integral})
                    n_ward += 1
                item['snapshots'].append(snapshot)
                n_snapshots += 1
        report_runs.append(item)
    refinements = []
    for (source,eta),row in summaries.items():
        for point in row['finite_k']:
            cutoff = point['K']
            gap = abs(ledger_lookup[(source,'fine',eta,cutoff)]-ledger_lookup[(source,'coarse',eta,cutoff)])
            require(gap <= gates['ward_refinement'], 'Registered ledger refinement gate')
            close(np.array([gap]),np.array([point['ward']['refinement_difference']]),'Ward refinement record')
            refinements.append({'source':source,'eta':eta,'K':cutoff,'ledger_refinement':gap})
    return {'status':'passed','classification':'post_run_saved_data_audit',
            'provenance':{'freeze_commit':FREEZE,'registration_sha256':REGISTRATION,
                          'mode_results_sha256':sha(modes_path),'audit_source_sha256':sha(__file__)},
            'runtime':{'python':platform.python_version(),'numpy':np.__version__,
                       'longdouble_nmant':int(np.finfo(LD).nmant),'python_optimization':sys.flags.optimize},
            'counts':{'mode_snapshots':n_snapshots,'reconstructed_quantity_sums':n_sums,
                      'registered_ward_endpoints':n_ward,'ward_refinements':len(refinements),
                      'panel_polynomial_moments':n_moments},
            'maximums':maximums,'runs':report_runs,'ward_refinements':refinements,
            'post_run_arithmetic_allowance':float(ARITHMETIC),
            'limits':['No new source sampling, mode evolution, or time/momentum-response experiment was performed.',
                      'Phase-free operators use true evolved u,w and are independent of the producer complex-phase multiplication.',
                      'Complete subtraction archives are hash-bound and separately audited by the frozen validator; this second audit reuses them.',
                      'Only six mode snapshots per run are stored. Complete histories store direct integrated stresses, not all intermediate modes.',
                      'All-even-time ledger residuals outside registered observations are descriptive post-run diagnostics, not new scientific gates.',
                      'Agreement and refinement estimate integration accuracy; only the separately derived omitted UV band has an analytic enclosure.',
                      'Postpulse linear stress is coherent response, not canonical occupation energy, particle yield or heating.']}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--checkpoint',required=True,type=Path)
    p.add_argument('--modes',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    args=p.parse_args()
    require(not args.output.exists(),'Refusing to overwrite post-run audit')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    report={'status':'started','audit_source_sha256':sha(__file__)}
    try:
        report=audit(args.checkpoint.resolve(),args.modes.resolve())
    except Exception as exc:
        report.update({'status':'failed','exception':str(exc),'traceback':traceback.format_exc()})
        raise
    finally:
        args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':report['status'],'counts':report['counts'],'maximums':report['maximums']}))


if __name__=='__main__':
    main()
