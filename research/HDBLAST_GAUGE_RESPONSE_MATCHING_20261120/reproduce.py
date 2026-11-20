"""Reproduce the photon response matching calculations on the calibrated closed descents."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
BASE = "dbdf3fe1e62981d3b5a92081a0dc9eb14029e21e"
PINS_SHA = "794114aa93db455df2e0f97bb3f8018b73774479a49c5708af9fca5570dd6f77"
PRODUCTS = {name: f"{name}/{name.upper()}.json" for name in ("matched", "probe", "band", "completion", "reach", "compressed")}
SCRIPTS = {name: f"{name}/{name}.py" for name in PRODUCTS}


def require(value, message):
    if not value:
        raise ValueError(message)


def read(path):
    return json.loads(path.read_bytes())


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def verify_files():
    manifest = HERE/'SHA256SUMS'
    seen = set()
    for line in manifest.read_text().splitlines():
        checksum, name = line.split('  ', 1)
        path = HERE/name
        require(not Path(name).is_absolute() and Path(name).as_posix() == name
                and path == path.resolve() and path.is_relative_to(HERE)
                and path != manifest and path not in seen, 'Invalid inventory path')
        require(digest(path) == checksum, f'Hash differs for {name}')
        seen.add(path)
    actual = {p for p in HERE.rglob('*') if p.is_file() and p != manifest
              and not (p.parent.name == '__pycache__' and p.suffix == '.pyc')}
    require(actual == seen, 'Archive coverage differs')
    pin_file = HERE/'SOURCES.json'
    require(digest(pin_file) == PINS_SHA, 'Unexpected source manifest')
    pins = read(pin_file)
    require(pins['base'] == BASE and pins['previous_doi'] == '10.5281/zenodo.24382645',
            'Historical revision differs')
    prior = set()
    for record in pins['sources']:
        path = HERE/record['bundled']
        require(path in seen and path not in prior and path.stat().st_size == record['bytes']
                and digest(path) == record['sha256'], 'Historical input mismatch')
        prior.add(path)
    require(prior == {p for p in seen if p.is_relative_to(HERE/'prior')}, 'Historical coverage differs')
    print(f'Verified {len(seen)} files and {len(prior)} historical objects.', flush=True)


def close_number(actual, expected, path):
    # Absolute floors apply to independently measured cancellation/solver diagnostics.
    # Physical endpoints and spectral amplitudes retain a relative comparison.
    product, field = path[0], path[-1]
    absolute = 0.
    if product in ('matched','reach'):
        absolute = {'scalar_residual': 4e-9, 'constraint': 4e-9,
                    'purity_residual': 3e-10, 'work_residual': 1e-13}.get(field,0.)
        if 'free_intervals' in path or 'strong_ungated_intervals' in path:
            absolute = {'constraint_error': 6e-12,'parent_balance_error': 3e-13,
                        'normalized_comoving_density_error': 2e-15,'independent_maximum_difference': 4e-11,
                        'feedback_first_minimum_difference': 4e-12}.get(field,0.)
            if product=='reach' and field in ('constraint_error','independent_maximum_difference'):
                absolute=6e-9
            if path[-2:]==('stop','h'):
                absolute=1e-14
        if 'independent_coordinate_differences' in path or 'independent_trajectory_error' in path:
            absolute = 2e-10
    elif product=='compressed':
        absolute = {'max_relative_moment_error': 1e-13,'relative_error': 1e-14,
                    'constraint': 3e-11,'scalar_residual': 3e-11,'work_residual': 2e-14,
                    'purity_residual': 1e-13}.get(field,0.)
        if 'archived_differences' in path or 'first_minimum_fractional_errors' in path:
            absolute = 2e-12
        elif 'independent_n8_coordinate_differences' in path:
            absolute = 2e-11
    elif product=='probe':
        absolute = {'outgoing_from_two_local_terms': 2e-10,'relative_outgoing_born_error': 1e-7,
                    'post_switch_fractional_energy_drift': 2e-11,'relative_l2_force_error': 2e-14}.get(field,0.)
    elif product=='completion':
        absolute = {'maximum_purity_error': 1e-16,'unweighted_historical_occupation_error': 2e-13,
                    'resolution_absolute_D_error': 2e-8,'signed_numerical_duality_gap_D': 1e-10}.get(field,0.)
    elif product=='band':
        absolute = {'clock_error': 1e-15, 'wronskian_error': 2e-14,'abs_work_balance_error': 9e-7,
                    'absolute_local_integrated_work_2': 4e-9,'absolute_local_integrated_work_4': 1e-10,
                    'archived_band_D_error': 2e-10,'integrated_exact_work': 9e-7,
                    'after_fourth_order_rms_over_exact': 1e-12}.get(field,0.)
        if '385-768' in path and 'delta_0.1_stretch_4' in path and field in ('D_exact','exit_D'):
            absolute = 8e-18
    return math.isclose(actual, expected, rel_tol=2e-7, abs_tol=absolute)


def compare(actual, expected, path=()):
    label = '/'.join(path)
    if isinstance(expected, dict):
        require(isinstance(actual,dict) and actual.keys()==expected.keys(), f'Different keys at {label}')
        for key in expected:
            compare(actual[key],expected[key],(*path,key))
    elif isinstance(expected, list):
        require(isinstance(actual,list) and len(actual)==len(expected), f'Different size at {label}')
        for i, (a,e) in enumerate(zip(actual,expected)):
            compare(a,e,(*path,str(i)))
    elif isinstance(expected, float):
        require(type(actual) in (int,float) and math.isfinite(actual) and math.isfinite(expected)
                and close_number(actual, expected, path),
                f'{label}: expected {expected!r}, obtained {actual!r}')
    else:
        require(type(actual) is type(expected) and actual==expected, f'Scalar differs at {label}')


def equal(a, b, label, *, rel=5e-12, absolute=0.):
    require(math.isclose(a, b, rel_tol=rel, abs_tol=absolute), f'{label}: {a} differs from {b}')


def check_probe(probe):
    ref = probe['reference']
    pole = 2*(ref+1)
    equal(probe['first_added_pair_frequency'], pole, 'First added pole')
    cal = (probe['matching_frequencies'][1]/2)**2
    for frequency, entry in probe['frequencies'].items():
        x = (float(frequency)/2)**2
        require(0 < x < (pole/2)**2, 'Frequency outside the subthreshold strip')
        for cut, row in entry['cutoffs'].items():
            modes = range(ref+1, int(cut)+1)
            weight = lambda l: l-1/l
            equal(row['chi'], math.fsum(weight(l)/(1-x/l**2) for l in modes), 'Static-shell resolvent')
            equal(row['static_coefficient'], math.fsum(weight(l) for l in modes), 'Zeroth dispersive moment')
            equal(row['derivative_coefficient'], math.fsum(weight(l)/l**2 for l in modes), 'First dispersive moment')
            slope = math.fsum(weight(l)*x**2/(l**2*(l**2-x)) for l in modes)
            secant = math.fsum(weight(l)*x*(x-cal)/((l*l-x)*(l*l-cal)) for l in modes)
            equal(row['after_value_and_slope'], slope, 'Slope remainder')
            equal(row['after_two_frequencies'], secant, 'Secant remainder')
            for partial, full, bound in [
                (slope,'infinite_after_value_and_slope','uncomputed_slope_tail_bound'),
                (secant,'infinite_after_two_frequencies','uncomputed_two_frequency_tail_bound')]:
                require(abs(entry[full]-partial) <= row[bound]+1e-9, 'Subtracted tail missed the exact infinite value')
    pulse = probe['compact_pulse']['cutoffs']
    outgoing = 'outgoing_born_a4_density_per_amplitude_squared'
    last = pulse[max(pulse, key=int)][outgoing]
    for row in pulse.values():
        energy = row[outgoing]
        require(energy > 0 and 0 <= last-energy <= row['omitted_outgoing_born_tail_upper'], 'Born energy tail')
        equal(row['outgoing_from_retarded_work'], energy, 'Pulse work', absolute=1e-7)
        require(abs(row['outgoing_from_two_local_terms']) < 1e-7, 'Local pulse endpoint work')
        equal(row['peak_free_signal_after']/row['peak_retarded_signal_during'],
              row['maximum_tail_signal_over_peak'], 'Free pulse response')


def check_state(record, margin, min_key):
    end = record['end']
    a, h, z, v = (end[k] for k in ('a','h','z','v'))
    photons, scalar = end['photon_density'], end['scalar_density']
    require(a>0 and photons >= 0 and scalar + photons > 0, 'Invalid endpoint densities')
    require(margin-1e-8 <= record[min_key] <= z+1e-8, 'Event or minimum violates the reserved inertia')
    equal(h*h + 1/a**2, photons+scalar, 'Friedmann constraint', absolute=8e-7)
    equal(end['scalar_pressure']+scalar, (2*z+2*(1-z)/3)*v*v/6, 'Pressure from the matched action')
    equal(photons/(photons+scalar), end['photon_fraction'], 'Physical photon fraction')
    equal(a**4*photons, end['comoving_energy'], 'Comoving energy')
    equal(end['reference_band_comoving'] + end['added_band_comoving'], end['comoving_energy'], 'Shell energy decomposition')
    equal(end['reference_band_comoving']/a**4/(photons+scalar), end['reference_band_fraction'], 'Low-band photon fraction')
    if record['status'] == 'kinetic_reserve':
        equal(z, margin, 'Terminal inertia', absolute=1e-7)
    elif record['status'] == 'first_minimum':
        equal(end['progress'], 1., 'Terminal field minimum')
    else:
        require(record['status'] == 'first_turnaround' and abs(h)<1e-7, 'Unexpected terminal state')


def check_coast(coast, live, invariant):
    for point in ('start','stop','peak','first_minimum'):
        state = coast[point]
        if state is None:
            continue
        a = state['a']
        radiation = state['photon_density']
        equal(state['fraction_of_matter'], radiation/(radiation+state['parent_density']), 'Free photon share')
        equal(radiation*a**4, invariant, 'Free-streaming photon redshift', absolute=2e-8)
        equal(state['h']**2+1/a**2, radiation+state['parent_density'], 'Free Friedmann constraint', absolute=2e-8)
    if coast['first_minimum'] is not None and live['status']=='first_minimum':
        equal(coast['first_minimum']['fraction_of_matter'], live['end']['photon_fraction'], 'Independent continuing solution', absolute=2e-7)
    intervals = coast['half_intervals']
    for sub in intervals:
        equal(sub['end_x']-sub['begin_x'], sub['proper_duration'], 'Time above one half')
        equal(math.log(sub['end_a']/sub['begin_a']), sub['delta_N'], 'Expansion above one half')
    equal(sum(sub['delta_N'] for sub in intervals), coast['total_efolds_at_or_above_half'], 'Integrated expansion')
    equal(sum(sub['proper_duration'] for sub in intervals), coast['total_proper_duration_at_or_above_half'], 'Integrated duration')


def check_live(matched, reach, gauge, previous, probe):
    ref, margin = matched['reference_cutoff'], matched['reserve']
    require(ref == reach['reference_cutoff'] == probe['reference'] and margin == reach['reserve'], 'Different calibration')
    for key, archive in gauge['branches'].items():
        points = archive['points']
        dq = points['first_minimum']['q'] - points['launch']['q']
        width = points['switch_end']['q'] - points['switch_start']['q']
        branch = matched['branches'][key]
        cap = reach['branches'][key]
        equal(abs(width), cap['absolute_field_width'], 'Historical switching width')
        def coefficient(cut):
            raw = math.fsum(l-1/l for l in range(ref+1, cut+1))
            probe_rows = probe['frequencies']['8']['cutoffs']
            if str(cut) in probe_rows:
                equal(raw, probe_rows[str(cut)]['static_coefficient'], 'Live versus fixed-band moment')
            return raw * archive['H0_over_Mpl']**2/(2*math.pi**2)
        def gradient(progress, delta):
            t = (progress-.2)/.6
            return delta*math.exp(-1/(t*(1-t)))/(.6*dq*gauge['bump_normalization']) if 0<t<1 else 0.
        for group in ('runs', 'uncalibrated'):
            for run_id, row in branch[group].items():
                cut, delta = row['cutoff'], row['delta']
                c = coefficient(cut) if group == 'runs' else 0.
                end = row['end']
                check_state(row, margin, 'minimum_bare_inertia')
                equal(c, row['added_inertia_coefficient'], 'Live inertia coefficient')
                equal(1.-c*(gradient(end['progress'],delta)/end['a'])**2, end['z'], 'Live physical inertia', absolute=3e-12)
                if (group == 'uncalibrated' or cut==ref) and run_id in previous[key]['runs']:
                    past = previous[key]['runs'][run_id]['end']
                    for field in ('x','a','h','v','photon_density'):
                        equal(end[field], past[field], 'Historical branch at '+run_id, absolute=2e-8)
        for run_id, coast in branch['free_intervals'].items():
            live = branch['runs'][run_id]
            switch = live['switch_end']
            equal(coast['start']['x'],switch['x'],'Start of switch-free evolution')
            check_coast(coast, live, switch['comoving_energy'])
        for cutoff, coast in cap['strong_ungated_intervals'].items():
            live = cap['tapered_profile_runs'][cutoff+'/3']
            check_coast(coast, live, live['end']['comoving_energy'])
        for cutoff, row in cap['cutoffs'].items():
            c = coefficient(int(cutoff))
            equal(c, row['c'], 'Capacity from the same shell action')
            scale = math.sqrt((1-margin)/c)
            upper = scale*abs(width)*cap['mean_reference_radius']
            equal(row['upper_contrast_on_reference'], upper, 'Integrated contrast ceiling')
            for name in ('original','tapered'):
                equal(row[name+'_contrast_limit'], upper*cap[name+'_capacity_fraction'], name+' capacity')
            radius_end = 3./(scale*abs(width))
            radius_mid = abs(gradient(.5,3.))/scale
            for point, name, radius in [('end','any_profile',radius_end),('midpoint','original',radius_mid)]:
                equal(row[f'strong_minimum_{point}_radius_{name}'],radius,'Necessary expansion radius')
                time_bound = math.acosh(max(1.,radius/points['launch']['a']))/math.sqrt(points['launch']['parent_over_V0'])
                equal(row[f'strong_earliest_{point}_{name}'], time_bound, 'Earliest admissible expansion')
        for row in cap['tapered_profile_runs'].values():
            check_state(row,margin,'minimum_z')


def check_band(band, gauge, probe):
    for name, match in band['regulator_matches'].items():
        if name.startswith('two_moment_'):
            other = band['regulator_matches']['two_moment_plus']
            for field in ('W2','W4'):
                equal(match[field],other[field],'Common two-moment window')
            continue
        k = str(match['sharp_maximum'])
        static = probe['frequencies']['8']['cutoffs'][k]['static_coefficient']
        equal(match['sharp_W2'], static/(4*math.pi**2), 'Band and probe static weights')
        equal(match['tapered_W2'], match['sharp_W2'], 'Smooth regulator calibration')
    for label, branch in band['branches'].items():
        spectrum = gauge['branches'][label]['primary_spectrum']
        for case in branch['cases'].values():
            rows = case['bands']
            for name,row in rows.items():
                equal(row['points']['1.0']['D_exact'], row['exit_D'], 'Endpoint occupation')
                equal(row['exit_D'], row['integrated_exact_work'], 'Band work',
                      absolute=max(2e-9,1e-8*row['peak_D_sampled']))
                require(max(row['absolute_local_integrated_work_2'],row['absolute_local_integrated_work_4'])<3e-5,
                        'Unexpected local outgoing work')
                if case['delta_log_I']==gauge['primary_delta'] and case['conformal_stretch']==1 and name[0].isdigit():
                    lo,hi = map(int,name.split('-'))
                    historical_D = math.fsum(l*(l*l-1)*spectrum['occupation'][l-2]/math.pi**2 for l in range(lo,hi+1))
                    equal(row['exit_D'], historical_D,'Historical fixed-band particle energy',absolute=2e-7)
            for point in ('0.25','0.5','0.75','1.0'):
                for field in rows['49-768']['points'][point]:
                    parts = [rows[k]['points'][point][field] for k in ('49-384','385-768')]
                    equal(math.fsum(parts), rows['49-768']['points'][point][field], 'Independent additive shell windows',absolute=1e-9)


def check_completion(data, gauge):
    ref, last = data['common_low_band'], data['highest_shell']
    modes = range(2,last+1)
    targets = data['high_band_matched_coefficients']
    def survival(cuts, probabilities, l):
        return math.fsum(p for k,p in zip(cuts,probabilities,strict=True) if k>=l)
    def moment(cuts, probabilities, j):
        return math.fsum(survival(cuts,probabilities,l)*l*(l*l-1)/l**(2*j) for l in range(ref+1,last+1))
    def weights(cuts, probabilities):
        require(all(ref <= k <= last for k in cuts) and all(0<=p<=1 for p in probabilities), 'Invalid cutoff weights')
        equal(sum(probabilities),1.,'Cutoff probabilities',absolute=1e-10)
    for i, name in enumerate(('A','B')):
        cuts, prob = data['cuts'], data['mixture_'+name]
        weights(cuts,prob)
        for j,target in enumerate(targets,1):
            equal(moment(cuts,prob,j),target,'Matched cutoff moments')
        equal(moment(cuts,prob,3),data['high_band_third_moments'][i],'Third moment')
        for bound in data['high_band_dispersive_bounds']:
            x = (bound['conformal_frequency']/2)**2
            value = math.fsum(survival(cuts,prob,l)*(l-1/l)/(1-x/l**2) for l in range(ref+1,last+1))/targets[0]
            equal(value,bound['normalized_'+name],'Dispersive profile')
            require(bound['lower'] <= value <= bound['upper'],'Positive response outside analytic envelope')
    for label, branch in data['branches'].items():
        archive = gauge['branches'][label]
        spectrum = archive['primary_spectrum']
        born = {l:l*(l*l-1)*spectrum['occupation_born'][l-2]/math.pi**2 for l in modes}
        low = math.fsum(l*(l*l-1)*spectrum['occupation'][l-2]/math.pi**2 for l in range(2,ref+1))
        def born_energy(cuts,probabilities,start):
            return math.fsum(survival(cuts,probabilities,l)*born[l] for l in range(start,last+1))
        for name,row in branch['profiles'].items():
            prob = data['mixture_'+name]
            equal(row['D_Born'], born_energy(data['cuts'],prob,2), 'Historical total Born output')
            equal(row['high_band_D_Born'], born_energy(data['cuts'],prob,ref+1), 'Historical high-band Born output')
            equal(row['D']-row['high_band_D'], low, 'Common low-band output', absolute=2e-7)
            equal(row['photon_over_original_parent'], row['D']*archive['points']['first_minimum']['photon_to_parent_per_D'],
                  'Fixed-trajectory photon share')
        for number, limits in branch['born_moment_ranges'].items():
            for row in limits.values():
                cuts,prob = row['cuts'],row['probabilities']
                weights(cuts,prob)
                for j in range(1,int(number)+1):
                    equal(moment(cuts,prob,j), targets[j-1], 'Feasible extremal moments', rel=2e-9)
                equal(row['total_D_Born'],born_energy(cuts,prob,2), 'Extremal Born output',absolute=1e-7)
                require(abs(row['signed_numerical_duality_gap_D'])<2e-7,'Relaxed dual certificate')
        for field,ratio in [('D','ratio_total_A_over_B'), ('high_band_D','ratio_high_A_over_B')]:
            equal(branch['profiles']['A'][field]/branch['profiles']['B'][field],branch[ratio],'Completion output ratio')


def check_compressed(data, probe, previous, matched):
    ref, top = data['reference_cutoff'],data['target_cutoff']
    require(ref==matched['reference_cutoff'] and top>ref, 'Reduced band differs from matched reference')
    modes = range(ref+1,top+1)
    for order,rule in data['rules'].items():
        labels,weights=rule['effective_frequencies'],rule['inertia_weights']
        require(len(labels)==len(weights)==int(order) and sorted(labels)==labels
                and min(labels)>=ref+1-1e-10 and max(labels)<=top+1e-10 and min(weights)>0,
                'Positive Gaussian frequencies or weights left their band')
        first = ref+1
        for j in range(2*int(order)):
            exact = math.fsum((l-1/l)*(first/l)**(2*j) for l in modes)
            approx = math.fsum(w*(first/l)**(2*j) for l,w in zip(labels,weights,strict=True))
            equal(exact,approx,'Finite positive quadrature moment',rel=1e-11)
        for omega, row in rule['static_frequency_response'].items():
            x = (float(omega)/2)**2
            if omega in probe['frequencies']:
                equal(probe['frequencies'][omega]['cutoffs'][str(top)]['chi'],row['exact'],'Independent retarded response')
            exact = math.fsum((l-1/l)/(1-x/l**2) for l in modes)
            approx = math.fsum(w/(1-x/l**2) for l,w in zip(labels,weights,strict=True))
            equal(exact,row['exact'],'Full susceptibility')
            equal(approx,row['approximate'],'Positive quadrature susceptibility')
    for label, cases in data['branches'].items():
        for delta, row in cases.items():
            run_id=f'{top}/{delta}'
            historical = previous[label]['runs'][run_id]['end']
            full = row['full']['end']
            equal(full['photon_fraction'],historical['fraction_of_matter'],'Full-band photon fraction',absolute=3e-8)
            equal(full['comoving_energy'],historical['comoving_excess'],'Full-band comoving energy',absolute=3e-8)
            if run_id in matched['branches'][label]['uncalibrated']:
                live = matched['branches'][label]['uncalibrated'][run_id]['end']
                for field in ('a','x','photon_fraction','added_band_comoving','reference_band_comoving'):
                    equal(live[field],full[field],'Unmatched reference for reduced dynamics',absolute=3e-8)
            for record in row['rules'].values():
                for stage in ('switch_end','end'):
                    state=record[stage]
                    total=state['photon_density']+state['scalar_density']
                    equal(state['h']**2+state['a']**-2,total,'Reduced Friedmann constraint',absolute=8e-7)
                    equal(state['photon_density']/total,state['photon_fraction'],'Reduced photon fraction')
                    equal(state['added_band_comoving']+state['reference_band_comoving'],state['comoving_energy'],'Reduced shell sum')
                equal(record['end']['comoving_energy'],record['switch_end']['comoving_energy'],'Reduced free-streaming energy',absolute=1e-8)
                for field,reported in record['first_minimum_fractional_errors'].items():
                    equal(reported,(record['end'][field]-full[field])/full[field],'Reported finite-band compression error',absolute=5e-15)


def physics(products):
    require(products.keys() == PRODUCTS.keys(), 'Missing numerical product')
    gauge = read(HERE/'prior/gauge/GAUGE.json')
    previous = read(HERE/'prior/feedback/FEEDBACK.json')['branches']
    check_probe(products['probe'])
    check_live(products['matched'], products['reach'], gauge, previous, products['probe'])
    check_band(products['band'], gauge, products['probe'])
    check_completion(products['completion'], gauge)
    check_compressed(products['compressed'], products['probe'], previous, products['matched'])


def regenerate(destination, jobs):
    destination = destination.resolve()
    require(not destination.is_relative_to(HERE),'Use a separate output directory')
    env = dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
    def run(name):
        output = destination/PRODUCTS[name]
        output.parent.mkdir(parents=True,exist_ok=True)
        require(not output.exists(),f'Output already exists: {output}')
        command = [sys.executable,'-I','-B',str(HERE/SCRIPTS[name]),'--output',str(output)]
        print('Recomputing '+SCRIPTS[name],flush=True)
        proc = subprocess.run(command,cwd=HERE,env=env,stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT,text=True)
        require(proc.returncode==0, f"{name} exited {proc.returncode}\n{proc.stdout[-4000:]}")
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        list(pool.map(run,PRODUCTS))
    data = {name:read(destination/path) for name,path in PRODUCTS.items()}
    for name,path in PRODUCTS.items():
        compare(data[name],read(HERE/path),(name,))
    physics(data)
    print(f'Verified {len(PRODUCTS)} numerical products and cross-product physical checks.',flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify',action='store_true',help='Repeat the numerical experiments')
    parser.add_argument('--jobs',type=int,default=2)
    parser.add_argument('--output-dir',type=Path,help='Keep regenerated files in this separate directory')
    args=parser.parse_args()
    require(args.jobs>=1 and (args.output_dir is None or args.verify),'Incompatible replay arguments')
    verify_files()
    if args.verify:
        if args.output_dir:
            regenerate(args.output_dir,args.jobs)
        else:
            with TemporaryDirectory(prefix='hdblast-suscept-') as tmp:
                regenerate(Path(tmp),args.jobs)
