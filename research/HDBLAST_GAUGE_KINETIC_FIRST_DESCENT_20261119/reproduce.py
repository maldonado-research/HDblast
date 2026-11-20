"""Reproduce the gauge-mode, radiation-fluid, analytic-control and bulk-origin calculations."""
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
BASE = "3309eb472d7fea083ea6bdf8eebfdfd276ee2f08"
PINS_SHA = "39903b9813055451d91d2c27f0e49a029cddaef143441e88c5ab4faf102ed729"
PRODUCTS = {"gauge": "gauge/GAUGE.json", "checks": "checks/CHECKS.json", "fluid": "fluid/FLUID.json",
            "origin": "origin/ORIGIN.json", "feedback": "feedback/FEEDBACK.json"}
SCRIPTS = {"gauge": "gauge/gauge.py", "checks": "checks/checks.py", "fluid": "fluid/fluid.py",
           "origin": "origin/origin.py", "feedback": "feedback/feedback.py"}


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
    require(pins['base'] == BASE and pins['previous_doi'] == '10.5281/zenodo.24381850',
            'Historical revision differs')
    prior = set()
    for record in pins['sources']:
        path = HERE/record['bundled']
        require(path in seen and path not in prior and path.stat().st_size == record['bytes']
                and digest(path) == record['sha256'], 'Historical input mismatch')
        prior.add(path)
    require(prior == {p for p in seen if p.is_relative_to(HERE/'prior')}, 'Historical coverage differs')
    print(f'Verified {len(seen)} files and {len(prior)} historical objects.', flush=True)


def numerical_floor(path):
    """Absolute replay tolerances for small diagnostics, in their reported units."""
    product, key = path[0], path[-1]
    if product == 'gauge' and key.endswith('_error'):
        return 2e-12
    if product == 'checks' and 'checks' in path[1:]:
        return 1e-9 if key == 'max_relative_occupation_error' else 5e-13
    if product == 'fluid':
        if 'balances' in path or 'absolute_discrepancies' in path or key=='unsourced_historical_maximum_difference':
            return 5e-10
        if key in {'h','fraction_rate','instantaneous_Q'}:
            return 2e-10
    if product == 'origin':
        if key == 'static_refinement_max_absolute' or path[-2:] == ('static_refined','phi_b'):
            return 4e-9
        if key.endswith(('_error','_max_shift')):
            return 5e-12
    if product == 'feedback':
        if 'free_intervals' in path and (key == 'h' or key.endswith(('_error','_difference'))):
            return 5e-12
        if key == 'purity_residual':
            return 3e-9
        if key in {'constraint','source_residual','parent_residual'}:
            return 2e-10
        if 'switch_start' in path:
            # Event interpolation can sample a tiny nonzero rate at the flat endpoint.
            return 1e-26
    return 0.


def close_number(actual, expected, path):
    if path[0] == 'gauge' and ('selected_occupations' in path or
            any(key in path for key in ('occupation','occupation_born'))):
        # The integrators evolve amplitudes. Near a node, a small absolute
        # amplitude error can produce a large relative change in its square.
        if actual < 0 or expected < 0:
            return False
        if '0.0' in path:
            return actual == expected == 0
        return math.isclose(math.sqrt(actual),math.sqrt(expected),rel_tol=1e-7,abs_tol=8e-14)
    return math.isclose(actual,expected,rel_tol=2e-7,abs_tol=numerical_floor(path))


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


def verify_gauge(gauge, historical):
    """Reconstruct a partial energy sum and the outgoing-to-parent energy ratios."""
    require(gauge['branches'].keys() == historical.keys(), 'Different rolling backgrounds')
    delta = str(gauge['primary_delta'])
    for label, branch in gauge['branches'].items():
        points = branch['points']
        expected_crossing = historical[label]['runs']['inverse_f_seed_0']['first_minimum_crossing']
        for field in ('x', 'a', 'q', 'h'):
            require(abs(points['first_minimum'][field]-expected_crossing[field]) < 2e-8,
                    'Gauge background missed an archived crossing')
        scan = branch['scan']
        require(scan['0.0']['D'] == 0 and scan[delta]['D'] > 0,
                'Photon energy in the vacuum or primary preparation differs')
        require(math.isclose(scan[delta]['D'], scan[str(-float(delta))]['D'], rel_tol=1e-9),
                'Opposite gauge steps yield different energy')
        spectrum = branch['primary_spectrum']
        first, last = spectrum['first_shell'], spectrum['last_shell']
        require(first == 2 and last == gauge['maximum_shell'], 'Primary spectrum misses a shell')
        shells = range(first, last+1)
        for field in ('occupation','occupation_born'):
            require(len(spectrum[field]) == len(shells) and all(v >= 0 for v in spectrum[field]),
                    'Incomplete spectrum or negative photon occupation')
        total = math.fsum(k*(k*k-1)*v/math.pi**2 for k,v in zip(shells,spectrum['occupation']))
        require(math.isclose(total, scan[delta]['D'], rel_tol=2e-12),
                'Transverse degeneracy or partial energy sum differs')
        for row in scan.values():
            partial = [v for _, v in sorted(row['D_by_cutoff'].items(), key=lambda item: int(item[0]))]
            require(all(0 <= small <= large for small, large in zip(partial, partial[1:]))
                    and math.isclose(partial[-1], row['D'], rel_tol=2e-13),
                    'Truncated outgoing energy is inconsistent')
            for stage, ratio in row['photon_to_parent'].items():
                point = points[stage]
                predicted = row['D'] * branch['H0_over_Mpl']**2 / (3*point['a']**4*point['parent_over_V0'])
                require(math.isclose(predicted, ratio, rel_tol=2e-12), 'Photon energy normalization differs')


def verify_benchmarks(checks):
    sech = checks['sech']
    amplitude = math.sinh(math.pi*sech['g'])**2
    def occupation(k):
        return amplitude/math.cosh(math.pi*sech['tau']*k)**2
    for shell in sech['samples']:
        expect = occupation(shell['k'])
        require(math.isclose(shell['exact'], expect, rel_tol=5e-13)
                and math.isclose(shell['canonical_ode'], expect, rel_tol=1e-7),
                'Analytic scattering benchmark disagrees')
    weighted = math.fsum(k*(k*k-1)*occupation(k) for k in range(2,sech['highest_k']+1))/math.pi**2
    require(math.isclose(weighted,sech['full_sphere_a4_density'],rel_tol=2e-12),
            'Exact spectrum missed a transverse state degeneracy')
    sums = checks['compact']['partial_sums']
    for lo,hi in zip(sums,sums[1:]):
        require(0 <= hi['sum_a4_density']-lo['sum_a4_density'] <= lo['rigorous_tail_bound_rounded']
                and 0 <= hi['rigorous_tail_bound_rounded'] < lo['rigorous_tail_bound_rounded'],
                'Compact analytic energy tail failed to enclose the partial sums')


def verify_fluid(fluid, gauge):
    require(fluid['branches'].keys() == gauge['branches'].keys(), 'Fluid backgrounds differ')
    for label, branch in fluid['branches'].items():
        runs = branch['scan']
        zero = runs['0']['first_minimum']
        target = gauge['branches'][label]['points']['first_minimum']
        require(zero['r'] == 0 and abs(zero['x']-target['x']) < 2e-8, 'Unsourced fluid missed the gauge roll')
        for which, thresholds in branch['half_thresholds'].items():
            require(len(thresholds) == 2 and thresholds[0]['gamma'] < thresholds[1]['gamma'],
                    'Fluid equality thresholds overlap')
            for row in thresholds:
                lower, upper = row['bracket']
                require(lower < row['gamma'] < upper and abs(row['peak']['fraction']-.5) < 5e-7,
                        'A fluid half-energy threshold missed equality')
                peak = row['peak']
                require(abs(peak['total_comoving_input_over_ai4']-peak['required_comoving_input_for_half_over_ai4']) < 1e-7,
                        'Comoving radiation input missed the half-energy threshold')
                if which == 'maximum_before_turnaround':
                    require(peak['instantaneous_Q'] == 0 and abs(peak['parent_equation_of_state']-1/3) < 2e-7,
                            'Full-expansion equality peak missed the source-free stationary condition')
                left, right = runs[str(lower)][which], runs[str(upper)][which]
                require((left['fraction']-.5)*(right['fraction']-.5) < 0,'Invalid fluid equality bracket')
        for run in runs.values():
            for stage in ('first_minimum','first_turnaround','first_half','maximum_before_turnaround'):
                row = run[stage]
                if row is None:
                    continue
                require(abs(row['r']/(row['scalar_density']+row['r'])-row['fraction']) < 2e-10,
                        'Fluid radiation fraction differs from the recorded energies')
                if stage=='first_half':
                    require(abs(row['fraction']-.5) < 2e-7, 'Fluid equality event failed')


def verify_feedback(feedback, gauge):
    for label, branch in feedback['branches'].items():
        reference = gauge['branches'][label]
        for step in (0., .1, .5):
            passive = branch['controls'][f'passive/{step:g}']['end']
            fixed_sum = reference['scan'][str(step)]['D_by_cutoff']['64']
            expect = reference['H0_over_Mpl']**2*fixed_sum/3
            require(math.isclose(passive['comoving_excess'],expect,rel_tol=2e-6,abs_tol=1e-18),
                    'Independent fixed-background gauge and passive discrete photon solvers disagree')
        for key, interval in branch['free_intervals'].items():
            start = branch['runs'][key]['stages']['switch_end']
            if start is None:
                require(interval['status'] == 'switch_not_reached', 'Unexpected coasting interval')
                continue
            require(abs(start['x']-interval['start']['x']) < 1e-9, 'Free radiation continuation missed switch exit')
            half_intervals = interval['half_intervals']
            crossings = interval['half_crossings']
            require(crossings == sorted(set(crossings)) and all(start['x'] < t < interval['stop']['x'] for t in crossings),
                    'Coasting equality roots are unordered or outside the freely propagating segment')
            edges = [t for row in half_intervals for t in (row['begin_x'],row['end_x'])]
            require(all(any(abs(t-edge) < 1e-9 for edge in edges) for t in crossings),
                    'Coasting equality root is missing from the dominance intervals')
            if key == '32/1.535' and label == 'phi_m=0.5':
                require(len(crossings) == 2, 'The supplementary free-radiation example missed a crossing')
            duration = math.fsum(r['end_x']-r['begin_x'] for r in half_intervals)
            efolds = math.fsum(math.log(r['end_a']/r['begin_a']) for r in half_intervals)
            require(abs(interval['total_efolds_at_or_above_half']-efolds) < 1e-9 and
                    abs(interval['total_proper_duration_at_or_above_half']-duration) < 1e-9,
                    'Radiation-dominance e-folds or duration differ from interval endpoints')
        envelope = branch['finite_band_envelope']
        variation = envelope['total_variation_log_I']
        ai = reference['points']['launch']['a']
        scale = reference['H0_over_Mpl']**2/(3*math.pi**2*ai**2)
        def photon_weight(limit):
            return sum(k*(k*k-1) for k in range(2,limit+1))
        cutoff = envelope['first_shell_cutoff_not_excluding_half']
        require(scale*photon_weight(cutoff-1)*math.sinh(variation)**2 < .5 <=
                scale*photon_weight(cutoff)*math.sinh(variation)**2,
                'Photon-band equality ceiling has the wrong integer threshold')
        for cap in envelope['rows']:
            n = math.floor(cap['initial_frequency_cap_over_Mpl']*ai/reference['H0_over_Mpl'])
            require(n == cap['maximum_shell'] and math.isclose(cap['fraction_ceiling'],
                    scale*photon_weight(n)*math.sinh(variation)**2, rel_tol=1e-12),
                    'Physical cutoff differs from the discrete frequency bound')
        for run in branch['runs'].values():
            require(run['sampled_peak_fraction'] <= run['fraction_ceiling_from_total_variation']+2e-9,
                    'Finite photon fraction exceeds its integrated squeezing bound')
            row = run['end']
            expected = row['photon_density']/(row['parent_density']+row['photon_density'])
            require(math.isclose(row['fraction_of_matter'],expected,rel_tol=1e-12,abs_tol=1e-21),
                    'Finite-system photon fraction differs from its densities')
            require(abs(row['a']**4*row['photon_density']-row['integrated_source_comoving']) < 4e-8,
                    'Finite-mode comoving source missed the energy change')


def verify_origin(origin, gauge):
    bulk = origin['infinite_tuned_leaf']
    require(bulk['vacuum_zero_mode_constant'] and not bulk['minimal_photon_normalizable']
            and bulk['continuum_lower_edge_mass_squared'] == 0,'Unexpected infinite-leaf photon spectrum')
    static = origin['static_cone']
    require(math.isclose(static['two_copy_constant_gauge_norm_times_g5_squared'],
                         2*static['cone_to_shell_proper_length'], rel_tol=1e-13),
            'Static two-sided Maxwell normalization differs from the proper interval')
    for key, branch in origin['branches'].items():
        response = branch['response']
        bound = math.cosh(2*abs(response['log_I_change']))
        require(math.isclose(bound,response['added_energy_over_adiabatic_upper_bound'],rel_tol=2e-12),
                'Localized gauge variation bound differs')
        spectrum = response['added_energy']
        require(spectrum[0]['L'] == 2 and all(1-1e-10 <= r['energy_over_adiabatic'] <= bound+2e-10 for r in spectrum),
                'A localized positive mode violates its finite-interval energy bound')
        time = gauge['branches'][key]['points']['first_minimum']['x']
        require(abs(branch['response']['crossing_x']-time) < 2e-8, 'Origin calculation used a different release')


def physics(products):
    require(products.keys() == PRODUCTS.keys(), 'Missing numerical product')
    gauge = products['gauge']
    verify_gauge(gauge,read(HERE/'prior/homogeneous/results.json')['branches'])
    verify_benchmarks(products['checks'])
    verify_fluid(products['fluid'],gauge)
    verify_feedback(products['feedback'],gauge)
    verify_origin(products['origin'],gauge)


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
            with TemporaryDirectory(prefix='hdblast-gauge-') as tmp:
                regenerate(Path(tmp),args.jobs)
