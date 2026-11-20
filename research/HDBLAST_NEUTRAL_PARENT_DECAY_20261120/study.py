"""Verify or replay the closed-wall neutral-parent decay experiment."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
PRIOR_DIGEST = '5dbed28e158584403a4cd238d392aab15c00480522cf489d1ee0a36f1d550f90'
PRODUCTS = ('coupled/DECAY', 'coupled/COHORTS', 'survival/SURVIVAL', 'spectrum/SPECTRUM',
            'consistency/CONSISTENCY', 'expansion/FATE', 'kinetic/KINETIC')


def require(test, message):
    if not test:
        raise ValueError(message)


def verify_files():
    manifest = ROOT/'SHA256SUMS'
    expected = {}
    for line in manifest.read_text().splitlines():
        digest, name = line.split('  ', 1)
        path = Path(name)
        require(not path.is_absolute() and path.as_posix()==name and '..' not in path.parts
                and name not in expected and not (ROOT/path).is_symlink(), 'Invalid manifest entry '+name)
        expected[name] = digest
    present = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file()
               and p!=manifest and '__pycache__' not in p.parts}
    require(set(expected)==present, 'Unexpected file inventory '+str(set(expected)^present))
    for name, digest in expected.items():
        with (ROOT/name).open('rb') as f:
            require(hashlib.file_digest(f,'sha256').hexdigest()==digest, 'Checksum mismatch '+name)
    require(expected.get('prior/CHARGED.json')==PRIOR_DIGEST, 'Historical input differs from #187')
    print(f'Verified {len(expected)} files including the pinned previous study.', flush=True)


def compare(value, expected, path):
    if isinstance(expected, dict):
        require(isinstance(value,dict) and value.keys()==expected.keys(), 'Changed keys at '+path)
        for key, item in expected.items():
            compare(value[key], item, path+'/'+key)
    elif isinstance(expected, list):
        require(isinstance(value,list) and len(value)==len(expected), 'Changed length at '+path)
        for i, (actual, reference) in enumerate(zip(value,expected)):
            compare(actual, reference, path+'/'+str(i))
    elif type(expected) is float:
        # The collisionless fluid's signed drift at maximum internal share is small
        # and sensitive to optimization at the maximum. The heating control is independent.
        absolute = 2e-9 if path=='coupled/DECAY/runs/no_explicit_drag/peak_daughter_velocity' else 2e-10
        require(type(value) in (int,float) and math.isfinite(value)
                and math.isclose(value,expected,rel_tol=5e-7,abs_tol=absolute),
                f'{path}: found {value}, expected {expected}')
    else:
        require(type(value) is type(expected) and value==expected, 'Changed value at '+path)


def physical_checks(products):
    c = products['coupled/DECAY']
    runs = c['runs']
    p = runs['primary']
    require(len(runs)==len(c['parameters']) and 0<p['bounce_tau']<p['exit_complete_tau']<p['threshold_tau']
            <p['half_decayed_tau']<p['ninety_decayed_tau']<p['maximum_tau'] and p['end_H']<0,
            'Expansion, threshold and birth events out of order')
    require(all(p[key]>0 for key in ('min_lapse','min_radial_gradient_norm','min_bulk_mass'))
            and p['friedmann_error']<1e-8 and p['gas_balance_error']<1e-8
            and p['parent_balance_error']<1e-13 and p['particle_number_error']<1e-10,
            'Constraint, energy or particle budget failed')
    require(.91 < p['maximum_intrinsic_fraction'] < .92 and .007 < p['ninety_intrinsic_expanding_efolds'] < .0085
            and p['end_radiation']==0 and p['peak_daughter_work_a3']>100000*p['peak_injected_kinetic_a3'],
            'Charged matter heating or its source changed')
    for key in ('subthreshold_parent','scalar_mass_only','threshold_missed','vacuum'):
        require(runs[key]['end_daughter_fraction']==0, 'Unexpected daughters in '+key)
    require(runs['no_charge']['maximum_intrinsic_fraction']<1e-5 and
            .04<runs['no_explicit_drag']['maximum_intrinsic_fraction']<.05 and
            all(runs[n]['maximum_intrinsic_fraction']<.72 for n in ('dose_reduced','dose_increased')),
            'Charge, drag or dose comparison failed')
    require(abs(runs['medium']['maximum_intrinsic_fraction']-p['maximum_intrinsic_fraction'])<1e-5,
            'Radial refinement changed peak heating')
    rows = [dict(zip(c['trace_columns'],row,strict=True)) for row in p['trace']]
    require(all(-1e-12<=row['f']<=1+1e-12 and row['parent']>=-1e-14 and row['rad']==0 for row in rows)
            and max(abs(row['f']+math.expm1(-row['optical_depth'])) for row in rows)<1e-7,
            'Interpolated trace violates decay or particle count')
    coh = products['coupled/COHORTS']
    snapshot = next(r for r in coh['comparison'] if r['t']==12)
    require(abs(snapshot['exact_packets_on_given_field']['mean_current_per_charge_density']
                -snapshot['fluid_current_per_charge_density'])>.7
            and coh['largest_current_change_halving_time_resolution']<.01,
            'Independent current reconstruction or refinement failed')
    s = products['survival/SURVIVAL']
    require(s['history_sha256']==PRIOR_DIGEST and s['ordinary_charge_normalization_controls'][1]['weak_yukawa_fraction_upper']<4e-10
            and abs(s['finite_time']['cases']['128.0']['rolling']['ratio_to_local']-1)<.01
            and abs(s['primary']['100']['first_arrivals']['0.5']-p['half_decayed_tau'])<.01,
            'Neutral-background survival and finite-time controls failed')
    sp = products['spectrum/SPECTRUM']
    mag = sp['magnetic_history']
    require(mag['maximal_allowed_q1_momentum']<mag['bounded_q1_momentum']/1e7
            and sp['checks']['matrix_absolute_error']<1e-12
            and abs(sp['checks']['many_shell_flat_test']['compact_over_flat_FGR']-1)<1e-5,
            'Closed spinor spectrum or continuum normalization failed')
    edge = sp['field_driven_edge']
    require(edge['checks']['static_transition']<1e-20 and
            edge['cases']['kinematic_floor']['long']['finite_instantaneous_probability']>.999,
            'Single-mode electric transition controls failed')
    con = products['consistency/CONSISTENCY']
    require(math.isclose(con['evaluation_point']['theta'],p['peak_daughter_theta'])
            and con['rows']['visible_coupling_mu100']['thermal_swave_drag_ceiling_times_L']<.04
            and abs(con['calibrated_additional_elastic_kernel']['fitted_nu_times_L']-6)<1e-7,
            'Microscopic rate checks failed')
    kinetic = products['kinetic/KINETIC']
    collisionless = kinetic['runs']['reference']
    require(kinetic['zero_charge_analytic_field_error']<1e-8 and
            abs(kinetic['runs']['zero']['turn_net_field_loss_fraction'])<1e-8 and
            .90<collisionless['turn_net_field_loss_fraction']<.92 and
            abs(collisionless['turn_current_work_over_initial_field']-collisionless['turn_net_field_loss_fraction'])<1e-7,
            'Collisionless Maxwell energy and zero-charge controls failed')
    require(collisionless['returned_current_work_over_initial_field']>.5
            and collisionless['terminal_spectrum']['parallel_share_of_pressure_in_energy_frame']>.99
            and collisionless['max_normalized_total_balance']<2e-5
            and abs(collisionless['turn_field_share_of_added']-kinetic['runs']['more_births']['turn_field_share_of_added'])<1e-4,
            'Kinetic pressure distribution, returned work or refinement changed')
    tail = products['expansion/FATE']['scalar_threshold']
    require(tail['intrinsic_share']>=.9 and abs(tail['friedmann_h2_residual'])<1e-7
            and .25<tail['dust']['critical_vacuum_share']<.26
            and .41<tail['radiation']['critical_vacuum_share']<.42
            and 1.98<tail['radiation']['radius_targets_without_vacuum']['double']['factor_over_measured_anchor']<2.01,
            'Late-time expansion bounds failed')


def run_all(out):
    out = out.resolve()
    require(not out.is_relative_to(ROOT) and (not out.exists() or (out.is_dir() and not any(out.iterdir()))),
            'Use a new or empty output directory outside the study')
    out.mkdir(parents=True, exist_ok=True)
    targets = {name: out/(name+'.json') for name in PRODUCTS}
    for path in targets.values():
        path.parent.mkdir(exist_ok=True)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
    def run(script, name, *args):
        print('Recomputing '+name, flush=True)
        subprocess.run([sys.executable,'-I','-B',str(ROOT/script),
                        '--output',str(targets[name]),*map(str,args)], check=True, env=env)
    with ThreadPoolExecutor(max_workers=3) as pool:
        independent = [
            pool.submit(run,'survival/local_decay.py','survival/SURVIVAL'),
            pool.submit(run,'spectrum/spectrum.py','spectrum/SPECTRUM'),
            pool.submit(run,'coupled/cohorts.py','coupled/COHORTS')]
        run('coupled/decay.py','coupled/DECAY','--jobs','2')
        primary = targets['coupled/DECAY']
        dependent = [
            pool.submit(run,'expansion/fate.py','expansion/FATE','--decay',primary),
            pool.submit(run,'consistency/consistency.py','consistency/CONSISTENCY','--dynamics',primary),
            pool.submit(run,'kinetic/vlasov.py','kinetic/KINETIC','--source',primary)]
        for future in independent+dependent:
            future.result()
    products = {name:json.loads(path.read_bytes()) for name,path in targets.items()}
    physical_checks(products)
    for name, result in products.items():
        compare(result, json.loads((ROOT/(name+'.json')).read_bytes()), name)
    print('All seven numerical products, the historical input and physical controls passed.', flush=True)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    parser.add_argument('--output-dir', type=Path)
    args = parser.parse_args()
    if args.output_dir and not args.verify:
        parser.error('--output-dir requires --verify')
    verify_files()
    if args.verify:
        if args.output_dir:
            run_all(args.output_dir)
        else:
            with tempfile.TemporaryDirectory(prefix='hdblast-decay-replay-') as tmp:
                run_all(Path(tmp))
