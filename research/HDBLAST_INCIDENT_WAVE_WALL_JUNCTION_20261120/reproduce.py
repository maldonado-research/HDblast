"""Reproduce the finite-energy incident waves and gravitating wall tests."""
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
BASE = "7fd4e3066499ed0e5264ec20a1bd97b634fa40ae"
PINS_SHA = "65427f3b6c3d6c70abcf91dcaf6abf3df1afc7a00705603207ae8d3e1a339981"
PRODUCTS = {name: f"{name}/{name.upper()}.json" for name in ("junction", "scattering", "horizon", "response")}
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
    require(pins['base'] == BASE and pins['previous_doi'] == '10.5281/zenodo.24385102',
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
    name, key = path[0], path[-1]
    if name=='junction':
        if 'max_errors' in path and 'radiation_exit' in path:
            # A volume-weighted residual can grow with a^3 or a^4. Its script
            # checks that scale directly. Compare small residuals by absolute size.
            floor = 2e-7 if key in ('comoving_first_law_absolute','radiation_comoving_transfer_absolute') else 1e-9
            return math.isclose(actual,expected,rel_tol=.01,abs_tol=floor)
        if key in ('max_friedmann_residual','max_force_free_scalar_error'):
            return 0<=expected<1e-7 and 0<=actual<1e-7
        if 'errors' in path or 'independent_expanding_solution_errors' in path:
            return 0 <= actual < 1e-8 and 0 <= expected < 1e-8
        # Turning points and vacuum densities are zero up to solver roundoff.
        return math.isclose(actual,expected,rel_tol=2e-7,abs_tol=2e-12)
    if name=='horizon':
        if 'proper_time_geodesic_controls' in path:
            return 0 <= actual < 5e-10 and 0 <= expected < 5e-10
        tolerance = 1e-12 if key in ('wave_metric_rr','measured_field_exponent','measured_tidal_exponent') else 0.
        return math.isclose(actual,expected,rel_tol=2e-7,abs_tol=tolerance)
    if name=='scattering' and 'branches' in path:
        if 'exact_semidiscrete_errors' in path:
            return 0 <= actual < 4e-5 and 0 <= expected < 4e-5
        if key=='energy_maximum_drift':
            return 0 <= actual < 1e-6 and 0 <= expected < 1e-6
        # Unit energy and unit wall width. Printed series have 12 decimal places.
        return math.isclose(actual,expected,rel_tol=2e-7,abs_tol=2e-10)
    if name=='response':
        if key in ('max_coordinate_junction_error','max_surface_energy_identity_error'):
            return 0 <= actual < 1e-8 and 0 <= expected < 1e-8
        if key=='far_field_amplitude':
            return 0 <= actual < 1e-12 and 0 <= expected < 1e-12
        if key in ('relative_grid_shift','relative_map_change','relative_time_change',
                   'peak_relative_balance_error','peak_absolute_moment_difference'):
            return math.isclose(actual,expected,rel_tol=2e-6,abs_tol=2e-11)
        if 'projected_energy'==key or 'moment_source_energy'==key:
            # The weak moving-wall projections exceed 1e-11. At a fixed wall
            # their analytic value is zero and the numerical result can underflow.
            return math.isclose(actual,expected,rel_tol=2e-7,abs_tol=1e-20)
        return math.isclose(actual,expected,rel_tol=2e-7,abs_tol=2e-14)
    return math.isclose(actual,expected,rel_tol=2e-7,abs_tol=0.)


def compare(actual, expected, path=(), differences=None):
    label = '/'.join(path)
    if isinstance(expected, dict):
        require(isinstance(actual,dict) and actual.keys()==expected.keys(), f'Different keys at {label}')
        for key in expected:
            compare(actual[key],expected[key],(*path,key),differences)
    elif isinstance(expected, list):
        require(isinstance(actual,list) and len(actual)==len(expected), f'Different size at {label}')
        for i, (a,e) in enumerate(zip(actual,expected)):
            compare(a,e,(*path,str(i)),differences)
    elif isinstance(expected, float):
        good = (type(actual) in (int,float) and math.isfinite(actual) and math.isfinite(expected)
                and close_number(actual, expected, path))
        if not good:
            message = f'{label}: expected {expected!r}, obtained {actual!r}'
            if differences is None:
                raise ValueError(message)
            differences.append(message)
    else:
        require(type(actual) is type(expected) and actual==expected, f'Scalar differs at {label}')


def equal(a, b, label, *, rel=5e-12, absolute=0.):
    require(math.isclose(a, b, rel_tol=rel, abs_tol=absolute), f'{label}: {a} differs from {b}')


def check_radiation_exit(exit_data):
    runs = exit_data['runs']
    require(runs['baseline']['e_folds_radiation']>0 and
            runs['zero_kick']['e_folds_radiation']==runs['zero_damping']['e_folds_radiation']==0,
            'Source or damping controls failed')
    for label,row in runs.items():
        intervals = row.get('expanding_radiation_intervals',[])
        total = 0.
        for interval in intervals:
            start,end=interval['start'],interval['end']
            require(end['t']>start['t'] and end['a']>start['a'] and
                    min(start['radiation_fraction'],end['radiation_fraction']) >= .9-1e-8,
                    'Inconsistent radiation interval '+label)
            amount = math.log(end['a']/start['a'])
            equal(interval['e_folds'],amount,label)
            total+=amount
        equal(row['e_folds_radiation'],total,label)
        for events in row.get('events',{}).values():
            for snap in events:
                a,e,h=snap['a'],snap['e'],snap['H']
                equal(h*h+1/a**2,e*(2+e)+snap['mu']/a**4,label,absolute=1e-7)
                equal(snap['radiation_fraction'],snap['radiation']/e,label)
    threshold = exit_data['plateau_exit_threshold']
    lo,hi = threshold['bracket']
    require(lo<threshold['kick']<hi, 'Unbracketed plateau exit boundary')


def check_scalar_capture(capture):
    primary = capture['cases']['primary']
    require(primary['pulse_end']['H']<0 and primary['terminal']['H']>0, 'Primary scalar wall failed to reverse')
    for control in ('same_potential_no_pulse','same_pulse_zero_potential'):
        require(not capture['cases'][control]['turn'], 'A scalar control reproduced the bounce')
    for label, row in capture['cases'].items():
        for snapshot in ('pulse_end','plateau_entry','bounce','terminal'):
            if snapshot not in row:
                continue
            snap = row[snapshot]
            a,h,mu,e,p = (snap[k] for k in ('a','H','mu','e','p'))
            equal(h*h+1/a**2,e*(2+e)+mu/a**4,label,absolute=1e-8)
            equal(e+p,snap['pi']**2,label,absolute=1e-12)
        energy, geom = row['energy'], row['domain']
        equal(energy['two_side_wall_deposition'],energy['terminal_internal_energy']+
              energy['terminal_pressure_work'],label,absolute=1e-8)
        equal(energy['mu_loss_each_side'],4-row['terminal']['mu'],label,absolute=1e-12)
        require(geom['min_mu']>0 and geom['min_f_at_source']>0 and geom['min_source_power']>=-1e-12,
                'Negative mass, flux or emitter lapse '+label)
        if row['turn']:
            turn = row['bounce']
            integral = row['plateau_first_integral']
            equal(turn['H'],0.,label,absolute=1e-10)
            equal(turn['a'],integral['outer_root'],label,absolute=1e-8)
            require(turn['t']>row['plateau_entry']['t']>=capture['duration'] and
                    turn['a']>capture['source_radius'] and integral['acceleration_at_bounce']>0,
                    'Invalid source or turning point '+label)


def check_junction(data):
    capture = data['scalar_capture']
    check_scalar_capture(capture)
    primary = capture['cases']['primary']
    validation = data['scalar_validation']
    equal(validation['bounce_radius'],primary['bounce']['a'],'independent scalar minimum',absolute=5e-9)
    equal(validation['bounce_time'],primary['bounce']['t'],'independent bounce time',absolute=5e-9)
    equal(validation['wall_deposited_energy'],primary['energy']['two_side_wall_deposition'],
          'independent energy deposited',absolute=5e-9)
    check_radiation_exit(data['radiation_exit'])
    prescribed = data['prescribed_source']
    prescribed_exit = prescribed['baseline']['radiation_exit']
    require(prescribed['no_pulse']['bounce'] is None and prescribed['no_potential']['bounce'] is None,
            'An independent mass-profile control bounced')
    for label,row in prescribed.items():
        equal(row['final_bulk_mass'],4*(1-row['mass_loss_fraction']), 'prescribed mass '+label)
        require(row['deposited_proper_energy'] <= row['proper_emitter_energy']+1e-12,
                'Deposited energy exceeds the emitter budget '+label)
    value = sum(math.log(w['end']['a']/w['start']['a']) for w in prescribed_exit['windows'])
    equal(value, prescribed_exit['e_folds'], 'prescribed pulse expansion interval')
    require(value>0 and bool(prescribed_exit['events']['maximum']),'Missing prescribed-source recollapse')
    for tag, case in data['cases'].items():
        start, end = case['initial'], case['end']
        for snap in (start, *case['turns'], end):
            radius, hubble, density, mass = (snap[k] for k in ('a','H','e','mass'))
            equal(hubble**2 + 1/radius**2, density*(2+density)+mass/radius**4, tag)
            equal(snap['B'], 1+density+hubble, tag)
            equal(snap['comoving_radiation'], density*radius**4, tag, absolute=2e-14)
        for event in case['turns']:
            require(abs(event['H']) < 1e-11, f'Nonzero H at turnaround {tag}')
        require(case['minimum_wall_clearance_from_source'] > 0 and
                case['minimum_source_f'] > 0 and case['minimum_sampled_matter'] > -1e-10,
                f'Invalid exterior or matter density {tag}')
        energy = case['proper_energy']
        equal(energy['wall_proper_heat'], energy['wall_internal_energy_change'] +
              energy['wall_expansion_work'], tag, absolute=2e-10)
        transfer = end['comoving_radiation'] - start['comoving_radiation']
        equal(transfer, case['independent']['deposited_comoving'], tag, absolute=1e-8)
        loss = start['mass']-end['mass']
        if loss:
            equal(2*transfer/loss, case['normalized_comoving_transfer'], tag)
            equal(energy['wall_proper_heat']/energy['emitter_proper_energy'], energy['heat_over_source'], tag)
        else:
            require(abs(transfer)<1e-10 and case['normalized_comoving_transfer'] is None, tag)
    for tag, row in data['asymmetric'].items():
        h,a,e = row['end_H'],row['end_radius'],row['end_e']
        left,right = row['end_radial_curvatures']
        equal((left+right)/2, 1+e, tag, absolute=1e-9)
        equal(row['end_harmonic_mean'], 2/(1/left+1/right), tag)
        mass1,mass2 = row['end_masses']
        contrast = (mass1-mass2)/(4*a**4*(1+e))
        equal(row['end_asymmetry_squared'], contrast**2, tag, absolute=2e-16)
        equal(h*h, e*(2+e)-1/a**2+(mass1+mass2)/(2*a**4)+contrast**2,
              tag, absolute=1e-9)
        require(row['maximum_sampled_acceleration_over_radius']<=0, 'Positive acceleration '+tag)
    symmetric = data['cases']['expanding/0.5/0.2']['end']
    two_sided = data['asymmetric']['two_sided_half']
    for k in ('H','e'):
        equal(symmetric[k],two_sided['end_'+k],'symmetric reduction',absolute=2e-10)
    pressure = data['normal_pressure']
    for geometry in ('double_null','vaidya'):
        for row in pressure[geometry]:
            equal(row['hessian_UU_over_radius'], row['predicted_hessian_over_radius'],
                  'spherical Hessian law', absolute=2e-12)
            if geometry=='vaidya':
                equal(row['kappa2_Tnn'], row['null_dust_kappa2_Tnn'], 'Vaidya pressure', absolute=2e-12)
    for stream in pressure['counterstreaming']:
        equal(stream['Tnn'], stream['incident']+stream['reflected'], 'counterstream pressure')
        equal(stream['net_flux'], stream['incident']-stream['reflected'], 'counterstream flux')
    unequal = pressure['asymmetric_pressure']
    radial = unequal['radial_curvatures']
    harmonic = 2/sum(1/x for x in radial)
    weighted = sum((m+t)/x for m,t,x in zip(unequal['mass_over_radius_fourth'],
                   unequal['kappa2_Tnn_over_three'],radial,strict=True))*harmonic/2
    equal(weighted,unequal['Qbar'],'weighted normal pressure')
    exact = unequal['exact_critical_pressure']
    equal((1-2*unequal['density']-3*exact)*harmonic,1+weighted,'asymmetric pressure threshold')
    require(exact <= unequal['necessary_ceiling'], 'Asymmetric pressure exceeds ceiling')
    threshold = pressure['pressure_threshold']
    vacuum, radiation = pressure['vacuum_bounce'], pressure['radiation_at_rest']
    require(vacuum['pressure']<threshold['critical_pressure']<radiation['pressure'],
            'Vacuum and radiation should straddle the rest-point threshold')
    equal(vacuum['acceleration_over_radius'], .25, 'vacuum bounce')
    require(radiation['acceleration_over_radius'] < 0, 'Radiation acceleration at rest')
    for sample in vacuum['samples']:
        equal(sample['radius'], 2*math.cosh(sample['proper_time']/2.), 'vacuum orbit')
    impulse = data['impulse']
    widths = sorted((float(k),v) for k,v in impulse['smooth'].items())
    error = [abs(row['critical_mass_fraction']-impulse['critical_fraction_limit']) for _,row in widths]
    require(error == sorted(error) and error[0] < .005*impulse['critical_fraction_limit'],
            'Smooth threshold failed to approach the impulse result')


def check_horizon(data):
    for tail in data['wave_tails'].values():
        p = tail['field_amplitude_power']
        equal(tail['measured_tidal_exponent'], 2*p-4, 'tidal slope', absolute=1e-5)
        equal(tail['measured_field_exponent'], p-3, 'Maxwell slope', absolute=1e-5)
        for row in tail['by_depth']:
            equal(row['signal_ricci'], row['weighted_maxwell_screen_squared'], 'Einstein-Maxwell projection')
    for row in data['proper_time_geodesic_controls'].values():
        require(max(row.values())<5e-10, 'Conserved momentum or proper-time geodesic failed')
    for chart in data['chart_controls'].values():
        for row in chart:
            equal(row['det_ratio'], 1., 'wave determinant', absolute=2e-12)
    for row in data['transverse_localization']:
        equal(row['rms_radius']*row['normalized_residual'], row['analytic_product'], 'transverse norm')
        require(row['product']>=1-5e-12, 'Transverse localization bound failed')
    sandwich = data['sandwich']
    equal(sandwich['exit_proper_time']+sandwich['proper_time_remaining_at_exit'],
          data['horizon_proper_time'], 'proper time to horizon')
    require(all(x == 0 for x in sandwich['wave_at_later_depths']), 'Compact wave failed to switch off')


def check_reflection(reflection):
    # A direct Fourier time series checks the sign, phase convention and flux
    # Jacobian in the spectral derivative calculation.
    import importlib.util
    import numpy as np
    spec = importlib.util.spec_from_file_location('reflection_check', HERE/'scattering/phase.py')
    solver = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(solver)
    packet = next(r for r in reflection['packets'] if r['center_mass']==1.)
    center, width, k = packet['center_mass'], packet['half_width'], reflection['k']
    low,high = np.hypot(k,[center-width,center+width])
    bins = 96
    spacing = (high-low)/bins
    omega = low+(np.arange(bins)+.5)*spacing
    mass = np.sqrt(omega**2-k**2)
    amp = np.sqrt(omega/mass)*np.exp(-.5/(1-((mass-center)/width)**2))
    times = np.fft.fftfreq(16*bins,d=spacing/(2*np.pi))
    centroids = []
    for selected in (False, True):
        wave = solver.jost_reflection(mass,8.,selected)
        require(np.max(np.abs(np.abs(wave)**2-1))<1e-12, 'Incident and reflected flux differ')
        flux = abs(np.fft.fft(amp*wave, len(times)))**2
        centroids.append(float(flux@times/flux.sum()))
    equal(centroids[1]-centroids[0], packet['mean_relative_delay'], 'direct arrival time', absolute=2e-6)
    for row in reflection['rows']:
        equal(row['relative_delay'], row['independent_relative_delay'], 'independent group delay', absolute=2e-6)
        delta = math.remainder(row['relative_phase_radians']-row['independent_phase_radians'],2*math.pi)
        require(abs(delta)<3e-7, 'Independent phase extraction differs')


def check_wave(data):
    from scipy.integrate import simpson
    for law, branch in data['branches'].items():
        for kind, row in branch.items():
            spec, obs = row['spec'], row['checks']
            equal(spec['pole_energy'], (spec['moment_v']**2+row['spec']['k']**2*spec['moment_q']**2)
                  /(2*spec['half_norm']), 'incident massless energy', absolute=1e-12)
            if kind == 'balanced':
                require(spec['pole_energy']<1e-20, f'Unbalanced incident seed {law}')
            require(obs['energy_maximum_drift'] < 1e-6 and obs['initial_outward_momentum']<0,
                    f'Incident wave or energy drift wrong for {law} {kind}')
            impulse = float(simpson(row['wall_pressure'],x=row['time']))
            equal(impulse, obs['integrated_wall_pressure'], 'wall impulse')
            for q,v,p in zip(row['wall_field'],row['wall_velocity'],row['wall_pressure'],strict=True):
                equal(p, .5*(v*v-data['k']**2*q*q), 'wall normal stress', absolute=1e-11)
            change = obs['final_outward_momentum']-obs['initial_outward_momentum']
            mismatch = change-impulse-obs['integrated_background_force']+obs['integrated_right_pressure']
            equal(mismatch/obs['initial_energy'], obs['momentum_balance_error'], 'radial momentum', absolute=5e-12)
            require(abs(mismatch/obs['initial_energy']) < .008, 'Radial momentum failed to balance')
    for law, bound in data['support_bounds'].items():
        tails = bound['exterior_norm_fraction']
        equal(tails['6']-tails['10'], bound['packet_interval_fraction'], 'supported norm fraction')
        equal(bound['one_copy_norm'], data['branches'][law]['raw']['spec']['half_norm'],
              'normalization of fixed mode', absolute=1e-12)
        for kind, row in data['branches'][law].items():
            require(row['spec']['pole_energy'] <= bound['packet_interval_fraction']+1e-12,
                    f'The support bound failed for {law} {kind}')
    ancestor = read(HERE/'prior/BOUNDS.json')['leaves']['2']['0.0']['single_copy_norm']
    equal(ancestor, data['support_bounds']['wall']['one_copy_norm'], 'archived central normalization')
    check_reflection(data['reflection'])


def check_response(data, waves):
    for law, branch in data['branches'].items():
        for row in (*branch['runs'].values(), *branch['controls'].values()):
            final, energy = row['final'], row['incident']
            states = row.get('times', [final])
            for state in states:
                error = abs(state['field_energy']+state['wall_mechanical_energy']+
                            state['dissipated_energy']-energy)/energy
                require(error <= row['peak_relative_balance_error'] + 3e-12, 'Inconsistent energy accounting')
                require(0 <= state['projected_energy'] <= state['field_energy'] and
                        state['dissipated_energy'] >= 0, 'Passive wall has invalid energy')
            require(row['peak_relative_balance_error'] < 3e-4 and row['far_field_amplitude'] < 1e-12,
                    'A wall run lost energy or reached the far boundary')
            equal(final['wall_mechanical_energy'],
                  .5*(final['wall_speed']**2+row['frequency']**2*final['wall_position']**2),
                  'wall mechanical energy', absolute=1e-24)
            if row['fixed']:
                require(final['dissipated_energy']==0 and row['peak_absolute_displacement']==0,
                        'Fixed wall absorbed work')
            if row['damping']==0:
                require(final['dissipated_energy']==0, 'Undamped wall dissipated energy')
            require(row['peak_wall_speed']<1 and waves['branches'][law]['balanced']['spec']['k']==1.,
                    'Wall speed or incident mode differs')
        fixed = branch['runs']['0.025/fixed']
        weak = branch['controls']['weak']
        for quantity,coefficient,power,ratio in (
                ('dissipated_energy','linear_dissipated_coefficient',2,'weak_dissipation_over_leading'),
                ('projected_energy','linear_projected_coefficient',3,'weak_projection_over_leading')):
            predicted = fixed[coefficient]*weak['incident']**power
            equal(weak['final'][quantity]/predicted, branch['tests'][ratio], 'predicted scaling')
        require(all(abs(v-1)<.025 for v in branch['tests'].values()), 'Weak-field limit failed')
    compatibility = data['compatibility']
    for row in compatibility['examples']:
        for key in ('numerical_positive_integral','numerical_negative_magnitude'):
            equal(row[key],row['exact_each_integrated_magnitude'],'integrated surface enthalpy',absolute=1e-11)
        equal(row['scaled_affine_null_energy'],row['scaled_affine_null_energy_by_parts'],
              'averaged null energy on the wall', absolute=1e-11)
        require((row['scaled_affine_null_energy'] < 0) == (row['displacement'] != 0),
                'Wrong averaged null-energy sign')
        require(row['peak_normal_scalar_ratio']>=1, 'Invalid scalar junction factor')


def physics(products):
    require(products.keys() == PRODUCTS.keys(), 'Missing numerical product')
    check_junction(products['junction'])
    check_horizon(products['horizon'])
    check_wave(products['scattering'])
    check_response(products['response'], products['scattering'])


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
    physics(data)
    differences = []
    for name,path in PRODUCTS.items():
        compare(data[name],read(HERE/path),(name,),differences)
    require(not differences,'Numerical replay differences:\n'+'\n'.join(differences))
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
            with TemporaryDirectory(prefix='hdblast-incident-') as tmp:
                regenerate(Path(tmp),args.jobs)
