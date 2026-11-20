"""Reproduce the curved Maxwell and scalar wall calculations."""
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
BASE = "315e1edf085d601a93a616947db75d0929f9f8fe"
PINS_SHA = "3bb4656670beac57ab8932413f8fb1afb7c636a6a506b5c7904bf1d731bbbdde"
PRODUCTS = {name: f"{name}/{name.upper()}.json" for name in ("fixed", "identities", "gravity", "loading")}
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
    require(pins['base'] == BASE and pins['previous_doi'] == '10.5281/zenodo.24386416',
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
    if name=='loading':
        if 'forward' in path and key=='tube_mass_error':
            # This is the peak difference of two independent mass integrals.
            # The primary runs also satisfy the absolute balance bound in FORWARD.md.
            if path[-2] in ('fluid_primary','maxwell_primary') and max(actual,expected)>=1.8e-6:
                return False
            return min(actual,expected)>=0 and math.isclose(actual,expected,rel_tol=2e-7,abs_tol=1e-10)
        if 'errors' in path:
            # The respective scripts enforce dimensional bounds on conservation
            # residuals. Their tiny floating-point values need not coincide.
            return 0<=actual<2e-6 and 0<=expected<2e-6
        if 'refinement_errors' in path:
            return 0<=actual<5e-6 and 0<=expected<5e-6
        if key in ('max_friedmann_error','friedmann_error','fluid_balance_error',
                   'surface_balance_error','sheet_balance_error'):
            return 0<=actual<2e-6 and 0<=expected<2e-6
        if 'minima' in path or key in ('minimum_mass','post_exit_min_velocity_during_expansion'):
            # These are extrema sampled on adaptively selected ODE steps.
            return math.isclose(actual,expected,rel_tol=2e-5,abs_tol=1e-7)
        if key=='net_force_per_step_slope':
            return math.isclose(actual,expected,rel_tol=3e-6,abs_tol=2e-8)
        if key in ('I','O','X','Y','Bmag','Kprime','pi','flux'):
            return math.isclose(actual,expected,rel_tol=2e-7,abs_tol=2e-9)
    if name=='gravity':
        # Quantities use unit AdS radius. Roundoff in constraint residuals
        # can exceed their central value; grid convergence is checked separately.
        return math.isclose(actual,expected,rel_tol=1e-6,abs_tol=3e-10)
    if name=='fixed':
        if key=='refined_mesh_error':
            return 0 <= expected < 1e-8 and 0 <= actual < 1e-8
        if key=='coarse_mesh_error':
            return math.isclose(actual,expected,rel_tol=.002,abs_tol=1e-10)
        return math.isclose(actual,expected,rel_tol=2e-7,abs_tol=2e-10)
    if name=='identities' and 'collar' in path:
        if key in ('max_null_constraint','normalized_constraint','max_mass_flux_residual'):
            return math.isclose(actual,expected,rel_tol=1e-4,abs_tol=2e-10)
        if key=='vacuum_enclosed_mass_error':
            return 0<=actual<1e-7 and 0<=expected<1e-7
    if name=='identities' and 'max_absolute_residuals' in path:
        return 0 <= expected < 1e-7 and 0 <= actual < 1e-7
    return math.isclose(actual,expected,rel_tol=2e-7,abs_tol=2e-12)


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


def check_fixed(data):
    for label,row in data['cases'].items():
        equal(row['source_radius']+row['gap'],row['radius'],label)
        for stage in ('end_of_port_drive','direct_front_after_end','readout'):
            s = row[stage]
            incident = s['incident_left']+s['incident_right']
            outgoing = s['outgoing_left']+s['outgoing_right']
            retained = sum(s[k] for k in ('scalar','surface_electric','surface_magnetic',
                                          'bulk_electric','bulk_radial_gradient','bulk_angular'))
            equal(incident-outgoing-retained,s['balance'],label,absolute=6e-12)
            require(abs(s['balance'])<1e-7 and retained>=-1e-12, 'Invalid normalized wave energy '+label)
            equal(s['scalar_velocity']**2/2,s['scalar'],label,absolute=2e-12)
            equal(s['K']*s['wall_velocity']**2/2,s['surface_electric'],label,absolute=1e-11)
            equal(2*s['K']*(s['wall_q']/row['radius'])**2,s['surface_magnetic'],label,absolute=1e-11)
            equal(s['electric_impulse']-s['magnetic_impulse'],s['scalar_velocity'],label,absolute=2e-9)
        final=row['readout']
        equal(final['incident_left']+final['incident_right'],0. if row['channel']=='none' else 1.,
              'incident dose '+label,absolute=1e-8)
        if row['channel']=='one':
            require(final['scalar']<=.5+1e-8, 'Single-sided pulse violates parity bound')
        if row['channel'] in ('odd','none') or row['g']==0:
            require(final['scalar']<1e-12, 'Control acquired scalar energy '+label)
    far = data['cases']['R4096']['readout']['scalar']
    equal(abs(far-data['flat_limit']['scalar']),data['scalar_difference_at_R4096'],'flat limit',
          absolute=2e-12)
    helm = data['helmholtz']
    for point in helm['points']:
        equal(point['one_side_reflectance']+point['one_side_transmittance'],1.,'two-port unitarity')
        require(point['refined_mesh_error']<1e-8,'Wrong Bessel amplitude or wall condition')
    for packet in helm['pulse_transmission']:
        direct = packet['outgoing_right']/packet['incident']
        remaining = packet['leftover_energy']/packet['incident']
        require(remaining < 3e-7 and direct-3e-7<=packet['spectral_fraction']<=direct+remaining+3e-7,
                'Time-domain transmitted energy differs from the spectrum')
        equal(packet['spectral_fraction'],packet['shorter_transform_fraction'],
              'incident Fourier resolution',absolute=1e-7)


def check_identities(data):
    require(data['boosted_checks']>=data['sampled_radii']*data['sampled_times'] and
            data['comoving_sheet_kinetic_coefficient_over_radius_required']<0,
            'Reference solution or static matching condition failed')
    collar = data['collar']
    finest = collar['launch_fine']
    require(finest['minimum_mass_on_strip']>0 and finest['minimum_F_on_strip']>0 and
            finest['normalized_constraint']<.001, 'Invalid reconstructed near-wall strip')
    for row in collar['levels']['4']['late_samples'].values():
        radius = row['r']
        equal(radius**2*(1+radius**2-row['F']), row['M'], 'spherical mass on reconstructed tube',
              absolute=2e-11)
    require(collar['vacuum_enclosed_mass_error'] < 1e-7, 'Vacuum mass changed')
    for name,residual in data['max_absolute_residuals'].items():
        require(0<=residual<1e-7,'Geometric residual outside stated numerical scale '+name)


def check_boundary(boundary):
    runs = boundary['runs']
    primary = runs['primary']
    for label,row in runs.items():
        states = [row['final'], *row.get('trajectory',[])]
        for events in row['events'].values():
            states.extend(events)
        for st in states:
            a,h,mu,e,p = (st[k] for k in ('a','H','mu','e','p'))
            equal(h*h+1/a**2, e*(2+e)+mu/a**4,'wall Friedmann '+label,absolute=2e-6)
            equal(st['I']+st['O'],st['X'],'characteristic velocity',absolute=1e-12)
            equal(st['O']-st['I'],st['Y'],'characteristic derivative',absolute=1e-12)
            equal(st['Bmag'],2*st['A']/a**2,'magnetic curvature term',absolute=1e-12)
            equal(st['sheet'],.5*st['K']*(st['X']**2+st['Bmag']**2),'Maxwell sheet energy',absolute=1e-12)
            # The scalar carries e+p=pi². The sheet and the material radiation
            # both have isotropic pressure equal to one third of their energy.
            equal(e+p,st['pi']**2+4*(st['sheet']+st['rad'])/3,'surface enthalpy',absolute=1e-12)
            normal=st['I']**2+st['O']**2+2*st['odd_squared']-.5*st['Bmag']**2
            equal(normal,st['Tnn'],'independent normal pressure',absolute=1e-13)
            equal(st['flux'],st['I']**2-st['O']**2,'wall flux',absolute=1e-13)
            equal(st['acc_over_a'],(1+e)*(1-2*e-3*p)-1-mu/a**4-2*normal,
                  'radial acceleration',absolute=1e-11)
        equal(row['incident_dose'],row['config']['dose'],'supplied dose '+label,absolute=1e-8)
        for epoch in row['epochs']:
            require(epoch['end']>epoch['start'] and epoch['end_radius']>epoch['start_radius'],
                    'Empty or contracting epoch '+label)
            equal(epoch['efolds'],math.log(epoch['end_radius']/epoch['start_radius']),'expansion '+label)
    require(primary['epochs'] and all(not runs[k]['epochs'] for k in
            ('no_pulse','zero_coupling','no_potential','no_damping')),'Control classification failed')
    jump = boundary['thin_limit']['jump']
    equal(jump['incident']-jump['returned'],jump['e'],'thin limit energy',absolute=2e-8)
    equal(jump['H'],-jump['incident']-jump['returned'],'thin limit impulse')
    equal(jump['efficiency'],jump['e']/jump['incident'],'thin limit efficiency')
    equal(jump['mu'],16*(jump['H']**2+.25-jump['e']*(2+jump['e'])),'thin limit closed constraint')
    old = read(HERE/'prior/JUNCTION.json')['radiation_exit']
    equal(primary['config']['damping'],old['damping'],'inherited dissipation')


def check_gravity(data):
    for label,row in data['cases'].items():
        s = row['final']
        a,h,mu,e,p=(s[k] for k in ('a','H','wall_mass','energy','pressure'))
        equal(h*a,s['v'],'wall proper speed '+label)
        equal(e,s['scalar_density']+s['surface_em_density'],'wall energy '+label)
        equal(s['scalar_density'],.5*s['pi']**2+s['potential'],'scalar source '+label)
        equal(s['surface_em_density'],.5*s['k']*(s['wall_electric_squared']+
              s['wall_magnetic_squared']),'surface field energy '+label,absolute=1e-14)
        equal(s['wall_magnetic_squared'],4*s['wall_field']**2/a**4,'curved wall harmonic '+label)
        equal(s['wall_acceleration_over_radius'],(1+e)*(1-2*e-3*p)-1-mu/a**4-2*s['Tnn'],
              'temporal junction '+label,absolute=1e-12)
        equal(s['wall_input']-s['wall_return'],s['sheet_proper_flux'],'wall fluence '+label,absolute=2e-8)
        equal(a**3*e+s['sheet_expansion_work']-s['sheet_proper_flux'],s['sheet_balance'],
              'proper wall work '+label,absolute=1e-13)
        equal(mu-row['initial_mass']-s['mass_ledger'],s['wall_mass_residual'],
              'wall gravitational mass ledger '+label,absolute=1e-13)
        require(min(row['minimum_source_mass'],row['minimum_f'],row['minimum_characteristic_speed'])>0,
                'Nonpositive source mass, radial metric or wave speed '+label)
    for k in ('vacuum','constant_K','rising_constant_K'):
        require(data['cases'][k]['final']['scalar_density']==0, 'Unforced control gained energy '+k)
    vacuum=data['cases']['vacuum']['final']
    equal(vacuum['a'],math.sqrt(4-vacuum['t']**2),'vacuum trajectory',absolute=2e-10)
    for branch in ('loading','rising_potential'):
        refined=data['cases'][branch+'_320']
        for metric,ratios in refined['constraint_reduction_factors'].items():
            errors=[data['cases'][f'{branch}_{n}']['errors'][metric] for n in (80,160,320)]
            require(min(ratios)>3.,'Radial constraint failed mesh convergence '+branch)
            for observed,a,b in zip(ratios,errors[:-1],errors[1:],strict=True):
                equal(observed,a/b, 'constraint convergence '+branch)


def check_return_field(data):
    for label,row in data['cases'].items():
        require(row['max_fluid_density']==0 and row['minimum_mass']>0,
                'Conservative Maxwell scenario acquired a separate fluid or lost regular mass')
        for s in (row['peak'],row['bounce'],row['terminal']):
            a=s['a']
            equal(s['H']**2+1/a**2,s['e']*(2+s['e'])+s['mu']/a**4,'conservative closure',absolute=2e-6)
            equal(s['sheet_fraction'],s['sheet']/s['e'],'Maxwell fraction')
            equal(s['e'],s['sheet']+.5*s['pi']**2+s['V'],'energy components',absolute=1e-12)
        for epoch in row['expanding_epochs']:
            start,end=epoch['start'],epoch['end']
            equal(epoch['efolds'],math.log(end['a']/start['a']),'Maxwell-dominated expansion')
            require(min(start['sheet_fraction'],end['sheet_fraction'])>=.9-1e-8 and
                    epoch['efolds']>0,'Invalid expanding Maxwell window')
        for exchange in row.get('exchange',{}).values():
            transfer = exchange['scalar_to_sheet']
            equal(exchange['scalar_internal_change']+exchange['scalar_pressure_work'],-transfer,
                  'scalar release',absolute=5e-7)
            equal(exchange['sheet_internal_change']+exchange['sheet_pressure_work']+
                  exchange['outgoing_to_both_copies']-exchange['bulk_input'],transfer,
                  'Maxwell energy gain',absolute=5e-7)
        seed=row.get('crossing_seed_bound')
        if seed:
            require(seed['realized_K_ratio'] >= seed['required_K_ratio']*(1-1e-8) and
                    seed['measured_seed_density'] >= seed['minimum_seed_density']*(1-1e-8),
                    'Monotone postpulse seed threshold failed')
    require(data['cases']['primary']['expanding_epochs'] and
            not data['cases']['no_late_coupling']['expanding_epochs'], 'Late transfer control failed')


def check_forward(data):
    runs=data['runs']
    columns=data['trace_columns']
    for label,row in runs.items():
        require(min(row['min_lapse'],row['min_bulk_mass'],row['min_radial_gradient_norm'])>0,
                'Nonpositive metric, mass or characteristic frame '+label)
        a=row['end_radius']
        e=row['end_density']
        equal(a**4*row['end_sheet_energy_density'],row['end_comoving_sheet_energy'],
              'forward output scale '+label,absolute=1e-12)
        equal(a**4*row['end_radiation'],row['fluid_comoving_radiation_produced'],
              'comoving radiation transfer '+label,absolute=1e-6)
        equal(row['end_H']**2+1/a**2,e*(2+e)+row['end_mu']/a**4,
              'finite-source Friedmann equation '+label,absolute=2e-6)
        equal(row['wall_incident_proper_energy']-row['wall_outgoing_proper_energy'],
              a**3*e+row['wall_pressure_work'],'finite-source surface work '+label,absolute=2e-6)
        equal(row['comoving_wave_to_sheet']+row['comoving_scalar_to_sheet'],row['end_comoving_sheet_energy'],
              'Maxwell transfer ledger '+label,absolute=2e-6)
        for values in row.get('time_series',[]):
            s=dict(zip(columns,values,strict=True))
            equal(s['e'],s['pi']**2/2+s['V']+s['rad']+s['sheet'],'collapsed energy components')
            equal(s['sheet'],.5*s['K']*(s['X']**2+s['Bmag']**2),'reconstructed surface amplitude',
                  absolute=2e-12)
    for control in ('maxwell_no_return_coupling','maxwell_weaker','maxwell_stronger'):
        require(runs[control]['peak_exit_sheet_fraction']<.9,'Control unexpectedly reached Maxwell domination')
    wave=runs['maxwell_primary']
    require(wave['end_radiation']==0 and wave['postexit_em90_intervals']==1 and
            wave['longest_em90_expanding_efolds']>0 and wave['end_H']<0 and
            wave['at_peak_comoving_late_scalar_to_sheet']>0,'Missing conservative Maxwell epoch')
    peak=wave['at_peak_sheet_electric_density']+wave['at_peak_sheet_magnetic_density']
    equal(peak/wave['at_peak_sheet_total_density'],wave['peak_exit_sheet_fraction'],'peak Maxwell fraction')
    equal(peak*wave['at_peak_sheet_radius']**4,wave['at_peak_comoving_sheet_energy'],'peak energy')
    require(runs['fluid_primary']['radiation_expanding_efolds']>0,'Missing material-radiation comparison')
    vacuum=runs['vacuum']
    equal(vacuum['end_radius'],math.sqrt(4-vacuum['end_tau']**2),'finite-source vacuum control',
          absolute=1e-10)


def physics(products):
    require(products.keys() == PRODUCTS.keys(), 'Missing numerical product')
    check_fixed(products['fixed'])
    check_identities(products['identities'])
    check_boundary(products['loading']['boundary'])
    check_gravity(products['gravity'])
    check_return_field(products['loading']['return_field'])
    check_forward(products['loading']['forward'])


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
            with TemporaryDirectory(prefix='hdblast-curvature-') as tmp:
                regenerate(Path(tmp),args.jobs)
