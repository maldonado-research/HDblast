"""Reproduce the Euclidean response, timelike spectrum and pulse tests on the static tuned leaves."""
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
BASE = "1475602d67399699bf4461d72c7a8c1742c8ab08"
PINS_SHA = "c01d9fa9465a5de16e62803804d93f5e8adcad16edd4aaf44a04b0191f9b4f9b"
PRODUCTS = {name: f"{name}/{name.upper()}.json" for name in ("euclidean", "retarded", "bounds", "pulse", "forces")}
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
    require(pins['base'] == BASE and pins['previous_doi'] == '10.5281/zenodo.24383815',
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
    # Solver residuals and mesh-difference diagnostics are compared in their own units.
    # Tiny spectral weights and pulse tails retain a relative accuracy requirement.
    # QUADPACK's own error estimate has the units of the integrated coefficient.
    # check_static enforces the absolute/relative tolerance requested for the integral.
    if path[0]=='euclidean' and path[-1]=='first_correction_quadrature_error_estimate':
        return actual >= 0 and expected >= 0
    # This bound independently checked the incomplete beta norm by direct quadrature.
    if path[0]=='bounds' and path[-1]=='check_norm_absolute_error':
        return 0 <= actual <= 5e-14 and 0 <= expected <= 5e-14
    floors = {
        'relative_bracket_and_integration_difference': 3e-12, 'refinement_absolute_in_scaled_fraction': 2e-11,
        'refinement_relative': 5e-10, 'independent_relative': 5e-10,
        'integration_order_difference': 5e-15, 'quadrature_relative': 2e-12,
        'pole_error': 2e-11, 'zero_mass_squared': 2e-10, 'relative_difference': 4e-11,
        'independent_box_relative': 5e-11, 'relative_extrapolation_change': 8e-11,
        'relative_finest_change': 2e-11, 'relative_box_change': 2e-11,
        'relative_larger_box_change': 2e-11, 'raw_relative_box_change': 2e-11,
        'raw_relative_larger_box_change': 2e-11, 'spectral_mesh_relative': 2e-12,
        'radius_quadrature_refinement': 2e-12, 'independent_max_relative': 2e-9,
        'refinement_max_relative': 2e-9, 'current_maximum_error_above_half': 3e-10,
        'spectral_current_maximum_relative_above_half': 3e-10}
    return math.isclose(actual, expected, rel_tol=2e-7, abs_tol=floors.get(path[-1],0.))


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


def varied_energies(phi_b, momentum, intervals=600):
    """Upper and lower admittance estimates from independent compact trial fields."""
    import numpy as np
    from scipy.linalg import solve_banded
    # The primal potential and the dual flux vanish beyond the last node.
    y_b = -math.atanh(phi_b)
    edges = np.linspace(y_b, y_b+12, intervals+1)
    height = (edges[-1]-edges[0])/intervals
    nodes, weights = np.polynomial.legendre.leggauss(5)
    right = (1+nodes)/2
    left = 1-right
    y = edges[:-1,None]+height*right
    def warp_primitive(q):
        return q + 2*np.log(np.cosh(q))/3 + np.tanh(q)**2/6
    warp = -(warp_primitive(y)-warp_primitive(y_b))/3
    dielectric = np.cosh(y)**-4
    p = dielectric*np.exp(2*warp)
    q = dielectric*momentum**2
    measure = height*weights/2

    def functional(stiffness, mass):
        derivative = (stiffness @ measure)/height**2
        diagonal_left = derivative + mass @ (left**2*measure)
        diagonal_right = derivative + mass @ (right**2*measure)
        offdiagonal = -derivative + mass @ (left*right*measure)
        bands = np.zeros((3,intervals-1))
        bands[1] = diagonal_right[:-1]+diagonal_left[1:]
        bands[0,1:] = offdiagonal[1:-1]
        bands[2,:-1] = offdiagonal[1:-1]
        boundary = np.zeros(intervals-1)
        boundary[0] = -offdiagonal[0]
        values = np.r_[1,solve_banded((1,1),bands,boundary),0]
        slope = np.diff(values)/height
        trial = values[:-1,None]*left+values[1:,None]*right
        return float(np.sum((stiffness*slope[:,None]**2+mass*trial**2)*measure))
    return 1/functional(1/q,1/p), functional(p,q)


def check_static(euclidean, retarded):
    grid_fields = {float(v) for v in euclidean['leaves']}
    require(grid_fields == {-.5,0.,.5}, 'Different static leaves')
    require(grid_fields == {float(v) for v in retarded['branches']}, 'Different timelike leaves')
    for field, leaf in euclidean['leaves'].items():
        phi = float(field)
        branch = retarded['branches'][f'{phi:g}']
        norm = (1+phi)**2*(2-phi)/3
        moment = leaf['first_dtn_correction']
        require(leaf['first_correction_quadrature_error_estimate'] <= max(2e-14,2e-13*abs(moment)),
                'Estimated quadrature error exceeds the integrated coefficient tolerance')
        equal(leaf['one_copy_norm'], norm, 'Zero-mode normalization')
        equal(leaf['independent_norm_integral'], norm, 'Independent zero-mode norm')
        equal(branch['half_norm'],norm, 'Timelike normalization')
        equal(branch['pole_residue'],1/norm, 'Boundary pole residue')
        equal(leaf['finite_propagator_after_zero_mode_at_zero'],moment/norm**2,'First finite response')
        for label, row in leaf['momenta'].items():
            q = float(label)
            kernel = row['dtn']
            require(0 < kernel < norm*q*q, 'Nonpositive pole or continuum fraction')
            equal(row['brane_propagator']*kernel, 1., 'Dirichlet and Neumann reciprocity')
            equal(row['extra_response_fraction'],1-kernel/(norm*q*q),'Nonzero-mode fraction')
            bound = q*q*moment/(norm+q*q*moment)
            equal(row['spectral_upper_bound'],bound,'Sharp trial-function response bound')
            require(row['extra_response_fraction'] < bound+3e-11, 'Trial-function bound fails')
            # Independent compact primal and dual finite-element trial fields. Checking
            # intermediate q avoids loss of resolution in either cancellation regime.
            if .029 < q <= 3:
                lower, upper = varied_energies(phi,q,intervals=2400)
                require(lower <= kernel <= upper and upper/lower-1 < 1e-4,
                        f'Complementary variational brackets fail at {phi}, {q}')
        psi_power = branch['threshold_power']
        equal(psi_power, 6.2, 'Warped Maxwell threshold power')
        for row in branch['samples']:
            mass, density = row['mass'], row['density']
            require(mass > 0 and density > 0, 'Outgoing spectral measure must be positive')
            equal(row['green_imaginary_from_flux'], math.pi*density/(2*mass), 'Retarded spectral normalization')
            equal(row['threshold_ratio'],density/(branch['threshold_coefficient']*mass**psi_power),'Threshold response')
        for label, transformed in branch['continuum_stieltjes'].items():
            q = float(label)
            direct = (1/leaf['momenta'][str(q)]['dtn'] - 1/(norm*q*q)) if q else moment/norm**2
            equal(transformed['euclidean_continuum'],direct,'Pole-subtracted static response',rel=3e-10)
            equal(transformed['timelike_transform'],direct,'Independent timelike dispersion integral',rel=4e-6)
        for current in branch['gaussian_current'].values():
            ratio = current['continuum_over_pole']
            require(ratio>0,'Gaussian wall current must excite the continuum')
            equal(current['continuum_fraction'],ratio/(1+ratio),'Gaussian work normalization')
        selected = euclidean['selected_weight_leaves'][field]
        alternate = retarded['alternate_weight'][f'{phi:g}']
        equal(alternate['half_norm'], selected['one_copy_norm'], 'Alternate weight zero-mode residue')
        equal(alternate['threshold_power'], 1, 'Logarithmic-response threshold')
        log_coeff = selected['response_fraction_logarithmic_coefficient']
        equal(alternate['threshold_coefficient']*selected['one_copy_norm'], log_coeff,
              'Threshold density from static logarithmic response')
        require(abs(alternate['samples'][0]['threshold_ratio']-1)<.001,'Alternate threshold not recovered')
        first = branch['samples'][0]
        require(abs(first['threshold_ratio']-1)<.01, 'Long-wave threshold coefficient')


def check_bounds(euclidean, bounds):
    require(set(bounds['leaves']['2']) == set(euclidean['leaves']), 'Missing weight-two leaves')
    expected_contrasts = read(HERE/'prior/origin/ORIGIN.json')['branches']
    for branch, change in bounds['contrasts']['2'].items():
        equal(change, expected_contrasts[branch]['response']['log_I_change'], 'Previously calibrated contrast',absolute=5e-8)
    equal(bounds['positive_contrast_upper_bound'],math.log(2)/2,'Symmetric-weight upper bound')
    for field, row in bounds['null_lift']['branches'].items():
        endpoint = row['endpoints']['first_minimum']
        phi = round(endpoint['phi_b'], 12)
        alternate = euclidean['selected_weight_leaves'][str(phi)]
        equal(endpoint['norm'],alternate['global_one_copy_norm'],'Geometric weight in common global normalization')
        equal(endpoint['gravity_norm'],2*alternate['one_copy_norm'],'Gravity and localized vector normalization')
        equal(endpoint['marginal_relative_loss_coefficient'],alternate['response_fraction_logarithmic_coefficient'],'Logarithmic fraction coefficient')
    for field, leaf in bounds['leaves']['2'].items():
        response = euclidean['leaves'][field]
        equal(leaf['single_copy_norm'],response['one_copy_norm'],'Static normalizations')
        equal(leaf['j_over_norm'],response['correction_over_norm'],'Euclidean inverse spectral moment')
        equal(leaf['one_percent_kernel_sufficient_q']**2*leaf['j_over_norm'],1/99,'One percent variational condition')
        for key, row in leaf['points'].items():
            if key in response['momenta']:
                actual = response['momenta'][key]['extra_response_fraction']
                require(row['relative_loss_lower']-3e-11 <= actual <= row['relative_loss_upper']+3e-11,
                        'Independent Euclidean response not within scalar-weight bound')


def check_forces(euclidean, retarded, forces):
    require(set(forces['leaves']) == set(retarded['branches']), 'Static force leaves differ')
    for field, leaf in forces['leaves'].items():
        phi = float(field)
        spacelike = euclidean['leaves'][str(phi)]
        spectrum = retarded['branches'][field]
        equal(leaf['one_copy_norm'],spacelike['one_copy_norm'],'Potential-field normalization')
        moment = leaf['finite_spectral_moment']
        equal(moment,spacelike['correction_over_norm'],'Potential inverse-mass moment')
        exponent = spectrum['threshold_power']+1
        leading = leaf['one_copy_norm']*spectrum['threshold_coefficient']*math.gamma(exponent)
        potential = leaf['potential']
        force = leaf['force']
        equal(potential['asymptotic_coefficient'],leading,'Long-distance potential coefficient')
        equal(force['asymptotic_coefficient'],(1+exponent)*leading,'Long-distance force coefficient')
        radii = sorted(float(key) for key in potential['samples'])
        u = [potential['samples'][f'{r:g}']['excess_over_zero_mode'] for r in radii]
        f = [force['samples'][f'{r:g}']['excess_over_zero_mode'] for r in radii]
        for radius, value, derivative in zip(radii,u,f):
            require(0 < value < derivative,'Positive Yukawa contribution to force and potential')
            require(value*radius**2 <= moment*4*math.exp(-2),'Zero-frequency bound on potential')
            x = 1+math.sqrt(3)
            require(derivative*radius**2 <= moment*x*x*(1+x)*math.exp(-x),'Zero-frequency bound on force')
        # Every positive spectral measure gives a log-convex Laplace transform.
        # For each interior tabulated separation, adjacent secants bound d log U/dR.
        slopes = [math.log(right/left)/(r-l) for l,r,left,right in zip(radii,radii[1:],u,u[1:])]
        for i in range(1,len(radii)-1):
            log_gradient = (1-f[i]/u[i])/radii[i]
            require(slopes[i-1] <= log_gradient <= slopes[i],'Force conflicts with potential slope')


def check_pulses(euclidean, retarded, pulse, selected):
    trajectories = pulse['selected_branches' if selected else 'branches']
    waves = retarded['alternate_weight' if selected else 'branches']
    spaces = euclidean['selected_weight_leaves' if selected else 'leaves']
    require(trajectories.keys()==waves.keys(), 'Moving-wall leaves differ')
    for field, row in trajectories.items():
        phi = float(field)
        space = spaces[str(phi)]
        spectral = waves[field]
        power = spectral['threshold_power']
        norm = row['half_norm']
        g = row['bulk_weight_at_wall']
        equal(norm,space['one_copy_norm'],'Moving-wall zero norm')
        equal(g,space['wall_weight' if selected else 'shell_weight'],'Moving-wall Maxwell weight')
        if selected:
            require(row['inverse_mass_squared_moment'] is None and row['tail_moment_bound']>0,
                    'Marginal throat uses a resolvent bound')
        else:
            equal(row['inverse_mass_squared_moment'],space['finite_propagator_after_zero_mode_at_zero'],'Moving-wall inverse moment')
        equal(row['fixed_k_slow_duration_power'],-(power+1)/2,'Slow motion energy exponent')
        independent = row['independent_box']
        estimates = independent['estimates']
        for key, case in row['cases'].items():
            k, width = case['k'], case['width']
            energy = case['energy_coefficient']
            factor = g**2/(4*norm)
            slow = factor*math.pi*spectral['threshold_coefficient']/2*math.gamma((power+5)/4)
            equal(case['slow_leading'],slow*(2*k/width)**((power+1)/2),'Moving-wall large-duration asymptote')
            require(energy>0 and 0<=case['omitted_high_mass_bound']<energy*1e-8,'Uncontrolled pulse tail')
            if 'static_potential_excess' in case:
                equal(case['static_separation']**2,width/(2*k),'Slow wall and static probe scales')
                equal(case['pulse_to_potential'],energy/case['static_potential_excess'],'Slow-pulse and force amplitudes')
                analytic = (g/norm)**2*math.pi*math.gamma((power+5)/4)/(8*math.gamma(power+1))
                equal(case['pulse_to_potential_asymptote'],analytic,'Pulse and static asymptotic relation')
                if width == 4096:
                    equal(case['pulse_to_potential'],analytic,'Long-duration pulse-to-potential limit',rel=.03)
                    equal(energy,case['slow_leading'],'Measured slow-wall law',rel=.03)
            if key in estimates:
                box = estimates[key]['energy']
                require(abs(energy/box-1)<2e-6,'Outgoing and finite-volume pulse energies differ')
                equal(case['independent_box_relative'],abs(energy/box-1),'Reported box discrepancy',absolute=2e-13)
                column = 'corrected_energy' if selected else 'energy'
                coarse = independent['grids'][2]['entries'][key][column]
                fine = independent['grids'][3]['entries'][key][column]
                equal(box,(4*fine-coarse)/3,'Independent continuum extrapolation')
        for grid in independent['grids']:
            require(abs(grid['zero_mass_squared'])<1e-8 and abs(grid['zero_overlap']-1)<1e-8,'Massless mode in finite box')
            if grid['time_check']:
                tc = grid['time_check']
                direct = grid['entries'][f"{tc['k']:g}/{tc['width']:g}"]['energy']
                equal(tc['energy'],direct,'Explicit time-domain forced-work check',rel=2e-8)


def physics(products):
    require(products.keys() == PRODUCTS.keys(), 'Missing numerical product')
    check_static(products['euclidean'], products['retarded'])
    check_bounds(products['euclidean'], products['bounds'])
    check_forces(products['euclidean'], products['retarded'], products['forces'])
    for selected in (False, True):
        check_pulses(products['euclidean'], products['retarded'], products['pulse'], selected)


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
            with TemporaryDirectory(prefix='hdblast-bulk-') as tmp:
                regenerate(Path(tmp),args.jobs)
