#!/usr/bin/env python3
"""Fresh acceptance audit of saved scientific source/evidence, stdlib only.

This is separate acceptance evidence, not an additional registered scientific
replay process. No producer, upstream program or other checker is imported.
The retained continuous waveform theorem remains a conditional analytic input.
"""
import argparse
import csv
from fractions import Fraction as Q
import hashlib
import json
from math import comb
from pathlib import Path
import subprocess
import sys

PIN_SHA = '2507bbb4f279ec146fcea7342bbd17322f3ded6eb86eba346aeb3bd985765481'
CSV_SHA = 'd7ac199c89109c360c5da52b6d8693f4d34950f314fc42a5bde35a39d71859c5'
FIELDS = ('y','yprime','u','utau','uy','v','vy','Sp','SK','EN','dSp','dSK','dEN')
SOURCES = ('Sp','SK','EN','dSp','dSK','dEN')
CONTROLS = {
    'omit_endpoint_material_square': 'central_endpoint_signs',
    'wrong_EN_angular_factor': 'angular_norm_factors',
    'wrong_ratio_denominator': 'strict_norm_ratio',
    'drop_proper_clock_radius': 'proper_time_derivative_reconstruction',
    'stale_Sol_report_identity': 'saved_checker_evidence',
    'source_sign_is_invariant': 'scientific_scope',
}


class Rejected(Exception):
    def __init__(self, group, detail):
        self.group, self.detail = group, detail
        super().__init__(group+': '+detail)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bounds(x):
    return tuple(Q(v) for v in x)


def product(x, y):
    p = [a*b for a in x for b in y]
    return min(p), max(p)


def square(x):
    return (Q(0) if x[0] <= 0 <= x[1] else min(a*a for a in x),
            max(a*a for a in x))


def moment(a, b):
    """Mean of n1^(2a)n2^(2b) on S3 from its beta moments."""
    result = Q(1)
    for count in (a, b):
        for j in range(count):
            result *= Q(2*j+1, 2)
    for j in range(a+b):
        result /= 2+j
    return result


def audit(root, injection=None):
    groups, details = [], {}

    def need(group, ok, detail):
        if not ok:
            raise Rejected(group, detail)

    def finish(group, detail):
        groups.append(group)
        details[group] = detail

    g = 'immutable_input_inventory'
    pinpath = root/'inputs/INPUT_PINS.json'
    need(g, sha(pinpath) == PIN_SHA, 'final eighteen-input inventory changed')
    pins = json.loads(pinpath.read_text())
    need(g, len(pins['files']) == 18 and len({x['path'] for x in pins['files']}) == 18,
         'exactly eighteen distinct retained inputs required')
    for entry in pins['files']:
        rel = Path(entry['path'])
        need(g, not rel.is_absolute() and rel.parts[0] == 'inputs' and '..' not in rel.parts,
             'noncanonical input path')
        path = root/rel
        need(g, path.stat().st_size == entry['bytes'] and sha(path) == entry['sha256'], str(rel))
        origin = entry['original_public_source']
        need(g, len(origin['commit']) == 40 and origin['url'].startswith(
             'https://github.com/'+origin['repository']+'/blob/'+origin['commit']+'/'),
             'immutable public origin missing: '+str(rel))
    finish(g, {'files':18, 'pins_sha256':PIN_SHA, 'public_origin_metadata_checked':True,
               'remote_blobs_refetched':False})

    report_path = root/'computed/COORDINATE_SHELL_SOURCE_ENCLOSURE.json'
    report = json.loads(report_path.read_text())
    exact = report['exact']
    g = 'source_and_mode_identity'
    need(g, report['status'] == 'PASS' and report['producer_sha256'] == sha(root/'compute_coordinate_sources.py')
         and report['input_pins_sha256'] == PIN_SHA and report['runtime']['optimization'] == 0,
         'producer/report identity')
    for rel in ('computed/COORDINATE_SHELL_SOURCE_ENCLOSURE_OPTIMIZED.json',
                'validation/optimized-generator/COORDINATE_SHELL_SOURCE_ENCLOSURE.json'):
        other = json.loads((root/rel).read_text())
        need(g, other['runtime']['optimization'] == 1, 'optimized report mode')
        a, b = dict(report), dict(other)
        a['runtime'], b['runtime'] = dict(a['runtime']), dict(b['runtime'])
        a['runtime'].pop('optimization'); b['runtime'].pop('optimization')
        need(g, a == b, 'complete report disagreement except optimization: '+rel)
    table_path = root/'computed/COORDINATE_SHELL_SOURCE_CELLS.csv'
    meta = report['cell_data']
    need(g, sha(table_path) == CSV_SHA == meta['sha256'] and table_path.stat().st_size == meta['bytes'],
         'fixed scientific table identity')
    need(g, (root/'validation/optimized-generator/COORDINATE_SHELL_SOURCE_CELLS.csv').read_bytes()
         == table_path.read_bytes(), 'complete mode CSV disagreement')
    finish(g, {'csv_sha256':CSV_SHA, 'report_sha256':sha(report_path),
               'complete_JSON_compared_except_individual_mode_field':True})

    g = 'whole_cell_table_and_global_hulls'
    den, grid = meta['value_denominator'], meta['time_grid_denominator']
    need(g, den == 2**144 and grid == 4096 and meta['fields'] == list(FIELDS), 'grid/field schema')
    rows = list(csv.DictReader(table_path.read_text().splitlines()))
    need(g, len(rows) == meta['rows'] == 4353, 'complete continuous cell count')
    last, data = Q(0), []
    for row in rows:
        a, b = Q(int(row['left_grid']),grid), Q(int(row['right_grid']),grid)
        need(g, a == last < b, 'gap, overlap or nonpositive cell')
        values = {k: (Q(int(row[k+'_lo']),den),Q(int(row[k+'_hi']),den)) for k in FIELDS}
        need(g, all(x[0] <= x[1] for x in values.values()), 'reversed field interval')
        data.append((a,b,values)); last = b
    need(g, last == 3, 'observation window')
    for k in FIELDS:
        hull = min(v[k][0] for _,_,v in data), max(v[k][1] for _,_,v in data)
        need(g, hull == bounds(exact['global_envelopes'][k]), 'global hull '+k)
    finish(g, {'continuous_cells':4353, 'authenticated_scientific_endpoints':4353*13*2})

    g = 'angular_norm_factors'
    y4 = sum(Q((-1)**j)*comb(4,j)*moment(4-j,j) for j in range(5))
    z2sq = 2*moment(2,0)+2*moment(1,1)-2*moment(1,0)+Q(1,4)
    h4sq = y4-Q(1,36)-z2sq/4
    need(g, (y4,z2sq,h4sq) == (Q(3,40),Q(1,12),Q(19,720)), 'independent S3 moments')
    pi, rho, ell = bounds(exact['pi']), bounds(exact['rho']), Q(exact['ell'])
    expected = tuple(Q(3,20)*pi[i]**2*rho[i]**4*ell for i in (0,1))
    enfactor = Q(3,20) if injection == 'wrong_EN_angular_factor' else 2*h4sq
    expected_en = tuple(enfactor*pi[i]**2*rho[i]**4*ell for i in (0,1))
    need(g, bounds(exact['L2_prefactor']) == expected and bounds(exact['EN_L2_prefactor']) == expected_en,
         'complete Y squared and retained degree-four H4 factors differ')
    finish(g, {'Y_squared_factor_over_pi_squared':'3/20','H4_factor_over_pi_squared':'19/360',
               'proper_measure_radial_power':4})

    g = 'saved_source_norm_quadratures'
    for k in SOURCES:
        totals = [Q(0),Q(0)]
        for a,b,values in data:
            sq = square(values[k]); clock_hi = 1+2*(ell*b)**2
            totals[0] += (b-a)*sq[0]
            totals[1] += (b-a)*clock_hi*sq[1]
        row = exact['norms'][k]
        need(g, tuple(totals) == bounds(row['weighted_x_squared_integral']), 'all-cell integral '+k)
        factor = expected_en if k in ('EN','dEN') else expected
        need(g, row['angular_field'] == ('H4' if k in ('EN','dEN') else 'Y^2')
             and bounds(row['squared_L2']) == product(factor,totals), 'norm factor '+k)
    z = (3*ell)**2
    need(g, z < Q(1,100) and 4-18*z-2*z**3 > 0,
         'cosh-cubed polynomial upper bound on full domain')
    finish(g, {'all_exported_source_and_correction_squared_norms':6,
               'full_continuous_cell_quadratures':True})

    old = json.loads((root/'inputs/prior/INERTIAL_WAVEFORM_ENCLOSURE.json').read_text())
    shell = json.loads((root/'inputs/prior/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json').read_text())
    box = lambda name: bounds(shell['boxes'][name]['exact'])
    a,p,sigma,r2 = map(box,('a_b','p_b','sigma_second_b','H_squared'))
    mu = bounds(old['shell_parameter_intervals']['mu'])
    g = 'independent_central_metric_positive'
    A,B = bounds(old['zero_existence_brackets'][0])
    need(g, (A,B) == (Q('1.49998'),Q('1.49999')) and ell == Q(1,10**6), 'registered state')
    oldden = old['arithmetic']['value_denominator']; oldgrid = old['cell_data']['time_grid_denominator']
    history = list(csv.DictReader((root/'inputs/prior/INERTIAL_WAVEFORM_CELLS.csv').read_text().splitlines()))
    at_A, Hnorm, BH = [Q(0),Q(0)], Q(0), Q(0)
    for row in history:
        left,right = Q(row['left_grid'])/oldgrid,Q(row['right_grid'])/oldgrid
        H = (Q(row['H_lo'])/oldden,Q(row['H_hi'])/oldden)
        absH = max(map(abs,H)); Hnorm += absH*(right-left)
        BH += absH*((3-left)**2-(3-right)**2)/2
        endpoint = min(A,right)
        if left < endpoint:
            w = (endpoint-left)*(A-(left+endpoint)/2)
            for i in (0,1): at_A[i] += H[i]*w
    M, q = Q(5,2)*ell**2, Q(9,800000000000)
    coarse = ell/mu[0]*BH/(1-q)
    # Different central proof from the producer: evaluate F(A) by a direct
    # integral and use its global Lipschitz constant across the short bracket.
    fA = product((ell/mu[1],ell/mu[0]),at_A)
    fprime_abs = ell/mu[0]*Hnorm
    lower_y = fA[0]-fprime_abs*(B-A)-q*coarse
    lower_u = r2[0]*(1-(ell*B)**2)*lower_y
    need(g, M*9/2 == q < 1 and lower_y > 0 and lower_u > 0
         and a[0] > 0 and p[1] < 0 and sigma[0] > 0, 'uniform positive metric on closed bracket')
    U = r2[1]*coarse
    finish(g, {'u_lower_exact':str(lower_u), 'u_lower_display':float(lower_u),
               'global_abs_u_upper_display':float(U),
               'method':'direct F(A) integration, global F-prime Lipschitz bound, coarse Volterra residual'})

    g = 'proper_time_derivative_reconstruction'
    prefix = [Q(0),Q(0)]
    force = ell/mu[1],ell/mu[0]
    proper_rho = box('rho_b')
    ep_coarse = 3*M*coarse
    material_I2 = [Q(0),Q(0)]
    for inherited,(left,right,values) in zip(history,data):
        width = right-left
        H = Q(inherited['H_lo'])/oldden,Q(inherited['H_hi'])/oldden
        fp = product(force,(prefix[0]+min(Q(0),H[0]*width),
                            prefix[1]+max(Q(0),H[1]*width)))
        # y bounded by the independently rebuilt absolute-H resolvent.
        friction_abs = Q(3,2)*ell**2*right*coarse
        core = fp[0]-ep_coarse-friction_abs,fp[1]+ep_coarse+friction_abs
        radius = (Q(1),Q(1)) if injection == 'drop_proper_clock_radius' else (1/proper_rho[1],1/proper_rho[0])
        prefactor = product(product(r2,(1-(ell*right)**2,Q(1))),radius)
        derivative = product((prefactor[0]/ell,prefactor[1]/ell),core)
        need(g, derivative[0] <= values['utau'][0] <= values['utau'][1] <= derivative[1],
             'all-cell proper-time derivative lies outside independent coarse reconstruction')
        for i in (0,1): prefix[i] += H[i]*width
        v = Q(inherited['ell_S_lo'])/oldden/ell,Q(inherited['ell_S_hi'])/oldden/ell
        v2 = square(v)
        raw_material = -2*v2[1],-2*v2[0]
        lscaled,hscaled = raw_material[0]*den,raw_material[1]*den
        material = Q(lscaled.numerator//lscaled.denominator,den),Q(-((-hscaled.numerator)//hscaled.denominator),den)
        material2 = square(material)
        material_I2[0] += width*material2[0]
        material_I2[1] += width*(1+2*(ell*right)**2)*material2[1]
    finish(g, {'all_cells':4353,'chain_rule':'rho_b^(-3) C / ell times (yprime-(3/2)ell*tanh(ell*x)*y)',
               'comparison':'reported utau boxes contained in independent coarser absolute-H derivative boxes'})

    g = 'central_endpoint_signs'
    witnesses = old['sign_witnesses']; endpoint_results = {}
    for name,x in (('central_left',A),('central_right',B)):
        w = witnesses[name]; v = tuple(Q(z)/ell for z in w['ell_S'])
        need(g, Q(w['x']) == x and v[0]*v[1] > 0, 'nonzero scalar endpoint '+name)
        vmin,vmax = min(map(abs,v)),max(map(abs,v))
        spmax = -4*p[0]*U**2+2*sigma[1]*U*vmax
        if injection != 'omit_endpoint_material_square': spmax -= 2*vmin**2
        skmin = sigma[0]*vmin**2/6-12*a[1]*U**2
        need(g, spmax < 0 < skmin, 'endpoint signs '+name)
        endpoint_results[name] = {'Sp_upper_exact':str(spmax),'SK_lower_exact':str(skmin)}
    need(g, Q(witnesses['central_left']['H'][0]) > 0 > Q(witnesses['central_right']['H'][1])
         and Q(witnesses['central_left']['ell_S'][1]) < 0 < Q(witnesses['central_right']['ell_S'][0]),
         'retained actual scalar sign change')
    need(g, -4*p[1]*lower_u**2 > 0 and -12*a[0]*lower_u**2 < 0,
         'opposite source signs at every interior scalar zero')
    finish(g, {'endpoint_bounds':endpoint_results,
               'inference':'Each coefficient has one zero in (A,x*) and one in (x*,B), by continuity; these intervals are disjoint.',
               'uniqueness_claimed':False})

    g = 'strict_norm_ratio'
    num = bounds(exact['norms']['dSp']['squared_L2'])
    denom = bounds(exact['norms']['Sp_material']['squared_L2'])
    need(g, denom == product(expected,material_I2), 'material denominator independently reconstructed from retained scalar cells')
    ratio = num[1]/denom[1 if injection == 'wrong_ratio_denominator' else 0]
    need(g, denom[0] > 0 and ratio == Q(exact['metric_correction_squared_norm_ratio_upper']['Sp'])
         and ratio < Q(9,10**13)**2, 'interval denominator and strict unsquared bound')
    finish(g, {'squared_ratio_upper_exact':str(ratio),'unsquared_strict_upper':'9/10000000000000'})

    g = 'scientific_scope'
    scope = dict(report['scope'])
    if injection == 'source_sign_is_invariant': scope['coordinate_source_sign_is_gauge_invariant'] = True
    need(g, scope['conditional_on_registered_analytic_inputs'] is True
         and all(scope[k] is False for k in ('full_U2_evaluated','nonlinear_remainder_established',
                    'coordinate_source_sign_is_gauge_invariant','cosmological_evidence'))
         and scope['EN_scope'] == 'Retained degree-four reduced Robin forcing only.',
         'conditional fixed-representative coefficient scope')
    finish(g, {'retained_PDE_proof_reproduced':False,'physical_force_or_U2_inferred':False,
               'upstream_support_inverse_constraint_gauge_corner_hypotheses_remain':True})

    g = 'saved_checker_evidence'
    evidence = []
    for directory, script, base, key in (
        ('astra','check_coordinate_sources.py','ASTRA_COORDINATE','checker_sha256'),
        ('sol','check_coordinate_shell_sources.py','SOL','checker_sha256'),
        ('referee','check_full_source_formulation.py','FULL_SOURCE','source_script_sha256')):
        codehash = sha(root/directory/script)
        for mode in ('NORMAL','OPTIMIZED'):
            rel = directory+'/'+base+'_'+mode+'.json'
            receipt = json.loads((root/rel).read_text())
            need(g, receipt['status'] == 'PASS' and receipt[key] == codehash, 'stale source receipt '+rel)
            if directory == 'sol':
                auth = dict(receipt['primary_authentication'])
                if injection == 'stale_Sol_report_identity': auth['report_sha256'] = '0'*64
                need(g, auth['report_sha256'] == sha(report_path) and auth['csv_sha256'] == CSV_SHA
                     and auth['scientific_endpoints_authenticated'] == 113178, 'Sol complete primary identity '+rel)
            evidence.append({'path':rel,'sha256':sha(root/rel),'checker_sha256':codehash})
    finish(g, evidence)
    reviewed_docs = ('README.md','COORDINATE_SHELL_SOURCES.md','CLAIMS.json','REPRODUCE.md',
        'astra/INDEPENDENT_COORDINATE_SOURCE_PROOF.md','sol/COORDINATE_SOURCE_AUDIT.md',
        'referee/SCIENTIFIC_SCOPE_AND_ASSUMPTION_AUDIT.md',
        'referee/FULL_BULK_SOURCE_AND_INTERIOR_NEXT_STEP.md','literature/LITERATURE_CHECK.md',
        'literature/LITERATURE_LEDGER.json')
    return {'schema':'hdblast-fresh-coordinate-acceptance-v1','status':'PASS',
            'acceptance_evidence_only':True,'registered_scientific_replay_process':False,
            'runtime':{'optimization':sys.flags.optimize,'python':sys.version.split()[0]},
            'checker_sha256':sha(Path(__file__)),'root':str(root),'check_groups':groups,
            'check_group_count':len(groups),'details':details,
            'reviewed_document_sha256':{rel:sha(root/rel) for rel in reviewed_docs}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--output',type=Path)
    parser.add_argument('--inject-error',choices=sorted(CONTROLS))
    parser.add_argument('--run-controls',action='store_true')
    args = parser.parse_args()
    try:
        if args.run_controls:
            results = []
            for optimize in (False,True):
                for name,expected in CONTROLS.items():
                    command = [sys.executable,'-I','-B']+(['-O'] if optimize else [])+[
                        str(Path(__file__).resolve()),'--root',str(args.root.resolve()),'--inject-error',name]
                    run = subprocess.run(command,capture_output=True,text=True)
                    receipt = json.loads(run.stdout)
                    if run.returncode != 1 or receipt.get('failed_group') != expected:
                        raise Rejected('control_diagnostic',name+': '+run.stdout+run.stderr)
                    results.append({'control':name,'optimized':optimize,'exit_code':1,
                                    'failed_group':expected,'status':'EXPECTED_REJECTION'})
            result = {'schema':'hdblast-fresh-acceptance-controls-v1','status':'PASS',
                      'checker_sha256':sha(Path(__file__)),'controls':results,
                      'rejected_acceptance_control_processes':len(results),
                      'registered_scientific_replay_processes_added':0}
        else:
            result = audit(args.root.resolve(),args.inject_error)
        code = 0
    except Rejected as exc:
        result = {'status':'FAIL','failed_group':exc.group,'detail':exc.detail,'injected_error':args.inject_error}
        code = 1
    raw = json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.output: args.output.write_text(raw)
    print(raw,end='')
    return code


if __name__ == '__main__':
    raise SystemExit(main())
