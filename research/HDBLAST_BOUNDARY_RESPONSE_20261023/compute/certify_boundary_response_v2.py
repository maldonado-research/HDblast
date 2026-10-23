#!/usr/bin/env python3
"""One conditional actual-source, shell-reaching auxiliary mean enclosure.

Run with python -I -S -B (optionally -O).  All input trees are authenticated
before importing their numerical code, and again after all evaluations.
Exact clipped areas multiply whole-box source intervals; no point quadrature,
comparator substitution, evaluated feedback waveform, or sign premise is used.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path, PurePosixPath
import sys
import time

WORK = Path('/workspace/hdblast-research-work')
OUTPUT_ROOT = WORK/'boundary-response-20261023/compute'
OLD = WORK/'geometric-bulk-source-20261016/capsule-source'
NEW = WORK/'forced-response-20261023/package-work/attempt-2/source'
FOLLOWON = WORK/'forced-response-20261023'
OLD_SHA = '3ff1f95c151d2da6a1ab5e9b54718c2bb325c10d9c62a7ed9f97ef2bb36f0bb2'
NEW_SHA = '72924328a4e79cc77f0da90999b8dd56956bfdebd5e614e15c3c0498c3a8ae76'
PINS = {
 'followon-exact-geometry/exact_geometry.py': '992dd82e3b48dc7ea29e7be90957daf4ba1ae44e8f629ed09db92fa240685b72',
 'followon-boundary-formulation/NEXT_BOUNDARY_CALCULATION.md': '903699df349ed1299252aa68ddec8b4e38f81730ae6418a0550da85222c8b5e2',
 'followon-boundary-formulation/ONE_POINT_FEEDBACK_AUDIT.md': '4a58981b9cd8ca0add3548a3fd01185a76cb791e06b0ccc1f17fece1bef3885e',
}
A, B, D = Q(3,2), Q(1,4), Q(5,4)
SIGMA, M = Q(203,1000000), Q(1,2000)
GLOBAL_F, GLOBAL_G = Q(17387513156573), Q(3851400)
E_TRACE = Q(14866325944167915,3992782)


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def unique_object(pairs):
    out = {}
    for key, value in pairs:
        need(key not in out, 'duplicate JSON key: '+key)
        out[key] = value
    return out


def canonical(name):
    pure = PurePosixPath(name)
    need(name == pure.as_posix() and not pure.is_absolute() and
         name not in ('', '.') and '..' not in pure.parts and '\\' not in name,
         'noncanonical manifest path')
    return name


def authenticate_tree(root, manifest_sha, count):
    need(root.is_dir() and not root.is_symlink(), 'source root regular directory')
    manifest_path = root/'MANIFEST.json'
    need(manifest_path.is_file() and not manifest_path.is_symlink(), 'regular manifest')
    raw = manifest_path.read_bytes()
    need(sha(raw) == manifest_sha, 'manifest external pin')
    manifest = json.loads(raw, object_pairs_hook=unique_object)
    expected = {}
    for row in manifest['files']:
        name = canonical(row['path'])
        need(name not in expected and name != 'MANIFEST.json', 'duplicate member')
        expected[name] = row
    actual = set()
    total = len(raw)
    for path in root.rglob('*'):
        need(not path.is_symlink() and (path.is_file() or path.is_dir()), 'nonregular member')
        if path.is_file():
            name = path.relative_to(root).as_posix()
            actual.add(name)
            if name != 'MANIFEST.json':
                need(name in expected, 'unexpected member: '+name)
                data = path.read_bytes()
                row = expected[name]
                need(len(data) == row['bytes'] and sha(data) == row['sha256'], 'member pin: '+name)
                total += len(data)
    need(actual == set(expected) | {'MANIFEST.json'} and len(actual) == count, 'exact source closure')
    return {'manifest_sha256': manifest_sha, 'files_including_manifest': count,
            'total_input_bytes': total, 'complete_source_authenticated': True}


def authenticate_inputs():
    out = {'geometric169': authenticate_tree(OLD, OLD_SHA, 169),
           'forced49': authenticate_tree(NEW, NEW_SHA, 49), 'followon': {}}
    for name, pin in PINS.items():
        path = FOLLOWON/name
        need(path.is_file() and not path.is_symlink(), 'regular followon input')
        data = path.read_bytes()
        need(sha(data) == pin, 'followon input pin: '+name)
        out['followon'][name] = {'bytes':len(data),'sha256':pin}
    return out


def import_path(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def row(value):
    return [str(value.lo), str(value.hi)]


def add_interval(left, right):
    return (left[0]+right[0], left[1]+right[1])


def scale_interval(value, scale):
    lo, hi = value
    return (scale*lo,scale*hi) if scale >= 0 else (scale*hi,scale*lo)


def interval_strings(value):
    return [str(v) for v in value]


def polygon_cell(geometry, polygon, base, label):
    cell = geometry._cell(polygon, base.mesh_index, base.region, False)
    if cell is None:
        return None
    return cell, label


def build_active_cells(geometry, n):
    mesh = geometry.build_geometry(A, B, A/n, first_moments=False)
    active, zero = [], []
    for base in mesh.dminus + mesh.strip:
        outside = geometry._clip(base.polygon,Q(1),Q(1),Q(1))
        if len(outside) >= 3:
            zarea, _ = geometry._measure(outside,False)
            if zarea:
                zero.append((base.region,zarea))
        # Keep xi=s+u>=1.  The exact standing zero wedge removes the rest.
        polygon = geometry._clip(base.polygon,Q(-1),Q(-1),Q(-1))
        # Split at eta=s-u=0. Error eligibility still uses the whole bbox.
        for ns,nu,label in ((Q(1),Q(-1),'eta_nonpositive'),
                            (Q(-1),Q(1),'eta_nonnegative')):
            entry = polygon_cell(geometry, geometry._clip(polygon,ns,nu,Q(0)), base,label)
            if entry is not None:
                active.append(entry)
    areas = {region: sum((c.area for c,_ in active if c.region==region),Q(0))
             for region in ('Dminus','strip')}
    zero_areas = {region:sum((area for r,area in zero if r==region),Q(0))
                  for region in ('Dminus','strip')}
    need(areas['Dminus']==Q(9,32) and areas['strip']==Q(11,16), 'active clipped areas')
    need(zero_areas=={'Dminus':Q(1,2),'strip':Q(0)}, 'zero wedge area')
    need(areas['Dminus']+areas['strip']/2==Q(5,8), 'active Robin kernel weight')
    return active, areas, zero_areas


def encoded_interval(coefficients, terms, m, s, bg):
    values = dict(m)
    values.update(s)
    values.update(bg)
    return coefficients.evaluate_encoded(terms,values,{})


def boundary_partition(n):
    # [0,1] is exactly zero. Positive history uses n equal intervals plus
    # the inherited 1/64 bump-tail transition, all with exact endpoints.
    points = {Q(1),D,Q(65,64)} | {Q(1)+i*(D-1)/n for i in range(n+1)}
    points = sorted(v for v in points if 1<=v<=D)
    return tuple(zip(points,points[1:]))


def run(n):
    need(n in (4,8,16), 'bounded initial/refinement partitions only')
    need(sys.flags.isolated and sys.flags.no_site and sys.flags.dont_write_bytecode,
         'use python -I -S -B')
    own_pin = sha(Path(__file__).read_bytes())
    authentication = authenticate_inputs()  # MUST precede numeric imports.
    geometry = import_path('boundary_exact_geometry',FOLLOWON/'followon-exact-geometry/exact_geometry.py')
    sys.path.insert(0,str(OLD/'numeric'))
    sys.path.insert(0,str(OLD/'projection'))
    import actual_first_jets as jets
    import fast_dyadic
    import source_coefficients as coefficients
    import shell_source_values as shell
    import angular_generic as angular
    need(tuple(angular.inner(s,angular.SCALARS[0]) for s in angular.SCALARS)==(Q(1),Q(0),Q(0)),
         'yy0 equals probability-sphere mean')
    accepted = json.loads((OLD/'numeric/ACTUAL_SOURCES_NORMAL.json').read_text(), object_pairs_hook=unique_object)
    need(Q(accepted['global_cover']['summary']['positive_bundle_bulk_forcing_H_bound']['F_H_upper'])==GLOBAL_F,
         'retained actual full23 forcing cap')
    data, errors, _ = jets.inputs()
    exact_maps = coefficients.load_maps(OLD/'projection/GEOMETRIC_COEFFICIENT_MAPS.json')
    F_terms = exact_maps['maps']['metric_forcing']['yy'][0]
    QC_terms = exact_maps['maps']['QC']['y'][0]
    active, areas, zero_areas = build_active_cells(geometry,n)
    convert = lambda d: {key:fast_dyadic.D(value) if hasattr(value,'lo') and hasattr(value,'hi') else value
                         for key,value in d.items()}
    totals = {'Dminus':(Q(0),Q(0)),'strip':(Q(0),Q(0))}
    bulk_rows = []
    started = time.monotonic()
    for index,(cell,label) in enumerate(active):
        slo,shi,ulo,uhi = cell.bbox
        need(0<=slo<=shi<=A and 0<=ulo<=uhi and shi+uhi<8, 'whole oracle box domain')
        m,s,bg = jets.actual_jet(jets.I(slo,shi),jets.I(ulo,uhi),mode='actual')
        value = encoded_interval(coefficients,F_terms,convert(m),convert(s),convert(bg))
        contribution = scale_interval((value.lo,value.hi),cell.area*cell.lambda_wyy_weight)
        totals[cell.region] = add_interval(totals[cell.region],contribution)
        bulk_rows.append({'index':index,'mesh_index':list(cell.mesh_index),'region':cell.region,
            'characteristic_split':label,'polygon':[[str(v) for v in p] for p in cell.polygon],
            'bbox_s':[str(slo),str(shi)],'bbox_u':[str(ulo),str(uhi)],
            'exact_area':str(cell.area),'wave_kernel_weight':str(cell.lambda_wyy_weight),
            'bbox_exterior_error_eligible':cell.bbox_exterior_error_eligible,
            'actual_error_regime':'exterior' if ulo>=shi else 'global',
            'actual_source_inputs':{'m':jets.serialize(m),'s':jets.serialize(s),'background':jets.serialize(bg)},
            'actual_Fyy_probability_mean_source_interval':row(value),
            'signed_bulk_contribution_interval':interval_strings(contribution)})
        print(json.dumps({'event':'bulk_cell','n':n,'cell':index+1,'cells':len(active),
                          'wall_seconds':round(time.monotonic()-started,3)}),flush=True)
    shell_rows = []
    shell_G_integral = (Q(0),Q(0))
    for index,(slo,shi) in enumerate(boundary_partition(n)):
        m,s,bg = jets.actual_jet(jets.I(slo,shi),jets.I(0),mode='actual')
        # The same actual boundary jet supplies every term, including the
        # varying-metric-lowered QC.  Never normally differentiate eliminated
        # physical shell rows.  QK is the ONE-copy OUTWARD affine Hessian.
        eta0 = jets.I(*data['shell']['boxes']['eta_b']['exact'])
        raw = shell.raw_fields(m['A'],m['D'],m['E'],m['L'],m['M'],s['s'],bg['a'],bg['p'],4*(1+eta0))
        qk = shell.projected(raw)
        qc = encoded_interval(coefficients,QC_terms,convert(m),convert(s),convert(bg))
        qk00,qkab = fast_dyadic.D(qk['00'][0]),fast_dyadic.D(qk['ab'][0])
        bracket = -qk00+3*qkab+qc
        gyy = 2*fast_dyadic.D(bg['ell'])*fast_dyadic.D(bg['rho'])*bracket
        gcontribution = scale_interval((gyy.lo,gyy.hi),shi-slo)
        shell_G_integral = add_interval(shell_G_integral,gcontribution)
        shell_rows.append({'index':index,'s':[str(slo),str(shi)],'r':['0','0'],
            'length':str(shi-slo),'actual_error_regime':'global',
            'actual_source_inputs':{'m':jets.serialize(m),'s':jets.serialize(s),'background':jets.serialize(bg)},
            'ell':str(bg['ell']),'rho':row(bg['rho']),
            'one_copy_outward_QK00_0':row(qk00),'one_copy_outward_QKab_0':row(qkab),
            'varying_metric_lowered_QC_y_0':row(qc),
            'lorentzian_trace_plus_QC_interval':row(bracket),'actual_Gyy_mean_interval':row(gyy),
            'G_integral_contribution_interval':interval_strings(gcontribution),
            'signed_negative_LG_contribution_interval':interval_strings(scale_interval(gcontribution,Q(-1)))})
        print(json.dumps({'event':'boundary_cell','n':n,'cell':index+1,'cells':len(boundary_partition(n)),
                          'wall_seconds':round(time.monotonic()-started,3)}),flush=True)
    bulk = add_interval(totals['Dminus'],totals['strip'])
    lg = scale_interval(shell_G_integral,Q(-1))
    s_interval = add_interval(bulk,lg)
    qd = A*SIGMA+3*M
    b_domain = (A*GLOBAL_F+GLOBAL_G)/(1-qd)
    e_bulk = Q(5,8)*SIGMA*b_domain
    e_robin = M*(D-1)**2*b_domain/2
    need(qd==Q(3609,2000000) and b_domain==Q(52162547172519000000,1996391), 'reviewed feedback norm arithmetic')
    need(e_bulk+e_robin==E_TRACE, 'reviewed full feedback allowance arithmetic')
    w_interval = (s_interval[0]-E_TRACE,s_interval[1]+E_TRACE)
    post = authenticate_inputs()
    need(post==authentication,'all input bytes unchanged after numeric evaluation')
    need(own_pin==sha(Path(__file__).read_bytes()),'producer unchanged during evaluation')
    return {'status':'CONDITIONAL_ACTUAL_BOUNDARY_AUXILIARY_MEAN_ENCLOSURE',
        'producer_sha256':own_pin,'input_authentication_before_numeric_import':authentication,
        'input_authentication_after_evaluation':post,'input_bytes_unchanged':True,
        'domain':{'a':str(A),'b':str(B),'d':str(D),'backward_closed':True,
            'shell_history':['0',str(D)],'standing_zero_wedge':'s+u<1; boundary s<1',
            'shell_history_exact_zero_interval':['0','1']},
        'partition':{'n':n,'h':str(A/n),'active_bulk_cells':len(active),'boundary_cells':len(shell_rows),
            'active_area':{k:str(v) for k,v in areas.items()},
            'zero_wedge_area':{k:str(v) for k,v in zero_areas.items()},
            'active_Robin_wave_kernel_weighted_area':'5/8',
            'boundary_positive_history_length':'1/4',
            'bbox_exterior_error_eligible_cells':sum(c.bbox_exterior_error_eligible for c,_ in active)},
        'actual_uncertainty':{'retained_actual_master_C3_error_budget':errors['master_error_budget_by_scaled_total_order'],
            'retained_actual_gauge_C3_error_caps':errors['actual_minus_gauge_comparison_mixed_scaled_derivative_caps'],
            'exterior_improvement_only_when_bbox_u_lo_ge_s_hi':True,
            'inherited_actual_uncertainty_in_source_intervals_once':True,
            'within_box_variation_and_outward_arithmetic_in_source_intervals':True,
            'separate_point_quadrature_error_required':False,
            'exact_source_map_terms_evaluated':{'metric_forcing_yy0':len(F_terms),'QC_y0':len(QC_terms)},
            'all23_feedback_couplings_in_analytic_allowance':True},
        'actual_bulk_source_cell_enclosures':bulk_rows,'actual_boundary_source_interval_enclosures':shell_rows,
        'signed_integral_ledger':{'bulk_Dminus_weight_1':interval_strings(totals['Dminus']),
            'bulk_strip_weight_1_over_2':interval_strings(totals['strip']),
            'KF_mean':interval_strings(bulk),'G_history_integral_mean':interval_strings(shell_G_integral),
            'LG_negative_history_mean':interval_strings(lg),'S_equals_KF_plus_LG_mean':interval_strings(s_interval)},
        'full_feedback_ledger':{'A_operator_bound':'1/10000','B_operator_bound':'1/10000',
            'D0_operator_bound':'1/1000000','sigma':str(SIGMA),'M_operator_bound':str(M),
            'F_full23_H_bound':str(GLOBAL_F),'G_full23_H_bound':str(GLOBAL_G),
            'q_D':str(qd),'W_X_domain_bound':str(b_domain),'bulk_value_allowance':str(e_bulk),
            'Robin_trace_value_allowance':str(e_robin),'total_point_mean_allowance':str(E_TRACE),
            'feedback_scope':'all23 bulk couplings and material Robin feedback; analytic bound, no evaluated feedback waveform',
            'norm_conversion':'yy probability mean has H-dual norm 1; no extra factor 3 or 17'},
        'actual_Wyy_probability_mean_interval':interval_strings(w_interval),
        'strict_positive':w_interval[0]>0,'strict_negative':w_interval[1]<0,
        'contains_zero':w_interval[0]<=0<=w_interval[1],
        'normalization':{'bulk':'retained metric_forcing=-2 ell^2 rho^2 QEred; no extra factorial two',
            'folded_kernel':'Robin yy: weight 1 on Dminus and 1/2 on Dplus minus Dminus',
            'shell':'Gyy=2 ell rho [-QK00_0+3 QKab_0+QC_y_0]; one-copy outward QK and varying-metric-lowered QC',
            'boundary_lift':'LG=-integral_0^d G(s) ds'},
        'limits':['Conditional on retained actual first-solution and coefficient premises.',
            'Auxiliary forced increment W at one point, not an evaluated global waveform or standalone physical total response.',
            'No physical force, cosmological-origin, or novelty claim.']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n',type=int,default=4)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--timing-output',type=Path,required=True)
    args = parser.parse_args()
    for path in (args.output,args.timing_output):
        need(path.is_absolute() and path.parent.resolve()==OUTPUT_ROOT.resolve() and not path.exists(),
             'new output directly in dedicated compute directory only')
    started = time.monotonic()
    result = run(args.n)
    runtime = time.monotonic()-started
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    timing = {'n':args.n,'python_optimized':bool(sys.flags.optimize),'actual_wall_seconds':runtime,
              'report_sha256':sha(args.output.read_bytes()),'producer_sha256':result['producer_sha256']}
    args.timing_output.write_text(json.dumps(timing,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'event':'completed','timing':timing,'S':result['signed_integral_ledger']['S_equals_KF_plus_LG_mean'],
                      'W':result['actual_Wyy_probability_mean_interval'],'contains_zero':result['contains_zero']},sort_keys=True),flush=True)


if __name__ == '__main__':
    main()
