#!/usr/bin/env python3
"""Own outward I/J and explicit all23 formula, only nominal absolute jets.

The immutable arithmetic helper is authenticated before import. No producer,
initial-core or Context method is imported or used here.
"""
import argparse
import hashlib
import importlib.util
import json
import math
import sys
from fractions import Fraction as Q
from pathlib import Path

sys.dont_write_bytecode = True
PACK = 'research/HDBLAST_EXACT_PACKET_INITIAL_20261031'


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def enc(x):
    if isinstance(x, Q): return str(x)
    if isinstance(x, dict): return {k: enc(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [enc(v) for v in x]
    return x


def authenticate(root, pins):
    need(root.is_dir() and not root.is_symlink(), 'SOURCE_ROOT_TYPE')
    out = {}
    for rel, pin in pins.items():
        name = Path(rel)
        need(not name.is_absolute() and '..' not in name.parts, 'SOURCE_PATH')
        p = root/name
        need(p.is_file(), 'SOURCE_FILE:'+rel)
        for pp in (p, *p.parents):
            if pp == root: break
            need(not pp.is_symlink(), 'SOURCE_SYMLINK:'+rel)
        need(sha(p) == pin, 'SOURCE_HASH:'+rel)
        out[rel] = pin
    return out


def campaign(m, coefficients, baseline, bits):
    m.BITS = bits
    I, J = m.I, m.J
    ell = m.ELL
    weights = m.GRAM
    limits = [Q(baseline['candidate_initial_absolute_H_caps'][k]) for k in ('u','u_r','u_rr')]
    limits[2] /= 2

    def norm(row):
        if all(x == 0 for x in row): return Q(0)
        return m.sqrt_bounds([I(x) for x in row])[1]

    def evaluate(lo, hi):
        if hi <= 1:
            return [[I(0) for _ in range(23)] for _ in range(3)]
        r = I(lo, hi)
        rz = J([r, 1/ell]+[0]*4)
        z = ell*rz
        bg = {k:m.poly(cs,z) for k,cs in coefficients.items()}
        if lo == 1:
            s = 1/((hi-1)*(2-hi))
            need(s >= Q(4096,63), 'TAIL_DOMAIN')
            ex = m.expminus(s-4)
            qjet = []
            for n in range(6):
                need(s >= 2*n, 'TAIL_MONOTONICITY')
                cap = s**(2*n)*ex.hi
                qjet.append((I(0,cap) if n == 0 else I(-cap,cap))/(ell**n*math.factorial(n)))
            q = J(qjet)
        else:
            den = (rz-1)*(2-rz)
            den.c[0] = Q(1,4)-(r-Q(3,2))**2
            q = m.expjetminus(1/den-4)
        e = m.expjetminus(ell*rz)
        displacement, unused_velocity = m.registered(bg,q,e)
        need(len(displacement) == 23 and all(v.n >= 2 for v in displacement), 'ALL23_JET_ORDER')
        return [[math.factorial(j)*ell**j*v.c[j] for v in displacement] for j in range(3)]

    seeds = [(Q(0),Q(1)),(Q(1),Q(65,64))]+[(Q(k,64),Q(k+1,64)) for k in range(65,112)]
    stack = [(lo,hi,0,i,None,None) for i,(lo,hi) in reversed(list(enumerate(seeds)))]
    leaves, splits = [], []
    while stack:
        lo,hi,depth,seed,parent,parent_id = stack.pop()
        jets = evaluate(lo,hi)
        raw = [[v.mag() for v in row] for row in jets]
        caps = raw if parent is None else [[min(x,y) for x,y in zip(row,old)] for row,old in zip(raw,parent)]
        hh = [norm(row) for row in caps]
        missed = any(x > y for x,y in zip(hh, limits))
        if missed and depth < 16 and len(seeds)+len(splits) < 4096:
            mid = (lo+hi)/2
            split_id = len(splits)
            splits.append({'r':[lo,hi], 'depth':depth, 'seed':seed, 'parent_split':parent_id,
                           'midpoint':mid, 'raw_component_caps':raw,'inherited_component_caps':caps})
            stack.append((mid,hi,depth+1,seed,caps,split_id))
            stack.append((lo,mid,depth+1,seed,caps,split_id))
        else:
            leaves.append({'r':[lo,hi], 'depth':depth, 'seed':seed,'parent_split':parent_id,
                           'raw_component_caps':raw,'component_caps':caps,'absolute_H_caps':hh,
                           'threshold_satisfied':not missed})
        if (len(leaves)+len(splits))%16 == 0:
            print(json.dumps({'bits':bits,'leaves':len(leaves),'splits':len(splits)}),flush=True)
    need(leaves[0]['r'][0] == 0 and leaves[-1]['r'][1] == Q(7,4) and
         all(a['r'][1] == b['r'][0] for a,b in zip(leaves,leaves[1:])), 'FULL_CLOSED_COVERAGE')
    need(len(leaves) == len(seeds)+len(splits) <= 4096 and max(v['depth'] for v in leaves) <= 16, 'FROZEN_BUDGET')
    caps = [max(row['absolute_H_caps'][j] for row in leaves) for j in range(3)]
    endpoint = evaluate(Q(7,4), Q(7,4))
    den = {}
    z = ell*J([I(0,Q(7,4)),1/ell]+[0]*4)
    for k in ('rho','p'):
        v = m.poly(coefficients[k], z).c[0]
        need(v.lo > 0 or v.hi < 0, 'NOMINAL_UNIFORM_DENOMINATOR:'+k)
        den[k] = v.row()
    return {'bits':bits, 'labels':m.LABELS,'Gram':weights,'seeds':seeds,'leaves':leaves,'splits':splits,
            'leaf_count':len(leaves),'split_count':len(splits),'nominal_denominators':den,
            'complete_absolute_H_caps':dict(zip(('u','u_r','u_rr'),caps)),
            'baseline_absolute_H_caps':{k:baseline['candidate_initial_absolute_H_caps'][k] for k in ('u','u_r','u_rr')},
            'strict_no_worse':all(x <= Q(baseline['candidate_initial_absolute_H_caps'][k]) for x,k in zip(caps,('u','u_r','u_rr'))),
            'rr_factor_two_gate':caps[2] <= limits[2],
            'all_leaf_thresholds_satisfied':all(v['threshold_satisfied'] for v in leaves),
            'endpoint_rows':[{'label':label,'Gram':g,'interior_closed_jet':[endpoint[j][i].row() for j in range(3)],
                              'quadratic_exterior_minus_interior_trace':[0,0,0]}
                             for i,(label,g) in enumerate(zip(m.LABELS,weights))],
            'exact_zero_region':['0','1'],'flat_endpoint_all_jets_zero':True,
            'outer_extension':'For each completed component at R=7/4 use its exact value, r and rr traces in quadratic continuation. Closed physical footprint only; no C3 claim.',
            'tail_certificate':m.tail_proof()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--source-root',required=True,type=Path)
    ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--freeze-sha256',required=True)
    args = ap.parse_args()
    here = Path(__file__).resolve().parent.parent
    freeze_path = here/'FREEZE.json'
    need(sha(freeze_path) == args.freeze_sha256, 'FROZEN_CONFIG_HASH')
    frozen = json.loads(freeze_path.read_text())
    for rel, pin in frozen['own_sha256'].items():
        need(sha(here/rel) == pin, 'OWN_FILE_HASH:'+rel)
    pins_path = here/'source-pins.json'
    need(sha(pins_path) == frozen['source_map_sha256'], 'SOURCE_MAP_HASH')
    pins = json.loads(pins_path.read_text())
    before = authenticate(args.source_root,pins)
    need(not args.output.exists(), 'FRESH_OUTPUT_REQUIRED')
    helper = here/'scripts/pinned_independent_exact.py'
    spec = importlib.util.spec_from_file_location('authenticated_own_independent_arithmetic',helper)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    bg = json.loads((args.source_root/'research/HDBLAST_INPUT_PRECISION_20261023/initial/RECENTERED_CANDIDATE_NORMAL.json').read_text())
    cs = {k:list(map(Q,v)) for k,v in bg['representation_background_coefficients_ascending'].items()}
    old = json.loads((args.source_root/PACK/'reference/INITIAL.json').read_text())
    runs = [campaign(m,cs,old['campaigns'][0],128)]
    if not (runs[0]['rr_factor_two_gate'] and runs[0]['strict_no_worse']):
        runs.append(campaign(m,cs,old['campaigns'][0],192))
    after = authenticate(args.source_root,pins)
    need(after == before and sha(pins_path) == frozen['source_map_sha256'], 'SOURCE_POSTAUTH')
    need(sha(freeze_path) == args.freeze_sha256, 'FREEZE_UNCHANGED')
    for rel,pin in frozen['own_sha256'].items():
        need(sha(here/rel) == pin, 'OWN_POSTAUTH:'+rel)
    report = {'schema':'HDBLAST_SAME_CANDIDATE_INDEPENDENT_ABSOLUTE_V1','candidate':'W1_exact_packet',
              'status':'PASS_INDEPENDENT_ABSOLUTE_REFINEMENT' if all(x['rr_factor_two_gate'] and x['strict_no_worse'] for x in runs) else 'INCONCLUSIVE',
              'freeze_sha256':args.freeze_sha256,'selected_source_map_sha256':frozen['source_map_sha256'],
              'selected_source_count':len(before),'source_preauth_before_helper_import':True,'source_postauth_equal':True,
              'own_sha256':frozen['own_sha256'],'campaigns':runs,
              'initial_error_velocity_target_source_operator_changed':False,
              'formula':'Own immutable outward I/J arithmetic and explicit all23 registered formulas; nominal polynomial background, exact q/E analytic intervals, completed components before continuation.',
              'scope':'Conditional absolute nominal initial jets only. Lower upper certificates do not establish smaller true error. Existing conditional exact initial/precontact results retain hypotheses; nonlinear shell contact separate.'}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(enc(report),sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'caps':[{k:float(v) for k,v in r['complete_absolute_H_caps'].items()} for r in runs]},sort_keys=True))


if __name__ == '__main__': main()
