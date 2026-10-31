#!/usr/bin/env python3
"""Independent exact-substitution review of the actual third-jet adapter.

The oracle substitutes x=x0+t/ell and a nonaffine polynomial r=r(y)
directly into generic raw polynomials. It does not differentiate the
producer's radial conversion or replay its Taylor-jet arithmetic.
"""
import argparse
import ast
from fractions import Fraction as Q
from functools import reduce
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path

GEOM = 'HDBLAST_GEOMETRIC_BULK_SOURCES_20261016'
MANIFEST_PIN = '3ff1f95c151d2da6a1ab5e9b54718c2bb325c10d9c62a7ed9f97ef2bb36f0bb2'
SOURCE_PINS = {
 'HDBLAST_CAUSAL_SOURCE_REGULARITY_20261030/reference/MASTER4.json': 'd04f67e2e1863cca423a4b0835bef0a739bf5c471cf1d3e2f6730677d970490e',
 'HDBLAST_CAUSAL_SOURCE_REGULARITY_20261030/checks/check_master_fourth.py': '6962e8c9bdca416ad072029df1d84ae666467e11656f0078f3e5ba4033659ed0',
 'HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/reference/GAUGE4.json': '94579087d39d9c849e9264d40bd1d71d64892a668fe3e2e08aade9ca4296b873',
 'HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/reference/FX.json': '63ed8e3a7c0544db93de0d3bcbdba4652eb4d9a6bcd8f434327e785cd9986e60',
 'HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/checks/check_actual_gauge_fourth.py': '9d88529967c6da55684779a67c828aa15062052ccd8307e188c7211963cf305b',
 'HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/checks/check_original_bulk_fx.py': '30def8af98798c18a5f48d73294945da5522d5c3fb11269b1d8226cde70676ee',
}
IMPORT_NAMES = {'argparse', 'hashlib', 'importlib.util', 'json', 'math', 'fractions', 'pathlib'}
EXECUTED = ('numeric/actual_first_jets.py', 'numeric/interval_core.py',
            'numeric/exterior_error_caps.py', 'numeric/check_inherited_errors.py')

def need(ok, message):
    if not ok:
        raise ValueError(message)

def sha(raw):
    return hashlib.sha256(raw).hexdigest()

def source_snapshot(repo):
    research = repo / 'research'
    cap = research / GEOM
    raw = (cap / 'MANIFEST.json').read_bytes()
    need(sha(raw) == MANIFEST_PIN, 'independent geometric manifest pin')
    out = {GEOM + '/MANIFEST.json': sha(raw)}
    seen = {'MANIFEST.json'}
    for row in json.loads(raw)['files']:
        rel = Path(row['path'])
        need(not rel.is_absolute() and '..' not in rel.parts and str(rel) not in seen,
             'independent safe manifest member')
        seen.add(str(rel))
        p = cap / rel
        need(p.is_file() and not p.is_symlink(), 'independent regular member')
        data = p.read_bytes()
        need(len(data) == row['bytes'] and sha(data) == row['sha256'],
             'independent geometric member authentication: ' + str(rel))
        out[GEOM + '/' + str(rel)] = sha(data)
    need({str(p.relative_to(cap)) for p in cap.rglob('*') if p.is_file()} == seen,
         'independent closed source inventory')
    need(not any(p.is_symlink() for p in cap.rglob('*')), 'independent no source symlinks')
    for rel, pin in SOURCE_PINS.items():
        p = research / rel
        need(p.is_file() and not p.is_symlink() and sha(p.read_bytes()) == pin,
             'independent selected source authentication: ' + rel)
        out[rel] = pin
    return out

def import_audit(path):
    tree = ast.parse(path.read_text())
    imports = []
    dynamic = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.extend(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            imports.append(node.module)
        elif isinstance(node, ast.Call):
            name = ast.unparse(node.func)
            need(name not in {'exec', 'eval', '__import__', 'compile'},
                 'no unchecked executable evaluation')
            if name.endswith('spec_from_file_location'):
                dynamic.append(ast.unparse(node.args[1]))
    need(set(imports) <= IMPORT_NAMES, 'permitted standard-library imports')
    return {'imports': sorted(set(imports)), 'dynamic_paths': sorted(dynamic)}

# Sparse direct bivariate polynomial algebra, truncated only at physical degree3.
def plus(p, q):
    out = dict(p)
    for key, value in q.items():
        out[key] = out.get(key, Q(0)) + value
    return {key: value for key, value in out.items() if value}

def times(p, q):
    out = {}
    for (a, b), v in p.items():
        for (c, d), w in q.items():
            key = a+c, b+d
            if sum(key) <= 3:
                out[key] = out.get(key, Q(0)) + v*w
    return {key: value for key, value in out.items() if value}

def power(p, n):
    out = {(0, 0): Q(1)}
    for unused in range(n):
        out = times(out, p)
    return out

def inverse(p):
    """Finite geometric inverse, independent of producer coefficient recursion."""
    base = p[(0,0)]
    h = {key:-value/base for key,value in p.items() if key != (0,0)}
    out = {}
    for n in range(4):
        out = plus(out,{key:value/base for key,value in power(h,n).items()})
    return out

def scaled(value,p):
    return {key:value*c for key,c in p.items() if value*c}

def derivative(p,axis):
    out = {}
    for key,value in p.items():
        if key[axis]:
            k = list(key); k[axis] -= 1
            out[tuple(k)] = key[axis]*value
    return out

def substitute(raw, xpoly, rpoly):
    out = {}
    for (i, j), value in raw.items():
        term = times(power(xpoly, i), power(rpoly, j))
        out = plus(out, {key: value*c for key, c in term.items()})
    return out

class RawJet:
    def __init__(self, polynomial, x0, r0):
        self.p = polynomial
        self.x0, self.r0 = x0, r0
    def derivative(self, i, j):
        return sum((value*math.prod(range(m-i+1, m+1))*math.prod(range(n-j+1, n+1))
                    *self.x0**(m-i)*self.r0**(n-j)
                    for (m, n), value in self.p.items() if m >= i and n >= j), Q(0))
    def val(self):
        return self.derivative(0, 0)

class Value:
    def __init__(self, value): self.value = value
    def val(self): return self.value

def exact_chain_checks(module):
    rows = []
    # Each nonaffine r(y) has nonzero second and third derivatives. The
    # local rho(y)=-1/(ell*r'(y)) is positive near0; its logarithmic
    # derivatives give a,ay without imposing any producer recurrence.
    fixtures = (
      (Q(2,3), Q(1,3), Q(2,5), Q(-2), Q(3), Q(5)),
      (Q(3,5), Q(-2,7), Q(4,9), Q(-3,2), Q(-7,3), Q(11,5)),
      (Q(5,4), Q(3,8), Q(-1,6), Q(-5,3), Q(2,7), Q(-13,4)),
    )
    for index, (ell, x0, r0, r1, r2, r3) in enumerate(fixtures):
        raw = {(i,j): Q((-1)**(i+j+index)*(2+7*i+11*j+3*i*j), 1+i+j)
               for i in range(6) for j in range(6-i)}
        physical = substitute(raw, {(0,0):x0, (0,1):1/ell},
                              {(0,0):r0, (1,0):r1, (2,0):r2/2, (3,0):r3/6})
        rho = -1/(ell*r1)
        a = -r2/r1
        ay = (r2*r2-r3*r1)/(r1*r1)
        actual = module.physical_through_three(RawJet(raw,x0,r0),Value(rho),Value(a),ay,ell)
        values = {}
        for n in range(4):
            for ny in range(n+1):
                nt = n-ny
                suffix = 'y'*ny+'t'*nt
                expected = physical.get((ny,nt), Q(0))*math.factorial(ny)*math.factorial(nt)
                need(actual[suffix] == expected,
                     'direct substitution witness '+str(index)+' '+suffix+': residual '+str(actual[suffix]-expected))
                values[suffix] = str(expected)
                # Direct partial differentiation of the polynomial in all
                # word orders independently certifies the commuting ordinary
                # coordinate slots used by the adapter.
                for word in set(itertools.permutations(suffix)):
                    p = dict(physical)
                    for letter in word:
                        axis = 0 if letter == 'y' else 1
                        changed = {}
                        for key, value in p.items():
                            if key[axis]:
                                k = list(key); k[axis] -= 1
                                changed[tuple(k)] = value*key[axis]
                        p = changed
                    need(p.get((0,0),Q(0)) == expected, 'ordinary derivative commutation')
        rows.append({'fixture':index,'physical_derivatives':values})
    return rows

def serialize(values):
    return {key: value.row() if hasattr(value,'row') else str(value)
            for key, value in sorted(values.items())}

def raw_composition_check(adapter):
    """Independent sparse products/geometric inverses of the registered Eq7."""
    ell = adapter.ell
    names = ('rho','p','a','R','B','k','c','H','vartheta','W','f','b','d')
    polynomials = {}
    for n,name in enumerate(names):
        degree = 4 if name in ('W','f','b','d') else 3
        # Nonzero entries in all mixed directions exercise product cross terms.
        polynomials[name] = {(i,j):Q((-1)**(i+j+n)*(n+2+3*i+5*j),2**(i+j+3))
                             for i in range(degree+1) for j in range(degree+1-i)}
    polynomials['rho'][(0,0)] = Q(2)
    polynomials['c'][(0,0)] = Q(3,2)
    g = polynomials
    add = lambda *ps: reduce(plus,ps,{})
    mul = lambda *ps: reduce(times,ps,{(0,0):Q(1)})
    invrho,invc = inverse(g['rho']),inverse(g['c'])
    u = mul(g['R'],g['vartheta'],g['W'])
    v = mul(g['B'],g['vartheta'],add(scaled(1/ell,derivative(g['W'],1)),mul(g['k'],g['W'])))
    expected = {
      'A':add(scaled(-4,u),scaled(-2/ell,mul(derivative(g['f'],1),invrho))),
      'B':add(scaled(1/ell,mul(derivative(g['f'],0),invrho)),scaled(1/ell,derivative(g['b'],1))),
      'C':add(mul(g['f'],invrho,invc),scaled(-1/ell,mul(g['c'],derivative(g['d'],1)))),
      'D':scaled(-2,add(u,mul(g['a'],g['f']),scaled(1/ell,derivative(g['b'],0)))),
      'E':add(scaled(1/ell,mul(g['c'],derivative(g['d'],0))),scaled(-1,mul(g['b'],invc))),
      'L':scaled(2,add(u,mul(g['a'],g['f']),mul(g['H'],g['b']))),
      'M':scaled(2,g['d']),
      's':add(v,mul(g['p'],g['f'])),
    }
    jets = {name:adapter.T(p) for name,p in g.items()}
    metric,scalar = adapter.old.frame_amplitudes(
      jets['W'],[jets[name] for name in ('f','b','d')],
      *[jets[name] for name in ('rho','p','a','R','B','k','c','H','vartheta')])
    actual = {**metric,'s':scalar}
    checks = 0
    for name,p in expected.items():
        for i in range(4):
            for j in range(4-i):
                exact = p.get((i,j),Q(0))*math.factorial(i)*math.factorial(j)
                value = actual[name].derivative(i,j)
                need(value.lo <= exact <= value.hi,
                     'independent composite derivative '+name+' '+str((i,j)))
                checks += 1
    return {'raw_composite_derivatives_enclosed':checks,
            'raw_gauge_and_master_degree':4,'coefficient_degree':3,
            'inverse_oracle':'finite geometric series, without producer inverse recursion'}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo-root', type=Path, required=True)
    parser.add_argument('--adapter', type=Path, required=True)
    parser.add_argument('--adapter-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    need(not args.output.exists(), 'fresh output required')
    initial = source_snapshot(args.repo_root)
    need(sha(args.adapter.read_bytes()) == args.adapter_sha256, 'adapter pinned before first import')
    audits = {'adapter': import_audit(args.adapter)}
    for rel in EXECUTED:
        audits[rel] = import_audit(args.repo_root/'research'/GEOM/rel)
    spec = importlib.util.spec_from_file_location('independently_pinned_actual_third', args.adapter)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    exact = exact_chain_checks(module)
    adapter = module.Adapter(args.repo_root)
    composite = raw_composition_check(adapter)
    nr = args.repo_root/'research'/GEOM/'numeric/inputs'
    master = json.loads((nr/'HIGHER_JETS_NORMAL.json').read_text())['exact']['interior']
    m4 = json.loads((args.repo_root/'research/HDBLAST_CAUSAL_SOURCE_REGULARITY_20261030/reference/MASTER4.json').read_text())
    g4 = json.loads((args.repo_root/'research/HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/reference/GAUGE4.json').read_text())
    expected_master = [Q(master['U']),Q(master['first']),max(map(Q,master['second'].values())),
                       max(map(Q,master['third'].values())),max(map(Q,m4['exact']['all_scaled_order_four_mixed_caps'].values()))]
    expected_gauge = {tuple(map(int,key.split(','))):list(map(Q,value))
                      for key,value in g4['actual_mixed_scaled_derivative_caps'].items()}
    need(adapter.master_caps == expected_master, 'faithful authenticated master0through4 cap join')
    need(adapter.gauge_caps == expected_gauge, 'faithful authenticated gauge0through4 cap join')
    I = adapter.I
    # Cell, shell one-sided trace, shell time interval, exterior, and final
    # characteristic point all exercise the actual retained adapter.
    cases = (
      ('interior',I(Q(3,2),Q(3,2)+Q(1,4096)),I(Q(1,4),Q(1,4)+Q(1,4096)),False),
      ('shell_point',I(Q(3,2)),I(0),False),
      ('shell_segment',I(1,2),I(0),False),
      ('exterior',I(Q(1,4)),I(Q(5,4)),False),
      ('causal_vertex',I(3),I(5),False),
      ('clipped_causal',I(2,3),I(5,6),True),
    )
    case_rows = []
    expected_metric = {name+suffix for name in 'ABCDELM'
                       for suffix in ('','y','t','yy','yt','tt','yyy','yyt','ytt','ttt')}
    expected_scalar = {'s'+suffix for suffix in ('','y','t','yy','yt','tt','yyy','yyt','ytt','ttt')}
    for label,x,r,clip in cases:
        m,s,b = adapter.actual_jet(x,r,clip_causal=clip)
        oldm,olds,oldb = adapter.old.actual_jet(x,r,clip_causal=clip)
        need(set(m) == expected_metric and set(s) == expected_scalar, 'complete actual coherent interface')
        need(serialize({key:m[key] for key in oldm}) == serialize(oldm), 'retained actual metric lower2 identical')
        need(serialize({key:s[key] for key in olds}) == serialize(olds), 'retained actual scalar lower2 identical')
        need(serialize(b) == serialize(oldb), 'retained same actual background')
        need(all(v.lo <= 0 <= v.hi for k,v in {**m,**s}.items()
                 if k.endswith(('yyy','yyt','ytt','ttt'))), 'uniform actual magnitude enclosures')
        case_rows.append({'label':label,'metric_slots':len(m),'scalar_slots':len(s),
                          'retained_lower2_identical':True})
    rejected = []
    for label,x,r in (('negative_time',I(-1),I(1)), ('excess_time',I(4),I(1)),
                      ('negative_depth',I(1),I(-1)), ('outside_triangle',I(3),I(6)),
                      ('unclipped_crossing',I(2,3),I(5,6))):
        try:
            adapter.actual_jet(x,r)
        except ValueError:
            rejected.append(label)
        else:
            raise ValueError('invalid domain accepted: '+label)
    sentinel, coeff = adapter._global_third(17)
    need(serialize(sentinel) == serialize(adapter.global_third), 'unused coefficient4 sentinel invariance')
    need(coeff == adapter.coefficient_boxes, 'coefficient through3 sentinel invariance')
    need(len(sentinel) == 32 and all(len(v) == 10 for v in coeff.values()), 'separate coefficient3 ledger')
    final = source_snapshot(args.repo_root)
    need(final == initial, 'full source closure unchanged after execution')
    need(sha(args.adapter.read_bytes()) == args.adapter_sha256, 'adapter unchanged after execution')
    out = {'status':'PASS_INDEPENDENT_ACTUAL_THIRD_JET_REVIEW',
           'adapter_sha256':args.adapter_sha256,
           'geometric_manifest_sha256':MANIFEST_PIN,
           'geometric_members':len(initial)-len(SOURCE_PINS)-1,
           'selected_higher_jet_source_sha256':SOURCE_PINS,
           'import_audit':audits,'direct_substitution_oracles':exact,
           'independent_raw_composition_oracle':composite,
           'actual_cases':case_rows,'invalid_domain_rejected':rejected,
           'uniform_third_jets':serialize(sentinel),
           'coefficient4_sentinel_invariant':True,
           'scope':{'actual_lower2_retained':True,'actual_uniform_third_magnitude':True,
                    'physical_one_sided_shell':True,'differentiate_selected_intervals':False,
                    'numerical_third_derivative_error_certified':False,
                    'higher_PDE_smoothness_and_cap_proofs_reestablished':False}}
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(out['status'])

if __name__ == '__main__':
    main()
