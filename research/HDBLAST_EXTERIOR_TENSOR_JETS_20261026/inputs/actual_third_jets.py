#!/usr/bin/env python3
"""Certified full TT forcing and its scaled source-coordinate radial derivative.

The source's only physical second-derivative inputs are Myy, Myt, Mtt.
Since M=2d, Fr requires gauge C3 and master C2, already certified for the
actual first solution.  No derivative of an independently chosen interval
is taken: Taylor coefficients are actual smooth-solution derivative bounds.
"""
import argparse
import hashlib
import importlib.util
import json
import math
import sys
from fractions import Fraction as F
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
TT_REL = Path('research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026')
TT_MANIFEST = 'a35f84119e84bf2ea8a712ded86ff79977d3ee91ac25f91ccbd202d379c859ac'
TT_SOURCE = 'bdcae6445126027611a4ff8e3677ec9c6d05743d1a73bd4f932aa58f42b0b9ca'
EXPECTED_SYMBOLS = ('A','At','Ay','B','Bt','By','C','Ct','Cy','D','Dt','Dy',
                    'E','Et','Ey','L','Lt','Ly','M','Mt','Mtt','My','Myt','Myy',
                    'Vp','a','ay','c','ell','rho','s','theta')
DECLARATION = {'source_point_probes': 3, 'source_cell_probes': 2,
               'mesh_search': False, 'arithmetic': 'exact rationals, outward grid 2^112',
               'acceptance': 'authenticated C3 derivative errors and full 78-term algebraic AD'}
PROTECTED_LABELS = ('HDBLAST_COUPLED_MEAN_PILOT_20261025',
    'HDBLAST_INPUT_PRECISION_20261023', 'HDBLAST_COMMON_GAUGE_FORCING_20261016',
    'HDBLAST_NUMERICAL_RESPONSE_INPUTS_20261016', 'HDBLAST_GEOMETRIC_BULK_SOURCES_20261016',
    'HDBLAST_FORCED_RESPONSE_BOUNDS_20261023', 'HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023',
    'HDBLAST_FULL_ANGULAR_RESPONSE_20261026', 'HDBLAST_SHELL_FEEDBACK_PRECISION_20261025')


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def guard_output(repo_root, output):
    """Reject protected destinations before authentication or evaluation work."""
    output = Path(output)
    need(not output.is_symlink(), 'OUTPUT_SYMLINK')
    need(not output.exists(), 'FRESH_OUTPUT_REQUIRED')
    target = output.resolve()
    research = Path(repo_root).resolve()/'research'
    for label in PROTECTED_LABELS:
        root = (research/label).resolve()
        need(target != root and root not in target.parents, 'OUTPUT_INSIDE_FROZEN_CAPSULE:'+label)
    # Scratch inputs remain external.  Once packaged under research/<capsule>,
    # protect that complete capsule even when the caller supplies another repo.
    for parent in Path(__file__).resolve().parents:
        if parent.parent.name == 'research':
            need(target != parent and parent not in target.parents, 'OUTPUT_INSIDE_PACKAGED_CAPSULE')
            break


def authenticate_tt(repo):
    """Authenticate the complete extra frozen TT closure before importing it."""
    root = Path(repo).resolve()/TT_REL
    need(not root.is_symlink(), 'TT_ROOT_SYMLINK')
    manifest = root/'MANIFEST.json'
    need(sha(manifest) == TT_MANIFEST, 'TT_MANIFEST_PIN')
    seen = {'MANIFEST.json'}
    for row in json.loads(manifest.read_text())['files']:
        rel = Path(row['path'])
        need(not rel.is_absolute() and '..' not in rel.parts and str(rel) not in seen, 'TT_MEMBER_PATH')
        seen.add(str(rel))
        path = root/rel
        need(path.is_file() and not path.is_symlink(), 'TT_MEMBER_TYPE')
        need(path.stat().st_size == row['bytes'] and sha(path) == row['sha256'], 'TT_MEMBER_BYTES:'+str(rel))
    for path in root.rglob('*'):
        need(not path.is_symlink(), 'TT_TREE_SYMLINK')
        if path.is_file():
            need(str(path.relative_to(root)) in seen, 'TT_UNLISTED_FILE')
    source = root/'tensor/tt_response.py'
    need(sha(source) == TT_SOURCE, 'TT_SOURCE_PIN')
    return {'manifest_sha256': TT_MANIFEST, 'files_including_manifest': len(seen)}


class Dual:
    """First source-r derivative in any exact outward interval scalar ring."""
    __slots__ = ('v', 'd')

    def __init__(self, value, derivative=0):
        self.v, self.d = value, derivative

    @staticmethod
    def cv(value):
        return value if isinstance(value, Dual) else Dual(value, 0)

    def __add__(self, other):
        other = Dual.cv(other)
        return Dual(self.v+other.v, self.d+other.d)
    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.v, -self.d)

    def __sub__(self, other):
        return self+-Dual.cv(other)

    def __rsub__(self, other):
        return Dual.cv(other)+-self

    def __mul__(self, other):
        other = Dual.cv(other)
        return Dual(self.v*other.v, self.d*other.v+self.v*other.d)
    __rmul__ = __mul__

    def inv(self):
        return Dual(1/self.v, -self.d/(self.v**2))

    def __truediv__(self, other):
        return self*Dual.cv(other).inv()

    def __rtruediv__(self, other):
        return Dual.cv(other)*self.inv()

    def __pow__(self, power):
        need(isinstance(power, int), 'INTEGER_POWER')
        if power == 0:
            return Dual(self.v**0, 0)
        return Dual(self.v**power, power*(self.v**(power-1))*self.d)


def derivative_order_audit(terms):
    symbols = tuple(sorted({name for term in terms for name in term['powers']}))
    need(len(terms) == 78 and symbols == tuple(sorted(EXPECTED_SYMBOLS)), 'FULL_TT_DEPENDENCY_SIGNATURE')
    rows = {}
    max_master = max_gauge = 0
    for name in symbols:
        if name[0] in 'ABCDELM' and name[1:] in ('','y','t','yy','yt','tt'):
            suffix = name[1:]
            physical = len(suffix)
            gauge_base = 0 if name[0] in 'LM' else 1
            # L depends on W without differentiation; M is independent of W.
            master_base = None if name[0] in 'BCEM' else 0
            gauge = gauge_base+physical+1
            master = None if master_base is None else master_base+physical+1
        elif name == 's':
            physical, gauge, master = 0, 1, 2
        else:
            continue
        need(gauge <= 3 and (master is None or master <= 2), 'C3_SUFFICIENCY:'+name)
        rows[name] = {'physical_derivative_order': physical,
                      'gauge_order_after_source_r': gauge,
                      'master_order_after_source_r': master}
        max_gauge = max(max_gauge, gauge)
        if master is not None:
            max_master = max(max_master, master)
    seconds = [term for term in terms if any(name in ('Myy','Myt','Mtt') for name in term['powers'])]
    need(len(seconds) == 3, 'THREE_TT_PRINCIPAL_TERMS')
    return {'full_terms': len(terms), 'symbols': list(symbols), 'per_dynamic_symbol': rows,
            'required_actual_master_total_order': max_master,
            'required_actual_gauge_total_order': max_gauge,
            'only_second_physical_derivative_symbols': ['Mtt','Myt','Myy'],
            'second_physical_derivative_terms': seconds,
            'C4_derivative_premise_required': False,
            'reason': 'A through E have reconstruction order at most1, L and M order0; source suffix order is at most1 except M order2; M=2d.'}


def actual_frame_jets(ctx, x, r):
    """Degree3 input jets; reconstructed A..E/scalar valid through degree2."""
    j = ctx.jets
    I, T = j.I, j.T
    x = I(x.lo, x.hi) if hasattr(x, 'lo') else I.cv(x)
    r = I(r.lo, r.hi) if hasattr(r, 'lo') else I.cv(r)
    need(0 <= x.lo <= x.hi <= 3 and 0 <= r.lo <= r.hi and x.hi+r.hi <= 8, 'ACTUAL_JET_DOMAIN')
    need(r.lo >= x.hi, 'BACKWARD_CLOSED_EXTERIOR_BOX')
    data, errors, pins = j.inputs()
    W, P = j.comparisons(x, r, data)
    U = [P, T(), T()]
    ew, eu = j.exterior_caps()
    need(len(ew) == 4, 'ACTUAL_C3_MASTER_CERTIFICATE')
    for i, k in j.KEYS:
        den = math.factorial(i)*math.factorial(k)
        W.c[i,k] += I(-ew[i+k], ew[i+k])/den
        for channel in range(3):
            cap = eu[f'{i},{k}'][channel]
            U[channel].c[i,k] += I(-cap, cap)/den
    rr, eta, pp, aa, hh, RR, BB, kk, _, _ = j.core.background(r, data, 'none')
    rho, p, a, R, B, k = map(j.radial, (rr, pp, aa, RR, BB, kk))
    c, H, vartheta = j.temporal(x)
    amps, scalar = j.frame_amplitudes(W, U, rho, p, a, R, B, k, c, H, vartheta)
    bg = {'rho':rho, 'c':c, 'theta':H, 'a':a,
          'ay':-1/rho**2-p**2/3, 'p':p,
          'py':p*(j.radial(hh)-4*a),
          'Vp':j.core.poly(j.core.UP[1], j.radial(eta)),
          'Vpp':j.core.poly(j.core.UP[2], j.radial(eta)),
          'Vppp':j.core.poly(j.core.UP[3], j.radial(eta)),
          'ell':T.cv(j.ELL)}
    return amps, scalar, bg


def physical_expression(amp, suffix, rho, a, ell):
    L = 1/(ell*rho)
    if suffix == '':
        return amp
    if suffix == 'y':
        return -L*amp.dr()
    if suffix == 't':
        return amp.dx()/ell
    if suffix == 'yy':
        return L**2*amp.dr().dr()+a*L*amp.dr()
    if suffix == 'yt':
        return -L*amp.dx().dr()/ell
    if suffix == 'tt':
        return amp.dx().dx()/ell**2
    raise ValueError('PHYSICAL_SUFFIX:'+suffix)


def source_variables(ctx, x, r):
    amps, scalar, bg = actual_frame_jets(ctx, x, r)
    names = sorted({name for term in ctx.terms[5] for name in term['powers']})
    raw = {}
    for name in names:
        if name in bg:
            expression = bg[name]
        elif name == 's':
            expression = scalar
        else:
            suffix = name[1:]
            reconstruction = 0 if name[0] in 'LM' else 1
            need(reconstruction+len(suffix)+1 <= 3, 'NO_UNCERTIFIED_TAYLOR_COEFFICIENT:'+name)
            expression = physical_expression(amps[name[0]], suffix, bg['rho'], bg['a'], ctx.ell)
        raw[name] = Dual(expression.val(), expression.dr().val())
    return raw


def evaluate(ctx, x, r, audit=False):
    """Return {'F','Fr'} on a whole exterior rectangle, with all TT terms."""
    raw = source_variables(ctx, x, r)
    converted = {name: Dual(ctx.dyadic.D(value.v), ctx.dyadic.D(value.d)) for name,value in raw.items()}
    out = ctx.coeff.evaluate_encoded(ctx.terms[5], converted)
    need(isinstance(out, Dual), 'FULL_SOURCE_RADIAL_AD')
    if audit:
        reference = ctx.coeff.evaluate_encoded(ctx.terms[5], raw)
        need(out.v.row() == reference.v.row() and out.d.row() == reference.d.row(), 'DUAL_DYADIC_REFERENCE_ENDPOINT_IDENTITY')
        m, s, b = ctx.jets.actual_jet(x, r, mode='actual')
        old_raw = {**m, **s, **b}
        for name in raw:
            expected = old_raw[name]
            expected = expected.row() if hasattr(expected, 'row') else ctx.I(expected).row()
            need(raw[name].v.row() == expected, 'VALUE_ORACLE_ENDPOINT_IDENTITY:'+name)
        old_source = ctx.coeff.evaluate_encoded(ctx.terms[5], old_raw)
        need(out.v.row() == old_source.row(), 'ACTUAL_FULL_SOURCE_ENDPOINT_IDENTITY')
    return {'F': out.v, 'Fr': out.d}


class Context:
    def __init__(self, repo_root):
        self.repo = Path(repo_root).resolve()
        self.tt_before = authenticate_tt(self.repo)
        old = load('authenticated_tt_source_for_radial_jets', self.repo/TT_REL/'tensor/tt_response.py')
        self.base = old.Context(self.repo)
        for name in ('I','ell','jets','data','dyadic','coeff','terms','maps','maps_path','roots','before'):
            setattr(self, name, getattr(self.base, name))
        self.calls = 0
        self.dependency_audit = derivative_order_audit(self.terms[5])
        self.certificate = {'schema':'HDBLAST_ACTUAL_FULL_TT_SOURCE_RADIAL_DERIVATIVE_V1',
            'status':'CERTIFIED_ACTUAL_FULL_TT_SOURCE_AND_RADIAL_DERIVATIVE_ORACLE',
            'source_module_sha256':sha(Path(__file__)),
            'frozen_TT_authentication':self.tt_before, 'seven_predecessor_authentication':self.before,
            'source_map_sha256':sha(self.maps_path),
            'actual_C3_error_sha256':sha(self.roots['G']/'numeric/inputs/APPROXIMATION_ERRORS_NORMAL.json'),
            'actual_C3_independent_review_sha256':sha(self.roots['G']/'numeric/inputs/C3_ERROR_FINAL_NORMAL.json'),
            'exterior_C3_algebra_sha256':sha(self.roots['G']/'numeric/exterior_error_caps.py'),
            'dependency_order_audit':self.dependency_audit,
            'budget_declaration':DECLARATION,
            'domain':'0<=x<=3, r>=x, x+r<=8; oracle requires whole rectangles r.lo>=x.hi',
            'requested_cone':'0<=x<=1/4, 1+x<=r<=3/2-x',
            'radial_coordinate':'Fr=partial_r F at fixed scaled time x; partial_r=-ell*rho*partial_y',
            'proof':'Actual normalized Taylor coefficients are comparison derivatives plus authenticated differentiated C3 actual errors. Algebraic product, inverse and chain rules propagate them. The collected full TT source uses gauge derivatives through3 and master through2 after one radial derivative. The Dual ring differentiates all 78 terms and every used background variable. No C3 interval is differentiated and no actual error is omitted.',
            'premises':'Authenticated retained smooth actual background/master/gauge solution, differentiated integral/reflection identities, actual C3 error certificates and formal all-point source projection remain conditional scientific premises.',
            'all_actual_error_caps_retained':True,
            'C4_required':False, 'finite_difference_used':False,
            'arithmetic':{'outward_grid':str(self.jets.core.GRID), 'exact_rational_decisions':True}}

    def forcing_with_radial(self, x, r, audit=False):
        out = evaluate(self, x, r, audit)
        self.calls += 1
        return out['F'], out['Fr']

    def forcing(self, x, r, audit=False):
        return self.forcing_with_radial(x, r, audit)[0]

    def finish(self):
        self.base.finish()
        need(authenticate_tt(self.repo) == self.tt_before, 'FROZEN_TT_CLOSURE_CHANGED')
        need(sha(Path(__file__)) == self.certificate['source_module_sha256'], 'SOURCE_MODULE_CHANGED')


def encode(value):
    if hasattr(value, 'row'):
        return value.row()
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key:encode(item) for key,item in value.items()}
    if isinstance(value, (tuple,list)):
        return [encode(item) for item in value]
    return value


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    guard_output(args.repo_root, args.output)
    ctx = Context(args.repo_root)
    rows = []
    for label,x,r,w in [('initial_packet_point',F(0),F(5,4),None),
                        ('middle_cone_point',F(1,8),F(5,4),None),
                        ('apex_point',F(1,4),F(5,4),None),
                        ('initial_left_N4_cell',F(0),F(1),F(1,16)),
                        ('middle_N4_cell',F(1,8),F(19,16),F(1,16))]:
        xx = ctx.I(x) if w is None else ctx.I(x,x+w)
        rr = ctx.I(r) if w is None else ctx.I(r,r+w)
        f,fr = ctx.forcing_with_radial(xx, rr, audit=True)
        rows.append({'label':label,'x':xx,'r':rr,'F':f,'Fr':fr,
                     'display_only_F':[float(f.lo),float(f.hi)],
                     'display_only_Fr':[float(fr.lo),float(fr.hi)],
                     'full_source_and_fraction_dyadic_endpoint_audits':True})
        print(json.dumps({'probe':label,'Fr_display_only':[float(fr.lo),float(fr.hi)]}),flush=True)
    ctx.finish()
    report = {**ctx.certificate,'probes':rows,'oracle_calls':ctx.calls,
              'all_frozen_closures_reauthenticated_after_run':True}
    with args.output.open('x') as handle:
        json.dump(encode(report),handle,indent=2,sort_keys=True)
        handle.write('\n')
    print(report['status'])


if __name__ == '__main__':
    main()
