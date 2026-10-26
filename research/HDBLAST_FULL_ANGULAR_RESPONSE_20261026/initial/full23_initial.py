#!/usr/bin/env python3
"""Bounded full23 initial candidate on the retained r <= 7/4 footprint.

The continuous rational displacement is expression-defined.  P1 tables only
approximate displacement/radial derivative in feedback and define velocity.
Actual-minus-nominal errors use the retained three-channel paired algebra.
Every project import follows authentication of all seven complete roots.
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
LABELS = {
    'current': 'HDBLAST_COUPLED_MEAN_PILOT_20261025',
    'I': 'HDBLAST_INPUT_PRECISION_20261023',
    'P': 'HDBLAST_COMMON_GAUGE_FORCING_20261016',
    'N': 'HDBLAST_NUMERICAL_RESPONSE_INPUTS_20261016',
    'G': 'HDBLAST_GEOMETRIC_BULK_SOURCES_20261016',
    'C': 'HDBLAST_FORCED_RESPONSE_BOUNDS_20261023',
    'T': 'HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023',
}
PINS = {
    'current': 'b932157f308c9912fa83ee06990f1132a4e75a498aafd66eb99dabcc30b256de',
    'I': 'c2e7bd40e94422c5fd98820941231c0da63a8392b14a5a14c4340eabdde37f61',
    'P': '8bf173293b5cecf7f78e138ac5004abd9c7b1821ef29addb7ecfbc2a6409b4f6',
    'N': '928f3e4e58d222a3c7e1274babf9e94930ec501b74cf6ba04209f97c6bbded44',
    'G': '3ff1f95c151d2da6a1ab5e9b54718c2bb325c10d9c62a7ed9f97ef2bb36f0bb2',
    'C': '72924328a4e79cc77f0da90999b8dd56956bfdebd5e614e15c3c0498c3a8ae76',
    'T': 'cffb213a69f7f9a27932e990ffdecb7569d1007e5e5baeda708d8dc39d0cce72',
}
COEFFICIENT_PIN = '739240283c2bee53afa559c4ce71c1e3356cf2d0ad4e4b0f18b9621c077e20d3'
KINDS = ('u', 'du', 'v')
CHANNELS = ('inherited_shell_Taylor_coefficient_difference',
            'rational_packet_cutoff_and_E_approximation',
            'validated_fifth_ODE_derivative_Taylor_remainder_outward_arithmetic')


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    return value


def authenticate(roots):
    result = {}
    for key, root in roots.items():
        need(not root.is_symlink(), 'ROOT_SYMLINK:' + key)
        manifest = root/'MANIFEST.json'
        need(sha(manifest) == PINS[key], 'MANIFEST:' + key)
        seen = {'MANIFEST.json'}
        for entry in json.loads(manifest.read_text())['files']:
            rel = Path(entry['path'])
            need(not rel.is_absolute() and '..' not in rel.parts and str(rel) not in seen,
                 'MEMBER_PATH:' + key)
            seen.add(str(rel))
            path = root/rel
            need(path.is_file() and not path.is_symlink(), 'MEMBER_TYPE:' + str(path))
            need(path.stat().st_size == entry['bytes'] and sha(path) == entry['sha256'],
                 'MEMBER_BYTES:' + str(path))
        for path in root.rglob('*'):
            need(not path.is_symlink(), 'TREE_SYMLINK:' + str(path))
            if path.is_file():
                need(str(path.relative_to(root)) in seen, 'UNLISTED_FILE:' + str(path))
        result[key] = {'manifest_sha256': PINS[key], 'files_including_manifest': len(seen)}
    return result


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sqrt_upper(value):
    need(value >= 0, 'NEGATIVE_SQUARE')
    if not value:
        return F(0)
    grid = 1 << 112
    return F(math.isqrt((value*grid*grid).__ceil__())+1, grid)


class Context:
    def __init__(self, repo_root):
        self.repo = Path(repo_root).resolve()
        self.roots = {key: self.repo/'research'/name for key, name in LABELS.items()}
        self.before = authenticate(self.roots)
        self.init = load_module('angular23_interval', self.roots['N']/'numeric/evaluate_initial_profiles.py')
        self.op = load_module('angular23_operator', self.roots['G']/'projection/inputs/linear_operator.py')
        self.adapter = load_module('angular23_adapter', self.roots['G']/'projection/inputs/source_adapters.py')
        self.I, self.J, self.ell = self.init.I, self.init.J, self.init.ELL
        need(self.ell == F(1, 10**6) and self.init.GRID == 1 << 112, 'ARITHMETIC_CONVENTION')
        self.frozen = json.loads((self.roots['I']/'initial/RECENTERED_CANDIDATE_NORMAL.json').read_text())
        self.frozen_profiles = json.loads((self.roots['current']/'results/PROFILES.json').read_text())
        cs = self.frozen['representation_background_coefficients_ascending']
        digest = hashlib.sha256(json.dumps(cs, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
        need(digest == COEFFICIENT_PIN == self.frozen['representation_background_coefficients_sha256'],
             'FIXED_BACKGROUND_COEFFICIENTS')
        need(self.frozen['full_footprint_covered'] and self.frozen['passes_1e28'], 'FROZEN_INITIAL_GATE')
        self.cs = {key: list(map(F, value)) for key, value in cs.items()}
        self.data, self.input_pins = self.init.load()
        self.weights = self.op.gram()
        self.mean_ids = [self.op.IDX[name+':0'] for name in ('K', 'L', 'N', 'J', 'w')]
        need(self.mean_ids == [0, 3, 6, 9, 12], 'MEAN_ORDER')
        self.power = 1 << 40
        self.P = self.make_pair_type()

    def reauthenticate(self):
        need(authenticate(self.roots) == self.before, 'UPSTREAM_INPUT_CHANGED')
        return self.before

    def constant(self, value, order=6):
        return self.J([value]+[0]*order)

    def polynomial(self, values, z):
        result = self.constant(0, z.n)
        for value in reversed(values):
            result = result*z+value
        return result

    def square40(self, value):
        for _ in range(40):
            value = value*value
        return value

    def nominal_background(self, r, order=6):
        rr = self.J([r, 1/self.ell]+[0]*(order-1))
        z = self.ell*rr
        bg = {key: self.polynomial(values, z) for key, values in self.cs.items()}
        return bg, rr, z

    def _registered(self, bg, E, packet):
        rho, eta, p, a, alpha = [bg[key] for key in ('rho', 'eta', 'p', 'a', 'alpha')]
        h = self.init.poly(self.init.UP[1], eta)/p
        k = rho*(F(5, 2)*a-h)
        ir, ip = 1/rho, 1/p
        dy = lambda value: -ir*value.dz()
        u = alpha/rho**2*packet
        v = 3*alpha/(rho**3*p)*(packet.dz()+k*packet)
        vy, ut = dy(v), E*u
        d = 3*ir*ip
        vt = E*(v-d*u)
        U2 = self.init.poly(self.init.UP[2], eta)
        aa = (-2*vy**2-8*p*u*vy+8*p*(h+2*a)*u*v
              +(16*p*p-192*ir**2)*u*u-(2*U2+F(8, 3)*p*p)*v*v
              -ir**2*(12*ut*ut+2*vt*vt))
        dd = ir**2*(36*u*u-2*v*v)
        sectors = [aa/6+F(4, 3)*dd, aa/2+2*dd, aa-4*dd]
        fs = [-rho**2*sectors[j]/den for j, den in enumerate((-6, 10, 42))]
        # Retain the frozen mean expression evaluation order as well as its formula.
        fs[0] = rho**2*(aa/6+F(4, 3)*dd)/6
        fy = [dy(value) for value in fs]
        K = E*ir/3*(18*u*u-v*v+d*u*v)
        Pi = E*ir*ip*(d*u*vy+(d*(h-ir)+10*p)*u*v
                       +(30*a-6*ir)*u*u+(2*ir-a)*v*v)
        zero = packet*0
        # source_adapters.initial_frame_data has exactly these formulas.
        # Fmom=3K and (Gmom-Fmom_y)/p=Pi, as in the retained explicit jet.
        Y2 = self.adapter.Y2
        values = {
            'yy': (zero,)*3, 'y0': (zero,)*3, 'ya': (zero,)*2,
            '00': tuple(-2*u*u*c for c in Y2), '0a': (zero,)*2,
            'ab': tuple(fs)+(zero,)*4,
            'scalar': tuple(-F(3, 2)*ip*value for value in fy),
        }
        y0 = tuple(self.ell*(rho*(-26*u*dy(u)-2*a*u*u)*c
                             -F(3, 2)*rho*fy[j]-3*rho*a*fs[j])
                   for j, c in enumerate(Y2))
        angular = self.adapter.vector_project(22*u*u)
        velocities = {
            'yy': tuple(self.ell*(2*rho*K-8*u*ut)*c for c in Y2),
            'y0': y0, 'ya': (zero,)*2,
            '00': tuple(self.ell*(-8*rho*K+44*u*ut)*c for c in Y2),
            '0a': (self.ell*(angular[0]-fs[1]/2), self.ell*(angular[1]-fs[2]/2)),
            'ab': tuple(self.ell*(2*rho*K+4*u*ut)*c for c in Y2)+(zero,)*4,
            'scalar': tuple(self.ell*(rho*Pi+2*u*vt)*c for c in Y2),
        }
        displacement = self.op.pack_fields(values)
        velocity = self.op.pack_fields(velocities)
        # Identical frozen five-mean formulas; these slots define the old means.
        mean_u = [zero, zero, -u*u/3, fs[0], -F(3, 2)*ip*fy[0]]
        mean_v = [self.ell*(2*rho*K-8*u*ut)/6,
                  self.ell*rho*((-26*u*dy(u)-2*a*u*u)/6-F(3, 2)*fy[0]-3*a*fs[0]),
                  self.ell*(-8*rho*K+44*u*ut)/6,
                  self.ell*(2*rho*K+4*u*ut)/6,
                  self.ell*(rho*Pi+2*u*vt)/6]
        for j, i in enumerate(self.mean_ids):
            displacement[i], velocity[i] = mean_u[j], mean_v[j]
        return displacement, velocity

    def nominal_profiles(self, r, order=6):
        r = self.I.cv(r)
        need(0 <= r.lo <= r.hi <= F(7, 4), 'NOMINAL_OUTSIDE_RETAINED_FOOTPRINT')
        if r.hi <= F(65, 64):
            return ([self.constant(0, order) for _ in range(23)],
                    [self.constant(0, order) for _ in range(23)])
        need(r.lo >= F(65, 64) and not (r.lo < F(66, 64) < r.hi), 'CUTOFF_BRANCH_CROSSING')
        bg, rr, z = self.nominal_background(r, order)
        dr = rr-F(3, 2)
        denominator = F(1, 4)-dr*dr
        denominator.c[0] = F(1, 4)-(r-F(3, 2))**2
        t = 1/denominator-4
        packet = self.square40(1/(1+t/self.power))
        if r.lo < F(66, 64):
            tt = 64*rr-65
            packet = packet*(126*tt**5-420*tt**6+540*tt**7-315*tt**8+70*tt**9)
        E = self.square40(1/(1+z/self.power))
        return self._registered(bg, E, packet)

    def profiles(self, r):
        u, v = self.nominal_profiles(r)
        values = [[item.val() for item in u],
                  [(self.ell*item.dz()).val() for item in u],
                  [item.val() for item in v]]
        seconds = [[2*self.ell**2*item.c[2] for item in u],
                   [6*self.ell**3*item.c[3] for item in u],
                   [2*self.ell**2*item.c[2] for item in v]]
        return values, seconds

    def actual_background(self, r):
        z = self.ell*r
        values = {}
        for key in ('rho', 'eta', 'p', 'a', 'alpha'):
            value = self.I(0)
            for coefficient in reversed(self.frozen['shell_Taylor_coefficients'][key]):
                value = value*z+self.I(*coefficient)
            values[key] = value+self.I(*self.frozen['uniform_fifth_normalized_ODE_derivative'][key])*z**5
        return values

    def make_pair_type(self):
        context = self
        class P:
            def __init__(self, a, b, d=None):
                self.a, self.b = a, b
                self.d = d if d is not None else [context.constant(0, a.n) for _ in range(3)]
            @property
            def n(self): return self.a.n
            def cv(self, value):
                return value if isinstance(value, P) else P(context.constant(value, self.n), context.constant(value, self.n))
            def __add__(self, value):
                value = self.cv(value)
                return P(self.a+value.a, self.b+value.b, [x+y for x, y in zip(self.d, value.d)])
            __radd__ = __add__
            def __neg__(self): return P(-self.a, -self.b, [-x for x in self.d])
            def __sub__(self, value): return self+-self.cv(value)
            def __rsub__(self, value): return self.cv(value)+-self
            def __mul__(self, value):
                value = self.cv(value)
                return P(self.a*value.a, self.b*value.b,
                         [x*value.b+self.a*y for x, y in zip(self.d, value.d)])
            __rmul__ = __mul__
            def inv(self):
                ia, ib = self.a.inv(), self.b.inv()
                return P(ia, ib, [-ia*x*ib for x in self.d])
            def __truediv__(self, value): return self*self.cv(value).inv()
            def __rtruediv__(self, value): return self.cv(value)*self.inv()
            def __pow__(self, power):
                if power < 0: return self.inv()**(-power)
                result = self.cv(1)
                for _ in range(power): result = result*self
                return result
            def dz(self): return P(self.a.dz(), self.b.dz(), [x.dz() for x in self.d])
        return P

    def paired_background(self, r):
        z = self.ell*self.J([r, 1/self.ell, 0, 0, 0])
        bg = {}
        for key in ('rho', 'eta', 'p', 'a', 'alpha'):
            nominal = self.polynomial(self.cs[key], z)
            coefficients = [self.I(*entry) for entry in self.frozen['shell_Taylor_coefficients'][key]]
            fifth = self.I(*self.frozen['uniform_fifth_normalized_ODE_derivative'][key])
            branch, remainder = [], []
            for n in range(5):
                branch.append(sum((math.comb(j, n)*(coefficients[j]-self.cs[key][j])*z.c[0]**(j-n)
                                   for j in range(n, 5)), self.I(0)))
                remainder.append(math.comb(5, n)*fifth*z.c[0]**(5-n))
            differences = [self.J(branch), self.constant(0, 4), self.J(remainder)]
            bg[key] = self.P(nominal+differences[0]+differences[2], nominal, differences)
        return bg

    def paired_kernel(self, t):
        b = 1+t/self.power
        rational = self.square40(1/b)
        t0, b0 = t.val(), 1+t.val()/self.power
        R = rational.val()
        diffs, rising = [], F(1)
        for k in range(5):
            if k: rising *= 1+F(k-1, self.power)
            coefficients = [F(math.comb(k, j), self.power**j) for j in range(k+1)]
            coefficients[0] -= rising
            pdiff = self.init.poly(coefficients, t0)
            D0 = self.I(-(R*t0**2/(2*self.power)).hi, 0)
            diffs.append((-1)**k*(D0+R*pdiff/(b0**k)))
        dt = self.J([self.I(0)]+t.c[1:])
        difference, actual, power_dt = self.constant(0, 4), self.constant(0, 4), self.constant(1, 4)
        E = self.init.expbox(t0)
        for k in range(5):
            difference = difference+power_dt*diffs[k]/math.factorial(k)
            actual = actual+power_dt*((-1)**k*E)/math.factorial(k)
            power_dt = power_dt*dt
        return self.P(actual, rational, [self.constant(0, 4), difference, self.constant(0, 4)])

    def paired_packet(self, r):
        rr = self.J([r, 1/self.ell, 0, 0, 0])
        if r.hi <= F(65, 64):
            actual = self.init.bump(r, self.data, True)
            return self.P(actual, self.constant(0, 4), [self.constant(0, 4), actual, self.constant(0, 4)])
        need(r.lo >= F(65, 64) and not (r.lo < F(66, 64) < r.hi), 'PAIRED_CUTOFF_BRANCH_CROSSING')
        dr = rr-F(3, 2)
        denominator = F(1, 4)-dr*dr
        denominator.c[0] = F(1, 4)-(r-F(3, 2))**2
        q = self.paired_kernel(1/denominator-4)
        if r.lo < F(66, 64):
            tt = 64*rr-65
            S = 126*tt**5-420*tt**6+540*tt**7-315*tt**8+70*tt**9
            q = q*self.P(self.constant(1, 4), S,
                         [self.constant(0, 4), self.constant(1, 4)-S, self.constant(0, 4)])
        oracle = self.init.bump(r, self.data, False)
        q.a = self.J([self.I(max(x.lo, y.lo), min(x.hi, y.hi)) for x, y in zip(q.a.c, oracle.c)])
        return q

    def paired_profiles(self, r):
        r = self.I.cv(r)
        need(1 <= r.lo <= r.hi <= F(7, 4), 'PAIRED_OUTSIDE_RETAINED_FOOTPRINT')
        bg = self.paired_background(r)
        E = self.paired_kernel(self.ell*self.J([r, 1/self.ell, 0, 0, 0]))
        u, v = self._registered(bg, E, self.paired_packet(r))
        return [u, [self.ell*item.dz() for item in u], v]

    def norm(self, caps):
        return sqrt_upper(sum(w*c*c for w, c in zip(self.weights, caps)))


def produce(context, output_dir, cells=256):
    need(cells == 256, 'REGISTERED_FIXED_MESH_CAP_256')
    output_dir = Path(output_dir).resolve()
    need(not output_dir.exists(), 'OUTPUT_DIRECTORY_MUST_BE_NEW')
    need(all(root not in output_dir.parents for root in context.roots.values()), 'OUTPUT_INSIDE_FROZEN_ROOT')
    I = context.I
    nodes = [F(0)]+[F(k, cells) for k in range(65*cells//64, 7*cells//4+1)]
    need(nodes == list(map(F, context.frozen_profiles['nodes'])), 'FROZEN_SELECTED_MESH')
    node_cache = {}
    frozen_node_map = {F(r): j for j, r in enumerate(context.frozen_profiles['nodes'])}
    mean_evaluation_overlap_checks = 0
    for r in nodes:
        values, _ = context.profiles(I(r))
        centers = [[(box.lo+box.hi)/2 for box in row] for row in values]
        for kind, name in enumerate(KINDS):
            for j, i in enumerate(context.mean_ids):
                retained = F(context.frozen_profiles['profiles'][name][frozen_node_map[r]][j])
                # The rational retained center may be a last-bit rounding distance
                # outside a differently ordered expression enclosure.  Reuse it
                # and charge the exact maximum endpoint distance below.
                centers[kind][i] = retained
                mean_evaluation_overlap_checks += 1
        radii = [[max(abs(box.lo-center), abs(box.hi-center)) for box, center in zip(row, c)]
                 for row, c in zip(values, centers)]
        node_cache[r] = centers, radii
    interpolation_caps = [[F(0)]*23 for _ in range(3)]
    nominal_caps = [[F(0)]*23 for _ in range(3)]
    derivative_caps = [[F(0)]*23 for _ in range(3)]
    interpolation_rows = []
    for left, right in zip(nodes, nodes[1:]):
        if right <= F(65, 64):
            caps = [[F(0)]*23 for _ in range(3)]
            seconds = [[I(0)]*23 for _ in range(3)]
            values = [[I(0)]*23 for _ in range(3)]
        else:
            values, seconds = context.profiles(I(left, right))
            caps = [[max(node_cache[left][1][kind][i], node_cache[right][1][kind][i])
                     +(right-left)**2*seconds[kind][i].mag()/8 for i in range(23)]
                    for kind in range(3)]
        for kind in range(3):
            for i in range(23):
                interpolation_caps[kind][i] = max(interpolation_caps[kind][i], caps[kind][i])
                nominal_caps[kind][i] = max(nominal_caps[kind][i], values[kind][i].mag())
                derivative_caps[kind][i] = max(derivative_caps[kind][i], seconds[kind][i].mag())
        interpolation_rows.append({'r': [left, right], 'component_error_caps': caps,
                                   'H_error_caps': [context.norm(row) for row in caps],
                                   'second_derivative_component_caps': [[box.mag() for box in row] for row in seconds]})
    channels = [[[F(0)]*23 for _ in range(3)] for _ in range(3)]
    error_rows, mean_dominated_checks = [], 0
    old_rows = context.frozen['rows']
    for k in range(cells, 7*cells//4):
        left, right = F(k, cells), F(k+1, cells)
        paired = context.paired_profiles(I(left, right))
        errors = [[[item.d[c].val().mag() for item in row] for row in paired] for c in range(3)]
        for c in range(3):
            for kind in range(3):
                channels[c][kind] = [max(x, y) for x, y in zip(channels[c][kind], errors[c][kind])]
        old = next(row for row in old_rows if F(row['r'][0]) <= left and right <= F(row['r'][1]))
        for c in range(3):
            for kind in range(3):
                for j, i in enumerate(context.mean_ids):
                    # Existing wider-cell certificates independently cap the same
                    # mathematical five mean differences, even if expression
                    # ordering changes the final outward arithmetic bit.
                    retained = F(old['channel_component_error_caps'][c][kind][j])
                    if errors[c][kind][i] <= retained:
                        mean_dominated_checks += 1
        error_rows.append({'r': [left, right], 'analytic_flat_endpoint_tail': right <= F(65, 64),
                           'channel_component_error_caps': errors})
        if (k-cells) % 16 == 0:
            print(json.dumps({'full23_initial_cell': k, 'fixed_mesh': cells}), flush=True)
    channel_norms = [[context.norm(row) for row in channel] for channel in channels]
    total_errors = [sum((channel_norms[c][kind] for c in range(3)), F(0)) for kind in range(3)]
    interpolation_norms = [context.norm(row) for row in interpolation_caps]
    parity = [(-1 if name.startswith(('L:', 'v:')) else 1) for name in context.op.LABELS]
    report = {
        'status': 'PASS_CONTINUOUS_FULL23_INITIAL_CANDIDATE_AND_ERROR_LEDGER',
        'scope': 'Retained exact branch/operator premises; full23 initial profiles only, on r in [0,7/4].',
        'labels': context.op.LABELS, 'gram': context.weights, 'ell': context.ell,
        'fixed_mesh_denominator': cells, 'candidate_nodes': len(nodes),
        'continuous_error_cells': len(error_rows),
        'footprint': {'r': ['0', '7/4'], 'time': '0', 'whole_sphere': True,
                      'actual_zero_on': ['0', '1'], 'nominal_zero_on': ['0', '65/64'],
                      'cutoff_transition': ['65/64', '66/64'],
                      'no_claim_beyond_retained_footprint': True},
        'candidate_definition': 'Exact expression-defined rational u; P1 u/du for feedback only; P1 v defines the free-wave velocity.',
        'candidate_background_coefficients_sha256': COEFFICIENT_PIN,
        'channels': CHANNELS,
        'channel_component_error_caps': channels,
        'channel_H_caps': channel_norms,
        'total_H_caps': total_errors,
        'total_H_caps_display': list(map(float, total_errors)),
        'P1_component_error_caps': interpolation_caps,
        'P1_H_error_caps': interpolation_norms,
        'P1_H_error_caps_display': list(map(float, interpolation_norms)),
        'nominal_component_caps': nominal_caps,
        'nominal_H_caps': [context.norm(row) for row in nominal_caps],
        'nominal_H_caps_display': [float(context.norm(row)) for row in nominal_caps],
        'nominal_second_derivative_component_caps': derivative_caps,
        'nominal_second_derivative_H_caps': [context.norm(row) for row in derivative_caps],
        'derivative_convention': {'u': 'U', 'du': 'partial_r U', 'v': 'partial_x U; x=t/ell',
                                  'second_derivatives': ['partial_r^2 U', 'partial_r^3 U', 'partial_r^2 V']},
        'mean_consistency': {'retained_mean_node_centers_reused': mean_evaluation_overlap_checks,
                             'mean_rows_and_columns': context.mean_ids,
                             'new_error_entries_dominated_by_frozen_wider_cells': mean_dominated_checks,
                             'entries_compared': len(error_rows)*3*3*5,
                             'changed_expression_rounding_charged_in_endpoint_radii': True},
        'regularity': 'Actual data smoothly vanish at r=1. Rational cutoff polynomial matches derivatives through4 at65/64,66/64; U and V are C1 on the retained footprint. Their higher one-sided derivative bounds are used only on cells inside one analytic branch.',
        'arithmetic': 'All decisions exact rational; 112-bit outward intervals; expression R=(1+t/2^40)^(-2^40), evaluated by40 squarings. Three telescoping difference channels avoid subtracting whole profile ranges.',
        'authentication_before': context.before,
        'complete_closures_authenticated_before_and_after': True,
        'interpolation_cells': interpolation_rows, 'rows': error_rows,
        'limits': ['No full spacetime response has been computed in this producer.',
                   'No physical origin, nonlinear family, or finite amplitude conclusion.',
                   'Retained shell uncertainty is enclosed, not reduced.'],
    }
    profiles = {
        'schema': 'HDBLAST_FULL23_INITIAL_P1_PROFILES_V1', 'labels': context.op.LABELS,
        'gram': context.weights, 'ell': context.ell,
        'candidate_coefficients_sha256': COEFFICIENT_PIN,
        'nodes': nodes,
        'profiles': {name: [node_cache[r][0][kind] for r in nodes] for kind, name in enumerate(KINDS)},
        'parity': {'u': parity, 'du': [-p for p in parity], 'v': parity},
        'initial_displacement': 'exact registered full23 rational displacement; P1 u/du only inside feedback',
    }
    context.reauthenticate()
    output_dir.mkdir(parents=True)
    for filename, value in (('INITIAL_REPORT.json', report), ('INITIAL_PROFILES.json', profiles)):
        with (output_dir/filename).open('x') as handle:
            json.dump(encode(value), handle, indent=2, sort_keys=True)
            handle.write('\n')
    print(json.dumps({key: report[key] for key in ('status', 'total_H_caps_display', 'P1_H_error_caps_display', 'nominal_H_caps_display')}), flush=True)
    return report, profiles


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    produce(Context(args.repo_root), args.output_dir)


if __name__ == '__main__':
    main()
