#!/usr/bin/env python3
"""Authenticated interval ranges of the fixed W0 residual on alpha<=9/8.

All intervals are pairs of Fraction endpoints.  Every quadrature panel encloses
the whole panel, rather than sampled values.  Actual affine F and G are retained
as fixed, integral-defined functions enclosed by their authenticated uniform
caps; neither source is evaluated at a point here.  The finite Robin series is
always supplemented by its analytic tail.  No Einstein residual is claimed.
"""
from bisect import bisect_right
from fractions import Fraction as F
from pathlib import Path
import hashlib
import importlib.util
import json
import math
import sys

sys.dont_write_bytecode = True

ZERO = (F(0), F(0))
GRID = 1 << 120
LIMIT = F(9, 8)


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def pin(path, digest):
    raw = Path(path).read_bytes()
    need(hashlib.sha256(raw).hexdigest() == digest, 'INPUT_HASH:' + str(path))
    return raw


def interval(z):
    if hasattr(z, 'lo'):
        return F(z.lo), F(z.hi)
    if isinstance(z, (tuple, list)):
        need(len(z) == 2, 'interval shape')
        lo, hi = map(F, z)
    else:
        lo = hi = F(z)
    need(lo <= hi, 'ordered interval')
    return lo, hi


def rounded(z):
    """Outward dyadic rounding controls denominator growth, with exact zeros."""
    lo, hi = z
    return (F((lo.numerator * GRID) // lo.denominator, GRID),
            F(-((-hi.numerator * GRID) // hi.denominator), GRID))


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def neg(a):
    return -a[1], -a[0]


def sub(a, b):
    return add(a, neg(b))


def mul(a, b):
    if a == ZERO or b == ZERO:
        return ZERO
    if a[0] == a[1]:
        return scale(b, a[0])
    if b[0] == b[1]:
        return scale(a, b[0])
    z = (a[0]*b[0], a[0]*b[1], a[1]*b[0], a[1]*b[1])
    return min(z), max(z)


def scale(a, c):
    if not c or a == ZERO:
        return ZERO
    return (a[0]*c, a[1]*c) if c > 0 else (a[1]*c, a[0]*c)


def widen(a, e):
    return a[0]-e, a[1]+e


def hull(*rows):
    return min(z[0] for z in rows), max(z[1] for z in rows)


def vadd(a, b):
    return [rounded(add(x, y)) for x, y in zip(a, b)]


def vsub(a, b):
    return [rounded(sub(x, y)) for x, y in zip(a, b)]


def vscale(a, c):
    return [rounded(scale(x, c)) for x in a]


def vmul(a, c):
    return [rounded(mul(x, c)) for x in a]


def vtotal(rows, dim=23):
    out = [ZERO]*dim
    for row in rows:
        out = vadd(out, row)
    return out


def mv(matrix, vector):
    out = []
    for row in matrix:
        z = ZERO
        for c, v in zip(row, vector):
            if c != ZERO and v != ZERO:
                z = add(z, mul(c, v))
        out.append(rounded(z))
    return out


def frozen_mv(matrix, vector):
    out = []
    for row in matrix:
        z = ZERO
        for c, v in zip(row, vector):
            if c and v != ZERO:
                z = add(z, scale(v, c))
        out.append(rounded(z))
    return out


def sqrt_upper(q):
    need(q >= 0, 'nonnegative norm square')
    grid = 1 << 64
    k = math.isqrt(q.numerator*grid*grid//q.denominator)
    if k*k*q.denominator < q.numerator*grid*grid:
        k += 1
    return F(k, grid)


PINS = {
    'driver': 'efa6a6c2a9c7db596fe024246684f760b807e727fd2ddac8e8126182ba941a73',
    'initial': 'bbde14bab0f52dc4371ee016aba48dd6e91874f610d09c197a942f8c8f2e1fad',
    'profiles': '00165ede35939b768150a9cb9b04cd0b01f0ec6be6d54c8dbcf1b37b90959fdc',
    'response': 'ce1a07ea7409253e679e84090ae9ffbcfa9a53dbfe8c620600f2b022ae3431fd',
    'F': '45a5550d2e4647bdfb84f9c4a36bf9861ef500a1de89a76f978c5afe9c51a60f',
    'Fmap': '34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615',
    'G': 'fec70c6c305ceb55b3c1b9917b4cbf5abb25c403e0399770877d2f3dc5b2bf2e',
    'local': '1832fba309ceaf0a138bbddc0aaf14602f81796bc1a7aa8cd938a74470e1b899',
}


class InteriorEvaluator:
    def __init__(self, root, primitive_panels=128, contact_capsule=None,
                 bulk_capsule=None):
        self.root = Path(root).resolve()
        need(1 <= primitive_panels <= 4096, 'primitive panel count')
        self.primitive_panels = primitive_panels
        base = self.root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026'
        contact = Path(contact_capsule) if contact_capsule else self.root/'research/HDBLAST_SMOOTH_SHELL_CONTACT_20261030'
        bulk = Path(bulk_capsule) if bulk_capsule else self.root/'research/HDBLAST_ACTUAL_BULK_UPDATE_20261030'
        self.paths = {
            'driver': base/'compute/full23_response.py',
            'initial': base/'initial/results/INITIAL_REPORT.json',
            'profiles': base/'compute/response-normal/RESPONSE_PROFILES.json',
            'response': base/'compute/response-normal/RESPONSE_REPORT.json',
            'F': self.root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/ACTUAL_SOURCES_NORMAL.json',
            'Fmap': self.root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/GEOMETRIC_COEFFICIENT_MAPS.json',
            'G': contact/'reference/AFFINE_SOURCE.json',
            'local': bulk/'reference/LOCAL.json',
        }
        data = {k: pin(p, PINS[k]) for k, p in self.paths.items()}
        initial, profile, report, source, Gsource, local = [json.loads(data[k]) for k in
            ('initial', 'profiles', 'response', 'F', 'G', 'local')]
        spec = importlib.util.spec_from_file_location('interior_authenticated_response', self.paths['driver'])
        self.module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.module)
        self.response = self.module.Response(self.root)
        self.ctx = self.response.context
        self.I = self.response.I
        self.labels = list(self.response.labels)
        self.weights = list(self.response.weights)
        self.parity = list(map(int, profile['parity']['u']))
        self.P = self.parity
        need(len(self.labels) == 23 and self.labels == profile['labels'] == initial['labels'] == report['labels'], 'registered labels')
        need(self.weights == list(map(F, profile['gram'])) == list(map(F, initial['gram'])) and min(self.weights) > 0, 'registered Gram')
        need(self.parity == [-1 if z.startswith(('L:', 'v:')) else 1 for z in self.labels], 'mixed reflection parity')
        need(initial['status'] == 'PASS_CONTINUOUS_FULL23_INITIAL_CANDIDATE_AND_ERROR_LEDGER', 'initial scope')
        need(report['status'] == 'PASS_BOUNDED_FULL23_RESPONSE_X_ERROR', 'response scope')
        self.frozen = {key: [list(map(F, row)) for row in profile['matrices'][key]] for key in ('A', 'B', 'D', 'M')}
        need(all(profile['matrices'][key] == report['constant_matrices'][key] for key in self.frozen), 'frozen matrices agreement')
        self.plus = [[a+b for a, b in zip(ar, br)] for ar, br in zip(self.frozen['A'], self.frozen['B'])]
        self.minus = [[b-a for a, b in zip(ar, br)] for ar, br in zip(self.frozen['A'], self.frozen['B'])]
        self.M = [[interval(z) for z in row] for row in report['actual_coefficient_boxes']['M']]
        shellM = self.response.assemble(self.I(0), self.I(0))['M']
        need(self.M == [[interval(z) for z in row] for row in shellM], 'actual constant shell M')
        self.m = F(report['constant_matrix_H_operator_upper']['M']) + F(report['coefficient_error_H_operator_upper']['M'])
        need(0 < self.m*(LIMIT-1) < 1, 'local Volterra ratio')
        for i, row in enumerate(self.M):
            for j, z in enumerate(row):
                need(self.parity[i] > 0 and self.parity[j] > 0 or z == ZERO, 'actual Dirichlet Robin zero')
        self.factors = [sqrt_upper(1/w) for w in self.weights]
        self.Hcap = F(local['auxiliary_trace_h_local_H_cap'])
        need(local['status'] == 'PASS' and local['point'] == ['9/8', '0'], 'local uniform trace cap')
        self.G = [interval(row['interval']) for row in Gsource['rows']]
        self.Gcap = [max(map(abs, z)) for z in self.G]
        self.GH = F(Gsource['source_H_norm_strict_integer_upper'])
        need(Gsource['uniform_domain'] == '0<=x<=3,r=0' and [r['label'] for r in Gsource['rows']] == self.labels, 'actual shell source domain/labels')
        need(all(z == ZERO for z, p in zip(self.G, self.parity) if p < 0), 'actual zero Dirichlet source')
        self._source_caps(source, data['Fmap'])
        self.du, self.dv = F(1, 2**30), F(1, 2048)
        self.u_collar = F(65, 64)-self.du
        self.v_collar = F(65, 64)-self.dv
        self.collar = min(self.u_collar, self.v_collar)
        M2 = list(map(F, initial['nominal_second_derivative_component_caps'][0]))
        need(len(M2) == 23 and min(M2) >= 0, 'second derivative caps')
        self.alpha = [z*self.du**2/18 for z in M2]
        self.beta = [35*z*self.du/128 for z in M2]
        nodes = list(map(F, profile['nodes']))
        vv = [list(map(F, row)) for row in profile['profiles']['v']]
        need(all(z == 0 for row in vv[:2] for z in row) and nodes[1]-max(self.du, self.dv) > 1, 'flat smoothed initial collar')
        slopes = [[(b-a)/(hi-lo) for a, b in zip(va, vb)] for lo, hi, va, vb in zip(nodes, nodes[1:], vv, vv[1:])]
        need(not any(slopes[0]), 'zero first velocity slope')
        self.jumps = [(nodes[k], [b-a for a, b in zip(slopes[k-1], slopes[k])]) for k in range(1, len(slopes)) if nodes[k]-self.dv < LIMIT]
        self._profile_cache = {}
        self._initial_cache = {}
        self._ridge_cache = {}
        self._lift_cache = {}
        self._coefficient_cache = {}
        self._build_primitives()
        self.calls = 0

    def _source_caps(self, source, map_bytes):
        need(source['source_map'] == {'bytes': len(map_bytes), 'sha256': PINS['Fmap'], 'status': 'PASS_FORMAL_ALL_POINT_GEOMETRIC_SOURCE_PROJECTION'}, 'complete source map')
        need(source['status'] == 'PASS_CONTINUOUS_ACTUAL_EVOLVED_GEOMETRIC_SOURCE_CELL_ENCLOSURES', 'actual affine F scope')
        need(source['scope']['actual_evolved_common_gauge_bulk_sources_evaluated'] is True, 'common gauge source')
        need(source['domain']['causal_triangle'] == '0<=x<=3,r>=0,x+r<=8', 'uniform F domain')
        need(source['normalization'] == 'Second parameter derivatives along the affine actual-first-pair path; metric_forcing=-2 ell^2 rho^2 QEred; scalar_forcing=ell^2 rho^2 Qphi. QC is the actual-metric-lowered quadratic gauge Hessian.', 'forcing normalization')
        summary = source['global_cover']['summary']
        caps = summary['absolute_coefficient_caps']
        hb = summary['positive_bundle_bulk_forcing_H_bound']
        mapping = [('yy', 'metric_forcing/yy', ('K:0', 'K:2', 'K:4'), 1), ('y0', 'metric_forcing/y0', ('L:0', 'L:2', 'L:4'), 2), ('00', 'metric_forcing/00', ('N:0', 'N:2', 'N:4'), 1), ('ya', 'metric_forcing/ya', ('v:2', 'v:4'), 2), ('0a', 'metric_forcing/0a', ('u:2', 'u:4'), 2), ('ab', 'metric_forcing/ab', ('J:0', 'J:2', 'J:4', 'H:2', 'H:4', 'T:+', 'T:-'), 1), ('scalar', 'scalar_forcing', ('w:0', 'w:2', 'w:4'), 1)]
        self.Fcap = [F(0)]*23
        square = F(0)
        count = 0
        fields = {}
        for field, prefix, labels, multiplicity in mapping:
            block = hb['block_Gram_certificates'][field]
            gg = [list(map(F, row)) for row in block['gram']]
            cc = [F(caps[prefix+'/'+str(j)]) for j in range(len(labels))]
            need(block['multiplicity'] == multiplicity and len(gg) == len(labels) and all(len(row) == len(labels) for row in gg) and min(cc) >= 0, 'source Gram shape')
            fields[field] = cc
            for i, label in enumerate(labels):
                self.Fcap[self.labels.index(label)] = cc[i]
                for j in range(len(labels)):
                    need(multiplicity*gg[i][j] == (self.weights[self.labels.index(label)] if i == j else 0), 'source/current Gram identity')
                    count += 1
            bsquare = multiplicity*sum(abs(gg[i][j])*cc[i]*cc[j] for i in range(len(labels)) for j in range(len(labels)))
            need(bsquare == F(block['squared_norm_upper']), 'source block norm')
            square += bsquare
        need(self.ctx.op.pack_fields(fields) == self.Fcap and count == 93, 'complete source packing')
        need(square == F(hb['squared_upper']) == sum(w*c*c for w, c in zip(self.weights, self.Fcap)), 'complete H source norm')
        self.FH = F(hb['F_H_upper'])
        need(self.FH == 17387513156573 and (self.FH-1)**2 <= square < self.FH**2, 'genuine pointwise F norm')
        self.affine_F = [(-c, c) for c in self.Fcap]

    def _profiles(self, z):
        z = interval(z)
        if z in self._profile_cache:
            return self._profile_cache[z]
        need(0 <= z[0] <= z[1] <= LIMIT, 'positive profile domain')
        if z[1] <= F(65, 64):
            out = ([ZERO]*23, [ZERO]*23)
        elif getattr(self, '_profile_grid_ready', False) and z[0] < z[1]:
            # A broad interval fed directly into the rational cutoff expression
            # can have enormous dependency inflation.  Hull the authenticated
            # fine-grid whole-panel ranges instead, including every intersected
            # panel; this remains an enclosure of the complete query interval.
            k0 = min(len(self.grid)-2, max(0, bisect_right(self.grid, z[0])-1))
            k1 = min(len(self.grid)-2, max(0, bisect_right(self.grid, z[1])-1))
            rows = self.panel_profiles[k0:k1+1]
            out = tuple([[hull(*(row[k][i] for row in rows)) for i in range(23)] for k in range(2)])
        else:
            cuts = sorted(set([z[0], z[1]] + [c for c in (F(65, 64), F(66, 64)) if z[0] < c < z[1]]))
            rows = []
            if len(cuts) == 1:
                cuts.append(cuts[0])
            for lo, hi in zip(cuts, cuts[1:]):
                values, _ = self.ctx.profiles(self.I(lo, hi))
                rows.append(([rounded(interval(v)) for v in values[0]], [rounded(interval(v)) for v in values[1]]))
            out = tuple([[hull(*(row[k][i] for row in rows)) for i in range(23)] for k in range(2)])
        self._profile_cache[z] = out
        return out

    def _build_primitives(self):
        n = self.primitive_panels
        first, middle = F(65, 64), F(66, 64)
        self.grid = sorted(set([F(0), first] + [lo+(hi-lo)*k/n for lo, hi in ((first, middle), (middle, LIMIT)) for k in range(n+1)]))
        self._profile_grid_ready = False
        self.panel_u = []
        self.panel_profiles = []
        self.prefix = [[[ZERO]*23 for _ in range(6)]]
        for left, right in zip(self.grid, self.grid[1:]):
            profiles = self._profiles((left, right))
            self.panel_profiles.append(profiles)
            urange = profiles[0]
            self.panel_u.append(urange)
            prev = self.prefix[-1]
            step = right-left
            nxt = [[ZERO]*23]
            for j in range(1, 6):
                rows = [vscale(prev[k], step**(j-k)/math.factorial(j-k)) for k in range(1, j+1)]
                rows.append(vscale(urange, step**j/math.factorial(j)))
                nxt.append(vtotal(rows))
            self.prefix.append(nxt)
        self._profile_grid_ready = True

    def _nominal_U_piece(self, lo, hi, k):
        left = self.grid[k]
        d = (lo-left, hi-left)
        need(0 <= d[0] <= d[1] <= self.grid[k+1]-left, 'local primitive piece')
        powers = [(F(1), F(1))] + [(d[0]**j/math.factorial(j), d[1]**j/math.factorial(j)) for j in range(1, 6)]
        prev = self.prefix[k]
        out = [[ZERO]*23]
        for j in range(1, 6):
            rows = [vmul(prev[h], powers[j-h]) for h in range(1, j+1)]
            rows.append(vmul(self.panel_u[k], powers[j]))
            out.append(vtotal(rows))
        return out

    def _positive_initial(self, z):
        lo, hi = z
        if hi <= self.collar:
            return {'U': [[ZERO]*23 for _ in range(6)], 'p': [ZERO]*23, 'V': [[ZERO]*23 for _ in range(7)]}
        unom, pnom = self._profiles(z)
        U = [[ZERO]*23 for _ in range(6)]
        U[0] = [rounded(widen(v, a)) for v, a in zip(unom, self.alpha)]
        p = [rounded(widen(v, b)) for v, b in zip(pnom, self.beta)]
        if hi <= self.u_collar:
            U[0] = p = [ZERO]*23
        k0 = min(len(self.grid)-2, max(0, bisect_right(self.grid, lo)-1))
        k1 = min(len(self.grid)-2, max(0, bisect_right(self.grid, hi)-1))
        pieces = [self._nominal_U_piece(max(lo, self.grid[k]), min(hi, self.grid[k+1]), k) for k in range(k0, k1+1)]
        for j in range(1, 6):
            rem_mass = max(F(0), hi-self.u_collar)**j/math.factorial(j)
            U[j] = [rounded(widen(hull(*(part[j][i] for part in pieces)), self.alpha[i]*rem_mass)) for i in range(23)]
        V = [self._velocity(z, j) for j in range(7)]
        return {'U': U, 'p': p, 'V': V}

    def initial(self, z):
        z = interval(z)
        need(-LIMIT <= z[0] <= z[1] <= LIMIT, 'initial range domain')
        if z in self._initial_cache:
            return self._initial_cache[z]
        if z[0] >= 0:
            out = self._positive_initial(z)
        elif z[1] <= 0:
            positive = self._positive_initial((-z[1], -z[0]))
            out = {'U': [[scale(v, (-1)**j*p) for v, p in zip(row, self.parity)] for j, row in enumerate(positive['U'])],
                   'V': [[scale(v, (-1)**j*p) for v, p in zip(row, self.parity)] for j, row in enumerate(positive['V'])],
                   'p': [scale(v, -p) for v, p in zip(positive['p'], self.parity)]}
        else:
            left, right = self.initial((z[0], 0)), self.initial((0, z[1]))
            out = {key: [[hull(a, b) for a, b in zip(ar, br)] for ar, br in zip(left[key], right[key])] for key in ('U', 'V')}
            out['p'] = [hull(a, b) for a, b in zip(left['p'], right['p'])]
        self._initial_cache[z] = out
        return out

    @staticmethod
    def _moment(k):
        return F(0) if k % 2 else F(35, 32)*sum(F(2*c, 2*j+k+1) for j, c in enumerate((1, -3, 3, -1)))

    def _hinge(self, z, n):
        if z <= -self.dv:
            return F(0)
        if z >= self.dv:
            return sum(F(math.comb(n, k), math.factorial(n))*z**(n-k)*(-self.dv)**k*self._moment(k) for k in range(n+1))
        s, out = z/self.dv, F(0)
        for j, c in enumerate((1, -3, 3, -1)):
            for k in range(n+1):
                power = 2*j+k
                out += c*math.comb(n, k)*s**(n-k)*(-1)**k*(s**(power+1)-(-1)**(power+1))/F(power+1)
        return F(35, 32)*self.dv**n*out/math.factorial(n)

    def _velocity(self, z, order):
        out = [ZERO]*23
        for knot, jump in self.jumps:
            if z[1] <= knot-self.dv:
                break
            # The scalar convolved normalized hinge is nondecreasing (its
            # derivative is a nonnegative convolved positive-part polynomial).
            h = (self._hinge(z[0]-knot, order+1), self._hinge(z[1]-knot, order+1))
            need(h[0] <= h[1], 'monotone exact velocity hinge')
            out = vadd(out, [scale(h, c) for c in jump])
        return out

    def project(self, vector, sign=1):
        return [v if p == sign else ZERO for v, p in zip(vector, self.parity)]

    def reflect(self, vector):
        return [scale(v, p) for v, p in zip(vector, self.parity)]

    def ridge(self, z, order, which):
        z = interval(z)
        key = z, order, which
        if key in self._ridge_cache:
            return self._ridge_cache[key]
        need(0 <= order <= 5 and which in ('F', 'G'), 'ridge order')
        data = self.initial(z)
        deriv = data['p'] if order == 0 else data['U'][order-1]
        if which == 'F':
            out = vscale(vadd(frozen_mv(self.plus, vadd(deriv, data['V'][order])), frozen_mv(self.frozen['D'], vadd(data['U'][order], data['V'][order+1]))), F(1, 2))
        else:
            out = vscale(vadd(frozen_mv(self.minus, vsub(deriv, data['V'][order])), frozen_mv(self.frozen['D'], vsub(data['U'][order], data['V'][order+1]))), F(1, 2))
        self._ridge_cache[key] = out
        return out

    def free_wave(self, a, q):
        ia, iq = self.initial(a), self.initial(q)
        return {'value': vscale(vadd(vadd(ia['U'][0], iq['U'][0]), vsub(ia['V'][1], iq['V'][1])), F(1, 2)),
                'dx': vscale(vadd(vsub(ia['p'], iq['p']), vadd(ia['V'][0], iq['V'][0])), F(1, 2)),
                'dr': vscale(vadd(vadd(ia['p'], iq['p']), vsub(ia['V'][0], iq['V'][0])), F(1, 2))}

    def frozen_kernel(self, x, r, a, q):
        f = lambda z, j=0: self.ridge(z, j, 'F')
        g = lambda z, j=0: self.ridge(z, j, 'G')
        fa1, fa2, ga1, ga2 = f(a, 1), f(a, 2), g(a, 1), g(a, 2)
        p_edge = vadd(vmul(f(a, 0), x), vscale(vsub(ga1, g(q, 1)), F(1, 2)))
        if -self.collar <= q[0] <= q[1] <= self.collar:
            # This branch also covers q straddling zero without a case split.
            value = vscale(vadd(vmul(fa1, x), vscale(vsub(ga2, fa2), F(1, 2))), F(1, 2))
            m_edge = vscale(fa1, F(1, 2))
        elif q[0] >= 0:
            value = vscale(vsub(vadd(vsub(vmul(fa1, x), vscale(vsub(fa2, f(q, 2)), F(1, 2))), vscale(vsub(ga2, g(q, 2)), F(1, 2))), vmul(g(q, 1), x)), F(1, 2))
            m_edge = vadd(vscale(vsub(fa1, f(q, 1)), F(1, 2)), vmul(g(q, 0), x))
        elif q[1] <= 0:
            t = neg(q)
            physical = vsub(vadd(vsub(vmul(fa1, x), vscale(vadd(fa2, f(t, 2)), F(1, 2))), vscale(vadd(ga2, g(q, 2)), F(1, 2))), vmul(g(q, 1), r))
            reflected = vadd(vsub(vmul(f(t, 1), t), f(t, 2)), vscale(vadd(g(t, 2), g(q, 2)), F(1, 2)))
            value = vscale(vadd(physical, self.reflect(reflected)), F(1, 2))
            reflected_m = vadd(vmul(f(t, 0), t), vscale(vsub(g(t, 1), g(q, 1)), F(1, 2)))
            m_edge = vadd(vadd(vscale(vsub(fa1, f(t, 1)), F(1, 2)), vmul(g(q, 0), r)), self.reflect(reflected_m))
        else:
            # Broad caller boxes may cross zero and an active collar.  Branch
            # restrictions are safe even while x,r retain their original boxes.
            left = self.frozen_kernel(x, r, a, (q[0], F(0)))
            right = self.frozen_kernel(x, r, a, (F(0), q[1]))
            return {k: [hull(l, h) for l, h in zip(left[k], right[k])] for k in left}
        return {'value': value, 'dx': vscale(vadd(p_edge, m_edge), F(1, 2)), 'dr': vscale(vsub(p_edge, m_edge), F(1, 2))}

    def h_lift(self, t):
        t = interval(t)
        if t in self._lift_cache:
            return self._lift_cache[t]
        if t[1] <= 1:
            out = {'value': [ZERO]*23, 'load': [ZERO]*23, 'trace': [ZERO]*23, 'homogeneous_tail_H': F(0), 'source_tail_H': F(0)}
            self._lift_cache[t] = out
            return out
        # The support premise gives zero through1.  Source caps alone do not.
        tl = (max(F(1), t[0]), t[1])
        need(tl[1] <= LIMIT, 'boundary history domain')
        initial = self.initial(tl)
        h = []
        for j in range(4):
            H = vadd(vsub(vmul(self.ridge(tl, j+1, 'F'), tl), vscale(self.ridge(tl, j+2, 'F'), j+1)), vscale(vadd(self.ridge(tl, j+2, 'G'), vscale(self.ridge(neg(tl), j+2, 'G'), (-1)**j)), F(1, 2)))
            h.append(self.project(vadd(vadd(initial['U'][j], initial['V'][j+1]), H)))
        terms = []
        for j in range(1, 4):
            term = h[j]
            for _ in range(j):
                term = mv(self.M, term)
            terms.append(vscale(term, (-1)**j))
        lift = vtotal(terms)
        d = tl[1]-1
        ratio = self.m*d
        htail = self.Hcap*ratio**4/math.factorial(4)/(1-ratio/5)
        gtail = self.m*d*d*self.GH/2/(1-ratio/3)
        for i in range(23):
            nonzero = any(z != ZERO for z in self.M[i])
            rem = (htail+gtail)*self.factors[i] if nonzero else F(0)
            lift[i] = rounded(widen(lift[i], d*self.Gcap[i]+rem))
        trace = vadd(h[0], lift)
        load = vadd(mv(self.M, trace), self.G)
        if t[0] <= 1:
            lift = [hull(z, ZERO) for z in lift]
            trace = [hull(z, ZERO) for z in trace]
            load = [hull(z, ZERO) for z in load]
        out = {'value': lift, 'load': load, 'trace': trace, 'homogeneous_tail_H': htail, 'source_tail_H': gtail}
        self._lift_cache[t] = out
        return out

    def _evaluate(self, x, r, a, q):
        self.calls += 1
        if a[1] <= 1:
            zjets = {key: [ZERO]*23 for key in ('value', 'dx', 'dr')}
            return {'r0_lower': [ZERO]*23, 'affine_F': [ZERO]*23, 'residual': [ZERO]*23, 'Bs': zjets, 'delta_W0': zjets}
        need(x[1] <= F(3, 2) and r[0] >= max(F(0), x[1]-F(5, 4))
             and x[1]+r[1] <= F(7, 4),
             'coefficient rectangle inside inherited retained footprint; split the caller box')
        Bs = self.free_wave(a, q)
        kernel = self.frozen_kernel(x, r, a, q)
        lift = self.h_lift(neg(q))
        delta = {'value': vadd(kernel['value'], lift['value']), 'dx': vsub(kernel['dx'], lift['load']), 'dr': vadd(kernel['dr'], lift['load'])}
        ckey = x, r
        if ckey not in self._coefficient_cache:
            assembled = self.response.assemble(self.I(*x), self.I(*r))
            self._coefficient_cache[ckey] = {k: [[rounded(interval(z)) for z in row] for row in assembled[k]] for k in ('A', 'B', 'D')}
        actual = self._coefficient_cache[ckey]
        rows = []
        for key, jet in (('A', 'dx'), ('B', 'dr'), ('D', 'value')):
            difference = [[sub(v, (nom, nom)) for v, nom in zip(ar, fr)] for ar, fr in zip(actual[key], self.frozen[key])]
            rows.extend((mv(difference, Bs[jet]), mv(actual[key], delta[jet])))
        lower = vtotal(rows)
        return {'r0_lower': lower, 'affine_F': self.affine_F, 'residual': vadd(lower, self.affine_F), 'Bs': Bs, 'delta_W0': delta}

    def evaluate(self, alpha, theta):
        """Whole α×θ box, x=α(1−θ)/2,r=α(1+θ)/2; no samples."""
        a, th = interval(alpha), interval(theta)
        need(0 <= a[0] <= a[1] <= LIMIT and -1 <= th[0] <= th[1] <= 1, 'physical folded domain')
        x = scale(mul(a, sub((F(1), F(1)), th)), F(1, 2))
        r = scale(mul(a, add((F(1), F(1)), th)), F(1, 2))
        # Correlation is ignored only outward; coefficient boxes remain physical.
        x = max(F(0), x[0]), min(LIMIT, x[1])
        r = max(F(0), r[0]), min(LIMIT, r[1])
        q = mul(a, th)
        return self._evaluate(x, r, a, q)

    def evaluate_xr(self, x, r):
        x, r = interval(x), interval(r)
        need(x[0] >= 0 and r[0] >= 0 and x[1]+r[1] <= LIMIT, 'physical rectangular domain')
        return self._evaluate(x, r, add(x, r), sub(r, x))

    def authenticate(self):
        self.ctx.reauthenticate()
        need(self.module.authenticate_local() == self.response.local_before, 'local dependency readback')
        for key, path in self.paths.items():
            pin(path, PINS[key])
        return True

    def metadata(self):
        return {'input_sha256': dict(PINS), 'primitive_panels_per_branch': self.primitive_panels,
                'primitive_method': 'outward exact-rational whole-panel nominal ranges and positive recursive convolution quadrature, plus uniform smoothing remainders',
                'velocity_method': 'exact compact-kernel smoothed P1 hinge primitives with monotone endpoint range enclosures',
                'rounding': 'every returned algebraic row rounded outward to 120 binary fractional bits',
                'source_method': 'same fixed actual functions; authenticated uniform F/G component caps, no source point evaluation',
                'robin_method': 'terms1..3 cancellation-free lift plus homogeneous factorial tail and full actual G convolution enclosure',
                'affine_F_H_uniform_cap': str(self.FH), 'actual_M_H_operator_upper': str(self.m),
                'support': 'inherited zero wedge alpha<=1, independently of source magnitude caps',
                'scope': 'conditional registered reduced W0 residual and first jets on alpha<=9/8; no full Einstein residual, nonlinear solution or physical observable'}
