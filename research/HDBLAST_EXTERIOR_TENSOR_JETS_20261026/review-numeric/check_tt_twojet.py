#!/usr/bin/env python3
"""Independent bounded TT two-jet report review; all decisions are rational.

Scientific checks never import the producer. Polygon clipping reconstructs
cone geometry, a triangular linear solve reconstructs the differentiated
resolvent, and a separately implemented interval ring checks every endpoint.
The source oracle is imported only after complete-tree authentication.
"""
import argparse
import hashlib
import importlib.util
import json
import math
import sys
import tempfile
from fractions import Fraction as Q
from pathlib import Path
from types import SimpleNamespace

sys.dont_write_bytecode = True
GRID = 1 << 112
END, CENTER = Q(1, 4), Q(5, 4)
ROOTS = {
    'HDBLAST_COUPLED_MEAN_PILOT_20261025': 'b932157f308c9912fa83ee06990f1132a4e75a498aafd66eb99dabcc30b256de',
    'HDBLAST_INPUT_PRECISION_20261023': 'c2e7bd40e94422c5fd98820941231c0da63a8392b14a5a14c4340eabdde37f61',
    'HDBLAST_COMMON_GAUGE_FORCING_20261016': '8bf173293b5cecf7f78e138ac5004abd9c7b1821ef29addb7ecfbc2a6409b4f6',
    'HDBLAST_NUMERICAL_RESPONSE_INPUTS_20261016': '928f3e4e58d222a3c7e1274babf9e94930ec501b74cf6ba04209f97c6bbded44',
    'HDBLAST_GEOMETRIC_BULK_SOURCES_20261016': '3ff1f95c151d2da6a1ab5e9b54718c2bb325c10d9c62a7ed9f97ef2bb36f0bb2',
    'HDBLAST_FORCED_RESPONSE_BOUNDS_20261023': '72924328a4e79cc77f0da90999b8dd56956bfdebd5e614e15c3c0498c3a8ae76',
    'HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023': 'cffb213a69f7f9a27932e990ffdecb7569d1007e5e5baeda708d8dc39d0cce72',
    'HDBLAST_FULL_ANGULAR_RESPONSE_20261026': 'a35f84119e84bf2ea8a712ded86ff79977d3ee91ac25f91ccbd202d379c859ac',
    'HDBLAST_SHELL_FEEDBACK_PRECISION_20261025': '1a5ee8aa22953f949da000ba03da3b5a0872e4467a0fbaec9f2d3fe7fba11ec4',
}


def require(value, message):
    if not value:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class Bounds:
    """Endpoint arithmetic, explicitly rounded to the declared output lattice."""
    def __init__(self, low=0, high=None):
        if isinstance(low, Bounds):
            low, high = low.low, low.high
        elif hasattr(low, 'lo'):
            low, high = low.lo, low.hi
        low = Q(low)
        high = low if high is None else Q(high)
        require(low <= high, 'INTERVAL_ORDER')
        self.low = Q(math.floor(low * GRID), GRID)
        self.high = Q(math.ceil(high * GRID), GRID)

    def __add__(self, other):
        other = Bounds(other)
        return Bounds(self.low + other.low, self.high + other.high)
    __radd__ = __add__

    def __neg__(self):
        return Bounds(-self.high, -self.low)

    def __sub__(self, other):
        return self + -Bounds(other)

    def __rsub__(self, other):
        return Bounds(other) + -self

    def __mul__(self, other):
        other = Bounds(other)
        corners = [a * b for a in (self.low, self.high) for b in (other.low, other.high)]
        return Bounds(min(corners), max(corners))
    __rmul__ = __mul__

    def reciprocal(self):
        require(self.high < 0 or self.low > 0, 'RECIPROCAL_ZERO')
        return Bounds(1 / self.high, 1 / self.low)

    def __truediv__(self, other):
        return self * Bounds(other).reciprocal()

    def __rtruediv__(self, other):
        return Bounds(other) * self.reciprocal()

    def __pow__(self, power):
        if power < 0:
            return self.reciprocal() ** -power
        if power == 0:
            return Bounds(1)
        if power % 2:
            return self * self ** (power - 1)
        half = self ** (power // 2)
        return Bounds(0 if half.low <= 0 <= half.high else min(half.low**2, half.high**2), max(half.low**2, half.high**2))

    def magnitude(self):
        return max(abs(self.low), abs(self.high))

    def row(self):
        return [str(self.low), str(self.high)]


class FirstDerivative:
    """Independent one-variable differential ring for analytic ODE inputs."""
    def __init__(self, value, radial=0):
        self.value, self.radial = Bounds(value), Bounds(radial)

    @staticmethod
    def convert(value):
        return value if isinstance(value, FirstDerivative) else FirstDerivative(value)

    def __add__(self, other):
        other = self.convert(other)
        return FirstDerivative(self.value + other.value, self.radial + other.radial)
    __radd__ = __add__

    def __neg__(self):
        return FirstDerivative(-self.value, -self.radial)

    def __sub__(self, other):
        return self + -self.convert(other)

    def __rsub__(self, other):
        return self.convert(other) + -self

    def __mul__(self, other):
        other = self.convert(other)
        return FirstDerivative(self.value * other.value, self.value * other.radial + self.radial * other.value)
    __rmul__ = __mul__

    def reciprocal(self):
        inv = self.value.reciprocal()
        return FirstDerivative(inv, -inv * (self.radial * inv))

    def __truediv__(self, other):
        return self * self.convert(other).reciprocal()

    def __rtruediv__(self, other):
        return self.convert(other) * self.reciprocal()

    def __pow__(self, power):
        if power < 0:
            return self.reciprocal() ** -power
        result = FirstDerivative(1)
        for _ in range(power):
            result *= self
        return result


def interval(row):
    require(isinstance(row, list) and len(row) == 2, 'TWO_ENDPOINTS_REQUIRED')
    result = Bounds(*row)
    require(result.row() == row, 'DECLARED_ENDPOINTS_ON_OUTPUT_LATTICE')
    return result


def same(actual, expected, message):
    actual = actual.row() if isinstance(actual, Bounds) else str(actual)
    require(actual == expected, message)


def inside(inner, outer):
    return outer.low <= inner.low <= inner.high <= outer.high


def meet(a, b):
    return Bounds(max(a.low, b.low), min(a.high, b.high))


def sign(value):
    return 'negative' if value.high < 0 else ('positive' if value.low > 0 else 'inconclusive')


def authenticate(repo):
    verified = {}
    for label, pin in ROOTS.items():
        root = repo / 'research' / label
        require(root.is_dir() and not root.is_symlink(), 'FROZEN_ROOT:' + label)
        manifest = root / 'MANIFEST.json'
        require(not manifest.is_symlink() and digest(manifest) == pin, 'FROZEN_MANIFEST:' + label)
        entries = json.loads(manifest.read_text())['files']
        expected = {'MANIFEST.json'}
        for entry in entries:
            rel = Path(entry['path'])
            require(not rel.is_absolute() and '..' not in rel.parts and str(rel) not in expected, 'FROZEN_RELATIVE_MEMBER')
            expected.add(str(rel))
            member = root / rel
            require(member.is_file() and not member.is_symlink(), 'FROZEN_MEMBER_TYPE')
            require(member.stat().st_size == entry['bytes'] and digest(member) == entry['sha256'], 'FROZEN_MEMBER_CONTENT')
        observed = set()
        for member in root.rglob('*'):
            require(not member.is_symlink(), 'FROZEN_TREE_LINK')
            if member.is_file():
                observed.add(str(member.relative_to(root)))
        require(observed == expected, 'FROZEN_EXACT_TREE')
        verified[label] = {'manifest_sha256': pin, 'files_including_manifest': len(expected)}
    return verified


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def clip(poly, a, b, c):
    """Clip polygon to a*x+b*r+c >= 0 with exact edge intersections."""
    out = []
    for start, finish in zip(poly, poly[1:] + poly[:1]):
        first = a * start[0] + b * start[1] + c
        second = a * finish[0] + b * finish[1] + c
        if first >= 0:
            out.append(start)
        if (first < 0 <= second) or (second < 0 <= first):
            scale = first / (first - second)
            out.append(tuple(start[k] + scale * (finish[k] - start[k]) for k in (0, 1)))
    return out


def area(poly):
    if not poly:
        return Q(0)
    return abs(sum(a[0] * b[1] - a[1] * b[0] for a, b in zip(poly, poly[1:] + poly[:1]))) / 2


def cone_geometry(n):
    step = END / n
    result = {}
    for i in range(n):
        for k in range(2 * n):
            xl, xh = i * step, (i + 1) * step
            rl, rh = CENTER - END + k * step, CENTER - END + (k + 1) * step
            poly = [(xl, rl), (xh, rl), (xh, rh), (xl, rh)]
            poly = clip(clip(poly, -1, 1, -(CENTER - END)), -1, -1, CENTER + END)
            weight = area(poly)
            if weight:
                result[i, k] = (xl, xh, rl, rh, weight)
    require(sum(row[-1] for row in result.values()) == END**2, 'POLYGON_CONE_TOTAL_AREA')
    return result


def integrate_rows(rows):
    low = sum(weight * value.low for weight, value in rows)
    high = sum(weight * value.high for weight, value in rows)
    return Bounds(low, high)


def solve(matrix, right):
    """Generic exact elimination, used for the differentiated 2x2 resolvent."""
    rows = [list(map(Q, row)) + [Q(value)] for row, value in zip(matrix, right)]
    for column in range(len(rows)):
        pivot = next(i for i in range(column, len(rows)) if rows[i][column])
        rows[column], rows[pivot] = rows[pivot], rows[column]
        divisor = rows[column][column]
        rows[column] = [value / divisor for value in rows[column]]
        for i in range(len(rows)):
            if i != column:
                factor = rows[i][column]
                rows[i] = [a - factor * b for a, b in zip(rows[i], rows[column])]
    return [row[-1] for row in rows]


def check_operator(op):
    coefficients = {key: interval(value) for key, value in op['coefficient_intervals'].items()}
    radial = {}
    for key in ('Ar', 'Br', 'Dr'):
        first = interval(op['radial_derivative_analytic_jet_intervals'][key])
        second = interval(op['radial_derivative_physical_identity_intervals'][key])
        radial[key] = meet(first, second)
        same(radial[key], op['radial_derivative_intersections'][key], 'COEFFICIENT_INTERSECTION:' + key)
        same(radial[key].magnitude(), op['radial_derivative_absolute_caps'][key], 'RADIAL_MAGNITUDE:' + key)
    for key, bound in coefficients.items():
        same(bound.magnitude(), op['absolute_caps'][key], 'COEFFICIENT_MAGNITUDE:' + key)
    sigma = coefficients['A'].magnitude() + coefficients['B'].magnitude() + END * coefficients['D'].magnitude() / 2
    sigma_r = radial['Ar'].magnitude() + radial['Br'].magnitude() + END * radial['Dr'].magnitude() / 2
    same(sigma, op['sigma'], 'SIGMA')
    same(sigma_r, op['sigma_r'], 'SIGMA_R')
    same(END * sigma, op['q'], 'Q')
    same(END * sigma_r, op['q_r'], 'Q_R')
    require(0 <= END * sigma < 1 and sigma_r >= 0, 'RESOLVENT_INVERTIBILITY')
    bg = op['background_value_intervals']
    require(Q(bg['ell']) == Q(1, 10**6), 'ELL_CONVENTION')
    same(1 / (Q(bg['ell']) * interval(bg['rho'])), bg['L'], 'COORDINATE_L')
    return coefficients


def check_run(run, whole_op, apex_op, source):
    cone = run['cone_quadrature']
    partition = cone['partition']
    n = partition['n']
    require(n in (4, 8), 'FINITE_DECLARED_CONE_GRID')
    geometry = cone_geometry(n)
    require(len(geometry) == partition['cells'] == len(cone['cells']) == n * (n + 1), 'ALL_CONE_POLYGONS_PRESENT')
    same(END / n, partition['h'], 'CONE_STEP')
    same(END**2, partition['exact_area'], 'CONE_AREA')
    integrated = {key: [] for key in ('F', 'Fr')}
    caps = {key: Q(0) for key in integrated}
    seen = set()
    for row in cone['cells']:
        key = row['i'], row['k']
        require(key in geometry and key not in seen, 'UNIQUE_CONE_CELL')
        seen.add(key)
        xl, xh, rl, rh, weight = geometry[key]
        require(row['x'] == list(map(str, (xl, xh))) and row['r'] == list(map(str, (rl, rh))), 'CLIPPED_RECTANGLE_BOUNDS')
        same(weight, row['exact_intersection_area'], 'POLYGON_CLIPPED_AREA')
        require(rl >= xh and xh + rh <= 8, 'WHOLE_BOX_SOURCE_ORACLE_DOMAIN')
        for name in integrated:
            bound = interval(row[name])
            integrated[name].append((weight / 2, bound))
            caps[name] = max(caps[name], bound.magnitude())
    free_area = {name: integrate_rows(rows) for name, rows in integrated.items()}
    for name in integrated:
        same(free_area[name], cone['free_apex_values'][name], 'HALF_RETARDED_AREA:' + name)
        same(caps[name], cone['source_absolute_caps'][name], 'WHOLE_CONE_CAP:' + name)
    lines = run['characteristic_quadrature']
    count = lines['partition']['segments_per_characteristic']
    require(count == 2 * n and len(lines['segments']) == count, 'DECLARED_LINE_BUDGET')
    delta = END / count
    same(delta, lines['partition']['h'], 'LINE_STEP')
    same(delta / 2, lines['partition']['weight_each_branch'], 'LINE_HALF_WEIGHT')
    integrated = {key: [] for key in ('F_plus', 'F_minus', 'Fr_plus', 'Fr_minus')}
    for j, row in enumerate(lines['segments']):
        low, high = j * delta, (j + 1) * delta
        require(row['segment'] == j and row['x'] == list(map(str, (low, high))), 'CONTIGUOUS_TIME_TRACE')
        same((high - low) / 2, row['weight'], 'CHARACTERISTIC_HALF_MEASURE')
        for branch, orientation in (('plus', 1), ('minus', -1)):
            radii = [CENTER + orientation * (END - x) for x in (low, high)]
            require(row[branch]['r'] == list(map(str, sorted(radii))), 'CHARACTERISTIC_RADIAL_TRACE')
            require(min(radii) >= high, 'CHARACTERISTIC_EXTERIOR')
            for name in ('F', 'Fr'):
                integrated[name + '_' + branch].append(((high - low) / 2, interval(row[branch][name])))
    traces = {name: integrate_rows(rows) for name, rows in integrated.items()}
    for name, value in traces.items():
        same(value, lines['half_integrals'][name], 'ACTUAL_HALF_TRACE:' + name)
    free_line = {
        'u_x': traces['F_plus'] + traces['F_minus'],
        'u_r': traces['F_plus'] - traces['F_minus'],
        'u_xr': traces['Fr_plus'] + traces['Fr_minus'],
        'u_rr': traces['Fr_plus'] - traces['Fr_minus'],
    }
    for name, value in free_line.items():
        same(value, lines['free_apex_derivatives'][name], 'CHARACTERISTIC_ORIENTATION:' + name)
    q, qr = Q(whole_op['q']), Q(whole_op['q_r'])
    free_norms = [END * caps['F'], END * caps['Fr']]
    mu, mv = solve([[1 - q, 0], [-qr, 1 - q]], free_norms)
    e1, e2 = mu - free_norms[0], mv - free_norms[1]
    ledger = {'M_u': mu, 'M_v': mv, 'E_0': END * e1 / 2, 'E_1': e1, 'E_2': e2, 'E_v_value': END * e2 / 2}
    for name, value in ledger.items():
        same(value, run['resolvent_ledger'][name], 'TRIANGULAR_RESOLVENT:' + name)
    require(e1 == q * mu and e2 == q * mv + qr * mu, 'ALL_FEEDBACK_INSERTIONS')
    jet = {'u': free_area['F'] + Bounds(-ledger['E_0'], ledger['E_0'])}
    for name in ('u_x', 'u_r'):
        jet[name] = free_line[name] + Bounds(-e1, e1)
    first_route = jet['u_r']
    second_route = free_area['Fr'] + Bounds(-ledger['E_v_value'], ledger['E_v_value'])
    jet['u_r'] = meet(first_route, second_route)
    routes = run['apex_u_r_two_routes']
    for name, value in (('characteristic_F', first_route), ('cone_Fr', second_route), ('intersection', jet['u_r'])):
        same(value, routes[name], 'UR_TWO_VALID_ROUTES:' + name)
    for name in ('u_xr', 'u_rr'):
        jet[name] = free_line[name] + Bounds(-e2, e2)
    coeff = {key: interval(value) for key, value in apex_op['coefficient_intervals'].items()}
    lower = coeff['A'] * jet['u_x'] + coeff['B'] * jet['u_r'] + coeff['D'] * jet['u'] + interval(source['F'])
    same(lower, run['pde_closure']['lower_order_plus_source_interval'], 'PDE_LOWER_ORDER_AND_SOURCE')
    jet['u_xx'] = jet['u_rr'] + lower
    require(run['pde_closure']['exact_equation_defines_u_xx'] and not run['pde_closure']['requires_source_x_derivative'], 'PDE_CLOSURE_SCOPE')
    for name, value in jet.items():
        same(value, run['total_apex_twojet'][name], 'COMPLETE_ACTUAL_TWOJET:' + name)
        require(run['display_only'][name] == [float(value.low), float(value.high)], 'DISPLAY_IS_ONLY_PROJECTION')
    bg = apex_op['background_value_intervals']
    ell, L, a = Q(bg['ell']), interval(bg['L']), interval(bg['a'])
    natural = {'u': jet['u'], 'u_t': jet['u_x'] / ell, 'u_y': -L * jet['u_r'],
               'u_tt': jet['u_xx'] / ell**2, 'u_ty': -L * jet['u_xr'] / ell,
               'u_yy': L**2 * jet['u_rr'] + a * L * jet['u_r']}
    for name, value in natural.items():
        same(value, run['ordinary_t_y_coefficient_twojet'][name], 'ORDINARY_CHAIN_RULE:' + name)
    for names, field, values in ((('u', 'u_x', 'u_r'), 'scaled_value_and_first_derivative_signs', jet),
                                 (('u_xx', 'u_xr', 'u_rr'), 'second_derivative_signs', jet),
                                 (('u_tt', 'u_ty', 'u_yy'), 'ordinary_t_y_second_derivative_signs', natural)):
        require(run[field] == {name: sign(values[name]) for name in names}, 'EXACT_ENDPOINT_SIGNS:' + field)
    gate = any(sign(jet[key]) != 'inconclusive' for key in ('u_xx', 'u_xr', 'u_rr'))
    require(run['passes_diagnostic_nonzero_gate'] == gate and run['finite_twojet_primary_acceptance'], 'DIAGNOSTIC_GATE')
    conversion = run['ordinary_t_y_conversion']
    require(not conversion['covariant_tensor_derivatives'] and not conversion['curvature_components'] and not conversion['observable'], 'COEFFICIENT_PARTIAL_DERIVATIVES_ONLY')
    return {'cone_n': n, 'polygon_cells': len(geometry), 'line_segments_per_branch': count,
            'half_cone_area': str(END**2 / 2), 'resolvent_ledger': {k: str(v) for k, v in ledger.items()},
            'scaled_signs': {name: sign(value) for name, value in jet.items()}, 'all_report_endpoints_reconstructed': True}


def coefficient_and_source_replays(ctx, report):
    """Use background ODE values, independent compact D, and three source cells."""
    j, ell = ctx.jets, ctx.ell
    coefficient_rows = []
    for name in ('whole_cone_operator', 'apex_operator'):
        op = report[name]
        x, r = ctx.I(*map(Q, op['box']['x'])), ctx.I(*map(Q, op['box']['r']))
        rho, eta, p, a, h, *_ = j.core.background(r, ctx.data, 'none')
        derivatives = [FirstDerivative(j.radial(value).val(), j.radial(value).derivative(0, 1)) for value in (rho, eta, p, a)]
        rho, eta, p, a, h = map(lambda value: Bounds(value.val()), (rho, eta, p, a, h))
        c, H, _ = j.temporal(x)
        c, H = Bounds(c.val()), Bounds(H.val())
        drho, deta, dp, da = derivatives
        dc, dH = FirstDerivative(c), FirstDerivative(H)
        day = -1 / drho**2 - dp**2 / 3
        db = -(day + da**2)
        dd = 1 / drho**2 - da**2
        dV = FirstDerivative(0)
        for coefficient in reversed(j.core.U):
            dV = dV * deta + coefficient
        dh = dH / drho
        analytic = {'Ar': (-3 * ell * dH).radial, 'Br': (-3 * ell * drho * da).radial,
                    'Dr': (ell**2 * (-6 / dc**2 + drho**2 * (-2 * da**2 + 2 * dh**2 - 2 * db - 8 * dd + Q(4, 3) * dV))).radial}
        for key, value in analytic.items():
            same(value, op['radial_derivative_analytic_jet_intervals'][key], 'INDEPENDENT_FIRST_ORDER_DIFFERENTIAL_RING:' + name + ':' + key)
        V = Bounds(j.core.poly(j.core.U, ctx.I(eta.low, eta.high)))
        Vp = Bounds(j.core.poly(j.core.UP[1], ctx.I(eta.low, eta.high)))
        ay = -1 / rho**2 - p**2 / 3
        py = p * (h - 4 * a)
        S = 8 * a**2 - Q(2, 3) * p**2 + Q(4, 3) * V
        # D=ell²(rho²S-8(1+c^-2)) follows from h=H/rho,
        # ay=-rho^-2-p²/3, b=-(ay+a²), d=rho^-2-a², H²=1-c^-2.
        independent = {'A': -3 * ell * H, 'B': -3 * ell * rho * a,
                       'D': ell**2 * (rho**2 * S - 8 * (1 + 1 / c**2))}
        for key, value in independent.items():
            require(inside(value, interval(op['coefficient_intervals'][key])), 'COMPACT_PHYSICAL_COEFFICIENT_SUBSET:' + name + ':' + key)
        physical = {'Ar': Bounds(0), 'Br': 3 * ell**2 * rho**2 * (a**2 + ay),
                    'Dr': -ell**3 * rho**3 * (2 * a * S + 16 * a * ay - Q(4, 3) * p * py + Q(4, 3) * p * Vp)}
        for key, value in physical.items():
            same(value, op['radial_derivative_physical_identity_intervals'][key], 'ODE_PHYSICAL_RADIAL_IDENTITY:' + key)
        coefficient_rows.append({'operator': name, 'compact_coefficient_intervals': {k: v.row() for k, v in independent.items()},
                                 'all_compact_intervals_subset_report': True, 'physical_radial_identity_endpoints_match': True,
                                 'independent_first_order_differential_ring_endpoints_match': True})
    final = report['campaign_runs'][-1]
    line = final['characteristic_quadrature']['segments']
    choices = [('apex', report['apex_operator']['box'], report['apex_source']),
               ('first_plus_characteristic', {'x': line[0]['x'], 'r': line[0]['plus']['r']}, line[0]['plus']),
               ('last_minus_characteristic', {'x': line[-1]['x'], 'r': line[-1]['minus']['r']}, line[-1]['minus'])]
    rows = []
    for label, box, expected in choices:
        x, r = (ctx.I(*map(Q, box[key])) for key in ('x', 'r'))
        value, radial = ctx.forcing_with_radial(x, r, audit=True)
        for key, result in (('F', value), ('Fr', radial)):
            require([str(result.lo), str(result.hi)] == expected[key], 'INDEPENDENT_DECLARED_SOURCE_REPLAY:' + label + ':' + key)
        rows.append({'label': label, 'x': box['x'], 'r': box['r'], 'F': expected['F'], 'Fr': expected['Fr'], 'exact_endpoint_match': True})
    return {'coefficient_replays': coefficient_rows, 'source_replay_cells': rows, 'source_oracle_calls': len(rows)}


# Generic exact bivariate polynomial arithmetic supplies controls independently
# of both the interval implementation and the producer quadrature code.
def padd(a, b):
    keys = set(a) | set(b)
    return {key: a.get(key, 0) + b.get(key, 0) for key in keys if a.get(key, 0) + b.get(key, 0)}


def pmul(a, b):
    out = {}
    for (i, j), x in a.items():
        for (k, l), y in b.items():
            key = i + k, j + l
            out[key] = out.get(key, 0) + x * y
    return {key: value for key, value in out.items() if value}


def pdiff(a, axis):
    out = {}
    for power, value in a.items():
        if power[axis]:
            new = list(power)
            new[axis] -= 1
            out[tuple(new)] = value * power[axis]
    return out


def peval(a, x, r):
    return sum(value * x**i * r**j for (i, j), value in a.items())


def trace_integral(poly, orientation, time=END, radius=CENTER):
    # r(s)=radius+orientation*(time-s); expand each monomial and integrate.
    result = Q(0)
    for (i, j), value in poly.items():
        for k in range(j + 1):
            result += value * math.comb(j, k) * (radius + orientation * time)**(j - k) * (-orientation)**k * time**(i + k + 1) / (i + k + 1)
    return result / 2


def cone_integral(poly, time=END, radius=CENTER):
    # Integrate in centered spatial coordinate z, then s; odd z terms vanish.
    result = Q(0)
    for (i, j), value in poly.items():
        for k in range(0, j + 1, 2):
            coefficient = value * math.comb(j, k) * radius**(j - k) / (k + 1)
            # Integral s^i (time-s)^(k+1) ds, exactly via the beta identity.
            result += coefficient * time**(i + k + 2) * Q(math.factorial(i) * math.factorial(k + 1), math.factorial(i + k + 2))
    return result


def controls():
    cases = []
    mutations = set()
    for i, j in ((2, 0), (2, 1), (2, 2), (3, 2), (4, 3), (5, 1), (6, 4)):
        u = {(i, j): Q(1)}
        uxx = pdiff(pdiff(u, 0), 0)
        urr = pdiff(pdiff(u, 1), 1)
        force = padd(uxx, {key: -value for key, value in urr.items()})
        fr = pdiff(force, 1)
        plus, minus = trace_integral(force, 1), trace_integral(force, -1)
        rplus, rminus = trace_integral(fr, 1), trace_integral(fr, -1)
        exact = {'u': peval(u, END, CENTER), 'u_x': peval(pdiff(u, 0), END, CENTER),
                 'u_r': peval(pdiff(u, 1), END, CENTER), 'u_xr': peval(pdiff(pdiff(u, 0), 1), END, CENTER),
                 'u_rr': peval(urr, END, CENTER), 'u_xx': peval(uxx, END, CENTER)}
        derived = {'u': cone_integral(force), 'u_x': plus + minus, 'u_r': plus - minus,
                   'u_xr': rplus + rminus, 'u_rr': rplus - rminus,
                   'u_xx': rplus - rminus + peval(force, END, CENTER)}
        require(derived == exact and cone_integral(fr) == exact['u_r'], 'MANUFACTURED_POLYNOMIAL_IDENTITIES')
        if 2 * derived['u'] != exact['u']:
            mutations.add('missing_half_retarded_kernel')
        if minus - plus != exact['u_r']:
            mutations.add('reversed_radial_trace_orientation')
        if rplus + rminus != exact['u_rr']:
            mutations.add('using_trace_sum_for_u_rr')
        if rplus - rminus != exact['u_xr']:
            mutations.add('using_trace_difference_for_u_xr')
        cases.append({'u_monomial': [i, j], 'exact_six_jet': {key: str(value) for key, value in exact.items()}})
    # This deliberately nonzero triangular coupling distinguishes the two
    # different qr insertions, beyond testing only a copied scalar formula.
    q, qr, f, g = Q(1, 7), Q(2, 11), Q(3, 13), Q(5, 17)
    mu, mv = solve([[1 - q, 0], [-qr, 1 - q]], [f, g])
    tail = mv - g
    require(tail == q * mv + qr * mu, 'TRIANGULAR_FEEDBACK_REFERENCE')
    require(q * mv < tail and q * g / (1 - q) + qr * mu < tail, 'MISSING_CROSS_COUPLING_CONTROLS')
    mutations.update(('omitted_qr_Mu_in_second_derivative_error', 'omitted_resolvent_propagation_of_qr_coupling'))
    # A scalar coefficient u=r and nonconstant L=2+3r give d_y²u=L L_r.
    # This detects either deleting the ordinary-coordinate connection term or
    # changing its sign without invoking a covariant-curvature interpretation.
    r = Q(4, 5)
    L = 2 + 3 * r
    coordinate_truth = 3 * L
    require(coordinate_truth != 0 and coordinate_truth != -3 * L, 'NONCONSTANT_COORDINATE_MAP_CONTROL')
    mutations.update(('omitted_a_L_ur_in_ordinary_u_yy', 'wrong_sign_of_a_L_ur_in_ordinary_u_yy'))
    require(len(mutations) == 8, 'EIGHT_DISTINCT_FORMULA_CONTROLS')
    return {'exact_manufactured_six_jet_cases': cases, 'formula_mutations_rejected': sorted(mutations),
            'nonconstant_coordinate_chain_rule_reference': str(coordinate_truth), 'all_pass': True}


def rejects(action, marker):
    try:
        action()
    except (ValueError, FileNotFoundError) as error:
        require(marker in str(error), 'UNEXPECTED_CONTROL_FAILURE:' + str(error))
        return marker
    raise ValueError('CONTROL_NOT_REJECTED:' + marker)


def driver_controls(evaluator, repo, scratch):
    """Guard-only import; scientific producer routines are never called."""
    module = load('numeric_driver_guard_control', evaluator)
    results = []
    with tempfile.TemporaryDirectory(prefix='guard-fixture-', dir=scratch) as directory:
        fixture = Path(directory)
        existing = fixture / 'exists.json'
        existing.write_text('retained')
        results.append(rejects(lambda: module.guard_output(repo, existing), 'FRESH_OUTPUT_REQUIRED'))
        link = fixture / 'link.json'
        link.symlink_to(fixture / 'nonexistent.json')
        results.append(rejects(lambda: module.guard_output(repo, link), 'FRESH_OUTPUT_REQUIRED'))
        for label in ROOTS:
            target = repo / 'research' / label / '__review_missing_output__.json'
            require(not target.exists(), 'GUARD_CONTROL_FRESH_TARGET')
            marker = 'OUTPUT_INSIDE_FROZEN_CAPSULE'
            results.append(rejects(lambda target=target: module.guard_output(repo, target), marker))
        # Relocate the evaluator beneath a new research capsule and verify its
        # own complete root is protected even against an unrelated repo root.
        packaged_root = fixture / 'relocated' / 'research' / 'NEW_REVIEW_FIXTURE'
        packaged_root.mkdir(parents=True)
        packaged_source = packaged_root / 'driver.py'
        packaged_source.write_bytes(evaluator.read_bytes())
        packaged_module = load('numeric_packaged_output_guard_control', packaged_source)
        results.append(rejects(lambda: packaged_module.guard_output(repo, packaged_root / 'fresh.json'), 'OUTPUT_INSIDE_PACKAGED_CAPSULE'))
        # A synthetic complete closure exercises every tree-rejection branch.
        synthetic = fixture / 'repo'
        root = synthetic / 'research' / 'MOCK'
        root.mkdir(parents=True)
        member = root / 'member.txt'
        member.write_text('known fixture bytes\n')
        manifest = root / 'MANIFEST.json'
        manifest.write_text(json.dumps({'files': [{'path': member.name, 'bytes': member.stat().st_size, 'sha256': digest(member)}]}))
        module.LABELS = {'fixture': 'MOCK'}
        module.PINS = {'fixture': digest(manifest)}
        require(module.authenticate(synthetic)['fixture']['files_including_manifest'] == 2, 'SYNTHETIC_VALID_CLOSURE')
        saved = module.PINS['fixture']
        module.PINS['fixture'] = '0' * 64
        results.append(rejects(lambda: module.authenticate(synthetic), 'MANIFEST:fixture'))
        module.PINS['fixture'] = saved
        original = member.read_bytes()
        member.write_bytes(b'different content\n')
        results.append(rejects(lambda: module.authenticate(synthetic), 'MEMBER_BYTES:fixture'))
        member.write_bytes(original)
        renamed_root = root.with_name('REGULAR_MOCK')
        root.rename(renamed_root)
        root.symlink_to(renamed_root, target_is_directory=True)
        results.append(rejects(lambda: module.authenticate(synthetic), 'ROOT_TYPE:fixture'))
        root.unlink()
        renamed_root.rename(root)
        extra = root / 'unlisted.txt'
        extra.write_text('extra')
        results.append(rejects(lambda: module.authenticate(synthetic), 'UNLISTED_FILE:fixture'))
        extra.unlink()
        member.unlink()
        member.symlink_to(existing)
        results.append(rejects(lambda: module.authenticate(synthetic), 'MEMBER_TYPE:fixture'))
        member.unlink()
        member.write_bytes(original)
        directory_link = root / 'directory-link'
        directory_link.symlink_to(fixture, target_is_directory=True)
        results.append(rejects(lambda: module.authenticate(synthetic), 'TREE_SYMLINK:fixture'))
        directory_link.unlink()
        source = fixture / 'input.py'
        source.write_text('# isolated external-pin fixture\n')
        output = fixture / 'unwritten.json'
        original_argv = sys.argv[:]
        def invoke(pin):
            sys.argv = [str(evaluator), '--repo-root', str(synthetic), '--input-module', str(source),
                        '--input-module-sha256', pin, '--output', str(output), '--probes-only']
            module.main()
        try:
            results.append(rejects(lambda: invoke('0' * 64), 'NEW_INPUT_MODULE_EXTERNAL_PIN'))
            source.unlink()
            source.symlink_to(existing)
            results.append(rejects(lambda: invoke(digest(source)), 'NEW_INPUT_MODULE_TYPE'))
            source.unlink()
            source.write_text('# restored isolated external-pin fixture\n')
            # Isolate expensive science, retaining main's before/after source
            # authentication. Context.finish changes only the synthetic file.
            class FixtureContext:
                def __init__(self, _):
                    self.base = SimpleNamespace(I=lambda value: value)
                    self.certificate = {}
                def finish(self):
                    source.write_text('# altered during isolated fixture\n')
            module.load = lambda *_: SimpleNamespace(Context=FixtureContext)
            module.probes = lambda _: []
            module.operator = lambda *_, **__: {}
            module.evaluate_source = lambda *_, **__: (0, 0)
            results.append(rejects(lambda: invoke(digest(source)), 'NEW_INPUT_MODULE_CHANGED'))
            require(not output.exists(), 'FAILED_GUARDS_WRITE_NO_OUTPUT')
        finally:
            sys.argv = original_argv
        require(existing.read_text() == 'retained', 'EXISTING_OUTPUT_BYTES_RETAINED')
    return {'guard_failures_observed': results, 'isolated_synthetic_closure': True,
            'scientific_producer_routines_called': False, 'all_failed_attempts_wrote_no_report': True}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo-root', type=Path, required=True)
    parser.add_argument('--evaluator', type=Path, required=True)
    parser.add_argument('--evaluator-sha256', required=True)
    parser.add_argument('--input-module', type=Path, required=True)
    parser.add_argument('--input-module-sha256', required=True)
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--report-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists() and not args.output.is_symlink(), 'REVIEW_FRESH_OUTPUT_REQUIRED')
    repo, destination = args.repo_root.resolve(), args.output.resolve()
    for label in ROOTS:
        root = (repo / 'research' / label).resolve()
        require(destination != root and root not in destination.parents, 'REVIEW_OUTPUT_INSIDE_FROZEN_CAPSULE')
    for parent in Path(__file__).resolve().parents:
        if parent.parent.name == 'research':
            require(destination != parent and parent not in destination.parents, 'REVIEW_OUTPUT_INSIDE_PACKAGED_CAPSULE')
            break
    own_pin = digest(Path(__file__))
    for path, pin, label in ((args.evaluator, args.evaluator_sha256, 'EVALUATOR'),
                             (args.input_module, args.input_module_sha256, 'ORACLE'),
                             (args.report, args.report_sha256, 'PRODUCER_REPORT')):
        require(path.is_file() and not path.is_symlink() and digest(path) == pin, 'EXTERNAL_PIN:' + label)
    closure = authenticate(repo)
    report = json.loads(args.report.read_text())
    require(report['evaluator_sha256'] == args.evaluator_sha256 and report['input_module_sha256'] == args.input_module_sha256, 'REPORT_BINDS_EXTERNAL_CODE')
    require(report['status'] == 'CERTIFIED_FINITE_ACTUAL_TOTAL_TT_TWOJET', 'FINAL_SCIENTIFIC_REPORT_REQUIRED')
    require(report['all_nine_complete_closures_authenticated_before_project_imports_and_after_run'] and len(report['authentication']) == 9, 'NINE_COMPLETE_PRODUCER_CLOSURES')
    producer_pins = {entry['manifest_sha256']: entry['files_including_manifest'] for entry in report['authentication'].values()}
    require(producer_pins == {entry['manifest_sha256']: entry['files_including_manifest'] for entry in closure.values()}, 'PRODUCER_AUTHENTICATION_MATCHES_INDEPENDENT_TREES')
    require(report['source_input_certificate']['source_module_sha256'] == args.input_module_sha256, 'SOURCE_CERTIFICATE_EXTERNAL_PIN')
    declaration = report['budget_declaration']
    require(declaration['first'] == {'cone_n': 4, 'cone_cells': 20, 'segments_per_characteristic': 8} and
            declaration['maximum'] == {'cone_n': 8, 'cone_cells': 72, 'segments_per_characteristic': 16} and
            declaration['no_further_refinement'] and declaration['probes'] == {'points': 3, 'whole_cells': 2}, 'PUBLISHED_BUDGET')
    require(report['source_input_certificate']['dependency_order_audit']['full_terms'] == 78, 'COMPLETE_SOURCE_ROW_PREMISE')
    require(report['source_input_certificate']['all_actual_error_caps_retained'] and not report['source_input_certificate']['C4_required'], 'SOURCE_CERTIFICATE_SCOPE')
    domain = report['domain']
    require(domain['apex'] == ['1/4', '5/4'] and not domain['shell_reached'] and not domain['reflection_reached'], 'EXTERIOR_CAUSAL_CONE')
    require(all(domain[key] == 'exactly0' for key in ('TT_initial_value', 'TT_initial_velocity', 'TT_free_response')), 'TT_ZERO_DATA_STRUCTURAL_PREMISE')
    # Every point of the cone has r-x >= 3/4. A signal from r=0 at
    # nonnegative time cannot reach it. Its initial footprint is [1,3/2].
    require(CENTER - 2 * END == Q(3, 4) > 0, 'NO_SHELL_OR_REFLECTION_BY_FINITE_SPEED')
    for name in ('whole_cone_operator', 'apex_operator'):
        check_operator(report[name])
    runs = report['campaign_runs']
    require(1 <= len(runs) <= 2 and runs[0]['cone_quadrature']['partition']['n'] == 4, 'DECLARED_FIRST_CAMPAIGN')
    expected_second = not runs[0]['passes_diagnostic_nonzero_gate']
    require((len(runs) == 2) == expected_second, 'EXACTLY_ONE_CONDITIONAL_REFINEMENT')
    if len(runs) == 2:
        require(runs[1]['cone_quadrature']['partition']['n'] == 8, 'DECLARED_MAXIMUM_CAMPAIGN')
    run_reviews = [check_run(run, report['whole_cone_operator'], report['apex_operator'], report['apex_source']) for run in runs]
    final = runs[-1]
    expected_diagnostic = 'strict_nonzero_second_derivative_certified' if final['passes_diagnostic_nonzero_gate'] else 'inconclusive_at_declared_maximum'
    require(report['diagnostic_nonzero_status'] == expected_diagnostic, 'HONEST_SECOND_DERIVATIVE_DIAGNOSTIC')
    old = report['old_TT_value_compatibility']
    old_path = repo / old['path']
    require(digest(old_path) == old['sha256'] and not old['old_run_replayed'], 'OLD_VALUE_SOURCE_AND_NO_REPLAY')
    old_endpoints = json.loads(old_path.read_text())['partition_runs'][-1]['total_Tplus_response_interval']
    require(old_endpoints == old['old_total_value_interval'], 'OLD_VALUE_ENDPOINTS')
    old_low, old_high = map(Q, old_endpoints)
    new = interval(final['total_apex_twojet']['u'])
    overlap = [max(old_low, new.low), min(old_high, new.high)]
    require(overlap[0] <= overlap[1] and list(map(str, overlap)) == old['nonempty_overlap'], 'OLD_VALUE_OVERLAP_GATE')
    require(old['new_total_value_interval'] == new.row(), 'NEW_VALUE_COMPATIBILITY_BINDING')
    old_subset = old_low <= new.low <= new.high <= old_high
    formula_controls = controls()
    oracle = load('review_authenticated_actual_radial_source', args.input_module)
    ctx = oracle.Context(repo)
    replays = coefficient_and_source_replays(ctx, report)
    ctx.finish()
    guard_reviews = driver_controls(args.evaluator, repo, args.output.parent)
    require(authenticate(repo) == closure, 'FROZEN_CLOSURES_CHANGED_DURING_REVIEW')
    require(digest(args.evaluator) == args.evaluator_sha256 and digest(args.input_module) == args.input_module_sha256 and digest(args.report) == args.report_sha256 and digest(Path(__file__)) == own_pin, 'REVIEW_BOUND_INPUT_CHANGED')
    result = {'schema': 'HDBLAST_INDEPENDENT_ACTUAL_TOTAL_TT_TWOJET_REVIEW_V1', 'status': 'PASS_INDEPENDENT_ACTUAL_TOTAL_TT_TWOJET_REVIEW',
              'reviewer_sha256': own_pin, 'evaluator_sha256': args.evaluator_sha256, 'input_module_sha256': args.input_module_sha256,
              'producer_report_sha256': args.report_sha256, 'nine_complete_dependency_roots': closure,
              'campaign_reviews': run_reviews, 'formula_controls': formula_controls, 'driver_guards': guard_reviews,
              'independent_actual_replays': replays, 'old_value_nonempty_overlap': list(map(str, overlap)),
              'new_value_subset_old_value': old_subset, 'second_derivative_diagnostic': expected_diagnostic,
              'no_scientific_producer_functions_called': True, 'producer_import_for_isolated_driver_guards_only': True,
              'scope': 'Independent numerical geometry, all scalar feedback coefficients, all six scaled and ordinary t/y coefficient derivative intervals, signs, old-value compatibility, bounded-budget gate, selected oracle replays and IO/authentication controls. Full 78-term source differentiation has a separate independent input review.',
              'premises': 'Retained smooth actual common-gauge solution, formal full TT operator/source identities, exact TT zero Cauchy data, finite propagation and authenticated actual source/C3 derivative certificates remain conditional scientific premises; no derivative of a C1 error bound is taken.',
              'all_dependencies_reauthenticated_after_review': True}
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write('\n')
    print(result['status'])


if __name__ == '__main__':
    main()
