#!/usr/bin/env python3
"""Frozen P1 homogeneous preparation-response enclosure, no oracle calls."""
import hashlib
import importlib.util
import json
import math
import sys
from fractions import Fraction as F
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ELL = F(1, 10**6)
T = F(5, 4)
DELTA = T - 1
Q_RESOLVENT = F(1201, 3200000)
ETA = F(1, 10**11)
GRAM = (1, 2, 1, 3, 1)
SOURCE_REL = Path('research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016')
CAP_REL = Path('research/HDBLAST_PREPARED_SHELL_VOLUME_20261101')


def need(ok, witness):
    if not ok:
        raise ValueError(witness)


def digest(p):
    raw = p.read_bytes()
    return {'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()}


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class I:
    """Exact endpoint interval arithmetic for aggregation of rounded raw boxes."""
    def __init__(self, lo=0, hi=None):
        self.lo = F(lo)
        self.hi = self.lo if hi is None else F(hi)
        need(self.lo <= self.hi, 'INTERVAL_ORDER')

    @staticmethod
    def cv(x):
        return x if isinstance(x, I) else I(x)

    def __add__(self, b):
        b = self.cv(b)
        return I(self.lo + b.lo, self.hi + b.hi)
    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, b):
        return self + -self.cv(b)

    def __rsub__(self, b):
        return self.cv(b) + -self

    def __mul__(self, b):
        b = self.cv(b)
        vals = [a*c for a in (self.lo, self.hi) for c in (b.lo, b.hi)]
        return I(min(vals), max(vals))
    __rmul__ = __mul__

    def __truediv__(self, b):
        b = self.cv(b)
        need(b.lo > 0 or b.hi < 0, 'ZERO_DENOMINATOR')
        return self * I(1/b.hi, 1/b.lo)

    def __pow__(self, n):
        need(isinstance(n, int) and n >= 0, 'NONNEGATIVE_POWER')
        if not n:
            return I(1)
        if n % 2 == 0:
            return I(0 if self.lo <= 0 <= self.hi else min(self.lo**n, self.hi**n),
                     max(self.lo**n, self.hi**n))
        return self * self**(n-1)

    def mag(self):
        return max(abs(self.lo), abs(self.hi))

    def row(self):
        return [str(self.lo), str(self.hi)]


def interval(row):
    return I(*row['exact'])


def enc(x):
    if isinstance(x, I) or hasattr(x, 'row'):
        return {'exact': x.row(), 'display_only': [float(x.lo), float(x.hi)]}
    if isinstance(x, F):
        return {'exact': str(x), 'display_only': float(x)}
    if isinstance(x, dict):
        return {k: enc(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [enc(v) for v in x]
    return x


def sqrtup(v):
    grid = 1 << 112
    scaled = v * grid * grid
    k = math.isqrt(scaled.numerator // scaled.denominator)
    if F(k*k) < scaled:
        k += 1
    return F(k, grid)


def gram_upper(component_caps):
    return sqrtup(sum((F(w)*x*x for w, x in zip(GRAM, component_caps)), F(0)))


def initial_profile(m, bg, radius, tail, bump_data, rho_b, py):
    """All outputs are fresh whole-cell outward boxes for this direction."""
    rho, p, a = [m.I(*bg[k]['exact']) for k in ('rho', 'p', 'a')]
    qjet = m.bump(m.I(*radius), {'bump': bump_data}, tail)
    # m.bump stores z Taylor coefficients; convert to dimensionless r derivatives.
    f = [qjet.c[n] * ELL**n * math.factorial(n) for n in range(4)]
    q0, q1, q2, q3 = f
    b = ELL*rho*a
    d = ELL*rho*py/p
    ay = -1/rho**2-p*p/3
    R = (rho_b/rho)**3
    c = -3*rho_b**3/(2*ELL*rho**4*p)
    H = q2+4*b*q1+3*ELL**2*q0
    h = -R*q1
    hr = -R*(q2+3*b*q1)
    vl = -R*(q2+4*b*q1)
    s = c*H
    sr = c*(q3+4*b*q2+3*ELL**2*q1+(4*b+d)*H
            -4*ELL**2*rho**2*(a*a+ay)*q1)
    chi = ELL*rho_b**3*q0
    hy = -hr/(ELL*rho)
    zeta = 3*(hy-a*h)/(2*p*p)+9*chi/(2*rho**4*p*p)
    # Separate local-flux adapter route must overlap the compact r formula.
    scalar_adapter = -p*zeta
    need(max(s.lo, scalar_adapter.lo) <= min(s.hi, scalar_adapter.hi),
         'SCALAR_LOCAL_FLUX_OVERLAP')
    need(vl.lo <= (hr+b*h).hi and (hr+b*h).lo <= vl.hi,
         'MIXED_VELOCITY_OVERLAP')
    zero = m.I(0)
    return enc({'background': {'rho': rho, 'p': p, 'a': a},
                'Q_r_jets': f,
                'auxiliary': {'b': b, 'd': d, 'a_y': ay, 'R': R, 'c': c,
                              'chi': chi, 'h_y': hy, 'zeta': zeta,
                              'scalar_from_local_flux': scalar_adapter},
                'profiles': {'h': h, 'h_r': hr, 'V_L': vl, 's': s, 's_r': sr},
                'U0': [h, zero, zero, h, s],
                'U0_r': [hr, zero, zero, hr, sr],
                'V0': [zero, vl, zero, zero, zero]})


def generate(repo, campaign_freeze_sha256):
    # Caller authenticates every local and source member before this import.
    source = repo/SOURCE_REL
    m = load_module(source/'numeric/interval_core.py', 'preparation_interval_core')
    cap = repo/CAP_REL
    old = json.loads((cap/'initial/references/RAW_BOXES.json').read_text())
    bump = json.loads((source/'numeric/inputs/BUMP_DERIVATIVE_BOUNDS.json').read_text())
    static = json.loads((source/'numeric/inputs/ACTUAL_ANNULUS_STATIC_JETS.json').read_text())
    shell = json.loads((source/'numeric/inputs/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json').read_text())['boxes']
    rho_b = m.I(*shell['rho_b']['exact'])
    py = m.I(*static['y_jets']['p'][1]['exact'])
    selected = [row for row in old['profile_cells'] if F(row['r'][1]) <= T]
    expected = [(F(1), F(65, 64), True)] + [
        (F(k, 1024), F(k+1, 1024), False) for k in range(1040, 1280)]
    need([(F(r['r'][0]), F(r['r'][1]), r['tail']) for r in selected] == expected,
         'EXACT_REGISTERED_CONTINUOUS_COVER')
    rows = []
    for oldrow in selected:
        row = {'r': oldrow['r'], 'tail': oldrow['tail']}
        row.update(initial_profile(m, oldrow['profiles'], oldrow['r'], oldrow['tail'],
                                   bump, rho_b, py))
        rows.append(row)
    oldtarget = next(x for x in old['targets'] if F(x['r']) == T)
    target = initial_profile(m, oldtarget['profiles'], [str(T), str(T)], False,
                             bump, rho_b, py)
    shared = {'rho_b': enc(rho_b), 'p_y_global': enc(py),
              'shell': {k: shell[k] for k in ('sigma_b', 'sigma_prime_b',
                         'sigma_second_b', 'U_prime_b')}}
    witness = {'schema': 'HDBLAST_P1_PREPARATION_RAW_WITNESS_V1',
               'campaign_freeze_sha256': campaign_freeze_sha256,
               'shared_inputs': shared, 'cells': rows, 'target': target,
               'calls': {'background': 0, 'evolved_actual_first_field': 0,
                         'bump_whole_cell': len(rows), 'bump_target': 1},
               'zero_extension': {'r_le_1': True, 'r_ge_2': True,
                                 'reason': 'Analytic compactness, exact moment derivative identity'}}
    return witness


def aggregate(witness, base, inject='none'):
    rows = witness['cells']
    target = witness['target']
    shell = witness['shared_inputs']['shell']
    rho_b = interval(witness['shared_inputs']['rho_b'])
    maxima = [[F(0)]*5 for _ in range(3)]
    jv = I(0)
    jt = I(0)
    av = F(0)
    at = F(0)
    for row in rows:
        lo, hi = map(F, row['r'])
        for caps, name in zip(maxima, ('U0', 'U0_r', 'V0')):
            for k, value in enumerate(row[name]):
                caps[k] = max(caps[k], interval(value).mag())
        scalar = interval(row['profiles']['s'])
        if inject == 'drop_conditioned_scalar':
            scalar = I(0)
        wv = T*(hi-lo)-(hi*hi-lo*lo)/2
        wt = hi-lo
        need(wv > 0 and wt > 0, 'POSITIVE_HISTORY_CELL_WEIGHTS')
        jv += wv*scalar
        jt += wt*scalar
        av += wv*scalar.mag()
        at += wt*scalar.mag()
    d0, dr, dv = map(gram_upper, maxima)
    free = max((d0+T*dv)/3, dr+dv)
    if inject == 'skip_initial_norm':
        free = F(0)
    error = Q_RESOLVENT*free/(1-Q_RESOLVENT)
    if inject == 'skip_resolvent_feedback':
        error = F(0)
    if inject == 'add_unchanged_first_forcing':
        error = (Q_RESOLVENT*free+T*F(17387513156573)+F(3851400))/(1-Q_RESOLVENT)
    m2 = interval(target['profiles']['s'])+I(-DELTA*error, DELTA*error)
    if inject == 'drop_conditioned_scalar':
        m2 = I(-DELTA*error, DELTA*error)
    if inject == 'skip_target_feedback':
        m2 = interval(target['profiles']['s'])
    ev = DELTA**3*error/6+ETA*av
    et = DELTA**2*error/2+ETA*at
    if inject == 'skip_kernel_comparison':
        ev = DELTA**3*error/6
        et = DELTA**2*error/2
    hv = jv+I(-ev, ev)
    ht = jt+I(-et, et)
    pick = lambda k: I(*shell[k]['exact'])
    sigma, sp, spp, vp = [pick(k) for k in ('sigma_b', 'sigma_prime_b',
                                           'sigma_second_b', 'U_prime_b')]
    fp = 2*vp-sp*spp/2+F(2, 3)*sigma*sp
    need(fp.hi < 0, 'STRICT_NEGATIVE_FIXED_SHELL_FPRIME')
    rhob = I(rho_b.lo, rho_b.hi)
    volume = (rhob*ELL)**2*fp*hv/2
    expansion = rhob*ELL*fp*ht/2
    if inject == 'drop_proper_time_factors':
        volume = fp*hv/2
        expansion = fp*ht/2
    postprocess = {'status': 'UNAVAILABLE_TRIGGER_NOT_MET'}
    if volume.hi < 0 and expansion.hi < 0:
        basev = interval(base['proper_v2'])
        baset = interval(base['proper_Theta2'])
        threshold = max((-basev.lo)/(-volume.hi), (-baset.lo)/(-expansion.hi))
        power = 0
        L = F(1)
        while L <= threshold:
            power += 1
            L *= 10
        lam = L if inject == 'reverse_lambda_sign' else -L
        newv = basev+lam*volume
        newt = baset+lam*expansion
        postprocess = {'status': 'PASS_FINITE_LAMBDA_SIMULTANEOUS_POSITIVE_SECOND_COEFFICIENTS'
                       if newv.lo > 0 and newt.lo > 0 else 'FAIL_FIXED_ALGEBRAIC_WITNESS',
                       'base_volume': basev, 'base_expansion': baset,
                       'threshold': threshold, 'smallest_power_of_ten_n': power,
                       'L': L, 'lambda': lam, 'new_volume_second_coefficient': newv,
                       'new_expansion_second_coefficient': newt,
                       'both_strict_positive': newv.lo > 0 and newt.lo > 0,
                       'interpretation': 'Formal second coefficients for one fixed finite lambda; no effective epsilon radius.'}
    excludes = lambda z: z.lo > 0 or z.hi < 0
    gates = {'M2_nonzero': excludes(m2), 'volume_history_nonzero': excludes(hv),
             'expansion_history_nonzero': excludes(ht), 'Fprime_negative': fp.hi < 0,
             'volume_slope_nonzero': excludes(volume),
             'expansion_slope_nonzero': excludes(expansion)}
    result = enc({'status': 'PASS_CONDITIONAL_NONZERO_PREPARATION_SLOPES'
                            if all(gates.values()) else 'INCONCLUSIVE_FROZEN_EXPERIMENT',
                  'campaign': 'ONE_DIRECTION_P1_N1024_T5_4', 'target': T, 'ell': ELL,
                  'oracle_calls': witness['calls'],
                  'initial_component_caps': dict(zip(('U0', 'U0_r', 'V0'), maxima)),
                  'initial_Gram_norm_caps': {'D0': d0, 'Dr': dr, 'Dv': dv},
                  'N_h': free, 'q': Q_RESOLVENT, 'E_h': error,
                  'free_scalar_target': interval(target['profiles']['s']),
                  'M2_slope': m2, 'J_volume': jv, 'A_volume': av,
                  'J_expansion': jt, 'A_expansion': at,
                  'volume_history_error': ev, 'expansion_history_error': et,
                  'kernel_volume_history': hv, 'kernel_expansion_history': ht,
                  'Fprime': fp, 'volume_slope': volume, 'expansion_slope': expansion,
                  'finite_lambda_witness': postprocess,
                  'gates': gates,
                  'premises': ['Authenticated shared frozen background boxes and operator/gauge/existence premises.',
                               'Fresh profile and C1 bounds for the registered direction; old prepared profiles are not reused as profile bounds.',
                               'Full coupled homogeneous five-mean resolvent; unchanged first forcing cancels only in the difference.',
                               'No new independent PDE observation, nonlinear epsilon radius, physical radiation or cosmological prediction.']})
    return result
