"""Exact envelope for the actual canonical first seed (pinned producer)."""
import importlib.util
from fractions import Fraction as F
from math import factorial, prod
from pathlib import Path

def evaluate(claims, check, producer_path=None):
    if producer_path is None:
        producer_path = Path(__file__).resolve().parents[1] / 'inputs/PINNED_SOURCE_JETS.py'
    spec = importlib.util.spec_from_file_location('pinned_seed_majorants', producer_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    check('seed-pinned-scale-and-orders', module.N == 16 and module.ELL == F(1, 10 ** 6) and (module.KAPPA == F(1, 10 ** 8)))
    check('seed-pinned-potential', module.VCOEFF == [F(-2, 27), 0, F(14, 9), F(50, 27), F(-1, 6), F(-4, 9), F(-2, 27)])
    bg = module.background()
    bump, _ = module.bump_bounds()
    q_bump, _ = module.compose_bump(bump, bg['b'])
    source = module.source_jets(bg, q_bump)
    profiles = {'u': (source['u'], 14, 86), 'v': (source['v'], 14, 97), 'tau_profile': (F(3, 2) * source['w_velocity'] * bg['b'], 13, 77), 'omega_profile': (source['vt'] * bg['b'], 13, 88), 'a': (bg['a'], 14, 0), 'b': (bg['b'], 14, 0)}
    for key, (jet, order, exponent) in profiles.items():
        ceiling = F(10 ** exponent)
        for n in range(order + 1):
            check(f'seed-radial-{key}-{n}', factorial(n) * jet.c[n] / module.ELL ** n <= ceiling)
    term_factor = 4 * 6 * 2 ** 15 * prod((4 * (5 + 4 * j) for j in range(14)))
    check('seed-all-slot-term-envelope', term_factor * 10 ** 97 < 10 ** int(claims['raw_tensor_exponent']))
    x = F(4, 3)
    n = 16
    exp_upper = sum((x ** k / factorial(k) for k in range(n + 1)), F(0))
    exp_upper += x ** (n + 1) / factorial(n + 1) / (1 - x / F(n + 2))
    check('seed-inverse-weight-bound', exp_upper < 4)
    check('seed-q-lower', F(1, 10 ** 16) * F(3, 4) ** 2 > F(1, 10 ** 17))
    check('seed-volume-square-bound', 2 * 4 ** 2 * 131 ** 4 * F(1, 10 ** 6) < 100 ** 2)
    norm_bound = (2 * 15 * 4 * 10 ** 34 + 2 * 14 * 4 * 10 ** 17) * 100 * 10 ** int(claims['raw_tensor_exponent'])
    check('seed-four-field-sum-bound', norm_bound < 10 ** int(claims['seed_norm_exponent']))
    check('seed-common-envelope', int(claims['seed_norm_exponent']) < int(claims['B_seed_exponent']))
    return {'radial_profile_count': len(profiles), 'radial_derivative_count': sum((o + 1 for _, o, _ in profiles.values())), 'seed_norm_exponent': claims['seed_norm_exponent'], 'B_seed_exponent': claims['B_seed_exponent']}
