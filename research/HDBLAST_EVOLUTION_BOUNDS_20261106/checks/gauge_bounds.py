"""Exact proper-y radial majorants and explicit unstraightened norm bridge."""
from fractions import Fraction as F
from math import factorial
import importlib.util


def evaluate(root, data, check):
    def need(condition, name):
        check(name, condition)
    engine = root/'inputs/derive_source_jets_06e.py'
    spec = importlib.util.spec_from_file_location('pinned_first_source_jets', engine)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    need(m.N >= 16, 'sixteen source Taylor orders retained')
    need(m.ELL == F(data['ell']), 'retained proper-y source normalization')
    j = m.background()
    bump, _ = m.bump_bounds()
    Q, _ = m.compose_bump(bump, j['b'])
    src = m.source_jets(j, Q)
    rho, a, u = j['rho'], j['a'], src['u']
    radial_S = src['w_velocity']*j['b']

    def bound(jet, order):
        return max(factorial(n)*jet.c[n]/m.ELL**n for n in range(order+1))

    need(bound(u, 15) < 10**data['u_derivative_decade'], 'proper-y u jets through fifteen')
    need(bound(radial_S, 15) < 10**data['S_derivative_decade'], 'proper-y S profile jets through fifteen')
    multipliers = {'rho': rho, 'rho2': rho**2, 'rho3': rho**3,
                   'rho2a': rho**2*a, 'a': a}
    profiles = {'rho2u2': rho**2*u*u,
                'rhou2': rho*u*u,
                'rho3uS': rho**3*u*radial_S,
                'rhouS': rho*u*radial_S,
                'rho2uuy': rho**2*u*u.dy(),
                'rho2au2': rho**2*a*u*u}
    ceilings = {}
    for name, jet in multipliers.items():
        value = bound(jet, 14)
        need(value < 10**data['fixed_radial_decade'], f'radial fixed multiplier {name} through fourteen')
        ceilings[name] = m.decade_ceiling(value)
    for name, jet in profiles.items():
        value = bound(jet, 14)
        need(value < 10**data['known_quadratic_radial_decade'], f'radial known quadratic {name} through fourteen')
        ceilings[name] = m.decade_ceiling(value)
    for stage in range(14):
        branches = stage+1+1+3*(2+stage)
        need(branches <= 60 and branches < 1000,
             f'all derivative-slot branches stage {stage}')
    bfix_exact = 1000**14*10**(14*data['fixed_radial_decade']+data['fixed_radial_decade']+200)
    need(bfix_exact < 10**data['fixed_multiplier_C14_decade'], 'full geometric fixed multiplier C14 budget')
    volume_upper = m.F(40*131**4, 10**6)
    need(volume_upper < 10**6, 'full geometric annulus volume budget')
    need(2**28*15 < 10**12, 'covariant Leibniz factor through fourteen')
    bknown = 4*10**6*1000**14*10**(14*data['fixed_radial_decade']+data['known_quadratic_radial_decade']+200)
    need(bknown < 10**data['known_quadratic_geometric_decade'], 'full geometric known quadratic budget')
    need(1+8*10**6*10**data['fixed_multiplier_C14_decade'] < 10**data['linear_completion_decade'], 'linear completion coefficient budget')
    need((1+8*10**6*10**data['fixed_multiplier_C14_decade'])*10**data['physical_I2_input_decade']+bknown < 10**data['unstraightened_wave_Cauchy_per_copy_decade'],
         'unstraightened wave Cauchy second norm below 10^3641')
    need(2*10**data['unstraightened_wave_Cauchy_per_copy_decade'] < 10**data['rounded_one_or_two_copy_decade'], 'two-copy rounded wave Cauchy bound')
    return {'radial_derivative_decade_ceilings_through_14': ceilings,
            'fixed_multiplier_C14_bound': '10^'+str(data['fixed_multiplier_C14_decade']),
            'known_quadratic_geometric_bound': '10^'+str(data['known_quadratic_geometric_decade']),
            'physical_I2_input_bound': '10^'+str(data['physical_I2_input_decade']),
            'unstraightened_wave_Cauchy_bound_per_copy': '10^'+str(data['unstraightened_wave_Cauchy_per_copy_decade']),
            'rounded_bound_one_or_two_copies': '10^'+str(data['rounded_one_or_two_copy_decade']),
            'annulus_volume_rational_upper': str(volume_upper)}
