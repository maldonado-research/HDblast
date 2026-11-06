#!/usr/bin/env python3
"""Exact scalar-square and weighted-curvature-source bounds from retained inputs.

No PDE producer is imported or replayed. Analytic premises are in SCALAR_SOURCE_BOUNDS.md.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from math import isqrt
from pathlib import Path


def need(condition, message):
    if not condition:
        raise ValueError(message)


PINS = {
    'HIGHER_JETS_NORMAL.json': '92b17fd71d78f83eb3f2c5ec8b9ea5f4e70365ec4d41aa4c0205551e7eaeb5df',
    'INERTIAL_WAVEFORM_ENCLOSURE.json': 'd2ffbcec5413580ff5703da58b062a9b83024431748fa40726555c09c2a968d8',
    'PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json': '8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8',
    'PHYSICAL_SHELL_RESPONSE.md': 'a78f2239893fdd0e5ff653d4d5ce5cd9d3244c57c794ed2d3bb43bd4f5e1484b',
    'REGISTERED_EXACT_FAMILY_BINDING.md': '69bfe846cc4fda8ec6a981969fbdd12fff7862b5ca7cf31a19d696b0a839165a',
    'HIGHER_JETS.md': '4c6453049f0d91e9223f06b8c881cd6bd42164d84473bf5fa5a5357bfbe35c8a',
    'PAST_FUTURE_PROOF.md': 'ddfe9ad88773133f8f79758802724058453fb3b0c75e13aafb892ea960c05c10',
    'QUADRATIC_AREA_BOUND.md': '56cfb7b4d2817df3c874a8a204cc3cf5931119b9e32fe06c822f38d1c554ceb4',
}


WITNESS_Q_NUMERATOR = 5


def sqrt_upper(a, bits=96):
    """Upper dyadic enclosure, certified by squaring, with no floating point."""
    need(a >= 0, 'nonnegative square root argument')
    k = isqrt((a.numerator << (2*bits)) // a.denominator)
    r = F(k, 1 << bits)
    if r*r < a:
        r += F(1, 1 << bits)
    need(r*r >= a and (r-F(1, 1 << bits))**2 < a if a else r == 0,
         'outward dyadic square root')
    return r


def trim(a):
    a = list(map(F, a))
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def add(a, b):
    return trim([(a[i] if i < len(a) else 0)+(b[i] if i < len(b) else 0)
                 for i in range(max(len(a), len(b)))])


def scale(a, c):
    return trim([c*x for x in a])


def mul(a, b):
    ans = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            ans[i+j] += x*y
    return trim(ans)


def derivative(a):
    return trim([i*a[i] for i in range(1, len(a))] or [0])


def value(a, x):
    result = F(0)
    for c in reversed(a):
        result = result*x+c
    return result


def polynomial_checks(delta, c):
    """Differentiate the full V,sigma model rather than the reduced P'' formula."""
    W = [1, -1, 0, F(1, 3)]
    Wp = derivative(W)
    V = add(scale(mul(Wp, Wp), F(1, 2)), scale(mul(W, W), F(-2, 3)))
    sigma = add(scale(W, 2), [delta, delta*c])
    sp = derivative(sigma)
    P = add(add(scale(mul(sp, sp), F(1, 16)), scale(V, F(-1, 2))),
            scale(mul(sigma, sigma), F(-1, 12)))
    reduced = add(scale(add(scale(Wp, c/F(4)),
                              scale(mul([1, c], W), F(-1, 3))), delta),
                  scale(add([c*c/F(16)],
                            scale(mul([1, c], [1, c]), F(-1, 12))), delta*delta))
    expected = [delta*F(7, 6)*c-delta*delta*c*c/6,
                -delta*F(2, 3), -delta*F(4, 3)*c]
    need(P == reduced, 'full polynomial P cancellation')
    Ppp = derivative(derivative(P))
    need(Ppp == expected, 'independent full-model Psecond differentiation')
    return Ppp


def source_joins(jets, shell, waveform):
    need(jets['status'] == 'PASS', 'higher jets status')
    need(shell['status'] == 'PASS_EXACT_RATIONAL_CONSEQUENCES_OF_RETAINED_ANALYTIC_THEOREM',
         'conditional actual shell status')
    need(shell['parameters_exact'] == {'delta': '1/1000', 'c': '96714321425692/161839258930731'},
         'same actual model parameters')
    for key, name in [('inputs/prior/INERTIAL_WAVEFORM_ENCLOSURE.json', 'INERTIAL_WAVEFORM_ENCLOSURE.json'),
                      ('inputs/prior/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json', 'PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json'),
                      ('inputs/reference/PHYSICAL_SHELL_RESPONSE.md', 'PHYSICAL_SHELL_RESPONSE.md')]:
        need(jets['input_pins'][key] == PINS[name], 'higher-report retained join '+name)
    ell = F(jets['source_bound_adapter']['ell'])
    need(ell == F(1, 10**6) == F(jets['exact']['interior']['ell']) == F(waveform['state']['ell']),
         'same registered ell and scaled derivative units')
    need(waveform['state']['delta'] == '1/1000'
         and waveform['state']['initial_velocity'] == 'exp(-z) Q(z/ell)',
         'same nonzero actual initial velocity')
    b = jets['exact']['boundary']
    mu = tuple(map(F, b['mu']))
    need(mu == tuple(map(F, waveform['shell_parameter_intervals']['mu'])), 'same unscaled mu')
    need(tuple(map(F, b['gamma'])) == tuple(map(F, waveform['shell_parameter_intervals']['gamma'])),
         'same unscaled gamma')
    actual_mu = tuple(map(F, shell['boxes']['mu']['exact']))
    need(0 < mu[0] <= actual_mu[0] <= actual_mu[1] <= mu[1] < F(1, 400),
         'actual shell mu in retained report')
    need(F(b['ell_over_mu']) == ell/mu[0], 'ell/mu is not Lambda')
    need(mu[0]/ell > 1000, 'positive dynamic endpoint mass')
    rho = tuple(map(F, shell['boxes']['rho_b']['exact']))
    p = tuple(map(F, shell['boxes']['p_b']['exact']))
    need(129 <= rho[0] <= rho[1] <= 131, 'actual rho box')
    need(F(-7, 50000) <= p[0] <= p[1] <= F(-3, 25000), 'actual negative p away from zero')
    phi = tuple(1+F(x) for x in shell['boxes']['eta_b']['exact'])
    need(phi == tuple(F(x)/4 for x in shell['boxes']['sigma_second_b']['exact']),
         'actual phi agrees with sigma_second=4phi')
    need(0 < phi[0] <= phi[1] < 1, 'actual near-plus-one branch, no phi=1 substitution')
    return ell, mu, rho, p, phi


def graph_bounds(ell, mu_lo, h, velocity_difference=2):
    need(velocity_difference == 2, 'actual past-minus-future velocity factor two')
    # Norm of g=ell exp(-ell r)Q on its length-one support.
    # Complete basic, half-graph and graph-1 energies retain endpoint mass/spring.
    eps = F(35, 4)*ell**4*h*h/(5500*ell*ell)
    need(0 <= eps < 1, 'A0+omega positive on shortened interval')
    G0 = velocity_difference*ell
    Ghalf = sqrt_upper(4*ell**2*((16+ell)**2+103000*ell**2)/(1-eps))
    G1 = velocity_difference*ell*(256+32*ell+103001*ell**2)
    need(Ghalf < 33*ell and G1 < 513*ell, 'complete graph energy norm ceilings')
    endpoint = sqrt_upper(ell/mu_lo)
    dy1 = endpoint*G0
    dy0 = (h-1)*dy1
    Omega = F(35, 4)*ell**4*h*h
    Omega1 = F(35, 2)*ell**4*h
    dy2 = endpoint*Ghalf+Omega*dy0
    dy3 = endpoint*G1+Omega*dy1+Omega1*dy0
    return dict(epsilon=eps, endpoint_trace_factor=endpoint, G0=G0, Ghalf=Ghalf, G1=G1,
                omega_abs=Omega, omega_s_abs=Omega1,
                delta_y=dy0, delta_yprime=dy1, delta_ysecond=dy2, delta_ythird=dy3)


def algebra_checks():
    # Independent fixtures for the cancellation and its differentiated form.
    # Variables are algebraic jet values of an actual endpoint solution; all
    # omega_s terms must cancel from the differentiated reconstruction.
    fixtures = [(F(2, 7), F(-3, 11), F(5, 13), F(-7, 17), F(11, 19), F(-13, 23)),
                (F(-3, 5), F(2, 9), F(-4, 7), F(6, 13), F(-8, 19), F(10, 29))]
    for y, yp, Ab, Avb, om, omp in fixtures:
        e, alpha, alphap = F(1, 10**6), F(5, 2), F(0)
        ypp, yppp = -Ab-om*y, -Avb-om*yp-omp*y
        a, ap = alpha+om/e**2, omp/e**3
        direct = ypp/e**2+a*y
        canceled = -Ab/e**2+alpha*y
        direct_p = yppp/e**3+a*yp/e+ap*y
        canceled_p = -Avb/e**3+alpha*yp/e
        need(direct == canceled and direct_p == canceled_p,
             'complete endpoint reconstruction and omega_s cancellation')
    # This exact spectral datum makes the initial velocity factor observable.
    # For z''+a z=0, data (0,-2g), z_s(0)=-2g and A z_s(0)=-2ag.
    a, g = F(5), F(3, 7)
    need(abs(a*(-2*g)) == 2*abs(a*g), 'nonzero velocity graph datum fixture')
    return 5


def weighted_integral_bound(waveform, ell, T, deltaS, Smax, Cmin, Khi):
    """Cauchy bound from one actual scalar witness and its quiet interval."""
    witness = waveform['sign_witnesses']['rising']
    need(F(witness['x']) == F(WITNESS_Q_NUMERATOR, 4), 'same actual interior scalar witness time')
    actual_upper = F(witness['ell_S'][1])/ell+deltaS
    need(actual_upper < -37376, 'actual scalar witness is divided by ell and transported')
    a = F(WITNESS_Q_NUMERATOR, 4)*ell
    need(0 < ell < a < T and T*T/2 < 1, 'Cauchy subinterval and cosh comparison')
    # cosh(T)<=sum_(k>=0)(T^2/2)^k=(1-T^2/2)^(-1).
    H = (1-T*T/2)**-2
    positive = F(27, 40)*(T-a)/(a-ell)*F(37376)**2/H
    need(F(27, 40)*(T-a)/(a-ell) == F(27, 32), 'weighted Cauchy exact coefficient')
    negative = H*(-Cmin)*Smax*Smax*(T-ell)**2/2
    lower = positive-negative
    upper = H*Khi*(T-ell)**2/2
    need(lower > 1170000000 and upper < 1040000000000,
         'strict nonzero known scalar curvature integral')
    return dict(witness_q=a, actual_past_witness_upper=actual_upper,
                cosh_squared_upper=H, derivative_square_integral_lower=positive,
                negative_potential_integral_abs_upper=negative,
                integral_W_K_R_lower=lower, integral_W_K_R_upper=upper,
                readable_integral_strict_box=[F(1170000000), F(1040000000000)],
                scalar_Raychaudhuri_particular_upper=-lower)


def evaluate(root):
    root = Path(root)
    pinfile = root/'SOURCE_PINS.json'
    need(pinfile.is_file() and not pinfile.is_symlink(), 'regular source metadata')
    pinmeta = json.loads(pinfile.read_text())
    base = '5ed540606f7ad5096a8baa9749cb13d753a551ef'
    need(pinmeta['ref'] == base and pinmeta['repository'] == 'https://github.com/maldonado-research/HDblast',
         'immutable public source base')
    rows = pinmeta['files']
    need(len(rows) == len(PINS) and {r['local'] for r in rows} == {'inputs/'+n for n in PINS},
         'complete public source metadata')
    for name, digest in PINS.items():
        path = root/'inputs'/name
        need(path.is_file() and not path.is_symlink(), 'regular source input '+name)
        data = path.read_bytes()
        need(hashlib.sha256(data).hexdigest() == digest, 'source bytes '+name)
        row = next(r for r in rows if r['local'] == 'inputs/'+name)
        need(row['sha256'] == digest and row['bytes'] == len(data) and row['ref'] == base,
             'source metadata digest/length/ref '+name)
        blob = hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
        need(row['git_blob_sha1'] == blob, 'source metadata Git blob '+name)
        need(row['url'] == 'https://github.com/maldonado-research/HDblast/blob/'+base+'/'+row['path']
             and row['raw_url'] == 'https://raw.githubusercontent.com/maldonado-research/HDblast/'+base+'/'+row['path'],
             'immutable source URLs '+name)
    jets = json.loads((root/'inputs/HIGHER_JETS_NORMAL.json').read_text())
    shell = json.loads((root/'inputs/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json').read_text())
    waveform = json.loads((root/'inputs/INERTIAL_WAVEFORM_ENCLOSURE.json').read_text())
    ell, mu, rho, p, phi = source_joins(jets, shell, waveform)
    h, T = F(25, 16), F(25, 16)*ell
    b = jets['exact']['boundary']
    future = dict(y=F(b['y']), yprime=F(b['yprime']), ysecond=F(b['ysecond']),
                  ythird=F(jets['exact']['interior']['ythird']))
    need(future['ysecond'] == F(jets['exact']['interior']['ysecond']), 'same actual boundary second jet')
    need(jets['exact']['interior']['scaled_domain']['x_max'] == '3' and h < 3,
         'short interval contained in actual future domain')
    graph = graph_bounds(ell, mu[0], h)
    past = {key: val+graph['delta_'+key] for key, val in future.items()}
    gain = 3*mu[1]/(rho[0]**3*(-p[1]))
    afuture = future['ysecond']/ell**2+F(5, 2)*future['y']
    Sfuture = gain*afuture
    Sqfuture = gain*(future['ythird']/ell**3+F(5, 2)*future['yprime']/ell
                    +F(35, 2)*T*future['y']+F(3, 2)*T*afuture)
    delta_inner = graph['endpoint_trace_factor']*graph['Ghalf']/ell**2+F(5, 2)*graph['delta_y']
    deltaS = gain*delta_inner
    deltaSq = gain*(graph['endpoint_trace_factor']*graph['G1']/ell**3
                   +F(5, 2)*graph['delta_yprime']/ell+F(3, 2)*T*delta_inner)
    Smax, Sqmax = Sfuture+deltaS, Sqfuture+deltaSq
    # A separate uncanceled, transferred-y''' route must enclose the result.
    apast = past['ysecond']/ell**2+F(5, 2)*past['y']
    directS = gain*apast
    directSq = gain*(past['ythird']/ell**3+F(5, 2)*past['yprime']/ell
                    +F(35, 2)*T*past['y']+F(3, 2)*T*apast)
    need(Smax <= directS and Sqmax <= directSq, 'independent direct past-jet envelope dominates cancellation')
    delta, c = F(shell['parameters_exact']['delta']), F(shell['parameters_exact']['c'])
    Ppp = polynomial_checks(delta, c)
    need(value(derivative(Ppp), phi[0]) < 0 and value(derivative(Ppp), phi[1]) < 0,
         'Psecond decreasing on positive actual phi interval')
    Plo, Phi = value(Ppp, phi[1]), value(Ppp, phi[0])
    need(Plo <= Phi < 0, 'actual negative Psecond enclosure')
    # The physical background response gives Pprime=-p/(2rho*mu).
    scalar2_coefficient = [-rho[0]*p[1]/(2*mu[1]), -rho[1]*p[0]/(2*mu[0])]
    need(3 < scalar2_coefficient[0] <= scalar2_coefficient[1] < 5,
         'actual nonzero second-scalar coefficient between three and five')
    # sech^2 q is between 1-T^2 and 1 since |tanh q|<=q.
    Cmin = -F(1, 5)+F(9, 20)*rho[1]**2*Plo
    Cmax = -F(1, 5)*(1-T*T)+F(9, 20)*rho[0]**2*Phi
    need(Cmin <= Cmax < 0, 'entire S-squared coefficient is negative')
    Klo, Khi = Cmin*Smax**2, F(27, 40)*Sqmax**2
    need(-Klo < Khi, 'upper source magnitude dominates negative lower magnitude')
    readable = dict(S_abs=F(51500), Sq_abs=F(3107000000000),
                    Psecond_lower=F(-766143, 10**9), Psecond_upper=F(-766128, 10**9),
                    K_lower=F(-16000000000), K_upper=F(6520000000000000000000000))
    need(Smax < readable['S_abs'] and Sqmax < readable['Sq_abs'], 'strict readable first-scalar bounds')
    need(readable['Psecond_lower'] < Plo <= Phi < readable['Psecond_upper'], 'readable actual Psecond bounds')
    need(readable['K_lower'] < Klo <= Khi < readable['K_upper'], 'strict readable weighted source enclosure')
    integral = weighted_integral_bound(waveform, ell, T, deltaS, Smax, Cmin, Khi)
    # This is a necessary condition, not an evaluated unknown-field bound.
    # The full-family identity and quiet cylinder are separate analytic premises.
    KD_upper = F(5357, 10**19)
    CR_upper = F(8614, 10**16)
    quiet = ell/4
    CR = integral['cosh_squared_upper']*(T-quiet)**2/2
    need(CR < CR_upper, 'full unknown source uses ell/4 quiet interval')
    threshold = (integral['integral_W_K_R_lower']-KD_upper)/CR_upper
    need(threshold > 1350000000000000000000, 'necessary compensation exceeds 1.35e21')
    compensation = dict(condition='IF the full family and Gaussian cylinder exist through T and D2/D0>=0',
                        unknown='Z=rho_b^2 Pprime_b mean_w(psi2)-mean_w(E2)',
                        full_quiet_q=quiet, K_D_strict_upper=KD_upper,
                        unknown_Green_multiplier=CR, unknown_Green_multiplier_strict_upper=CR_upper,
                        Z_sup_abs_necessary_strict_lower=threshold,
                        readable_Z_sup_abs_necessary_strict_lower=F(1350000000000000000000),
                        unknown_projection_evaluated=False)
    exact = dict(ell=ell, q_max=T, scaled_q_max=h, rho_box=rho, p_box=p, mu_box=mu,
                 phi_box=phi, physical_gain_upper=gain, future_scaled_boundary=future,
                 complete_graph_transfer=graph, actual_past_scaled_boundary=past,
                 S_future_abs=Sfuture, Sq_future_abs=Sqfuture,
                 S_past_future_difference_abs=deltaS, Sq_past_future_difference_abs=deltaSq,
                 S_actual_abs=Smax, Sq_actual_abs=Sqmax,
                 direct_past_S_abs=directS, direct_past_Sq_abs=directSq,
                 Psecond_actual_box=[Plo, Phi], rho_squared_Pprime_actual_box=scalar2_coefficient,
                 S_squared_coefficient_box=[Cmin, Cmax],
                 K_R_known_box=[Klo, Khi], readable_strict_bounds=readable,
                 weighted_Green_integral=integral, conditional_compensation=compensation)
    checks = ['eight public input byte identities and immutable Git source metadata', 'actual state, branch, shell and unscaled-mu joins',
              'complete positive endpoint basic/half-graph/graph-one transfer on 25ell/16',
              'actual future and transferred past boundary C3 jets',
              'physical reconstruction and differentiated endpoint cancellation',
              'full polynomial model P and Psecond identities', 'actual phi monotonic Psecond enclosure',
              'actual nonzero rho_squared_Pprime coefficient enclosure',
              'continuous weighted first-scalar source enclosure',
              'actual past witness transport and nonzero weighted Cauchy integral',
              'conditional compensation threshold with full ell/4 quiet interval']
    algebra_checks()
    result = dict(schema='hdblast-projected-curvature-source-v1', status='PASS',
                named_checks=checks, named_check_count=len(checks), exact=exact,
                display_only={k: float(exact[k]) for k in ['S_future_abs', 'Sq_future_abs',
                    'S_past_future_difference_abs', 'Sq_past_future_difference_abs', 'S_actual_abs', 'Sq_actual_abs']}
                    | {'Psecond_actual_box': [float(Plo), float(Phi)],
                       'rho_squared_Pprime_actual_box': [float(v) for v in scalar2_coefficient],
                       'S_squared_coefficient_box': [float(Cmin), float(Cmax)],
                       'K_R_known_box': [float(Klo), float(Khi)],
                       'weighted_Green_integral': {k: float(v) for k, v in integral.items()
                                                   if isinstance(v, F)},
                       'conditional_Z_sup_abs_necessary_lower': float(threshold)},
                input_pins=PINS,
                scope=dict(conditional_on_retained_actual_smooth_development=True,
                           all_bounds_continuous_on_q_interval=True,
                           actual_nonzero_velocity_retained=True,
                           nonlinear_target_cylinder_existence_assumed_not_proved=True,
                           full_R2_or_Cr_evaluated=False, psi2_or_electric_Weyl2_evaluated=False,
                           uniform_amplitude_radius_evaluated=False, old_PDE_campaign_replayed=False))
    return json.loads(json.dumps(result, default=serial))


def serial(obj):
    if isinstance(obj, F):
        return str(obj)
    raise TypeError(type(obj).__name__)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = evaluate(args.root)
    rendered = json.dumps(result, default=serial, indent=2, sort_keys=True)+'\n'
    if args.output:
        args.output.write_text(rendered)
    else:
        print(rendered, end='')


if __name__ == '__main__':
    main()
