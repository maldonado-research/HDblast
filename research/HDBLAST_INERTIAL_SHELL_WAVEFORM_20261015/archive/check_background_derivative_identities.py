"""Formal coefficient identities only; not an analytic/PDE proof certificate."""
import json
from pathlib import Path
import sympy as s

r, a, p, h, u1, u2, u3, w2, mu, t = s.symbols('rho a p h U1 U2 U3 Wsecond mu t', nonzero=True)
Q = 2*h*h-10*a*h+s.Rational(35, 4)*a*a-u2-p*p/6
flow = {r: -r*r*a, a: 1/r+r*p*p/3, p: -r*p*(h-4*a),
        h: r*(h*h-4*a*h-u2), u2: -r*p*u3}
direct = s.expand(sum(s.diff(r*r*Q, x)*v for x, v in flow.items()))
factored_Qz = ((4*h-10*a)*flow[h]+(-10*h+s.Rational(35, 2)*a)*flow[a]
                +r*p*u3-p*flow[p]/3)
factored = -2*r**3*a*Q+r*r*factored_Qz
expanded = (r**3*(4*h**3-30*a*h*h+60*a*a*h-s.Rational(35, 2)*a**3
                   +(12*a-4*h)*u2+(s.Rational(29, 6)*a-3*h)*p*p+p*u3)
            +r*(s.Rational(35, 2)*a-10*h))
gain_ratio_derivative = ((-r*p*u2)*p-u1*(-r*(u1-4*a*p)))/(p*p)
gamma = r*(h-s.Rational(5, 2)*a)
inverse_mu = r*(h-4*a+w2)
gap = r*(w2-s.Rational(3, 2)*a)
checks = {
    'W_z_factored_equals_total_derivative': s.cancel(direct-factored) == 0,
    'W_z_expanded_independent_referee_formula': s.cancel(direct-expanded) == 0,
    'h_ratio_derivative': s.cancel(gain_ratio_derivative.subs(u1, h*p)-flow[h]) == 0,
    'inverse_mu_minus_gamma_cancellation': s.cancel(inverse_mu-gamma-gap) == 0,
    'inverse_mu_minus_beta0_cancellation': s.cancel(inverse_mu-(gamma+5*mu/2)-(gap-5*mu/2)) == 0,
    'W_t_exact': s.simplify(s.diff(-s.Rational(35, 4)*s.tanh(t)**2, t)
                          +s.Rational(35, 2)*s.tanh(t)*s.sech(t)**2) == 0,
}
controls = {
    'omit_U_third_term_rejected': s.cancel(direct-(expanded-r**3*p*u3)) != 0,
    'flip_inward_derivative_sign_rejected': s.cancel(direct+expanded) != 0,
    'wrong_shell_half_factor_rejected': s.cancel(inverse_mu-gamma-r*(w2-3*a)) != 0,
    'omit_beta0_inertia_shift_rejected': s.cancel(inverse_mu-(gamma+5*mu/2)-gap) != 0,
}
if not all(checks.values()) or not all(controls.values()):
    raise ValueError('Formal derivative identity check failed')
receipt = {'status': 'PASS_FORMAL_COEFFICIENT_IDENTITIES_ONLY', 'sympy_version': s.__version__,
           'checks': checks, 'wrong_formula_controls': controls,
           'scope': 'Formal identities and deliberately wrong coefficients; no analytic branch/PDE proof certification.'}
Path(__file__).with_name('BACKGROUND_DERIVATIVE_IDENTITY_CHECK.json').write_text(json.dumps(receipt, indent=2)+'\n')
print(json.dumps(receipt, indent=2))
