"""Exact tensor factorization and registered potential-sign checks.

This checks the algebra used by the companion proof; functional analysis and
the physical transverse-traceless reduction require their stated arguments.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sympy as s


def run():
    y = s.symbols('y', positive=True)
    rho = s.Function('rho')(y)
    h = s.Function('h')(y)
    a, p, U, ap = s.symbols('a p U ap', real=True)
    mu = s.symbols('mu', real=True)
    checks = {}

    def zero(name, expression):
        reduced = s.simplify(expression)
        checks[name] = reduced == 0
        if reduced != 0:
            raise ValueError(f'{name}: residual {reduced}')

    dz = lambda f: rho * s.diff(f, y)
    psi = rho**s.Rational(3, 2) * h
    S = s.Rational(3, 2) * s.diff(rho, y)
    Vplus = S**2 + dz(S)
    Vminus = S**2 - dz(S)
    D = lambda f: dz(f) - S * f
    Dstar = lambda f: -dz(f) - S * f
    zero('TT_equation_to_Schrodinger',
         (-dz(dz(psi)) + Vplus*psi - mu*psi)
         + rho**s.Rational(7, 2)*(s.diff(h,y,2)
              + 4*s.diff(rho,y)*s.diff(h,y)/rho + mu*h/rho**2))
    zero('DstarD_factorization', Dstar(D(psi)) + dz(dz(psi)) - Vplus*psi)
    zero('DDstar_factorization', D(Dstar(psi)) + dz(dz(psi)) - Vminus*psi)
    zero('Neumann_to_D_boundary', D(psi)-rho**s.Rational(5,2)*s.diff(h,y))
    zero('massless_graviton_kernel', D(rho**s.Rational(3,2)))
    zero('partner_formal_kernel', Dstar(rho**s.Rational(-3,2)))

    Vmin_y = Vminus.subs({s.diff(rho,y,2):rho*(a*a+ap), s.diff(rho,y):rho*a})
    Vmin_bg = s.expand(Vmin_y.subs(ap,-1/rho**2-p*p/3))
    Vmin_constraint = s.expand(Vmin_bg.subs(a*a,1/rho**2+p*p/12-U/6))
    expected = s.Rational(9,4)+rho**2*(s.Rational(9,16)*p*p-U/8)
    zero('partner_potential_on_background', Vmin_constraint-expected)
    eta=s.symbols('eta', real=True)
    W=s.Rational(1,3)+eta**2+eta**3/3
    potential=s.Poly(s.expand(s.diff(W,eta)**2/2-s.Rational(2,3)*W**2),eta)
    r=s.Rational(1,100)
    upper=potential.nth(0)+sum(abs(potential.nth(i))*r**i
                               for i in range(1,potential.degree()+1))
    if not upper < 0:
        raise ValueError('Potential not certified negative on polynomial radius')
    checks['registered_U_negative_entire_radius']=True

    wrong_controls={
        'wrong_scalar_gradient_coefficient':
            s.simplify(expected-(s.Rational(9,4)+rho**2*(s.Rational(7,16)*p*p-U/8)))!=0,
        'wrong_potential_sign':
            s.simplify(expected-(s.Rational(9,4)+rho**2*(s.Rational(9,16)*p*p+U/8)))!=0,
        'wrong_factor_order':s.simplify(Vplus-Vminus)!=0,
        'wrong_even_boundary':s.simplify(Dstar(psi)-rho**s.Rational(5,2)*s.diff(h,y))!=0,
    }
    if not all(wrong_controls.values()):
        raise ValueError('A deliberately wrong formula was not rejected')
    return {'status':'PASS_TENSOR_ALGEBRA_AND_EXACT_POTENTIAL_SIGN',
            'checks':checks,'check_count':len(checks),
            'wrong_formula_controls':wrong_controls,
            'potential_polynomial':str(potential.as_expr()),
            'radius':str(r),'U_upper_bound':str(upper),
            'U_upper_bound_decimal':str(s.N(upper,20)),
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'scope':'Algebra and exact polynomial sign only; operator domains, spectral mapping and physical mode classification are proved/reviewed separately.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',required=True)
    args=parser.parse_args();result=run()
    Path(args.output).write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(result['status'],result['check_count'],'identities/sign checks;',
          len(result['wrong_formula_controls']),'wrong formulas rejected')
