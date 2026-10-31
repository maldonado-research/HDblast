#!/usr/bin/env python3
"""Independent exact screen-map reduction and unit audit; no evolved fields."""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
SOURCE_PINS = {
 'research/HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026/theory/DERIVATION.md': 'f27efc433a20689e97ea45e2c0cfd0918323318a4a1cdc4e798166de56e0172b',
 'research/HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026/theory/screen_coefficient.py': '359f4771a27b490614c91a44266d156d1696abc35f72795b68564703e9fe49bd',
 'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/actual_first_jets.py': 'a68fa8025d2170efca6c3e85d443910a796d0c30a9b7dd17ed988cd9e00bbb00',
 'research/HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026/review/INDEPENDENT_SCREEN_AUDIT.md': '068d40a4fe8664053962749f02846805fc52c753fc2af9edcea2b22a29df87a3',
}
EPSILON_TERM_ENABLED = True
COMPACT_B_SIGN = 1
LINEAR_PARAMETER_FACTOR = Q(1)


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def authenticate(root):
    for rel, pin in SOURCE_PINS.items():
        p = root / rel
        need(p.is_file() and not p.is_symlink(), 'SOURCE_TYPE:' + rel)
        need(all(not q.is_symlink() for q in p.absolute().parents), 'SOURCE_ANCESTOR_SYMLINK')
        need(sha(p) == pin, 'SOURCE_HASH:' + rel)


def expanded_F(m, b):
    """14-term Laurent expansion, separately entered without factorized calls."""
    D, M, L, E, C, B = (m[k] for k in ('D', 'M', 'L', 'E', 'C', 'B'))
    Mt, Lt, Ct, Ey, My = (m[k] for k in ('Mt', 'Lt', 'Ct', 'Ey', 'My'))
    rho, c, theta = (b[k] for k in ('rho', 'c', 'theta'))
    R = rho * c
    return (D*Mt/(4*rho) + M*Mt/rho - L*Mt/(4*rho) - M*Lt/(2*rho)
            - E*D/(2*R) - E*L/(2*R) + C*Ct/(2*rho) - theta*C*C/(2*rho)
            - C*Ey/2 + COMPACT_B_SIGN*B*C/(2*R) - R*My*Ct/(4*rho)
            + R*theta*C*My/(4*rho) + R*My*Ey/4 - COMPACT_B_SIGN*B*My/4)


def four_contributions(m, b):
    """Own unsimplified ADM/Codazzi sum, including screen-volume variation."""
    D, M, L, E, C, B = (m[k] for k in ('D', 'M', 'L', 'E', 'C', 'B'))
    Mt, Lt, Ct, Ey, My = (m[k] for k in ('Mt', 'Lt', 'Ct', 'Ey', 'My'))
    rho, c, theta, a = (b[k] for k in ('rho', 'c', 'theta', 'a'))
    R = rho * c
    H = theta/rho
    small_l = Lt/(2*rho) + H*D/2 - a*B
    small_m = Mt/(2*rho) - E/R
    beta = small_m + H*M
    z = C/R - a*M - My/2
    eta = R*(Ct-theta*C)/(2*rho) - R*Ey/2 - B/2
    return {'curl_K2': D*beta/2 + B*z - E*(L-2*M)/R,
            'screen_connection': -(L/2-M)*small_m,
            'epsilon': M*(-small_l+small_m) if EPSILON_TERM_ENABLED else 0,
            'radial': (C/R-My/2)*eta}


def own_linear(ux, b):
    # Parameter derivative k, ordinary t clock, then proper observer derivative.
    ut = ux/b['ell']
    return LINEAR_PARAMETER_FACTOR * 3*ut/(2*b['rho']**2*b['c'])


def symbolic_audit():
    import sympy as s
    names = 'B C D E L M Mt Lt Ct Ey My rho c theta a ell ux'
    z = dict(zip(names.split(), s.symbols(names)))
    m = {k:z[k] for k in ('B','C','D','E','L','M','Mt','Lt','Ct','Ey','My')}
    b = {k:z[k] for k in ('rho','c','theta','a','ell')}
    need(s.cancel(sum(four_contributions(m,b).values()) - expanded_F(m,b)) == 0,
         'SYMBOLIC_FOUR_CONTRIBUTION_IDENTITY')
    need(s.diff(expanded_F(m,b),z['a']) == 0, 'CANCELLED_RADIAL_BACKGROUND')
    need(s.cancel(own_linear(z['ux'],b)-3*z['ux']/(2*z['ell']*z['rho']**2*z['c'])) == 0,
         'SECOND_PARAMETER_AND_TIME_UNITS')
    # The conversion is an operator identity, not differentiation of boxes.
    rho, ell, a, jr, jrr = tuple(z[k] for k in ('rho','ell','a')) + tuple(s.symbols('jr jrr'))
    rho_r = -ell*rho**2*a
    L = 1/(ell*rho)
    # D_y(-L*j_r)=L*(L_r*j_r+L*j_rr).
    yy_by_chain_rule = L**2*jrr + L*s.diff(L,rho)*rho_r*jr
    need(s.cancel(yy_by_chain_rule - (jrr/(ell**2*rho**2)+a*jr/(ell*rho))) == 0,
         'PHYSICAL_YY_CONNECTION')
    return {'four_contribution_equals_expanded_F':True,
            'a_cancels_exactly':True,
            'physical_yy_chain_rule':True,
            'first_metric_derivative_order':1,
            'second_response_derivative_order':1,
            'expanded_term_count':14}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--source-root',required=True,type=Path)
    ap.add_argument('--output',required=True,type=Path)
    args = ap.parse_args()
    need(not args.output.exists() and not args.output.is_symlink(), 'FRESH_REPORT')
    authenticate(args.source_root)
    own_pin = sha(__file__)
    rel = 'research/HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026/theory/screen_coefficient.py'
    spec = importlib.util.spec_from_file_location('authenticated_old_screen_map',args.source_root/rel)
    canonical = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(canonical)
    symbolic = symbolic_audit()
    keys = ('A','B','C','D','E','L','M','Mt','Lt','Ct','Ey','My')
    cases = []
    for i in range(1,33):
        m = {k:Q(((i+2)*(j+3)*7)%31-15, (i+3*j)%9+1) for j,k in enumerate(keys)}
        b = {'rho':Q(3*i+7,5),'c':Q(i+11,7),'theta':Q(i-15,37),
             'a':Q(i+2,23),'ell':Q(1,1000000)}
        ux = Q(3*i-13,37)
        expected = expanded_F(m,b)
        need(expected == sum(four_contributions(m,b).values()), 'EXACT_FOUR_CONTRIBUTION_CASE')
        got = canonical.coefficients(m,b,ux)
        need(got['F'] == expected, 'CANONICAL_MAP_EXPANDED_CASE')
        need(got['linear'] == own_linear(ux,b), 'CANONICAL_LINEAR_CASE')
        need(got['quadratic'] == 4*expected/(b['rho']*b['c']), 'PARAMETER_QUADRATIC_FACTOR')
        need(got['minus'] == -got['plus'], 'CURL_ORIENTATION')
        cases.append({'case':i,'F':str(expected),'plus':str(got['plus'])})
    zero = {k:Q(0) for k in keys}
    b = {'rho':Q(1),'c':Q(1),'theta':Q(0),'a':Q(7,3),'ell':Q(1,1000000)}
    controls = [({'D':1,'Mt':1},Q(1),'lapse_time_jet'),
                ({'E':1,'D':1},Q(-2),'angular_shift_lapse'),
                ({'C':1,'B':1},Q(2),'radial_time_shift'),
                ({'M':1,'Lt':1},Q(-2),'epsilon_volume'),
                ({'B':1,'My':1},Q(-1),'radial_shift_Hessian_jet')]
    control_results = []
    for updates, expected, name in controls:
        m = dict(zero);m.update({k:Q(v) for k,v in updates.items()})
        got = 4*expanded_F(m,b)
        need(got == expected, 'NONZERO_CONTROL:'+name)
        need(4*sum(four_contributions(m,b).values()) == expected, 'GEOMETRIC_CONTROL:'+name)
        control_results.append({'name':name,'quadratic_exact':str(got)})
    need(own_linear(Q(1,1000000), b) == Q(3,2), 'PHYSICAL_CLOCK_CONTROL')
    authenticate(args.source_root)
    need(sha(__file__) == own_pin, 'CHECKER_CHANGED')
    report = {'status':'PASS_PRIVATE_INDEPENDENT_SCREEN_FORMULA_AUDIT',
              'source_pins':SOURCE_PINS,'checker_sha256':own_pin,
              'authenticated_before_import_and_after':True,
              'imports_main_scientific_producers':False,
              'symbolic_audit':symbolic,'deterministic_rational_cases':cases,
              'nonzero_exact_controls':control_results,
              'conventions':{'parameter':'g=g0+lambda*q+lambda^2*k/2; coefficient is H_second_derivative',
                             'clock':'ut=ux/ell; e0u=ux/(ell*rho)',
                             'radial':'partial_y=-(ell*rho)^-1 partial_r',
                             'derivatives':'Ct,Ey,Lt,Mt,My are ordinary physical coordinate derivatives',
                             'orientation':'outward ambient sphere EPS_0123=+1; reverse orientation reverses h',
                             'observer':'future t-slice normal; radial normal is positive spatial y normal',
                             'screen':'coordinate t,y surfaces and positive symmetric inverse-square-root frame'},
              'required_first_field_slots':['B','C','D','E','L','M','Ct','Ey','Lt','Mt','My'],
              'scope':'Exact formula/unit audit only; no evolved postcontact value, sign, gauge-invariance, nonlinear existence or new protected-point result.'}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:
        f.write(json.dumps(report,sort_keys=True,indent=2)+'\n')
    print(report['status'])


if __name__ == '__main__':
    main()
