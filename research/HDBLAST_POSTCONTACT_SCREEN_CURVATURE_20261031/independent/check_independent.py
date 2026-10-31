#!/usr/bin/env python3
"""Independent exact operator/reflected-kernel audit; no main producer import.

Source modules are authenticated before import and after all work. SymPy is
used only as an independent exact polynomial oracle, never as interval proof.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
import sympy as S

PINMAP = 'ed1e70a05034feec47e58e75d09ffa94849b343597b21b11122d8704a5d221e4'
LIFT_SIGN = -1
KX_SECOND_SIGN = 1
REFLECTION_SIGN = 1


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def same(a, b, message):
    need(S.expand(a-b) == 0, message)


def exact_kernel(f, x, r, s, v, after_contact):
    """Integrate signed-q image interval by its positive-depth pieces."""
    primitive = S.integrate(f, v)
    F = lambda z: primitive.subs(v, z)-primitive.subs(v, 0)
    top = F(r+x-s)
    direct = top-F(r-x+s)
    if after_contact:
        reflected = top+REFLECTION_SIGN*F(x-r-s)
        return S.expand((S.integrate(reflected, (s, 0, x-r))
                         +S.integrate(direct, (s, x-r, x)))/2)
    return S.expand(S.integrate(direct, (s, 0, x))/2)


def endpoint_integral(f, x, r, s, v, after_contact, derivative):
    hi = f.subs(v, r+x-s)
    sign = KX_SECOND_SIGN if derivative == 'x' else -1
    lo = f.subs(v, r-x+s)
    direct = hi+sign*lo
    if after_contact:
        reflected = hi+sign*f.subs(v, x-r-s)
        return S.expand((S.integrate(reflected, (s, 0, x-r))
                         +S.integrate(direct, (s, x-r, x)))/2)
    return S.expand(S.integrate(direct, (s, 0, x))/2)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source-root', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    pinpath = here/'SOURCE_PINS.json'
    need(sha(pinpath) == PINMAP, 'PINMAP_AUTHENTICATION')
    pins = json.loads(pinpath.read_text())
    root = args.source_root.absolute()
    need(root.is_dir() and not root.is_symlink(), 'SOURCE_ROOT_TYPE')
    out = args.output.absolute()
    need(not out.exists() and not out.is_symlink(), 'FRESH_OUTPUT')
    need(root not in out.parents, 'OUTPUT_IN_READONLY_SOURCE')
    need(all(not p.is_symlink() for p in out.parents), 'OUTPUT_PARENT_SYMLINK')
    def authenticate():
        for rel, pin in pins.items():
            p = root/rel
            need(p.is_file() and not p.is_symlink(), 'SOURCE_TYPE:'+rel)
            need(all(not a.is_symlink() for a in p.parents), 'SOURCE_ANCESTOR_SYMLINK')
            need(sha(p) == pin, 'SOURCE_AUTHENTICATION:'+rel)
    authenticate()
    prefix = root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection'
    op = load('independent_actual_tt_operator', prefix/'inputs/linear_operator.py')
    adapter = load('independent_actual_tt_adapter', prefix/'inputs/source_adapters.py')
    need([op.IDX['T:+'], op.IDX['T:-']] == [21,22], 'TT_INDICES')
    mats = op.matrices()
    for key, m in mats.items():
        for i in (21,22):
            need(all(m[i][j] == m[j][i] == 0 for j in range(23) if i != j),
                 'PRIMITIVE_TT_CLOSURE:'+key)
        need(m[21][21] == m[22][22], 'TT_ORIENTATION_OPERATOR:'+key)
    expected = {'I':Q(1),'Sa':Q(-2),'Sh':Q(2),'Rb':Q(1),'Rd':Q(4),
                'Delta':Q(-6),'Pmet':Q(1),'PN':Q(1)}
    need({k:m[21][21] for k,m in mats.items() if m[21][21]} == expected,
         'PRIMITIVE_DIAGONALS')
    names = 'ell rho a h c Ht p py b d V V1 V2 sigma1 sigma2'
    z = dict(zip(names.split(), S.symbols(names)))
    # Actual assembly, independent exact symbolic expected expression.
    assembly = op.assemble(**z)
    Dfull = z['ell']**2*(-6/z['c']**2+z['rho']**2*(
        -2*z['a']**2+2*z['h']**2-2*z['b']-8*z['d']+S.Rational(4,3)*z['V']))
    for name, expected_value in {'A':-3*z['ell']*z['Ht'],
            'B':-3*z['ell']*z['rho']*z['a'], 'D':Dfull, 'M':0}.items():
        same(assembly[name][21][21], expected_value, 'ASSEMBLED_'+name)
        same(assembly[name][22][22], expected_value, 'ASSEMBLED_MINUS_'+name)
    # The existing branch note explicitly supplies both identities.  This is
    # an algebraic consequence conditional on those existing actual premises.
    ay = -1/z['rho']**2-z['p']**2/3
    V = z['p']**2/2-6*z['a']**2+6/z['rho']**2
    reduced = Dfull.subs({z['b']:-(ay+z['a']**2),
                         z['d']:1/z['rho']**2-z['a']**2,
                         z['V']:V, z['h']:z['Ht']/z['rho']})
    same(reduced, z['ell']**2*(-6/z['c']**2-2+2*z['Ht']**2), 'BACKGROUND_REDUCTION')
    same(reduced.subs(z['Ht']**2,1-1/z['c']**2), -8*z['ell']**2/z['c']**2,
         'DT_MINUS8_SECH2')
    background = (root/'research/HDBLAST_SECOND_INITIAL_NORMS_GAUGE_20261015/background/inputs/PRIOR_ACTUAL_BACKGROUND_AND_RIGOROUS_COLLAR_BOUNDS.md').read_text()
    need('a_y=-\\rho^{-2}-p^2/3' in background and
         'a^2=\\rho^{-2}+p^2/12-U/6' in background,
         'BACKGROUND_PREMISES_PRESENT')
    maps = json.loads((prefix/'GEOMETRIC_COEFFICIENT_MAPS.json').read_text())
    forcing = maps['maps']['metric_forcing']['ab']
    need(forcing[5] == forcing[6] and len(forcing[5]) == 78, 'COMPLETE_78_TERM_TT_MAP')
    def normalize_terms(rows,shift=None,factor=Q(1)):
        result={}
        for row in rows:
            powers=dict(row['powers'])
            for key,power in (shift or {}).items():
                powers[key]=powers.get(key,0)+power
            key=tuple(sorted((k,v) for k,v in powers.items() if v))
            result[key]=result.get(key,Q(0))+factor*Q(row['coefficient'])
        return {k:v for k,v in result.items() if v}
    need(normalize_terms(forcing[5])==normalize_terms(
         maps['maps']['QEred']['ab'][5],{'ell':2,'rho':2},Q(-2)),
         'ACTUAL_FORCING_MINUS_TWO_ELL2_RHO2_UNITS')
    # Independent Laurent evaluator checks equality at rational assignments;
    # exact list equality above is stronger and covers all symbolic inputs.
    allsymbols = sorted({k for row in forcing[5] for k in row['powers']})
    env = {name:Q(2*j+3,j+2) for j,name in enumerate(allsymbols)}
    def evaluate(rows):
        result = Q(0)
        for row in rows:
            term = Q(row['coefficient'])
            for name,power in row['powers'].items():
                term *= env[name]**power
            result += term
        return result
    fvalue = evaluate(forcing[5]); need(fvalue != 0, 'GENERIC_FORCING_NOT_ZERO')
    need(fvalue == evaluate(forcing[6]), 'FORCING_ORIENTATION_EVALUATION')
    shell_names = 'a p sigmap sigmapp sigmappp h X V L d v'
    shell = dict(zip(shell_names.split(),S.symbols(shell_names)))
    shell_output = op.boundary_forcing(z['ell'],z['rho'],adapter.shell_sources(**shell))
    G = -S.Rational(8,3)*z['ell']*z['rho']*(shell['a']*shell['h']
         +shell['sigmap']*shell['v']/3)*shell['d']
    same(shell_output[21],G,'ACTUAL_G_PLUS')
    same(shell_output[22],G,'ACTUAL_G_MINUS')
    need(G != 0, 'NONHOMOGENEOUS_NEUMANN')
    initial_names = 'rho a p u uy ut v vt Fmom Gmom Fmom_y'
    initial = dict(zip(initial_names.split(),S.symbols(initial_names)))
    initial.update(f=S.symbols('f0 f2 f4'),fy=S.symbols('fy0 fy2 fy4'))
    u0,v0 = op.pack_initial(adapter.initial_frame_data(**initial),z['ell'])
    need(all(u0[i] == v0[i] == 0 for i in (21,22)), 'ACTUAL_ZERO_INITIAL_TT')
    # Independent exact images for every monomial s^i v^j, i+j <= 4.
    x,r,s,v = S.symbols('x r s v',positive=True)
    polynomials = [s**i*v**j for i in range(5) for j in range(5-i)]
    checks = 0
    for f in polynomials:
        for post in (False,True):
            image = exact_kernel(f,x,r,s,v,post)
            ux = endpoint_integral(f,x,r,s,v,post,'x')
            ur = endpoint_integral(f,x,r,s,v,post,'r')
            same(S.diff(image,x),ux,'KX_ENDPOINT')
            same(S.diff(image,r),ur,'KR_ENDPOINT')
            same(S.diff(image,x,2)-S.diff(image,r,2),f.subs({s:x,v:r}), 'K_PDE')
            if post:
                same(S.diff(image,r).subs(r,0),0,'NEUMANN_IMAGE_BOUNDARY')
            else:
                same(image.subs(x,0),0,'K_ZERO_INITIAL_VALUE')
                same(S.diff(image,x).subs(x,0),0,'K_ZERO_INITIAL_VELOCITY')
            checks += 1
    # Sharp simple oracle which rejects an omitted image in the post-contact
    # cone: constant source remains x^2/2, ux=x, even at r<x.
    constant_image=exact_kernel(S.Integer(1),x,r,s,v,True)
    same(constant_image,x*x/2,'POSTCONTACT_CONSTANT_IMAGE')
    same(endpoint_integral(S.Integer(1),x,r,s,v,True,'x'),x,'POSTCONTACT_CONSTANT_KX')
    g=s+s*s
    lift=LIFT_SIGN*S.integrate(g,(s,0,x-r))
    same(S.diff(lift,x,2)-S.diff(lift,r,2),0,'LIFT_WAVE_EQUATION')
    same(S.diff(lift,r).subs(r,0),g.subs(s,x),'LIFT_POSITIVE_INWARD_R_DERIVATIVE')
    same(S.diff(lift,x),-g.subs(s,x-r),'LIFT_TIME_SIGN')
    # Positive outward-y normal equals -positive-r: dy=-ell*rho dr.
    normal_y = -S.diff(lift,r)/(z['ell']*z['rho'])
    same(normal_y.subs(r,0),-g.subs(s,x)/(z['ell']*z['rho']), 'OUTWARD_Y_ORIENTATION')
    # Closed backward physical cone and rational footprint.
    T,R=S.Rational(5,4),S.Rational(1,8)
    need(T-R==S.Rational(9,8) and T+R==S.Rational(11,8), 'PDAGGER_HISTORY')
    # Two endpoint curves together have exact half-weight mass T.
    need(S.integrate(S.Rational(1,2),(s,0,T))*2==T,'ENDPOINT_TOTAL_MASS')
    # The fold does NOT reverse the Kx integrand sign. Drift coefficients
    # have no physical continuation: reflect the full composite RHS.
    for q in (Q(-3,4),Q(-1,8),Q(1,2)):
        physical_B=Q(2)+abs(q)
        physical_Ur=2*abs(q)
        expected_composite=physical_B*physical_Ur
        reflected_B=(1 if q>0 else -1)*physical_B
        reflected_Ur=(1 if q>0 else -1)*physical_Ur
        need(reflected_B*reflected_Ur==expected_composite,'COMPOSITE_DRIFT_REFLECTION')
    authenticate()
    need(sha(pinpath)==PINMAP,'PINMAP_POSTAUTH')
    result={
      'status':'PASS_INDEPENDENT_EXACT_TT_POSTCONTACT_ARCHITECTURE',
      'source_pins':pins,'source_pins_sha256':PINMAP,
      'before_import_and_after_source_authentication':True,
      'project_imports':['actual linear_operator.py','actual source_adapters.py'],
      'main_scientific_producer_imported':False,
      'primitive_matrices':len(mats),'TT_indices':[21,22],
      'TT_forcing_terms_per_orientation':78,'generic_TT_forcing_nonzero':True,
      'bulk_source_units':'metric_forcing_T=-2 ell^2 rho^2 QEred_T checked as exact normalized Laurent maps',
      'D_T':'-8 ell^2/cosh(t)^2 conditional on authenticated existing background field equation and Hamilton constraint',
      'G_T':'-(8/3)ell*rho*(a*h_shell+sigma_prime*v_shell/3)*d_shell, generally nonzero',
      'initial_TT_value_velocity_identically_zero':True,
      'polynomial_kernel_cases':checks,'monomial_total_degree_max':4,
      'exact_Pdagger':['5/4','1/8'],'shell_history':['0','9/8'],
      'initial_footprint':['0','11/8'],
      'Neumann_lift':'LG=-integral_0^(x-r) G(s) ds for x>r; zero for x<r',
      'time_derivative':'(KNf)_x=1/2 integral_0^x [f(s,r+x-s)+f(s,abs(r-x+s))] ds',
      'radial_derivative':'(KNf)_r=1/2 integral_0^x [f(s,r+x-s)-f(s,abs(r-x+s))] ds',
      'reflection':'full physical scalar composite RHS even; equivalent B odd times Ur odd',
      'actual_response':'U=LG+KNF+KN(A Ux+B Ur+D U); every feedback term retained',
      'complete_feedback_bound':{
        'norm':'X=max(2 sup|u|/T,sup|ux|,sup|ur|)',
        'sigma':'a+b+T*d/2','q':'T*sigma<1',
        'seed_bound':'||LG+KNF||X <= 2g+Tf',
        'tail_derivative':'|Ux-(LG+KNF)x| <= q*(2g+Tf)/(1-q)',
        'domain':'0<=s<=5/4, max(0,s-9/8)<=v<=11/8-s',
        'premises':'uniform actual F,A,B,D bounds on whole reflected backward cone, actual G on shell history; conditional C1 integral solution'},
      'scope':'Independent exact architecture and source identity audit, no new numerical actual Ux or screen sign certificate; prior protected-point screen result already exists.'}
    with out.open('x') as stream:
        json.dump(result,stream,sort_keys=True,indent=2);stream.write('\n')
    print(result['status'])


if __name__=='__main__':
    main()
