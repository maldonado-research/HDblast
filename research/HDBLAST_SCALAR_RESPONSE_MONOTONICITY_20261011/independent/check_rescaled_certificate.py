#!/usr/bin/env python3
"""Independent SymPy reconstruction of the October 11 scalar certificate.

The producer is never imported or executed. Only its certificate and source
bytes (for provenance hashing) are read. Every acceptance decision is exact.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path

import sympy as S

R = S.Rational
x = S.Symbol('x')
HERE = Path(__file__).resolve().parent
WORK = HERE.parent
DEFAULT_CERTIFICATE = WORK / 'bounds/RESCALED_BOUND_CERTIFICATE.json'
DEFAULT_PRODUCER = WORK / 'bounds/evaluate_rescaled_bounds.py'
DEFAULT_PROOF = WORK / 'analytic/JACOBI_ENDPOINT_MONOTONICITY.md'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def equal(actual, expected, label):
    if R(actual) != expected:
        raise ValueError(f'{label}: exact arithmetic mismatch')


def truth(actual, expected, label):
    if actual is not bool(expected):
        raise ValueError(f'{label}: Boolean mismatch')


def bound(poly, radius):
    poly = S.Poly(poly, x)
    return sum(abs(poly.nth(i)) * radius**i for i in range(poly.degree() + 1)) if not poly.is_zero else R(0)


def lower(poly, radius):
    poly = S.Poly(poly, x)
    return poly.nth(0) - sum(abs(poly.nth(i)) * radius**i for i in range(1, poly.degree() + 1)) if not poly.is_zero else R(0)


def reconstruct(row, model='registered'):
    p = row['parameters']
    coupling = 2/R('1.0357712571566784') - R(4, 3)
    actual_parameters = dict(k=R(1,9), w=R(2), cubic_W=R(2),
                             f0=1+coupling, f1=coupling, f2=R(0))
    if model in ('binary64_coupling','binary64_stored'):
        # Fixed independently inspected exact coefficients, not an arithmetic
        # claim equating either model with the registered exact-decimal model.
        actual_parameters['f1']=R(21026002552427,35184372088832)
        actual_parameters['f0']=1+actual_parameters['f1']
        if model=='binary64_stored':
            actual_parameters['k']=R(2001599834386887,18014398509481984)
            actual_parameters['f0']=R(56210374641259,35184372088832)
    elif model!='registered':
        raise ValueError('Unknown independently specified model')
    for key, value in actual_parameters.items():
        equal(p[key], value, 'parameters.'+key)
    k, w = actual_parameters['k'], actual_parameters['w']
    cubic = actual_parameters['cubic_W']
    f0, f1, f2 = actual_parameters['f0'], actual_parameters['f1'], actual_parameters['f2']
    radius, tube, delta = R(row['r']), R(row['b']), R(row['delta_star'])
    if min(radius, tube, delta) <= 0:
        raise ValueError('Nonpositive input radius, tube, or interval')
    reported = row['exact_constants']
    s, L = R(reported['s']), R(reported['L'])
    if s <= 0 or L < abs(f1)/w:
        raise ValueError('Invalid norm scale or scalar root bound')
    W = 3*k + w*x*x/2 + cubic*x**3/6
    f = f0+f1*x+f2*x*x/2
    Wp, fp = S.diff(W,x), S.diff(f,x)
    U = S.expand(Wp**2/2 - R(2,3)*W**2)
    A = S.expand(W*f/9 - Wp*fp/12)
    B_fun = S.expand(f*f/36 - fp*fp/48)
    lam = w-4*k
    m2 = w*lam
    nu = lam/2
    mu = nu*(nu+4*k)
    gap = m2-mu
    zeta_lower = mu*R(15,16)/(R(16,3)*k+nu)
    J_upper = 2/k+1/(2*zeta_lower)
    M, B = max(R(1),lam/s), 2*max(R(1),lam/s)
    C_G = max(1/gap,(1+m2/gap)/(4*k*s))
    M2, M3 = bound(S.diff(U,x,2),radius),bound(S.diff(U,x,3),radius)
    L1 = M3/2
    C_F = s*s/4+M2/12
    epsilon = C_F*B*B*tube*tube*J_upper/(2*k)
    C_N = L1*B*B + 4*C_F*s*B**3*tube/k
    L_N = 2*L1*B*tube + 20*C_F*s*B*B*tube*tube/k
    C_P = s*C_G*B*(2*L1*B+20*C_F*s*B*B*tube/k)
    C_D = s*C_G*C_N
    C_a = 4*C_F*B*B/k
    C_L_upper = 8*lam + R(28,3)*k*lam/(w-2*k)
    C_H_accepted = 4*C_L_upper/(9*k*k)
    C_H_riccati = 16*lam/(9*k*(w-3*k))
    C_H_convex_riccati = lam/(k*(w-3*k))
    endpoint_mode = reported['endpoint_mode']
    C_H = {'accepted':C_H_accepted,'riccati':C_H_riccati,
           'convex_riccati':C_H_convex_riccati}[endpoint_mode]
    W2, W3 = bound(S.diff(W,x,2),radius),bound(S.diff(W,x,3),radius)
    F1, F2 = bound(fp,radius),bound(S.diff(f,x,2),radius)
    F_abs, W_abs = bound(f,radius),bound(W,radius)
    a1 = k*f0/3
    C_K = L*bound(S.diff(A,x),radius)+bound(B_fun,radius)
    a_bar = R(4,3)*k+C_F*B*B*tube*tube/k
    a_min = k-C_F*B*B*tube*tube/k
    C_Tp = M2+4*a_bar*s*B+s*B*(lam+C_P*tube)
    C_Tbeta = 2*C_Tp/w
    C_slope_retired = s*B*C_a*L*L + ((C_a+W2/3)*L+F1/6)*C_Tbeta*L
    C_v = 8*C_F*B*B*J_upper/(3*k)
    C_Q_coarse = 2*C_v/k
    C_Q = C_v/a_min
    d0 = 2*(w-2*k)
    alpha = -f1/(2*d0)
    K_eta = (C_H*(3*a1/2)*L+(C_D+W3/2)*L*L+F2*L/2)/d0
    K_H = abs(S.Poly(A,x).nth(1))*K_eta + bound(S.diff(A,x,2),radius)*L*L/2
    K_H += bound(S.diff(B_fun,x),radius)*L
    suppression_q = k*f1*f1/(24*(w-2*k))
    expected = {name:val for name,val in locals().copy().items()
                if isinstance(val, S.Rational)}
    # Compare every scalar mathematical constant shared by the published schema.
    compared = []
    for name,value in expected.items():
        if name in reported and not isinstance(reported[name],list):
            equal(reported[name],value,'constant.'+name)
            compared.append(name)
    for name,poly in [('W',W),('f',f),('Wp',Wp),('fp',fp),('U',U),('A',A),('B_fun',B_fun)]:
        recorded = sum(R(coefficient)*x**i for i,coefficient in enumerate(reported[name]))
        if S.expand(recorded-poly) != 0:
            raise ValueError('polynomial.'+name+': symbolic mismatch')
    return locals().copy(), compared


def inspect_gate(item,lhs,rhs,strict):
    equal(item['lhs'],lhs,item['name']+'.lhs')
    equal(item['rhs'],rhs,item['name']+'.rhs')
    truth(item['strict'],strict,item['name']+'.strict')
    truth(item['pass_check'],lhs<rhs if strict else lhs<=rhs,item['name']+'.pass')
    equal(item['exact_margin'],rhs-lhs,item['name']+'.margin')


def verify_row(row, model='registered'):
    c, compared = reconstruct(row,model)
    k,B,b,r,d = c['k'],c['B'],c['tube'],c['radius'],c['delta']
    theta = c['C_Q']*(c['C_Tbeta']+c['s']*B)*c['L']**2*d*d
    gates = {
      'field_tube_inside_polynomial_radius':(B*b,r,False),
      'metric_volterra_contraction':(c['epsilon'],R(1,4),False),
      'scalar_contraction':(c['C_G']*c['L_N'],R(1,2),False),
      'scalar_ball_preservation':(c['C_G']*c['C_N']*b,c['M'],False),
      'positive_metric_derivative':(c['C_F']*B*B*b*b/k,k/2,False),
      'scalar_junction_derivative_feasibility':((c['C_P']+c['W3'])*b,c['w']/2,True),
      'scalar_root_inside_boundary_tube':(c['L']*d,b,False),
      'scalar_junction_derivative':((c['C_P']+c['W3'])*b+d*c['F2']/2,c['w']/2,False),
      'positive_reduced_curvature':(c['C_K']*d,c['a1']/2,False),
      'metric_residual_positive_at_Y':(3*c['a1']*d/2,R(9,32)*k*k,True),
      'positive_shell_tension':(c['W2']*b*b/6+d*c['F_abs']/6,k/2,False),
      'strict_downward_shell_crossing_factored_H2':(theta,R(1),True),
      'convex_bulk_potential_for_sharp_endpoint_bound':(R(0),lower(S.diff(c['U'],x,2),r),False),
      'convex_shell_tension_for_structural_monotonicity':(R(0),lower(2*S.diff(c['W'],x,2),r)+min(R(0),d*c['f2']),False),
    }
    seen = set()
    for item in row['checks']:
        name = item['name']
        if name not in gates or name in seen:
            raise ValueError('Unknown or duplicate gate: '+name)
        inspect_gate(item,*gates[name]); seen.add(name)
    required = set(gates)-{'strict_downward_shell_crossing_factored_H2',
                            'convex_bulk_potential_for_sharp_endpoint_bound',
                            'convex_shell_tension_for_structural_monotonicity'}
    if not required <= seen:
        raise ValueError('Missing construction gate')
    truth(row['pass_all_sufficient_gates'],all(i['pass_check'] for i in row['checks']),
          'all construction gates')
    # Independently establish the automatic convexity actually used by proof.
    if c['C_G']*c['L_N']<=R(1,2):
        if not c['M3']*B*b<=c['gap']/2<c['m2']:
            raise ValueError('Automatic convexity implication failed')
    E = min(c['L'],abs(c['alpha'])+c['K_eta']*d)
    pa,pb = S.Poly(c['A'],x),S.Poly(c['B_fun'],x)
    Kpoly = abs(pa.nth(1))*c['K_eta']
    Kpoly += sum(abs(pa.nth(i))*E**i*d**(i-2) for i in range(2,pa.degree()+1))
    Kpoly += sum(abs(pb.nth(i))*E**i*d**(i-1) for i in range(1,pb.degree()+1))
    Kbest = min(c['K_H'],Kpoly)
    equal(row['E_beta_over_delta'],E,'refined E')
    equal(row['K_H_polynomial'],Kpoly,'refined curvature remainder')
    equal(row['K_H_best'],Kbest,'selected curvature remainder')
    equal(row['absolute_errors_at_delta_star']['scalar'],c['K_eta']*d*d,'scalar error')
    equal(row['absolute_errors_at_delta_star']['curvature'],Kbest*d**3,'curvature error')
    sign = row['strict_curvature_suppression']
    margin = c['suppression_q']-Kbest*d
    equal(sign['q'],c['suppression_q'],'suppression q')
    equal(sign['exact_margin'],margin,'suppression margin')
    truth(sign['pass_check'],c['suppression_q']>0 and margin>0,'strict suppression')
    if 'retired_slope_gate' in row:
        inspect_gate(row['retired_slope_gate'],c['C_slope_retired']*d,c['a1']/4,False)
    if 'redundant_factored_slope_gate' in row:
        inspect_gate(row['redundant_factored_slope_gate'],theta,R(1),True)
    if 'additional_convexity_bounds' in row:
        diagnostic=row['additional_convexity_bounds']
        U_lower=lower(S.diff(c['U'],x,2),r)
        sigma_lower=lower(2*S.diff(c['W'],x,2),r)+min(R(0),d*c['f2'])
        equal(diagnostic['U_second_derivative_lower'],U_lower,'U convexity diagnostic')
        equal(diagnostic['sigma_second_derivative_lower'],sigma_lower,'sigma convexity diagnostic')
        truth(diagnostic['pass_both'],U_lower>=0 and sigma_lower>=0,'convexity diagnostics')
    return dict(compared_scalar_constants=len(compared),constant_names=sorted(compared),
                polynomials=7,construction_gates=len(seen),K_eta=c['K_eta'],K_H=Kbest,
                suppression_margin=margin,delta=d,theta=theta)


def verify_document(doc):
    rows = [('primary_uniform_certificate',doc['primary_uniform_certificate'])]
    rows += [(f'sample_interval_certificates[{i}]',row)
             for i,row in enumerate(doc['sample_interval_certificates'])]
    for name in ('conservative_without_bulk_convexity','conservative_metric_ratio_16over9',
                 'legacy_loose_root_bound_comparison'):
        if name in doc: rows.append((name,doc[name]))
    if 'conservative_riccati_factor_16_over_9' in doc:
        rows.append(('conservative_riccati_factor_16_over_9',doc['conservative_riccati_factor_16_over_9']))
    reports = {name:verify_row(row) for name,row in rows}
    variants={'historical_binary64_coupling_as_exact_rational':'binary64_coupling',
              'historical_stored_binary64_coefficients_as_exact_rationals':'binary64_stored'}
    for name,model in variants.items():
        if name in doc:
            reports[name]=verify_row(doc[name]['certificate'],model)
    primary = reports['primary_uniform_certificate']
    equal(primary['delta'],R(1,500),'primary interval')
    if not doc['primary_uniform_certificate']['pass_all_sufficient_gates']:
        raise ValueError('Primary certificate fails a required gate')
    if primary['suppression_margin']<=0:
        raise ValueError('Primary strict curvature suppression not proved')
    for row in doc['sample_errors_using_uniform_constants']:
        d = R(row['delta'])
        equal(row['scalar'],primary['K_eta']*d*d,'uniform scalar sample bound')
        equal(row['curvature'],primary['K_H']*d**3,'uniform curvature sample bound')
    for label,value,upper in [('K_eta',primary['K_eta'],R('3.196')),
                              ('K_H',primary['K_H'],R('0.254'))]:
        if not value<upper:
            raise ValueError('Readable upper rounding fails: '+label)
    if not primary['suppression_margin']>R('0.000423'):
        raise ValueError('Readable strict sign lower bound fails')
    if 'readable_strict_bounds' in doc:
        expected=[(primary['K_eta'],R('3.196'),True),
                  (primary['K_H'],R('0.254'),True),
                  (R('0.000423'),primary['suppression_margin'],True)]
        if len(doc['readable_strict_bounds'])!=3:
            raise ValueError('Unexpected readable-bound inventory')
        for item,values in zip(doc['readable_strict_bounds'],expected):
            inspect_gate(item,*values)
    if 'conservative_strict_suppression_subinterval' in doc:
        sub=doc['conservative_strict_suppression_subinterval']
        equal(sub['delta_star'],R(1,500),'strict sign full interval')
        equal(sub['K_H'],primary['K_H'],'strict sign remainder')
        equal(sub['exact_positive_margin'],primary['suppression_margin'],'strict sign full margin')
        equal(sub['q'],R(doc['primary_uniform_certificate']['strict_curvature_suppression']['q']),
              'strict sign q')
    # With the old norm's fixed exact certificate this ratio is purely rational.
    old_delta = R(161839258930731,1482649651195872051200)
    return reports,primary['delta']/old_delta


def mutation_controls(doc):
    tests = []
    def changed(label,edit):
        altered = copy.deepcopy(doc); edit(altered)
        try: verify_document(altered)
        except (ValueError,KeyError,TypeError): tests.append(label); return
        raise RuntimeError('Corrupt certificate accepted: '+label)
    def primary(d): return d['primary_uniform_certificate']
    changed('halved_scalar_error_constant',lambda d:primary(d)['exact_constants'].__setitem__('K_eta','1/1'))
    changed('omitted_derivative_scale_factor',lambda d:primary(d)['exact_constants'].__setitem__('C_D','1/1'))
    changed('falsified_polynomial_coefficient',lambda d:primary(d)['exact_constants']['U'].__setitem__(2,'0/1'))
    changed('falsified_positive_sign_margin',lambda d:primary(d)['strict_curvature_suppression'].__setitem__('exact_margin','1/1'))
    changed('missing_construction_gate',lambda d:primary(d)['checks'].pop(2))
    changed('understated_remainder',lambda d:primary(d).__setitem__('K_H_best','0/1'))
    return tests


def pack(value):
    if isinstance(value,S.Rational): return f'{value.p}/{value.q}'
    if isinstance(value,dict): return {k:pack(v) for k,v in value.items()}
    if isinstance(value,list): return [pack(v) for v in value]
    return value


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--certificate',type=Path,default=DEFAULT_CERTIFICATE)
    parser.add_argument('--producer',type=Path,default=DEFAULT_PRODUCER)
    parser.add_argument('--proof',type=Path,default=DEFAULT_PROOF)
    parser.add_argument('--require-proof-pin',action='store_true')
    parser.add_argument('--output',type=Path,default=HERE/'INDEPENDENT_EXACT_CHECK.json')
    args=parser.parse_args()
    doc=json.loads(args.certificate.read_text())
    if doc['producer_sha256']!=digest(args.producer):
        raise ValueError('Producer source digest differs from certificate')
    binding=doc.get('proof_binding',{})
    if args.require_proof_pin and binding.get('runtime_proof_pin_checked') is not True:
        raise ValueError('A verified runtime proof pin is required')
    if binding.get('runtime_proof_pin_checked') is True and binding.get('proof_sha256')!=digest(args.proof):
        raise ValueError('Certificate analytic proof digest mismatch')
    reports,gain=verify_document(doc)
    controls=mutation_controls(doc)
    result=dict(status='PASS_INDEPENDENT_EXACT_ARITHMETIC',
                arithmetic='SymPy Rational and independent symbolic polynomial expansion; no producer import or execution.',
                certificate_sha256=digest(args.certificate),producer_sha256=digest(args.producer),
                proof_sha256=digest(args.proof),checker_sha256=digest(Path(__file__)),
                reports=reports,certificate_gain_over_october8=gain,
                readable_bounds_checked=dict(K_eta_upper='3.196',K_H_upper='0.254',
                                             suppression_margin_lower='0.000423'),
                corrupt_certificate_rejections=controls,
                scope='Exact sufficient inequalities and polynomial remainder arithmetic only. Analytic proof review is separately recorded; no ODE integration, stability or novelty certification.')
    args.output.write_text(json.dumps(pack(result),indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],rows_checked=len(reports),
                         corrupt_certificate_rejections=len(controls),
                         constants_per_primary=reports['primary_uniform_certificate']['compared_scalar_constants'])))


if __name__=='__main__':main()
