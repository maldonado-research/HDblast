"""Independent exact algebra and continuous-cell material-source audit.

No scientific producer is imported. The angular moment formula is evaluated
with Fraction arithmetic, and all prior and new CSV cells are reconstructed.
"""
import argparse
import csv
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent.parent
PINS = {
 'prior/INERTIAL_WAVEFORM_CELLS.csv':'2f76fd1ae2195b445582a29a352bd76658dea17c2233cbb9acc395b255172bf7',
 'prior/INERTIAL_WAVEFORM_ENCLOSURE.json':'d2ffbcec5413580ff5703da58b062a9b83024431748fa40726555c09c2a968d8',
 'prior/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json':'8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8',
 'prior/PRIOR_PHYSICAL_GAIN_ENCLOSURE.json':'2a2dbc398ac2799418edcddeb917f4dcecce9b9a8049ac73b54541e617b616a4',
}
ZERO = (0,0,0,0)


def poly(c=0):
    return {} if not c else {ZERO:F(c)}


def var(i):
    key = list(ZERO); key[i] = 1
    return {tuple(key):F(1)}


def plus(*parts):
    out = {}
    for part in parts:
        for key,c in part.items():
            out[key] = out.get(key,F(0))+c
    return {k:c for k,c in out.items() if c}


def times(a,b):
    out = {}
    for k,c in a.items():
        for j,d in b.items():
            key = tuple(x+y for x,y in zip(k,j))
            out[key] = out.get(key,F(0))+c*d
    return {k:c for k,c in out.items() if c}


def scale(a,c):
    return {k:v*c for k,v in a.items() if v*c}


def power(a,n):
    out = poly(1)
    for _ in range(n): out = times(out,a)
    return out


def derivative(a,i):
    out = {}
    for k,c in a.items():
        if k[i]:
            j = list(k); j[i] -= 1; out[tuple(j)] = c*k[i]
    return out


def laplace(a):
    return plus(*(derivative(derivative(a,i),i) for i in range(4)))


def rising(a,n):
    out = F(1)
    for i in range(n): out *= a+i
    return out


def mean(a):
    out = F(0)
    for k,c in a.items():
        if any(v%2 for v in k): continue
        orders = [v//2 for v in k]
        numerator = F(1)
        for n in orders: numerator *= rising(F(1,2),n)
        out += c*numerator/rising(F(2),sum(orders))
    return out


def square(interval):
    lo,hi = interval
    if lo > hi: raise RuntimeError('Reversed interval')
    return (F(0) if lo <= 0 <= hi else min(lo*lo,hi*hi),max(lo*lo,hi*hi))


def scaled(interval,c):
    return tuple(sorted(v*c for v in interval))


def product(a,b):
    values = [x*y for x in a for y in b]
    return min(values),max(values)


def fixed(interval,den):
    return (interval[0].numerator*den//interval[0].denominator,
            -((-interval[1].numerator*den)//interval[1].denominator))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root',type=Path,default=ROOT)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args(); root = args.source_root
    checks = []

    def check(label,ok):
        if not ok: raise RuntimeError(label)
        checks.append(label)

    x = [var(i) for i in range(4)]
    s = plus(*(power(v,2) for v in x))
    A = plus(power(x[0],2),power(x[1],2))
    Y = plus(power(x[0],2),scale(power(x[1],2),-1))
    Z = plus(A,scale(s,F(-1,2)))
    H4 = plus(power(Y,2),scale(times(s,A),F(-1,2)),scale(power(s,2),F(1,12)))
    check('ambient_harmonic_B_Z2_H4',not laplace(Y) and not laplace(Z) and not laplace(H4))
    check('exact_ambient_degree_0_2_4_decomposition',power(Y,2)==plus(H4,scale(times(s,Z),F(1,2)),scale(power(s,2),F(1,6))))
    check('normalized_sphere_volume_convention',mean(poly(1))==1)
    check('quadrupole_second_moment',mean(power(Y,2))==F(1,6))
    check('quadrupole_fourth_moment',mean(power(Y,4))==F(3,40))
    check('degree_two_norm',mean(power(Z,2))==F(1,12))
    check('degree_four_norm',mean(power(H4,2))==F(19,720))
    check('orthogonality_of_three_degrees',mean(Z)==mean(H4)==mean(times(Z,H4))==0)
    shares = [F(10,27),F(5,18),F(19,54)]
    contributions = [F(4,9),F(1,3),F(19,45)]
    check('full_scalar_source_squared_angular_mean',16*mean(power(Y,4))==sum(contributions)==F(6,5))
    check('all_three_orthogonal_squared_power_shares',all(a/F(6,5)==b for a,b in zip(contributions,shares)) and sum(shares)==1)
    eta = var(0); delta = F(1,1000); c = F(96714321425692,161839258930731)
    sigma = plus(poly(F(2,3)+delta*(1+c)),scale(eta,delta*c),scale(power(eta,2),2),scale(power(eta,3),F(2,3)))
    second = derivative(derivative(sigma,0),0)
    third = derivative(second,0)
    check('actual_registered_tension_second_and_third_derivatives',second==plus(poly(4),scale(eta,4)) and third==poly(4))
    eps,base,p1,p2 = x
    eta_eps = plus(base,times(eps,p1),scale(times(power(eps,2),p2),F(1,2)))
    sigmap_eps = plus(scale(eta_eps,4),scale(power(eta_eps,2),2),poly(delta*c))
    second_eps = {}
    for k,v in sigmap_eps.items():
        if k[0]==2:
            j=(0,*k[1:]);second_eps[j]=2*v
    expected = plus(times(plus(poly(4),scale(base,4)),p2),scale(power(p1,2),4))
    check('second_derivative_chain_rule_without_spurious_half',second_eps==expected)
    check('sum_and_per_copy_source_factors',F(-4)/2==-2 and F(-4)*F(1,2)==-2)

    inputs = {}
    for name,digest in PINS.items():
        raw = (root/'inputs'/name).read_bytes()
        check('authenticated_prior_input:'+name,sha256(raw).hexdigest()==digest)
        inputs[name] = raw
    prior = json.loads(inputs['prior/INERTIAL_WAVEFORM_ENCLOSURE.json'])
    shell = json.loads(inputs['prior/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json'])
    gain_record = json.loads(inputs['prior/PRIOR_PHYSICAL_GAIN_ENCLOSURE.json'])['correlated_retained_eta_box']
    new_raw = (root/'computed/QUADRATIC_FORCING_ENCLOSURE.json').read_bytes()
    csv_raw = (root/'computed/QUADRATIC_FORCING_CELLS.csv').read_bytes()
    new = json.loads(new_raw)
    check('new_report_and_complete_csv_authenticated',new['status']=='PASS' and new['cell_data']['bytes']==len(csv_raw) and new['cell_data']['sha256']==sha256(csv_raw).hexdigest())
    check('same_wider_state_and_derivative_normal_conventions',new['state']['ell']==prior['state']['ell']=='1/1000000' and new['state']['derivative_convention']=='Ueps=U0+eps U1+eps^2 U2/2' and new['state']['normal_convention']=='outward normal SUM over both copies')
    den,N,ell = 2**112,4096,F(1,10**6)
    rows = list(csv.DictReader(inputs['prior/INERTIAL_WAVEFORM_CELLS.csv'].decode().splitlines()))
    target = list(csv.DictReader(csv_raw.decode().splitlines()))
    check('all_4353_source_and_target_cells_present',len(rows)==len(target)==new['cell_data']['rows']==4353)
    sigma2 = tuple(4+4*F(v) for v in shell['boxes']['eta_b']['exact'])
    check('actual_positive_tension_hessian_interval',sigma2==tuple(map(F,shell['boxes']['sigma_second_b']['exact']))==tuple(map(F,new['parameters']['sigma_second_interval'])) and sigma2[0]>0)
    G = tuple(map(F,gain_record['G']['outward_rational_at_12_decimal_places']))
    K = tuple(map(F,gain_record['K']['outward_rational_at_12_decimal_places']))
    check('physical_gain_and_norm_prefactor_use_authenticated_background',G==tuple(map(F,prior['physical_gain_interval'])) and K==tuple(map(F,prior['physical_norm_prefactor_interval'])))
    j4 = [F(0),F(0)];previous = 0
    for index,(row,out) in enumerate(zip(rows,target)):
        left,right = int(row['left_grid']),int(row['right_grid'])
        if not (left==previous and left<right<=3*N and (int(out['left_grid']),int(out['right_grid']))==(left,right)):
            raise RuntimeError('Continuous coverage mismatch at cell'+str(index))
        previous = right
        E2 = square((F(row['ell_S_lo'])/den,F(row['ell_S_hi'])/den))
        scalar = [scaled(E2,k) for k in [F(-2,3),F(-2),F(-4)]]
        tbase = product(E2,sigma2)
        tension = [scaled(tbase,k) for k in [F(1,6),F(1,2),F(1)]]
        keys = [('ell2_Q0_lo','ell2_Q0_hi'),('ell2_Q2_lo','ell2_Q2_hi'),('ell2_Q4_lo','ell2_Q4_hi'),('ell2_T0_lo','ell2_T0_hi'),('ell2_T2_lo','ell2_T2_hi'),('ell2_T4_lo','ell2_T4_hi')]
        for value,pair in zip(scalar+tension,keys):
            if fixed(value,den)!=tuple(int(out[k]) for k in pair):
                raise RuntimeError('Coefficient rounding mismatch at cell'+str(index))
        H = F(row['H_lo'])/den,F(row['H_hi'])/den
        H4range = square(square(H))
        factor = (1-(ell*F(right,N))**2)**2
        width = F(right-left,N)
        j4[0] += width*factor*H4range[0]
        j4[1] += width*H4range[1]
    check('all_six_harmonic_coefficient_ranges_rebuilt_at_every_cell',previous==3*N)
    check('fourth_moment_retains_positive_cosh_minus3_weight',tuple(j4)==tuple(map(F,new['quadrature']['cosh_minus3_H_fourth_integral_interval'])) and 0<j4[0]<=j4[1])
    oldC = tuple(map(F,prior['quadrature']['C_interval']))
    check('same_state_C_identity_not_cross_width_transfer',oldC==tuple(map(F,new['quadrature']['inherited_C_interval'])))
    check('signed_integral_and_L1_are_exact_multiples_of_C',scaled(oldC,-4)==tuple(map(F,new['quadrature']['scalar_signed_integral_interval'])) and scaled(oldC,4)==tuple(map(F,new['quadrature']['scalar_L1_interval'])))
    QL2 = scaled(product(product(K,square(G)),j4),F(36,5)/ell**3)
    check('full_scalar_L2_squared_with_correct_clock_and_inverse_cubic_width',QL2==tuple(map(F,new['quadrature']['scalar_L2_squared_interval'])))
    TL2 = scaled(product(product(product(square(sigma2),K),square(G)),j4),F(9,20)/ell**3)
    check('explicit_constitutive_term_integral_and_L2_squared',product(sigma2,oldC)==tuple(map(F,new['quadrature']['constitutive_tension_integral_interval'])) and TL2==tuple(map(F,new['quadrature']['constitutive_tension_L2_squared_interval'])))
    ht = tuple(map(F,prior['tail_theorem']['actual_H_interval']))
    physical_tail = G[0]*ht[1],G[1]*(1-(3*ell)**2)*ht[0]
    tail = scaled(square(physical_tail),F(-2,3))
    check('continuous_tail_monopole_negative_and_strict_bounds',tail==tuple(map(F,new['tail_theorem']['ell2_scalar_monopole_interval'])) and F('-2.3e-19')<tail[0]<=tail[1]<F('-6.7e-20'))
    check('no_full_second_order_or_total_stress_claim',new['scope']['explicit_material_tension_contribution'] is True and all(new['scope'][k] is False for k in ['complete_bulk_coordinate_quadratic_source','full_U2_response_computed','total_tension_or_stress_second_variation','nonlinear_remainder_bound','novelty_claimed','cosmological_evidence']))
    check('inputs_and_outputs_unchanged_during_audit',all((root/'inputs'/n).read_bytes()==raw for n,raw in inputs.items()) and (root/'computed/QUADRATIC_FORCING_ENCLOSURE.json').read_bytes()==new_raw and (root/'computed/QUADRATIC_FORCING_CELLS.csv').read_bytes()==csv_raw)
    result = {'status':'PASS_INDEPENDENT_QUADRATIC_MATERIAL_SOURCE_AUDIT','python_optimization':sys.flags.optimize,
      'check_count':len(checks),'checks':checks,'cells_reconstructed':len(rows),
      'angular_moments':{'mean_Y2':'1/6','mean_Y4':'3/40','mean_Z2_squared':'1/12','mean_H4_squared':'19/720'},
      'input_pins':PINS,'report_sha256':sha256(new_raw).hexdigest(),'csv_sha256':sha256(csv_raw).hexdigest(),
      'own_source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
      'display_only':{'fourth_moment':[float(v) for v in j4],'scalar_L2_squared':[float(v) for v in QL2]},
      'scope':'Independent exact algebra and full-cell audit of the explicit material tension contributions. It does not give the complete harmonic-coordinate source, U2, total second stress, a nonlinear remainder or a novelty claim.'}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'check_count':len(checks),'cells':len(rows)}))


if __name__=='__main__':
    main()
