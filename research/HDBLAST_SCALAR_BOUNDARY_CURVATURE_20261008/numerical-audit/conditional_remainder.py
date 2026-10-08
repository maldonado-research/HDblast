"""Exact rational curvature-remainder bound, conditional on a scalar-shift bound.

This module does NOT prove the boundary-value solution exists or derive K_eta.
Its output is an implication: if |eta-a*delta| <= K_eta*delta**2 for an exact
solution and 0 < delta <= delta_max, then |H2-P2(delta)| <= C_H*delta**3.
All polynomial arithmetic and bound constants below use fractions.Fraction.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path


def add(a, b):
    return [ (a[i] if i < len(a) else Q(0)) +
             (b[i] if i < len(b) else Q(0))
             for i in range(max(len(a),len(b))) ]


def scale(a, c):
    return [x*c for x in a]


def multiply(a, b):
    out = [Q(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            out[i+j] += x*y
    return out


def derivative(a):
    return [i*a[i] for i in range(1,len(a))]


def value(a, x):
    out = Q(0)
    for c in reversed(a):
        out = out*x+c
    return out


def curvature_polynomials(p):
    k,w,v,f0,f1,f2 = [p[x] for x in ('k','w','v','f0','f1','f2')]
    W,f = [3*k,Q(0),w/2,v/6], [f0,f1,f2/2]
    Wp,fp = derivative(W),derivative(f)
    A = add(scale(multiply(W,f),Q(1,9)),scale(multiply(Wp,fp),-Q(1,12)))
    B = add(scale(multiply(f,f),Q(1,36)),scale(multiply(fp,fp),-Q(1,48)))
    return A,B


def bound(p, delta_max, K_eta):
    if any(type(p[x]) is not Q for x in ('k','w','v','f0','f1','f2')):
        raise TypeError('Exact Fraction model coefficients required')
    if type(delta_max) is not Q or type(K_eta) is not Q:
        raise TypeError('Exact Fraction bound inputs required')
    k,w,f0,f1 = (p[x] for x in ('k','w','f0','f1'))
    if not (k > 0 and w > 4*k and f0 > 0 and delta_max > 0 and K_eta >= 0):
        raise ValueError('Outside the stated positive-curvature nondegenerate domain')
    a = -f1/(4*(w-2*k))
    E = abs(a)+K_eta*delta_max
    A,B = curvature_polynomials(p)
    C = abs(A[1])*K_eta
    C += sum(abs(A[n])*E**n*delta_max**(n-2) for n in range(2,len(A)))
    C += sum(abs(B[n])*E**n*delta_max**(n-1) for n in range(1,len(B)))
    c2 = f0*f0/36-k*f1*f1/(24*(w-2*k))
    if A[0] != k*f0/3 or A[1]*a+B[0] != c2:
        raise ArithmeticError('Response identity failed in exact arithmetic')
    return dict(a=a,E=E,A=A,B=B,C_H=C,coefficient_delta=k*f0/3,coefficient_delta2=c2)


def exact_text(x):
    if type(x) is Q:
        return str(x.numerator)+'/'+str(x.denominator)
    if isinstance(x,dict):
        return {k:exact_text(v) for k,v in x.items()}
    if isinstance(x,list):
        return [exact_text(v) for v in x]
    return x


def controls(models):
    """Check proof transcription and reject invalid inputs; no ODE claim."""
    failures = []
    covered_points = 0
    feedback_omission_detected = 0
    for p in models:
        D,K = Q(1,500),Q(1)
        b=bound(p,D,K)
        # Exact points including both error-envelope endpoints and its center.
        # The accompanying triangle-inequality proof gives continuum coverage.
        for divisor in (1,2,4,16,1000):
            d = D/divisor
            for signed_fraction in (Q(-1),-Q(1,2),Q(0),Q(1,2),Q(1)):
                eta=b['a']*d+signed_fraction*K*d*d
                h=d*value(b['A'],eta)+d*d*value(b['B'],eta)
                pred=b['coefficient_delta']*d+b['coefficient_delta2']*d*d
                if abs(h-pred) > b['C_H']*d**3:
                    failures.append('exact envelope point escaped')
                covered_points += 1
        if p['f1']:
            d=Q(1,10**6)
            b0=bound(p,d,Q(0))
            eta=b0['a']*d
            h=d*value(b0['A'],eta)+d*d*value(b0['B'],eta)
            wrong=p['k']*p['f0']*d/3+p['f0']**2*d*d/36
            if abs(h-wrong) <= b0['C_H']*d**3:
                failures.append('omitted scalar-feedback coefficient went undetected')
            else:
                feedback_omission_detected += 1
    invalids=[(dict(models[0],k=Q(0)),Q(1,500),Q(1)),
              (dict(models[0],w=4*models[0]['k']),Q(1,500),Q(1)),
              (models[0],Q(0),Q(1)),(models[0],Q(1,500),Q(-1))]
    invalid_rejections = 0
    for p,d,k in invalids:
        try:
            bound(p,d,k)
        except ValueError:
            invalid_rejections += 1
        else:
            failures.append('invalid domain accepted')
    if failures:
        raise ArithmeticError('; '.join(failures))
    return dict(exact_polynomial_envelope_points=covered_points,
                feedback_omission_negative_controls=feedback_omission_detected,
                invalid_domain_rejections=invalid_rejections,
                scope='Implementation controls. Continuum implication is proved by the displayed polynomial triangle inequality, not point sampling.')


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--replay-dir',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--delta-max',default='0.002')
    ap.add_argument('--conditional-scalar-error-bound',default='1')
    args=ap.parse_args()
    D,K=Q(args.delta_max),Q(args.conditional_scalar_error_bound)
    models=[]
    for name in ('family-benchmark-results.json','matched-pair-results.json'):
        data=json.loads((args.replay_dir/name).read_text())
        for case in data['cases']:
            p={key:Q.from_float(float(case['parameters'][key])) for key in ('k','w','v','f0','f1','f2')}
            models.append((case['parameters']['name'],p))
    results=[]
    for name,p in models:
        b=bound(p,D,K)
        results.append(dict(name=name,exact_model=exact_text(p),
                            bound=exact_text(b),display_C_H=float(b['C_H']),
                            exact_max_abs_remainder_at_delta_max=exact_text(b['C_H']*D**3)))
    out=dict(status='PASS_EXACT_CONDITIONAL_POLYNOMIAL_REMAINDER',
             source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             assumptions=dict(delta_max=exact_text(D),K_eta=exact_text(K),
                              scalar_shift_bound='An exact branch must independently satisfy |eta_b-a*delta| <= K_eta*delta^2 for every 0<delta<=delta_max.',
                              exact_boundary_identity='H2=delta*A(eta_b)+delta^2*B(eta_b)'),
             numerical_example_status='K_eta is an input assumption. No output here establishes it for the physical boundary-value problem.',
             rationalization='Each recorded binary64 model coefficient is interpreted as its exact rational value; these are explicit nearby polynomial models.',
             proof='Let E=|a|+K_eta*delta_max. Then |eta_b|<=E*delta. Write A=A0+A1*eta+sum(n>=2)A_n*eta^n and B=B0+sum(n>=1)B_n*eta^n. Subtract P2=A0*delta+(A1*a+B0)*delta^2. Bound each term by the triangle inequality to get C_H=|A1|*K_eta+sum(n>=2)|A_n|*E^n*delta_max^(n-2)+sum(n>=1)|B_n|*E^n*delta_max^(n-1).',
             controls=controls([p for _,p in models]),models=results,
             limitations=['No branch existence, local uniqueness, stability, ODE enclosure or novelty is established by this certificate.',
                          'Decimal display_C_H is illustrative; the accompanying rational value is the exact certified conditional constant.'])
    args.out.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(status=out['status'],models=len(results),controls=out['controls'],conditional=True)))


if __name__=='__main__':
    main()
