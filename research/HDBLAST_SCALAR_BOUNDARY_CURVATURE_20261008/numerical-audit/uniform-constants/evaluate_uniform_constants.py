"""Evaluate sufficient inequalities of the uniform branch proof candidate.

All calculations are exact rational arithmetic. Transcendental constants are
replaced by proved rational outer bounds. Passing these inequalities is
conditional on the independently reviewed validity of the cited analytic proof;
this program alone does not prove its operator estimates.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path


def add(a,b):
    return [(a[i] if i<len(a) else Q(0))+(b[i] if i<len(b) else Q(0)) for i in range(max(len(a),len(b)))]


def scale(a,c): return [x*c for x in a]


def mul(a,b):
    out=[Q(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): out[i+j]+=x*y
    return out


def deriv(a): return [i*a[i] for i in range(1,len(a))] or [Q(0)]


def sup_abs(a,r): return sum(abs(x)*r**i for i,x in enumerate(a))


def pack(x):
    if type(x) is Q: return str(x.numerator)+'/'+str(x.denominator)
    if isinstance(x,dict): return {k:pack(v) for k,v in x.items()}
    if isinstance(x,list): return [pack(v) for v in x]
    return x


def constants(p,r,b):
    k,w,v,f0,f1,f2=(p[x] for x in ('k','w','v','f0','f1','f2'))
    if not (k>0 and w>4*k and f0>0 and r>0 and b>0):
        raise ValueError('Invalid theorem domain')
    W,f=[3*k,Q(0),w/2,v/6],[f0,f1,f2/2]
    Wp,fp=deriv(W),deriv(f)
    U=add(scale(mul(Wp,Wp),Q(1,2)),scale(mul(W,W),-Q(2,3)))
    A=add(scale(mul(W,f),Q(1,9)),scale(mul(Wp,fp),-Q(1,12)))
    BB=add(scale(mul(f,f),Q(1,36)),scale(mul(fp,fp),-Q(1,48)))
    m2=w*(w-4*k); lam=w-4*k; nu=lam/2
    mu=nu*(nu+4*k); gap=m2-mu
    zeta_lower=mu*Q(15,16)/(Q(16,3)*k+nu)
    J_upper=2/k+1/(2*zeta_lower)
    C_G=max(1/gap,(1+m2/gap)/(4*k*k))
    M=max(Q(1),lam/k); ball=2*M
    M2=sup_abs(deriv(deriv(U)),r); M3=sup_abs(deriv(deriv(deriv(U))),r)
    L1=M3/2; C_F=k*k/4+M2/12
    epsilon=C_F*ball**2*b*b*J_upper/(2*k)
    C_N=L1*ball**2+4*C_F*ball**3*b
    L_N=2*L1*ball*b+20*C_F*ball**2*b*b
    C_P=k*C_G*ball*(2*L1*ball+20*C_F*ball**2*b)
    C_D=k*C_G*C_N; C_a=4*C_F*ball**2/k
    C_L_upper=8*lam+Q(28,3)*k*lam/(w-2*k)
    C_H_DN=4*C_L_upper/(9*k*k)
    W2=sup_abs(deriv(Wp),r); W3=sup_abs(deriv(deriv(Wp)),r)
    F1=sup_abs(fp,r); F2=sup_abs(deriv(fp),r)
    F_abs=sup_abs(f,r); W_abs=sup_abs(W,r)
    L=(abs(f1)+1)/w; a1=k*f0/3
    C_K=L*sup_abs(deriv(A),r)+sup_abs(BB,r)
    a_bar=Q(4,3)*k+C_F*ball**2*b*b/k
    C_Tp=M2+4*a_bar*k*ball+k*ball*(lam+C_P*b)
    C_Tbeta=2*C_Tp/w
    C_slope=k*ball*C_a*L**2+((C_a+W2/3)*L+F1/6)*C_Tbeta*L
    d0=2*(w-2*k); alpha=-f1/(2*d0)
    K_eta=(C_H_DN*(3*a1/2)*L+(C_D+W3/2)*L**2+F2*L/2)/d0
    K_H=abs(A[1])*K_eta+sup_abs(deriv(deriv(A)),r)*L**2/2+sup_abs(deriv(BB),r)*L
    return locals()


def stage_one(c):
    return [
        ('field_tube_inside_polynomial_radius',c['ball']*c['b'],c['r'],False),
        ('metric_volterra_contraction',c['epsilon'],Q(1,4),False),
        ('scalar_contraction',c['C_G']*c['L_N'],Q(1,2),False),
        ('scalar_ball_preservation',c['C_G']*c['C_N']*c['b'],c['M'],False),
        ('positive_metric_derivative',c['C_F']*c['ball']**2*c['b']**2/c['k'],c['k']/2,False),
        ('scalar_junction_derivative_feasibility',(c['C_P']+c['W3'])*c['b'],c['w']/2,True),
    ]


def stage_two(c,d):
    return [
        ('scalar_root_inside_boundary_tube',c['L']*d,c['b'],False),
        ('scalar_junction_derivative',(c['C_P']+c['W3'])*c['b']+d*c['F2']/2,c['w']/2,False),
        ('positive_reduced_curvature',c['C_K']*d,c['a1']/2,False),
        ('metric_residual_positive_at_Y',3*c['a1']*d/2,Q(9,32)*c['k']**2,True),
        ('positive_shell_tension',c['W2']*c['b']**2/6+d*c['F_abs']/6,c['k']/2,False),
        ('strict_downward_shell_crossing',c['C_slope']*d,c['a1']/4,False),
    ]


def check(rows): return all(lhs<rhs if strict else lhs<=rhs for _,lhs,rhs,strict in rows)


def represent_checks(rows):
    return [dict(name=n,lhs=pack(lhs),rhs=pack(rhs),strict=strict,pass_check=lhs<rhs if strict else lhs<=rhs,
                 exact_margin=pack(rhs-lhs)) for n,lhs,rhs,strict in rows]


def evaluate(name,p,r=Q(1,100)):
    # A deterministic conservative dyadic search, not an optimality claim.
    M=max(Q(1),(p['w']-4*p['k'])/p['k'])
    b=r/(2*M)
    for n in range(512):
        c=constants(p,r,b)
        if check(stage_one(c)): break
        b/=2
    else: raise ArithmeticError('No first-stage constants found')
    d=min(Q(1,1000),b/c['L'])
    for j in range(512):
        if check(stage_two(c,d)): break
        d/=2
    else: raise ArithmeticError('No second-stage constants found')
    E=abs(c['alpha'])+c['K_eta']*d
    K_poly=abs(c['A'][1])*c['K_eta']
    K_poly+=sum(abs(c['A'][i])*E**i*d**(i-2) for i in range(2,len(c['A'])))
    K_poly+=sum(abs(c['BB'][i])*E**i*d**(i-1) for i in range(1,len(c['BB'])))
    best_K_H=min(c['K_H'],K_poly)
    q=p['k']*p['f1']**2/(24*(p['w']-2*p['k']))
    suppression_margin=q-best_K_H*d
    if not (q>0 and suppression_margin>0):
        raise ArithmeticError('Selected examples do not establish strict curvature suppression')
    if check(stage_two(c,Q(1))): raise ArithmeticError('Large-detuning negative control did not fail')
    bad=constants(p,r,r)
    if check(stage_one(bad)): raise ArithmeticError('Oversized tube negative control did not fail')
    if not check(stage_one(c)+stage_two(c,d)): raise ArithmeticError('Certificate check failed')
    keep=('m2','lam','nu','mu','gap','zeta_lower','J_upper','C_G','M','ball','M2','M3','L1','C_F','epsilon',
          'C_N','L_N','C_P','C_D','C_a','C_L_upper','C_H_DN','W2','W3','F1','F2','F_abs','W_abs','L','a1','C_K',
          'a_bar','C_Tp','C_Tbeta','C_slope','d0','alpha','K_eta','K_H')
    return dict(name=name,parameters=pack(p),r=pack(r),b=pack(b),delta_star=pack(d),
                exact_constants=pack({x:c[x] for x in keep}),polynomial_coefficients=pack({x:c[x] for x in ('W','f','U','A','BB')}),
                checks=represent_checks(stage_one(c)+stage_two(c,d)),
                K_H_polynomial_refinement=pack(K_poly),
                exact_bound_at_delta_star=dict(scalar=pack(c['K_eta']*d*d),curvature=pack(best_K_H*d**3)),
                strict_curvature_suppression=dict(q=pack(q),best_K_H=pack(best_K_H),
                    lower_magnitude_coefficient=pack(suppression_margin),upper_magnitude_coefficient=pack(q+best_K_H*d),
                    exact_positive_margin=pack(suppression_margin),pass_check=q>0 and suppression_margin>0,
                    inequality='-(q+K*delta_star)*delta^2 <= H2-H2_metric <= -(q-K*delta_star)*delta^2 < 0 for every 0<delta<=delta_star.',
                    comparator='H2_metric=k*f0*delta/3+f0^2*delta^2/36 is the exact constant-scalar branch of the separate constant-tension model f(eta)=f0 with the same bulk W. For f1!=0 it generally violates the scalar junction of the scalar-dependent tension model.'),
                display_only=dict(b=float(b),delta_star=float(d),K_eta=float(c['K_eta']),K_H_taylor=float(c['K_H']),
                                  K_H_polynomial=float(K_poly),suppression_margin=float(suppression_margin),
                                  registered_sample_min_delta_over_delta_star=float(Q(1,2000)/d) if name.startswith('registered_') else None),
                search_halvings=dict(first_stage=n,second_stage=j),
                negative_controls=dict(delta_one_rejected=True,oversized_tube_rejected=True),
                uploaded_detuning_range_covered=d>=Q(1,500) if name.startswith('registered_') else None)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--theorem',type=Path,required=True)
    ap.add_argument('--theorem-sha256',required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    proof_bytes=args.theorem.read_bytes()
    proof_sha=hashlib.sha256(proof_bytes).hexdigest()
    if proof_sha != args.theorem_sha256: raise ValueError('Proof candidate pin mismatch')
    c=2/Q('1.0357712571566784')-Q(4,3)
    models=[('registered_exact_decimal_model',dict(k=Q(1,9),w=Q(2),v=Q(2),f0=1+c,f1=c,f2=Q(0))),
            ('simple_polynomial_model',dict(k=Q(1),w=Q(5),v=Q(0),f0=Q(1),f1=Q(1),f2=Q(0)))]
    results=[evaluate(n,p) for n,p in models]
    out=dict(status='PASS_EXACT_SUFFICIENT_INEQUALITIES_CONDITIONAL_ON_ANALYTIC_PROOF',
             source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             theorem_sha256=proof_sha,theorem_review_status='PENDING_INDEPENDENT_ANALYTIC_REVIEW_AT_COMPUTATION',
             arithmetic='All inputs, polynomial derivative bounds, constants, inequalities and output bounds use exact fractions.Fraction. Display fields are floating-point approximations only.',
             transcendental_bounds=['7<exp(2)<8: its degree-4 Taylor partial sum is 7; the remaining positive series is at most (4/15)/(1-1/3)=2/5, so exp(2)<7.4<8.',
                 'coth(1)=(exp(2)+1)/(exp(2)-1)<4/3.',
                 '1/(1-exp(-2))<7/6.',
                 'C0/k=4*coth(1)+nu/k>4, and exp(4)>16, hence 1-exp(-C0/k)>15/16.',
                 'C0<=16*k/3+nu, so zeta>=mu*(15/16)/(16*k/3+nu).',
                 'sinh(1)^2=(exp(2)+exp(-2)-2)/4<7/4<2, hence H0(1/k)^2>k^2/2.'],
             derivative_bound='For polynomial P(x)=sum p_i*x^i on |x|<=r, sup|P|<=sum|p_i|*r^i, applied after exact polynomial differentiation.',
             models=results,
             model_scope='The registered model here interprets the declared decimal 1.0357712571566784 and k=1/9 as exact rationals. It is a specified polynomial model, not a claim that binary64 shooting has zero parameter-rounding error.',
             limitations=['The program checks arithmetic sufficient inequalities, not the analytic operator estimates supplying them.',
                          'The registered-model range is smaller than every uploaded detuning. The simple model was not in the uploaded sweep.',
                          'No dynamical stability, physical observation or external novelty follows.'])
    args.out.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(status=out['status'],models=[dict(name=x['name'],**x['display_only']) for x in results]),indent=2))


if __name__=='__main__': main()
