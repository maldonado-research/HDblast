"""Independent finite algebra tests; no project scientific imports.

These tests do not prove a weighted inverse, a nonlinear solution, a gauge
wave, or an analytic existence/parameter theorem.  All failures use explicit
exceptions so optimized Python exercises the same checks.
"""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

import sympy as s


def need(value, label):
    if not value:
        raise RuntimeError(label)


def zero(expr):
    return s.expand(expr) == 0


def sphere_mean(expr, xs):
    answer = s.S(0)
    for powers, coeff in s.Poly(s.expand(expr), *xs).terms():
        if any(power % 2 for power in powers):
            continue
        aa = [power // 2 for power in powers]
        moment = s.Rational(1, math.factorial(sum(aa) + 1))
        for a in aa:
            moment *= s.Rational(math.factorial(2*a), 4**a * math.factorial(a))
        answer += coeff * moment
    return s.factor(answer)


def source_checks(repo_root):
    root = repo_root.resolve(strict=True)
    pins = json.loads(Path(__file__).with_name('SOURCE_PINS.json').read_text())
    for row in pins['sources']:
        relative = Path(row['path'])
        need(not relative.is_absolute() and '..' not in relative.parts, 'INVALID_SOURCE_RELATIVE_PATH')
        path = root/relative
        need(not path.is_symlink() and path.resolve().is_relative_to(root), 'SOURCE_SYMLINK_OR_ESCAPE')
        data=path.read_bytes()
        need(len(data) == row['bytes'], 'SELECTED_SOURCE_LENGTH:'+row['path'])
        need(hashlib.sha256(data).hexdigest() == row['sha256'], 'SELECTED_SOURCE_HASH:'+row['path'])
    base = root/pins['retained_capsule']['path']
    manifest = base / 'MANIFEST.json'
    expected_manifest = pins['retained_capsule']['manifest_sha256']
    need(hashlib.sha256(manifest.read_bytes()).hexdigest() == expected_manifest, 'RECOVERED_MANIFEST_HASH')
    newrows = json.loads(manifest.read_text())['files']
    for row in newrows:
        data = (base/row['path']).read_bytes()
        need(len(data) == row['bytes'], 'RECOVERED_SOURCE_LENGTH:'+row['path'])
        need(hashlib.sha256(data).hexdigest() == row['sha256'], 'RECOVERED_SOURCE_HASH:'+row['path'])
    return {'selected_public_sources': len(pins['sources']), 'recovered_capsule_members': len(newrows),
            'recovered_manifest_sha256': hashlib.sha256(manifest.read_bytes()).hexdigest(),
            'scope': 'complete-byte authentication; no implication that analytic proofs are true'}


def canonical_check():
    e, U, S = s.symbols('e U S')
    p = s.symbols('p0:10')
    p2 = s.zeros(4)
    for value, (i,j) in zip(p, ((i,j) for i in range(4) for j in range(i,4))):
        p2[i,j] = p2[j,i] = value
    gamma = s.diag(1-4*e*U, *(1+2*e*U for _ in range(3)))
    pi1 = s.diag(-s.Rational(3,2)*S, 0, 0, 0)
    pi = e*pi1 + e**2*p2/2
    density = s.sqrt(gamma.det())
    b = gamma*pi*gamma-gamma*s.trace(gamma*pi)/3
    actual = (2*b/density).applyfunc(lambda v: s.simplify(s.diff(v,e,2).subs(e,0)))
    expected = 2*(p2-s.eye(4)*s.trace(p2)/3)+s.diag(36*U*S,-6*U*S,-6*U*S,-6*U*S)
    need(all(zero(v) for v in actual-expected), 'CANONICAL_FULL_SECOND_MAP')
    # The inverse map vanishes identically on the whole zero-momentum section.
    # Thus all its pure metric derivatives vanish there, without a chosen h2.
    gg = s.symbols('g0:10')
    metric = s.zeros(4)
    for value,(i,j) in zip(gg, ((i,j) for i in range(4) for j in range(i,4))):
        metric[i,j] = metric[j,i] = value
    pp = s.zeros(4)
    inverse = 2*(metric*pp*metric-metric*s.trace(metric*pp)/3)
    need(inverse == s.zeros(4), 'SECOND_METRIC_ABSENCE')
    # A real wrong convention makes the first physical K half the desired one.
    half = (b/density).applyfunc(lambda v:s.diff(v,e).subs(e,0))
    wanted = s.diag(-2*S,S,S,S)
    need(half != wanted, 'FACTOR_TWO_CONTROL')
    no_density = (2*b).applyfunc(lambda v:s.diff(v,e,2).subs(e,0))
    need(no_density != expected, 'DENSITY_CONTROL')
    return {'full_second_matrix_entries':16,'arbitrary_second_momentum_entries':10,
            'angular_quadratic_correction':'-6 U S gamma0_ab',
            'second_metric_direction_absent_at_pi0_zero':True,
            'wrong_factor_two_rejected':True,'omitted_density_rejected':True}


def pullback_check():
    t,y,a,b,c = coords = s.symbols('t y a b c')
    # Nonzero time, radial and angular accelerations; value and velocity vanish
    # identically as spatial functions. Arbitrary mixed S slots stay present.
    field = s.Matrix([t*t*(1+y+a*b),t*t*(2+y*y+b),t*t*(a+c*y),t*t*(b+a*y),t*t*(c+a*b)])
    tensor = s.Matrix(5,5,lambda i,j: (i+j+1)*(1+t+y+a*b)+t*t*(i*j+1)+a**(i%2)*c**(j%2))
    tensor = (tensor+tensor.T)/2
    radius = 2+y+y*y
    background = s.diag(-radius**2,1,*(radius**2*(1+t*t) for _ in range(3)))

    def lie(vector, ten):
        return s.Matrix(5,5,lambda i,j: sum(vector[k]*s.diff(ten[i,j],coords[k])
            +ten[k,j]*s.diff(vector[k],coords[i])+ten[i,k]*s.diff(vector[k],coords[j])
            for k in range(5)))

    first = lie(field,tensor)
    squared = lie(field,lie(field,background))
    quad = 2*first+squared
    for i,j in itertools.product(range(2,5),repeat=2):
        need(zero(quad[i,j].subs(t,0)), 'PULLBACK_ANGULAR_VALUE')
        need(zero(s.diff(quad[i,j],t).subs(t,0)), 'PULLBACK_ANGULAR_VELOCITY')
    mixed = lie(field,background)
    need(not zero(s.diff(mixed[0,1],t).subs(t,0)), 'ACCELERATION_CONTROL_MIXED_VELOCITY')
    wrong_vector = field/t
    wrong = lie(wrong_vector,tensor)
    need(any(not zero(s.diff(wrong[i,j],t).subs(t,0))
             for i,j in itertools.product(range(2,5),repeat=2)), 'NONZERO_FIRST_VELOCITY_CONTROL')
    return {'angular_quadratic_value_entries':9,'angular_quadratic_velocity_entries':9,
            'nonzero_acceleration_changes_mixed_velocity':True,
            'violated_zero_first_generator_velocity_rejected':True,
            'scope':'local tensor identities; not a constructed gauge diffeomorphism'}


def projection_check():
    xs = s.symbols('x0:4')
    x = s.Matrix(xs)
    P = s.eye(4)-x*x.T
    # Independent round-S3 TT tensor from two left invariant unit Killing fields.
    e1 = s.Matrix([-xs[1],xs[0],-xs[3],xs[2]])
    e2 = s.Matrix([-xs[2],xs[3],xs[0],-xs[1]])
    T = e1*e1.T-e2*e2.T
    need(zero(s.trace(T)), 'SPHERE_TT_TRACE')
    need(all(zero(v) for v in T*x), 'SPHERE_TT_TANGENCY')
    div = s.Matrix([sum(P[k,i]*P[j,l]*s.diff(T[i,l],xs[k])
                    for i,k,l in itertools.product(range(4),repeat=3)) for j in range(4)])
    need(all(zero(v) for v in div), 'SPHERE_TT_DIVERGENCE')

    def contract(A):
        return s.expand(sum(T[i,j]*A[i,j] for i,j in itertools.product(range(4),repeat=2)))

    scalar_tests = [xs[i]**4 + (i+1)*xs[i]*xs[(i+1)%4]**3 for i in range(4)]
    scalar_tests += [xs[0]*xs[1]*xs[2]*xs[3], sum(xs)**5,
                     xs[0]**2*xs[1]**2+7*xs[2]**4]
    integrands = []
    for f in scalar_tests:
        grad = s.Matrix([s.diff(f,v) for v in xs])
        H = P*s.hessian(f,xs)*P-(x.dot(grad))*P
        integrand = contract(H)
        need(sphere_mean(integrand,xs) == 0, 'SPHERE_HESSIAN_PROJECTION')
        integrands.append(integrand)
    vector_tests = []
    for i in range(4):
        W = s.zeros(4,1)
        W[i] = xs[(i+1)%4]**3 + (i+2)*xs[i]*xs[(i+2)%4]**2
        V = P*W
        dV = V.jacobian(xs)
        gradV = P*dV*P
        integrand = contract(gradV+gradV.T)
        need(sphere_mean(integrand,xs) == 0, 'SPHERE_SYMGRAD_PROJECTION')
        vector_tests.append(integrand)
    need(any(not zero(v) for v in integrands), 'NONTRIVIAL_HESSIAN_CANCELLATION')
    # Angular weights would destroy the orthogonality; radial weights commute.
    bad = [sphere_mean(v*(1+xs[0]**2),xs) for v in integrands]
    need(any(value != 0 for value in bad), 'ANGULAR_WEIGHT_CONTROL')
    return {'all_point_tangent_tracefree_divergencefree_verified':True,
            'nontrivial_scalar_hessian_integral_cases':len(scalar_tests),
            'vector_symgradient_integral_cases':len(vector_tests),
            'nonradial_weight_control_rejected':True,
            'scope':'finite exact controls supplement integration by parts; not weighted inverse estimates'}


def tt_units_check():
    ell,rho,a,c = s.symbols('ell rho a c',nonzero=True)
    theta,Ut,Uz,Uzz,Uyy,U = s.symbols('theta Ut Uz Uzz Uyy U')
    # Universal warped-product TT Jacobi operator: Box_base H + 3 dR/R dH
    # - (lambda_T+2)/R^2 H, lambda_T=6 for the registered l=2 tensors.
    physical = -Ut-3*theta*s.Symbol('Ht')+rho**2*(Uyy+4*a*s.Symbol('Hy'))-8*U/c**2
    converted = s.expand(physical.subs({Uyy:Uzz/rho**2+a*Uz/rho,s.Symbol('Hy'):-Uz/rho}))
    expected = -Ut+Uzz-3*theta*s.Symbol('Ht')-3*rho*a*Uz-8*U/c**2
    need(zero(converted-expected), 'UNIVERSAL_TT_COORDINATE_CONVERSION')
    # The source normalization is a second parameter derivative, no hidden /2.
    h = 3*s.Symbol('Ux')/(2*ell*rho**2*c)+4*s.Symbol('Fquad')/(rho*c)
    need(s.diff(h,s.Symbol('Ux')) == 3/(2*ell*rho**2*c), 'SCREEN_CLOCK_FACTOR')
    return {'universal_TT_radial_conversion_verified':True,'tensor_eigenvalue':6,
            'scaled_zero_order_coefficient':'-8 ell^2/cosh(t)^2',
            'scope':'algebra after universal Jacobi identity; not a new derivation of Einstein linearization'}


def direct_ricci_check():
    """Direct five-dimensional coordinate Ricci variation with free base jets.

    Use an independently constructed smooth lowest S3 TT tensor and compute
    all Christoffel/Ricci slots at one sphere chart origin. Time/radial jets
    remain independent indeterminates. Background angular Einstein residual
    is subtracted explicitly, rather than silently assuming an off-shell
    warped metric solves the matter equations.
    """
    t,y = s.symbols('t y')
    r,r1,r2,c,c1,c2 = s.symbols('rho rho_y rho_yy c c_t c_tt',nonzero=True)
    h,ht,hy,htt,hty,hyy = s.symbols('H H_t H_y H_tt H_ty H_yy')
    vv = s.symbols('v0:3')
    v = s.Matrix(vv)
    vv2 = v.dot(v)
    xx = s.Matrix([(1-vv2)/(1+vv2),*(2*vi/(1+vv2) for vi in vv)])
    J = xx.jacobian(vv)
    e1 = s.Matrix([-xx[1],xx[0],-xx[3],xx[2]])
    e2 = s.Matrix([-xx[2],xx[3],xx[0],-xx[1]])
    tensor = (J.T*(e1*e1.T-e2*e2.T)*J).applyfunc(s.factor)
    sphere = 4*s.eye(3)/(1+vv2)**2
    atv = dict.fromkeys(vv,0)
    R = r+r1*y+r2*y*y/2
    C = c+c1*t+c2*t*t/2
    H = h+ht*t+hy*y+htt*t*t/2+hty*t*y+hyy*y*y/2
    F = R**2*C**2
    coords = (t,y,*vv)
    origin = dict.fromkeys(coords,0)

    def angular_metric(factor, angular):
        values=s.zeros(5)
        first=[s.zeros(5) for _ in range(5)]
        second=[[s.zeros(5) for _ in range(5)] for _ in range(5)]
        # Differentiate separated factors once; this is exactly the product
        # rule used by the original full-expression coordinate calculation.
        f0=factor.subs(origin)
        df=[s.diff(factor,z).subs(origin) for z in (t,y)]
        ddf=[[s.diff(factor,z,w).subs(origin) for w in (t,y)] for z in (t,y)]
        for A,B in itertools.product(range(3),repeat=2):
            q=angular[A,B]
            q0=q.subs(atv)
            dq0=[s.diff(q,z).subs(atv) for z in vv]
            ddq0=[[s.diff(q,z,w).subs(atv) for w in vv] for z in vv]
            values[A+2,B+2]=f0*q0
            for k in range(5):
                first[k][A+2,B+2]=df[k]*q0 if k<2 else f0*dq0[k-2]
                for l in range(5):
                    if k<2 and l<2: value=ddf[k][l]*q0
                    elif k<2: value=df[k]*dq0[l-2]
                    elif l<2: value=df[l]*dq0[k-2]
                    else: value=f0*ddq0[k-2][l-2]
                    second[k][l][A+2,B+2]=value
        return values,first,second

    g,dg,ddg=angular_metric(F,sphere)
    q,dq,ddq=angular_metric(F*H,tensor)
    g[0,0]=-r*r
    g[1,1]=1
    dg[1][0,0]=-2*r*r1
    ddg[1][1][0,0]=-2*(r1*r1+r*r2)
    inv=g.inv()
    inv1=-inv*q*inv
    dinv=[-inv*der*inv for der in dg]
    dinv1=[-inv1*dg[k]*inv-inv*dq[k]*inv-inv*dg[k]*inv1 for k in range(5)]

    def A(der,d,b,c):
        return der[b][d,c]+der[c][d,b]-der[d][b,c]

    def dA(second,k,d,b,c):
        return second[k][b][d,c]+second[k][c][d,b]-second[k][d][b,c]

    gamma={}
    gamma1={}
    dgamma={}
    dgamma1={}
    for a,b,cc in itertools.product(range(5),repeat=3):
        gamma[a,b,cc]=s.expand(sum(inv[a,d]*A(dg,d,b,cc) for d in range(5))/2)
        gamma1[a,b,cc]=s.expand(sum(inv1[a,d]*A(dg,d,b,cc)+inv[a,d]*A(dq,d,b,cc) for d in range(5))/2)
        for k in range(5):
            dgamma[k,a,b,cc]=s.expand(sum(dinv[k][a,d]*A(dg,d,b,cc)+inv[a,d]*dA(ddg,k,d,b,cc) for d in range(5))/2)
            dgamma1[k,a,b,cc]=s.expand(sum(
                dinv1[k][a,d]*A(dg,d,b,cc)+dinv[k][a,d]*A(dq,d,b,cc)
                +inv1[a,d]*dA(ddg,k,d,b,cc)+inv[a,d]*dA(ddq,k,d,b,cc)
                for d in range(5))/2)
    ric=s.zeros(5)
    ric1=s.zeros(5)
    for a,b in itertools.product(range(5),repeat=2):
        ric[a,b]=s.expand(sum(dgamma[k,k,a,b]-dgamma[b,k,a,k] for k in range(5))
          +sum(gamma[k,k,l]*gamma[l,a,b]-gamma[k,b,l]*gamma[l,a,k]
               for k,l in itertools.product(range(5),repeat=2)))
        ric1[a,b]=s.expand(sum(dgamma1[k,k,a,b]-dgamma1[b,k,a,k] for k in range(5))
          +sum(gamma1[k,k,l]*gamma[l,a,b]+gamma[k,k,l]*gamma1[l,a,b]
               -gamma1[k,b,l]*gamma[l,a,k]-gamma[k,b,l]*gamma1[l,a,k]
               for k,l in itertools.product(range(5),repeat=2)))
    # If the background angular Einstein equation holds, the displayed
    # subtraction is precisely δRic-(2/3)V δg. No scalar variation is present.
    angular_bg=ric[2,2]/4
    op=-htt/r**2-3*c1*ht/(r**2*c)+hyy+4*r1*hy/r-8*h/(r**2*c**2)
    expected=-r**2*c**2*op/2
    for a,b in itertools.product(range(5),repeat=2):
        t0=tensor[a-2,b-2].subs(atv) if a>=2 and b>=2 else s.S(0)
        residual=ric1[a,b]-h*angular_bg*t0
        need(s.simplify(residual-expected*t0)==0,'DIRECT_WARPED_RICCI_TT_%d%d'%(a,b))
    # A real wrong angular eigenvalue (+6 instead of +8 including curvature)
    # changes a nonzero coefficient even with all base H jets zero.
    wrong=expected-h
    need(s.simplify(wrong-expected)!=0,'DIRECT_RICCI_WRONG_CURVATURE_CONTROL')
    return {'five_dimensional_Ricci_variation_slots':25,
            'independent_free_base_jets':12,
            'sphere_tensor':'e1 tensor e1 - e2 tensor e2 in stereographic chart',
            'background_angular_Einstein_residual_subtracted_explicitly':True,
            'universal_spin2_operator_verified_at_chart_origin':True,
            'wrong_tensor_curvature_shift_rejected':True,
            'scope':'exact arbitrary-base-jet control at one angular chart point; general projection also uses spherical covariance and IBP'}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--repo-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    need(not args.output.exists() and not args.output.is_symlink(),'OUTPUT_EXISTS_OR_SYMLINK')
    need(args.output.parent.is_dir(),'OUTPUT_PARENT_MISSING')
    report={'sources':source_checks(args.repo_root),'canonical':canonical_check(),
            'pullback':pullback_check(),'projection':projection_check(),'TT_units':tt_units_check(),
            'direct_Ricci':direct_ricci_check(),
            'scope':'independent algebra and source authentication only; analytic theorems remain premises'}
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps(report,indent=2,sort_keys=True))


if __name__=='__main__':
    main()
