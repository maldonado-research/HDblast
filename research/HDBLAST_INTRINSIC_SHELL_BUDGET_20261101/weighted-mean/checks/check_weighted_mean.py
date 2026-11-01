#!/usr/bin/env python3
"""Exact monopole algebra; authenticates repository sources before SymPy import."""
from pathlib import Path
import argparse
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]


def authenticate(repo_root, pins):
    for entry in pins['source_files']:
        path = (repo_root / entry['path']).resolve()
        if not path.is_relative_to(repo_root):
            raise ValueError('Source path escapes repository: ' + entry['path'])
        raw = path.read_bytes()
        if hashlib.sha256(raw).hexdigest() != entry['sha256'] or len(raw) != entry['bytes']:
            raise ValueError('Changed source: ' + entry['path'])


def evaluate(repo_root, pins):
    import sympy as S
    y = S.symbols('y', real=True)
    rho, a, p, u, v, s, j, f, xi, N, A, B, W = [S.Function(k)(y) for k in
        ('rho', 'a', 'p', 'u', 'v', 's', 'j', 'f', 'xi', 'N', 'A', 'B', 'W')]
    Vp, Vpp = S.symbols('Vp Vpp')
    V = p**2 / 2 - 6*a**2 + 6/rho**2
    D = lambda z: S.diff(z, y)
    rules = {D(rho): a*rho, D(a): -1/rho**2-p**2/3,
             D(p): Vp-4*a*p}

    def reduce(z):
        z = S.expand(z)
        for unused in range(4):
            z = S.expand(z.doit().subs(rules, simultaneous=True))
        return S.factor(z)

    checks = []
    controls = []

    def check(name, difference):
        residual = reduce(difference)
        if residual != 0:
            raise ValueError(f'{name}: {residual}')
        checks.append(name)

    def reject(name, difference):
        residual = reduce(difference)
        if residual == 0:
            raise ValueError(f'Wrong-formula control was not rejected: {name}')
        controls.append(name)

    def Pg(Ar, Bt):
        return (-3*D(D(Bt))/2 + 3*a*D(Ar)/2 - 6*a*D(Bt)
                + (3*D(a)+6*a*a+p*p/2)*Ar - 3*Bt/rho**2)

    def Jr(n):
        return (-3*a*D(n)+(3*D(a)+3*a*a+p*p)*n)/2

    def Jt(n):
        return (-D(D(n))-2*a*D(n)+(D(a)+3*a*a-2/rho**2)*n)/2

    def mom(Tr, Ut):
        return -2*(D(Tr)+3*a*(Tr-Ut))

    def run():
        # Linearized scalar curvature, independently built from divergence and Ricci.
        divh = D(A)+3*a*(A-B)
        divdiv = D(divh)+3*a*divh
        laptrace = D(D(A+3*B))+3*a*D(A+3*B)
        rich = -3*(D(a)+a*a)*A+3*(2/rho**2-D(a)-3*a*a)*B
        check('radial_metric_operator', (divdiv-laptrace-rich+p*p*A)/2-Pg(A,B))
        # Formal adjoints via Euler-Lagrange coefficients in the radial volume.
        lag = rho**3*N*Pg(A,B)
        adjA = S.diff(lag,A)-D(S.diff(lag,D(A)))
        adjB = S.diff(lag,B)-D(S.diff(lag,D(B)))+D(D(S.diff(lag,D(D(B)))))
        check('radial_adjoint', adjA/rho**3-Jr(N))
        check('tangential_adjoint_multiplicity', adjB/(3*rho**3)-Jt(N))
        check('radial_adjoint_background_simplification', Jr(N)-(-3*a*D(N)+3*(a*a-rho**-2)*N)/2)
        check('momentum_weighted_operator', mom(2*W*D(xi),2*W*a*xi)
              +4*(D(W*D(xi))+3*a*W*(D(xi)-a*xi)))
        # Weighted weak momentum Euler-Lagrange derivative.
        lagM = 2*rho**3*W*(D(xi)**2+3*a*a*xi*xi)
        eulerM = S.diff(lagM,xi)-D(S.diff(lagM,D(xi)))
        check('momentum_weak_form', eulerM/rho**3-mom(2*W*D(xi),2*W*a*xi))

        # Canonical affine Hessian, reconstructed from exact gamma/pi/Pi maps.
        e, Y = S.symbols('epsilon Y')
        gm = S.diag(1-4*e*u*Y,1+2*e*u*Y,1+2*e*u*Y,1+2*e*u*Y)
        dens = S.sqrt(gm.det())
        Tr, Ut = S.symbols('Tr Ut')
        tau = S.diag(-3*e*s*Y/2+e*e*Tr/2,e*e*Ut/2,e*e*Ut/2,e*e*Ut/2)
        K = 2*(gm*tau*gm-gm*S.trace(gm*tau)/3)/dens
        K2 = K.diff(e,2).subs(e,0)
        Pi = e*j*Y/dens
        Pi2 = S.diff(Pi,e,2).subs(e,0)
        check('canonical_radial_K2', K2[0,0]-(4*Tr/3-2*Ut+36*u*s*Y*Y))
        check('canonical_tangential_K2', K2[1,1]-(-2*Tr/3-6*u*s*Y*Y))
        check('physical_Pi2', Pi2+2*u*j*Y*Y)
        scalar_t = rho*(1+e*u*Y)*Pi
        check('coordinate_scalar_second_velocity_zero', S.diff(scalar_t,e,2).subs(e,0))
        vel = 2*rho*(1+e*u*Y)*K
        vel2 = vel.diff(e,2).subs(e,0)
        mean = lambda expr: S.expand(expr).subs(Y*Y,S.Rational(1,6))
        check('mean_metric_velocity_radial', mean(vel2[0,0])-2*rho*(4*Tr/3-2*Ut+16*u*s/3))
        check('mean_metric_velocity_tangent', mean(vel2[1,1])-2*rho*(-2*Tr/3-2*u*s/3))
        check('mean_scalar_normal_velocity', mean(Pi2)+u*j/3)
        gtt2 = S.diff(-rho*rho*(1+e*u*Y)**2,e,2).subs(e,0)
        check('mean_initial_lapse_metric', mean(gtt2)+rho*rho*u*u/3)

        # Direct density-momentum products are independent of physical K Hessians.
        radial_direct = -24*u*(D(s)+3*a*s)-24*s*D(u)+12*s*D(u)+2*j*D(v)
        QM = -4*u*D(s)-12*a*u*s-2*s*D(u)+j*D(v)/3
        check('mean_canonical_momentum_source', radial_direct/6-QM)
        check('momentum_source_first_constraint_reduction',
              (QM-(-2*s*D(u)+j*(D(v)+4*p*u)/3)).subs(D(s),-p*j/3-3*a*s))
        # A second derivation converts physical source S3 to canonical variables.
        qphys = 12*u*D(s)-36*a*u*s-6*s*D(u)-2*j*D(v)
        Kindr, Kindt = 36*u*s, -6*u*s
        induced = -3*D(Kindt)+3*a*(Kindr-Kindt)+2*p*u*j
        check('physical_to_canonical_radial_source',
              (-qphys-induced-radial_direct).subs(D(s),-p*j/3-3*a*s))
        check('physical_to_canonical_angular_source', -(18*u*s-v*j-24*u*s)-(6*u*s+j*v))

        AH = (-24*D(u)**2-288*a*u*D(u)-24*u*D(D(u))
              -(192*D(a)+384*a*a+32*p*p)*u*u
              -336*u*u/rho**2-12*s*s-2*j*j
              -2*D(v)**2-16*p*u*D(v)-2*Vpp*v*v)
        DH = (36*u*u-2*v*v)/rho**2
        QH = (-2*D(u)**2-24*a*u*D(u)-2*u*D(D(u))
              -(16*D(a)+32*a*a+S.Rational(8,3)*p*p)*u*u-4*u*u/rho**2
              -s*s-(j*j+D(v)**2+Vpp*v*v)/6
              -S.Rational(4,3)*p*u*D(v)-S.Rational(4,3)*v*v/rho**2)
        check('Hamiltonian_normalization_and_sphere_average', AH/12+2*DH/3-QH)
        # Direct kinetic expansion checks the half-constraint normalization.
        kin = 2*(S.trace(gm*tau)**2/3-S.trace(gm*tau*gm*tau))/gm.det()-Pi**2/2
        k2 = S.diff(kin,e,2).subs(e,0)
        check('direct_canonical_kinetic_Hessian_mean', mean(k2)+s*s+j*j/6)

        # Independent full Hessian from the radial lapse / angular conformal geometry.
        GY = S.symbols('gradient_Y_squared')
        L = S.sqrt(1-4*e*u*Y)
        blogY = e*u/(1+2*e*u*Y)
        lapblog = -8*Y*blogY+GY*S.diff(blogY,Y)
        lapL = -8*Y*S.diff(L,Y)+GY*S.diff(L,Y,2)
        radial_k = (a+e*D(u)*Y/(1+2*e*u*Y))/L
        sphere_inverse = 1/(rho*rho*(1+2*e*u*Y))
        R4 = (sphere_inverse*(6-4*lapblog-2*blogY*blogY*GY)
              -12*radial_k*radial_k-6*D(radial_k)/L
              -2*sphere_inverse*(lapL+blogY*S.diff(L,Y)*GY)/L)
        R2 = S.diff(R4,e,2).subs(e,0)
        R2target = ((-24*D(u)**2-288*a*u*D(u)-24*u*D(D(u))
                     -(192*D(a)+384*a*a)*u*u-336*u*u/rho**2)*Y*Y
                    +36*u*u*GY/rho**2)
        check('independent_full_curvature_Hessian', R2-R2target)
        gradient = (p+e*D(v)*Y)**2/(1-4*e*u*Y)+e*e*v*v*GY*sphere_inverse
        H2direct = R2/2+k2-S.diff(gradient,e,2).subs(e,0)/2-Vpp*v*v*Y*Y
        check('independent_full_Hamiltonian_Hessian_mean',
              mean(H2direct).subs(GY,S.Rational(4,3))-QH)
        n04, n02n12, n02 = S.Rational(1,8), S.Rational(1,24), S.Rational(1,4)
        check('sphere_quadrupole_square_moment',2*n04-2*n02n12-S.Rational(1,6))
        check('sphere_gradient_square_moment',4*(2*n02-(2*n04-2*n02n12))-S.Rational(4,3))

        # Geometry and source-only moment identities.
        check('areal_radius_adjoint_radial', Jr(rho)+3/(2*rho))
        check('areal_radius_adjoint_tangential', Jt(rho)+1/rho)
        check('conformal_metric_operator', Pg(f,f)+3*(D(D(f))+3*a*D(f))/2+V*f)
        gaugeP = Pg(2*D(xi),2*a*xi)
        check('scalar_fixed_radial_gauge_operator', gaugeP-D(rho**4*p*p*xi)/rho**4)
        # Twice integrating C[f] with weight rho^4 gives this Euler coefficient.
        Clag = rho**4*(-3*(D(D(f))+3*a*D(f))/2-V*f)
        coefficient = S.diff(Clag,f)-D(S.diff(Clag,D(f)))+D(D(S.diff(Clag,D(D(f)))))
        check('compact_kernel_integral_condition', coefficient+S.Rational(9,2)*rho*rho)
        # After its defining first-order equation, the constructed kernel vanishes.
        q = S.Function('q')(y)
        replacement = {D(q): -rho**4*Pg(f,f)}
        xi_kernel = q/(rho**4*p*p)
        check('constructed_kernel_equation',
              (Pg(f,f)+D(rho**4*p*p*xi_kernel)/rho**4).subs(replacement))

        square_form = ((4/rho**2+S.Rational(8,3)*p*p)*u*u-s*s-j*j/6
                       -(D(v)+4*p*u)**2/6-(Vpp/6+4/(3*rho**2))*v*v)
        flux = -2*rho**4*u*D(u)-8*rho**4*a*u*u
        check('integrated_source_complete_square', QH-square_form-D(flux)/rho**4)
        check('support_primitive_relation',
              (D(rho*rho*u)+rho*rho*p*v/3).subs(D(u),-p*v/3-2*a*u))
        source6 = next(x for x in pins['source_files'] if x['key']=='S6')
        bounds = json.loads((repo_root/source6['path']).read_text())['background_intervals']
        rlo,rhi = map(S.Rational,bounds['rho']['exact'])
        plo,phi = map(S.Rational,bounds['p']['exact'])
        elo,ehi = map(S.Rational,bounds['eta']['exact'])
        pmax,emax = max(abs(plo),abs(phi)),max(abs(elo),abs(ehi))
        ell = S.Rational(1,1000000)
        vpplower = S.Rational(28,9)-S.Rational(100,9)*emax-2*emax**2-S.Rational(80,9)*emax**3-S.Rational(20,9)*emax**4
        epspos = (4/rlo**2+S.Rational(8,3)*pmax*pmax)*(rhi*ell)**2*pmax*pmax/18
        if not (vpplower>3 and epspos<S.Rational(1,100)):
            raise ValueError('Retained exact bounds do not establish the strict moment sign')
        checks.extend(['retained_potential_second_derivative_above_three','strict_moment_positive_part_below_one_percent'])

        reject('wrong_Hamiltonian_factor_two', AH/6+4*DH/3-QH)
        reject('wrong_momentum_source_sign', -radial_direct/6-QM)
        reject('omit_density_Pi_conversion', 2*rho*u*j*Y*Y)
        reject('wrong_adjoint_angular_multiplicity', adjB/rho**3-Jt(N))
        reject('wrong_moment_radial_factor', Jr(rho)+1/rho)
        reject('wrong_compact_kernel_integral_factor', coefficient+3*rho*rho)
        reject('wrong_metric_second_velocity_cross_factor', mean(2*rho*K2[1,1])+4*rho*u*s/3
               -2*rho*(-2*Tr/3-2*u*s/3))
        reject('wrong_square_completion_sign',QH-square_form+D(flux)/rho**4)
        return {'status':'PASS','checks':checks,'check_count':len(checks),
                'wrong_formula_controls':controls,'control_count':len(controls),
                'source_commit':pins['commit'],'sympy':S.__version__,
                'moment_bound':{'Vpp_rational_lower':str(vpplower),
                                'positive_part_fraction':str(epspos),
                                'source_integral_factor_upper':'-49/100',
                                'metric_moment_factor_upper':'-49/150'},
                'scope':'Exact symbolic identities only; no weighted inverse, numerical selected profile, evolution or physical claim.'}


    return run()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', required=True, type=Path,
                        help='Repository containing every path in SOURCE_PINS.json')
    parser.add_argument('--output', type=Path,
                        help='Fresh JSON report path; default writes canonical JSON to stdout')
    args = parser.parse_args()
    repo_root = args.repo_root.resolve(strict=True)
    if not repo_root.is_dir():
        raise ValueError('--repo-root must be a directory')
    if args.output is not None:
        output_path = args.output.resolve()
        if output_path.is_relative_to(ROOT) or output_path.is_relative_to(repo_root):
            raise ValueError('--output must be outside the code and source roots')
    if args.output is not None and args.output.exists():
        raise FileExistsError('Refusing to replace existing report: ' + str(args.output))
    pin_path = ROOT / 'SOURCE_PINS.json'
    pin_bytes = pin_path.read_bytes()
    pins = json.loads(pin_bytes)
    authenticate(repo_root, pins)  # Must precede the SymPy import in evaluate.
    report = evaluate(repo_root, pins)
    if pin_path.read_bytes() != pin_bytes:
        raise ValueError('SOURCE_PINS.json changed during execution')
    authenticate(repo_root, pins)
    raw = json.dumps(report, indent=2) + '\n'
    if args.output is None:
        print(raw, end='')
    else:
        with args.output.open('x', encoding='utf-8') as stream:
            stream.write(raw)
        print(json.dumps({'status': report['status'],
                          'check_count': report['check_count'],
                          'control_count': report['control_count']}))


if __name__ == '__main__':
    main()
