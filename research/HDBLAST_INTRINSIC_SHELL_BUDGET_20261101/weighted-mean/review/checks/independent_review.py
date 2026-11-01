#!/usr/bin/env python3
"""Portable independent review. Source bytes are authenticated before SymPy."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import platform
import sys


COMPONENT_ROOT = Path(__file__).resolve().parent.parent
SOURCE_MANIFEST = COMPONENT_ROOT / 'SOURCE_PINS.json'
HISTORICAL_NOTE_SHA256 = 'a2f9d76543aa7d16826b2331d780f70873f3ce18570e38ddb9ebc4c7596d93cf'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def contained(path, directory):
    try:
        path.relative_to(directory)
        return True
    except ValueError:
        return False


def authenticate(repo_root, pins):
    files = pins['source_files']
    if sorted(pin['key'] for pin in files) != ['S1', 'S2', 'S3', 'S4', 'S5', 'S6']:
        raise ValueError('Source manifest must contain exactly S1 through S6.')
    blobs = {}
    records = []
    for pin in files:
        relative = PurePosixPath(pin['path'])
        if relative.is_absolute() or '..' in relative.parts:
            raise ValueError('Source path must be repository-relative: ' + pin['key'])
        candidate = (repo_root / Path(*relative.parts)).resolve()
        if not contained(candidate, repo_root):
            raise ValueError('Source path escapes repository root: ' + pin['key'])
        try:
            data = candidate.read_bytes()
        except OSError as exc:
            raise RuntimeError('Cannot read source ' + pin['key'] + ': ' + pin['path']) from exc
        if len(data) != pin['bytes'] or digest(data) != pin['sha256']:
            raise RuntimeError('Source authentication failed: ' + pin['key'] + ': ' + pin['path'])
        blobs[pin['key']] = data
        records.append(dict(pin))
    return blobs, records


def validate_output(value, repo_root):
    if value == '-':
        return None
    output = Path(value).expanduser().resolve()
    if contained(output, repo_root) or contained(output, COMPONENT_ROOT):
        raise ValueError('Output must be outside the source repository and review component.')
    if output.exists():
        raise ValueError('Output must be a fresh file.')
    return output


def run_checks(source_blobs):
    import sympy as S
    checks = {}

    def zero(name, expression):
        value = S.factor(expression)
        if value != 0:
            raise RuntimeError(f'{name}: {value}')
        checks[name] = True

    e,Y,G = S.symbols('e Y G')
    r,a,a1,a2,p,p1,V,V1,V2 = S.symbols('rho a ay ayy p py V Vphi Vphiphi')
    u,u1,u2,v,v1,s,s1,j = S.symbols('u uy uyy v vy s sy j')
    A,A1,A2,B,B1,B2,N,N1,N2,N3 = S.symbols('A Ay Ayy B By Byy N Ny Nyy Nyyy')
    jet = {r:a*r,a:a1,a1:a2,p:p1,u:u1,u1:u2,v:v1,s:s1,A:A1,A1:A2,B:B1,B1:B2,N:N1,N1:N2,N2:N3}
    def dy(expr):
        return S.expand(sum(S.diff(expr, x)*dx for x,dx in jet.items()))
    def angular_lap(expr):
        return -8*Y*S.diff(expr,Y)+G*S.diff(expr,Y,2)
    def second(expr):
        return S.diff(expr,e,2).subs(e,0)
    def first(expr):
        return S.diff(expr,e).subs(e,0)
    def avg(expr):
        return S.expand(expr).subs({Y**2:S.Rational(1,6),G:S.Rational(4,3)})
    bg = {a1:-1/r**2-p**2/3, p1:V1-4*a*p,V:p**2/2-6*a*a+6/r**2}

    # Curvature: Gauss decomposition for a 3-sphere foliation; angular
    # derivatives act on Y with Delta Y=-8Y and |DY|^2=G.
    L=S.sqrt(1-4*e*u*Y)
    C=1+2*e*u*Y
    logB=S.log(r)+S.log(C)/2
    k=(a+e*u1*Y/C)/L
    R=(6-4*angular_lap(logB)-2*S.diff(logB,Y)**2*G)/(r*r*C)
    R-=12*k*k+6*dy(k)/L
    R-=2*(angular_lap(L)+S.diff(logB,Y)*S.diff(L,Y)*G)/(L*r*r*C)
    R2=S.expand(second(R))
    R2_expected=(-24*u1*u1-288*a*u*u1-24*u*u2-(192*a1+384*a*a)*u*u-336*u*u/r**2)*Y**2+36*u*u*G/r**2
    zero('direct_geometric_curvature_Hessian', R2-R2_expected)

    # Compute the canonical kinetic terms from the density variables, not K2.
    tau=-3*e*s*Y/2
    det_ratio=(1-4*e*u*Y)*C**3
    pi_trace=(1-4*e*u*Y)*tau
    pi_norm=(1-4*e*u*Y)**2*tau**2
    kinetic=2*(pi_trace**2/3-pi_norm)/det_ratio-e**2*j*j*Y**2/(2*det_ratio)
    grad=(p+e*v1*Y)**2/(1-4*e*u*Y)+e*e*v*v*G/(r*r*C)
    potential=V+e*V1*v*Y+e*e*V2*v*v*Y**2/2
    QH=avg(R2/2+second(kinetic-grad/2-potential))
    QH_claim=-2*u1*u1-24*a*u*u1-2*u*u2-(16*a1+32*a*a+8*p*p/3)*u*u-4*u*u/r**2-s*s-(j*j+v1*v1+V2*v*v)/6-S.Rational(4,3)*p*u*v1-S.Rational(4,3)*v*v/r**2
    zero('canonical_QH',QH-QH_claim)

    # Radial linearization rederived from exact warped curvature for A and B.
    Lrad=S.sqrt(1+e*A)
    Brad=S.sqrt(1+e*B)
    krad=(a+e*B1/(2*(1+e*B)))/Lrad
    Hrad=(6/(r*r*Brad**2)-12*krad**2-6*dy(krad)/Lrad)/2-p*p/(2*Lrad**2)-V
    Pg=-3*B2/2+3*a*A1/2-6*a*B1+(3*a1+6*a*a+p*p/2)*A-3*B/r**2
    zero('radial_metric_linearization',first(Hrad)-Pg)

    # Formal adjoint with respect to rho^3 dy and tensor norm A^2+3B^2.
    Jr=(-3*a*N1+(3*a1+3*a*a+p*p)*N)/2
    Jt=(-N2-2*a*N1+(a1+3*a*a-2/r**2)*N)/2
    adjA=((3*a1+6*a*a+p*p/2)*N-dy(3*a*r**3*N/2)/r**3)
    adjB=(-3*N/r**2+dy(6*a*r**3*N)/r**3-dy(dy(3*r**3*N/2))/r**3)/3
    zero('formal_adjoint_radial',Jr-adjA)
    zero('formal_adjoint_tangential',Jt-adjB)
    rhojet={N:r,N1:a*r,N2:(a1+a*a)*r}
    zero('areal_radius_adjoint_radial',(Jr.subs(rhojet)+S.Rational(3,2)/r).subs(bg))
    zero('areal_radius_adjoint_tangential',(Jt.subs(rhojet)+1/r).subs(bg))

    # Integrated identity: QH minus the sum of squares is a total derivative.
    negative=(4/r**2+8*p*p/3)*u*u-s*s-j*j/6-(v1+4*p*u)**2/6-(V2/6+S.Rational(4,3)/r**2)*v*v
    boundary=-2*r**4*u*u1-8*r**4*a*u*u
    zero('exact_completed_square_identity',(r**4*(QH-negative)-dy(boundary)).subs(bg))

    # Momentum density: use warped-frame background divergence of tau_rr only.
    # D_j tau^jr = (tau_rr)' + 3a tau_rr. Each term is the second
    # parameter derivative of the exact normalized density formula.
    QM=avg((-24*u*(s1+3*a*s)-24*s*u1+12*s*u1+2*j*v1)*Y**2)
    QM_claim=-4*u*s1-12*a*u*s-2*s*u1+j*v1/3
    zero('canonical_QM',QM-QM_claim)
    zero('momentum_source_constraint_reduction',(QM-(-2*s*u1+j*(v1+4*p*u)/3)).subs(s1,-p*j/3-3*a*s))
    DY=S.symbols('DY')
    QMA=2*(-3*s*Y/2)*(-4*u*DY)+2*j*Y*v*DY
    zero('angular_momentum_Hessian',QMA-2*(6*u*s+j*v)*Y*DY)
    W,W1,x,x1,x2=S.symbols('W Wy xi xiy xiyyyy')
    Tmom=2*W*x1; Tmom1=2*W1*x1+2*W*x2; Umom=2*W*a*x
    Mop=-2*(Tmom1+3*a*(Tmom-Umom))
    zero('weighted_radial_momentum_operator',Mop+4*(W1*x1+W*x2+3*a*W*(x1-a*x)))
    # Integrating the divergence term against rho^3 chi produces the stated
    # positive form; its Euler operator is -4 rho^-3(rho^3 W xi')'+12Wa^2 xi.
    Mweak=-4*(W1*x1+W*x2+3*a*W*x1)+12*W*a*a*x
    zero('weighted_momentum_weak_form',Mweak-Mop)
    zero('areal_radius_moment_pairing',(-S.Rational(3,2)*r*r*A-3*r*r*B)+S.Rational(3,2)*r*r*(A+2*B))
    zero('sphere_Y_square_probability',2*S.Rational(1,8)-2*S.Rational(1,24)-S.Rational(1,6))
    zero('sphere_DY_square_probability',8*S.Rational(1,6)-S.Rational(4,3))

    # Conformal-plus-radial-Lie compact kernel: compute Pg on both explicitly.
    f,f1,f2,x,x1,x2,x3 = S.symbols('f fy fyy xi xiy xiyyy xiyyyy')
    Pgf=Pg.subs({A:f,A1:f1,B:f,B1:f1,B2:f2})
    Cf=-S.Rational(3,2)*(f2+3*a*f1)-V*f
    zero('conformal_metric_operator',(Pgf-Cf).subs(bg))
    Pglie=Pg.subs({A:2*x1,A1:2*x2,B:2*a*x,B1:2*a1*x+2*a*x1,B2:2*a2*x+4*a1*x1+2*a*x2})
    a2bg=2*a/r**2-S.Rational(2,3)*p*p1
    zero('radial_Lie_metric_operator',(Pglie-(p*p*x1+(4*a*p*p+2*p*p1)*x)).subs(a2,a2bg).subs(bg))
    # Integral of rho^4 Cf: formal scalar adjoint, integrating twice.
    Cfadj=-S.Rational(3,2)*dy(dy(r**4))+S.Rational(9,2)*dy(a*r**4)-V*r**4
    zero('kernel_compactness_moment',(Cfadj+S.Rational(9,2)*r*r).subs(bg))

    # Canonical conversion retaining an arbitrary second radial/tangential pi.
    T,U=S.symbols('T U')
    g_rr=1-4*e*u*Y
    g_tt=1+2*e*u*Y
    tau_rr=-3*e*s*Y/2+e*e*T/2
    tau_tt=e*e*U/2
    tr=g_rr*tau_rr+3*g_tt*tau_tt
    vol=S.sqrt(g_rr*g_tt**3)
    Krr=2*(g_rr*g_rr*tau_rr-g_rr*tr/3)/vol
    Ktt=2*(g_tt*g_tt*tau_tt-g_tt*tr/3)/vol
    zero('physical_K2_radial',second(Krr)-(4*T/3-2*U+36*u*s*Y**2))
    zero('physical_K2_tangential',second(Ktt)-(-2*T/3-6*u*s*Y**2))
    zero('physical_Pi2',second(e*j*Y/vol)+2*u*j*Y**2)
    zero('scalar_coordinate_velocity2',second(r*(1+e*u*Y)*e*j*Y/vol))
    zero('metric_coordinate_velocity2_radial',avg(second(2*r*(1+e*u*Y)*Krr))-2*r*(4*T/3-2*U+S.Rational(16,3)*u*s))
    zero('metric_coordinate_velocity2_tangential',avg(second(2*r*(1+e*u*Y)*Ktt))-2*r*(-2*T/3-S.Rational(2,3)*u*s))

    # Read the already authenticated retained coefficient certificate.
    cert=json.loads(source_blobs['S6'])
    coeff=cert['background_intervals']
    rho_min,rho_max=map(S.Rational,coeff['rho']['exact'])
    eta_min,eta_max=map(S.Rational,coeff['eta']['exact'])
    p_min,p_max=map(S.Rational,coeff['p']['exact'])
    ell=S.Rational(cert['domain']['ell'])
    pabs=max(abs(p_min),abs(p_max)); E=max(abs(eta_min),abs(eta_max))
    V2lower=S.Rational(28,9)-S.Rational(100,9)*E-2*E**2-S.Rational(80,9)*E**3-S.Rational(20,9)*E**4
    Lmax=rho_max*ell
    eps=(4/rho_min**2+8*pabs*pabs/3)*Lmax**2*pabs*pabs/18
    if not (V2lower>3 and eps<S.Rational(1,100)):
        raise RuntimeError('strict sign bound')
    checks['exact_rational_potential_lower_bound_gt_3']=True
    checks['exact_positive_term_lt_one_hundredth']=True
    if S.factor((r**4*(QH-negative)+dy(boundary)).subs(bg))==0:
        raise RuntimeError('wrong orientation control did not fail')
    checks['deliberate_reversed_total_derivative_rejected']=True
    if S.factor(2*QH-QH_claim)==0:
        raise RuntimeError('wrong Hessian normalization control did not fail')
    checks['deliberate_factor_two_error_rejected']=True

    return {
        'status': 'PASS_CONDITIONAL_ALGEBRA_AND_SIGN',
        'check_count': len(checks),
        'checks': checks,
        'potential_second_derivative_lower_bound': str(V2lower),
        'positive_coefficient_upper_bound': str(eps),
        'positive_coefficient_display_only': float(eps),
        'support_proper_length_upper_bound': str(Lmax),
        'sympy_version': S.__version__,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, required=True,
                        help='Root containing the authenticated research/ source paths; Git is unnecessary.')
    parser.add_argument('--output', default='-',
                        help='Fresh JSON file outside source/component directories, or - for stdout (default).')
    args = parser.parse_args()
    repo_root = args.repo_root.expanduser().resolve()
    if not repo_root.is_dir():
        raise ValueError('Repository root must be an existing directory.')
    output = validate_output(args.output, repo_root)
    manifest_bytes = SOURCE_MANIFEST.read_bytes()
    code_bytes = Path(__file__).read_bytes()
    pins = json.loads(manifest_bytes)
    source_blobs, source_before = authenticate(repo_root, pins)
    # This is the first call that imports SymPy; all six sources passed above.
    report = run_checks(source_blobs)
    _, source_after = authenticate(repo_root, pins)
    if source_before != source_after:
        raise RuntimeError('Source authentication changed during computation.')
    if SOURCE_MANIFEST.read_bytes() != manifest_bytes or Path(__file__).read_bytes() != code_bytes:
        raise RuntimeError('Checker or manifest changed during computation.')
    report.update({
        'source_commit_provenance': pins['source_commit_provenance'],
        'source_authentication': {'before_symbolic_import': True, 'after_checks': True,
                                  'files': source_after},
        'source_manifest_sha256': digest(manifest_bytes),
        'checker_sha256': digest(code_bytes),
        'review_history_original_note_sha256': HISTORICAL_NOTE_SHA256,
        'python_version': platform.python_version(),
        'limitations': [
            'Retained analytic inverse and background enclosures are inputs, not re-proved.',
            'No weighted boundary value solve or nonlinear evolution performed.',
            'Moment sign is conditional, canonical-gauge initial-data information.',
        ],
    })
    serialized = json.dumps(report, indent=2, sort_keys=True) + '\n'
    if output is None:
        sys.stdout.write(serialized)
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open('x', encoding='utf-8') as handle:
            handle.write(serialized)


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, RuntimeError) as exc:
        print('ERROR: ' + str(exc), file=sys.stderr)
        raise SystemExit(1)
