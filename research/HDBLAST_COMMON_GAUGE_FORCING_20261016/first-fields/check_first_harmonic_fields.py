#!/usr/bin/env python3
"""Unpublished exact finite algebra checks; no actual-field norm enclosure."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile


ERRORS = {
    "clock_sign": "metric component tt",
    "missing_c_squared": "metric component yi",
    "wrong_inward_sign": "physical inward mixed conversion",
    "half_connection_hessian": "connection second-parameter derivative",
    "missing_scalar_connection": "full scalar affine Hessian",
}


def require(ok, msg):
    if not ok:
        raise ValueError(msg)


def eq(a, b, msg):
    require(a == b, msg)


def matrix(n, v=F(0)):
    return [[v for _ in range(n)] for _ in range(n)]


def inv(a):
    n = len(a)
    aug = [list(row) + [F(i == j) for j in range(n)]
           for i, row in enumerate(a)]
    for j in range(n):
        k = next((i for i in range(j, n) if aug[i][j]), None)
        require(k is not None, "nondegenerate generic metric")
        aug[j], aug[k] = aug[k], aug[j]
        scale = aug[j][j]
        aug[j] = [v / scale for v in aug[j]]
        for i in range(n):
            if i != j:
                scale = aug[i][j]
                aug[i] = [v - scale * w for v, w in zip(aug[i], aug[j])]
    return [row[n:] for row in aug]


def mm(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def neg(a):
    return [[-v for v in row] for row in a]


def rationals(seed):
    rng = random.Random(seed)
    return lambda: F(rng.randint(-7, 7), rng.randint(1, 7))


def symmetric(n, get):
    a = matrix(n)
    for i in range(n):
        for j in range(i, n):
            a[i][j] = a[j][i] = get()
    return a


def component_checks(mutation):
    count = 0
    for seed in (61016, 181016, 291016):
        get = rationals(seed)
        rho, c, H, a = F(130), F(5, 3), F(4, 5), F(1, 9)
        u, v, p, Y = [get() for _ in range(4)]
        f, b, d = [get() for _ in range(3)]
        fy, by, dy, ft, bt, dt = [get() for _ in range(6)]
        yg = [get() for _ in range(3)]
        yh = symmetric(3, get)
        m = [[get() for _ in range(3)] for _ in range(3)]
        kap = [[sum(m[k][i] * m[k][j] for k in range(3)) + F(i == j)
                for j in range(3)] for i in range(3)]
        ki = inv(kap)
        dk = [symmetric(3, get) for _ in range(3)]
        dki = [neg(mm(mm(ki, dd), ki)) for dd in dk]
        grad = [sum(ki[i][j] * yg[j] for j in range(3)) for i in range(3)]
        dgrad = [[sum(dki[k][i][j] * yg[j] + ki[i][j] * yh[k][j]
                       for j in range(3)) for i in range(3)] for k in range(3)]
        gam = [[[sum(ki[k][l] * (dk[i][j][l] + dk[j][i][l] - dk[l][i][j])
                         for l in range(3)) / 2
                 for j in range(3)] for i in range(3)] for k in range(3)]
        hess = [[yh[i][j] - sum(gam[k][i][j] * yg[k] for k in range(3))
                 for j in range(3)] for i in range(3)]
        g = matrix(5)
        g[0][0], g[1][1] = F(1), -rho**2
        dg = [matrix(5) for _ in range(5)]
        dg[0][1][1] = -2 * a * rho**2
        for i in range(3):
            for j in range(3):
                g[i+2][j+2] = rho**2 * c**2 * kap[i][j]
                dg[0][i+2][j+2] = 2 * a * g[i+2][j+2]
                dg[1][i+2][j+2] = 2 * H * g[i+2][j+2]
                for k in range(3):
                    dg[k+2][i+2][j+2] = rho**2 * c**2 * dk[k][i][j]
        xi = [f*Y, b*Y] + [d*w for w in grad]
        dxi = [[fy*Y, by*Y] + [dy*w for w in grad],
               [ft*Y, bt*Y] + [dt*w for w in grad]]
        dxi += [[f*yg[k], b*yg[k]] + [d*w for w in dgrad[k]] for k in range(3)]
        long = [[(-4*u*Y if i == j == 0 else
                  2*u*Y*g[i][j] if i > 0 and j > 0 else F(0))
                 for j in range(5)] for i in range(5)]
        direct = [[long[i][j] + sum(xi[k]*dg[k][i][j]
                      + g[k][j]*dxi[i][k] + g[i][k]*dxi[j][k]
                      for k in range(5)) for j in range(5)] for i in range(5)]
        expected = matrix(5)
        expected[0][0] = (-4*u+2*fy)*Y
        expected[0][1] = expected[1][0] = (ft-rho**2*by)*Y
        expected[1][1] = -2*rho**2*(u+a*f+bt)*Y
        if mutation == "clock_sign":
            expected[1][1] = -2*rho**2*(u+a*f-bt)*Y
        for i in range(3):
            c2 = F(1) if mutation == "missing_c_squared" else c**2
            expected[0][i+2] = expected[i+2][0] = (f+rho**2*c2*dy)*yg[i]
            expected[1][i+2] = expected[i+2][1] = rho**2*(c**2*dt-b)*yg[i]
            for j in range(3):
                expected[i+2][j+2] = 2*rho**2*c**2*((u+a*f+H*b)*Y*kap[i][j]
                                                          + d*hess[i][j])
        for i in range(5):
            for j in range(i, 5):
                name = "tt" if i == j == 1 else "yi" if i == 0 and j >= 2 else str((i,j))
                eq(direct[i][j], expected[i][j], "metric component " + name)
                count += 1
        eq(v*Y + xi[0]*p, (v+p*f)*Y, "scalar pullback component")
        count += 1
        # The background-frame conversion is diagonal, with timelike basis
        # vector rho^-1 partial_t; no signature sign is inserted in evaluation.
        scale = [F(1), 1/rho] + [1/(rho*c)]*3
        eq(direct[0][1]*scale[1], (ft/rho-rho*by)*Y, "frame y0 normalization")
        eq(direct[1][2]*scale[1]*scale[2], (c*dt-b/c)*yg[0], "frame 0a normalization")
        count += 2
    return count


def radial_checks(mutation):
    # Independent exact evaluation on rho(y)=2+y+y^2 and F(r)=r^3+2r^2-r.
    # y coordinate is converted by r_y=-(ell*rho)^-1. The forward automatic
    # derivatives below use chain-rule jets, not the displayed radial formulas.
    ell, rho, rhoy, rhoyy = F(1, 1000000), F(7, 3), F(4, 3), F(2)
    a = rhoy/rho
    ay = rhoyy/rho-a*a
    r = F(5, 4)
    L = 1/(ell*rho)
    ry = -L
    ryy = rhoy/(ell*rho**2)
    ryyy = rhoyy/(ell*rho**2)-2*rhoy**2/(ell*rho**3)
    Fr, Frr, Frrr = 3*r*r+4*r-1, 6*r+4, F(6)
    direct = [ry*Fr, ryy*Fr+ry*ry*Frr,
              ryyy*Fr+3*ry*ryy*Frr+ry**3*Frrr]
    stated = [-L*Fr, L*L*Frr+a*L*Fr,
              -L**3*Frrr-3*a*L*L*Frr+(ay-a*a)*L*Fr]
    for n in range(3):
        eq(direct[n], stated[n], "proper radial conversion order " + str(n+1))
    fx, br = F(3, 7), F(4, 9)
    by = (1 if mutation == "wrong_inward_sign" else -1)*br/(ell*rho)
    eq(fx/ell-rho**2*by, (fx+rho*br)/ell, "physical inward mixed conversion")
    return 4


def connection_and_scalar_checks(mutation):
    count = 0
    n = 5
    gi = [[F(0) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        gi[i][i] = F(-1 if i == 1 else 1)
    for seed in (31016, 41016, 51016):
        get = rationals(seed)
        h = symmetric(n, get)
        dh = [symmetric(n, get) for _ in range(n)]
        hi = mm(mm(gi, h), gi)
        hu = mm(gi, h)
        inv2 = mm(mm(gi, h), mm(mm(gi, h), gi))
        P = [[[sum(gi[k][l]*(dh[i][j][l]+dh[j][i][l]-dh[l][i][j])
                         for l in range(n))/2 for j in range(n)]
                         for i in range(n)] for k in range(n)]
        J = [[[-2*sum(hu[k][l]*P[l][i][j] for l in range(n))
                 for j in range(n)] for i in range(n)] for k in range(n)]
        # Direct inverse-metric series times the first metric derivative.
        gamma2coef = [[[-sum(hi[k][l]*(dh[i][j][l]+dh[j][i][l]-dh[l][i][j])
                         for l in range(n))/2 for j in range(n)]
                         for i in range(n)] for k in range(n)]
        factor = 1 if mutation == "half_connection_hessian" else 2
        for k in range(n):
            for i in range(n):
                for j in range(i,n):
                    eq(J[k][i][j], factor*gamma2coef[k][i][j],
                       "connection second-parameter derivative")
                    count += 1
        phi0g, chig = [get() for _ in range(n)], [get() for _ in range(n)]
        phi0h, chih = symmetric(n,get), symmetric(n,get)
        chi, V3 = get(), get()
        C = [sum(gi[i][j]*P[k][i][j] for i in range(n) for j in range(n))
             for k in range(n)]
        # Expand g(eps)^ab (partial_ab phi(eps)-Gamma(eps)^c_ab partial_c phi(eps)).
        direct = 2*sum(inv2[i][j]*phi0h[i][j]-hi[i][j]*(chih[i][j]
                    -sum(P[k][i][j]*phi0g[k] for k in range(n)))
                    +gi[i][j]*(-sum(P[k][i][j]*chig[k]
                    +gamma2coef[k][i][j]*phi0g[k] for k in range(n)))
                    for i in range(n) for j in range(n))-V3*chi**2
        connection_factor = 0 if mutation == "missing_scalar_connection" else 2
        stated = 2*sum(inv2[i][j]*phi0h[i][j]-hi[i][j]*chih[i][j]
                       for i in range(n) for j in range(n))
        stated += connection_factor*sum(hi[i][j]*P[k][i][j]*phi0g[k]
                       for i in range(n) for j in range(n) for k in range(n))
        stated += 2*sum(hu[k][l]*C[l]*phi0g[k] for k in range(n) for l in range(n))
        stated -= 2*sum(C[k]*chig[k] for k in range(n)) + V3*chi**2
        eq(direct, stated, "full scalar affine Hessian")
        count += 1
        qc_direct = [sum(gi[i][j]*J[k][i][j]-2*hi[i][j]*P[k][i][j]
                       for i in range(n) for j in range(n)) for k in range(n)]
        qc_stated = [-2*sum(hu[k][l]*C[l] for l in range(n))
                     -2*sum(hi[i][j]*P[k][i][j] for i in range(n) for j in range(n))
                     for k in range(n)]
        eq(qc_direct, qc_stated, "quadratic relative harmonic vector")
        count += 1
        # Flat-background two-jets give an independent inverse/Christoffel
        # epsilon expansion of every Ricci Hessian component. The covariant
        # formula in the note applies to curved backgrounds by connection
        # difference; these finite witnesses check its algebraic coefficients.
        ddh = [[None for _ in range(n)] for _ in range(n)]
        for aa in range(n):
            for bb in range(aa,n):
                ddh[aa][bb] = ddh[bb][aa] = symmetric(n,get)
        dP = [[[[sum(gi[k][l]*(ddh[aa][i][j][l]+ddh[aa][j][i][l]
                            -ddh[aa][l][i][j]) for l in range(n))/2
                  for j in range(n)] for i in range(n)] for k in range(n)]
                  for aa in range(n)]
        dhu = [mm(gi,dd) for dd in dh]
        dJ = [[[[-2*sum(dhu[aa][k][l]*P[l][i][j]+hu[k][l]*dP[aa][l][i][j]
                         for l in range(n)) for j in range(n)]
                         for i in range(n)] for k in range(n)] for aa in range(n)]
        dInv1 = [neg(mm(mm(gi,dd),gi)) for dd in dh]
        dGamma2coef = [[[[sum(dInv1[aa][k][l]*(dh[i][j][l]+dh[j][i][l]-dh[l][i][j])
                         -hi[k][l]*(ddh[aa][i][j][l]+ddh[aa][j][i][l]-ddh[aa][l][i][j])
                         for l in range(n))/2 for j in range(n)]
                         for i in range(n)] for k in range(n)] for aa in range(n)]
        V1, V2 = get(), get()
        for i in range(n):
            for j in range(i,n):
                pp = sum(P[k][k][l]*P[l][i][j]-P[k][j][l]*P[l][i][k]
                         for k in range(n) for l in range(n))
                direct_ricci = 2*(sum(dGamma2coef[k][k][i][j]-dGamma2coef[j][k][i][k]
                                  for k in range(n))+pp)
                stated_ricci = sum(dJ[k][k][i][j]-dJ[j][k][i][k] for k in range(n))+2*pp
                eq(direct_ricci, stated_ricci, "full Ricci affine Hessian")
                direct_matter = -2*chig[i]*chig[j]-F(4,3)*(V2*chi**2*gi[i][j]/2+V1*chi*h[i][j])
                stated_matter = -2*chig[i]*chig[j]-F(2,3)*(V2*chi**2*gi[i][j]+2*V1*chi*h[i][j])
                eq(direct_matter, stated_matter, "trace-reversed potential affine Hessian")
                count += 2
        # A nontrivial exactly trace-free first connection tests the harmonic
        # simplification independently.
        Pt = [[list(row) for row in upper] for upper in P]
        for k in range(n):
            Pt[k][0][0] = -sum(gi[i][i]*Pt[k][i][i] for i in range(1,n))
        Ct = [sum(gi[i][j]*Pt[k][i][j] for i in range(n) for j in range(n))
              for k in range(n)]
        eq(Ct, [F(0)]*n, "nontrivial first harmonic connection trace")
        scalar_trace_terms = 2*sum(hu[k][l]*Ct[l]*phi0g[k] for k in range(n) for l in range(n))
        scalar_trace_terms -= 2*sum(Ct[k]*chig[k] for k in range(n))
        eq(scalar_trace_terms, F(0), "harmonic scalar trace cancellation")
        count += 2
    return count


def checks(mutation):
    return {
        "component_and_frame_identities": component_checks(mutation),
        "physical_scaled_conversion_identities": radial_checks(mutation),
        "connection_scalar_and_gauge_identities": connection_and_scalar_checks(mutation),
    }


def controls():
    result = []
    with tempfile.TemporaryDirectory(prefix="first-fields-controls-") as td:
        for mutation, reason in ERRORS.items():
            dest = Path(td) / (mutation + ".json")
            cmd = [sys.executable]
            if not __debug__:
                cmd.append("-O")
            cmd += [str(Path(__file__).resolve()), "--inject-error", mutation,
                    "--output", str(dest)]
            p = subprocess.run(cmd, text=True, capture_output=True)
            require(p.returncode != 0 and reason in p.stderr and not dest.exists(),
                    "mutation-specific rejection failed: " + mutation)
            result.append({"mutation": mutation, "expected_reason": reason,
                           "rejected_without_output": True})
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--run-controls", action="store_true")
    ap.add_argument("--inject-error", choices=list(ERRORS))
    args = ap.parse_args()
    require(not args.output.exists(), "refuse to overwrite output")
    counts = checks(args.inject_error)
    report = {
        "status": "PASS_UNPUBLISHED_FINITE_ALGEBRA_ONLY",
        "scope": {
            "actual_registered_field_evaluated": False,
            "actual_transformed_source_enclosed": False,
            "second_order_evolution_solved": False,
            "finite_generic_jet_witnesses_supplement_analytic_derivation": True,
        },
        "checks": counts,
        "conventions": {
            "first_transformation": "h_prime=h+Lie_xi(g0); chi_prime=chi+Lie_xi(phi0)",
            "raw_vector": "xi_y_contravariant=fY; xi_t_contravariant=bY; xi_angular=d_gradY",
            "radial": "dy=-rho_dz; y_derivative=-r_derivative/(ell*rho)",
            "quadratic": "second_parameter_derivative; epsilon_squared_coefficient_is_half",
        },
        "controls": controls() if args.run_controls else [],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True)+"\n")


if __name__ == "__main__":
    try:
        main()
    except ValueError as exc:
        sys.exit("REJECT: " + str(exc))
