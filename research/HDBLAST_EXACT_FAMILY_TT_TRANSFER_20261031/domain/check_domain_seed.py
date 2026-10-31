#!/usr/bin/env python3
"""Pinned arithmetic/algebra audit; does not certify the analytic PDE theorem."""
from pathlib import Path, PurePosixPath
from fractions import Fraction as F
import argparse
import hashlib
import json
import sys


def require(value, message):
    if not value:
        raise RuntimeError(message)


def checkout_path(repo_root, relative):
    rel = PurePosixPath(relative)
    require(not rel.is_absolute() and '..' not in rel.parts and rel.as_posix() == relative,
            'UNSAFE_CHECKOUT_RELATIVE_PATH')
    path = (repo_root / Path(*rel.parts)).resolve(strict=True)
    require(path.is_relative_to(repo_root) and path.is_file(), 'SOURCE_OUTSIDE_CHECKOUT')
    return path


def authenticate(pins_path, expected, repo_root, overrides_path):
    raw = pins_path.read_bytes()
    require(hashlib.sha256(raw).hexdigest() == expected, 'PINS_AUTHENTICATION_FAILED')
    pins = json.loads(raw)
    require(pins['schema'] == 'hdblast-portable-domain-seed-pins-v2', 'SOURCE_MANIFEST_SCHEMA')
    repo_root = repo_root.resolve(strict=True)
    require(repo_root.is_dir(), 'REPO_ROOT_NOT_DIRECTORY')
    overrides = {} if overrides_path is None else json.loads(overrides_path.read_text())
    require(set(overrides) <= set(pins['sources']), 'UNKNOWN_OVERRIDE')
    data = {}
    for name, row in pins['sources'].items():
        bound_path = checkout_path(repo_root, row['path'])
        path = Path(overrides[name]) if name in overrides else bound_path
        blob = path.read_bytes()
        require(len(blob) == row['bytes'] and hashlib.sha256(blob).hexdigest() == row['sha256'],
                'SOURCE_AUTHENTICATION_FAILED:' + name)
        data[name] = blob
    # Bind to the independently supplied full-capsule external manifest anchor.
    external = '1758d6479eaa4ed3637ef8d0aaf9f61fcbf3d8b0518e5f95fe5ecdd6fd8b979a'
    require(hashlib.sha256(data['full_manifest']).hexdigest() == external,
            'FULL_MANIFEST_AUTHENTICATION_FAILED')
    manifest = json.loads(data['full_manifest'])
    require(len(manifest['files']) == 51, 'FULL_MANIFEST_MEMBER_COUNT')
    root = checkout_path(repo_root, pins['sources']['full_manifest']['path']).parent
    for entry in manifest['files']:
        capsule_relative = str(root.relative_to(repo_root) / entry['path'])
        path = checkout_path(repo_root, capsule_relative)
        blob = path.read_bytes()
        require(len(blob) == entry['bytes'] and hashlib.sha256(blob).hexdigest() == entry['sha256'],
                'FULL_MEMBER_AUTHENTICATION_FAILED:' + entry['path'])
    require(data['full_initial'] == data['initial_exact'], 'EXACT_INITIAL_BYTES_DIFFER')
    return data, len(pins['sources'])


def domain(data, case):
    bg = json.loads(data['background_bounds'])
    ell, background_t, background_z = map(F, (bg['domain']['ell'], bg['domain']['T'], bg['domain']['Z']))
    require(ell == F(1, 10**6), 'ACTUAL_ELL')
    require(bg['domain']['verified_enlargement'] == '[0,T] x [0,Z]', 'BACKGROUND_RECTANGLE')
    rho_lo, rho_hi = map(F, bg['background_intervals']['rho']['exact'])
    p_lo, p_hi = map(F, bg['background_intervals']['p']['exact'])
    require((rho_lo, rho_hi) == (F(129), F(131)) and rho_lo > 0, 'BACKGROUND_RHO')
    require((p_lo, p_hi) == (F(-7, 50000), F(-3, 25000)) and p_hi < 0, 'SCALAR_TRANSVERSALITY')
    packet = data['packet'].decode()
    actual_domain = data['first_actual_domain'].decode()
    require('Cq_\\ell' in packet and 'Ce^{-z}q_\\ell' in packet, 'ACTUAL_DISPLACEMENT_VELOCITY_ANCHORS')
    require('Y=n_1^2-n_2^2' in packet, 'ACTUAL_ANGULAR_ANCHOR')
    require('D={0<=x<=3, r>=0, x+r<=8}' in actual_domain, 'ACTUAL_FIRST_DOMAIN_ANCHOR')
    lo, hi = ell / 2, 5 * ell / 2
    support_lo, support_hi = ell, 2 * ell
    slab_t, collar_z, speed = 3 * ell / 2, 6 * ell, F(2)
    h_cap = ell / 4  # h itself is existential; no numerical flow radius is asserted.
    target_x, target_r = F(5, 4), F(1, 8)
    if case == 'causal_equality':
        speed = F(7, 3)
    elif case == 'packet_touches_annulus':
        lo = support_lo
    elif case == 'collar_touches_annulus':
        h_cap = lo
    elif case == 'target_outside_slab':
        target_x = F(2)
    elif case == 'first_domain_escape':
        collar_z = 7 * ell
    elif case != 'nominal':
        raise RuntimeError('UNKNOWN_CASE')
    require(0 < lo < support_lo <= support_hi < hi, 'STRICT_PACKET_IN_ANNULUS')
    require(0 < h_cap < lo, 'FIXED_COLLAR_IN_ZERO_INITIAL_REGION')
    require(speed > 1, 'OPEN_BACKGROUND_SPEED_NEIGHBORHOOD')
    margin = collar_z - hi - speed * slab_t
    require(margin > 0, 'STRICT_FULL_CORRECTION_CAUSAL_MARGIN')
    require(slab_t < background_t and collar_z < background_z, 'SLAB_IN_BACKGROUND_RECTANGLE')
    require(slab_t / ell < 3 and (slab_t + collar_z) / ell < 8, 'SLAB_IN_ACTUAL_FIRST_TRIANGLE')
    require(0 < target_x * ell < slab_t and 0 < target_r * ell < collar_z, 'TARGET_INTERIOR')
    require((target_r + speed * target_x) * ell < collar_z, 'TARGET_SPEED2_BACKWARD_HISTORY')
    require(F(131, 125000) < 9, 'COLLAR_ON_REGISTERED_REGULAR_BRANCH')
    return {
        'ell': str(ell), 'annulus_z_over_ell': ['1/2', '5/2'],
        'slab_t_over_ell': '3/2', 'Z_over_ell': '6', 'c_star': '2',
        'causal_margin_over_ell': str(margin / ell),
        'actual_first_triangle_margin_over_ell': str(8 - (slab_t + collar_z) / ell),
        'target_scaled': [str(target_x), str(target_r)],
        'target_speed2_initial_reach_over_ell': str(target_r + speed * target_x),
        'h': 'Choose 0<h<min(uniform local flow radius,ell/4); no numerical h is certified.',
    }


def seed_algebra():
    import sympy as s
    require(s.__version__ == '1.14.0', 'SYMPY_VERSION')
    R, b, rho, alpha, Q, Qz, k, p0 = s.symbols('R b rho alpha Q Qz k p0', nonzero=True)
    Omega = rho / R
    phis = p0 * rho / (2 * R**2)
    u = alpha * Q / rho**2
    ut = u * R / b  # exp(-z)=R/b, and ∂t=R∂ambientT on the common initial section.
    qmaster = Omega**2 * u / R**2
    pmaster = Omega**2 * ut / R**3
    Dz = lambda f: s.diff(f, R) * (-R) + s.diff(f, alpha) * k * alpha + s.diff(f, Q) * Qz
    D = lambda f: -Dz(f)
    v_expected = 3 * alpha / (rho**3 * p0) * (Qz + k * Q)
    vt_expected = 3 * alpha / (rho**3 * p0) * R / b * (Qz + (k - 1) * Q)
    scalar = -3 * (D(qmaster) + 4 * qmaster) / (2 * Omega**2 * phis)
    scalar_density = -3 * Omega * (D(pmaster) + 5 * pmaster) / (2 * phis)
    rows = {
        'ambient_profile_ratio': pmaster - qmaster / b,
        'metric_radial': -4 * R**2 * qmaster - (-4 * Omega**2 * u),
        'metric_angular': 2 * R**4 * qmaster - 2 * rho**2 * u,
        'extrinsic_radial': -2 * R**2 * pmaster / Omega - (-2 * Omega**2 * ut / rho),
        'extrinsic_angular': R**4 * pmaster / Omega - rho * ut,
        'scalar_displacement': scalar - v_expected,
        'scalar_density_velocity': scalar_density - Omega**4 * vt_expected / rho,
        'G': qmaster / R**2 - alpha * Q / R**6,
        'K': pmaster / R**2 - alpha * Q / (b * R**6),
    }
    for name, expression in rows.items():
        require(s.cancel(expression) == 0, 'SEED_IDENTITY_FAILED:' + name)
    wrong = {
        'drop_ambient_clock_factor': Omega**2 * ut / R**2 - pmaster,
        'drop_velocity_minus_one': scalar_density - Omega**4 * (3 * alpha / (rho**3 * p0) * R / b * (Qz + k * Q)) / rho,
        'drop_density_conversion': scalar_density - vt_expected / rho,
    }
    for name, expression in wrong.items():
        require(s.cancel(expression) != 0, 'WRONG_SEED_FORMULA_NOT_REJECTED:' + name)
    return {'identities': sorted(rows), 'wrong_formula_controls_rejected': sorted(wrong)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pins', type=Path, default=Path(__file__).with_name('SOURCE_PINS.json'))
    ap.add_argument('--pins-sha256', required=True)
    ap.add_argument('--repo-root', type=Path, required=True)
    ap.add_argument('--overrides', type=Path)
    ap.add_argument('--case', default='nominal')
    args = ap.parse_args()
    data, count = authenticate(args.pins, args.pins_sha256, args.repo_root, args.overrides)
    report = {
        'status': 'PASS_PINNED_DOMAIN_AND_COMMON_SEED_ALGEBRA',
        'sources_authenticated': count, 'full_theorem_members_authenticated': 51,
        'domain': domain(data, args.case), 'seed': seed_algebra(),
        'scope': 'Exact inequalities and common Cauchy seed identities only. Analytic weighted inverse, mixed evolution, finite propagation and uniqueness remain retained mathematical theorems; this is not a numerical certificate of them.',
    }
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(1)
