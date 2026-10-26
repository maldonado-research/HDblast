#!/usr/bin/env python3
"""Bounded actual total TT coefficient two-jet on an exterior backward cone.

The source module is authenticated against an explicit external CLI pin before
import. Every retained capsule is authenticated as a complete tree first.
All endpoints, quadrature weights, decisions and closure arithmetic are rational.
"""
import argparse
import hashlib
import importlib.util
import json
import sys
from fractions import Fraction as F
from pathlib import Path

sys.dont_write_bytecode = True
LABELS = {
    'current': 'HDBLAST_COUPLED_MEAN_PILOT_20261025',
    'I': 'HDBLAST_INPUT_PRECISION_20261023',
    'P': 'HDBLAST_COMMON_GAUGE_FORCING_20261016',
    'N': 'HDBLAST_NUMERICAL_RESPONSE_INPUTS_20261016',
    'G': 'HDBLAST_GEOMETRIC_BULK_SOURCES_20261016',
    'C': 'HDBLAST_FORCED_RESPONSE_BOUNDS_20261023',
    'T': 'HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023',
    'F': 'HDBLAST_FULL_ANGULAR_RESPONSE_20261026',
    'precision35': 'HDBLAST_SHELL_FEEDBACK_PRECISION_20261025',
}
PINS = {
    'current': 'b932157f308c9912fa83ee06990f1132a4e75a498aafd66eb99dabcc30b256de',
    'I': 'c2e7bd40e94422c5fd98820941231c0da63a8392b14a5a14c4340eabdde37f61',
    'P': '8bf173293b5cecf7f78e138ac5004abd9c7b1821ef29addb7ecfbc2a6409b4f6',
    'N': '928f3e4e58d222a3c7e1274babf9e94930ec501b74cf6ba04209f97c6bbded44',
    'G': '3ff1f95c151d2da6a1ab5e9b54718c2bb325c10d9c62a7ed9f97ef2bb36f0bb2',
    'C': '72924328a4e79cc77f0da90999b8dd56956bfdebd5e614e15c3c0498c3a8ae76',
    'T': 'cffb213a69f7f9a27932e990ffdecb7569d1007e5e5baeda708d8dc39d0cce72',
    'F': 'a35f84119e84bf2ea8a712ded86ff79977d3ee91ac25f91ccbd202d379c859ac',
    'precision35': '1a5ee8aa22953f949da000ba03da3b5a0872e4467a0fbaec9f2d3fe7fba11ec4',
}
T, R = F(1, 4), F(5, 4)
DECLARATION = {
    'first': {'cone_n': 4, 'cone_cells': 20, 'segments_per_characteristic': 8},
    'maximum': {'cone_n': 8, 'cone_cells': 72, 'segments_per_characteristic': 16},
    'refinement_gate': 'Exactly one maximum-budget refinement if no complete apex second derivative interval excludes zero in the first run.',
    'primary_acceptance': 'Finite outward actual total TT two-jet intervals, differentiated scalar resolvent and PDE closure, all source uncertainty retained.',
    'diagnostic_nonzero_gate': 'At least one of u_xx,u_xr,u_rr excludes zero; failure is recorded as inconclusive and does not invalidate the finite two-jet enclosure.',
    'no_further_refinement': True,
    'probes': {'points': 3, 'whole_cells': 2},
}


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def encode(value):
    if hasattr(value, 'lo') and hasattr(value, 'hi'):
        return [str(value.lo), str(value.hi)]
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    return value


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def authenticate(repo):
    """Independent complete-tree authentication before any project import."""
    result = {}
    for key, label in LABELS.items():
        root = repo / 'research' / label
        need(root.is_dir() and not root.is_symlink(), 'ROOT_TYPE:' + key)
        manifest = root / 'MANIFEST.json'
        need(not manifest.is_symlink() and sha(manifest) == PINS[key], 'MANIFEST:' + key)
        seen = {'MANIFEST.json'}
        for entry in json.loads(manifest.read_text())['files']:
            rel = Path(entry['path'])
            need(not rel.is_absolute() and '..' not in rel.parts and str(rel) not in seen, 'MEMBER_PATH:' + key)
            seen.add(str(rel))
            member = root / rel
            need(member.is_file() and not member.is_symlink(), 'MEMBER_TYPE:' + key + '/' + str(rel))
            need(member.stat().st_size == entry['bytes'] and sha(member) == entry['sha256'], 'MEMBER_BYTES:' + key + '/' + str(rel))
        for member in root.rglob('*'):
            need(not member.is_symlink(), 'TREE_SYMLINK:' + key)
            if member.is_file():
                need(str(member.relative_to(root)) in seen, 'UNLISTED_FILE:' + key + '/' + str(member.relative_to(root)))
        result[key] = {'manifest_sha256': PINS[key], 'files_including_manifest': len(seen)}
    return result


def guard_output(repo, output):
    need(not output.exists() and not output.is_symlink(), 'FRESH_OUTPUT_REQUIRED')
    target = output.resolve()
    for label in LABELS.values():
        root = (repo / 'research' / label).resolve()
        need(target != root and root not in target.parents, 'OUTPUT_INSIDE_FROZEN_CAPSULE')
    # Scratch evaluators remain external. Once packaged under research/<capsule>,
    # protect the entire containing capsule even with another --repo-root.
    for parent in Path(__file__).resolve().parents:
        if parent.parent.name == 'research':
            need(target != parent and parent not in target.parents, 'OUTPUT_INSIDE_PACKAGED_CAPSULE')
            break


def intersection(I, first, second, message):
    lo, hi = max(first.lo, second.lo), min(first.hi, second.hi)
    need(lo <= hi, message)
    return I(lo, hi)


def strict_sign(interval):
    return 'positive' if interval.lo > 0 else ('negative' if interval.hi < 0 else 'inconclusive')


def operator(ctx, x=None, r=None):
    """Analytic ODE/Taylor jets, plus a separate physical-identity audit."""
    base = ctx.base
    I, jets, ell = base.I, base.jets, base.ell
    x = I(0, T) if x is None else I.cv(x)
    r = I(R-T, R+T) if r is None else I.cv(r)
    rho0, eta0, p0, a0, h0, *_ = jets.core.background(r, base.data, 'none')
    rho, eta, p, a = map(jets.radial, (rho0, eta0, p0, a0))
    c, Ht, _ = jets.temporal(x)
    ay = -1/rho**2-p**2/3
    b = -(ay+a**2)
    d = 1/rho**2-a**2
    V = jets.core.poly(jets.core.U, eta)
    h = Ht/rho
    A = -3*ell*Ht
    B = -3*ell*rho*a
    D = ell**2*(-6/c**2+rho**2*(-2*a**2+2*h**2-2*b-8*d+F(4, 3)*V))
    values = {'A': A.val(), 'B': B.val(), 'D': D.val()}
    derivatives = {'Ar': A.derivative(0, 1), 'Br': B.derivative(0, 1), 'Dr': D.derivative(0, 1)}
    # These first-physical-derivative identities use a different algebraic form.
    rv, pv, av, ayv, Vv = rho.val(), p.val(), a.val(), ay.val(), V.val()
    Vp = jets.core.poly(jets.core.UP[1], eta.val())
    py = pv*(jets.radial(h0).val()-4*av)
    S = 8*av**2-F(2, 3)*pv**2+F(4, 3)*Vv
    independent = {
        'Ar': I(0),
        'Br': 3*ell**2*rv**2*(av**2+ayv),
        'Dr': -ell**3*rv**3*(2*av*S+16*av*ayv-F(4, 3)*pv*py+F(4, 3)*pv*Vp),
    }
    audited = {key: intersection(I, derivatives[key], independent[key], 'COEFFICIENT_RADIAL_AUDIT:' + key) for key in derivatives}
    caps = {key: value.mag() for key, value in values.items()}
    rcaps = {key: value.mag() for key, value in audited.items()}
    sigma = caps['A']+caps['B']+T*caps['D']/2
    sigma_r = rcaps['Ar']+rcaps['Br']+T*rcaps['Dr']/2
    q, qr = T*sigma, T*sigma_r
    need(0 <= q < 1 and qr >= 0, 'SCALAR_C1_RESOLVENT_GATE')
    return {
        'box': {'x': x, 'r': r},
        'background_value_intervals': {'rho': rv, 'a': av, 'ell': ell, 'L': 1/(ell*rv)},
        'coefficient_intervals': values,
        'radial_derivative_analytic_jet_intervals': derivatives,
        'radial_derivative_physical_identity_intervals': independent,
        'radial_derivative_intersections': audited,
        'absolute_caps': caps,
        'radial_derivative_absolute_caps': rcaps,
        'sigma': sigma, 'sigma_r': sigma_r, 'q': q, 'q_r': qr,
        'norm': 'X(w)=max(2 sup|w|/T,sup|w_x|,sup|w_r|), T=1/4',
        'radial_differentiation': 'Background normalized Taylor coefficients follow exact ODE recurrences, converted by ell^n. No numerical differentiation or differentiation of selected value intervals.',
        'differentiated_equation': 'v_xx-v_rr=A v_x+B v_r+D v+Fr+Ar ux+Br ur+Dr u, v=u_r; exact zero Cauchy data.',
    }


def evaluate_source(ctx, x, r, audit=False):
    value, radial = ctx.forcing_with_radial(x, r, audit=audit)
    need(value.lo <= value.hi and radial.lo <= radial.hi, 'SOURCE_INTERVAL_ORDER')
    # The source uses the equivalent integer-dyadic ring; normalize explicitly
    # before combining with analytic-background Fraction intervals.
    return ctx.base.I(value.lo, value.hi), ctx.base.I(radial.lo, radial.hi)


def probes(ctx):
    I = ctx.base.I
    rows = []
    cells = [
        ('initial_packet_point', F(0), R, None),
        ('middle_cone_point', F(1, 8), R, None),
        ('apex_point', T, R, None),
        ('initial_left_N4_cell', F(0), F(1), F(1, 16)),
        ('middle_N4_cell', F(1, 8), F(19, 16), F(1, 16)),
    ]
    for label, x, r, width in cells:
        xx = I(x) if width is None else I(x, x+width)
        rr = I(r) if width is None else I(r, r+width)
        value, radial = evaluate_source(ctx, xx, rr, audit=True)
        rows.append({'label': label, 'x': xx, 'r': rr, 'F': value, 'Fr': radial, 'F_sign': strict_sign(value), 'Fr_sign': strict_sign(radial)})
        print(json.dumps({'probe': label, 'F_display_only': [float(value.lo), float(value.hi)], 'Fr_display_only': [float(radial.lo), float(radial.hi)]}), flush=True)
    return rows


def cone_quadrature(ctx, n):
    need(n in (4, 8), 'DECLARED_CONE_BUDGET')
    I = ctx.base.I
    h, area_sum = T/n, F(0)
    accum = {'F': [F(0), F(0)], 'Fr': [F(0), F(0)]}
    caps = {'F': F(0), 'Fr': F(0)}
    rows = []
    for i in range(n):
        for k in range(i, 2*n-i):
            xl, xh = i*h, (i+1)*h
            rl, rh = R-T+k*h, R-T+(k+1)*h
            area = h*h/2 if k in (i, 2*n-i-1) else h*h
            need(0 <= xl < xh <= T and rl >= xh, 'EXTERIOR_CONE_CELL')
            value, radial = evaluate_source(ctx, I(xl, xh), I(rl, rh))
            for key, enclosure in (('F', value), ('Fr', radial)):
                accum[key][0] += area*enclosure.lo/2
                accum[key][1] += area*enclosure.hi/2
                caps[key] = max(caps[key], enclosure.mag())
            area_sum += area
            rows.append({'i': i, 'k': k, 'x': [xl, xh], 'r': [rl, rh], 'exact_intersection_area': area, 'F': value, 'Fr': radial})
        print(json.dumps({'cone_n': n, 'completed_rows': i+1, 'cells': len(rows)}), flush=True)
    need(len(rows) == n*(n+1) and area_sum == T*T, 'EXACT_CONE_COVER')
    return {
        'partition': {'n': n, 'h': h, 'cells': len(rows), 'exact_area': area_sum, 'cone': '|r-5/4|<=1/4-x', 'coverage': 'For i=0..N-1, k=i..2N-i-1. Boundary cells contribute exact triangles h²/2; other cells squares h². Disjoint interiors cover the entire cone.'},
        'cells': rows,
        'source_absolute_caps': caps,
        'free_apex_values': {key: I(*endpoints) for key, endpoints in accum.items()},
    }


def characteristic_quadrature(ctx, segments):
    need(segments in (8, 16), 'DECLARED_CHARACTERISTIC_BUDGET')
    I = ctx.base.I
    h = T/segments
    accum = {key: [F(0), F(0)] for key in ('F_plus', 'F_minus', 'Fr_plus', 'Fr_minus')}
    rows = []
    for j in range(segments):
        sl, sh = j*h, (j+1)*h
        row = {'segment': j, 'x': [sl, sh], 'weight': h/2}
        for branch, rlo, rhi in (('plus', R+T-sh, R+T-sl), ('minus', R-T+sl, R-T+sh)):
            need(rlo >= sh, 'EXTERIOR_CHARACTERISTIC_CELL')
            value, radial = evaluate_source(ctx, I(sl, sh), I(rlo, rhi))
            row[branch] = {'r': [rlo, rhi], 'F': value, 'Fr': radial}
            for key, enclosure in (('F_'+branch, value), ('Fr_'+branch, radial)):
                accum[key][0] += h*enclosure.lo/2
                accum[key][1] += h*enclosure.hi/2
        rows.append(row)
        print(json.dumps({'characteristic_segments': segments, 'completed_segments': j+1}), flush=True)
    integrals = {key: I(*endpoints) for key, endpoints in accum.items()}
    return {
        'partition': {'segments_per_characteristic': segments, 'h': h, 'characteristics': 'r_plus(s)=5/4+1/4-s; r_minus(s)=5/4-1/4+s; 0<=s<=1/4', 'weight_each_branch': h/2},
        'segments': rows,
        'half_integrals': integrals,
        'free_apex_derivatives': {
            'u_x': integrals['F_plus']+integrals['F_minus'],
            'u_r': integrals['F_plus']-integrals['F_minus'],
            'u_xr': integrals['Fr_plus']+integrals['Fr_minus'],
            'u_rr': integrals['Fr_plus']-integrals['Fr_minus'],
        },
    }


def campaign_run(ctx, n, segments, op, apex_op, apex_source):
    I = ctx.base.I
    cone = cone_quadrature(ctx, n)
    lines = characteristic_quadrature(ctx, segments)
    fcap, frcap = cone['source_absolute_caps']['F'], cone['source_absolute_caps']['Fr']
    q, qr = op['q'], op['q_r']
    mu = T*fcap/(1-q)
    mv = (T*frcap+qr*mu)/(1-q)
    e1, e2 = q*mu, q*mv+qr*mu
    e0, ev0 = T*e1/2, T*e2/2
    free = lines['free_apex_derivatives']
    jet = {'u': cone['free_apex_values']['F']+I(-e0, e0)}
    for key in ('u_x', 'u_r'):
        jet[key] = free[key]+I(-e1, e1)
    ur_line = jet['u_r']
    ur_area = cone['free_apex_values']['Fr']+I(-ev0, ev0)
    jet['u_r'] = intersection(I, ur_line, ur_area, 'APEX_UR_AREA_LINE_COMPATIBILITY')
    for key in ('u_xr', 'u_rr'):
        jet[key] = free[key]+I(-e2, e2)
    coeff = apex_op['coefficient_intervals']
    potential_forcing = coeff['A']*jet['u_x']+coeff['B']*jet['u_r']+coeff['D']*jet['u']+apex_source[0]
    jet['u_xx'] = jet['u_rr']+potential_forcing
    need(all(value.lo <= value.hi for value in jet.values()), 'FINITE_OUTWARD_TWOJET')
    # Ordinary coordinate derivatives of the registered frame coefficient.
    # These are neither covariant tensor derivatives nor curvature components.
    bg = apex_op['background_value_intervals']
    ell, L, a = bg['ell'], bg['L'], bg['a']
    natural = {
        'u': jet['u'],
        'u_t': jet['u_x']/ell,
        'u_y': -L*jet['u_r'],
        'u_tt': jet['u_xx']/ell**2,
        'u_ty': -L*jet['u_xr']/ell,
        'u_yy': L**2*jet['u_rr']+a*L*jet['u_r'],
    }
    signs = {key: strict_sign(jet[key]) for key in ('u_xx', 'u_xr', 'u_rr')}
    passes = any(sign != 'inconclusive' for sign in signs.values())
    display = {key: [float(value.lo), float(value.hi)] for key, value in jet.items()}
    result = {
        'cone_quadrature': cone,
        'characteristic_quadrature': lines,
        'resolvent_ledger': {'M_u': mu, 'M_v': mv, 'E_0': e0, 'E_1': e1, 'E_2': e2, 'E_v_value': ev0, 'proof': 'M_u=T Fcap/(1-q); M_v=(T Frcap+q_r M_u)/(1-q). E_1=q M_u; E_2=q M_v+q_r M_u; E_0=T E_1/2. All feedback and coefficient radial cross-couplings are summed through the resolvent.'},
        'apex_u_r_two_routes': {'characteristic_F': ur_line, 'cone_Fr': ur_area, 'intersection': jet['u_r']},
        'pde_closure': {'identity': 'u_xx=u_rr+A u_x+B u_r+D u+F', 'lower_order_plus_source_interval': potential_forcing, 'exact_equation_defines_u_xx': True, 'requires_source_x_derivative': False, 'independent_six_interval_residual_is_not_claimed_zero': True},
        'total_apex_twojet': jet,
        'ordinary_t_y_coefficient_twojet': natural,
        'ordinary_t_y_conversion': {'identities': 'u_t=ux/ell; u_y=-L ur; u_tt=uxx/ell²; u_ty=-L uxr/ell; u_yy=L² urr+a L ur; L=1/(ell rho), L_r=a', 'apex_actual_background': bg, 'covariant_tensor_derivatives': False, 'curvature_components': False, 'observable': False},
        'ordinary_t_y_second_derivative_signs': {key: strict_sign(natural[key]) for key in ('u_tt', 'u_ty', 'u_yy')},
        'scaled_value_and_first_derivative_signs': {key: strict_sign(jet[key]) for key in ('u', 'u_x', 'u_r')},
        'second_derivative_signs': signs,
        'passes_diagnostic_nonzero_gate': passes,
        'finite_twojet_primary_acceptance': True,
        'display_only': display,
    }
    print(json.dumps({'cone_n': n, 'segments': segments, 'total_twojet_display_only': display, 'second_derivative_signs': signs, 'passes_diagnostic_nonzero_gate': passes}), flush=True)
    return result


def old_value_comparison(repo, interval):
    path = repo/'research'/LABELS['F']/'tensor/TT_TOTAL_RESPONSE_NORMAL.json'
    old = json.loads(path.read_text())
    endpoints = old['partition_runs'][-1]['total_Tplus_response_interval']
    lo, hi = map(F, endpoints)
    overlap = [max(lo, interval.lo), min(hi, interval.hi)]
    need(overlap[0] <= overlap[1], 'OLD_TT_VALUE_INCOMPATIBILITY')
    return {'path': 'research/'+LABELS['F']+'/tensor/TT_TOTAL_RESPONSE_NORMAL.json', 'sha256': sha(path), 'old_total_value_interval': [lo, hi], 'new_total_value_interval': interval, 'nonempty_overlap': overlap, 'old_run_replayed': False}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo-root', type=Path, required=True)
    parser.add_argument('--input-module', type=Path, required=True)
    parser.add_argument('--input-module-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--probes-only', action='store_true')
    args = parser.parse_args()
    repo, module_path = args.repo_root.resolve(), args.input_module.resolve()
    guard_output(repo, args.output)
    own_pin = sha(Path(__file__))
    authenticated = authenticate(repo)
    need(not args.input_module.is_symlink() and module_path.is_file(), 'NEW_INPUT_MODULE_TYPE')
    need(sha(module_path) == args.input_module_sha256, 'NEW_INPUT_MODULE_EXTERNAL_PIN')
    source = load('authenticated_actual_TT_radial_source', module_path)
    ctx = source.Context(repo)
    sample_rows = probes(ctx)
    op = operator(ctx)
    apex_op = operator(ctx, ctx.base.I(T), ctx.base.I(R))
    apex_source = evaluate_source(ctx, ctx.base.I(T), ctx.base.I(R), audit=True)
    runs = []
    comparison = None
    if not args.probes_only:
        runs.append(campaign_run(ctx, 4, 8, op, apex_op, apex_source))
        if not runs[-1]['passes_diagnostic_nonzero_gate']:
            runs.append(campaign_run(ctx, 8, 16, op, apex_op, apex_source))
        comparison = old_value_comparison(repo, runs[-1]['total_apex_twojet']['u'])
    ctx.finish()
    need(authenticate(repo) == authenticated, 'FROZEN_CLOSURES_CHANGED')
    need(sha(module_path) == args.input_module_sha256, 'NEW_INPUT_MODULE_CHANGED')
    need(sha(Path(__file__)) == own_pin, 'EVALUATOR_CHANGED')
    status = 'PROBED_ONLY' if not runs else 'CERTIFIED_FINITE_ACTUAL_TOTAL_TT_TWOJET'
    diagnostic = 'not_run' if not runs else ('strict_nonzero_second_derivative_certified' if runs[-1]['passes_diagnostic_nonzero_gate'] else 'inconclusive_at_declared_maximum')
    certificate = ctx.certificate() if callable(ctx.certificate) else ctx.certificate
    result = {
        'schema': 'HDBLAST_BOUNDED_ACTUAL_TOTAL_TT_TWOJET_V1',
        'status': status,
        'diagnostic_nonzero_status': diagnostic,
        'evaluator_sha256': own_pin,
        'input_module_sha256': args.input_module_sha256,
        'authentication': authenticated,
        'all_nine_complete_closures_authenticated_before_project_imports_and_after_run': True,
        'new_input_module_authenticated_before_import_and_after_run': True,
        'budget_declaration': DECLARATION,
        'source_input_certificate': certificate,
        'domain': {'apex': [T, R], 'cone': '0<=x<=1/4, 1+x<=r<=3/2-x', 'shell_reached': False, 'reflection_reached': False, 'shell_lift': 'exactly0 by backward domain of dependence', 'TT_initial_value': 'exactly0', 'TT_initial_velocity': 'exactly0', 'TT_free_response': 'exactly0'},
        'probes': sample_rows,
        'whole_cone_operator': op,
        'apex_operator': apex_op,
        'apex_source': {'F': apex_source[0], 'Fr': apex_source[1]},
        'campaign_runs': runs,
        'old_TT_value_compatibility': comparison,
        'conventions': {'coefficient': 'tau_plus=tau_minus=u, multiplying Tplus+Tminus; not their sum and no half-source factor', 'derivatives': 'ordinary scaled x,r and converted ordinary t,y derivatives of registered radially parallel frame TT coefficient', 'mixed_partials': 'u_xr=u_rx by retained smooth actual solution', 'coordinate_conversion': 't=ell*x and partial_r=-ell*rho*partial_y; coefficient partial derivatives only, not covariant tensor two-jet, curvature or detector observable'},
        'limits': ['Conditional registered-gauge angular metric coefficient of the full total second response.', 'Actual complete TT source, its certified radial derivative, scalar feedback and coefficient derivative cross-couplings retained.', 'No numerical differentiation, source-as-exact-comparison substitution, old campaign replay, or frozen edits.', 'PDE compatibility is recorded through the exact closure identity; independent interval residual width is not an Einstein residual bound.', 'No curvature, gauge invariant force, nonlinear family, cosmology, origin, novelty or proper-motion claim.'],
    }
    with args.output.open('x') as handle:
        json.dump(encode(result), handle, indent=2, sort_keys=True)
        handle.write('\n')
    print(json.dumps({'status': status, 'diagnostic_nonzero_status': diagnostic, 'campaign_runs': len(runs)}), flush=True)


if __name__ == '__main__':
    main()
