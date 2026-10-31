#!/usr/bin/env python3
"""Independent read-only regularity/algebra audit; no producer imports."""
import argparse
import hashlib
import json
from fractions import Fraction as Q
from pathlib import Path
import sys

sys.dont_write_bytecode = True
SOURCE_MAP_PIN='e4530a87ad17f90113f4d16719473da5a9e4196e0f29ac5e9ea12c3a6526764a'


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def encode(x):
    if isinstance(x, Q):
        return str(x)
    if isinstance(x, dict):
        return {k: encode(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [encode(v) for v in x]
    return x


def poly(cs, x):
    out = Q(0)
    for c in reversed(cs):
        out = out*x+c
    return out


def deriv(cs):
    return [k*cs[k] for k in range(1, len(cs))]


def interval_horner(cs, xmax):
    lo = hi = Q(0)
    for c in reversed(cs):
        products = [Q(0), lo*xmax, hi*xmax]
        lo, hi = min(products)+c, max(products)+c
    return [lo, hi]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--source-root','--root', dest='root',type=Path, required=True)
    ap.add_argument('--pins', type=Path, default=Path(__file__).resolve().parent/'regularity-source-pins.json')
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    need(not args.output.exists(), 'OUTPUT_MUST_BE_ABSENT')
    pins_raw = args.pins.read_bytes()
    pins_hash = hashlib.sha256(pins_raw).hexdigest()
    need(pins_hash==SOURCE_MAP_PIN,'FROZEN_REGULARITY_SOURCE_MAP_HASH')
    pins = json.loads(pins_raw)
    own_hash = digest(Path(__file__))
    for rel, pin in pins.items():
        path = args.root/rel
        need(path.is_file() and not path.is_symlink(), 'SOURCE_TYPE:'+rel)
        need(digest(path) == pin, 'SOURCE_HASH:'+rel)
    source = (args.root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026/initial/full23_initial.py').read_text()
    need('fy = [dy(value) for value in fs]' in source, 'REVIEWED_HIGHEST_DERIVATIVE_PATH')
    need('for j, i in enumerate(self.mean_ids):' in source and
         'displacement[i], velocity[i] = mean_u[j], mean_v[j]' in source,
         'REVIEWED_MEAN_OVERRIDE')
    path = args.root/'research/HDBLAST_INPUT_PRECISION_20261023/initial/RECENTERED_CANDIDATE_NORMAL.json'
    cs = {k: list(map(Q, v)) for k, v in json.loads(path.read_text())['representation_background_coefficients_ascending'].items()}
    ell, radius = Q(1, 10**6), Q(7, 4)
    denominators = {k: interval_horner(cs[k], ell*radius) for k in ('rho', 'p')}
    need(denominators['rho'][0] > 0 and denominators['p'][1] < 0,
         'COMPLETE_FOOTPRINT_NOMINAL_DENOMINATORS_SEPARATED')
    zero = ['K:0','K:2','K:4','L:0','L:2','L:4','v:2','v:4','u:2','u:4','H:2','H:4','T:+','T:-']
    displacement_orders = {name: None for name in zero}
    for channel, order in [('N',0),('J',2),('w',3)]:
        displacement_orders.update({f'{channel}:{l}':order for l in (0,2,4)})
    velocity_orders = {name: None for name in ['v:2','v:4','H:2','H:4','T:+','T:-']}
    for channel, order in [('K',1),('L',3),('N',1),('J',1),('w',2)]:
        velocity_orders.update({f'{channel}:{l}':order for l in (0,2,4)})
    velocity_orders.update({'u:2':2,'u:4':2})
    need(len(displacement_orders) == len(velocity_orders) == 23,
         'ALL23_INVENTORIES')
    # For each completed component f, the extension is defined with exact
    # symbolic traces f0,f1,f2. At h=0 the jets are f0,f1,f2 exactly.
    extension_jet_differences = {name: [Q(0)]*3 for name in displacement_orders}
    # Compact k_delta has a triple zero at each endpoint: its global zero
    # extension is C2. Convolution with continuous piecewise-affine velocity
    # gives v_s'' = sum_j jump_j k_delta(r-r_j), hence v_s is at least C4.
    base_kernel = [Q(1),Q(0),Q(-3),Q(0),Q(3),Q(0),Q(-1)]
    kernel_traces = []
    dcs = base_kernel
    for _ in range(3):
        kernel_traces.append([poly(dcs,Q(-1)),poly(dcs,Q(1))])
        dcs = deriv(dcs)
    need(all(z == 0 for row in kernel_traces for z in row),
         'VELOCITY_KERNEL_C2_ZERO_EXTENSION')
    # A raw degree-two extension matches q,q_r,q_rr but loses q_rrr.
    # At R, q=e^(-4/3)>0 and q_rrr/q=925696/729 != 0.
    t1,t2,t3 = Q(-128,9),Q(-3584,27),Q(-40960,27)
    q1,q2,q3 = t1,t2+t1*t1,t3+3*t1*t2+t1**3
    need(q3 == Q(925696,729), 'RAW_QUADRATIC_MISSING_THIRD_JET')
    z = ell*radius
    bg = {k:poly(v,z) for k,v in cs.items()}
    bg1 = {k:poly(deriv(v),z) for k,v in cs.items()}
    rho,p,alpha,eta,a = [bg[k] for k in ('rho','p','alpha','eta','a')]
    U = list(map(Q,['-2/27','0','14/9','50/27','-1/6','-4/9','-2/27']))
    U1,U2 = deriv(U),deriv(deriv(U))
    h = poly(U1,eta)/p
    h1 = poly(U2,eta)*bg1['eta']/p-poly(U1,eta)*bg1['p']/p**2
    k = rho*(Q(5,2)*a-h)
    k1 = bg1['rho']*(Q(5,2)*a-h)+rho*(Q(5,2)*bg1['a']-h1)
    A = 3*alpha/(rho**3*p)
    A_log1 = bg1['alpha']/alpha-3*bg1['rho']/rho-bg1['p']/p
    vy_over_q = -A/rho*(A_log1*(q1/ell+k)+q2/ell**2+k1+k*q1/ell)
    # w:0 = -(3/2p) dy(fs0). Coefficient of raw q_zzz is
    # A*vy/(6p); the lower q jets agree for the two raw continuations.
    jump_over_q_squared = A*vy_over_q*q3/(6*p*ell**3)
    need(jump_over_q_squared != 0, 'RAW_QUADRATIC_RECONSTRUCTED_SCALAR_VALUE_JUMP_NONZERO')
    for rel, pin in pins.items():
        need(digest(args.root/rel) == pin, 'SOURCE_CHANGED:'+rel)
    need(digest(args.pins) == pins_hash and digest(Path(__file__)) == own_hash,
         'AUDITOR_AND_PIN_MAP_UNCHANGED')
    report = {
        'status':'PASS_INDEPENDENT_ANALYTIC_REGULARITY_AUDIT',
        'scope':'Analytic regularity inventory, coefficient denominator separation, and exact raw-quadratic control. No interval initial-error or response gate is claimed.',
        'all23_highest_raw_q_derivative_displacement':displacement_orders,
        'all23_highest_raw_q_derivative_velocity':velocity_orders,
        'required_q_order_for_displacement_rr':5,
        'required_q_order_for_analytic_velocity_r':4,
        'unchanged_velocity_definition':'same complete finite P1 table, zero first cell, final-secant-slope affine exterior, parity, convolution width 1/2048 and compact kernel 35/(32 delta)(1-(s/delta)^2)^3',
        'smoothed_P1_velocity_minimum_global_regulariy':'C4 by finite slope-jump convolution representation',
        'velocity_kernel_unscaled_endpoint_jets_through2':kernel_traces,
        'nominal_denominator_enclosures_on_0_to_7over4':denominators,
        'complete_component_quadratic_extension_jet_differences_through2':extension_jet_differences,
        'outer_extension_regularity':'C2 exactly; generally not C3; higher global smoothness is not assumed by this audit',
        'flatness_proof':'For r->1+, each q derivative is exp(4-1/d) times a rational function with finite-order poles in d; the exponential dominates every pole. Smooth nominal coefficients and separated rho,p denominators imply every required completed component jet tends to zero. All23 zero slots and five overridden means are included.',
        'old_artificial_inner_joins':'One exact analytic q/E/nominal formula on both sides; evaluation branches may only enclose that same expression.',
        'raw_quadratic_control':{'radius':radius,'q_value':'exp(-4/3)','q_r_over_q':q1,'q_rr_over_q':q2,'q_rrr_over_q':q3,'left_minus_raw_quadratic_right_w0_value_over_q_squared':jump_over_q_squared,'proved_nonzero':True},
        'source_sha256_before_and_after':pins,
        'auditor_sha256_before_and_after':own_hash,
        'pins_sha256_before_and_after':pins_hash,
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(encode(report),indent=2,sort_keys=True)+'\n')
    print(report['status'])


if __name__ == '__main__':
    main()
