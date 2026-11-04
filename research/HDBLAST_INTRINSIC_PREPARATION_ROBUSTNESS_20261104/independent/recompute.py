#!/usr/bin/env python3
"""Independent Fraction proof from two old byte-pinned JSON artifacts.

No producer imports, new oracle calls, sampling, or scientific dependencies.
"""
from fractions import Fraction as F
from hashlib import sha256
from math import factorial
from pathlib import Path
import argparse
import copy
import json
import sys

PINS = {
    'static': ('research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/inputs/ACTUAL_ANNULUS_STATIC_JETS.json',
               '3aa2c8d2e9ad652b72a34b4d8a6e6a6f14f72c32b8dcd1fbdc67556093265c22'),
    'response': ('research/HDBLAST_PREPARED_SHELL_VOLUME_20261101/response/references/FINAL_NORMAL.json',
                 '1a9ebffb6bee9332bc0fd0df0f5d1d9f7f2a796825c2a563e4c644974ba4d8a3'),
}


class Rejected(Exception):
    pass


def need(test, why):
    if not test:
        raise Rejected(why)


class Parser(argparse.ArgumentParser):
    def error(self, message):
        raise Rejected('invalid CLI arguments')


def read_sources(root, candidate=None):
    result = {}
    for key, (relative, expected) in PINS.items():
        source = candidate if key == 'response' and candidate is not None else root / relative
        payload = source.read_bytes()
        need(sha256(payload).hexdigest() == expected, 'source hash: ' + key)
        result[key] = json.loads(payload)
    return result


def interval(node):
    lo, hi = map(F, node['exact'])
    need(lo <= hi, 'reversed interval')
    return lo, hi


def contained(node, lo, hi, name):
    a, b = interval(node)
    need(lo <= a <= b <= hi, 'cap: ' + name)


def record(x):
    if isinstance(x, tuple):
        return {'exact': list(map(str, x)), 'display_only': list(map(float, x))}
    return {'exact': str(x), 'display_only': float(x)}


def compute(sources, n_claim=F(1, 2000)):
    st, old = sources['static'], sources['response']
    ell, T, delta = F(1, 10**6), F(5, 4), F(1, 4)
    rho_lo, rho_hi = F(129), F(131)
    amin, amax = F(111, 1000), F(14, 125)
    pmin, pmax = F(3, 25000), F(1, 7000)
    pymax, aymax = F(1, 4000), F(1, 10000)
    contained(st['y_jets']['rho'][0], rho_lo, rho_hi, 'rho')
    contained(st['y_jets']['a'][0], amin, amax, 'a')
    contained(st['y_jets']['a'][1], -aymax, aymax, 'a_y')
    contained(st['y_jets']['p'][0], -pmax, -pmin, 'p')
    contained(st['y_jets']['p'][1], -pymax, pymax, 'p_y')
    contained(old['first']['rho'], rho_lo, rho_hi, 'rho_b')
    contained(old['first']['Wprime'], -F(1, 1000), -F(1, 1200), 'Fprime')
    need(F(st['domain']['t']) == 0 and F(st['domain']['ell']) == ell, 'initial scale')
    need(F(st['domain']['z_min']) <= ell and F(st['domain']['z_max']) >= ell*T,
         'causal initial segment coverage')
    need(st['domain']['coverage'] == 'Every point of the original packet annulus; no samples.',
         'continuous domain declaration')
    need(st['normalization'] == 'coefficient j = D_z^j(value)/j! or D_y^j(value)/j!, as named',
         'first y-jet equals actual first derivative')
    q = F(1201, 3200000)
    need(F(old['q']['exact']) == q, 'contraction constant')
    b = ell*rho_hi*amax
    g = ell*rho_hi*pmax*pmax/(3*amin)
    # Max absolute row sum of [[3b,-1],[3ell^2,-b+g]].
    mu = 1+3*b
    need(3*ell*ell+b+g <= mu, 'continuous matrix norm')
    taylor_x = mu*delta
    need(0 < taylor_x < 8, 'exponential Taylor ratio')
    exp_upper = sum((taylor_x**j/factorial(j) for j in range(7)), F())
    exp_upper += taylor_x**7/factorial(7)/(1-taylor_x/8)
    need(exp_upper < F(4, 3), 'continuous Gronwall exponential')
    # |D|<=R, zero data at r=1: max(|z|,|h|)<=H R.
    H = F(2, 3)*g*delta*F(4, 3)
    Hr = 3*ell*ell*H+(b+g)*H+F(2, 3)*g
    beta = pmax/(3*amin)
    beta_r = ell*rho_hi*(pymax/(3*amin)+pmax*aymax/(3*amin*amin))
    S = beta*(F(3, 2)*H+1)
    Sr = beta_r*(F(3, 2)*H+1)+beta*(F(3, 2)*Hr+1)
    # V_L=h_r+b h, and H=diag(1,2,1,3,1).
    V = Hr+b*H
    D0, Dr, Dv = 2*H+S, 2*Hr+Sr, 2*V
    N = max((D0+T*Dv)/3, Dr+Dv)
    need(N <= n_claim, 'full five-channel N/R guard')
    Cv = (rho_hi*ell)**2*F(1, 2000)*delta**3/(6*(1-q))
    Ct = rho_hi*ell*F(1, 2000)*delta**2/(2*(1-q))
    radius = F(2*10**33)
    N_safe = radius*n_claim
    need(N_safe <= 10**30, 'complete added preparation norm')
    bv, bt = interval(old['proper_v2']), interval(old['proper_Theta2'])
    full_v = (bv[0]-Cv*N_safe, bv[1]+Cv*N_safe)
    full_t = (bt[0]-Ct*N_safe, bt[1]+Ct*N_safe)
    need(full_v[1] < 0 and full_t[1] < 0, 'both coefficient signs')
    return {
        'status': 'PASS_INDEPENDENT_GEOMETRIC_CLASS_BOUND',
        'method': 'Fresh exact Fraction arithmetic; two old pinned inputs; no producer import or numerical oracle.',
        'source_pins': {relative: digest for relative, digest in PINS.values()},
        'segment': ['1', '5/4'],
        'matrix': {'b_max': record(b), 'g_max': record(g), 'mu': record(mu),
                   'Taylor_degree': 6, 'exp_upper': record(exp_upper), 'exp_simple_cap': '4/3'},
        'profile_caps_per_R': {key: record(value) for key, value in
             [('h', H), ('z', H), ('h_r', Hr), ('beta', beta), ('beta_r', beta_r),
              ('s', S), ('s_r', Sr), ('V_L', V), ('D0', D0), ('Dr', Dr), ('Dv', Dv), ('N', N)]},
        'N_per_R_guard': str(n_claim), 'class_radius': str(radius), 'N_safe': str(N_safe),
        'homogeneous_operator': {'Cv': record(Cv), 'Ct': record(Ct), 'q': record(q)},
        'base_coefficients': {'v2': record(bv), 'Theta2': record(bt)},
        'class_coefficients': {'v2': record(full_v), 'Theta2': record(full_t)},
        'new_oracle_calls': 0,
        'scope': 'Conditional coefficient class on the restricted compact metric-mean freedom with other canonical data fixed. No physical class selection or uniform nonlinear epsilon radius.'
    }


def controls(sources):
    result = []
    for name, sector, branch, key, value in [
        ('p_contains_zero', 'static', 'y_jets', 'p', {'exact': ['-1/7000', '0']}),
        ('a_contains_zero', 'static', 'y_jets', 'a', {'exact': ['0', '14/125']}),
        ('p_y_exceeds_cap', 'static', 'y_jets', 'p_y', {'exact': ['-1', '1']}),
    ]:
        trial = copy.deepcopy(sources)
        if key == 'p_y':
            trial[sector][branch]['p'][1] = value
        else:
            trial[sector][branch][key][0] = value
        try:
            compute(trial)
        except Rejected as error:
            result.append({'mutant': name, 'status': 'REJECTED', 'reason': str(error)})
        else:
            raise Rejected('negative control accepted: '+name)
    for name, edit in [('q_zero', lambda x: x['response']['q'].update(exact='0')),
                       ('Fprime_exceeds_cap', lambda x: x['response']['first'].update(Wprime={'exact':['-1','0']}))]:
        trial = copy.deepcopy(sources)
        edit(trial)
        try:
            compute(trial)
        except Rejected as error:
            result.append({'mutant': name, 'status': 'REJECTED', 'reason': str(error)})
        else:
            raise Rejected('negative control accepted: '+name)
    try:
        compute(sources, F(1, 2500))
    except Rejected as error:
        result.append({'mutant': 'false_stronger_N_guard_1over2500', 'status': 'REJECTED', 'reason': str(error)})
    else:
        raise Rejected('false N guard accepted')
    return {'in_process_mutant_count': len(result), 'checks': result,
            'note': 'Parsed-source mutants are evaluated after the independent byte-pin read, to reach semantic checks.'}


def main():
    need(sys.flags.isolated and sys.flags.no_site and sys.dont_write_bytecode, 'required interpreter flags: -I -S -B')
    parser = Parser(add_help=False, allow_abbrev=False)
    parser.add_argument('--repo-root', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--controls', action='store_true')
    argv = sys.argv[1:]
    for flag in ['--repo-root','--candidate','--output','--controls']:
        need(sum(arg == flag or arg.startswith(flag+'=') for arg in argv) <= 1, 'duplicate CLI argument')
    args = parser.parse_args(argv)
    for path, label in [(args.repo_root, 'repo-root'), (args.candidate, 'candidate'), (args.output, 'output')]:
        need(path.is_absolute() and str(path) == str(path.resolve()), 'absolute normalized non-symlink '+label+' required')
    need(args.repo_root.is_dir() and args.candidate.is_file(), 'source paths must exist')
    need(args.output.parent.is_dir(), 'output parent missing')
    need(not args.output.exists() and not args.output.is_symlink(), 'output already exists')
    sources = read_sources(args.repo_root, args.candidate)
    report = compute(sources)
    if args.controls:
        report['negative_controls'] = controls(sources)
    payload = (json.dumps(report, indent=2, sort_keys=True)+'\n').encode()
    with args.output.open('xb') as stream:
        stream.write(payload)
    print(json.dumps({'status': report['status'], 'output_sha256': sha256(payload).hexdigest()}))


if __name__ == '__main__':
    try:
        main()
    except (Rejected, OSError, ValueError, KeyError, TypeError, ZeroDivisionError) as error:
        print(json.dumps({'status': 'REJECTED', 'reason': str(error)}))
        sys.exit(2)
