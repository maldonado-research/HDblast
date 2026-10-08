#!/usr/bin/env python3
"""Independent exact arithmetic on completed JSON/GZ exports only.

Imports no production numerical module, source builder, or binary-array codec.
"""
from pathlib import Path
from fractions import Fraction as Q
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
import gzip
import hashlib
import json

W = Path('/workspace/hdblast-research-work/continuation-state-comparison-20261007')
OUT = W / 'root/actual-normal-001'
REVIEW = W / 'actual-independent-review'
REG = '05e9943f0b4c8134252a2fecef7631ddba8bd398d18553d6e126fbf50fc3aacc'
GO = '27784552d9d7e10167066927d54f1757351b0f925b499e951b79849a44b62769'
FREEZE = 'bedcca7e86995da1230c12a20fae3755b31f4a94'


def require(condition, detail):
    if not condition:
        raise ValueError(detail)


def q(obj):
    if isinstance(obj, str):
        return Q(obj)
    return Q(int(obj['numerator']), int(obj['denominator']))


def canon(value):
    return f'{value.numerator}/{value.denominator}'


def decimal_outward(value, upper):
    with localcontext() as ctx:
        ctx.prec = 10
        ctx.rounding = ROUND_CEILING if upper else ROUND_FLOOR
        result = str(Decimal(value.numerator) / Decimal(value.denominator))
    require(Q(result) >= value if upper else Q(result) <= value,
            'Decimal direction changed')
    return result


def pin(path):
    raw = path.read_bytes()
    return {'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()}


def new_acc():
    return {'count': 0, 'previous_k': Q(0), 'max_U': Q(0), 'max_W': Q(0),
            'max_c': Q(0), 'weighted_U': Q(0), 'weighted_W': Q(0),
            'mu_sum': Q(0), 'weight_sum': Q(0), 'witnesses': {}}


def main():
    data = json.loads((OUT / 'DATA.json').read_text())
    strict = json.loads((REVIEW / 'STRICT_READBACK_NORMAL_001.json').read_text())
    require(strict['status'] == 'PASS_STANDALONE_EXACT_ALL_NODE_READBACK',
            'Completed authenticated strict readback required')
    require(strict['nodes_verified'] == 49152 and strict['prefix_cases_verified'] == 12,
            'Wrong strict readback universe')
    epsilon, pi = q(data['represented_epsilon']), q(data['represented_Pi'])
    require(epsilon == Q(3777893186295716171,37778931862957161709568), 'epsilon')
    require(pi == Q(14488038916154245685,4611686018427387904), 'Pi')
    radii = {part: q(value) for part, value in data['max_complete_target_export_L1_radius'].items()}
    summaries = {r['case']: r for r in data['rows']}
    capsules = tuple(f'{s}/{g}' for s in ('positive_B','signed_uB') for g in ('coarse','fine'))
    prefixes = {'coarse': {64:2048,128:4096,256:8192},
                'fine': {64:4096,128:8192,256:16384}}
    observed = {capsule: new_acc() for capsule in capsules}
    cases = []
    scale = Q(1, 1 << 96)
    count = 0
    with gzip.open(OUT / 'NODE_TARGETS.jsonl.gz', 'rt', encoding='ascii') as handle:
        for line in handle:
            row = json.loads(line)
            capsule = row['capsule']
            require(capsule in observed, 'Capsule identity')
            acc = observed[capsule]
            require(row['index'] == acc['count'], 'Capsule node index')
            k, weight = q(row['k']), q(row['weight'])
            require(acc['previous_k'] < k < 256 and weight > 0, 'Momentum/weight domain')
            grid = capsule.rsplit('/',1)[1]
            for cutoff, n in prefixes[grid].items():
                require((acc['count'] < n and k < cutoff) or
                        (acc['count'] >= n and k > cutoff), 'Strict prefix mask')
            saved = {'U': tuple(q(t)/epsilon for t in row['stored_raw_u']),
                     'W': tuple(q(t)/epsilon for t in row['stored_raw_w'])}
            mu = weight*k*k/(2*pi*pi)
            upper = {}
            for part in ('U','W'):
                target = row['target_dyadic96'][part]
                total = Q(0)
                radius = Q(0)
                for index, axis in enumerate(('real','imag')):
                    lo, hi = target[axis]
                    require(type(lo) is int and type(hi) is int and lo <= hi, 'Rectangle')
                    total += max(abs(saved[part][index]-lo*scale),
                                 abs(saved[part][index]-hi*scale))
                    radius += (hi-lo)*scale/2
                require(radius <= radii[part] <= Q(1,10**18), 'Complete target radius')
                upper[part] = total
                if total > acc['max_'+part]:
                    acc['max_'+part] = total
                    acc['witnesses'][part] = {
                        'index': row['index'], 'k': canon(k),
                        'rectangle_L1_upper': canon(total),
                        'node_complete_target_L1_radius': canon(radius)}
                acc['weighted_'+part] += mu*total
            c = saved['U'][0]-saved['W'][1]/(2*k)
            if abs(c) > acc['max_c']:
                acc['max_c'] = abs(c)
                acc['witnesses']['c'] = {'index': row['index'], 'k': canon(k), 'signed_c': canon(c)}
            acc['mu_sum'] += mu
            acc['weight_sum'] += weight
            acc['count'] += 1
            acc['previous_k'] = k
            count += 1
            reverse = {n: cutoff for cutoff,n in prefixes[grid].items()}
            if acc['count'] not in reverse:
                continue
            cutoff = reverse[acc['count']]
            case = f'{capsule}/{cutoff}'
            summary = summaries[case]
            require(summary['node_count'] == acc['count'], 'Prefix count')
            require(q(summary['last_k']) == k, 'Last exact k')
            for part in ('U','W'):
                require(q(summary['suprema']['delta_'+part+'_L1_upper']) == acc['max_'+part], 'Max '+part)
                require(q(summary['sums']['delta_'+part+'_weighted_L1_upper']) == acc['weighted_'+part], 'Weighted '+part)
            require(q(summary['suprema']['c_abs']) == acc['max_c'], 'Exact c maximum')
            require(q(summary['suprema']['first_order_Wronskian_defect_abs']) == 2*epsilon*acc['max_c'], 'First-order Wronskian')
            require(q(summary['suprema']['kA_L1_upper']) == acc['max_W']/2, 'Stable kA maximum')
            require(q(summary['sums']['kA_weighted_L1_upper']) == acc['weighted_W']/2, 'Stable kA weighted')
            for field in ('mu_sum','weight_sum'):
                require(q(summary['sums'][field]) == acc[field], field)
            intervals = {}
            for part in ('U','W'):
                hi = acc['max_'+part]
                lo = max(Q(0), hi-2*radii[part])
                intervals[part] = {
                    'exact_lower': canon(lo), 'exact_upper': canon(hi),
                    'decimal_lower_outward': decimal_outward(lo,False),
                    'decimal_upper_outward': decimal_outward(hi,True),
                    'raw_state_exact_lower': canon(epsilon*lo),
                    'raw_state_exact_upper': canon(epsilon*hi)}
            bound_keys = ('R_uniform_abs_upper','P_uniform_abs_upper','work_abs_upper',
                          'int_pressure_abs_upper','integral_absolute_pressure_abs_upper')
            cases.append({'case': case, 'node_count': acc['count'],
                          'max_normalized_L1_error_intervals': intervals,
                          'exact_max_abs_c': canon(acc['max_c']),
                          'exact_max_first_order_Wronskian_defect': canon(2*epsilon*acc['max_c']),
                          'witnesses': dict(acc['witnesses']),
                          'finite_scaled_stress_upper_bounds': {
                              key: {'exact': canon(q(summary['sums'][key])),
                                    'decimal_upper_outward': decimal_outward(q(summary['sums'][key]),True)}
                              for key in bound_keys}})
    require(count == 49152 and len(cases) == 12, 'Complete all-node/prefix universe')
    for capsule in capsules:
        require(observed[capsule]['count'] == prefixes[capsule.rsplit('/',1)[1]][256], 'Full capsule count')

    entry = json.loads((OUT / 'ENTRY_RECEIPT.json').read_text())
    execution = json.loads((OUT / 'EXECUTION.json').read_text())
    require(entry['status'] == 'PASS_AUTHENTICATED_STORED_COMPARISON', 'Entry status')
    require(entry['registration_sha256'] == REG and entry['public_go_sha256'] == GO and
            entry['freeze_commit'] == FREEZE and entry['fabricated_only'] is False, 'Entry provenance')
    require(entry['counter'] == data['counter'], 'Counter agreement')
    for name, expected in entry['files'].items():
        require(pin(OUT / name) == expected, 'Payload pin '+name)
    require(pin(OUT / 'ENTRY_RECEIPT.json') == execution['entry_receipt'], 'Custodian entry pin')
    require(execution['status'] == 'PASS_BOUNDED_REGISTERED_EXECUTION' and
            execution['exit_code'] == execution['wrapper_exit_code'] == execution['wait_status'] == 0 and
            execution['stop_reason'] is None and execution['monitor_error'] is None and
            execution['fabricated_only'] is False and execution['optimized'] is False, 'Custodian successful actual normal run')
    command = execution['command']
    require('-O' not in command and command[command.index('--registration-sha256')+1] == REG and
            command[command.index('--public-go-sha256')+1] == GO, 'Custodian command binding')
    require(pin(OUT / 'child.log')['sha256'] == execution['captured_child_log_sha256'] ==
            execution['observed_child_log_sha256'], 'Child log identity')

    result = {
        'status': 'PASS_INDEPENDENT_SERIALIZED_NORMS_AND_ACTUAL_RECEIPTS',
        'registration_sha256': REG, 'public_go_sha256': GO, 'freeze_commit': FREEZE,
        'scope': 'SERIALIZED_EXPORTS_ONLY; physical target truth retains registered analytic-proof/worker dependency',
        'nodes_independently_recomputed': count, 'prefix_cases_independently_recomputed': len(cases),
        'reviewer_original_array_reads': 0, 'reviewer_original_array_decodes': 0,
        'reviewer_physical_source_callbacks': 0, 'reviewer_target_evaluations': 0,
        'complete_target_L1_radii': {p: canon(r) for p,r in radii.items()},
        'lower_bound_rule': 'max(0,S_X-2R_X), in Cartesian complex L1 norm, with global complete target radius R_X',
        'cases': cases,
        'normal_output_pins': {name: pin(OUT / name) for name in sorted(entry['files'])+['ENTRY_RECEIPT.json','EXECUTION.json','child.log']},
        'first_order_defect_scope': '2 epsilon c; full finite-epsilon quadratic Wronskian term excluded',
        'finite_stress_units': 'inherited a_0^4/epsilon scaled linear stress with exact retained dk weights and represented Pi',
        'limitations': {'full_twelve_case_continuous_pressure_contacts': 'UNRESOLVED',
                        'later_saved_numerical_trajectory_error': 'NOT_CERTIFIED',
                        'metric_calibration': 'FAIL_UNCHANGED',
                        'higher_dimensional_origin': 'NOT_ESTABLISHED',
                        'external_novelty': 'NOT_ASSESSED'}}
    destination = REVIEW / 'INDEPENDENT_SERIALIZED_NORMS_NORMAL_001.json'
    destination.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status': result['status'], 'nodes': count, 'cases': len(cases),
                      'receipt': str(destination), **pin(destination)},sort_keys=True))


if __name__ == '__main__':
    main()
