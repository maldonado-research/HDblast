#!/usr/bin/env python3
"""Portable independent scientific comparator; explicit pinned local bootstrap."""
from fractions import Fraction as F
from hashlib import sha256
from math import factorial
from pathlib import Path
import argparse
import copy
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
LOCAL_SHA = '747ba7cbb2cab70911e0db1d23beef416b53611b917ee214f4a5f1d1a0f6140b'
BASELINE_SHA = 'e89d21d78395601cbf2af98b1fc127c10d1c6768a5cca489353d8774bcfd405a'
PROVENANCE_SHA = 'a220391f72e8eaf98cbe15cd319f320556fc24ff8d70bbb255669b50ccc0f6c3'
SOURCE_PINS = {
 'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/inputs/ACTUAL_ANNULUS_STATIC_JETS.json':'3aa2c8d2e9ad652b72a34b4d8a6e6a6f14f72c32b8dcd1fbdc67556093265c22',
 'research/HDBLAST_PREPARED_SHELL_VOLUME_20261101/response/references/FINAL_NORMAL.json':'1a9ebffb6bee9332bc0fd0df0f5d1d9f7f2a796825c2a563e4c644974ba4d8a3',
}

class BootstrapRejected(Exception): pass

def require(test, reason):
    if not test: raise BootstrapRejected(reason)

class Parser(argparse.ArgumentParser):
    def error(self, message): raise BootstrapRejected('invalid CLI arguments')

def pinned(path, digest, label):
    data = path.read_bytes()
    require(sha256(data).hexdigest() == digest, 'byte pin: '+label)
    return data

def bootstrap(repo):
    # Old source and local-input pins precede execution of any project module.
    for relative,digest in SOURCE_PINS.items(): pinned(repo/relative,digest,relative)
    baseline=json.loads(pinned(HERE/'inputs/FROZEN_SCIENTIFIC_BASELINE.json',BASELINE_SHA,'scientific baseline'))
    provenance=json.loads(pinned(HERE/'inputs/PORTABLE_PRODUCER_PROVENANCE.json',PROVENANCE_SHA,'portable provenance'))
    path=HERE/'recompute.py'
    pinned(path,LOCAL_SHA,'local independent recompute source')
    spec=importlib.util.spec_from_file_location('_hdblast_independent_recompute',path)
    require(spec is not None and spec.loader is not None,'local module spec')
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module,baseline,provenance

def exact(report, *path):
    node = report
    for key in path:
        node = node[key]
    return F(node['exact'])


def inspect(candidate, own, baseline=None, provenance=None):
    done = []
    def eq(path, expected):
        need(exact(candidate, *path) == expected, '/'.join(path))
        done.append('/'.join(path))
    own_caps = own['profile_caps_per_R']
    for key in ['D0', 'Dr', 'Dv', 'N']:
        eq(['Gram_norm_caps_per_R', key], F(own_caps[key]['exact']))
    for key in ['h', 'h_r', 's', 's_r', 'V_L']:
        eq(['conditioned_component_caps_per_R', key], F(own_caps[key]['exact']))
    eq(['conditioned_component_caps_per_R', 'k'], F(3, 2)*F(own_caps['beta']['exact']))
    eq(['conditioned_component_caps_per_R', 'abs_k_r'], F(3, 2)*F(own_caps['beta_r']['exact']))
    for label, own_key in [('H_per_R','h'),('Hr_per_R','h_r'),('Z_per_R','z')]:
        eq(['ode_caps',label], F(own_caps[own_key]['exact']))
    for label, own_key in [('bmax','b_max'),('gmax','g_max')]:
        eq(['ode_caps',label], F(own['matrix'][own_key]['exact']))
    mu = exact(candidate, 'ode_caps', 'matrix_inf_norm_guard')
    row1 = exact(candidate, 'ode_caps', 'matrix_row_1')
    row2 = exact(candidate, 'ode_caps', 'matrix_row_2')
    b, g, ell = F(own['matrix']['b_max']['exact']), F(own['matrix']['g_max']['exact']), F(1,10**6)
    need(row1 == 1+3*b and row2 == 3*ell*ell+b+g and mu >= max(row1,row2), 'matrix bound')
    done.append('continuous matrix max row sums')
    exponential = candidate['exponential_certificate']
    need(exponential['Taylor_degree'] == 6, 'Taylor degree')
    x = mu/4
    need(F(exponential['argument']['exact']) == x, 'Taylor argument')
    partial = sum((x**j/factorial(j) for j in range(7)), F())
    tail = x**7/factorial(7)
    ratio = x/8
    upper = partial+tail/(1-ratio)
    for key, value in [('partial_sum',partial),('first_tail_term',tail),('tail_ratio_upper',ratio),('exponential_upper',upper)]:
        need(F(exponential[key]['exact']) == value, 'Taylor '+key)
    need(upper < F(4,3), 'strict exponential guard')
    done.append('independent positive Taylor tail reconstruction')
    for key, name in [('Cv','volume'),('Ct','expansion')]:
        eq(['observable_operator_caps_per_N',name], F(own['homogeneous_operator'][key]['exact']))
    eq(['Gram_norm_caps_per_R','fixed_N_per_R_guard'], F(own['N_per_R_guard']))
    eq(['class','radius'], F(own['class_radius']))
    eq(['class_N','fixed_guard'], F(own['N_safe']))
    eq(['class_N','actual_derived_cap'], F(own['class_radius'])*F(own_caps['N']['exact']))
    need(candidate['class']['complete_increment_required'] is True, 'complete increment')
    need(candidate['class']['causal_initial_segment'] == own['segment'], 'causal segment')
    done.append('complete increment and causal scope')
    for key, name in [('Cv','volume'),('Ct','expansion')]:
        eq(['absolute_observable_difference_caps',name], F(own['N_safe'])*F(own['homogeneous_operator'][key]['exact']))
    for key in ['v2','Theta2']:
        need(interval(candidate['base_exact_intervals'][key]) == interval(own['base_coefficients'][key]), 'base '+key)
        need(interval(candidate['whole_class_exact_intervals'][key]) == interval(own['class_coefficients'][key]), 'whole class '+key)
        lo,hi = interval(candidate['whole_class_exact_intervals'][key])
        a,b = interval(candidate['frozen_outward_intervals'][key])
        need(a <= lo <= hi <= b < 0, 'outward final '+key)
        done.append('base and whole-class exact '+key)
    threshold = min(-interval(own['base_coefficients']['v2'])[1]/F(own['homogeneous_operator']['Cv']['exact']),
                    -interval(own['base_coefficients']['Theta2'])[1]/F(own['homogeneous_operator']['Ct']['exact']))
    eq(['old_exact_N_threshold','both'], threshold)
    need(F(own['N_safe']) < threshold, 'strict old threshold')
    if baseline is not None:
        excluded = {'authenticated_inputs','frozen_design','provenance','scope'}
        need(set(candidate) == set(baseline) | excluded, 'exact candidate top-level schema')
        projected = {key:value for key,value in candidate.items() if key not in excluded}
        need(projected == baseline, 'complete frozen scientific projection including all gates')
        need({key:candidate[key] for key in excluded} == provenance, 'exact declared portable provenance')
    return done


def mutants(candidate, own):
    outcomes = []
    for label,path in [
        ('omit_scalar_coefficient_derivative',['conditioned_component_caps_per_R','abs_k_r']),
        ('omit_mixed_metric_velocity',['conditioned_component_caps_per_R','V_L']),
        ('drop_angular_Gram_weight',['Gram_norm_caps_per_R','D0']),
        ('drop_radial_derivative',['Gram_norm_caps_per_R','Dr']),
        ('extra_half_in_volume_transfer',['observable_operator_caps_per_N','volume']),
        ('extra_half_in_expansion_transfer',['observable_operator_caps_per_N','expansion']),
        ('understate_matrix_row',['ode_caps','matrix_row_1']),
        ('understate_exponential_tail',['exponential_certificate','first_tail_term']),
        ('unscaled_direction_norm',['class_N','fixed_guard'])]:
        trial = copy.deepcopy(candidate)
        value = trial
        for key in path:
            value = value[key]
        value['exact'] = str(F(value['exact'])/2)
        try:
            inspect(trial, own)
        except Rejected as error:
            outcomes.append({'label':label,'status':'REJECTED','reason':str(error)})
        else:
            raise Rejected('mutant accepted: '+label)
    return outcomes


def main():
    require(sys.flags.isolated and sys.flags.no_site and sys.dont_write_bytecode,'required interpreter flags: -I -S -B')
    parser=Parser(add_help=False,allow_abbrev=False)
    parser.add_argument('--repo-root',type=Path,required=True)
    parser.add_argument('--candidate',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    argv=sys.argv[1:]
    for flag in ['--repo-root','--candidate','--output']:
        require(sum(arg == flag or arg.startswith(flag+'=') for arg in argv) <= 1,'duplicate CLI argument')
    args=parser.parse_args(argv)
    for path,label in [(args.repo_root,'repo-root'),(args.candidate,'candidate'),(args.output,'output')]:
        require(path.is_absolute() and str(path)==str(path.resolve()),'absolute normalized non-symlink '+label+' required')
    require(args.repo_root.is_dir() and args.candidate.is_file(),'source paths must exist')
    require(args.output.parent.is_dir(),'output parent missing')
    require(not args.output.exists() and not args.output.is_symlink(),'output already exists')
    core,baseline,provenance=bootstrap(args.repo_root)
    globals().update(need=core.need,interval=core.interval,Rejected=core.Rejected)
    own=core.compute(core.read_sources(args.repo_root))
    payload=args.candidate.read_bytes()
    candidate=json.loads(payload)
    checks=inspect(candidate,own,baseline,provenance)
    outcomes=mutants(candidate,own)
    report={
      'status':'PASS_INDEPENDENT_CANDIDATE_COMPARISON',
      'candidate_logical_role':'portable bounds producer report',
      'candidate_sha256':sha256(payload).hexdigest(),
      'method':'Independent two-source Fraction reconstruction; explicit SHA-pinned local module; no producer imports.',
      'accepted_check_count':len(checks),'accepted_checks':checks,
      'full_scientific_projection_matches_frozen_baseline':True,
      'exact_portable_provenance_schema_checked':True,
      'scientific_projection_excluded_top_level_keys':['authenticated_inputs','frozen_design','provenance','scope'],
      'in_process_report_mutants':len(outcomes),'report_mutant_outcomes':outcomes,
      'actual_CLI_negative_controls':0,
      'count_scope':'Nine in-process candidate mutants in this suite; actual CLI controls have a separate validation receipt.',
      'new_oracle_calls':0}
    out=(json.dumps(report,indent=2,sort_keys=True)+'\n').encode()
    with args.output.open('xb') as stream: stream.write(out)
    print(json.dumps({'status':report['status'],'output_sha256':sha256(out).hexdigest(),'accepted_check_count':len(checks),'report_mutants_rejected':len(outcomes)},sort_keys=True))

if __name__=='__main__':
    try: main()
    except Exception as error:
        print(json.dumps({'status':'REJECTED','reason':str(error)},sort_keys=True))
        sys.exit(2)
