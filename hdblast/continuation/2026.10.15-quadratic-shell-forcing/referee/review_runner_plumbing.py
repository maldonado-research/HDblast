"""Read-only helper and deterministic embedding review; no scientific execution."""
import hashlib
import json
from pathlib import Path
import pprint
import shutil
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent
PUB = ROOT / 'publication'


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def load(path, name):
    namespace = {'__name__': name, '__file__': str(path)}
    exec(compile(path.read_bytes(), str(path), 'exec'), namespace)
    return namespace


template = load(PUB / 'runner_template.py', 'referee_runner_template')
builder = load(PUB / 'build_replay_runner.py', 'referee_runner_builder')
checks = []


def passed(name):
    checks.append(name)


def rejects(name, action, expected):
    try:
        action()
    except RuntimeError as error:
        require(str(error) == expected, 'Unexpected failure for ' + name)
    else:
        raise RuntimeError('Required rejection missing: ' + name)
    passed(name)


scratch = HERE / '.runner-plumbing-focused-fixture'
require(not scratch.exists(), 'Refuse existing focused fixture')
scratch.mkdir()
try:
    work = scratch / 'work'
    outside = scratch / 'outside'
    work.mkdir()
    outside.mkdir()
    (work / 'fresh.json').write_text('{"status":"PASS"}\n')
    require(template['read_generated'](work / 'fresh.json', work) == b'{"status":"PASS"}\n',
            'Ordinary real generated file refused')
    passed('regular_fresh_generated_file_accepted')
    (outside / 'golden.json').write_text('{"status":"PASS"}\n')
    (work / 'stale-link.json').symlink_to(outside / 'golden.json')
    message = 'Fresh generated output must be a real regular file inside its workspace'
    rejects('output_symlink_rejected',
            lambda: template['read_generated'](work / 'stale-link.json', work), message)
    (work / 'linked-parent').symlink_to(outside, target_is_directory=True)
    rejects('output_ancestor_symlink_rejected',
            lambda: template['read_generated'](work / 'linked-parent' / 'golden.json', work), message)
    rejects('outside_regular_file_rejected',
            lambda: template['read_generated'](outside / 'golden.json', work), message)
    checker = {'name': 'focused', 'status_path': ['status'], 'status': 'PASS',
               'mode_fields': [], 'checks': None, 'embedded_checks': None,
               'check_names': [], 'embedded_control_names': [],
               'required_fields': [{'path': ['nested'], 'equals': {'answer': [1]}}]}
    require(template['positive_inventory'](checker, {'status': 'PASS', 'nested': {'answer': [1]}}, 'normal') == (0, 0),
            'Legitimate nested typed required field rejected')
    passed('nested_typed_required_field_accepted')
    rejects('nested_bool_integer_alias_rejected',
            lambda: template['positive_inventory'](checker, {'status': 'PASS', 'nested': {'answer': [True]}}, 'normal'),
            'Required positive output field differs')
    rejects('nested_integer_float_alias_rejected',
            lambda: template['positive_inventory'](checker, {'status': 'PASS', 'nested': {'answer': [1.0]}}, 'normal'),
            'Required positive output field differs')
    control = {'failure_format': 'traceback', 'exception': 'RuntimeError',
               'reason': 'named gate', 'diagnostic_prefix': 'REJECTED: named gate\n'}
    diagnostic = ('REJECTED: named gate\nTraceback (most recent call last):\n'
                  '  File "fixture.py", line 1, in <module>\n'
                  '    raise RuntimeError("named gate")\nRuntimeError: named gate\n')
    template['validate_failure_diagnostic'](control, diagnostic)
    passed('exact_configured_traceback_prefix_and_final_reason_accepted')
    bad_diagnostics = {
        'unexpected_traceback_leading_bytes_rejected': 'noise\n' + diagnostic,
        'missing_configured_traceback_prefix_rejected': diagnostic.split('\n', 1)[1],
        'wrong_configured_traceback_prefix_rejected': diagnostic.replace('REJECTED: named gate', 'REJECTED: other gate', 1),
        'unexpected_traceback_trailing_bytes_rejected': diagnostic + 'noise\n',
        'trailing_traceback_blank_line_rejected': diagnostic + '\n',
        'wrong_final_traceback_reason_rejected': diagnostic[:-len('RuntimeError: named gate\n')] + 'RuntimeError: other gate\n',
        'missing_final_traceback_newline_rejected': diagnostic[:-1],
    }
    for name, bad in bad_diagnostics.items():
        rejects(name, lambda value=bad: template['validate_failure_diagnostic'](control, value),
                'Traceback did not end with the exact named expected failure')
finally:
    shutil.rmtree(scratch)


class Template:
    pass


proxy = Template()
for name, value in template.items():
    setattr(proxy, name, value)
config_raw = (PUB / 'RUNNER_CONFIG.json').read_bytes()
config = json.loads(config_raw)
checkers, controls, roster, pins = builder['build'](config, ROOT, proxy)
roster_raw = (json.dumps(roster, indent=2) + '\n').encode()
require(roster_raw == (PUB / 'replay-tools/REPLAY_ROSTER.json').read_bytes(),
        'Actual roster differs from approved configuration and complete golden inventories')
passed('actual_roster_exactly_reconstructed_from_config_and_named_goldens')
constants = {'CAPSULE': builder['CAPSULE'], 'PYTHON_VERSION': config['python_version'],
             'DEPENDENCY_VERSIONS': config['dependency_versions'],
             'CONFIG_SHA256': hashlib.sha256(config_raw).hexdigest(),
             'CHECKERS': checkers, 'CONTROLS': controls}
source = (PUB / 'runner_template.py').read_text()
start = source.index('# BEGIN EMBEDDED APPROVED CONFIGURATION')
end = source.index('# END EMBEDDED APPROVED CONFIGURATION') + len('# END EMBEDDED APPROVED CONFIGURATION')
source = source[:start] + '# BEGIN EMBEDDED APPROVED CONFIGURATION\n' + '\n'.join(
    key + ' = ' + pprint.pformat(value, sort_dicts=False, width=100) for key, value in constants.items()
) + '\n# END EMBEDDED APPROVED CONFIGURATION' + source[end:]
require(source.encode() == (PUB / 'replay-tools/reproduce.py').read_bytes(),
        'Generated runner differs from reviewed template and approved embedded configuration')
passed('actual_runner_exactly_reconstructed_from_template_and_approved_config')
require(roster['counts'] == {'successful_checker_runs': 12,
                            'external_corruption_processes_rejected': 70,
                            'embedded_wrong_formula_controls_rejected': 16,
                            'finite_positive_check_groups_each_mode': 111}, 'Unexpected actual roster')
passed('counts_derived_from_actual_6_checkers_35_controls_111_names_and_8_embedded_names')
expected_modes = {'generator': [], 'root': [], 'astra': [(['python_optimized'], False, True)],
                  'archive': [(['python_optimized'], False, True)],
                  'sol': [(['runtime', 'optimization'], 0, 1)],
                  'referee': [(['python_optimization'], 0, 1)]}
for checker in checkers:
    fields = [(x['path'], x['normal'], x['optimized']) for x in checker['mode_fields']]
    require(fields == expected_modes[checker['name']], 'Unexpected scientific normalization field')
passed('all_normalized_fields_are_individually_typed_execution_metadata')
require(len(controls) == 35 and len({(c['checker'], c['name']) for c in controls}) == 35,
        'Duplicate or missing control identity')
for control in controls:
    if 'target' not in control:
        require(control['mechanism'] == 'formula comparison CLI injection',
                'CLI comparison controls described as stronger production mutation')
passed('all_35_controls_named_and_cli_comparison_mechanisms_described_precisely')
tracebacks = [control for control in controls if control['failure_format'] == 'traceback']
require(len(tracebacks) == 4 and all(control['checker'] == 'archive' and
        control['diagnostic_prefix'] == 'REJECTED: ' + control['reason'] + '\n'
        for control in tracebacks), 'Unexpected real traceback prefix contract')
passed('four_archive_traceback_controls_pin_exact_registered_diagnostic_prefix')

paths = ['publication/runner_template.py', 'publication/build_replay_runner.py',
         'publication/RUNNER_CONFIG.json', 'publication/replay-tools/reproduce.py',
         'publication/replay-tools/REPLAY_ROSTER.json', 'publication/RUNNER_PLUMBING_VALIDATION.json']
receipt = {'status': 'ACCEPT_REPLAY_HARDENING_AND_APPROVED_CONFIGURATION',
           'scope': 'Read-only source inspection, focused independent hardening regressions, complete approved-golden inventory inspection and deterministic build reconstruction. No actual scientific subprocess or frozen-archive fresh replay was run by this review.',
           'blockers': [], 'checks': checks, 'check_count': len(checks),
           'actual_roster_counts': roster['counts'],
           'source_pins': [{'path': name, 'bytes': len((ROOT / name).read_bytes()),
                            'sha256': hashlib.sha256((ROOT / name).read_bytes()).hexdigest()} for name in paths],
           'scientific_processes_executed': 0,
           'remaining_gate': 'Hash-specific frozen source inventory, external manifest authentication and actual fresh replay receipt still require review after source freeze.',
           'private_or_internal_research_added': False}
destination = HERE / 'RUNNER_PLUMBING_REVIEW.json'
destination.write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps({'status': receipt['status'], 'check_count': len(checks),
                  'review_sha256': hashlib.sha256(destination.read_bytes()).hexdigest(),
                  'counts': roster['counts']}, indent=2))
