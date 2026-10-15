"""Read-only independent audit of this checkpoint's frozen bytes and fresh evidence."""
import ast
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path, PurePosixPath
import stat
import zipfile

ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent
CAP = 'HDBLAST_QUADRATIC_SHELL_FORCING_20261015'
MANIFEST_PIN = '114d4698545e2f57e050b145d035af2a5b7fad73aa51791293e9666f2e02c1a6'
ARCHIVE_PIN = '71b550dcdd8649c8d5bea88ebaa5d00d7d5a7ffc40cbe7f4da6765694b75d2ed'
RUNNER_PIN = '0acb62ed9035b9ca90543811de009f4d075388a37778d331bd877b9cdfad991f'
ROSTER_PIN = '92b1c67f3622f0596e4e192c45751df602c69fbfd084a58f7929bb0ef4b0cb87'
REVIEW_PIN = '1c8c3d468e05fa5042db7af9785d41845974d068aad651ff05dbd0e45aec4f4a'
COUNTS = {'successful_checker_runs': 12, 'external_corruption_processes_rejected': 70,
          'embedded_wrong_formula_controls_rejected': 16,
          'finite_positive_check_groups_each_mode': 111}
checks = []


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def checked(name, condition=True):
    require(condition, name)
    checks.append(name)


def equal(left, right):
    if type(left) is not type(right):
        return False
    if isinstance(left, dict):
        return left.keys() == right.keys() and all(equal(left[k], right[k]) for k in left)
    if isinstance(left, list):
        return len(left) == len(right) and all(equal(a, b) for a, b in zip(left, right))
    return left == right


def at(data, path):
    for key in path:
        data = data[key]
    return data


def raw_file(path, base):
    require(path.is_file() and stat.S_ISREG(path.stat().st_mode) and
            not any(p.is_symlink() for p in [path, *path.parents]) and
            path.resolve().is_relative_to(base.resolve()), 'Real contained regular file required: ' + str(path))
    return path.read_bytes()


def pin(raw):
    return {'bytes': len(raw), 'sha256': sha256(raw).hexdigest()}


def match_record(record, raw):
    return type(record['bytes']) is int and record['bytes'] == len(raw) and record['sha256'] == sha256(raw).hexdigest()


frozen = ROOT / 'frozen-source'
source = frozen / CAP
manifest_raw = raw_file(source / 'MANIFEST.json', source)
archive_raw = raw_file(frozen / (CAP + '.zip'), frozen)
checked('external_manifest_pin_matches', sha256(manifest_raw).hexdigest() == MANIFEST_PIN)
checked('external_archive_pin_matches', sha256(archive_raw).hexdigest() == ARCHIVE_PIN)
manifest = json.loads(manifest_raw)
freeze = json.loads(raw_file(frozen / 'SOURCE_FREEZE_RECEIPT.json', frozen))
checked('freeze_identity_counts_and_external_records', manifest['schema'] == 'hdblast-frozen-source-v1' and
        manifest['capsule'] == CAP == freeze['capsule'] and freeze['status'] == 'PASS' and
        freeze['source_files'] == 51 and freeze['archive_members'] == 52 and
        match_record(freeze['manifest'], manifest_raw) and match_record(freeze['archive'], archive_raw))
inventory = {}
source_bytes = {}
for row in manifest['files']:
    name = row['path']
    relative = PurePosixPath(name)
    require(name and not relative.is_absolute() and '..' not in relative.parts and
            '\\' not in name and relative.as_posix() == name and
            name not in inventory and name != 'MANIFEST.json', 'Bad frozen source path')
    raw = raw_file(source / name, source)
    require(match_record(row, raw), 'Frozen source pin differs: ' + name)
    inventory[name] = row
    source_bytes[name] = raw
checked('all_51_unique_canonical_source_bytes_authenticated', len(inventory) == 51)
actual = set()
for path in source.rglob('*'):
    require(not path.is_symlink() and (path.is_dir() or stat.S_ISREG(path.stat().st_mode)),
            'Frozen symlink/special entry')
    if path.is_file():
        actual.add(path.relative_to(source).as_posix())
checked('full_tree_has_exact_52_manifest_and_source_files', actual == set(inventory) | {'MANIFEST.json'})
assembly = json.loads((ROOT / 'publication/SOURCE_ASSEMBLY_RECEIPT.json').read_bytes())
checked('source_matches_explicit_51_file_assembly', assembly['status'] == 'PASS' and
        assembly['source_files'] == 51 and {row['path']: row for row in assembly['files']} == inventory)
with zipfile.ZipFile(frozen / (CAP + '.zip')) as archive:
    names = archive.namelist()
    require(len(names) == len(set(names)) == 52 and
            set(names) == {CAP + '/' + name for name in actual}, 'ZIP exact inventory differs')
    require(archive.testzip() is None, 'ZIP CRC failure')
    for entry in archive.infolist():
        require(stat.S_IFMT(entry.external_attr >> 16) == stat.S_IFREG, 'Nonregular ZIP entry')
        name = entry.filename[len(CAP) + 1:]
        require(archive.read(entry) == (manifest_raw if name == 'MANIFEST.json' else source_bytes[name]),
                'ZIP and source bytes differ')
checked('all_52_regular_zip_members_crc_and_complete_bytes_match_source')
for name in actual:
    parts = set(PurePosixPath(name).parts)
    require(not parts & {'.git', '__pycache__', 'future_work', 'publication', 'private', 'internal',
                         'raw', 'raw-responses', 'scratch'} and not name.endswith(('.pyc', '.pyo')),
            'Out-of-scope source filename: ' + name)
checked('no_scratch_private_future_work_raw_response_or_repository_metadata_files')
checked('final_scientific_review_receipt_frozen_exactly', sha256(source_bytes['referee/FINAL_QUADRATIC_SOURCE_REVIEW.json']).hexdigest() == REVIEW_PIN)
scientific_review = json.loads(source_bytes['referee/FINAL_QUADRATIC_SOURCE_REVIEW.json'])
for name, row in scientific_review['pins'].items():
    require(name in source_bytes and match_record(row, source_bytes[name]), 'Prior accepted scientific review pin absent or changed: ' + name)
checked('all_13_final_scientific_review_pins_present_and_unchanged', len(scientific_review['pins']) == 13)
inputs = json.loads(source_bytes['inputs/INPUT_PINS.json'])
require(inputs['status'] == 'PASS' and len(inputs['files']) == 9, 'Input pin inventory differs')
prior_roots = {
    'HDBLAST_INERTIAL_SHELL_WAVEFORM_20261015': ROOT.parent / 'inertial-waveform-20261015/frozen-source/HDBLAST_INERTIAL_SHELL_WAVEFORM_20261015',
    'HDBLAST_FULL_SECOND_ORDER_CONTACT_20261014': ROOT / 'publication/hdblast-worktree/research/HDBLAST_FULL_SECOND_ORDER_CONTACT_20261014',
}
prior_manifests = {}
for name, prior in prior_roots.items():
    prior_raw = raw_file(prior / 'MANIFEST.json', prior)
    expected = inputs['external_inertial_manifest_sha256'] if 'INERTIAL' in name else inputs['external_full_contact_manifest_sha256']
    require(sha256(prior_raw).hexdigest() == expected, 'Previous external manifest differs')
    prior_manifests[name] = {row['path']: row for row in json.loads(prior_raw)['files']}
for row in inputs['files']:
    name = row['path']
    provenance = row['original_public_source']
    prior_cap = provenance['capsule']
    prior_name = provenance['path'].split('/' + prior_cap + '/', 1)[1]
    prior_row = prior_manifests[prior_cap][prior_name]
    require(match_record(row, source_bytes[name]) and match_record(prior_row, source_bytes[name]) and
            raw_file(prior_roots[prior_cap] / prior_name, prior_roots[prior_cap]) == source_bytes[name],
            'Pinned input provenance differs: ' + name)
checked('all_9_minimal_public_inputs_match_pinned_previous_manifests_and_source_bytes')
checked('runner_and_roster_match_reviewed_final_hashes',
        sha256(source_bytes['reproduce.py']).hexdigest() == RUNNER_PIN and
        sha256(source_bytes['REPLAY_ROSTER.json']).hexdigest() == ROSTER_PIN)
constants = {}
for node in ast.parse(source_bytes['reproduce.py']).body:
    if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
        name = node.targets[0].id
        if name in {'CAPSULE', 'PYTHON_VERSION', 'DEPENDENCY_VERSIONS', 'CONFIG_SHA256', 'CHECKERS', 'CONTROLS'}:
            constants[name] = ast.literal_eval(node.value)
checkers = constants['CHECKERS']
controls = constants['CONTROLS']
roster = json.loads(source_bytes['REPLAY_ROSTER.json'])
require(constants['CAPSULE'] == CAP and len(checkers) == 6 and len(controls) == 35, 'Embedded runner identity differs')
checked('frozen_roster_counts_equal_actual_named_embedded_runner_inventory', equal(roster['counts'], COUNTS) and
        sum(c['groups'] for c in checkers) == 111 and sum(c['embedded'] for c in checkers) == 8)
expected_positive = []
for checker in checkers:
    for name in [checker['script'], *checker['dependencies'], checker['golden'],
                 checker.get('optimized_golden', checker['golden'])]:
        require(name in inventory, 'Missing minimal checker dependency or golden: ' + name)
    for mode in ('normal', 'optimized'):
        name = checker.get('optimized_golden', checker['golden']) if mode == 'optimized' else checker['golden']
        golden = json.loads(source_bytes[name])
        require(at(golden, checker['status_path']) == checker['status'], 'Wrong positive golden status')
        for field in checker['mode_fields']:
            require(equal(at(golden, field['path']), field[mode]), 'Wrong typed mode metadata')
        for specification, expected in ((checker['checks'], checker['check_names']),
                                        (checker['embedded_checks'], checker['embedded_control_names'])):
            value = [] if specification is None else at(golden, specification['path'])
            if isinstance(value, dict):
                require(all(v is True for v in value.values()), 'Failed named map')
                value = list(value)
            require(equal(value, expected) and len(value) == len(set(value)), 'Named positive/embedded golden inventory differs')
    expected_positive.append({'name': checker['name'], 'script': checker['script'],
        'finite_positive_groups_each_mode': checker['groups'], 'embedded_controls_each_mode': checker['embedded'],
        'finite_positive_check_names': checker['check_names'], 'embedded_control_names': checker['embedded_control_names']})
expected_controls = [{'checker': c['checker'], 'name': c['name'], 'expected_reason': c['reason'], 'mechanism': c['mechanism']} for c in controls]
checked('all_checker_and_control_names_scripts_reasons_and_mechanisms_match_frozen_roster',
        equal(roster['positive_processes'], expected_positive) and equal(roster['external_controls'], expected_controls))
checked('all_checker_minimal_dependencies_and_complete_goldens_present')

receipt_path = ROOT / 'fresh-zip-replay/REPRODUCTION_RECEIPT.json'
if receipt_path.exists():
    replay_root = receipt_path.parent
    fresh_raw = raw_file(receipt_path, replay_root)
    fresh = json.loads(fresh_raw)
    extracted = ROOT / 'fresh-zip-source' / CAP
    require({path.relative_to(extracted).as_posix() for path in extracted.rglob('*') if path.is_file()} == actual,
            'Fresh extracted source inventory differs')
    for name in actual:
        require(raw_file(extracted / name, extracted) == (manifest_raw if name == 'MANIFEST.json' else source_bytes[name]),
                'Fresh extracted source bytes differ')
    authentication = json.loads(raw_file(ROOT / 'FRESH_ZIP_AUTHENTICATION.json', ROOT))
    require(authentication['status'] == 'PASS' and authentication['extracted_members'] == 52 and
            match_record(authentication['archive'], archive_raw) and match_record(authentication['manifest'], manifest_raw) and
            authentication['bootstrap_pins_before_execution'] == {
                name: sha256(manifest_raw if name == 'MANIFEST.json' else source_bytes[name]).hexdigest()
                for name in ('MANIFEST.json', 'reproduce.py', 'verify_manifest.py')},
            'External pre-execution ZIP extraction authentication differs')
    checked('all_52_fresh_extracted_source_bytes_and_external_bootstrap_authentication_match')
    require(fresh['status'] == 'PASS' and fresh['capsule'] == CAP and
            fresh['authentication']['manifest_sha256'] == MANIFEST_PIN and
            fresh['authentication']['files_verified'] == 51, 'Fresh receipt authentication differs')
    require(all(equal(fresh[key], value) for key, value in COUNTS.items()), 'Fresh receipt declared counts differ')
    for key, name in [('runner', 'reproduce.py'), ('integrity_verifier', 'verify_manifest.py'),
                      ('requirements', 'requirements.txt'), ('roster', 'REPLAY_ROSTER.json')]:
        require(equal(fresh[key], inventory[name]), 'Fresh receipt source pins differ')
    require(fresh['embedded_configuration_sha256'] == constants['CONFIG_SHA256'] and
            equal(fresh['dependency_versions'], constants['DEPENDENCY_VERSIONS']) and
            fresh['python']['version'].split()[0] == constants['PYTHON_VERSION'], 'Fresh runtime/configuration differs')
    checked('fresh_replay_exact_manifest_runtime_runner_verifier_requirements_roster_pins')
    expected_ids = [(mode, 'positive', c['name'], None) for mode in ('normal', 'optimized') for c in checkers]
    expected_ids += [(mode, 'external_corruption', c['checker'], c['name']) for mode in ('normal', 'optimized') for c in controls]
    steps = fresh['steps']
    actual_ids = [(row['mode'], row['kind'], row['checker'], row.get('control')) for row in steps]
    require(len(actual_ids) == len(set(actual_ids)) == 82 and set(actual_ids) == set(expected_ids), 'Fresh process identity inventory incomplete')
    checked('exact_82_process_identities_include_12_positives_and_70_external_controls')
    by_name = {c['name']: c for c in checkers}
    by_control = {(c['checker'], c['name']): c for c in controls}
    positive_data = {}
    preserved_files = set()
    for row in steps:
        checker = by_name[row['checker']]
        control = by_control.get((row['checker'], row.get('control')))
        work_prefix = Path(row['mode']) / row['checker'] / (row['control'] if control else 'positive')
        work = replay_root / work_prefix
        require(row['script'] == checker['script'] and
                row['flags'] == ['-I', '-B'] + (['-O'] if row['mode'] == 'optimized' else []) and
                row['arguments'] == (control.get('arguments', checker['arguments']) if control else checker['arguments']),
                'Actual fresh command identity differs')
        def read_record(record):
            raw = raw_file(replay_root / record['path'], replay_root)
            require(match_record(record, raw), 'Fresh record bytes/hash differ')
            require(Path(record['path']).is_relative_to(work_prefix), 'Process record crosses workspaces')
            preserved_files.add(record['path'])
            return raw
        diagnostic = read_record(row['log']).decode()
        output_raw = read_record(row['output'])
        data = json.loads(output_raw)
        expected_inputs = [checker['script'], *checker['dependencies']]
        require([copy['source']['path'] for copy in row['copied_inputs']] == expected_inputs, 'Fresh copied input inventory differs')
        for copy in row['copied_inputs']:
            name = copy['source']['path']
            require(equal(copy['source'], inventory[name]) and
                    copy['copy']['path'] == (work_prefix / name).as_posix() and
                    match_record(copy['copy'], source_bytes[name]), 'Authenticated initial copy snapshot differs')
            current = raw_file(work / name, work)
            if not control or name != control.get('target'):
                require(current == source_bytes[name], 'Untargeted copied source changed')
        if control:
            require(row['exit_code'] == row['expected_failure_exit_code'] == 1 and
                    row['expected_reason_detected'] is True and row['reported_failure'] == control['reason'] and
                    row['expected_reason'] == control['reason'] and row['control_mechanism'] == control['mechanism'],
                    'Fresh external failure metadata differs')
            require(checker['output'] is None or not ((work / checker['output']).exists() or
                    (work / checker['output']).is_symlink()), 'Rejected process left result output')
            require(data['status'] == 'FAIL' and data['failure'] == control['reason'] and
                    data['exit_code'] == 1 and data['passing_output_absent'] is True and
                    data['complete_diagnostic'] == diagnostic, 'Synthetic failure record differs from full log')
            if control['failure_format'] == 'fail_line':
                require(diagnostic == 'FAIL: ' + control['reason'] + '\n', 'Wrong named fail line')
            elif control['failure_format'] == 'json':
                failed = json.loads(diagnostic)
                require(at(failed, control.get('status_path', ['status'])) == control.get('status', 'FAIL') and
                        at(failed, control['failure_path']) == control['reason'], 'Wrong named JSON failure')
            else:
                require(diagnostic.startswith(control.get('diagnostic_prefix', '') + 'Traceback (most recent call last):\n') and
                        diagnostic.endswith(control['exception'] + ': ' + control['reason'] + '\n') and
                        diagnostic.splitlines()[-1] == control['exception'] + ': ' + control['reason'], 'Wrong exact traceback failure')
            if 'target' in control:
                name = control['target']
                original = source_bytes[name]
                if 'append' in control:
                    mutated = original + control['append'].encode()
                elif 'json_path' in control:
                    document = json.loads(original)
                    parent = document
                    for key in control['json_path'][:-1]:
                        parent = parent[key]
                    key = control['json_path'][-1]
                    if 'factor' in control:
                        value = parent[key]
                        parent[key] = [str(Fraction(v) * Fraction(control['factor'])) for v in value] if isinstance(value, list) else str(Fraction(value) * Fraction(control['factor']))
                    else:
                        parent[key] = control['value']
                    mutated = (json.dumps(document, indent=2) + '\n').encode()
                else:
                    require(original.decode().count(control['before']) == 1, 'Nonunique mutation')
                    mutated = original.decode().replace(control['before'], control['after'], 1).encode()
                require(mutated != original and raw_file(work / name, work) == mutated and
                        equal(row['mutation']['source'], inventory[name]) and
                        equal(row['mutation']['recipe'], control) and
                        row['mutation']['independent_checker_unchanged'] is True and
                        read_record(row['mutation']['mutated']) == mutated, 'Copied-input corruption provenance differs')
            else:
                require(row['mutation'] is None and control['mechanism'] == 'formula comparison CLI injection',
                        'CLI comparison control overstated as mutated production formula')
        else:
            require(row['exit_code'] == 0 and row['scientific_output_matches'] is True and
                    row['ignored_golden_fields'] == [], 'Fresh positive process failed/ignored golden fields')
            golden_name = checker.get('optimized_golden', checker['golden']) if row['mode'] == 'optimized' else checker['golden']
            require(equal(row['golden'], inventory[golden_name]) and equal(data, json.loads(source_bytes[golden_name])), 'Complete positive scientific golden differs')
            if checker['byte_equal']:
                require(output_raw == source_bytes[golden_name] and row['complete_generated_json_bytes_equal'] is True, 'Complete producer JSON byte equality failed')
            require(output_raw == raw_file(work / (checker['output'] or 'execution.log'), work), 'Fresh result snapshot differs from actual producer output')
            require(row['finite_positive_check_groups'] == checker['groups'] and
                    row['embedded_wrong_formula_controls'] == checker['embedded'] and
                    equal(row['finite_positive_check_names'], checker['check_names']) and
                    equal(row['embedded_control_names'], checker['embedded_control_names']), 'Fresh positive named inventory differs')
            require(len(row['additional_outputs']) == len(checker['additional_outputs']), 'Missing complete generated additional output')
            for additional, (name, golden) in zip(row['additional_outputs'], checker['additional_outputs']):
                require(read_record(additional['output']) == source_bytes[golden] and
                        equal(additional['golden'], inventory[golden]) and additional['complete_bytes_equal'] is True and
                        additional['output']['path'] == (work_prefix / name).as_posix(), 'Complete generated CSV bytes differ')
            positive_data[(row['mode'], checker['name'])] = data
    checked('all_82_commands_logs_results_initial_copy_snapshots_and_16_actual_mutated_inputs_authenticated')
    generated_artifacts = {item['path'] for row in steps for item in
                           [row['log'], row['output'], *[part['output'] for part in row['additional_outputs']]]}
    checked('all_166_distinct_generated_logs_results_and_additional_artifacts_authenticated',
            len(generated_artifacts) == 166 and generated_artifacts <= preserved_files)
    checked('all_70_exact_named_failures_preserve_independent_checker_and_have_no_passing_output')
    checked('all_12_complete_positive_json_goldens_match_with_2_full_generated_json_and_csv_pairs')
    expected_agreements = []
    for checker in checkers:
        pair = [json.loads(json.dumps(positive_data[(mode, checker['name'])])) for mode in ('normal', 'optimized')]
        for item, mode in zip(pair, ('normal', 'optimized')):
            for field in checker['mode_fields']:
                require(equal(at(item, field['path']), field[mode]), 'Fresh exact execution metadata differs')
                parent = item
                for key in field['path'][:-1]:
                    parent = parent[key]
                del parent[field['path'][-1]]
        require(equal(*pair), 'Scientific fields differ across Python modes')
        expected_agreements.append({'checker': checker['name'], 'scientific_json_equal': True,
                                    'only_normalized_fields': [field['path'] for field in checker['mode_fields']]})
    require(equal(fresh['mode_agreements'], expected_agreements) and
            equal(fresh['finite_positive_check_groups_by_mode'], {'normal': 111, 'optimized': 111}), 'Mode agreement inventories differ')
    checked('all_6_mode_agreements_normalize_only_pinned_typed_execution_metadata')
    checked('actual_positive_external_embedded_and_named_group_counts_reconstructed',
            sum(row['kind'] == 'positive' for row in steps) == 12 and
            sum(row['kind'] == 'external_corruption' for row in steps) == 70 and
            sum(row.get('embedded_wrong_formula_controls', 0) for row in steps) == 16 and
            all(sum(row.get('finite_positive_check_groups', 0) for row in steps if row['mode'] == mode) == 111
                for mode in ('normal', 'optimized')))
    require(raw_file(source / 'MANIFEST.json', source) == manifest_raw and
            raw_file(frozen / (CAP + '.zip'), frozen) == archive_raw and
            raw_file(receipt_path, replay_root) == fresh_raw, 'Pinned delivery changed during review')
    for name, raw in source_bytes.items():
        require(raw_file(source / name, source) == raw, 'Frozen science changed during full artifact review')
    checked('all_51_frozen_source_files_manifest_zip_and_fresh_receipt_immutable_across_review')
    fresh_pin = pin(fresh_raw)
    status = 'ACCEPT_FROZEN_SOURCE_AND_AUTHORITATIVE_FRESH_REPLAY'
else:
    fresh_pin = None
    status = 'ACCEPT_FROZEN_SOURCE_FRESH_REPLAY_PENDING'

review = {'schema': 'hdblast-quadratic-frozen-delivery-independent-review-v1', 'status': status,
          'blockers': [], 'manifest_sha256': MANIFEST_PIN, 'archive_sha256': ARCHIVE_PIN,
          'source_files': 51, 'archive_members': 52, 'approved_roster_counts': COUNTS,
          'scientific_review_sha256': REVIEW_PIN, 'runner_sha256': RUNNER_PIN, 'roster_sha256': ROSTER_PIN,
          'fresh_replay_receipt': fresh_pin, 'checks': checks, 'check_count': len(checks),
          'fresh_generated_artifacts_authenticated': 166 if fresh_pin else None,
          'scope': 'Independent read-only byte authentication and explicit current-checkpoint process-identity audit; no archived scientific program imported or executed, no Git/network or frozen-file mutation. This audit authenticates calculations and accepted source scope; it does not prove analytic PDE hypotheses, a quantitative full U2 response, an effective nonlinear remainder, a finite-amplitude cosmology or novelty.',
          'scope_exclusions_verified': ['future_work files', 'private/internal project files', 'raw web responses', 'scratch outputs', 'repository metadata', 'bytecode caches']}
destination = HERE / ('FROZEN_SOURCE_AND_FRESH_REPLAY_REVIEW.json' if fresh_pin else 'FROZEN_SOURCE_REVIEW.json')
destination.write_text(json.dumps(review, indent=2) + '\n')
print(json.dumps({'status': status, 'checks': len(checks), 'review_file': destination.name,
                  'review_sha256': sha256(destination.read_bytes()).hexdigest(),
                  'fresh_replay_receipt': fresh_pin}, indent=2))
