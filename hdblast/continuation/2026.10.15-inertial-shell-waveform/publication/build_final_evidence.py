"""Stage final additive publication evidence after all real merge/CI/readback receipts pass.

No Git mutation or network request occurs. Initial delivery files are copied
from the reviewed immutable Git commit and preserved exactly. Final evidence
PR readback stays outside this stage to avoid recursive publication receipts.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path, PurePosixPath
import re
import subprocess

PUBLIC = 'maldonado-research/HDblast'
SHARED = 'maldonado-research/maldonado-research.github.io'
CAPSULE = 'research/HDBLAST_INERTIAL_SHELL_WAVEFORM_20261015'
DELIVERY = 'hdblast/continuation/2026.10.15-inertial-shell-waveform'
ROSTER_PATH = 'REPLAY_ROSTER.json'
ROSTER_KEYS = ('successful_checker_runs', 'external_corruption_processes_rejected',
               'embedded_wrong_formula_controls_rejected', 'finite_positive_check_groups_each_mode')


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def canonical(name):
    require(isinstance(name, str), 'String relative path required')
    path = PurePosixPath(name)
    require(name and not path.is_absolute() and '..' not in path.parts
            and '\\' not in name and path.as_posix() == name, 'Canonical relative path required')
    return name


def read_file(path):
    require(path.is_file() and not any(p.is_symlink() for p in [path, *path.parents]),
            'Real regular evidence file required: '+str(path))
    return path.read_bytes()


def validate_roster(roster, replay):
    """Counts come from a separately approved frozen source file, not old checkpoints."""
    require(roster['schema'] == 'hdblast-replay-roster-v1'
            and roster['status'] == 'APPROVED_FOR_FROZEN_REPLAY'
            and roster['capsule'] == PurePosixPath(CAPSULE).name
            and roster['modes'] == ['normal', 'optimized'], 'Approved current frozen replay roster required')
    counts = roster['counts']
    require(set(counts) == set(ROSTER_KEYS)
            and all(type(v) is int and v >= 0 for v in counts.values())
            and counts['successful_checker_runs'] > 0
            and counts['finite_positive_check_groups_each_mode'] > 0, 'Valid explicit replay counts required')
    require(replay['status'] == 'PASS' and replay['capsule'] == roster['capsule']
            and all(replay[k] == v for k, v in counts.items()), 'Fresh replay differs from its approved frozen roster')
    positive = {r['name']: r for r in roster['positive_processes']}
    controls = {(r['checker'], r['name']): r for r in roster['external_controls']}
    require(positive and len(positive) == len(roster['positive_processes'])
            and len(controls) == len(roster['external_controls'])
            and all(r['checker'] in positive for r in controls.values()),
            'Unique approved checker/control identities required')
    expected_identities = {(name, 'positive', 'positive') for name in positive}
    expected_identities |= {(checker, name, 'external_corruption') for checker, name in controls}
    steps = replay['steps']
    require(all(r['mode'] in roster['modes'] and r['kind'] in ('positive', 'external_corruption')
                for r in steps), 'Unexpected replay mode or process kind')
    require(len({(r['mode'], r['checker'], r.get('control', 'positive')) for r in steps}) == len(steps),
            'Repeated replay process identity')
    actual = {'successful_checker_runs': sum(r['kind'] == 'positive' for r in steps),
              'external_corruption_processes_rejected': sum(r['kind'] == 'external_corruption' for r in steps),
              'embedded_wrong_formula_controls_rejected': sum(r.get('embedded_wrong_formula_controls', 0) for r in steps)}
    for key, value in actual.items():
        require(value == counts[key], 'Replay step inventory differs: '+key)
    for mode in roster['modes']:
        mode_steps = [r for r in steps if r['mode'] == mode]
        require(sum(r.get('finite_positive_check_groups', 0) for r in mode_steps)
                == counts['finite_positive_check_groups_each_mode'], 'Per-mode finite group inventory differs')
        require({(r['checker'], r.get('control', 'positive'), r['kind']) for r in mode_steps}
                == expected_identities, 'Each mode must cover the explicit approved process roster')
    for row in steps:
        require(row['checker'] in positive and row['script'] == positive[row['checker']]['script'],
                'Executed checker identity/script differs from approved roster')
        if row['kind'] == 'positive':
            require(row['exit_code'] == 0 and row['scientific_output_matches'] is True,
                    'Positive replay process did not reproduce its scientific output')
            approved = positive[row['checker']]
            require(row['finite_positive_check_groups'] == approved['finite_positive_groups_each_mode']
                    and row['embedded_wrong_formula_controls'] == approved['embedded_controls_each_mode'],
                    'Checker groups or embedded controls differ from approved roster')
        else:
            approved = controls.get((row['checker'], row['control']))
            require(approved is not None and row['expected_reason'] == row['reported_failure']
                    == approved['expected_reason'], 'Named rejection differs from its approved control')
            require(row['exit_code'] == row['expected_failure_exit_code'] == 1
                    and row['expected_reason_detected'] is True and bool(row['expected_reason']),
                    'Named external wrong-formula rejection required')
    return counts


def main():
    task = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', type=Path, default=task/'evidence-staging')
    parser.add_argument('--source-graph-audit', type=Path, default=task/'publication/COMMITTED_SOURCE_GRAPH_AUDIT.json')
    parser.add_argument('--readback-preflight', type=Path, default=task/'publication/READBACK_HELPER_PREFLIGHT.json')
    parser.add_argument('--include', action='append', default=[], help='Additional task-relative support/planning file, preserving that path')
    args = parser.parse_args()
    stage = args.stage.absolute()
    require(not stage.exists() and not any(p.is_symlink() for p in [stage, *stage.parents]), 'Fresh real evidence stage required')
    raw_inputs = {}
    def load(rel):
        raw = read_file(task/rel)
        raw_inputs[rel] = raw
        return json.loads(raw)
    prepared = load('PREPARED_COMMITS.json')
    frozen = load('frozen-source/SOURCE_FREEZE_RECEIPT.json')
    figures = load('publication/FIGURE_INTEGRATION_RECEIPT.json')
    pre = load('PRE_MERGE_CHECKS.json')
    merges = load('MERGE_RESULTS.json')
    post = load('POST_MERGE_WORKFLOW_RUNS.json')
    live = load('FINAL_LIVE_READBACK.json')
    replay = load('fresh-zip-replay/REPRODUCTION_RECEIPT.json')
    require(live['status'] == 'PASS', 'Completed real source/site/figure readback required')
    helper_path = task/'publication/readback_published.py'
    helper_raw = read_file(helper_path)
    require(live['helper_sha256'] == sha(helper_raw), 'Readback helper changed after execution')
    expected_inputs = ['PREPARED_COMMITS.json', 'frozen-source/SOURCE_FREEZE_RECEIPT.json',
                       'publication/FIGURE_INTEGRATION_RECEIPT.json', 'PRE_MERGE_CHECKS.json',
                       'MERGE_RESULTS.json', 'POST_MERGE_WORKFLOW_RUNS.json']
    require(set(live['input_receipts']) == {Path(p).name for p in expected_inputs}, 'Live readback input inventory differs')
    for rel in expected_inputs:
        raw, pin = raw_inputs[rel], live['input_receipts'][Path(rel).name]
        require(pin['bytes'] == len(raw) and pin['sha256'] == sha(raw), 'A readback input changed: '+rel)
    # Import only the currently hash-authenticated newly authored helper.
    # Its import has no HTTP or Git side effects.
    spec = importlib.util.spec_from_file_location('authenticated_waveform_readback', helper_path)
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    pins = helper.prepared_pins(prepared, require_pr=True)
    roots = {repo: Path(pin['worktree']).resolve() for repo, pin in pins.items()}
    require(not any(stage.resolve().is_relative_to(root) for root in roots.values())
            and not stage.resolve().is_relative_to((task/'frozen-source').resolve()),
            'Evidence stage must stay outside reviewed repositories and frozen source')
    targets, inventory, scientific_pins = helper.build_targets(roots, prepared, frozen, figures)
    snapshots = helper.validate_observations(pre, post, merges, prepared, roots, live['input_receipts'])
    require(live['inventory'] == inventory and live['pins'] == scientific_pins
            and live['target_count'] == len(targets) == len(live['targets']), 'Live inventory or source pins differ')
    by_url = {r['url']: r for r in live['targets']}
    require(len(by_url) == len(targets) and set(by_url) == {r['url'] for r in targets}, 'Complete live URL set required')
    for expected in targets:
        row = by_url[expected['url']]
        require(all(row[k] == expected[k] for k in ('repository', 'path', 'kind', 'expected_bytes', 'expected_sha256'))
                and row['http_status'] == 200 and row['byte_equal'] is True
                and row['bytes'] == expected['expected_bytes'] and row['sha256'] == expected['expected_sha256'],
                'Successful byte-equal live target required: '+expected['url'])
    source_root, source_head = roots[PUBLIC], pins[PUBLIC]['head_sha']
    roster_raw = subprocess.check_output(['git', '-C', str(source_root), 'show',
                                          source_head+':'+CAPSULE+'/'+ROSTER_PATH])
    source_manifest = json.loads(subprocess.check_output(['git', '-C', str(source_root), 'show',
                                                         source_head+':'+CAPSULE+'/MANIFEST.json']))
    roster_records = [r for r in source_manifest['files'] if r['path'] == ROSTER_PATH]
    require(len(roster_records) == 1 and roster_records[0]['bytes'] == len(roster_raw)
            and roster_records[0]['sha256'] == sha(roster_raw), 'Replay roster must be authenticated by frozen source')
    replay_counts = validate_roster(json.loads(roster_raw), replay)
    require(replay['authentication']['manifest_sha256'] == scientific_pins['manifest']
            and replay['runner']['sha256'] == scientific_pins['runner']
            and replay['integrity_verifier']['sha256'] == scientific_pins['verifier'],
            'Fresh replay belongs to a different frozen source, runner or verifier')
    graph_raw = read_file(args.source_graph_audit)
    graph = json.loads(graph_raw)
    require(graph['status'] == 'PASS' and graph['basis'] == 'immutable_commit'
            and graph['commit'] == pins[PUBLIC]['head_sha'] and graph['tree'] == pins[PUBLIC]['reviewed_tree']
            and graph['capsule_path'] == CAPSULE and graph['all_git_bytes_equal_frozen_manifest_and_zip'] is True
            and graph['source_files'] == frozen['source_files']
            and graph['tracked_source_and_manifest'] == graph['archive_members'] == frozen['archive_members']
            and graph['manifest_sha256'] == scientific_pins['manifest']
            and graph['archive_sha256'] == scientific_pins['archive'], 'Passing reviewed-commit complete source-graph audit required')
    preflight_raw = read_file(args.readback_preflight)
    preflight = json.loads(preflight_raw)
    require(preflight['status'] == 'PREPARED_NOT_EXECUTED' and preflight['inventory'] == inventory
            and preflight['pins'] == scientific_pins and preflight['helper_sha256'] == sha(helper_raw),
            'Matching local readback preparation receipt required')
    public_root, head = roots[PUBLIC], pins[PUBLIC]['head_sha']
    def committed(name):
        return subprocess.check_output(['git', '-C', str(public_root), 'show', head+':'+name])
    initial_raw = committed(DELIVERY+'/INITIAL_REPLAY_MANIFEST.json')
    require(sha(initial_raw) == scientific_pins['initial_delivery_manifest'], 'Initial-delivery manifest pin differs')
    initial = json.loads(initial_raw)
    initial_files = {'INITIAL_REPLAY_MANIFEST.json': initial_raw}
    for row in initial['files']:
        name = canonical(row['path'])
        require(name not in initial_files, 'Duplicate/reserved initial evidence path')
        raw = committed(DELIVERY+'/'+name)
        require(type(row['bytes']) is int and row['bytes'] == len(raw) and row['sha256'] == sha(raw),
                'Initial committed evidence differs')
        initial_files[name] = raw
    require(len(initial_files) == inventory['initial_delivery_files'], 'Initial delivery count differs')
    require(initial_files.get('REPRODUCTION_RECEIPT.json') == raw_inputs['fresh-zip-replay/REPRODUCTION_RECEIPT.json']
            and initial_files.get('SOURCE_FREEZE_RECEIPT.json') == raw_inputs['frozen-source/SOURCE_FREEZE_RECEIPT.json']
            and initial_files.get('FIGURE_INTEGRATION_RECEIPT.json') == raw_inputs['publication/FIGURE_INTEGRATION_RECEIPT.json'],
            'Published initial replay/freeze/figure receipts must match the authenticated current inputs')
    staged = dict(initial_files)
    def add(name, raw):
        canonical(name)
        if name in staged:
            require(staged[name] == raw, 'Existing initial evidence must not be altered: '+name)
        else:
            staged[name] = raw
    for rel in ('PREPARED_COMMITS.json', 'PRE_MERGE_CHECKS.json', 'MERGE_RESULTS.json',
                'POST_MERGE_WORKFLOW_RUNS.json', 'FINAL_LIVE_READBACK.json'):
        add(rel, raw_inputs[rel])
    add('publication/COMMITTED_SOURCE_GRAPH_AUDIT.json', graph_raw)
    add('publication/READBACK_HELPER_PREFLIGHT.json', preflight_raw)
    helper_names = ['monitor_premerge.py', 'monitor_postmerge.py', 'readback_published.py',
                    'build_final_evidence.py', 'readback_evidence.py', 'verify_git_inventory.py',
                    'validate_monitor_readback_tools.py', 'MONITOR_READBACK_TOOL_VALIDATION.json',
                    'MONITOR_READBACK_USAGE.md', 'validate_evidence_roster.py',
                    'EVIDENCE_ROSTER_TOOL_VALIDATION.json', 'REPLAY_ROSTER_CONTRACT.md',
                    'verify_publication_scope.py']
    for name in helper_names:
        add('publication/'+name, read_file(task/'publication'/name))
    helper_manifest = {'schema': 'hdblast-publication-support-source-v1',
                       'scope': 'Support tool/document hashes only; no remote observation is inferred.',
                       'files': [{'path': 'publication/'+name, 'bytes': len(staged['publication/'+name]),
                                  'sha256': sha(staged['publication/'+name])} for name in helper_names]}
    add('publication/SUPPORT_HELPER_MANIFEST.json', (json.dumps(helper_manifest, indent=2)+'\n').encode())
    for rel in args.include:
        canonical(rel)
        require('internal' not in PurePosixPath(rel).parts and 'private' not in PurePosixPath(rel).parts
                and (not rel.startswith('millennium/')
                     or rel == 'millennium/HDBLAST_METHODS_TRANSFER_PUBLIC_SAFE.md'),
                'Private audit or paper-body files cannot be added as publication evidence')
        source = task/rel
        require(source.resolve().is_relative_to(task.resolve()), 'Additional evidence outside task')
        add(rel, read_file(source))
    add('PUBLICATION_COMPLETION.md', (f'''# Completed research and website publication

This additive record supplements the preserved initial replay evidence.
The source and website PRs are merged at the reviewed trees; three premerge
scientific/regression workflows and five exact-merge workflows/deployments
passed. The complete source, initial evidence, nine summaries, four live HTML
pages and any declared figure assets were read back byte for byte. Actual inventories,
URLs, commit pins and counts are in `PUBLICATION_RECEIPT.json` and its inputs.

The computational scope is {replay_counts['successful_checker_runs']} positive processes,
{replay_counts['external_corruption_processes_rejected']} external corrupted-formula process
rejections, {replay_counts['embedded_wrong_formula_controls_rejected']} embedded wrong-formula
controls and {replay_counts['finite_positive_check_groups_each_mode']} finite positive groups
per mode, as specified in the authenticated source `REPLAY_ROSTER.json`.
The source describes a conditional inertial first-variation waveform enclosure
for its specified model state. These checks do not formally certify
the analytic PDE hypotheses or provide external peer review.

No effective nonlinear epsilon/remainder bound, finite-amplitude nonlinear
waveform, observational scale or higher-dimensional Big Bang origin is
established. No new Zenodo version or APS submission is part of this record.
The later evidence-only PR and its final readback remain separate local
records, avoiding a publication receipt that tries to certify its own upload.
''').encode())
    receipt = {'schema': 'hdblast-inertial-waveform-publication-v1',
               'status': 'RESEARCH_SITES_AND_INITIAL_REPLAY_PUBLISHED_VERIFIED',
               'observed_utc': datetime.now(timezone.utc).isoformat(),
               'research_and_site_pull_requests': merges, 'reviewed_merge_snapshots': snapshots,
               'source_freeze': frozen, 'source_pins': scientific_pins,
               'initial_delivery_files_preserved': len(initial_files),
               'initial_delivery_manifest_sha256': sha(initial_raw),
               'premerge_checks': {'status': 'PASS', 'reviewed_workflow_count': 3},
               'postmerge_workflows': post,
               'live_readback': {'status': 'PASS', 'target_count': len(targets), 'inventory': inventory,
                                 'receipt_sha256': sha(raw_inputs['FINAL_LIVE_READBACK.json'])},
               'fresh_replay': {'status': 'PASS', **replay_counts, 'frozen_roster_sha256': sha(roster_raw)},
               'analytic_scope': 'Conditional inertial first-variation waveform enclosure for the specified retained model state. Finite arithmetic and internal AI review do not formally certify the analytic PDE inputs or replace external peer review.',
               'limitations': ['No effective nonlinear epsilon or remainder bound.',
                               'No validated finite-amplitude nonlinear waveform.',
                               'No observational scale or cosmological-origin inference.'],
               'zenodo_version_created': False, 'aps_submission': False,
               'evidence_publication': 'This receipt records completed research/site publication. The later evidence-only PR merge and its final readback stay local separately.'}
    add('PUBLICATION_RECEIPT.json', (json.dumps(receipt, indent=2)+'\n').encode())
    require('PUBLICATION_EVIDENCE_MANIFEST.json' not in staged, 'Reserved final evidence manifest already present')
    manifest = {'schema': 'hdblast-inertial-waveform-final-evidence-manifest-v1',
                'initial_delivery_files_preserved': len(initial_files),
                'scope': 'Complete staged byte inventory; no analytic mathematics or future upload is certified.',
                'files': [{'path': name, 'bytes': len(raw), 'sha256': sha(raw)} for name, raw in sorted(staged.items())]}
    manifest_raw = (json.dumps(manifest, indent=2)+'\n').encode()
    staged['PUBLICATION_EVIDENCE_MANIFEST.json'] = manifest_raw
    stage.mkdir(parents=True)
    for name, raw in staged.items():
        target = stage/name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(raw)
    for name, raw in initial_files.items():
        require((stage/name).read_bytes() == raw, 'Initial evidence changed during staging')
    print(json.dumps({'status': 'PASS', 'initial_files_preserved': len(initial_files),
                      'final_evidence_files': len(staged), 'new_files': len(staged)-len(initial_files),
                      'manifest_sha256': sha(manifest_raw), 'stage': str(stage)}, indent=2))


if __name__ == '__main__':
    main()
