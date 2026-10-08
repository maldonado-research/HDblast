#!/usr/bin/env python3
"""Prepare original actual-science prose and statically scoped metadata checks only.

No controller is instantiated. The copied metadata-only validation AST contains
no inventory, journal, transport, state or publication calls.
"""
import ast
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import types

sys.dont_write_bytecode = True
WORK = Path('/workspace/hdblast-research-work/continuation-state-comparison-20261007')
OUT = WORK / 'publication-final/metadata'
BASE = Path('/workspace/hdblast-research-work/continuation-next-20261007/release-planning/FINAL_METADATA_MODERN.json')
BASE_SHA = '72415e84db788b8dec0c421012c1d264fc669e6f8e5d08d93e80f0075e3485ba'
CONTROLLER = WORK / 'publication-planning/new_edition_controller.py'
CONTROLLER_SHA = '27fd844cc903a445688a177339b7af633c0a3581f71eeedf7e20bb015bee36a7'

def pin(path):
    raw = path.read_bytes()
    return {'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()}

def require(condition, reason):
    if not condition:
        raise ValueError(reason)

def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False) + '\n')

def main():
    require(pin(BASE)['sha256'] == BASE_SHA, 'Prior complete metadata changed')
    require(pin(CONTROLLER)['sha256'] == CONTROLLER_SHA, 'Prospective controller changed')
    baseline = json.loads(BASE.read_text())
    result = copy.deepcopy(baseline)
    m = result['metadata']
    m['title'] = 'HDBLAST: Exact Stored-State Comparison with the Certified BD Prehistory Target'
    m['publication_date'] = '2026-10-07'
    m['version'] = '2026.10.07-stored-state-comparison'
    m['description'] = '\n'.join([
        '<p>This scientific checkpoint encloses the difference between exact represented incoming first-order states retained in four source/resolution capsules and the prescribed Bunch-Davies prehistory target. It covers all 49,152 capsule-node occurrences, 20 selected exact-reader arrays and 12 registered finite weighted prefixes at K=64,128,256. All registered stored-minus-target error rectangles are ENCLOSED. Capsule identities remain separate even when their momentum values coincide.</p>',
        '<p>The fixed reader preserves represented binary80 values as exact rationals, including momenta, quadrature weights and u_1/w_1. It uses U=u_1/epsilon_represented and W=w_1/epsilon_represented without projecting or renormalizing the saved state. The preregistered target engine uses global degree 832 and Arb precision 512. Every exported U rectangle and every exported W rectangle passes the complete L1-radius gate of at most 1e-18. This gate controls certified target uncertainty, not incoming mismatch or stress error.</p>',
        '<p>Certified positive lower and upper brackets resolve nonzero incoming mismatch. On the complete K=256 prefixes, coarse-capsule maximum normalized Cartesian complex L1 W errors are approximately 1.91e-14 at retained momenta near 255.86 and 255.90, while fine-capsule maxima are approximately 5.3-5.8e-16. These are decimal display approximations; exact outward brackets and maximizing-node witnesses in RESULTS.md are authoritative. These observations do not establish the cause of the errors or a convergence theorem.</p>',
        '<p>Exact finite-weight transport encloses signed anchor differences and phase-uniform linear density/pressure differences on eta in [-4.5,-3.5], using the represented measure dk*k^2/(2 Pi_represented^2) and the inherited a_0(eta)^4/epsilon_represented scaled R/P convention. Subsequent forcing remains active; the difference evolves homogeneously under identical forcing and complete contacts. These bounds concern the twelve fixed finite represented weighted prefixes. They do not establish continuum interpolation or quadrature, later saved numerical-trajectory accuracy, full source/contact error or ultraviolet completion.</p>',
        '<p>Actual normal and optimized registered executions pass, all eight stable scientific artifacts agree byte for byte, and independent saved-output review passes. Original target-node exports, analytic bounds, source certificates, execution receipts/logs, exact norm witnesses and original figures accompany the results. The scientific source was publicly frozen before execution at Git commit bedcca7e86995da1230c12a20fae3755b31f4a94. The portable delivery carries all eight original/compact input files at their original paths, original mathematical and literature synthesis, and a guarded standalone replay helper. Final package pins and the separate extraction-replay outcome are supplied by the completed immutable delivery record.</p>',
        '<p>The full continuum pressure/contact certificate remains UNRESOLVED, historical metric calibration FAIL, higher-dimensional Big Bang causation NOT_ESTABLISHED and external novelty NOT_ASSESSED. Internal AI-assisted checks are not external peer review or proof-assistant formalization. Literature access comprises 46 recorded TLS GET responses with status 200 and eight selected primary-source passage readings. This is bounded coverage, not an exhaustive search or full-paper/proof verification; private caches and third-party paper bodies are excluded.</p>',
        '<p>This is an additive edition in the existing concept family 10.5281/zenodo.17088132, following record 23225288. All 27 prior files are preserved and the portable comparison archive plus an original overview are the two intended additions. Final file counts, content hashes and delivery status come from the sealed release inventory and complete public content readback.</p>'
    ])
    m['additional_descriptions'][0]['description'] = '\n'.join([
        '<p>Semantic version 2026.10.07-stored-state-comparison has scientific checkpoint date October 7, 2026, Pacific time (2026-10-07). Actual preparation/execution journals use October 8, 2026 UTC timestamps; actual Zenodo publication time is recorded independently.</p>',
        '<p>The immutable scientific source freeze is bedcca7e86995da1230c12a20fae3755b31f4a94, with complete registration SHA256 05e9943f0b4c8134252a2fecef7631ddba8bd398d18553d6e126fbf50fc3aacc and independently verified public-readback GO SHA256 27784552d9d7e10167066927d54f1757351b0f925b499e951b79849a44b62769. Both actual normal/optimized study executions and independent serialized-output review pass; all eight stable scientific artifacts are byte-identical.</p>',
        '<p>ROOT_FINAL_PACKAGE_PINS: PENDING_ROOT_FINAL_SEAL_AND_PORTABLE_REPLAY. Replace this marker only with the actual archive filename/bytes/MD5/SHA256, payload-manifest/helper/replay-pins SHA256, final complete 29-file inventory and successful separately recorded portable extraction-replay evidence.</p>',
        '<p>All 27 inherited files and historical proof/failure evidence remain preserved. The earlier completed publisher state is not reused. The new edition requires complete fresh MD5/SHA256 content streams for all 29 inherited-plus-new objects before and after publication, including objects whose file identities change on version creation. Only client fresh guards and public readback are established; no server-side atomic publication compare-and-swap is claimed.</p>',
        '<p>The previously recorded automatic GitHub-to-Zenodo archiving setting remains historical OFF evidence; no fresh account-side configuration check is claimed. Only original source, results, explanatory synthesis and typed citations are delivered. The full continuum pressure/contact conclusion remains UNRESOLVED, metric calibration FAIL, higher-dimensional origin NOT_ESTABLISHED and external novelty NOT_ASSESSED.</p>',
        '<p>Historical Notes from the prior completed record 23225288 follow verbatim. Their previous versions, asset hashes, packaging failures and subsequent corrections describe that prior edition, separately from this current stored-state study and its forthcoming portable extraction replay.</p>',
        baseline['metadata']['additional_descriptions'][0]['description']
    ])
    m['related_identifiers'].append({'identifier': '10.5281/zenodo.23225288',
        'relation_type': {'id': 'references'}, 'resource_type': {'id': 'software'}, 'scheme': 'doi'})
    OUT.mkdir(parents=True, exist_ok=True)
    target = OUT / 'FINAL_METADATA_MODERN_CANDIDATE.json'
    write_json(target, result)

    spec = importlib.util.spec_from_file_location('metadata_static_controller', CONTROLLER)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    tree = ast.parse(CONTROLLER.read_bytes())
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'Controller')
    method = next(n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == 'validate_inputs')
    start = next(i for i, n in enumerate(method.body)
                 if isinstance(n, ast.Expr) and isinstance(n.value, ast.Call) and
                 ast.unparse(n.value.args[0]).startswith('self.metadata is not None'))
    selected = copy.deepcopy(method.body[start:])
    # The selected segment uses only metadata values and the pinned read-only baseline.
    forbidden_attributes = {'inventory', 'transport', 'state', 'save', 'write', 'event', 'capture', 'get', 'add', 'dir'}
    require(not any(isinstance(n, ast.Attribute) and n.attr in forbidden_attributes
                    and isinstance(n.value, ast.Name) and n.value.id == 'self'
                    for body in selected for n in ast.walk(body)), 'Metadata verifier includes state/network operation')
    fn = ast.FunctionDef(name='metadata_only',
        args=ast.arguments(posonlyargs=[], args=[ast.arg(arg='self')], vararg=None,
                           kwonlyargs=[], kw_defaults=[], kwarg=None, defaults=[]),
        body=selected, decorator_list=[])
    extracted = ast.fix_missing_locations(ast.Module(body=[fn], type_ignores=[]))
    scope = dict(vars(module))
    exec(compile(extracted, str(CONTROLLER) + ':metadata_only_AST', 'exec'), scope)
    published_path = WORK / 'publication-planning/read-only-observation/CURRENT_PRIOR.raw.json'
    published = json.loads(published_path.read_text())

    def check(value):
        scope['metadata_only'](types.SimpleNamespace(metadata=value, published_baseline=published))

    check(result)
    checks = [{'case': 'actual_candidate_complete_metadata_whitelist', 'expected': 'PASS', 'status': 'PASS'}]
    mutations = [
        ('creators_changed', lambda d: d['metadata']['creators'].clear()),
        ('license_changed', lambda d: d['metadata']['rights'].clear()),
        ('historical_notes_changed', lambda d: d['metadata']['additional_descriptions'][1].update(description='changed')),
        ('current_notes_type_changed', lambda d: d['metadata']['additional_descriptions'][0].update(type={'id':'other'})),
        ('old_related_reference_changed', lambda d: d['metadata']['related_identifiers'][0].update(identifier='changed')),
        ('unreviewed_reference_added', lambda d: d['metadata']['related_identifiers'].append({'identifier':'unreviewed'})),
        ('custom_fields_changed', lambda d: d['custom_fields'].clear()),
        ('access_changed', lambda d: d['access'].update(files='restricted')),
        ('metadata_field_removed', lambda d: d['metadata'].pop('copyright')),
        ('server_assigned_field_added', lambda d: d.update(id='999')),
        ('previous_semantic_version_reused', lambda d: d['metadata'].update(version=module.PRIOR_SEMANTIC_VERSION))
    ]
    for name, mutation in mutations:
        changed = copy.deepcopy(result)
        mutation(changed)
        try:
            check(changed)
        except module.Stop as e:
            checks.append({'case': name, 'expected': 'REJECT', 'status': 'PASS', 'rejection': str(e)})
        else:
            raise ValueError('Metadata boundary accepted ' + name)
    editable = {'title', 'description', 'additional_descriptions', 'publication_date', 'version', 'related_identifiers'}
    unchanged = {k: baseline['metadata'][k] == result['metadata'][k]
                 for k in baseline['metadata'] if k not in editable}
    require(all(unchanged.values()) and result['custom_fields'] == baseline['custom_fields']
            and result['access'] == baseline['access'], 'Literal inherited field preservation failed')
    notes_preserved = result['metadata']['additional_descriptions'][1:] == baseline['metadata']['additional_descriptions'][1:]
    refs_preserved = result['metadata']['related_identifiers'][:-1] == baseline['metadata']['related_identifiers']
    require(notes_preserved and refs_preserved, 'Literal historical Notes/references changed')
    proof = {
        'status': 'PASS_ACTUAL_METADATA_WHITELIST_WITH_DELIVERY_PINS_PENDING',
        'optimized': not __debug__, 'sealed_for_live_publication': False,
        'baseline': {'path': str(BASE), **pin(BASE)},
        'published_baseline': {'path': str(published_path), **pin(published_path)},
        'prospective_controller': {'path': str(CONTROLLER), **pin(CONTROLLER)},
        'candidate_metadata': {'path': str(target), **pin(target)},
        'overview': {'path': str(OUT / 'HDBLAST_STORED_STATE_COMPARISON_OVERVIEW_20261007.md'),
                     **pin(OUT / 'HDBLAST_STORED_STATE_COMPARISON_OVERVIEW_20261007.md')},
        'allowed_metadata_fields': sorted(editable),
        'changed_metadata_fields': sorted(k for k in baseline['metadata'] if baseline['metadata'][k] != result['metadata'][k]),
        'unchanged_metadata_fields_exact': unchanged,
        'custom_fields_exact': True, 'access_exact': True,
        'historical_notes_exact': notes_preserved, 'current_notes_type_and_language_exact': True,
        'prior_current_notes_preserved_verbatim_as_labeled_history': result['metadata']['additional_descriptions'][0]['description'].endswith(baseline['metadata']['additional_descriptions'][0]['description']),
        'all_old_related_references_exact': refs_preserved,
        'old_reference_count': len(baseline['metadata']['related_identifiers']),
        'new_reference_count': len(result['metadata']['related_identifiers']),
        'only_appended_reference': result['metadata']['related_identifiers'][-1],
        'checks': checks, 'passes': len(checks), 'skipped': 0,
        'controller_constructed': False, 'controller_state_reused': False,
        'verifier_scope': 'metadata_only_AST_extracted_from_pinned_validate_inputs',
        'network_actions': 0, 'repository_actions': 0, 'retained_array_decodes': 0, 'source_callbacks': 0,
        'root_final_package_pins': 'PENDING_ROOT_FINAL_SEAL_AND_PORTABLE_REPLAY',
        'new_publication_DOI': 'NOT_ASSIGNED',
        'future_asset_bytes_or_hashes_invented': False}
    mode = 'OPTIMIZED' if not __debug__ else 'NORMAL'
    write_json(OUT / ('METADATA_TRANSITION_' + mode + '.json'), proof)
    print(json.dumps({'status':proof['status'], 'mode':mode, 'checks':len(checks),
                      'metadata':pin(target), 'overview':proof['overview']}, indent=2))

if __name__ == '__main__':
    main()
