#!/usr/bin/env python3
"""Stdlib-only supplemental inspection of the approved unsealed refresh."""
import ast
import importlib.util
import json
from pathlib import Path
from datetime import datetime, timezone

W = Path('/workspace/hdblast-research-work/continuation-state-comparison-20261007')
OUT = W / 'replay-delivery/package-stage-independent-review'
spec = importlib.util.spec_from_file_location('original_static_byte_reviewer', OUT / 'review_unsealed_stage.py')
check = importlib.util.module_from_spec(spec)
spec.loader.exec_module(check)  # Original reviewer code only, outside P; no main().
P = check.P
STAGE = check.STAGE
GENERATOR_SHA = '87983dd395aad29d723a9c27f6287c4b138150544b9a9a8bcfc03070a375236f'
STAGE_SHA = 'f27160a4c97ae71fbe992d804e948486f98c45ee7c0af1477c3e666aa7c79b50'


def main():
    stage_pin = check.opaque_pin(STAGE)
    check.require(stage_pin['sha256'] == STAGE_SHA, 'refreshed staging inventory differs')
    stage = check.load_json(STAGE)
    check.require(stage['sealed'] is False and stage['status'] == 'UNSEALED_PORTABLE_PACKAGE_STAGED', 'refreshed stage is not unsealed')
    old = check.load_json(OUT / 'OBSERVED_INITIAL_SNAPSHOT.json')['files']
    expected = stage['files']
    figure_names = ['EXACT_TABLE.csv', 'FIGURE_RECEIPT.json', 'stored_state_bounds.pdf', 'stored_state_bounds.png']
    added_expected = ['figures/' + n for n in figure_names] + [check.C_REL + '/figures/' + n for n in figure_names] + ['figures/generate_figures.py']
    added = sorted(set(expected) - set(old))
    changed = sorted(n for n in set(expected) & set(old) if expected[n] != old[n])
    check.require(added == sorted(added_expected) and changed == ['README.md'] and not set(old) - set(expected), 'refresh has unexpected map changes')
    # Complete opaque reauthentication is inexpensive locally and confirms retained
    # baseline hashes; no archive or scientific JSON numeric content is parsed.
    observed, dirs = check.tree(P)
    check.require(observed == expected, 'refreshed P differs from external leaf inventory')
    check.require(check.summary(observed) == {'file_count': 236, 'total_leaf_bytes': 257974595}, 'refreshed counts differ')
    check.require(not (P / 'PAYLOAD_MANIFEST.json').exists() and not (P / 'REPLAY_PINS.json').exists(), 'final controls present before root seal')
    initial_receipt = check.load_json(OUT / 'INITIAL_STAGE_REVIEW_RECEIPT.json')
    archived_stage = W / 'package-candidate/STAGING_INVENTORY_ee321df46a117a1f.json'
    check.require(check.opaque_pin(archived_stage) == initial_receipt['staging_inventory'], 'original stage inventory was not preserved')
    readme_receipt = check.load_json(OUT / 'UPDATED_README_STATIC_REVIEW_RECEIPT.json')
    check.require(observed['README.md'] == readme_receipt['readme_source_pin'] == check.opaque_pin(W / 'replay-delivery/PACKAGE_README_UNSEALED.md'), 'refreshed README differs from accepted static text')
    generator = P / 'figures/generate_figures.py'
    check.require(observed['figures/generate_figures.py']['sha256'] == GENERATOR_SHA, 'approved generator identity differs')
    check.require(check.opaque_pin(W / 'root/generate_figures.py') == observed['figures/generate_figures.py'], 'approved source generator bytes differ')
    generator_ast = ast.parse(generator.read_text(encoding='utf-8'), filename='figures/generate_figures.py')
    imports = []
    for node in ast.walk(generator_ast):
        if isinstance(node, ast.Import):
            imports += [n.name for n in node.names]
        elif isinstance(node, ast.ImportFrom):
            imports.append(node.module)
    check.require(set(imports) == {'fractions', 'pathlib', 'csv', 'argparse', 'hashlib', 'json', 'matplotlib', 'matplotlib.pyplot'}, 'approved generator imports unexpected code')
    sources = {}
    for name in added:
        row = stage['copied_sources'][name]
        check.require(set(row) == {'bytes', 'sha256', 'source'}, 'refresh source pin shape differs')
        pin = {'bytes': row['bytes'], 'sha256': row['sha256']}
        check.require(observed[name] == pin == check.opaque_pin(Path(row['source'])), 'approved addition differs from pinned source: ' + name)
        sources[name] = row
    for name in figure_names:
        top_pin = observed['figures/' + name]
        check.require(top_pin == observed[check.C_REL + '/figures/' + name] == check.opaque_pin(W / 'figures-portable-002' / name), 'top/C/final-generation figure bytes differ: ' + name)
    receipt = check.load_json(P / 'figures/FIGURE_RECEIPT.json')
    check.require(receipt['status'] == 'PASS_PLOT_OF_CERTIFIED_JSON_SUMMARIES', 'figure status differs')
    check.require(receipt['generator_sha256'] == GENERATOR_SHA and receipt['input_DATA_sha256'] == observed['evidence/actual/normal/DATA.json']['sha256'], 'figure provenance links differ')
    check.require(type(receipt['source_callbacks']) is int and receipt['source_callbacks'] == 0 and type(receipt['original_array_decodes']) is int and receipt['original_array_decodes'] == 0, 'figure typed-zero counters differ')
    check.require(receipt['continuous_momentum_or_full_stress_certificate'] == 'UNRESOLVED', 'figure scope differs')
    check.require('a0^4/epsilon' in receipt['units'] and '[-9/2,-7/2]' in receipt['time_interval'], 'figure units/time scope differ')
    # Original generated PDF duplicated at two package paths, not a primary-paper body.
    pdf_paths = sorted(n for n in observed if n.lower().endswith('.pdf'))
    check.require(pdf_paths == sorted(['figures/stored_state_bounds.pdf', check.C_REL + '/figures/stored_state_bounds.pdf']), 'unexpected PDF body addition')
    c_source, _ = check.tree(check.SOURCE_C)
    c_copy = check.subset(observed, check.C_REL)
    check.require(c_copy == c_source, 'refreshed complete C copy differs from current repository source')
    check.require(check.opaque_pin(STAGE) == stage_pin, 'refreshed stage inventory changed during review')
    again, again_dirs = check.tree(P)
    check.require(again == observed and again_dirs == dirs, 'refreshed package changed during review')
    snapshot = {'schema_version': 1, 'status': 'INDEPENDENT_REFRESHED_UNSEALED_SNAPSHOT', 'sealed': False, 'package_path': str(P), 'files': observed, 'directories': dirs, **check.summary(observed)}
    check.json_write(OUT / 'OBSERVED_REFRESHED_SNAPSHOT.json', snapshot)
    report = {
        'schema_version': 1,
        'status': 'PASS_REFRESHED_UNSEALED_STAGE_ADDITIONS_AND_BYTE_CLOSURE_REVIEW',
        'sealed': False,
        'utc': datetime.now(timezone.utc).isoformat(),
        'package_path': str(P),
        'staging_inventory': stage_pin,
        'original_stage_inventory_preserved': check.opaque_pin(archived_stage),
        'initial_stage_receipt': check.opaque_pin(OUT / 'INITIAL_STAGE_REVIEW_RECEIPT.json'),
        'updated_readme_receipt': check.opaque_pin(OUT / 'UPDATED_README_STATIC_REVIEW_RECEIPT.json'),
        'observed_snapshot': check.opaque_pin(OUT / 'OBSERVED_REFRESHED_SNAPSHOT.json'),
        **check.summary(observed),
        'all_unchanged_initial_leaf_hashes_reverified': True,
        'added_files': sources,
        'modified_initial_files': changed,
        'removed_initial_files': [],
        'approved_generator': observed['figures/generate_figures.py'],
        'approved_generator_AST': 'PASS_STATIC_ONLY_NOT_EXECUTED',
        'figure_receipt': observed['figures/FIGURE_RECEIPT.json'],
        'figure_generator_and_actual_DATA_links': 'PASS_OPAQUE_HASHES',
        'figure_original_decode_and_source_callback_counters': 'PASS_STRICT_INTEGER_ZERO',
        'figure_exact_table_or_plot_values_numeric_review': 'NOT_PERFORMED_BY_THIS_REVIEWER',
        'top_and_C_figures_exact_identity': 'PASS_FOUR_FILES',
        'complete_C_current_source_identity': check.summary(c_copy),
        'generated_original_PDF_paths': pdf_paths,
        'README_scope_and_generic_commands': 'PASS_ACCEPTED_UPDATED_TEXT',
        'unsafe_paths_links_specials_empty_dirs_or_caches': 0,
        'RESULTS': 'PENDING_ROOT_SUPPLEMENTAL_REPORT',
        'final_payload_manifest_exists': False,
        'final_replay_pins_exists': False,
        'fresh_portable_replay': 'NOT_RUN_IN_THIS_REVIEW',
        'retained_array_decodes': 0,
        'scientific_output_numeric_reads': 0,
        'study_imports': 0,
        'source_callbacks': 0,
        'target_evaluations': 0,
        'generator_or_archived_script_executions': 0,
        'network_calls': 0,
        'repository_writes': 0,
        'package_writes': 0,
        'limitations': ['Supplemental review authenticates artifact bytes and static provenance; it does not independently recalculate figure/table values or visually evaluate the plots.', 'Root RESULT/supplemental reports, final seal, portable replay and publication remain outside this receipt.'],
        'review_script': check.opaque_pin(Path(__file__)),
        'shared_original_reviewer_script': check.opaque_pin(OUT / 'review_unsealed_stage.py'),
    }
    check.json_write(OUT / 'REFRESHED_STAGE_REVIEW_RECEIPT.json', report)
    print(json.dumps({'status': report['status'], 'files': len(observed), 'bytes': report['total_leaf_bytes'], 'added': len(added), 'changed': changed, 'current_C': check.summary(c_copy), 'receipt': check.opaque_pin(OUT / 'REFRESHED_STAGE_REVIEW_RECEIPT.json')}, sort_keys=True))


if __name__ == '__main__':
    main()
