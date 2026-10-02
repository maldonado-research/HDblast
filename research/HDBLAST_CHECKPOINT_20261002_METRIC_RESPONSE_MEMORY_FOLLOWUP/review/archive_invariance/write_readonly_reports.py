#!/usr/bin/env python3
"""Render completed, authenticated read-only audit reports. No physics."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def main():
    compare = json.loads((ROOT / 'EXACT_ARRAY_COMPARISON.json').read_text())
    padding = json.loads((ROOT / 'LONGDOUBLE_PADDING_AUDIT.json').read_text())
    audit = json.loads((ROOT / 'COMPLETED_CONTROL_AUDIT.json').read_text())
    bindings = audit['bindings']
    table = '\n'.join('| ' + r['setting'] + ' | ' + str(r['prior_member_count']) + ' | ' + str(sum(f['numeric'] for f in r['fields'])) + ' | ' + str(sum(f['text'] for f in r['fields'])) + ' | ' + str(sum(not f['array_storage_bytes_equal'] for f in r['fields'])) + ' |' for r in compare['archives'])
    archive_pins = '\n\n'.join(f"- `{r['setting']}` prior SHA-256: `{r['prior_sha256']}`\n- `{r['setting']}` streamed-followup SHA-256: `{r['current_sha256']}`" for r in compare['archives'])
    compare_md = f"""# Post-run archival invariance: completed metric memory followup

The storage-only followup preserved every archived positive-source value exactly. Both preserved `positive_B` coarse/fine archives have 631 members: 629 numeric arrays and two Unicode name arrays. All {compare['compared_member_pairs']:,} corresponding members ({compare['compared_elements']:,} elements) have identical dtypes, shapes and ndarray values. There are zero mismatched, missing or extra members and zero nonfinite array pairs. This is a post-run archival/reproducibility result; the registered metric calibration remains **FAIL**.

| Archive | Members | Numeric | Text | Storage-byte differences |
| --- | ---: | ---: | ---: | ---: |
{table}

The old and new ZIP file hashes differ. The array-level storage audit separately inspected every one of the {padding['audited_storage_unequal_arrays']:,} nonidentical byte strings on this little-endian x86-64 ABI, where each 16-byte long-double lane contains a 10-byte binary80 value and six unused padding bytes. All binary80 value bytes are identical. All {padding['unused_padding_bytes_different']:,} differing bytes occur exclusively in unused padding; no signed-zero/value-byte differences were found. Complex long-double values were checked as two consecutive lanes. ZIP hash equality and unused padding equality are not numerical equality gates.

The comparison loaded only one corresponding saved-array pair at a time. The largest pair occupied {compare['maximum_loaded_member_pair_bytes']:,} bytes; observed process peak RSS was {compare['comparison_process_peak_rss_kib']:,} KiB for exact member comparison and {padding['peak_rss_kib']:,} KiB for the padding audit. Python {compare['python']} and NumPy {compare['numpy']} were used. No source, mode equation, physical solver, quadrature or scientific validator was executed.

The prior run completed only the positive-source coarse/fine evolutions before a resource failure. Consequently, this comparison concerns those two historical archives. Signed-source archives in the current run have no corresponding old completed archives. No provenance metadata changes occur within the 631 NPZ members: both text arrays also match exactly. Producer JSON provenance and ZIP/storage hashes are recorded separately.

## Immutable pins and complete evidence

- Third public source freeze: `{bindings['public_freeze_commit']}`.
- Third registration SHA-256: `{bindings['registration_sha256']}`.
- Third independent manifest SHA-256: `{bindings['independent_manifest_sha256']}`.
- Exact member report SHA-256: `{sha(ROOT / 'EXACT_ARRAY_COMPARISON.json')}`.
- Padding report SHA-256: `{sha(ROOT / 'LONGDOUBLE_PADDING_AUDIT.json')}`.
- Prior producer JSON SHA-256: `{compare['prior_producer_report_sha256']}`.
- Third producer JSON SHA-256: `{compare['current_producer_report_sha256']}`.

{archive_pins}

`EXACT_ARRAY_COMPARISON.json` records every member, not only a summary. `LONGDOUBLE_PADDING_AUDIT.json` records every storage-unequal member, including binary80 and padding byte counts. Both reports pin their producing audit sources and retain the scientific **FAIL** status. The sources were prepared/executed outside the checkout after explicit confirmation that the four-run control had completed; these are post-run analyses, not prospective frozen producer inputs.

## Reproduce the read-only comparison

Use the pinned Python interpreter. Both commands refuse to overwrite an existing report. The completion receipt must report `PASS_EXPECTED_SCIENTIFIC_FAILURE_CONTROL` while retaining underlying scientific `FAIL`.

```bash
/workspace/hdblast-cloud-setup/venv-frw/bin/python compare_preserved_metric_archives.py \\
  --prior /workspace/hdblast-research-work/metric-followup-complete-replay-001/fresh/independent \\
  --current /workspace/hdblast-research-work/metric-memory-complete-control-001/replay/fresh/independent \\
  --completion-receipt /workspace/hdblast-research-work/metric-memory-complete-control-001/EXPECTED_SCIENTIFIC_FAILURE.json \\
  --output /tmp/EXACT_ARRAY_COMPARISON.json

/workspace/hdblast-cloud-setup/venv-frw/bin/python audit_longdouble_padding.py \\
  --comparison /tmp/EXACT_ARRAY_COMPARISON.json \\
  --output /tmp/LONGDOUBLE_PADDING_AUDIT.json
```

The required original/current run archives remain at the recorded paths or may be supplied through equivalent preserved directories. A copied exact report contains absolute archive paths; adjust those paths in a separate replay input if moving the archives, and record that adaptation rather than altering the original report.
"""
    failure_table = '\n'.join(f"| `{q}` | {value:.17g} |" for q, value in audit['maximum_observable_refinement_differences'].items())
    outcome_md = f"""# Completed memory followup: scientific failure retained

All four independent evolutions completed within the unchanged resource budgets, and the unchanged internal Ward gates then failed. The expected-failure wrapper successfully authenticated that outcome; the registered metric calibration remains **FAIL**, with cross-route calibration **NOT ESTABLISHED**. This read-only audit verified all 370 frozen replay-copy inputs against their registered SHA-256 and pinned Git-object blobs, verified all 32 captured fresh output hashes, and independently reconstructed the saved gate failure list. It did not query GitHub remote publication or recompute raw physical operators.

The checkpoint concerns homogeneous linear metric response at `H=1`, physical mass squared `r=2`, `xi=0`, fixed physical scalar and incoming BD state. It does not establish an actual shifted-root retarded matrix or coupled stability.

## Completed coverage and resources

- Frozen plan: 32 commands; attempted: 28; successful: 27, comprising 26 pure checks and the primary producer.
- Primary: all 12 source/time rows and 36 finite-cutoff points completed in {audit['primary_elapsed_seconds']:.12g} seconds.
- Independent: `positive_B` and `signed_uB`, each coarse/fine; four full evolutions, six observations and 18 finite-cutoff points per run; all 12 combined source/time rows saved.
- Independent producer elapsed time: {audit['resources']['elapsed_seconds']:.12g} seconds, below the 900-second budget.
- Independent peak RSS: {audit['resources']['peak_rss_kib']:,} KiB, below the 262,144-KiB budget.
- Failing command: `independent_modes`, ordinary exit 1; no timeout or resource-budget failure.
- Not reached after this failure: {', '.join('`' + c + '`' for c in audit['unattempted_commands'])}. No cross-route calibration, summary/report or plot acceptance is asserted.

The prior high-precision round exceeded the same memory budget at 306,088 KiB after completing its two positive-source runs. The current lower recorded peak accompanies completion of all four runs. It is an observed storage-engineering improvement, supported by unchanged numerical-source ASTs and exact archived-value invariance; this audit does not infer a measured allocation profile or a physical change from those resource records.

## Unchanged internal gate failures

Saved-data gate recomputation produces **59 failures**: **29 Ward endpoint** failures against `2e-6`, and **30 Ward refinement** failures against `1e-6`. There are no observable refinement failures. All saved observables are finite; canonical and physical Wronskian maxima are respectively {audit['maximum_wronskian_scaled']:.17g} and {audit['maximum_physical_wronskian_scaled']:.17g}, below the unchanged `1e-10` gate.

The maximum Ward endpoint difference is {audit['maximum_ward_endpoint_difference']:.17g}, including the historical witness `positive_B`, fine, `eta=-1.5`, `K=256` (about {audit['maximum_ward_endpoint_difference']/2e-6:.7g} times its gate). Maximum coarse/fine Ward-ledger difference is {audit['maximum_ward_refinement_difference']:.17g} (about {audit['maximum_ward_refinement_difference']/1e-6:.8g} times its gate). This establishes failure of the registered numerical Ward diagnostic. It does not establish physical instability, heating, an inconsistency of general relativity or evidence for the blast hypothesis. A raw-operator diagnostic and analytic Ward review remain separate work.

| Observable | Maximum saved coarse/fine difference |
| --- | ---: |
{failure_table}

Every one of the 36 combined source/time/cutoff gate points, every failed label, every verified input and every captured output pin is preserved in `COMPLETED_CONTROL_AUDIT.json`. The adjacent archival report documents exact equality of all 1,262 preserved positive-source member pairs and separates unused long-double padding differences from value equality.

## JSON consistency and audit history

The original post-run audit imposed `1e-18` absolute consistency when recomputing a saved Ward endpoint difference from separately serialized binary64 operands. That audit-only check failed: the maximum serialization discrepancy is about `1.054e-17`. Its initial source and failure receipt are retained as `audit_completed_memory_control.initial_too_strict.py.txt` and `READ_ONLY_AUDIT_INITIAL_FAILURE.json`. The final audit uses the frozen expected-failure wrapper's consistency rule (`1e-15` absolute, `1e-12` relative), records each residual and additionally requires that every scientific Ward gate classification agree before and after operand recomputation. No producer, physical value or scientific gate was changed. The retained calibration status remains **FAIL** throughout.

## Pins and reproduction

- Public source freeze: `{bindings['public_freeze_commit']}`.
- Registration SHA-256: `{bindings['registration_sha256']}`.
- Independent manifest SHA-256: `{bindings['independent_manifest_sha256']}`.
- Completed-control audit SHA-256: `{sha(ROOT / 'COMPLETED_CONTROL_AUDIT.json')}`.
- Expected-failure completion receipt SHA-256: `{audit['pins']['completion_receipt_sha256']}`.
- Underlying failed replay ledger SHA-256: `{audit['pins']['underlying_ledger_sha256']}`.
- Independent failed producer JSON SHA-256: `{audit['pins']['producer_results_sha256']}`.

```bash
/workspace/hdblast-cloud-setup/venv-frw/bin/python audit_completed_memory_control.py \\
  --control /workspace/hdblast-research-work/metric-memory-complete-control-001 \\
  --repository /workspace/HDblast \\
  --checkpoint-relative research/HDBLAST_CHECKPOINT_20261002_METRIC_RESPONSE_MEMORY_FOLLOWUP \\
  --output /tmp/COMPLETED_CONTROL_AUDIT.json
```

This audit only streams saved file and Git-object bytes and evaluates arithmetic on already saved JSON values. No new physical evaluation occurs.
"""
    (ROOT / 'POSTRUN_ARCHIVAL_INVARIANCE.md').write_text(compare_md)
    (ROOT / 'COMPLETED_MEMORY_CONTROL_FAILURE.md').write_text(outcome_md)
    execution = {'scope': 'Observed post-run read-only audit execution summaries; not raw stdout/stderr captures.',
                 'new_physical_evaluations': 0, 'commands': [
        {'source': 'compare_preserved_metric_archives.py', 'exit_code': 0, 'report': 'EXACT_ARRAY_COMPARISON.json', 'status': compare['status']},
        {'source': 'audit_longdouble_padding.py', 'exit_code': 0, 'report': 'LONGDOUBLE_PADDING_AUDIT.json', 'status': padding['status']},
        {'source': 'audit_completed_memory_control.py', 'exit_code': 0, 'report': 'COMPLETED_CONTROL_AUDIT.json', 'status': audit['status']}],
                 'initial_read_only_consistency_failure_preserved': 'READ_ONLY_AUDIT_INITIAL_FAILURE.json'}
    (ROOT / 'POSTRUN_EXECUTION_RECEIPT.json').write_text(json.dumps(execution, indent=2, sort_keys=True) + '\n')
    manifest = {'scope': 'Post-run read-only archival invariance and completed scientific-failure audit. Does not alter the 370-file prospective registration.',
                'scientific_status': 'FAIL', 'metric_calibration_passed': False, 'new_physical_evaluations': 0,
                'recorded_utc': datetime.now(timezone.utc).isoformat(),
                'files': {p.name: sha(p) for p in sorted(ROOT.iterdir()) if p.is_file() and p.name != 'MANIFEST.json'}}
    (ROOT / 'MANIFEST.json').write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': 'PASS_POSTRUN_READ_ONLY_PACKAGE', 'files': len(manifest['files']),
                      'manifest_sha256': sha(ROOT / 'MANIFEST.json'), 'directory': str(ROOT)}))

if __name__ == '__main__':
    main()
