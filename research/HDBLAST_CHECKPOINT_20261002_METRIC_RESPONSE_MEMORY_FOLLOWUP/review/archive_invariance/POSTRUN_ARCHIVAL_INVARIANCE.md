# Post-run archival invariance: completed metric memory followup

The storage-only followup preserved every archived positive-source value exactly. Both preserved `positive_B` coarse/fine archives have 631 members: 629 numeric arrays and two Unicode name arrays. All 1,262 corresponding members (14,151,140 elements) have identical dtypes, shapes and ndarray values. There are zero mismatched, missing or extra members and zero nonfinite array pairs. This is a post-run archival/reproducibility result; the registered metric calibration remains **FAIL**.

| Archive | Members | Numeric | Text | Storage-byte differences |
| --- | ---: | ---: | ---: | ---: |
| coarse | 631 | 629 | 2 | 561 |
| fine | 631 | 629 | 2 | 494 |

The old and new ZIP file hashes differ. The array-level storage audit separately inspected every one of the 1,055 nonidentical byte strings on this little-endian x86-64 ABI, where each 16-byte long-double lane contains a 10-byte binary80 value and six unused padding bytes. All binary80 value bytes are identical. All 15,872,188 differing bytes occur exclusively in unused padding; no signed-zero/value-byte differences were found. Complex long-double values were checked as two consecutive lanes. ZIP hash equality and unused padding equality are not numerical equality gates.

The comparison loaded only one corresponding saved-array pair at a time. The largest pair occupied 1,048,576 bytes; observed process peak RSS was 32,492 KiB for exact member comparison and 35,440 KiB for the padding audit. Python 3.12.14 and NumPy 2.2.6 were used. No source, mode equation, physical solver, quadrature or scientific validator was executed.

The prior run completed only the positive-source coarse/fine evolutions before a resource failure. Consequently, this comparison concerns those two historical archives. Signed-source archives in the current run have no corresponding old completed archives. No provenance metadata changes occur within the 631 NPZ members: both text arrays also match exactly. Producer JSON provenance and ZIP/storage hashes are recorded separately.

## Immutable pins and complete evidence

- Third public source freeze: `19fde76912af6e2f30d6f55e066b27d88cc7e34b`.
- Third registration SHA-256: `1c5bc9b21b34b1d036b7a7d12ccdb6766ac2ba7835ea04182dba868843866607`.
- Third independent manifest SHA-256: `71830eb8aaf31bccd2cb628efeb9735490e41a7d7e5b5bffc527177f7c70e075`.
- Exact member report SHA-256: `b13c55607a47e8795d6cf198866ffb70729cb36cb492f37991224f35bb46ec02`.
- Padding report SHA-256: `93ad07f04cb52400120d180abfd73b3e62935cb8075d13149e1719118962ba51`.
- Prior producer JSON SHA-256: `0acc6c57f63327fc9999d77649c37ecef025026a0e225ebc90b2efc4a322b2a0`.
- Third producer JSON SHA-256: `6ffd98ab0a18ba954cee7044d7d4f0581721eb8da6085c949be51f3785e16785`.

- `coarse` prior SHA-256: `0b366d473a672fa110a499759aa61fb8e879b36a7104f9049fee1c5147337102`
- `coarse` streamed-followup SHA-256: `978e7b985ea67644563649d6617884daafeda67046209e7f9ccabe79d5259ea6`

- `fine` prior SHA-256: `8bd8daac59d12cf0af5acfc7ff33185e0a18bd6910a27c7e30963fea1c89ac35`
- `fine` streamed-followup SHA-256: `9f76ff2b57a33b1341f6763a71d10bda8f22fa86fb2d5f7968dffbb716cc517d`

`EXACT_ARRAY_COMPARISON.json` records every member, not only a summary. `LONGDOUBLE_PADDING_AUDIT.json` records every storage-unequal member, including binary80 and padding byte counts. Both reports pin their producing audit sources and retain the scientific **FAIL** status. The sources were prepared/executed outside the checkout after explicit confirmation that the four-run control had completed; these are post-run analyses, not prospective frozen producer inputs.

## Reproduce the read-only comparison

Use the pinned Python interpreter. Both commands refuse to overwrite an existing report. The completion receipt must report `PASS_EXPECTED_SCIENTIFIC_FAILURE_CONTROL` while retaining underlying scientific `FAIL`.

```bash
/workspace/hdblast-cloud-setup/venv-frw/bin/python compare_preserved_metric_archives.py \
  --prior /workspace/hdblast-research-work/metric-followup-complete-replay-001/fresh/independent \
  --current /workspace/hdblast-research-work/metric-memory-complete-control-001/replay/fresh/independent \
  --completion-receipt /workspace/hdblast-research-work/metric-memory-complete-control-001/EXPECTED_SCIENTIFIC_FAILURE.json \
  --output /tmp/EXACT_ARRAY_COMPARISON.json

/workspace/hdblast-cloud-setup/venv-frw/bin/python audit_longdouble_padding.py \
  --comparison /tmp/EXACT_ARRAY_COMPARISON.json \
  --output /tmp/LONGDOUBLE_PADDING_AUDIT.json
```

The required original/current run archives remain at the recorded paths or may be supplied through equivalent preserved directories. A copied exact report contains absolute archive paths; adjust those paths in a separate replay input if moving the archives, and record that adaptation rather than altering the original report.
