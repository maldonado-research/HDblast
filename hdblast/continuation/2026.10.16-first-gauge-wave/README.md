# First relative-harmonic gauge wave evidence

The dedicated workflow authenticates the complete source at
`research/HDBLAST_FIRST_GAUGE_WAVE_20261016`, its matching ZIP, and a fresh
replay of its explicit scientific roster. It does not infer additional
mathematical claims from the replay count.

`PUBLICATION_PINS.json` records the reviewed frozen source: 67 files, a matching ZIP, and 19 registered replay processes. Fresh extraction reproduced every complete declared report in normal and optimized modes; 64 intentional error/overwrite-rejection runs failed as required. `FRESH_ZIP_REPRODUCTION.json`, `LOCAL_PUBLICATION_VERIFICATION.json`, and `VERIFICATION_SUMMARY.json` record the local results. The reviews are independent internal AI-assisted reviews, with the retained background/master assumptions explicit. The full second-order evolution remains unfinished.

No example or placeholder pin file is supplied. Required pin keys are:

| Key | Required value |
| --- | --- |
| `capsule` | `HDBLAST_FIRST_GAUGE_WAVE_20261016` |
| `manifest_sha256` | SHA-256 of the frozen `MANIFEST.json` |
| `archive_sha256` | SHA-256 of the published ZIP |
| `runner_sha256` | SHA-256 of the frozen `reproduce.py` |
| `archive_bytes` | Positive exact ZIP byte count |
| `source_files_including_manifest` | Positive exact complete source inventory count |
| `expected_replay_processes` | Positive exact count of expanded roster jobs/modes |
| `replay_timeout_seconds` | Optional positive overall timeout; defaults to 900 |

Hashes are complete lowercase hexadecimal strings. The verifier rejects a
missing pin file, zero hashes, malformed counts, incomplete source inventory,
changed files, roster paths outside that inventory, and ZIP/source mismatch
before executing the source. It requires the source runner's
`PASS_COMPLETE_REPORT_AGREEMENT` receipt, the expected process count, and
unchanged original and extracted source afterward.

From the repository root, after real pins are installed, run
`python -I -S -B hdblast/continuation/2026.10.16-first-gauge-wave/verify_and_reproduce.py`
with `--repo-root .`, `--pins` pointing to this directory's concrete pin file,
and `--output-dir` pointing to a new directory outside the repository.
The output contains `VERIFICATION_RECORD.json`, the fresh replay receipt and
checker reports, and runner stdout/stderr. The verifier uses only the Python
standard library and performs no network or Git operation.
