# Offline 29-file publication witness verification

This standard-library-only verifier authenticates saved evidence from independently supplied SHA256 pins. It makes no network request and does not establish a fresh remote state. A root-produced inventory or hash manifest is not independent public readback. Actual publication verification remains pending until the controller has captured and completely streamed all 29 files after publication and the root has sealed the corresponding witness packet.

The edition must preserve concept DOI `10.5281/zenodo.17088132`, inherit all 27 files and 448,387,915 bytes from record `23225288`, and add exactly two files. The final upload inventory uses the controller's actual key `expected_total_files_final` with integer value `29`, declares the completed scientific replay contract and contains exact filename, size, MD5 and SHA256 pins. All 29 assets, including inherited files, require the `FRESH_COMPLETE_CONTENT_STREAM` basis. This verifier allows zero inherited-stream reuses.

`VERIFIED_PUBLIC.json` must have status `PASS_COMPLETE_29_FILE_EDITION_READBACK`, the independently expected new record ID, correct concept family, completed publication state, exact count/byte total and the separately pinned inventory/metadata hashes. Every one of its 29 content rows must bind its SHA256 and immutable file/content-version IDs to the direct public `/files` capture. Both the direct published record and latest-version record must match complete editable metadata and inventory and expose the same immutable file identity. Zenodo's embedded public-record entries omit `version_id`; when that optional field is exposed it must agree too. The `/files` response and verification receipt must always provide both identities.

After the terminal unauthenticated direct public file-list GET in the controller journal, there must be exactly one verified complete stream for every file, with exact size, MD5 and SHA256. Content URLs must identify the independently expected record and encoded filename. A new-public file listing, matching checksum string or earlier dated receipt alone cannot substitute for those streams.

The existing `prior_full_stream` witness role retains the dated `PRIOR27_CONTENT_RECEIPT.json` under exact SHA256 `2664ec6a867d806e6146847fcf677d2e5c20061a400262acb8563d3ca552529c`. Its rows combine a prior complete public stream with the later GET-only immutable-identity check; this is historical provenance, not a repeated fresh content stream. It cannot satisfy any of the new 29 stream obligations.

Four earlier records are preserved: `23225288`, `23114217` and `22347452` in concept family `17088132`, plus companion `23111008` in concept family `22922927`. Each witness group must contain exactly `record`, `files`, `baseline` and `baseline_files`. The latter two are the unchanged raw dated captures; their eight SHA256 values are fixed in the verifier. Complete editable metadata, file membership, sizes, checksums and both dated immutable identities must be preserved. Embedded entries are also checked against their endpoint file identity. Supplying a newly substituted baseline and re-pinning its manifest does not replace the fixed historical baseline.

For prior `23225288`, dated raw baselines originate from `publication-planning/read-only-observation/CURRENT_PRIOR.raw.json` and `CURRENT_FILES.raw.json`. The three older records' raw baselines originate from `proof/preserved/{record_id}-current-record.json` and `...-current-files.json` in the previous frozen `2026.10.07-bd-prehistory` publication packet. The source packet's witness manifest remains unchanged.

Metadata comparison preserves exact JSON value types and required editable-key presence. It ignores only whitelisted RDM vocabulary decorations and compares decoded parsed HTML structure, attributes and character references, allowing document-edge ASCII whitespace. Every interior whitespace character, declaration and processing instruction is retained. Malformed closing tags are rejected. This removes broad block-layout assumptions, including cases where a frozen inline style makes a space visibly significant. Dictionary file-list keys must agree with nested filename fields.

The historical original literal Notes guard remains `FAIL_PRESERVED`; successful semantic HTML readback does not rewrite that historical failure. Cached genuine prior metadata still passes the stricter comparison against its separately frozen edition contract.

The first new verifier revision and its initially passing author controls are archived under `revisions/4a030499.../`. Independent review then demonstrated six accepted corruptions: a description string/list type collision, deleted inline interword space, inserted processing instruction or declaration, omitted empty custom fields, and a contradictory dictionary filename. A later preserved revision `35689ea2...` closed those six but still admitted styled-block space loss and malformed closing-tag deletion. Those findings and their source-pinned reproductions remain in `independent-review`. Current controls reject those cases without changing previous genuine publication receipts.

The current verifier is `verify_publication.py`, SHA256 `e2f8be2a1bddbffb928cbf1b6b6b5546d0e33558793aaf2d799fd62cd3f4f02b`. Manufactured controls use `test_manufactured_witnesses.py`, SHA256 `ee076eee032a7710f6fd2dfff606835384e1fed3e31fda8c78f3ec65c27081b8`, and the externally pinned seed manifest SHA256 `b778c6c6583365db2892a722eadb4719f8246572dc20e688e6e1b8417aee2db7`. They pass five positive checks and 117 rejection mutants in each of normal and optimized Python. The invented new record `90000001` and invented additions are explicitly manufactured controls; they are never publication assertions.

Independent final acceptance passes six positive checks and 60 rejection mutants in each mode, independently repeats the complete author suite, and verifies all eight dated baseline raw captures. The accepted receipt is `independent-review/INDEPENDENT_REVIEW_RECEIPT.json`, SHA256 `0d0fd056a3bc48003fb8ee4f1b6fea98ee8f52b7c5342e7cf9d3f6faa5ea9e68`, with no remaining blockers. This accepts offline witness handling; genuine new-publication readback is still a separate required step.

Authenticate the verifier source against its independent external hash before invoking it. After root completes the genuine saved witness packet, invoke it in isolated, standard-library-only Python using separately supplied expected values:

```sh
python -I -S -B verify_publication.py \
  --root /absolute/saved-witness-root \
  --witness-manifest /absolute/saved-witness-root/WITNESS_MANIFEST.json \
  --witness-manifest-sha256 "$WITNESS_MANIFEST_SHA256" \
  --expected-record-id "$EXPECTED_NEW_RECORD_ID" \
  --inventory-sha256 "$INVENTORY_SHA256" \
  --metadata-sha256 "$METADATA_SHA256" \
  --output /absolute/fresh-output/PUBLICATION_normal.json
```

Repeat with `-O` and a separate fresh output file. A successful actual saved-witness result is `PASS_OFFLINE_29_FILE_PUBLICATION_WITNESSES`; preparation controls cannot establish that actual result.

The scientific status field `ENCLOSED_RETAINED_NODES` names the release scope: 49,152 saved capsule-node occurrences, with conditional fixed-weight continuous-time error transport under identical subsequent forcing and contacts. The publication verifier does not itself reprove that mathematical or numerical result. The full twelve-case continuous stress/contact/momentum/UV certificate remains unresolved, metric calibration remains failed, external novelty is unassessed and a higher-dimensional Big Bang origin remains unestablished.
