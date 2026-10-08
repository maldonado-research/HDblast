# Independent publication verifier review

Accepted for the declared offline, saved-witness contract at verifier SHA256
`e2f8be2a1bddbffb928cbf1b6b6b5546d0e33558793aaf2d799fd62cd3f4f02b`.
No remaining blocker was identified within this review. This is a manufactured
control and source review, not a claim of a new publication or fresh remote
readback. The original verifier and all diagnosed revisions remain evidence;
their earlier receipts do not authorize use of those revisions.

The final independent controls passed **6 positive cases and 60 rejection
mutants**, separately under normal Python and `python -O`. A separate repeat of
the author's final suite passed **5 positive cases and 117 rejection mutants**
in each mode. The independent suite builds its own complete 29-file fixtures;
it does not import or reuse the author's fixture generator. It includes exact
JSON type preservation, file-ID aliases and conflicting aliases, optional
embedded version IDs, schema and metadata omissions, false fresh-stream bases,
terminal stream evidence, and coordinated substitutions of both dated and
current witnesses.

The final source requires the `PASS_COMPLETE_29_FILE_EDITION_READBACK` receipt,
27 inherited plus 2 new files, the inherited byte total `448387915`, all 29
`FRESH_COMPLETE_CONTENT_STREAM` basis rows, and all 29 terminal journal streams
with exact byte count, MD5 and SHA256. Each receipt's file and content-version
IDs must agree with `/files`. Public record and latest embeddings must expose
the matching file ID through `file_id` or `id`; exposed version IDs must agree.
The positive controls also establish that equal historical inherited IDs are
allowed when every file still has a fresh complete stream. Historical receipt
reuse does not satisfy the fresh-stream requirement.

All four preservation groups require `record`, `files`, `baseline` and
`baseline_files`. Both dated identities must agree with their current listing;
record embeddings remain linked to those identities. Eight fixed raw baseline
pins agree byte for byte with the original captures and with the externally
pinned seed manifest. The prior 27-file receipt retains its exact historical
pin and `PASS_PRIOR_COMPLETE_BYTES_WITH_FRESH_IMMUTABLE_IDENTITY` row status. It
is a historical complete-byte plus fresh-GET-identity composite, not permission
to reuse old SHA256 streams for the new 29-file edition.

Descriptions retain exact JSON types, all interior whitespace, declarations
and processing instructions. Parsed HTML structure, attributes and character
references may normalize; document-edge ASCII whitespace may normalize.
Malformed closing tags are rejected. Editable metadata, custom fields and
access require their complete keys. Dictionary filename keys must agree with
nested names, including agreement between `key` and `filename` when both are
present. Fixed dated baseline pins prevent a caller from replacing both
current and baseline history while merely recomputing a witness-manifest pin.

Earlier independent findings were reproduced against complete manufactured
witnesses before correction:

| Source revision | Independent outcome | Disposition |
| --- | --- | --- |
| `4a030499…` | 5 positives; 44 rejects; 6 accepted corruptions, in both modes | Superseded: description event-array type alias; deleted visible inline spacing; ignored declaration/PI; missing empty companion custom fields; contradictory filename map |
| `35689ea2…` | Expanded suite: 6 positives; 57 rejects; 3 accepted corruptions | Superseded: ignored empty end tag, ignored end-tag junk, lost spacing with inline-styled block tags |
| `e2f8be2a…` | 6 positives; 60 rejects; 0 accepted corruptions, in both modes | Accepted for this offline witness contract |

The first `CONTROLS_original_4a030499_normal.json` receipt included a creator
reorder that was a no-op because the seed has one creator. Its exact script is
preserved as `independent_controls_superseded_initial_noop.py`. The valid first
revision results are the `*_002` normal and optimized receipts, pinned to
`independent_controls_original_findings.py`, which use a meaningful creator-add
mutant. Each other historical script and source snapshot is retained with its
matching receipt. Intermediate passing receipts only cover the case set named
by their script; the expanded cases subsequently found the remaining gaps.

This review performed no live network calls, remote writes, publisher imports,
array decoding, source callbacks, target evaluations or physical calculations.
The verifier imports only standard-library modules. The controls use dated
record/file-list witnesses and fabricated additions; their record ID
`90000007` is explicitly manufactured. Repository and verifier source files
were not edited by the independent reviewer.

The output's `ENCLOSED_RETAINED_NODES` state-error scope is a separately pinned
release contract, not numerical reproof by this verifier. The full twelve-case
certificate remains `UNRESOLVED`, metric calibration `FAIL`, higher-dimensional
Big Bang origin `NOT_ESTABLISHED`, and external novelty `NOT_ASSESSED`.
Externally trusted witness pins authenticate saved evidence; this verifier
does not authenticate the present remote state. Actual publication acceptance
still requires genuine new-edition witnesses, separately trusted pins and
successful verification of all required fresh streams.

Evidence is indexed and SHA256-pinned by `INDEPENDENT_REVIEW_RECEIPT.json`.
