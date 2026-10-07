# Verified HDBLAST publication — 7 October 2026, Pacific

The existing same-family record [23114217](https://zenodo.org/records/23114217) is now published, independently verified on **7 October 2026, Pacific time**, with **23 files / 440,222,994 bytes**. Version DOI: **10.5281/zenodo.23114217**; concept DOI: **10.5281/zenodo.17088132**. Its semantic version **2026.10.03-verified-source-operator** and declared publication date **2026-10-03** retain the frozen packet's labels. Authenticated and direct public GETs returned HTTP 200 and published state; every attachment was independently streamed before publication and matched its frozen size, MD5 and SHA256.

The action completed existing record 23114217 in the established main family. [PUBLICATION_RESULT.json](PUBLICATION_RESULT.json) records one publish attempt, HTTP 202, followed by independent authenticated/public GET 200 witnesses and the unchanged guard's `PASS_VERIFIED_PUBLISHED_RECORD` against the reviewed request. It records no additional upload, standalone record or integration-setting change. The declared metadata date remains 3 October; the actual publication/readback occurred on 7 October.

## Metadata representation is separately reviewed

Saved Notes use exactly one `&rsquo;` entity for U+2019 and four LF separators between five paragraphs. Independent decoding proves the approved Unicode paragraph text, order and element/attribute structure are identical. The original frozen literal metadata comparison still **FAILS on Notes only** and its request, guard and failure receipt are preserved. A separately pinned reviewed request changes only that exact HTML serialization; the unchanged guard passes against it, including authenticated and direct public publication witnesses. Raw Notes HTML and its outside whitespace are not claimed byte-identical.

Only `/metadata/notes` differs between the frozen and reviewed JSON requests. The [reviewed request](REVIEWED_HTML_SERIALIZATION_METADATA_REQUEST.json) SHA256 is `66fea4d6a5666b441963e63372977221917dd3239335b660887fea8869821cb3`; the [semantic certificate](SEMANTIC_SERIALIZATION_CERTIFICATE.json) SHA256 is `5ab8bffd7efa67e57388f381dec894f50dcde20afa5b17ec9e475e2d1e1212ab`. All other metadata fields, owner, identity, family and twenty-three file pins keep the strict original checks. No general HTML normalization is allowed.

The [original frozen post-publication result](ORIGINAL_FROZEN_PUBLISHED_RECORD_GUARD.json) is `FAIL_SAVED_RECORD_VALIDATION`, with Notes as its only metadata mismatch. The [reviewed unchanged-guard result](PUBLISHED_RECORD_GUARD.json) is `PASS_VERIFIED_PUBLISHED_RECORD`. [Independent addendum verification](INDEPENDENT_LIVE_REVIEWED_PUBLISHED_CHECK.json) records `PASS_REVIEWED_HTML_SERIALIZATION_PUBLISHED_RECORD`. These are separate results, and the literal frozen gate is not relabelled as passed.

Primary offline controls pass 27 cases and independent review checks 91 cases, including rejected alternate entities, apostrophes, separators, metadata/file/identity/state changes and altered/missing public proof. [Controls](OFFLINE_HTML_SERIALIZATION_CONTROLS.json) and [root review](ROOT_NOTES_SERIALIZATION_REVIEW.json) establish the narrow verifier rules. Live publication is established by the separate witnesses and recorded action/readback, not those offline controls.

## Attachment identity and scientific scope

[All twenty-three streamed checks](ALL23_CONTENT_CHECKSUMS.json) verify sizes, MD5 and SHA256 directly against the original frozen manifest before publication. [Inventory confirmation](CONTENT_AND_INVENTORY_CONFIRMATION.json) shows remote identity and inventory stayed unchanged during those streams. The ten inherited files were reread in this completion, alongside the thirteen additions. The [post-publication content-identity receipt](PUBLISHED_CONTENT_IDENTITY_CONFIRMATION.json) confirms that all twenty-three public file IDs and version IDs match those fully streamed payloads, with unchanged sizes and MD5. A second post-publication payload download is not claimed.

The [authenticated witness](AUTHENTICATED_PUBLISHED_GET.json) and [direct public witness](DIRECT_PUBLIC_PUBLISHED_GET.json) are minimal credential-free captures. No token, request authorization value or upload URL is included. The frozen manifest pin remains `00f0ed73c27e846603680fdc52667392935e59bbe18ac942f30e8e66248937bb`. The [independent publication review](INDEPENDENT_PUBLICATION_RECEIPT.json) confirms main-family publication and the retained previous v24 22347452 inventory/version. The [fresh public companion check](INDEPENDENT_FRESH_PUBLIC_COMPANION_CHECK.json) confirms its complete public response remains identical; companion 23111008/concept 22922927 is preserved. The previous main edition's latest-version relationship can change after publication, so its complete response is not asserted unchanged.

This edition archives the frozen 3 October conditional source/operator certificate and preceding ledger diagnostics. The later 5 October finite-band proof and 7 October incoming-state theorem remain GitHub-only. The actual incoming state remains NOT_ENCLOSED, full twelve-case certificate UNRESOLVED, metric calibration FAIL and higher-dimensional Big Bang origin NOT_ESTABLISHED. Internal exact verification supplies no external peer review or novelty assessment.

Automatic GitHub-to-Zenodo archiving stays OFF. No integration setting was changed; its prior OFF observation is retained without claiming a fresh account-side inspection. Website deployment and final PR integration require their own verification.

## Reproduce the published metadata checks locally

The copied serialization verifier and frozen guard retain their original bytes. The portable wrapper changes only the packet input directory; all internal source/request/manifest pins and acceptance rules remain enforced. From the repository root, with a new output filename:

```sh
python -B hdblast/publication/2026.10.07-verified-publication/run_verification.py --saved-response hdblast/publication/2026.10.07-verified-publication/AUTHENTICATED_PUBLISHED_GET.json --public-response hdblast/publication/2026.10.07-verified-publication/DIRECT_PUBLIC_PUBLISHED_GET.json --output /tmp/hdblast-published-new.json
```

Require `PASS_REVIEWED_HTML_SERIALIZATION_PUBLISHED_RECORD`, with the original literal failure disclosed. Use honestly captured fresh witnesses for a new live-state check; these preserved captures establish the dated readback. The wrapper performs no network requests or remote mutations. The proof manifest authenticates this additive package; the old preparation packets remain immutable.
