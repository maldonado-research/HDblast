The verifier reads saved publication evidence only. It makes no network request, imports no numerical package or controller, launches no worker, and changes no remote record. A PASS means the externally authenticated saved witnesses agree. It does not establish a new live remote state, authenticate self-supplied pins, or rerun the scientific certificate.

Copy the original saved JSON bytes and complete controller JSONL journal into a new proof directory outside repository checkouts. Fill `WITNESS_MANIFEST_TEMPLATE.json` with the byte count and SHA256 of every referenced witness file. The `files` map must include all eight witness roles and all nine preserved-record roles; a file can serve multiple roles where the same original capture is appropriate. Relative paths must remain within the proof directory and must contain no symlink. Obtain and record the full manifest hash externally before running the verifier. Preserve the original captures alongside this additive receipt.

The inputs include the frozen 27-file upload inventory and complete modern editable metadata body; the new public record, independent public file-list capture, and latest-version capture; `VERIFIED_PUBLIC.json`; the journal showing the terminal unauthenticated new-record file-list GET followed by all fresh content streams; and the original 23-file complete-stream receipt. Each prior stream may be reused only when the public file UUID and content-version UUID both equal that original streamed identity. Every one of the four additions requires a fresh complete public stream. The receipt's SHA256, size, and MD5 must match the separately pinned upload inventory.

Supply full current record and file-list captures for preserved versions 23114217 and 22347452 in concept family 17088132, and companion 23111008 in family 22922927. Supply their separately authenticated original baseline record captures, including complete metadata and file inventories. A historical summary containing only an ID or file count is insufficient to establish preservation of metadata. The verifier compares complete editable metadata, file membership, sizes, MD5, and exposed immutable file identities. If a baseline was first captured during this release, disclose that observation date rather than describing it as an earlier observation.

Authenticate the verifier source by its separately reviewed SHA256 before execution. Run with the actual published ID obtained from the independently checked publication response/readback; no new ID is embedded in the source:

```bash
python -B verify_publication.py \
  --root /absolute/path/to/new-proof-directory \
  --witness-manifest /absolute/path/to/new-proof-directory/WITNESS_MANIFEST.json \
  --witness-manifest-sha256 EXTERNALLY_RECORDED_FULL_MANIFEST_SHA256 \
  --expected-record-id ACTUAL_VERIFIED_NEW_RECORD_ID \
  --inventory-sha256 a78c80160a48e0171e8bee61e0455aa9b283b8c20e4c17e6e461c5c5fdb45153 \
  --metadata-sha256 72415e84db788b8dec0c421012c1d264fc669e6f8e5d08d93e80f0075e3485ba \
  --output /absolute/path/to/new-proof-directory/OFFLINE_PUBLICATION_READBACK.json
```

The output must be a fresh file. The original prior complete-stream receipt is independently fixed at SHA256 `f3cb8d3c5280fe90307c98988bb2bf97de252bc2bebfbc20aef296612deb96ff`. The edition inventory totals 27 files and 448,387,915 bytes, including all 23 inherited files totaling 440,222,994 bytes. Modern vocabulary display decorations and equivalent HTML character references are normalized only at explicitly whitelisted metadata paths. All editable values and scalar types are retained.

`test_witnesses.py` provides fabricated witness controls, including a deliberately nonexistent test record ID. These tests assert no publication and perform no source calculation. They reject missing/duplicate files, incorrect size/MD5/SHA256 evidence, missing public-stream evidence, changed immutable identities, stale versions, changed metadata/ORCID/references, altered prior records, symlink/path escapes, changed witness bytes, and promotion of the original literal Notes guard. Normal and optimized receipts retain the same checks.

The prior literal Notes guard remains `FAIL_PRESERVED`, even when the new complete metadata contract and scoped scientific replay pass. This verifier does not promote the unresolved twelve-case certificate, original binary80 state-error claim, metric calibration, higher-dimensional Big Bang origin, or external novelty. The replay packaging diagnosis remains additive: the first clean replay preserved its 31-versus-32-check preparation-receipt failure, and the later successful replay used current proof receipts without changing the numerical source or gates.
