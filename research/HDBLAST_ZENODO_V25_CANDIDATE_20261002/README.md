# Reviewed HDBLAST v25 publication candidate

This folder prepares a concrete update to the main HDBLAST concept DOI **10.5281/zenodo.17088132**. It contains a metadata request, a complete file inventory, a publication overview and a standard-library helper that reconstructs and verifies the existing archives. It makes no network requests and runs no scientific calculations. The prepared version is **2026.10.02-v25**; a local candidate, a reserved draft DOI or an HTTP 200 response does not establish publication.

The authenticated latest main-family record is [22347452](https://zenodo.org/records/22347452), version 2026.09.05-v24. Its creator, record license and software upload type are preserved exactly. All ten inherited files remain unchanged and are labeled historical v24 material. The nine new attachments retain three failed metric-response experiments, a separate matched-stress positive control, the completed source-free ledger diagnostic, two reports, its fresh-ZIP replay receipt and the version overview. The earlier metric calibration remains **FAIL**.

The scalar-junction/static-branch companion has a separate concept DOI, **10.5281/zenodo.22922927**. Its known draft 23111008 must remain separate; it must not be repurposed as the main HDBLAST update. Earlier metric publication receipts targeted this companion. Their original evidence is retained, while this candidate explicitly selects the main family after an authenticated latest-version check.

The public GitHub integration repair was verified in the user's signed-in browser on October 2, 2026. No GitHub release was created to test it. Enabling public HDblast does not automatically attach a future GitHub-generated Zenodo record to the existing manually created main DOI family. The private archive repository cannot use the public GitHub integration.

## Prepare the nine files locally

Use a clone containing the source artifacts from commit `a8394d0127e58200cd4b9a87c14ae63a9dd69f02` and this candidate folder. The helper verifies every source against FILE_MANIFEST.json, including the three memory-package parts and two ledger-package parts. It writes a fresh directory outside the repository; the directory must not already exist. It copies small files and reconstructs archives without changing their bytes.

```bash
python research/HDBLAST_ZENODO_V25_CANDIDATE_20261002/prepare_candidate.py \
  --repository-root /path/to/HDblast \
  --output-dir /path/to/fresh-prepared-v25 \
  --expected-manifest-sha256 2ac71c9b14d459aa239bc7ed5e7ed04c058d50ed9af97c6781b987bb8a636144
```

The nine Zenodo attachments are listed in the generated PREPARATION_RECEIPT.json. METADATA_REQUEST.json supplies the API/UI metadata and is not an additional Zenodo attachment. The preparation receipt proves local byte integrity, not publication. No raw token is needed for this helper. The manifest excludes local absolute source paths; its relative paths and SHA256 values authenticate the public artifacts.

## Existing-draft browser fallback

The existing main-family draft is [23112891](https://zenodo.org/deposit/23112891). Reuse it for all subsequent attempts. It is authenticated, unsubmitted and belongs to concept 17088132. The new version was created from the latest published record 22347452. Do not create another duplicate draft. The separate companion draft 23111008 remains untouched. A single JSON metadata PUT returned HTTP 500; an authenticated read afterward confirmed the draft was unchanged, with ten inherited files and no new uploads. The candidate metadata was not applied and no publication was requested. [The public status receipt](PUBLICATION_STATUS.json) records these distinctions. Preserve the ten inherited v24 files; do not delete or replace them merely because their titles describe older research.

If authorized API upload transport still produces empty bodies or metadata errors, upload the nine prepared attachments in the signed-in browser to that same main-family draft. Enter the values in METADATA.json and read the saved metadata back. Compare every final remote filename, byte size and checksum against FILE_MANIFEST.json: ten historical files plus nine new files. Do not publish a clone containing only inherited files, empty metadata, missing new attachments or a mismatched DOI family.

Use the existing environment secret for API authentication, with its unchanged network placeholder in the Bearer header and destination zenodo.org. Do not print, log, commit or place the token in prompts, URLs or artifacts. This local helper handles no credentials. Browser upload uses the signed-in account and does not require adding another token.

## Scientific scope

The historical v24 certificate remains conditional on the premises of its fixed-background registered five-field linear model. It does not validate the later metric-response calibration or the cosmological mechanism. The completed twelve-case source-free diagnostic identifies time-integration error in eight cases under its fixed attribution test. It does not establish active-source conservation, complete momentum error, coupled dynamics, stability, heating, a hot Big Bang, extra dimensions or external mathematical novelty. Internal and AI-assisted review is not external peer review.

[Zenodo API documentation](https://developers.zenodo.org/) · [GitHub integration](https://help.zenodo.org/docs/github/enable-repository/) · [Managing versions](https://help.zenodo.org/docs/deposit/manage-versions/)

## Check saved metadata and uploaded files

The offline verifier accepts a raw legacy Zenodo deposition JSON object or a read-only receipt containing that object under `data`. It compares the explicitly selected record ID, main DOI family, saved metadata and all nineteen filenames, byte counts and MD5 checksums. Remote MD5 is compared with the local pinned MD5; no remote SHA256 is claimed.

```bash
python research/HDBLAST_ZENODO_V25_CANDIDATE_20261002/verify_saved_record.py \
  --saved-response /path/to/saved-main-draft.json \
  --expected-record-id 23112891 \
  --expected-manifest-sha256 2ac71c9b14d459aa239bc7ed5e7ed04c058d50ed9af97c6781b987bb8a636144 \
  --output /path/to/new-verification-report.json
```

The current draft fails this guard because it contains only the ten inherited files and its v25 metadata was not saved. This expected failure is recorded in CURRENT_DRAFT_NOT_READY.json. After complete verified uploads, a passing draft check establishes readiness only. To check a later saved published response, add `--published`; the verifier requires `submitted: true`, `state: done`, and a top-level published DOI matching the selected record ID. A metadata-only reserved DOI does not satisfy that check. The tool itself never publishes.

Synthetic guard tests and independent review receipts are included. The local preparation receipt records byte-identical normal and optimized preparation outcomes, separate from the scientific replay already archived.
