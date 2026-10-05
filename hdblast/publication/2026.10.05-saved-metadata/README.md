# Finish the existing HDBLAST draft

The descriptive metadata is saved in owned draft **23114217**, revision **18**.
Fresh authenticated readback confirms the intended title, creator, description,
version, license, subjects, references and related works. The unchanged exact
guard finds one difference: Notes has `ZIP's` where the approved text has `ZIP’s`.
The draft remains unpublished with ten unchanged inherited files; thirteen
approved additions are still missing. The copied preview alone did not establish
publication or a complete upload.

The main concept remains **10.5281/zenodo.17088132**, with published record
[22347452](https://zenodo.org/records/22347452). The separate static companion
[23111008](https://zenodo.org/records/23111008), concept
**10.5281/zenodo.22922927**, is preserved. The removed test record 23112891 stays
removed. Automatic GitHub archiving stays OFF.

## Correct the saved Notes

1. Open [existing draft 23114217](https://zenodo.org/uploads/23114217) in Safari,
   signed in as its owner.
2. Open the [reviewed correction as plain text](https://raw.githubusercontent.com/maldonado-research/HDblast/b15434c6c48ae1b602bba4b5dd550e215e258726/hdblast/publication/2026.10.05-saved-metadata/CORRECT_HDBLAST_NOTES.js)
   and copy all its text. Return to the Zenodo draft tab.
3. Choose **Develop → Show JavaScript Console**, paste and press **Return**.
   Return the resulting popup text to Codex.

After a verified save, click **OK** to reload the editor. This loads the saved
fields into the form. If a save stops, close that editor tab **without Save draft
or Publish**, return the message, and reopen after saved-state reconciliation.
A save may have applied even when its readback fails. Do not repeat an uncertain
write or use the earlier empty-field restoration over this saved metadata.

The correction checks the exact owner, concept, unpublished state and ten-file
inventory. It pins revision 18 and the complete reviewed metadata/custom-field
state, checks the unchanged approved body's SHA256, and permits at most one
conditional metadata PUT followed by fresh readback. A changed draft blocks the
write. Already matching metadata causes only a read and editor reload. It has no
upload, publication, deletion or new-record operation. Session credentials and
the CSRF cookie stay within Zenodo and never enter diagnostics or attempt markers.

**Actual Safari execution of this correction remains unverified.** Its 28 offline
cases and independent internal review are recorded in [OFFLINE_CHECKS.json](OFFLINE_CHECKS.json)
and [PROVENANCE.json](PROVENANCE.json). The frozen candidate and its saved-record
guard are unchanged; the apostrophe difference is recorded, not waived.

## Get the thirteen missing files

Open the [file-download page](https://maldonado-research.github.io/HDblast/zenodo-files/),
click **Prepare download**, keep the tab open, then click **Download ZIP** when
verification finishes. This creates **HDBLAST_ZENODO_UPLOAD_ADDITIONS_23114217.zip**,
about 425 MB, from files already published in the repository. It does not need a
token or a Zenodo login.

Extract this outer delivery ZIP on the Mac. Its `additions` folder contains the
thirteen intended attachments. Keep the research ZIPs inside that folder intact.
After Codex confirms a fresh saved-metadata guard pass, upload those thirteen
files to the existing draft. The outer delivery ZIP itself is not an attachment.
Keep the draft unpublished until fresh complete-file verification passes.

The download checks every transport piece, reconstructed file and the final ZIP.
Its exact size is **424,670,960 bytes**, SHA256:

```text
7b498e7aaa9cd201e2826f06032143e28252424b9e9723268a366fe378016d8a
```

The thirteen additions total **424,668,350 bytes**. Together with the unchanged
ten inherited files, the candidate is **23 files / 440,222,994 bytes**. The newer
5 October mathematical checkpoint stays GitHub-only, outside this frozen
3 October candidate.

Cloud release-asset transport returned `HTTP 400: Bad Content-Length`. The
prepared GitHub delivery release remains an unpublished draft with zero assets;
it is not the download route. No Zenodo file upload or publication occurred in
this recovery round. [DELIVERY_TRANSPORT_RESULT.json](DELIVERY_TRANSPORT_RESULT.json)
preserves the bounded diagnostic without credentials.

## Verification and reproduction

The [saved state](SAVED_STATE.json), [exact metadata check](SAVED_METADATA_CHECK.json)
and [character difference](NOTES_DIFFERENCE.json) record the authenticated
read-only evidence. The current guard result fails only Notes. Require a fresh
`PASS_EXACT_SAVED_DRAFT_METADATA` before uploads, `PASS_COMPLETE_SAVED_DRAFT`
before publication, and authenticated/public post-publication readback before
claiming publication. HTTP 200, a reserved DOI and a browser success popup are
not substitutes for those checks.

The [download provenance](DOWNLOAD_PROVENANCE.json) binds the five static files
and manifest. [Offline controls](DOWNLOAD_OFFLINE_CHECKS.json) pass 24 grouped
cases covering independent SHA vectors, CRC, ZIP framing, all mocked transports,
corruption, truncation, oversize responses, URL restrictions, cancellation and
download gating. The [full local byte pipeline](DOWNLOAD_LOCAL_BUILD.json)
reproduces the entire ZIP and compares its incremental SHA256 with Node crypto.
[Independent Python validation](DOWNLOAD_ZIP_VALIDATION.json) reads every ZIP
member, checks CRC and compares size, MD5 and SHA256 to the frozen candidate.
The [independent manifest comparison](DOWNLOAD_MANIFEST_REVIEW.json) confirms
all thirteen file pins and seventeen ordered source URLs.

[Seventeen live HEAD checks](DOWNLOAD_SOURCE_HEADS.json) establish availability,
exact URLs, CORS `*` and matching advertised sizes; they do not establish live
payload identity. Runtime downloads must check all bytes and hashes.
**Actual Safari download execution remains unverified.** No science replay or
new scientific claim is part of these delivery checks.

From the repository root, with Node 24 and Python 3.11 or later, use new receipt
filenames in an existing writable directory:

```sh
node hdblast/publication/2026.10.05-saved-metadata/test_notes_correction.cjs /tmp/hdblast-notes-new.json
node hdblast/publication/2026.10.05-saved-metadata/test_downloader.mjs /tmp/hdblast-download-new.json
python -B hdblast/publication/2026.10.05-saved-metadata/validate_download_zip.py /path/to/HDBLAST_ZENODO_UPLOAD_ADDITIONS_23114217.zip /tmp/hdblast-zip-new.json
```

The Python command requires the actual downloaded ZIP and runs with assertions
enabled. The [browser recovery workflow](../../../.github/workflows/hdblast_browser_recovery.yml)
repeats the two offline Node suites on relevant changes, with no live Zenodo
writes or large-payload downloads. This is verification automation; no perpetual
research service is running.
