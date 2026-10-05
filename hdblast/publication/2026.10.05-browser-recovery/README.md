# Restore the existing HDBLAST draft's metadata in the signed-in browser

This helper prepares a single conditional metadata save for owned draft
**23114217**, in concept **10.5281/zenodo.17088132**. It uses the exact reviewed
3 October candidate fields, including descriptions, creator, version, license,
44 subjects, seven references, related identifiers and historical custom fields.
It follows the official Invenio editor's same-origin session/CSRF protocol.

**Current status: prepared and offline checked; actual Safari execution is
unverified.** The earlier cloud saves returned HTTP 200 with empty readback.
This route runs directly in the user's signed-in Zenodo browser. Publication
remains pending saved-metadata and complete 23-file verification.

## One browser action

1. Open [draft 23114217](https://zenodo.org/uploads/23114217) in Safari, signed in
   as the owner. Keep that draft tab available.
2. Open the [script as plain text](https://raw.githubusercontent.com/maldonado-research/HDblast/19f28373b5dc87a384ebd6eebe22cbb8c07852b6/hdblast/publication/2026.10.05-browser-recovery/RESTORE_HDBLAST_METADATA.js)
   and copy all its text. Return to the Safari Zenodo draft tab.
3. Choose **Develop → Show JavaScript Console**, paste the copied script and
   press **Return**. The screenshot supplied by the user already shows Safari's
   Develop menu. If that command is unavailable, continue with the
   [manual field guide](../2026.10.03-verified-source-operator/browser/METADATA_FIELDS.md).

If the helper reports saved metadata, click **OK** to refresh the editor and
return to Codex with the result. Refreshing loads the saved fields into the form;
it prevents the old blank form from overwriting the recovered values. Codex must
then authenticate a fresh GET and run the existing saved-record guard before
any of the thirteen additions are uploaded.

If a write stops, close the Zenodo editor tab **without clicking Save draft or
Publish** and return the message to Codex. A save may have applied while the form
still shows the old fields. Reopen the editor after the saved state is reconciled.
A recorded write must be reconciled before another attempt. The helper records a small credential-free attempt marker
in this tab's session storage so another paste cannot repeat an uncertain save.
If the metadata already matches, it performs only a read and refreshes the form.

## Scope and checks

The helper checks the exact origin/editor path, owner, parent concept, unpublished
state, revision and all ten inherited filenames, byte sizes and MD5 checksums.
It verifies the embedded approved request's SHA256 before any network request.
Unreviewed existing metadata blocks replacement. One conditional PUT is followed
by a new GET; all approved metadata and custom fields must match, the inherited
inventory must remain intact and reserved identifiers must remain unchanged.
HTTP 200 alone is insufficient.

The only API URL is the existing draft's same-origin endpoint. Browser session
credentials stay within Zenodo. The CSRF cookie is used in the standard header
and excluded from every displayed or stored diagnostic. The script contains no
personal access token, file upload, publication, deletion, new-record operation
or external script loader. It does not change automatic GitHub archiving.

The ten existing files are retained. The reviewed full candidate has
**23 files / 440,222,994 bytes** and thirteen additions totaling 424,668,350 bytes.
The separate static companion is preserved. The newer 5 October mathematical
checkpoint remains GitHub-only and is not silently added to this candidate.

## Offline validation and provenance

The plain-text execution link above is pinned to the reviewed code commit.

Run `node test_browser_restore.cjs /path/to/a-new-receipt.json` from this directory,
using an existing writable parent directory and a new output filename.
The output file is created once; preserve existing receipts. The [23-case receipt](OFFLINE_CHECKS.json) reports synthetic
browser control flow only: **zero live requests or live metadata writes**.
It exercises successful/reordered readback, already-saved and partial states,
wrong origin/path/owner/family, published state, changed checksums/membership,
unreviewed fields, repeated attempts, blank HTTP 200 readback, conflicts, uncertain
writes, identifier changes, unexpected URLs and body tampering. It also checks
that the synthetic CSRF value never enters diagnostic output or session storage.

The [provenance](PROVENANCE.json) pins the exact approved body, frozen inherited
manifest, helper bytes and official client source blobs. These checks are not a
live browser test or a Zenodo publication receipt. The existing
[saved-record guard](../2026.10.03-verified-source-operator/verify_saved_record.py)
remains the independent metadata/full-file/publication gate.
