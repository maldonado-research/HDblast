# HDBLAST pending v25: complete local publication candidate

This derived candidate retains all ten historical v24 files and all nine
previously reviewed v25 additions, and adds the completed active-source ZIP,
measured report and external fresh replay receipt. The expected saved Zenodo
inventory is exactly **22 unique files**. The original source-free candidate and
its metadata remain preserved under `source_free_snapshot`; the historical
22347452/latest and 23112891/draft targets are snapshots and no longer current.
Retargeting requires authenticated same-family latest/new-draft GETs with exact
owner1386319, all ten inherited filename/size/MD5 pins, and unchanged creators
and license. The original
candidate folder is untouched. The prior overview remains a timestamped
source-free subset guide. Start with the new active-source report and this
candidate's metadata for the complete update.

The canonical main concept is **10.5281/zenodo.17088132**, latest published record
**23112891**, and reconciled unsubmitted draft **23114217**. The published record
has ten inherited files and no API `metadata.version`; semantic v25 is the
prepared candidate, not a verified complete publication. Automatic GitHub
archiving stays **OFF for all eight public repositories**. The separate static
companion DOI family22922927/published23111008 is preserved. Preparing files locally
does not publish a record, reserve a new main family or create a GitHub release.

The science and local files are ready; there is **NO_CURRENT_PUBLISH_GO**.
The user's latest instruction narrows removal to updates that were not verified.
Main record23112891 matches all ten reviewed inherited files but lacks the twelve
intended complete-v25 additions. That is an incomplete update, not invalidity of
the inherited science. Preserve the companion and every other verified study
unless its corresponding reviewed manifest demonstrates a failure. Two earlier
owner removal requests for23112891 and23111008 each returned HTTP500, with no
accepted-request response; no removal success was verified. The effect on the
main draft remains under review. Re-read current authenticated state before any
metadata, upload or publication action; an offline snapshot cannot establish
that a target remains current.

Inspect `BUILD_RECEIPT.json`, `FILE_MANIFEST.json` and `METADATA.json`. The builder
checks the final package bytes, its manifest, report and external fresh receipt.
Its completed-science assertion is taken from the separately reviewed pinned
completion evidence; it runs no model or array calculation. All prior FAILs and
the first execution-schema failure remain disclosed.

Prepare a concrete folder of the **twelve additions** outside the repository:

```bash
python prepare_candidate.py \
  --repository-root /path/to/HDblast \
  --expected-manifest-sha256 FULL_SHA256_FROM_BUILD_RECEIPT \
  --output-dir /path/to/new-upload-folder
```

This verifies/reconstructs the earlier multipart packages and copies the active
attachments. If the previous nine files already exist in a verified preparation
folder, add `--prepared-previous-dir /path/to/previous-nine-file-folder` to verify
and reuse those complete bytes instead of reconstructing multipart ZIPs. Every
reused filename, byte count, SHA256 and final MD5 is still checked. It does not
copy the ten inherited files, which must stay unchanged
in the existing Zenodo draft. `METADATA_REQUEST.json` and
`PREPARATION_RECEIPT.json` are local sidecars, not extra Zenodo attachments.
Do not upload an entire candidate directory blindly.

The public GitHub candidate omits the 32,195,260-byte complete ZIP to keep each
distribution file below24MiB. Reconstruct that ZIP from the registered public
package's two pinned parts with:

```bash
python research/HDBLAST_ACTIVE_SOURCE_LEDGER_PACKAGE_20261002/reassemble.py \
  --parts-manifest research/HDBLAST_ACTIVE_SOURCE_LEDGER_PACKAGE_20261002/PACKAGE.json \
  --expected-sha256 0cd3dd2da604444541c06e98167f537f769c0ce430108bf02354e5090f328c4f \
  --output /path/to/HDBLAST_CHECKPOINT_20261002_ACTIVE_SOURCE_LEDGER.zip \
  --receipt /path/to/reassembly-receipt.json
```

Then copy the verified ZIP to this candidate's
`attachments/HDBLAST_CHECKPOINT_20261002_ACTIVE_SOURCE_LEDGER.zip` before running
the preparation command. The private external staging candidate already holds
the complete ZIP. All other candidate files are copied without alteration.

If API transport still fails, open the existing draft at
<https://zenodo.org/deposit/23114217>, preserve all ten inherited filenames and
upload exactly the twelve additions listed in the manifest. Apply the complete
`METADATA.json` values through the form. Never create another main-family record,
delete inherited files, or publish the inherited-only clone. The prior unchanged
metadata PUT returned HTTP500; retry only after a concrete transport/schema
diagnosis, or use this existing-draft browser fallback. Reuse the configured
credential, never add its value to a file or URL.

Before publication, save the authenticated draft GET response locally and run:

```bash
python verify_saved_record.py \
  --saved-response /path/to/draft-get.json \
  --expected-record-id 23114217 \
  --expected-manifest-sha256 FULL_SHA256_FROM_BUILD_RECEIPT \
  --output /path/to/draft-check.json
```

Require `PASS_COMPLETE_SAVED_DRAFT`: correct family/draft, nonempty complete
metadata, and all22 exact filenames, byte sizes and MD5 checksums. HTTP200 or a
reserved draft DOI alone is insufficient. After the authorized publication,
repeat the saved-response check with `--published`, using the actual returned
record ID; require submitted=true, state=done, the corresponding record DOI and
all22 complete pins before recording publication success. This candidate's
helpers perform no network, upload or publish action.
