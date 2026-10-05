# Zenodo draft repair retry, 5 October 2026 UTC

The existing HDBLAST draft is **not ready to publish**. Three diagnosed metadata
restoration attempts returned HTTP 200, but fresh authenticated readbacks still
have empty descriptive metadata. The final observed draft revision is **12**.
The [diagnostics](DIAGNOSTICS.json), [final readback summary](FINAL_READBACK.json)
and [failed saved-metadata check](FINAL_METADATA_GUARD.json) record the result.
The responsible server or transport layer remains unresolved.

| Record | Verified state |
|---|---|
| Main concept [10.5281/zenodo.17088132](https://doi.org/10.5281/zenodo.17088132) | Latest published record remains [22347452](https://zenodo.org/records/22347452), with ten unchanged inherited files |
| Existing main-family draft [23114217](https://zenodo.org/uploads/23114217) | Owned, unsubmitted, ten inherited files, descriptive metadata missing |
| Separate static companion [23111008](https://zenodo.org/records/23111008) | Published; its original two-file inventory and metadata are preserved |
| Owner-removed 23112891 | Its archival tombstone remains; no further deletion is needed |

## What was checked

The existing `ZENODO_ACCESS_TOKEN` binding is saved for HTTPS `zenodo.org` and
authenticated ownership reads succeed. No replacement credential is needed.
The same reviewed metadata was sent through HTTP/1.1 chunked streaming, an
independent HTTPx JSON client and the installed Go-based native API client.
Conditional writes advanced revisions 9 to 12 while leaving metadata empty.
The compact HTTPx serialization has a different byte hash but the same JSON
fields. Official Zenodo/Invenio examples confirm the supported JSON structure;
TLS verification and the configured proxy were preserved.

All [thirteen local additions](LOCAL_FILES_READY.json) passed fresh filename,
byte-size, MD5 and SHA256 checks. They total **424,668,350 bytes**. The reviewed
complete candidate has **23 files / 440,222,994 bytes**, including the ten
inherited files already in the draft. No file upload, publication, deletion,
new record or test release was attempted during this retry. Automatic GitHub
archiving stays off.

## Complete the existing draft

1. Use a session that can control the owner's signed-in browser, or open the
   [existing draft](https://zenodo.org/uploads/23114217) manually. This cloud
   task can queue navigation but has no browser inspection or editing tools.
2. Restore the reviewed [browser fields](../2026.10.03-verified-source-operator/browser/METADATA_FIELDS.md),
   including the exact description, notes, creator, version, license,
   references and related identifiers, then **Save draft**. The
   [modern request](MODERN_METADATA_REQUEST.json) also preserves historical
   custom fields and the explicitly labeled v24 additional description.
   API JSON is a reference, not text to paste into a description field.
3. Require a fresh authenticated saved-metadata check to pass before uploading
   the [thirteen full files](../2026.10.03-verified-source-operator/browser/ZENODO_UPLOAD_FILES.md).
   Preserve the ten inherited files. Upload the full reconstructed ZIPs rather
   than their transport parts; the outer recovery bundle is not an attachment.
4. Require exact saved metadata and all 23 filenames, sizes and checksums before
   publishing under the existing authorization. Then verify submitted/done
   state and a direct public record readback in the same concept DOI family.

The already prepared complete recovery bundle remains available in the cloud
workspace to the continuing session. Public download and reconstruction links
for each addition are in the exact inventory above. Reuse this draft and token.
Another record or repeated unchanged API write would not solve the blocker.

The [support diagnostic](SUPPORT_DIAGNOSTIC.txt) is prepared but **has not been
sent**. It contains no authorization values. Platform/API transport repair is
the other recovery route; the current evidence cannot identify which layer
loses or ignores the submitted fields.

The [3 October candidate](../2026.10.03-verified-source-operator/README.md) remains
immutable. The [5 October mathematical checkpoint](../../../research/HDBLAST_CHECKPOINT_20261005_FINITE_BAND_WARD/RESULTS.md)
is available on GitHub and has not been silently added to this Zenodo candidate.
