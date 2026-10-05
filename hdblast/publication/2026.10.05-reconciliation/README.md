# Zenodo read-only reconciliation, 5 October 2026 UTC

The [status receipt](STATUS.json) records bounded read-only GET checks of the
existing draft and two published records. No metadata, file or publication
mutation occurred in this round.

- Main concept **10.5281/zenodo.17088132** resolves to published
  [22347452](https://zenodo.org/records/22347452), with ten inherited files.
- Existing linked draft [23114217](https://zenodo.org/uploads/23114217) is
  **unsubmitted**, remains in that concept family and holds ten inherited files.
  Its descriptive metadata is missing: the legacy readback includes only
  access and reserved-DOI fields, with no title, description, creator or version.
- The separate static companion [23111008](https://zenodo.org/records/23111008),
  concept **10.5281/zenodo.22922927**, remains published with its two original files.

During the 4 October UTC completion attempt, metadata saves returned HTTP 200
but actual saved metadata readback was empty. The attempted save removed the
draft's descriptive fields. Earlier HTTP 400 upload and HTTP 500 metadata
failures are preserved as diagnostics, rather than evidence of completion.
Published records and inherited files were preserved.

The next recovery step is to restore the reviewed fields through the account's
signed-in Zenodo editor and **Save draft**. Do not publish the incomplete draft.
The [prior reviewed browser fields and 23-file candidate](../2026.10.03-verified-source-operator/README.md)
remain immutable preparation evidence. After recovery, read back the saved
metadata and verify every intended filename, byte size and checksum before
publication. A reserved DOI or HTTP 200 response is not a publication receipt.

The [5 October mathematical checkpoint](../../../research/HDBLAST_CHECKPOINT_20261005_FINITE_BAND_WARD/RESULTS.md)
is published on GitHub only. It has not been uploaded to this draft, and does
not silently replace the prior 23-file candidate. Preserve the main DOI family,
reuse the existing draft and retain the owner's removed-record tombstone for
23112891. Automatic GitHub archiving remains off; no test release or duplicate
standalone record was created.
