# Deterministic memory-followup package and honest failure controls

Both package modules belong in the prospective checkpoint's `code/` directory
before registration. The builder never runs scientific code. It authenticates
the registration, independent manifest and post-freeze receipt, hashes the
complete curated payload including the receipt, and constructs the full ZIP
twice using fixed timestamps/permissions/order and bounded streaming. It
compares both complete byte streams and verifies every payload digest and CRC.

The explicit exclusion constants contain the payload manifest, logical ZIP,
ZIP SHA, part ledger and 64 fixed part names. Listing the entire fixed namespace
in advance avoids a part-count/hash cycle and keeps the unchanged replay's
literal exclusion verifier compatible. The maximum package is 4GiB. A ZIP
at most64MiB is tracked as one file; a larger ZIP remains outside the checkpoint
and is split into tracked64MiB parts, with a shorter final part. Every part
and the logical ZIP have SHA256 receipts in `PACKAGE_PARTS.json`. Final replay
or publication receipts belong outside the protected payload.

After results/evidence are curated, root can run:

```bash
python CHECKPOINT/code/build_package.py --checkpoint CHECKPOINT --external-dir FRESH_EXTERNAL_DIRECTORY
```

The standard-library `reassemble_package.py` is standalone. It verifies part
order, every part length and digest, the reconstructed full ZIP, exact payload
membership, all hashes/CRCs, manifest, receipt and frozen-input bindings. It
requires a fresh external destination. A caller can independently supply the
published full ZIP SHA:

```bash
python reassemble_package.py --parts-manifest CHECKPOINT/PACKAGE_PARTS.json --output EXTERNAL_FULL_ZIP --expected-sha256 PUBLISHED_ZIP_SHA256
```

The synthetic tests execute real single-file and split-file paths. A65MiB
incompressible fabricated payload exceeds the actual64MiB threshold; two
complete builds and a copied standalone helper produce identical full ZIP
bytes. The unchanged32-command driver's input verifier accepts the split
package. Corrupt parts, ordering, lengths, public bindings and digests fail.
No physical source, mode producer or quadrature is invoked.

The sole staged workflow replaces `hdblast_metric_response_memory_followup.yml`
when the final package and part ledger are published. It verifies and extracts
the stored third package, runs the preserved original quadrature-failure
control, then invokes the separately frozen memory-followup expected-scientific
failure control against the extracted standalone checkpoint. The latter keeps
the original32-command driver and all scientific gates unchanged. It must
complete both primary sources/all12 rows and all four full independent archives
within the original budgets, then reproduce the registered internal scientific
failure including the known fine positive_B Ward witness. An unsupported
failure remains fatal. Control success always retains the underlying science
status FAIL and metric calibration unestablished.

The workflow uses Python3.12, the five exact package pins, one numerical thread,
45minutes, read-only contents permission, sparse checkout and no persisted
credentials. Assembly and execution evidence is always uploaded for7days.
The second round's hardware-sensitive memory failure remains immutable stored
evidence; it is not rerun or advertised as a fresh second physical reproduction.
