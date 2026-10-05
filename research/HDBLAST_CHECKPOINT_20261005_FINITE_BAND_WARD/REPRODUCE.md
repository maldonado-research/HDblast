# Reproduce the finite-band mathematical checkpoint

Python 3.12 and SymPy 1.14.0 are the reference environment. The primary 41-check
verifier uses the standard library alone. Install the two pinned symbolic
dependencies for the other independent routes:

```sh
python -m pip install -r requirements-replay.txt
```

From this directory, run the sealed checkpoint replay with its manifest hash
from the standalone package's `ARCHIVE.json`. Use a fresh output directory
outside the checkpoint:

```sh
python -B execution/replay_checkpoint.py \
  --checkpoint . \
  --expected-manifest-sha256 MANIFEST_SHA256 \
  --output-directory /tmp/hdblast-finite-band-new-replay
```

The manifest hash is an external trust input, not supplied by the checkpoint
itself. The script validates the complete payload and rejects symlinks,
missing/extra files, size or checksum differences and an overlapping output
directory. It executes all three exact verifiers in normal and optimized
Python. Scientific receipts must reproduce the archived results; only the
explicit runtime version/optimization fields can differ. It validates the
whole payload again after execution and writes `REPLAY_RECEIPT.json` outside it.

The adjacent standalone package automates archive validation, extraction into
a fresh directory and replay:

```sh
python -B ../HDBLAST_FINITE_BAND_WARD_PACKAGE_20261005/replay_archive.py \
  --packet ../HDBLAST_FINITE_BAND_WARD_PACKAGE_20261005 \
  --output-directory /tmp/hdblast-finite-band-new-archive-replay
```

These are exact algebra and conditional bound checks. They import no project
numerical module and evaluate no physical source, mode or saved array. They
do not attempt the unresolved twelve-case stress/contact calculation.
