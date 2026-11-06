# Reproduce the finite high-order calculation

Python 3.12 or later is sufficient. The replay uses the standard library,
exact rational arithmetic, finite polynomial coefficient envelopes and
outward Taylor majorants. It imports no old private producer or report
oracle. It does not run the historical PDE simulations or evolve the
nonlinear Einstein–scalar family.

Obtain the complete capsule at the reviewed commit and the external
MANIFEST.json SHA256 pin from
`.github/workflows/weighted-regularity-20261106.yml` at that same commit.
A pin recomputed from untrusted inputs is not an authenticity anchor.
Set PIN to that literal 64-character value. From the repository root,
choose fresh report paths outside this closed capsule:

```sh
python -I -S -B research/HDBLAST_WEIGHTED_REGULARITY_20261106/checks/check_checkpoint.py --manifest-sha256 "$PIN" --output /tmp/weighted-regularity-normal.json
python -O -I -S -B research/HDBLAST_WEIGHTED_REGULARITY_20261106/checks/check_checkpoint.py --manifest-sha256 "$PIN" --output /tmp/weighted-regularity-optimized.json
cmp /tmp/weighted-regularity-normal.json /tmp/weighted-regularity-optimized.json
python -I -S -B research/HDBLAST_WEIGHTED_REGULARITY_20261106/checks/run_controls.py --manifest-sha256 "$PIN" --output-dir /tmp/weighted-regularity-controls-fresh
```

The controls directory must not exist before the command. Canonical
reports should be byte-identical in both modes. The control receipt
records its passing baseline and six intended failures: incorrect pin,
source-byte change, unlisted member, internally rehashed false scalar
coefficient ceiling, symlink, and a lexical alternative script path
attempting to write inside the capsule. The mathematical mutation must
reach its false inequality, rather than fail an unrelated pin check.

Read the proof notes alongside the report. Finite replay establishes
source identity and arithmetic conditional on the continuous arguments;
it does not prove the retained background/inverse theorem or external
originality. The global bootstrap, closure, covariant reconstruction and
physical Cauchy conversion received separate internal review. None of
these checks establishes a uniform amplitude interval, an evolved shell
observable, a cosmological prediction or observational confirmation.
