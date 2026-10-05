# Standalone finite-band Ward mathematical checkpoint

The [research result and limits](../HDBLAST_CHECKPOINT_20261005_FINITE_BAND_WARD/RESULTS.md)
derive exact retarded density/pressure kernels, a complete state/residual/contact
conservation-error identity and conditional continuous-band source-truncation
bounds. They also demonstrate why normalization, conservation and nine probes
cannot certify the incoming continuum state. The full twelve-case stress/contact
certificate remains **UNRESOLVED** and earlier metric calibration remains **FAIL**.

The [ZIP](HDBLAST_CHECKPOINT_20261005_FINITE_BAND_WARD.zip) has **57 regular-file
members**, **173,211 bytes**, SHA256
`dbd69f253b7a6d444253496ee769b2b5942ead17eb432a95a4c06ddfde669812`.
Its manifest SHA256 is
`562a38390e2c995e75b5992493095ec151d9c227e6fa5e00d4d9608f743d9e38`.
[`ARCHIVE.json`](ARCHIVE.json) records these external pins.

The archive contains the proofs, three separate exact verifiers, sealed receipts,
independent mathematical review, preparation history, seven inherited premise
copies and a bounded six-paper literature review. It contains no tokens, private
chat transcripts, third-party PDFs or physical arrays.

With Python 3.12, SymPy 1.14.0 and mpmath 1.3.0, run:

```sh
python -B replay_archive.py --packet . \
  --output-directory /tmp/hdblast-finite-band-clean-replay
```

The replay checks archive size/hash, membership, regular-file safety, CRC and
the externally pinned payload manifest, then runs all verifiers in ordinary
and optimized Python. The 41-check primary result must reproduce byte for byte;
the 55-check/12-control and 42-check/10-control symbolic results must reproduce
every scientific field. Explicit runtime fields may differ. Frozen payload
bytes must remain unchanged.

[`FRESH_REPRODUCTION.json`](FRESH_REPRODUCTION.json) records the completed clean
extraction and six successful verifier executions. No physical source callback,
saved-array decode or trajectory rerun occurred. These checks are internal
computational verification, not external peer review or proof-assistant
formalization. The model's physical viability and external novelty remain
unestablished. This package is available on GitHub; it is not a new Zenodo version.
