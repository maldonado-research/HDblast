# Final publication evidence

This directory records verification of the [frozen screen magnetic curvature checkpoint](https://github.com/maldonado-research/HDblast/tree/3c49f3377f2e1f471012f905668ca4e7a96fee56/research/HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026), merged through [HDblast PR 88](https://github.com/maldonado-research/HDblast/pull/88), and the shared project page merged through [PR 79](https://github.com/maldonado-research/maldonado-research.github.io/pull/79). It is added on a separate verification branch, leaving the merged scientific source and earlier manuscripts frozen.

Read [FINAL_HANDOFF.md](FINAL_HANDOFF.md) for the result, assumptions and next scientific questions. [FINAL_PUBLICATION_VERIFICATION.json](FINAL_PUBLICATION_VERIFICATION.json) reconciles the exact merged commits, all required CI, reproduction parity, preservation and complete public readbacks.

From this directory, verify the retained files without writing or rerunning the science:

```sh
python -I -S -B verify_evidence.py
python -I -S -B assemble_verification.py --check-existing
```

The first checks the complete inventory, byte lengths and SHA256 hashes against [EVIDENCE_MANIFEST.json](EVIDENCE_MANIFEST.json). The manifest excludes its own bytes to avoid self-reference; the immutable Git commit supplies an external content identity. The second recomputes the saved publication receipt from the retained evidence. These commands verify the saved record; they do not independently repeat its remote observations or prove a physical origin model.

[EVIDENCE_SELECTION.json](EVIDENCE_SELECTION.json) identifies the selected proof inputs. The top-level receipts include all 69 immutable changed files, three complete live HTML pages, and parsed pre/post native job evidence with 28 scientific hashes and 32 run identities. `independent-review/` retains the corresponding independent checks. `independent-package-review/` and `packaging-review/` retain source/authentication/replay provenance. Some records describe earlier dated preparation; final publication status is in the final receipt and postmerge review.

`proposals/` records the four additive public-text insertions, original tree metadata and frozen manuscript identities. No other shared project was edited. The evidence update itself only adds files beneath this directory, preserving all 11,978 leaves of the merged source tree exactly.

`follow-on/` is separate from the 60-file scientific capsule. It contains the standard symmetry obstruction to selecting a local natural timelike observer from the symmetric background alone, a proposed additional clock/transport datum, and an exact cancellation for the explicitly restricted same-harmonic scalar observer tilts. The transport enclosures, arbitrary-tilt extension and general gauge-invariance question remain open. These notes add no actual curvature-sign or nonlinear-solution claim to the frozen result.

Raw job/connector logs, downloaded paper bodies, private editorial/account/session history, workspace clones, temporary fixtures and bytecode are excluded. Full CI artifacts were not downloaded. Internal AI-assisted review is distinct from external specialist review and proof-assistant formalization. No new Zenodo version or journal submission was made.
