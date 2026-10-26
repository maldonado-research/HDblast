# Final publication verification

This directory retains the publication evidence for the frozen [exterior tensor derivative checkpoint](https://github.com/maldonado-research/HDblast/tree/f7265bbf495ee0be9852ff5174e117c59d8b418a/research/HDBLAST_EXTERIOR_TENSOR_JETS_20261026), merged through [PR 87](https://github.com/maldonado-research/HDblast/pull/87). The shared project-page update was merged through [PR 78](https://github.com/maldonado-research/maldonado-research.github.io/pull/78).

Read [FINAL_HANDOFF.md](FINAL_HANDOFF.md) for the result, assumptions, limitations and next calculation. [FINAL_PUBLICATION_VERIFICATION.json](FINAL_PUBLICATION_VERIFICATION.json) records the exact merged commits, successful CI coverage, reproduction parity and byte readbacks. The scientific source capsule is preserved on main; this evidence is added separately on a verification branch.

From this directory, these commands check the saved evidence without writing files or rerunning the scientific calculations:

```sh
python -I -S -B verify_evidence.py
python -I -S -B assemble_verification.py --check-existing
```

The first command checks the complete file inventory against [EVIDENCE_MANIFEST.json](EVIDENCE_MANIFEST.json), including lengths and SHA256 hashes. The second recomputes the publication receipt from the retained local/API evidence. These checks authenticate the retained record; they do not independently reproduce its original remote observations or prove the physical model. Immutable commit links and Git blob identities support separate remote comparison.

`EVIDENCE_SELECTION.json` records the selected inputs and their hashes. `independent-review-evidence/` retains the final successful API snapshots, canonical source trees and source readback. Intermediate polling snapshots, raw connector job responses, raw logs, downloaded paper bodies, private editorial history and workspace clones are excluded. Pre-merge and post-merge native-log proofs retain parsed scientific hashes and expected rejection identities; full CI artifacts were not downloaded.

The four `proposals/` metadata files and frozen-entry list document the reviewed additive public copy and preserved earlier work. These are historical preparation records; final status is given by `FINAL_POSTMERGE_REVIEW.json` and `FINAL_PUBLICATION_VERIFICATION.json`.

`packaging-review/` retains wrapper verification provenance. The separate `follow-on/` directory contains a qualitative feasibility study and toy angular controls for a possible curvature diagnostic. The complete quadratic and observer terms are uncomputed, and no physical curvature sign or radiation observable is certified there.

No new Zenodo version or journal submission was made. Internal AI-assisted review differs from external specialist review and proof-assistant formalization. This checkpoint provides conditional mathematical bounds and no evidence establishing a higher-dimensional origin of the Big Bang.
