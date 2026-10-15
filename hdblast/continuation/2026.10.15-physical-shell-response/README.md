# Physical shell scalar response: publication evidence

The [research capsule](../../../research/HDBLAST_PHYSICAL_SHELL_RESPONSE_20261015/README.md)
records a conditional qualitative nonzero physical shell scalar response for
one specified aligned quadrupole with both displacement and velocity nonzero.
Its proof retains the coupled metric–scalar reconstruction, exact initial-data
and finite-time evolution hypotheses, shell matching laws and controlled
positive-cone region. The intrinsic proper-time window observable has a
strictly positive leading quadratic coefficient for sufficiently small nonzero
amplitudes. That coefficient depends on the first physical tangent and not on
the second exact initial correction. This is a scalar response diagnostic,
not a gravitational energy or numerical waveform.

- [Source/ZIP pins and reviewed proof hashes](FROZEN_SOURCE_RECEIPT.json)
- [Merged PRs, reviewed trees and five successful postmerge workflows](PUBLICATION_RECEIPT.json)
- [Three passing checks on the reviewed public head; no configured shared-site checks](PRE_MERGE_CHECKS.json)
- [67 preservation checks and 9185 unchanged other old files](SITE_INTEGRATION_PREFLIGHT.json)
- [Four byte-pinned reference layout exceptions](SOURCE_LAYOUT_PREFLIGHT.json),
  authenticated by the preservation receipt and frozen scientific manifest.
- [39 existing workflows preserved](CI_INSTALL_PREFLIGHT.json)
- [Separate canonical-normalization clarification](CANONICAL_NORMALIZATION_CLARIFICATION.md),
  explaining the uniform momentum rescaling between retained conventions;
  [its pinned receipt](CANONICAL_NORMALIZATION_CLARIFICATION.json) authenticates
  the referenced source context. This clarification changes no frozen proof
  or source bytes and makes no new proof claim.
- [Authenticated fresh replay](fresh-replay/REPRODUCTION_RECEIPT.json): six positive
  processes, thirty-two external wrong-formula processes rejected, zero embedded
  controls, and forty-one finite positive check groups in each Python mode.
  These algebraic checks do not computer-certify the analytic proof.
- [182 byte-equal research, summary and live-page targets](FINAL_LIVE_READBACK.json)
- [Next validated-response calculation](NEXT_VALIDATED_RESPONSE_TARGET.md)

The 164-file source capsule, its manifest and matching ZIP,
the new workflow, nine current summaries, two initial delivery files and four
live pages were read back after successful deployments at the observed merge
commits. These receipts record dated observations; later intentional changes
to main or website content can change the resulting hashes.

With the source worktrees at their recorded reviewed heads and the merge
commit objects fetched locally, rerun the standard-library helper as follows:

```sh
python -I -B verify_published.py --hdblast-root /path/to/HDblast --shared-root /path/to/maldonado-research.github.io --merge-evidence MERGE_AND_DEPLOYMENT_PREREQUISITES.json --output /tmp/physical-shell-response-readback.json
```

The helper checks both observed merge/deployment prerequisites and local
reviewed-versus-merged tree equality before issuing bounded concurrent reads.
Its exact bytes are pinned by READBACK_HELPER_PREFLIGHT.json. The fresh replay
delivery contains 77 receipt/result/log files: 38 results/logs per execution
mode and one overall receipt. Each original record was authenticated before
copying; the temporary executable copies are excluded. FRESH_REPLAY_FILES.json
pins those delivered replay bytes. PUBLICATION_EVIDENCE_MANIFEST.json pins all
evidence files except itself. The final readback of the later evidence PR is
retained separately, avoiding a self-referential receipt.

No numerical waveform, coefficient value, effective amplitude bound or exact
nonlinear first-arrival time is established. Global stability, observational
agreement, external peer review, certified novelty and a higher-dimensional
origin of the Big Bang are not established. Earlier dated entries retain
their original scope. This checkpoint makes no new Zenodo version or APS
submission.
