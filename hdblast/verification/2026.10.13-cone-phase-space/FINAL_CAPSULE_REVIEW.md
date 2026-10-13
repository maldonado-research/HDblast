# Bounded final frozen-capsule review

13 October 2026. **PASS within the reviewed scope; no blocker found.**
This is internal static, provenance and replay-evidence review. The reviewer
did not rerun any scientific checker, the reproduction runner, typesetting,
remote publication or GitHub Actions. Independent file hashing and byte
comparisons were performed. Nothing inside the frozen capsule was changed.

The reviewed capsule is `portable/`, with manifest SHA-256
`6021915456815936a80c864d5e503b05725325a2cb15cbd882051fade0890194` and
runner SHA-256
`b5f80083ac29747723883af4fbd83aaacfc12dbc49f28724f5249e144881bcce`.

- The manifest covers exactly 63 other files and 1,084,540 bytes. Every
  digest and byte count matches, with no symbolic links or unmanifested
  files. All seven copied prior inputs match the separately pinned prior
  manifest.
- The final contracts are 79 identities and 23 wrong-formula controls per
  mode: covariant 20+4, domain 24+7, corrected tensor 23+6, independent 12+6.
  Their counts, script/source pins where recorded, and frozen receipts agree.
- All 19 authoring/public adaptation records match the actual files.
  Scientific checker content is unchanged by the documented path packaging;
  normalized authoring/public receipt contents agree. The nine literature
  mappings and public review references resolve to their recorded originals.
- The corrected tensor note and checker now include the full ambient
  near-null counterexample: `r=epsilon^3` and
  `xi_ell0^1=ell0(epsilon^3-epsilon^-1)`. The note distinguishes this
  divergence from bounded finite-time rays and does not infer an additional
  physical charged mode. Its final 23+6 receipt includes the three added
  identities and added boundedness rejection. This resolves the previously
  excluded correction for this bounded inspection.
- The manuscript Markdown, TeX and PDF exactly match the final mathematical,
  bibliography and build pins; the visual receipt pins the same PDF. Review
  mappings correctly distinguish the unchanged authoring manuscript from
  path-adapted supporting documents. The historical build-stage status is
  explicitly explained in `reviews/README.md`.

The completed `final-replay-001/REPRODUCTION_RECEIPT.json` reports
`PASS_REVIEWED_CONE_CURRENT_REPLAY` for Python 3.12.14, SymPy 1.14.0 and
mpmath 1.3.0. Inspection confirms that all eight normal/optimized output
JSONs are byte-identical to the frozen reviewed receipts and have the hashes
recorded in that replay receipt. The relocated capsule is byte-identical to
the frozen source. Both corrupted-input logs end with the intended source
integrity error; the altered-manifest log ends with the external-pin error;
both injection logs end with `injected_failure: residual 1`. Neither injection
created a success receipt. Thus the negative controls failed for the intended
reasons, rather than an unrelated dependency or invocation error.

The final CI candidate contains the matching external manifest and runner
pins with no placeholders. Its standard-library bootstrap authenticates
both before project execution; the authenticated runner verifies the entire
capsule before dependency installation. The pinned Python, narrow triggers,
fresh runner-temp replay destination, pipefail handling, and receipt/log
artifact paths remain consistent with the earlier static review. This is
review of the workflow source, not a claim that remote CI has run.

The result certifies consistency of this frozen package and the inspected
local replay evidence. It does not add global solution existence, physical
phase-space completeness, new literature priority, external peer review,
publication/readback verification, or a cosmological-origin result.
