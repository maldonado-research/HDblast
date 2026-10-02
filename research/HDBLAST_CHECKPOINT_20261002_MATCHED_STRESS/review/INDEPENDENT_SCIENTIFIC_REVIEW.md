# Internal scientific prose review

**Reviewer role:** author of the primary numerical route and post-run report
renderer. Despite this requested filename, this is an internal author review,
not independent implementation verification or external peer review. The
separate forced-mode route, tail-module review, and other internal reviews
supply distinct checks; this note does not substitute for them.

Reviewed `00_READ_FIRST.md`, `REPRODUCE.md`,
`review/CLAIMS_AND_NEXT_STEPS.md`, `outputs/NUMERICAL_RESULTS.md`, and the saved
summary/validator evidence. Input hashes are in the adjacent
`INPUT_SHA256.json`. This review only read existing files; it performed no
new source, quadrature, mode, stress, or coupled-physics calculation and
modified no checkout file.

No scientific claim requiring rejection was found in the reviewed prose.

- **Agreement and bounds are distinguished.** The density/pressure finite-K
  cross-route maxima, approximately 4.51e-17 and 6.24e-14, are numerical
  discrepancies, not total-error certificates. The larger K=256 stress UV
  bounds, approximately 2.51e-5 and 6.76e-4, concern only the omitted combined
  momentum band. Quadrature, mode refinement, arithmetic and Ward-ledger
  discretization retain their separate status. The report correctly includes
  derivative and local-contact tails rather than reusing the variance bound.
- **Control scope is accurate.** The twenty controls comprise five alternate
  operator/subtraction evaluations on saved modes, ten algebraic sensitivity
  diagnostics, and five rejected synthetic guards. They are not twenty new
  physical simulations. An improved-operator evaluation on unchanged modes
  is not a consistent changed-coupling theory. The detailed report also
  distinguishes nonzero sensitivity from operational gate rejection.
- **Physical interpretation remains bounded.** The reported density and
  pressure are first-order perturbations about the fixed incoming BD
  baseline. Their signs do not imply negative total energy, a thermal fluid,
  particle yield, or heating. The positive pulse has opposite reported
  pressure signs at the two post-pulse observations; the signed pulse does
  not share that pattern. No universal equation of state is inferred.
- **Mathematical dependencies remain explicit.** Both routes inherit the
  same minimal common-action prescription at H=1,r=2 and the same initial
  state. Their independent numerical agreement does not independently select
  that finite prescription or prove general off-reference, metric, bulk or
  shell dynamics. The Ward ledger is separately integrated from direct
  stresses; it does not define the tested density.
- **The next task is appropriately separate.** A prospective O(epsilon²)
  canonical excitation-energy calibration, compared across spectral,
  Bogoliubov and time-domain work routes with infrared/ultraviolet controls
  and the inherited reference/background drift, is a bounded follow-up.
  Its positive canonical occupation energy must not be substituted for the
  full second-order minimal stress or promoted to a heating/stability result.

Two presentation clarifications were sent to the root editor. First, describe
the opposite late pressure signs as a numerical result: the positive value
P(-1.5)=1.286269883e-4 is smaller than the conservative K=256 pressure-tail
bound, so that tail bound alone does not certify its sign. The existing
general statement that complete numerical errors are not certified is
correct. Second, the report's raw-hash table uses logical
`outputs/validation/CHECKS.json` labels, whereas the packaged check files are
`outputs/CHECKS.json`; the hashes agree, but matching the displayed paths
would simplify navigation.

`REPRODUCE.md` now describes eighteen commands, including the post-run
presentation steps being integrated. This prose review did not execute or
audit that evolving replay driver and makes no claim that a fresh or hosted
replay has passed. The document correctly requires actual execution receipts.
