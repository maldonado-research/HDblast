# Combined operator-audit status

13 October 2026. Read both reports together:

1. [The physical operator audit](PHYSICAL_OPERATOR_AUDIT_20261013.md)
   independently derives the generic scalar operator, full displaced-shell
   conditions, exceptional even harmonic classification, and local even
   kinetic coefficient.
2. [The odd kinetic completion](ODD_KINETIC_COMPLETION_20261013.md)
   **supersedes the unresolved local odd-action gap** described in Sections
   5 and 7 of the first-stage audit. It includes the moving integration
   domain, graph GHY term, Einstein tangential-divergence flux, and shell
   volume. Their cancellation establishes the full local odd kinetic
   coefficient.

At separated spectral roots, subject to the stated admissible gauge and
cone-flux assumptions,

\[
K_{\rm even}=18\Lambda\mathcal N_{\rm even},\qquad
K_{\rm odd}=18\Lambda\int w_{SL}|P_o|^2dy.
\]

The physical kinetic norm is not identified with the abstract
Sturm–Liouville norm. Global admissible-gauge coverage, the physical cone
symplectic domain, and completeness of the physical state space remain
outside these local action calculations. Exceptional harmonic subspaces
require their separate unreduced treatment. This package makes no
unrestricted no-ghost or nonlinear-stability claim.

## Reproduction

With Python 3 and SymPy 1.14.0, run from this directory:

```sh
python check_operator.py --output results/operator-check.json
python check_odd_kinetic.py --output results/odd-kinetic-check.json
```

The scripts also work from another working directory when invoked by path.
Isolated Python invocation (`python -I -B`, with or without `-O`) is supported
and was used for every portability replay.
Only Python, SymPy, this directory, and the bundled sibling `basis` inputs
are required. They import no historical derivation code and access no
network. The source inputs are pinned in [SOURCE_INPUTS.json](SOURCE_INPUTS.json).
Missing or modified pinned sources cause a failing exit status. There are
no assertion-dependent success conditions.

The operator check has **32 checks including three wrong-formula controls**;
the odd-action check has **11 checks including two wrong-formula controls**.
Source integrity is reported separately and is not added to these counts.
Receipts contain no elapsed runtime, generated timestamp, or absolute
location-dependent paths. Normal, optimized (`python -O`), and relocated
executions must produce byte-identical receipts.

The [portability review receipt](PORTABILITY_REVIEW.json) records those
verification results. [Original provenance](provenance/ORIGINAL_ARTIFACTS.json)
pins unmodified scratch sources, reports, and original receipts under
`provenance/original`. Those historical files may retain machine paths and
timing fields; they are provenance only and are not executed by these
portable checks. Historical September derivations under
`../basis/historical_operator` are exact source copies with their own pinned
provenance. No historical spectrum grid was copied or rerun.
