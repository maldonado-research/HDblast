# All-angular linear shell-contact control and exact boundary reduction

GitHub research checkpoint. Artifact identifiers use 15 October 2026 UTC;
this work began 14 October Pacific time. The results are conditional on the
registered Einstein–scalar background, physical copy-even shell laws and
the fixed controlled positive-cone collar stated in the proofs.

The previous contact construction produced the complete second Taylor
coefficient for one aligned quadrupole in a fixed finite angular space.
This checkpoint closes the specific **all-angular frozen-linear control**
target and proves a separate **exact local nonlinear boundary reduction**.
It does not construct an exact nonlinear family through shell contact.

| Result | Established scope | Remaining limit |
|---|---|---|
| All-angular linear contact estimate | Positive auxiliary Sobolev control and unique smooth evolution for compatible smooth data, with constants independent of angular frequency and of the full tensorial SO(4) cutoff | The operator is linearized at the registered background; smooth boundary forcing has a declared lifting regularity cost |
| Exact nonlinear boundary reduction | A local pointwise state flow converts all twelve physical and normal-gauge normal rows to homogeneous Neumann conditions; the four coordinate Dirichlet rows and common quasilinear wave principal symbol are preserved | Variable-coefficient energy estimates, compatible iteration, parameter regularity and a lifespan through prescribed contact still require proof |

## Why the full boundary matrix can be controlled

All fifteen metric fields and the scalar share the scalar wave principal
symbol in relative harmonic gauge after their invertible equation
normalizations. The full twelve-row Robin matrix is generally nonsymmetric
in a positive component norm. A bounded invertible collar matrix exponential
turns its homogeneous boundary conditions into Neumann conditions and adds
only bounded lower-order bulk terms. The four Dirichlet components are kept
separate by a fixed positive bundle metric and split-preserving connection.
The resulting homogeneous boundary flux vanishes, giving a genuine spatial
energy estimate with no angular cutoff in its constants.

Higher estimates explicitly track time/angular commutation, normal derivative
recovery and sufficient forcing norms. Inhomogeneous boundary data enter
through a smooth explicit lift. This is not a lossless estimate for arbitrary
`L²` Neumann boundary forcing. The positive norm is analytic control, not
physical gravitational energy, a ghost test or a stability theorem. The
same method gives all-angular uniqueness for the subsidiary constraint
system, with its required zero initial gauge and physical constraint data.

The separate nonlinear state flow straightens the exact physical and normal
harmonic laws without dropping a nonlinear remainder. It is local near the
background state and fixes the zero perturbation. Its preservation of the
common quasilinear principal symbol prepares a mixed boundary existence
proof; it does not replace that proof.

## Read the arguments

- [Main all-angular linear theorem](ANGULAR_UNIFORM_LINEAR_CONTACT.md)
- [Exact local nonlinear boundary reduction and remaining analytic steps](NONLINEAR_BOUNDARY_STRAIGHTENING.md)
- [Earlier complete finite-angular second-order contact coefficient](../HDBLAST_FULL_SECOND_ORDER_CONTACT_20261014/README.md)
- [Earlier exact localized initial-data family and short precontact evolution](../HDBLAST_EXACT_LOCALIZED_PRECONTACT_20261014/README.md)
- [Next quasilinear contact target](../../hdblast/continuation/2026.10.15-angular-uniform-contact/NEXT_QUASILINEAR_CONTACT_TARGET.md)

The earlier exact initial-family theorem remains available at each fixed
`k ≥ 12`, with metric in `H^{k+2}` and momentum in `H^{k+1}` under its stated
assumptions. Its finite-regularity parameter-smooth family and short
precontact evolution are retained as inputs. The next task is to carry
that family through the full quasilinear mixed boundary construction with
exact initial harmonic completion, an explicit parameter derivative-loss
budget, a proved lifespan through contact and exact subsidiary uniqueness.
No common spatial `C∞` family or identification with one particular earlier
second initial correction is asserted.

## Primary-source context

The bounded literature assessment checks primary-source assumptions against
this model. It is not an exhaustive web search or an external novelty
certificate. Useful comparison sources include:

- Fournodavlos and Smulevici, [The initial boundary value problem in General Relativity: the umbilic case](https://arxiv.org/abs/2104.08851v1): a four-dimensional vacuum constant-umbilic boundary problem, requiring a separate derivation for the present five-dimensional scalar model and field-dependent tension.
- Kreiss, Reula, Sarbach and Winicour, [Well-posed initial-boundary value problem for the harmonic Einstein equations using energy estimates](https://arxiv.org/abs/0707.4188v2): harmonic Einstein boundary energy methods; its chosen boundary conditions are not automatically the present physical shell laws.
- An and Anderson, [Well-Posed Geometric Boundary data in General Relativity, III: Conformal-mean curvature boundary data](https://arxiv.org/abs/2503.12599v3): a geometric boundary formulation whose data and hypotheses must be distinguished from this shell model.

These comparisons do not transfer a quasilinear theorem automatically.
The actual shell rows, collar transformations and constraint cancellation
are derived for the registered model in this checkpoint. The distributed
literature assessment and internal reviews retain the detailed source
locations and qualifications:

- [Primary-literature assessment](literature/ANGULAR_UNIFORM_CONTACT_LITERATURE.md)
- [Astra final acceptance of the main note](reviews/ASTRA_MAIN_NOTE_ACCEPTANCE.md)
- [Independent Sol final review of both pinned files](reviews/SOL_FINAL_TWO_FILE_REVIEW.md)
- [Independent boundary referee acceptance of the main note and companion](reviews/BOUNDARY_MAIN_NOTE_ACCEPTANCE.md)

The main note received three final internal AI-assisted reviews. Astra
authored the nonlinear companion; its author review is distinguished from
the two other agents' independent reviews of that companion.

## Verification

Three independently derived standalone algebra checkers run under ordinary
and optimized (`python -O`) execution: **six successful positive processes**.
The fresh composite runner also exercises eleven actual injected corruptions
per mode, requiring **twenty-two rejections** with exit code one and the
specified failure reason. Optimized execution must preserve the scientific
outputs; only Astra's explicit optimization metadata flag is ignored.

| Checker | Reported checks in each positive execution | Scope |
|---|---|---|
| Astra | 20 exact identity checks | Actual boundary coefficients, collar conjugation, principal and flux algebra |
| Sol | 8 validation groups, with embedded wrong-formula controls | Exact differentiation of the boundary laws and quadratic sources, noncommuting matrix derivative and glancing normalization |
| Boundary referee | 5 groups containing 7 explicit checks | Direct-metric matrix, trace conventions, harmonic/glancing rank and Robin conjugation |

These counts use different units and are **not additive**. The finite checks
support the specified exact identities and controls; they do not
computer-certify the analytic PDE estimates or prove an exact nonlinear
spacetime. The written analytic arguments receive separate mathematical
review.

Four controlling Markdown documents from earlier checkpoints are included
as unchanged contextual input snapshots. [Snapshot provenance](INPUT_SNAPSHOT_PROVENANCE.json)
pins their origins and hashes. The source capsule, matching standalone ZIP
and [external frozen-source receipt](../../hdblast/continuation/2026.10.15-angular-uniform-contact/FROZEN_SOURCE_RECEIPT.json)
record the file inventory and authentication pins. Outputs from a new replay
must remain outside the frozen source tree.

## Reproduction

Use an external environment with **Python 3.12.14** and **SymPy 1.14.0**
(`requirements.txt`). Do not put a virtual environment or generated outputs
inside this capsule: the integrity verifier rejects extra files. The
standalone ZIP contains the same source bytes as the capsule.

Obtain the expected manifest SHA-256 and the verifier/runner SHA-256 values
from the reviewed external
[frozen-source receipt](../../hdblast/continuation/2026.10.15-angular-uniform-contact/FROZEN_SOURCE_RECEIPT.json),
checked at the published reviewed commit. Before executing project code,
authenticate `MANIFEST.json`, `verify_manifest.py` and `reproduce.py` with
independent standard-library hashing or a trusted hashing tool against
those external values. Do not infer the expected pin from the files being
verified. For example, `sha256sum MANIFEST.json verify_manifest.py reproduce.py`
computes the three values to compare with that receipt.

In the commands below, replace `EXTERNAL_PIN` with the complete reviewed
external manifest SHA-256. Choose a new or empty absolute output directory
outside the capsule. From this capsule's directory:

```sh
python -I -S -B verify_manifest.py --expected EXTERNAL_PIN
python -I -B reproduce.py --manifest-sha256 EXTERNAL_PIN --output-dir /tmp/hdblast-angular-uniform-fresh-replay
```

The first command authenticates every source-file size and SHA-256 and
rejects missing, extra or symlinked files. The composite runner repeats
manifest authentication before importing the project checkers, copies each
checker to a fresh directory and executes the six positive processes and
twenty-two actual corruption controls. Its output directory contains the
fresh JSON/logs and `REPRODUCTION_RECEIPT.json`. Successful reproduction
establishes the declared algebraic checks and pinned outputs, with the
separate analytic-review scope described above.

## Limits

This is a GitHub-only mathematical checkpoint. Internal AI-assisted review
is not external peer review. Exact nonlinear contact existence, a common
spatial `C∞` perturbation family, global stability, external novelty and a
higher-dimensional origin of the Big Bang remain unestablished. No new
Zenodo version or APS submission accompanies this checkpoint. All earlier
dated research capsules and their original scope remain preserved.
