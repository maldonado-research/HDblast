# Conditional spectra of regular de Sitter shells

Ricardo Maldonado · Independent Researcher · 13 October 2026

This reproducible research checkpoint connects the previously certified static
scalar-boundary response to scalar spectral positivity and a tensor mass
threshold on the same regular Einstein–scalar shell family. The specified
operator domains, shell action and regular-cone assumptions are part of the
result. This is not evidence that a higher-dimensional blast caused the Big Bang.

The older project hilltop figure describes the original unstable shell branch.
This checkpoint concerns the distinct regular branch near `phi = +1` in the
same Einstein–scalar shell theory, using `eta = phi - 1`. It neither removes
the original branch's instability nor proves that evolution reaches this branch.

## Results and their scope

For every fixed detuning `0 < delta <= 0.002` in the registered exact-decimal
model, the even and odd scalar radial operators have their entire shifted
spectra bounded below by explicit positive background-dependent functionals.
The scalar essential spectrum begins at `mu^2 = 9/4`. The bound is not evaluated
numerically here, and no uniform lower bound as detuning tends to zero is proved.

The full local Einstein–Hilbert, Gibbons–Hawking–York and shell action gives
positive kinetic coefficients for the generic modes at spectral roots. The odd
calculation includes the moving integration domain and graph boundary. Constant
and translation-type harmonics are treated separately; those exceptions do not
remove all harmonics at their corresponding mass values.

In the transverse-traceless tensor sector, the only radial state below
`mu_T^2 = 9/4` is the even massless profile. A partner-operator proof covers the
nonzero spectra of both reflection parities and uses an exact negative bound
on the potential over the previously certified branch tube.

These are conditional linear spectral results. Global admissibility of the
gauge quotient, physical state-space completeness and possible cone symplectic
flux require separate work. Vector sectors, nonlinear stability, quantum decay,
cosmological evolution and observational evidence are not established. The
physical mass scale is set by the shell curvature, which tends to zero with
detuning; there is no uniform positive physical mass floor.

## Reading map

- [Ten-page technical note](manuscript/JACOBI_RESPONSE_TO_SPECTRAL_POSITIVITY.pdf).
- [Results and claim boundaries](RESULTS.md).
- [Scalar closed domains and full spectral gap](domain/SCALAR_SELF_ADJOINT_DOMAIN_AND_GAP.md).
- [Covariant operator and exceptional harmonics](operator/PHYSICAL_OPERATOR_AUDIT_20261013.md).
- [Completed moving-shell odd kinetic calculation](operator/ODD_KINETIC_COMPLETION_20261013.md).
- [Tensor partner proof](tensor/TENSOR_PARTNER_GAP.md).
- [Independent integration review](independent/INTEGRATION_REVIEW.md).
- [Closest prior work and literature search](literature/STABILITY_PRIOR_ART_REVIEW_20261013.md).
- [Unchanged static proof basis](basis/static/README.md) and [dated prospective bridge](basis/bridge/PROSPECTIVE_ZERO_MODE_STRATEGY.md).

The bridge and the first operator audit preserve their original status statements.
The later domain and odd-action reports explicitly resolve their identified local
gaps. Historical authoring scripts and receipts are provenance, not the portable
entry points. Third-party paper PDFs are not distributed in this capsule.

## Reproduction

Use Python 3.12.14 and the pinned dependencies:

```sh
python -m pip install -r requirements.txt
python -I -B reproduce.py --output /tmp/hdblast-regular-shell-replay
```

Choose a new output directory outside this capsule. The runner verifies the
complete manifest, executes the new checkers normally and under `python -O`,
compares their deterministic receipts with the archived results, replays the
unchanged static proof arithmetic, and checks that frozen source files remain
unchanged. Deliberately wrong formulas and corrupted-input controls must fail.
The new checkers cover 103 identities or exact signs and reject 30 deliberately
wrong formulas per mode. Static dependency and source-integrity controls are
reported separately. Symbolic checks validate specified algebra, not analytic
proof completeness.

The relevant prior literature already contains soft-boundary Sturm–Liouville
problems, modified inner products, ground-state identities and tensor
factorization. The present contribution is the explicit connection to this
certified family and its fully stated radial realizations; global mathematical
or phenomenological novelty is not claimed.

Substantial OpenAI Codex assistance was used for derivation, coding, literature
review and drafting, with separate AI-assisted derivations and checks. This is
internal verification, not external peer review. No funding or known competing
interests are declared. The previously approved APS submission package is
unchanged; this checkpoint is not a new journal submission.
