# HDBLAST continuation: registered-shell controls and compatible initial data

Research checkpoint, 22 September 2026, America/Los_Angeles. This is a new local research package; it has not been published or externally peer reviewed.

We made a concrete advance at the original registered detuning, δ=0.001: a replacement solver recovers the independently established scalar instability rate under grid refinement, and short evolutions reproduce that growth. We also derived an exact constraint-transport diagnostic, identified an initial-data obstruction, and constructed floating-point candidates designed to remove it. These results strengthen the calculation needed to test the hypothesis. They do not establish that a higher-dimensional blast created our universe.

## What was actually completed

| Result | Evidence | What it establishes |
|---|---|---|
| Registered instability recovered | Six completed spectral controls, including four registered resolutions and a larger-domain comparison | Numerical calibration of the expected physical mode |
| Short time evolutions | Six completed nonlinear-PDE runs, including opposite signs and a dissipation comparison | Recovery of linear growth in the small-amplitude regime |
| Exact constraint transport | Independent symbolic identities and deliberate wrong-formula controls | A diagnostic separating constraint propagation from physical instability |
| Initial-data obstruction and repair candidates | Independent algebra; one-bump and balanced two-bump boundary-value calculations | Why satisfying initial constraints alone is insufficient, and how to remove the identified boundary obstruction to numerical tolerance |
| Proper-clock transformation | 54 exact rational-polynomial checks | A local coordinate identity preserving shell observables and junction conditions |
| Source and literature audit | Relevant D-Blast 3 folders 146–151, scoped project memory, original archives, public publication metadata, primary literature | A documented starting point and specific next physics tests |

The main measured rates are:

| Calculation | Growth per initial Hubble time |
|---|---:|
| Previous independent registered result | 1.657193631245 |
| New finest selected eigenvalue | 1.657193663767 |
| New fine-grid time evolution, fitted over t=1 to 2.5 | 1.657194158956 |

The finest spectral difference from the prior rate is about 3.25×10⁻⁸. The prior rate was known when selecting the eigensolver target, so this is a validation test, not a blind prediction. It is not an interval-certified numerical error. See [full control tables](CONTROL_ANALYSIS.md), [machine-readable analysis](CONTROL_ANALYSIS.json), and [independent solver review](numerics_review/NEW_SOLVER_REVIEW.md).

## The unchanged model being tested

We use the registered, matter-free five-dimensional Einstein–scalar model:

\[
W(\phi)=1-\phi+\phi^3/3,\quad
U(\phi)=\tfrac12(\phi^2-1)^2-\tfrac23W^2,\quad
\sigma(\phi)=2W(\phi)+\delta(1+c\phi),
\]

with δ=0.001 and c=2/1.0357712571566784−4/3. The metric is
\(ds^2=e^{2B}(-dt^2+dz^2)+e^{2A}d\mathbf{x}_3^2\), with shell z=0 and bulk z<0. Here t is coordinate time; δ is detuning. The static background has A=t+lnρ(z), B=lnρ(z). No physical length or energy scale has been selected.

The older stable shell with a quadratic tension correction is a different model. Likewise, adding matter, higher-curvature terms, extra fields or a timelike extra dimension would create an explicitly new branch of the investigation.

## A useful exact result: where constraint errors travel

Let C_H and C_M denote the Hamiltonian and momentum constraint residuals, with the definitions in the independent derivation. For the continuum equations without numerical dissipation, we verified

\[
(\partial_t-\partial_z)\{e^{3A}(C_H+2C_M)\}=0,\qquad
(\partial_t+\partial_z)\{e^{3A}(C_H-2C_M)\}=0.
\]

Consequently a small unweighted outgoing constraint error can increase as it moves into a region with a smaller warp factor, even though its weighted characteristic quantity is conserved. Evaluation on the archived static profile gives a gain of approximately 10,939 at t=2.5 along the ray leaving the shell. This offers a diagnostic for large apparent errors; it does not prove the cause of every failed earlier run. Numerical dissipation adds explicit source terms, and the discrete scheme does not inherit exact continuum conservation automatically.

This is a useful derivation for this project, based on standard constraint propagation. External mathematical novelty has not been established. Definitions, derivation, negative controls, and the dissipation-source correction are in [the independent audit](numerics_review/INDEPENDENT_NUMERICS_AUDIT.md).

## A stronger way to prepare a disturbance

For an initial slice with A=B=lnρ, A_t=1 and B_t=φ_t=0, using dz=dy/ρ, the Hamiltonian constraint becomes

\[
\rho_{yy}=\frac{1-\rho_y^2}{\rho}-\frac{\rho\phi_y^2}{6}-\frac{\rho U}{3}.
\]

We parameterize the initial scalar profile using
\(\phi_{yy}+4(\rho_y/\rho)\phi_y-U_\phi=\epsilon b\), with a smooth compactly supported b away from the shell and cone. This equation specifies initial data only; no external forcing is added to the physical evolution.

Define

\[
D=1-\rho_y^2+\rho^2\left(\phi_y^2/12-U/6\right).
\]

Direct calculation gives

\[
(\rho^2D)_y=\epsilon\rho^4\phi_y b/6.
\]

One localized disturbance generally leaves D_b≠0 at the shell. Even when the initial constraints and initial junction conditions hold, the second time derivatives of the junction conditions then fail: the B residual is −12ρ_yD_b and the scalar residual is −4ρφ_yD_b. This is a concrete obstruction, not merely a numerical tuning problem.

Two separated bumps with opposite weights can instead be tuned to cancel the integrated defect while solving both unchanged junction conditions. Floating candidates have been computed for positive and negative disturbance amplitudes. If b vanishes in an open neighborhood of the shell and D_b=0 exactly, the defect equation implies D=0 there; the data then match a local static solution and satisfy its formal local time-compatibility conditions. In floating arithmetic, these quantities are not certified zero: they are generally nonzero, and a rounded zero is not an exact certificate. These are candidates with measured numerical tolerances, not an exact existence certificate or an evolved nonlinear solution. See the [construction and obstruction report](initial_data/INITIAL_DATA_AND_CORNER_OBSTRUCTION.md), [numerical results](initial_data/BALANCED_CONSTRAINT_SEED_RESULTS.json), and [independent compatibility review](numerics_review/SEED_COMPATIBILITY_REVIEW.md).

The algebraic Hamiltonian residual partly checks the equation used to construct the data. It must not be mistaken for an independent finite-difference validation. The separate defect ledger and tolerance comparisons provide additional checks. A further independent differentiation of the saved profiles passed 15 checks, with maximum sampled normalized Hamiltonian residual about 2.103×10⁻¹⁰. That test uses a proper-distance sample grid; a future PDE run must still test its actual interpolated initial data on the evolution grid.

## What the literature contributes

The closest older numerical precedent already studies unstable de Sitter braneworlds evolving toward lower curvature or brane collisions. Its topology has two branes and differs from this model. That overlap limits any claim that the broad two-outcome idea is new. [Martin et al., BraneCode](https://arxiv.org/abs/hep-th/0309001).

A useful correction to the prior interpretation is that scalar oscillations are not necessary for all reheating mechanisms. Instant preheating can produce particles during a rapid passage through a mass threshold. Our current calculation contains no such matter sector. A meaningful extension must specify its action, physical field normalization, modified junction conditions, energy transfer, backreaction, and thermalization. [Felder, Kofman and Linde](https://arxiv.org/abs/hep-ph/9812289).

Recent screened work includes inflation from regular higher-curvature bulk black holes, warm braneworld inflation, a timelike-extra-dimension bounce proposal, and a gapped bulk graviton continuum. Each makes assumptions absent from the registered model. Their results are research leads, not interchangeable fixes. Dates, primary links, applicability and exclusions are in [the literature review](literature/PRIMARY_LITERATURE_AND_PUBLICATION_REVIEW.md).

The local proper-clock construction also passed exact checks. It keeps the shell fixed and makes its new time coordinate equal proper time wherever the coordinate Jacobian remains positive. It has not yet been implemented as a stable evolution method, and cannot recover spacetime absent from failed late-time data. See [the proof](literature/PROPER_CLOCK_PROOF.md).

## What remains unresolved

1. The full nonlinear fate at δ=0.001. Current new runs remain in the small-displacement regime; the initial-data candidates have not been evolved. Earlier Chat 13 nonlinear results at δ=.1 and .03 cannot be silently transferred to .001.
2. Constraint-converged finite-amplitude initial data after mapping onto the evolution grid. The older shifted-tension seed fails constraint refinement; a linear eigenmode only satisfies the nonlinear junction to first order.
3. Resolution of the outgoing moving wall and continuation through conformal-lapse failure. A grid refined only near the initial shell is insufficient without demonstrating wall resolution throughout.
4. Actual visible matter production, thermalization, a sustained radiation era, acceptable residual vacuum energy, and controlled bulk losses.
5. Cosmological perturbations and observations that discriminate this model from alternatives. Homogeneous evolution alone supplies neither the measured primordial spectrum nor evidence for our universe's origin.

The next computational experiment should put the balanced candidates onto a grid with independently checked initial constraints, then compare amplitudes, resolutions and domain sizes while following the moving wall. The proper-clock gauge needs its own benchmark before use for late-time conclusions. The next physics extension should add one explicit matter-production channel and derive the coupled equations from its action. Mixing several new theories at once would obscure which assumption produces any apparent success.

## Files, provenance and limits

Start with this file, then CONTROL_ANALYSIS.md. The folders `solver/`, `runs/`, `initial_data/`, `numerics_review/`, `literature/`, `state/`, and `inputs/` retain code, raw data, reviews and reference snapshots. [Reproduction instructions](REPRODUCE.md) distinguish quick replay of checks from regeneration of the scientific runs.

The source audit covered selected relevant folders and documents, not every file on the computer or the entire web. All original research and synced project sources were treated as read-only. The source v2 solver was changing externally during inspection; frozen snapshots and observed hashes distinguish versions. The supplied Chat 13 ZIP passes CRC but omits five intermediate files listed by its own manifest; those files were found intact in the extracted source folder. This limitation is preserved in the audit.

The supplied [Zenodo record 22347452](https://zenodo.org/records/22347452) is v24, dated 5 September 2026, concerning a conditional linear response certificate. It does not publish or validate these later shell-evolution results. This package makes no claim of a verified Big Bang mechanism, a Nobel-worthy discovery, or established external mathematical novelty. It provides specific, reproducible progress toward testing the hypothesis.
