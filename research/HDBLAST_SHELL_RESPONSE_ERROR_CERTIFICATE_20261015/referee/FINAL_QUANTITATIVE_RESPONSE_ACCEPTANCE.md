# Final independent referee acceptance: conditional quantitative shell response

15 October 2026 UTC. Internal AI-assisted mathematical review, not external
peer review, formal proof certification, or a scientific novelty assessment.
Exact reviewed bytes are pinned in the accompanying JSON receipt.

**Decision: accept the conditional analytic result and its stated arithmetic
enclosure. No substantive blocker remains in the reviewed files.**

The primary note establishes, for the newly specified normalized aligned
quadrupole with `ell=1e-10` and observation coordinate time `3ell`,

\[
10^{15}<\mathcal C_{3\ell}<4\,10^{17},\qquad
\mathcal R_\varepsilon(\rho_b3\ell)
 =\varepsilon^2\mathcal C_{3\ell}+o(\varepsilon^2),
\]

conditional on the retained actual-background, localized exact-initial-data,
gauge, uniqueness, finite-time evolution, regularity, and causal-margin
hypotheses. The coefficient depends on the stated master normalization and
model units. The proof does not compute a waveform or give an effective
nonlinear parameter interval.

## Analytic checks

- The physical observable is the derivative trace
  `D=w_z-(C_b/2)w`, including the retained embedding/gauge identification.
  The proper-time coefficient is exactly
  `3*pi^2/(rho_b^2*p_b^2) * ||D||_L2^2`: the `cosh^3` weight cancels the
  squared time conjugation. No physical clock factor was dropped.
- A positive shifted energy cancels the physical shell kinetic flux. On
  `0<=z<=L-t`, the moving edge contributes
  `-[(e_t-e_z)^2+lambda e^2]/2`, so no artificial boundary condition is
  used. The growth constant and forcing factor in the primary certificate
  are conservative upper bounds for the independently derived energy law.
- Undifferentiated bulk energy alone does not control this derivative
  observation. The separate incoming characteristic identity and
  Minkowski estimate supply its `L2` normal trace. Boundary acceleration
  of the exact error is never assumed bounded by the bulk energy.
- The smooth odd/free Dirichlet comparator has the exact master initial
  data, zero shell value, trace `D_hat=q'+exp(-t)q`, bulk source `F=-W u`,
  and boundary source **`G=+D_hat`** in the error boundary convention.
  Its nonzero boundary residual is included, rather than imposed away.
- `|q'|<=16/ell` follows from `u^2 exp(4-u)<=16` for `u>=4`.
  On the rising half, `D_hat>=q'>=0`, giving the primary conservative
  squared-norm lower bound `2/ell`; a separate integration-by-parts check
  proves the stronger `4/ell`. No lower bound from a contradictory
  zero-trace extension is transferred to the actual solution.
- The retained exact background branch is the exact-decimal near-`phi=+1`
  branch at `delta=1/1000`. Its analytic endpoint estimates are applied to
  that branch. The inward ODE signs, positive shell inertia, nonzero `p`,
  first-exit continuation argument and derivative-free potential identity
  were checked. Historical floating profiles and the different S85 branch
  do not enter the enclosure.
- The energy-domain initial length `L=8ell` is distinguished from the full
  nonlinear correction annulus. The displayed sufficient geometry gives
  `L_corr+(5/4)T_evol=59ell/8<8ell`, with `T_1<T_evol`. The claim that
  nearby radial speeds are below `5/4` follows qualitatively from the
  retained parameter continuity after shrinking the parameter interval;
  no effective value of that interval is inferred.
- Coefficient uncertainty, initial derivative error, boundary residuals
  and the additional output-coefficient error are present in the general
  theorem. The present comparator's zero shell value makes its particular
  output-coefficient error vanish for every enclosed actual background.
- Sol's positive characteristic majorant, boundary-defect signs, explicit
  short-slab contraction and finite stepping were independently checked.
  Astra's energy/trace route, corner conditions, comparator algebra and
  optional small-depth criterion were independently checked. The two
  approaches avoid derivative loss through different explicit estimates.

## Arithmetic checks and independent corroboration

The primary `certify_response.py` uses rational input, integer-square-root
enclosures, a proved rational exponential majorant and Machin/alternating
arctangent bounds. Its source was reviewed. A read-only referee replay of
the pinned checker passed all 25 groups in both ordinary and optimized
isolated Python, with identical parsed results.

The independently written `independent_coarse_response_audit.py` imports
none of the primary implementation. Its 23 rational checks pass in both
modes and prove the same deliberately broad interval by cruder constants:
bulk forcing below 31, energy norm below 1532, derivative-trace error below
8, and elementary bounds `3<pi<22/7`. The separate receipt pins the exact
actual-background collar JSON used as its analytic input. The primary and
independent programs verify arithmetic consequences of mathematical input
inequalities; they do not machine-certify their PDE or background proofs.

## Required scope in publication

The result concerns one newly specified, very narrow smooth quadrupole in
fixed model units. It is not a numerical response for an earlier unspecified
width, a measured physical length/time, gravitational energy, or a cosmological
prediction. Large coefficients reflect the fixed bump normalization and do
not supply a finite nonlinear amplitude prediction. A resolved waveform,
practical-scale conditioning, an effective nonlinear Taylor remainder and
an observational connection remain open. No global stability, universal
observability, external peer review, certified novelty, or higher-dimensional
origin of the Big Bang follows from this checkpoint.
