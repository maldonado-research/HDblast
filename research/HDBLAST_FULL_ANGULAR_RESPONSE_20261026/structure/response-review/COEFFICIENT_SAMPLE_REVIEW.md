# Independent review of the full23 coefficient and error calculation

26 October 2026. No blocking mathematical or implementation issue was found
in the reviewed native full23 response. This is an internal bounded review of
the new code and its reported artifacts. All seven frozen dependency roots
were authenticated before imports and again afterward; no frozen producer was
run. The native response source and outputs were unchanged during the audit.

The accepted conditional bound is `E=1.1787044290072725e28` in display precision,
with its exact rational value retained in the reports. It certifies the
expression-defined continuous candidate on the whole stated causal cone.
The six numerical locations are evaluations of that candidate, not the source
of the global error bound.

## Time-dependent coefficients and interval bounds

The code uses the genuine time connection `h=tanh(t)/rho`, with `t=ell*x`.
It does not confuse this h with the initial-data abbreviation U'(eta)/p.
For `0≤tlo≤t≤thi`, its bounds

    tanh(t) in [tlo−tlo³/3, thi],
    sech(t) in [1−thi²/2, 1]

are valid. Indeed tanh(t)≤t, so `tanh'(t)=1−tanh²(t)≥1−t²`; integrating gives
the lower cubic bound, and monotonicity handles the interval lower endpoint.
For `f(t)=sech(t)−1+t²/2`, put s=sech(t) in (0,1]. Then
`f''(t)=1+s−2s³=(1−s)(2s²+2s+1)≥0`, while f(0)=f'(0)=0.
This proves the sech lower bound. The upper bound is elementary. The audited
tiny-time box also makes the displayed lower endpoint positive, as explicitly
checked by the native code. Outward interval rounding retains these bounds.

The frozen background equations give `a_y=−rho^−2−p²/3` and
`p_y+4ap=V'`. Hence `b=−(a_y+a²)=d+p²/3`, with d=rho^−2−a².
The scalar potential is the complete registered U(eta), with no factor 1/2.
The shell substitutions are `sigma'=-2p` and `sigma''=4(1+eta)`.

The independent checker extracts only the five coefficient assignments from
the new `assemble` function. In a separate exact Laurent algebra it compares
all 4*23*23 entries of A,B,D,M to the authenticated frozen full operator after
these explicit substitutions. Every entry agrees identically for independent
formal background parameters. This is stronger than agreement at sampled
coefficient values and covers Delta, Ca, Ch, Ephi and every retained
metric/scalar coupling. It also verifies the selected rough-angular spectral
cap 24 used in the separate intrinsic-angular-derivative lemma.

The new response constructor was then used solely to reconstruct the current
interval coefficient boxes and constant matrices. These agree exactly with
the saved native report. Its actual background enclosure uses the retained
shell Taylor coefficients and validated fifth-derivative remainder, whose
premises remain inherited. The audit does not independently re-prove that
background ODE certificate. All reused mean constants are compared to the
same actual boxes and charged through the same entrywise errors; no reuse is
treated as an exact actual coefficient identity.

## Matrix norms and free-solution bound

Let H=diag(g_i) be the positive angular Gram matrix. The induced operator norm
of M equals the Euclidean operator norm of `H^(1/2) M H^(-1/2)`. For each exact
invariant block it is bounded by

    sqrt(sum_(i,j in block) (g_i/g_j) M_ij²).

The five blocks are orthogonal, so the full operator norm is the maximum of
their norms. The same bound applies to coefficient errors using nonnegative
entrywise magnitude caps. The review verifies exact off-block zeros, all 40
reported block Frobenius upper bounds, every global maximum join, and coverage
of every coefficient-box difference by its entrywise cap.

If the nominal initial value, derivative and velocity caps are nu,ndu,nv,
and the velocity interpolation cap is pv, the free candidate has

    X0=max((nu+(3/2)(nv+pv))/3, ndu+nv+pv).

This follows directly from reflected d'Alembert endpoints and the velocity
integral of length at most 2x. Reflection is an H-isometry. It does not assume
that a previous mean-only full-response cap also bounds all 23 amplitudes.

With actual-minus-nominal caps eu,edu,ev, the code correctly takes
`efree=max((eu+(3/2)ev)/3,edu+ev)` and `e0=efree+pv+S`, where S is the full
retained bulk/shell source-omission bound. The velocity interpolation charge
pv suffices in X because the derivative bound dominates its value bound.

## Feedback defect and complete error

The same continuous P1 velocity is already part of B0 and B1. Only the exact
nominal displacement and its derivative are replaced independently inside
the feedback forcing. Consequently the profile-error terms are exactly

    df_profile=(nA+nB)*eps_du+nD*eps_u,
    dg_profile=nM*eps_u.

The coefficient-error terms use the zero-wedge bounds on B0:

    df_coefficient=(eA+eB+(3/4)eD)*X0,
    dg_coefficient=(1/4)eM*X0.

These give `e1=(3/2)df+dg`. The review recomputes every one of these terms
directly from the initial report and matrix-norm report, as well as
`Ucap=(X0+e0)/(1−q)` and `E=q²Ucap+(1+q)e0+e1`. All exact rational values
match. The retained full-H q and S, C1 support/regularity, finite Robin
mismatch and selected pointwise evaluation constant are separately reviewed
in `SCOPE_REVIEW.md`.

The TT slots of this coarse full23 candidate are identically zero, because
its initial data have zero TT coefficients and its autonomous TT blocks have
no coupling to the other 21 coordinates. The omitted actual TT forcing is
included in S. This is consistent with a separately resolved much smaller
nonzero actual TT signal; the global full23 error bound does not resolve that
signal or assert its absence.

## Sample evaluation checks

The review independently integrates each P1 velocity segment, rather than
using the native hinge primitives, and applies the prescribed odd/even
reflection. It joins those exact primitives to the authenticated nominal
displacement and derivative endpoint intervals. All 414 value/x-derivative/
r-derivative component checks at the six sample points pass. The saved
candidate intervals contain these independently formed intervals.

The feedback at all six locations is replayed with the separately reviewed
exact reflected kernel. Its universal polynomial and mutation tests are in
`shell/kernel-review/`; this review does not call that replay a second
independent derivation of the kernel. The additional 138 ordinary-time
velocity conversions are checked against division by ell. Every exact-response
component enclosure has at least the stated conservative radius
`3E/sqrt(g_i)` for values and `E/sqrt(g_i)` for the scaled first derivatives.
All three native output files agree between normal and optimized execution.

## Additional consequences and scope

`INTEGRAL_RESIDUAL_LEMMA.md` derives a complete whole-cone integral-residual
bound for the same candidate; its inverse bound remains below 10^29 and
does not add another tail. `ANGULAR_DERIVATIVE_LEMMA.md` gives the intrinsic
unit-sphere covariant derivative bound 255E/4 by equivariant finite-dimensional
evaluation and the exact angular spectral cap. It is separate from physical
spacetime covariant derivative conversions.

The result remains conditional on the retained exact branch, first solution,
complete-source and formal physical completion premises. Neither this C1
certificate nor its integral residual is a pointwise second-derivative
Einstein residual or a nonlinear existence result.

## Portable replay

Use a fresh output path outside the capsule. `--normal-results` defaults to
the retained native normal-results directory; pass both result paths to check
newly generated normal and optimized outputs:

    python -I -S -B structure/response-review/check_response_coefficients.py \
      --repo-root /path/to/frozen/repo --capsule-root /path/to/capsule \
      --normal-results /tmp/native-normal --optimized-results /tmp/native-optimized \
      --output /tmp/new-coefficient-review.json

`COEFFICIENT_SAMPLE_CHECKS.json` records the examined input hashes and exact
check counts. It contains no environment-specific absolute input path.
