# Independent numerical-input review

16 October 2026. Accepted within the conditional scope below. This is an
internal computational/scientific review, not external specialist peer review.
Only this new review directory was written; producer and frozen files were not
edited. The final producer and review artifact identities are in
`OUTPUT_PINS.json`.

The reviewed results are continuous bounds for selected initial profiles,
explicit first-master/raw-gauge comparison errors, a common affine background
approximation, and a refined inherited response majorant. They do not evaluate
the second response or its continuous numerical residual. The sharp selected
angular evaluation constant and invariant-space closure have separate reviews.

The review first read the actual frozen physical initial-jet, relative-harmonic
completion, first harmonic-pair, first-master reflected-derivative, first-gauge
wave, and full second-response constructions. Previously accepted smooth
actual background/master/gauge solutions and authenticated coefficient reports
remain conditional premises. This review does not re-prove the regular branch
or the original PDE existence theorems.

`check_independent_initial_algebra.py` uses standard-library sparse Laurent
polynomials with independent jet symbols. Eighteen comparisons are all-point
identities rather than point samples. They recover both initial gauge-velocity
orders directly from the exact ADM prescription, retain the second lapse
square, and independently recover the mixed metric velocity

    (g_ty)_2,t = rho^2[-26 U U_y - 3 f_y/2 - 2 a U^2 - 3 a f].

The first gauge acceleration is exactly
`xi^y_tt=-6 rho^2(U_y+4aU)`; the temporal and angular accelerations vanish.
It changes the first metric velocity in the mixed radial/time component while
preserving its physical initial data. The scalar compensator cancellation and
the actual three sector denominators `-6,10,42` are verified independently.

The producer's cancellation-preserving Hamiltonian and scalar-velocity
expressions also agree identically with the original physical-Cauchy-affine
sources. The independent reduction uses

    u_y=-2au-pv/3,  p_y=p(h-4a),  a_y=-rho^-2-p^2/3,
    E_y=E/rho,  d_y=d(3a-h),  d=3/(rho p),  v_t=E(v-du).

This detects the radial momentum term that distinguishes physical-Cauchy and
spacetime affine paths. It does not silently set the second scalar or momentum
data to zero. Normal and optimized Python give identical reports; eight
mutation controls detect missing factors, compensator signs, the physical
affine term, lapse square, mixed conversion and fourth packet derivative.

The fourth packet derivative is independently derived as

    Q''''=(s^8-20s^7+136s^6-384s^5+420s^4-120s^3) exp(4-s).

Both `s^8` plus and minus its polynomial have nonnegative coefficients after
`s=u+4`. Thus `|Q''''|<16777216/49` globally, with the flat zero extensions.
This independently corroborates the retained packet certificate. The initial
scalar spatial derivative reaches this fourth derivative through `f_yy`;
third packet derivatives alone would not suffice.

The initial evaluator's rational outward interval arithmetic, normalized
Taylor recurrences, polynomial evaluation, exponential enclosure and complete
498-cell/tail cover were inspected. Actual background values come from the
authenticated shell box and ODE image, and derivatives are calculated from the
actual ODE. The negative logarithmic alpha derivative and positive enclosure
are checked explicitly. The scalar and metric sector expressions use shared
actual functions; dependency loss widens enclosures without invalidating them.
The final caps use `<=`, correctly named `component_caps`. Unnormalized
`Z2,H4` coefficients are correctly named sector amplitudes.

The frame conversion is correct. The radial frame is parallel, so
`partial_r=-ell rho partial_y`. In the normal mixed metric derivative the
connection cancels the displayed `-2aU^2`; converting back to the actual
scaled partial time velocity restores that term. The angular factors
`1,1/2,7/12` and gradient factors `0,1,5/2` are valid over the full sphere.
A separate reviewer execution of the complete initial evaluator at 64 cells
passes in both Python modes with byte-identical reports. This is a coarser
whole-domain replay, not an independent reimplementation of all tight
512-partition interval arithmetic.

For the master comparison, its nonzero true initial velocity, dynamic shell
lift through third derivative, and reflected Duhamel derivatives are all
retained once. Reflection shell terms and flat-corner joins are explicitly
covered by the inherited differentiated integral identities. No derivative is
obtained from a C0 error bound.

For the gauge comparison, differentiating the prefactor `r` requires the
additional product term. The final conservative cap
`c[5 Q_k+k Q_(k-1)]` includes it. The true error equation retains homogeneous
physical mixed rows. Its differentiated equations retain coefficient
commutators, initial acceleration `R0`, its radial derivative, the third
initial jet `A0 R0+Rx0`, and the complete radial residual in third-radial
recovery. These points were checked analytically and with independently
manufactured polynomial fields, including the corresponding failure controls.

`check_independent_error_certificate.py` imports no producer program. It
recomputes all four master errors, six mixed source-defect bounds, six complete
three-channel gauge residuals, three differentiated resolvent vectors, and ten
mixed gauge errors with exact fractions. Its matrix inverse uses an adjugate,
independent of producer elimination. Final normal/optimized outputs are
identical. All 14 producer wrong-formula/input processes were also independently
re-executed and reject without producing reports.

The common affine background approximation defines one actual rational
function per state and derives all coefficients from those functions. Actual
and approximate derivatives are formed separately before subtraction. The
scaled C4 bounds are therefore derivative errors, not differentiated interval
selections. The final checker rejects an empty partition and verifies complete
collar coverage. Shell uncertainty and low affine order still dominate these
errors; the approximation has not been used to tighten the first-field errors.

The refined response join is accepted independently. A bound on every 5x5
metric component gives a positive Frobenius bound five times as large; adding
the scalar gives the three inherited initial H bounds. The mixed free-wave
estimate and retained source bounds then give exactly the reported
`X<1.03e33`. All 15 exact tail counts and their minimality checks agree with an
independent rational calculation. The manuscript's four displayed master-error
roundups and response products are strict conservative roundups. Selected
pointwise factors are conditional on the separate angular proof.

No unresolved numerical-input defect remains in the pinned artifacts. These
bounds certify neither common-gauge quadratic source approximation errors nor
a full spacetime waveform, invariant force, nonlinear remainder, or novelty.
