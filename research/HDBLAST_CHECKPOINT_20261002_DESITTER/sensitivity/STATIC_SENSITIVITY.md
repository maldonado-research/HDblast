# Resolved classical static source susceptibility

The inherited corrected +1 branch at delta=0.001 has a small **positive**
metric-junction derivative with respect to ell=log10(abs(phi_h-1)):

    E1_ell = 9.63825e-14, empirical absolute error scale 3.34e-18.

This resolves the sign and order of the component that the archived residual
finite difference could not resolve. Its archived value −6.9388939e-12 is
about −72 times the resolved result. The difference does not undermine the
archived nonsingularity conclusion; it does change its tiny inverse-matrix
entries. This calculation is a differential classical response to supplied
source values, with no physical quantum parameters or quantum source chosen.

## Registration, inputs, and derivation

`REGISTRATION.md` and `INPUT_PINS.json` were written and hashed at
2026-10-02T06:12:59.790789+00:00, before the first radial integration.
Four inherited files are pinned, including the actual September 27 background
shooting coordinates. All calculations hold

    delta=.001, c=.5975949350280132,
    e_h=-1.3366923651588084e-23, y_b=30.276680395885563.

The shifted variable e=phi−1 avoids rounding a cone displacement of order
1e-23 out of existence. No root search, physical parameter scan, or tuning was
performed. The fixed archived coordinates have residuals about
(2.78e-17,−2.11e-16) in the finest independent implementation; their estimated
Newton correction is (Delta ell,Delta y)=(-6.23e-13,4.69e-13).

The seven integrated variables are (R,e,w,R_ell,phi_ell,w_ell,I), where
I'=R R_ell. The regular initial series is

    R=y−U_h y³/36+(U_h²/4320−U_phi,h²/750)y⁵+O(y⁷),
    e=e_h+U_phi,h y²/10+U_phi,h(U_phiphi,h/280+U_h/630)y⁴+O(y⁶).

Initial variational values are analytic derivatives of this series, and the
initial integral includes its y⁵ and y⁷ terms. The verifier independently
substitutes the series into the acceleration and scalar equations with SymPy.

For z=R_ell/R, f=phi_ell, and v=z'=q_ell, varying the radial constraint gives

    2qv=−2z/R²+(w f'−U_phi f)/6.

Since q'=−1/R²−w²/3, define S=v+w f/3 and obtain exactly

    S'+4qS=−2z/R²,
    R_b⁴ S_b=−2 integral_0^yb R R_ell dy.

The regular cone contributes zero to R⁴S at y=0. The exact off-shell formula is

    E1_ell = −2 I_b/R_b⁴ − E2 phi_ell/3.

At the scalar junction E2=0 the second term vanishes. Thus no exact zero is
expected: a weighted interior response fixes the small derivative. This
identity was separately checked by the theory reviewer before execution and
is also checked symbolically in `verify_sensitivity.py`. Direct variational
subtraction and the integral agree to at most 1.39e-20 over every registered
run (5.04e-24 in the reference run). The off-shell correction was retained.

The moving-endpoint metric derivative is exactly
E1_y=−H²−w E2/3, reducing to −H² on shell. E2_y=w'+sigma_phiphi w/2=Bw.
Both forms and the differentiated endpoint constraint were checked.

## Numerical result and source convention

Let s=K rho and t=K J_phi denote independently supplied infinitesimal shell
source values. With residuals E1=q−sigma/6−s/6 and
E2=w+sigma_phi/2+t/2, the linear source solve is

    Delta u=J_vac^(-1) diag(1/6,−1/2) (s,t),   u=(ell,y_b).

The reference run gives

    phi_ell = −1.93527082528243e-4  [error scale 3.35e-14],
    (H²)_ell = 3.21211710044819e-13 [error scale 1.11e-22],

    J_vac = [[ 9.63825146e-14, −5.92401479432896e-5],
             [−6.87932216814216e-4, −4.64606629798131e-4]].

Its determinant is −4.07532063438e-8 and condition number 16.9364147.
The coordinate response matrix, with columns (s,t), is

    [[ 1900.08211002,  726.815792646],
     [−2813.40733050, 1.18251450e-6]].

The empirical error scale of the small dy_b/dt entry is 4.09e-11. Unlike the
archived numerical-floor estimate −8.51e-5, the resolved derivative is
positive and approximately 1.18e-6.

Including the moving shell endpoint gives the main observable result:

    Delta phi_b = −7.3550543e-6 s −0.140658540041 t + O(source²),
    Delta H_b²  =  0.0371257908752 s +2.1785725e-10 t + O(source²).

The respective empirical absolute error scales for the four coefficients are

    [[7.60e-13, 1.28e-13],
     [1.39e-16, 5.40e-16]].

These scales are ten times the largest observed full-run envelope, including
coarse Radau comparisons. They are conservative numerical diagnostics, not
rigorous interval bounds or a theorem about integration error. The small
current-to-curvature response is resolved relative to that diagnostic; the
nonlinear O(source²) terms have not been bounded.

An independent endpoint constraint derivative gives

    Delta H² = [sigma sigma_phi/18−sigma_phi sigma_phiphi/24+U_phi/6]
                 Delta phi + sigma s/18−sigma_phi t/24.

Its reconstructed coefficients agree with the radial-response calculation to
within 1.81e-16 across the run set, including the tiny current coefficient. The
leading detuning expansion cancels the current contribution at quadratic
order; the resolved nonzero coefficient concerns the exact finite-detuning
background and is consistent with higher-order terms.

At the inherited eta_b≈−8.40527e-5 and H_b²≈5.92401e-5, useful first-order
sufficient response estimates are

    |Delta phi_b|/|eta_b| <= 0.087506 |s| +1673.457 |t|,
    |Delta H_b²|/H_b²     <= 626.700 |s| +3.67753e-6 |t|.

Both must be much less than one for this local small-response interpretation.
No actual source values are supplied here, so no small-response pass is
asserted. These frozen supplied-source derivatives are the coefficients needed
for a formal leading-loop source calculation. Finite-amplitude source feedback
or phi-dependent density requires its full common-action source derivatives.

## Refinements, controls, and limitations

Nineteen registered runs completed: 13 DOP853 variational runs, four Radau
variational runs, and two independent complex-step background integrations.
The methods span rtol=1e-10 to 3e-14, cone starts y0=1e-2 to 1e-5, and maximum
steps .1 and .05. The finest fixed-tolerance cone-start envelopes are

    phi_ell: 1.44e-18;  (H²)_ell: 7.48e-27;  E1_ell: 1.77e-22.

Complex-step h=1e-15 and 1e-20 gives E1_ell≈9.638251656e-14, differing from the
reference by 2.00e-21. Its phi_ell and (H²)_ell differences are 1.92e-17 and
6.44e-26. All 11 registered numerical gates pass; all numerical run records,
including the deliberately coarser results, remain in `RUNS.json`.

Three registered wrong-formula controls are detected:

1. Dropping the scalar junction derivative gives E1_ell≈−8.43118e-9 and fails
   the weighted identity by over 8e-9.
2. Reversing the current sign produces a unit source-normalized scalar
   residual error.
3. Omitting the moving endpoint terms changes the density-induced scalar
   response by about 0.368 and loses the density-induced curvature response.

There were no integration failures. The first symbolic-verifier execution
failed because substitution was attempted before polynomial expansion; the
unexpanded expression concealed the matching product. The saved failed
verifier and `CHECK.log` preserve that implementation failure. Expanding the
expression before substitution fixes it; `CHECK_v2.log` records the passing
run. The calculation and numerical outputs were unchanged by that correction.

This evidence resolves classical first derivatives near one inherited regular
branch. It establishes neither a quantum-corrected BVP nor time evolution,
nonlinear attraction, quantum stability, reheating, or a physical prediction.

## Reproduction

The computation requires Python, NumPy, SciPy, and SymPy. To preserve archived
outputs, choose new output filenames on replay:

    python compute_sensitivity.py --repo /path/to/HDblast --output /tmp/sensitivity-runs.json
    python verify_sensitivity.py --runs /tmp/sensitivity-runs.json --output /tmp/sensitivity-checks.json

`compute_sensitivity.py` verifies pinned input and registration hashes before
integrating; the checker records the input run hash and its own source hash.
An explicit output file is required and its parent directory is created as
needed. The repository is detected from script/current-directory ancestors if
`--repo` is omitted. `PORTABILITY.json` records these mechanical CLI changes;
`registered_sources/` retains the exact executed sources. The failed verifier is
retained as evidence and is not part of the reproduction command.
