# Higher-jet research targets

26 October 2026. Conditional planning review of commit `218d58bb6ba3889d60bd816aecc076f35f474e0d`, not a new theorem or numerical certificate. No capsule was changed and no numerical campaign rerun. Repository-relative source abbreviations: `F/` = `research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026/`; `C/` = `research/HDBLAST_FORCED_RESPONSE_BOUNDS_20261023/`; `N/` = `research/HDBLAST_NUMERICAL_RESPONSE_INPUTS_20261016/`.

**Recommended next target:** certify second coordinate jets on a specified exterior domain, then extend to a physical-side C2 approximation on the full retained domain with a smoother numerical candidate. The regularity issue found here does **not** require a new physical seed.

## Distinguish actual data from the approximation

The current full23 certificate bounds values and scaled first x/r derivatives by `1.1787044291e28`. The separate total TT sign at `(1/4,5/4)` is not a second-derivative or physical-observable certificate (`F/DERIVATION.md:43–70,98–150,152–189`).

The numerical packet cutoff matches derivatives through order four; registered profiles contain packet derivatives through order three. Thus the rational initial profiles have guaranteed C1 regularity, and second derivatives can jump at the upper join. The review's exact jump witness uses general polynomial inputs; it does not establish a nonzero jump in every particular retained slot (`F/initial/INITIAL_DERIVATION.txt:22–27,49–64`; `F/initial/full23_initial.py:151–227`; `F/initial/review/check_exact_adapter.py:150–181`).

Independently, the free velocity is continuous P1, while feedback uses independent continuous P1 displacement and derivative tables. Velocity slope jumps can propagate along knot characteristics. Smoothing the packet alone does not fix those knots; a higher-degree table alone does not fix the expression-defined displacement cutoff (`F/structure/response-review/SCOPE_REVIEW.md:19–25`).

Actual data use the exponential flat packet, not this artificial cutoff (`N/numeric/evaluate_initial_profiles.py:106–117`; `F/initial/full23_initial.py:322–339`). The retained construction proves qualitative smoothness by characteristic/trace-ODE induction under its background and actual-first-field premises, explicitly distinguishing this from numerical high-jet bounds (`C/inputs/FULL_SECOND_RESPONSE_CONSTRUCTION.md:310–368`). Flat actual data give a zero initial shell collar and qualitative corner compatibility at every order (`F/shell/FULL23_SHELL_AND_SUBSIDIARY.md:182–208`). Calling the actual data “only C1” would therefore be incorrect.

## First deliverable: exterior second jets

A focused target is the total TT two-jet at the existing apex or on a smaller closed region within its backward cone, `1+x<=r<=3/2−x`, `0<=x<=1/4`. The cone avoids the shell, TT initial data vanish, and the exact TT equations decouple (`F/DERIVATION.md:152–176`). This removes shell derivative estimates and much full23 initial-jet bookkeeping. Its numerical viability remains untested.

Use a fresh smooth exterior approximation, or retain the old candidate only on audited smooth pieces. The latter requires excluding characteristics from **all relevant P1 knots**, not merely cutoff joins. The free-wave set includes `r+x=a` and `r−x=±a`; derive additional feedback/reflection break curves from the actual formulas and partition integrals at moving knot crossings. Avoiding `r=66/64` at an evaluation point is insufficient. Neither option obtains a C2 error by differentiating the existing C1 number.

A concrete route is to set `W=U_x` in

`U_xx−U_rr=A U_x+B U_r+D U+F`.

Then

`W_xx−W_rr=A W_x+B W_r+D W+A_x U_x+B_x U_r+D_x U+F_x`,

with `W(0)=V0` and `W_x(0)=U0_rr+A0 V0+B0 U0_r+D0 U0+F(0)`.

A certified C1 estimate for W supplies U_xx and U_xr; the original equation recovers U_rr. The highest-order operator is unchanged, but derivative-input errors and initial acceleration offsets need a new ledger and proved inverse bound. The old omitted-source cap does not bound F_x (`C/inputs/FULL_SECOND_RESPONSE_CONSTRUCTION.md:343–395`).

For full23 initial comparison, U0_rr needs packet derivatives through five and V0_r through four. Current paired primitives and the actual packet oracle stop at four (`F/initial/full23_initial.py:284–347`; `N/numeric/evaluate_initial_profiles.py:106–117`). Extend the algebra and flat-tail checks. Inspect the existing fifth ODE derivative enclosure before requiring a sixth: it may already bound the needed fifth background derivative. Additional reconstruction orders depend on the chosen interpolant.

## Full-domain target and shell gates

A sufficient replacement packet cutoff matches through order five. One candidate is the degree-eleven polynomial `S5(s)=2772 integral_0^s t^5(1−t)^5 dt`, giving C2 registered profiles after at most three packet derivatives. Verify its joins exactly and certify a new actual-minus-candidate ledger. This changes an approximation to the same actual data.

The free wave needs displacement C2 and velocity C1. Use a globally C1 velocity reconstruction with derivative-error control; refining P1 alone does not supply that regularity. Prove sufficient derivative continuity and approximation errors for feedback and lift inputs as well. Existing one-sided cell bounds do not establish matching derivatives across cells.

The differentiated Robin law is `PN W_r=M PN W+M_x PN U+G_x`. M is constant in the retained model, so M_x=0; G_x and derivative-feedback errors still need certification. For `Lg=−integral_0^max(x−r,0)g(s)ds`, continuous second derivatives across the corner characteristic require g C1 and g(0)=g'(0)=0. Preserve the open zero collar or charge new candidate corner defects. Retain all five Dirichlet and eighteen Robin rows and use physical-side traces; do not differentiate an odd reflected source as though its shell jump were smooth (`F/shell/FULL23_SHELL_AND_SUBSIDIARY.md:158–200`; `C/inputs/FULL_SECOND_RESPONSE_CONSTRUCTION.md:343–395`). “Global C2” here means the retained physical half-domain, not a C2 metric across a distributional thin shell.

## Residual or observable: define the endpoint first

A pointwise **natural** Einstein/scalar residual needs complete second-order source terms, second x/r jets, second angular and mixed angular-coordinate bounds, and physical conversions. Use `partial_t=ell^−1 partial_x` and `partial_y=−(ell rho)^−1 partial_r`, differentiating conversion coefficients and including frame/connection terms. The intrinsic first-angular lemma alone is insufficient. Treat shell junction residuals separately. With the same exact first field q, the natural residual error is the linearized natural operator applied to the second-field error; reconstructed q adds quadratic-source errors. A reduced wave residual alone omits gauge/subsidiary obligations (`F/shell/FULL23_SHELL_AND_SUBSIDIARY.md:212–261`; `F/structure/response-review/SCOPE_REVIEW.md:73–83`).

If retaining a piecewise candidate, an alternative is a precisely declared weak/interface residual with its test-function norm, derivative jumps and shell terms proved. Do not call it a pointwise C2 certificate.

Before an observable campaign, specify observer, frame, anchoring, units and gauge behavior, and derive its complete error map. Curvature/tidal observables normally need second jets and quadratic first-field/frame corrections; an integrated quantity may need less. The registered nonzero TT amplitude is not itself a gauge-invariant signal (`F/DERIVATION.md:178–189`).

Acceptance order: authenticate inputs; prove candidate regularity/corner joins; certify derivative-input errors; prove the differentiated response/error estimate; perform physical and angular conversions; certify only the declared residual or observable. Preserve the frozen capsule and put new calculations in a separate follow-on capsule. No physical-origin, finite-amplitude or observational conclusion follows from this plan.
