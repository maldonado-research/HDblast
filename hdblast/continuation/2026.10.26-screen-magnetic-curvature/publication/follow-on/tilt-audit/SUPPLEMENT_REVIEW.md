# Final independent review of the scalar-tilt supplement

**PASS within the stated scope.** No significant equation, implementation, or
scope fault was found in `SCALAR_OBSERVER_TILT_REDUCTION.md` and the independent
direct-coordinate fixture. This is internal mathematical/code review, not
external specialist review.

The supplement's equations agree with the independent derivation recorded in
`TILT_AUDIT.md`: the full relative-frame contraction is
`Delta H''=2(VP-WQ)A`, with the stated `P,Q`; the scalar radial/time boost
cancels. The exact correlated registered amplitudes give `P=Q=0` without field
equations. The text restricts this conclusion to the same-event, same-background,
standing-quadrupole scalar tilt class and does not infer arbitrary observer
independence, a completed clock construction, gauge invariance, or a new
actual-field sign.

The direct oracle was reviewed independently of its PASS output. Its quartic
stereographic expansions contain the derivatives needed for the full metric
two-jets, including the angular Hessian metric slot. The background is diagonal
at the selected point, so the implemented first inverse-metric variation is
exact. The Riemann ordering matches the stipulated curvature convention; the
five-dimensional Schouten subtraction uses `P=(Ric-R g/8)/3`. Background Weyl
is checked to vanish before conversion of its first variation with the
background frame, which justifies omitting first frame changes multiplying
`C0`. The ambient normal and triad have the required outward orientation, and
the factors `v=grad_coordinate(Y)/2`, `h=Hess_coordinate(Y)/4` correctly convert
the stereographic sphere metric at the chosen point. The slotwise relative
frame contraction retains observer and screen variations and the full second
parameter cross factor.

Independent isolated normal and optimized executions both pass and give equal
reports: three arbitrary populated scalar metric two-jets, four relative tilt
cases per fixture, and all nine screen components per case. Each fixture
explicitly rejects zero `P` or zero `Q`, preventing a generic-zero shortcut.
These finite checks support the all-point analytic derivation; the supplement
does not present them as an all-point proof by themselves.

Reviewed input SHA256 values:

- `SCALAR_OBSERVER_TILT_REDUCTION.md`:
  `193dd2bef26fbb2d94b3b6e922d0e9290f2a25402886ee3fc1002f782439e482`
- `tilt/check_direct_tilt.py`:
  `debe39d9f938621739002d3bfde6f01ae61a14a10b0ef91069aa5db45998ded9`

The independent replay evidence is `DIRECT_REVIEW_NORMAL.json`,
`DIRECT_REVIEW_OPTIMIZED.json`, and `SUPPLEMENT_REVIEW.json` in this directory.
Neither the frozen scientific capsule nor the completed anchoring note was
modified. No actual-field numerical evaluation was performed.
