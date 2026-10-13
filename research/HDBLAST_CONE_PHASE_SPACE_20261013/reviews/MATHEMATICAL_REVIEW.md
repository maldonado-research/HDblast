# Independent bounded mathematical review of the integrated manuscript

13 October 2026 UTC. Reviewer: cone-domain audit agent. This is an internal
mathematical review by a separate agent process, **not outside peer review**.
It reads the integrated Markdown and TeX and compares their mathematics with
the four technical notes. No checker tests or PDF typesetting were rerun.
The exact reviewed inputs and candidate disposition are in
`INTEGRATED_MANUSCRIPT_REVIEW.json`.

## Bounded review findings

The covariant potential uses raised covariant variations, so its sign is
consistent with the displayed Lee–Wald convention. The conformal nonlinear
family correctly gives the GHY/canonical equality in (6). The current curl
and momentum constraint give the `-2p` coefficient in (8), and the
tangential coefficient is `+6`. Mixed spectral values are retained in
(10)–(12); the radial derivative has `(Lambda1-Lambda2)` and the tangential
divergence has the opposite sign. There is no unjustified replacement by
a common-Lambda Wronskian for packets.

The even endpoint condition cancels mixed current without an added physical
shell kinetic field. The odd kinetic action/current includes two gravity
copies and one shell. Its divergence equals the positive increasing-y
shell flux; adding that current to the integrated bulk leaves the lower
cone term with the sign in (16). These conclusions retain the constrained
generic-sector scope and the prior graph-action input.

The Liouville Green form has the orientation in (18). The sufficient packet
criterion in (19) uses `F in C(I;H^2_x(L^2Omega))` and
`G in C(I;H^1_x(L^2Omega))`, so angular-integrated products vanish uniformly
on a fixed compact time interval. The manuscript correctly avoids a
pointwise derivative trace for arbitrary `H^1` form vectors and does not
extend the uniform estimate to intervals translated to infinite time.

The exact conformal radius gives the correct null Jacobian:
`det partial(y,t)/partial(u,v)=Omega/(2r)` and
`partial_y u=u/rho`. Thus (23) has the correct density factors and needs
both radial and tangential components. The coefficient of the per-copy
null current is `v^3/8`; the doubled coefficient in (24) is `1/4`.
The opposite outward lower-boundary orientation is explicitly stated.
The null calculation is conditional on actual on-shell solutions with
uniform differentiated characteristic expansions. Candidate leading data
are not presented as constructed shell-compatible solutions or as proof
of a global longitudinal gauge.

The exact test-scalar fields in (25) satisfy the Minkowski wave equation.
Compact plateau initial data give the local fields by finite propagation
speed, with finite global inertial energy. The fixed-time and moving-time
fluxes (26)–(27) have the correct measure and noncommuting limits. This
example is explicitly separate from the registered coupled gravity model.
The powers and cutoff costs in (27a) are correct; the lapse weighted radial
and inertial gradient energies have convergence thresholds `3/2` and `1`.
Their distinction does not establish inequivalence of completed gravity
phase spaces.

The tensor Weyl identity (28) follows from the warped Einstein slice and
is invariant because the background Weyl tensor vanishes. A bounded-frame
component divergence excludes a C2 vertex extension; a vanishing Lorentzian
curvature contraction cannot evade that argument. The radial cutoff costs
in (29) are O(epsilon^3), and inertial initial energy costs are
O(epsilon^2). Angular TT metric and extrinsic-curvature initial data remain
tracefree/divergence-free after radial cutoff. On the time-symmetric slice,
angular Ricci isotropy and unchanged scalar data preserve both constraints.
These are initial-data approximants, not separated zero-mass eigenmodes.
The manuscript correctly retains the massless mode in the completion.

The charge normalization (30), small-area vertex estimate (31), smooth
translation vector and the cancellation `a+F=-integral p^2/3` are correct.
No finite-null-cut charge theorem, charged exceptional mode, or global
gauge coverage follows from those bounded-jet and local statements. The
potential/corner convention, generic-sector restrictions, exceptional
harmonic qualifications, and unproved physical existence/completeness
questions are consistently retained.

## Supporting-note correction found during review

The original supporting gauge note said that the ell0 residual
`ell0[r c^A-(c_B X^B)X^A/r]` was continuous at the vertex. That is false
in the full ambient topology. Take `c=(0,1,0,...)`,
`X^1=epsilon`, `X^0=sqrt(epsilon^2-epsilon^6)`, so
`r=sqrt(eta_AB X^A X^B)=epsilon^3`. Its spatial component is
`ell0(epsilon^3-epsilon^-1)`, which diverges while `X -> 0` within the
spacelike cone wedge. Bounded finite-t rays hide this correlated approach.
The weaker original manuscript statement that this residual is not C1 was
already correct. The corrected supporting note and manuscript now remove
the continuity claim and include this explicit near-null counterexample.
Its exclusion from the smooth gauge class remains valid. The authoring
agent added affected exact algebra controls; this reviewer did not rerun them.

Three small clarity changes were requested: identify the momentum
constraint, rather than the first printed relation in (7), as the reduction
used in (8); state explicitly that the test-field energy example establishes
no gravity completion equivalence or inequivalence; and distinguish bounded
Cartesian metric components from O(r^-1) derivatives and velocities in the
tensor energy estimate. All three changes were applied and verified in both
the final Markdown and TeX. The JSON receipt pins the corrected candidate.

No mathematical error was found in the manuscript's current, shell balance,
conditional cone-flux coefficient, exact comparator, or completed-energy
and local gauge conclusions. Literature priority and primary-source scope
were handled by the separate literature audit; this review did not repeat
that search or independently verify every bibliographic assertion.
