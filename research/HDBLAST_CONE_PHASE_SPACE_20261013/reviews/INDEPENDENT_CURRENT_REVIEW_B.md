# Bounded mathematical review of the independent current derivation

13 October 2026 UTC. Reviewer: independent cone-domain audit agent. This
review reads the root note and script; it **does not rerun their tests**.
Exact hashes of the note, script and supplied reports are recorded in
`ROOT_CURRENT_REVIEW.json`. No frozen/source file was modified.

**Disposition: no local sign or normalization correction required.** This
is internal bounded mathematical review, not a global phase-space theorem.

The nonlinear tangential potential can be independently recovered using
`h^AB=2 delta psi g^AB+2(delta n-delta psi)n^A n^B`. Its tangential
divergence is `2 grad^mu delta psi-2(delta n-delta psi)grad^mu n`.
Subtracting the trace derivative, for trace `2 delta n+8 delta psi`, gives
the note's potential, including the sign of `(delta n-delta psi)D^mu n`.
The scalar contribution has the displayed negative sign.

The radial contraction similarly gives the displayed
`-rho^4 exp(-n+4psi)[4 delta psi_y+4(a+psi_y)(delta psi-delta n)]`.
Restricting the commuting family to `n=-2psi` gives the curl coefficients
`-12 W(Psi)-6p(Psi1 X2-Psi2 X1)`. The momentum constraint yields
`W(Psi)=-p(Psi1 X2-Psi2 X1)/3`, hence the reduced cross coefficient
is correctly `-2p`.

The first-order equations give
`9 alpha[Lambda1 P1 P2_y-Lambda2 P2 P1_y]`. The root script retains
mixed-eigenvalue terms in the bulk conservation check, as required for
packets. The even eigenparameter shell condition cancels this flux.
The odd doubled flux is `18 c_o(Lambda2-Lambda1)P1 P2`, equal to the
stated shell-current divergence. Integrating bulk conservation gives
`div bulk=j_y(cone)-j_y(shell)`; adding the shell divergence leaves the
cone term with the printed orientation. This is conditional on the prior
shell-action input and does not independently derive that action.

For `y^2=uv`, the density Jacobian is `1/(2y)` and
`partial_y u=u/y`, `partial_t u=-u`. Therefore
`J^u=[(u/y)j_y-u j_t]/(2y)` is correct. The two leading scalar pieces
give `(u+v)^3 W(F)/8`, and its null limit is `v^3 W(F)/8`. Dropping
the tangential piece gives half the answer. This agrees with the direct
scalar null current and the independent domain audit's master calculation.

The note makes the null result conditional on actual on-shell fields with
uniform differentiated expansions on `0<a<=v<=b`, and correctly avoids
claiming global shell-compatible packet existence. Those expansions give
`O(y^2)` errors in the integrated null density on this interval. An exact
conformal-flat formula needs the conformal radius `r`, but the `y`-based
null chart is also smooth when the background has smooth even coefficients:
`g=f(uv) eta+[1-f(uv)](v du+u dv)^2/(4uv)`, `f=rho^2/y^2`.
The finite curvature correction to `g_uu` on `u=0` does not change the
limiting null density `sqrt(-g)g^(uv)=v^3/8`.

The supplied report's identity and deliberately wrong formula counts are
recorded in the JSON receipt. Their successful runs check algebra, while
this review independently checks the potential contraction and signs.
Neither supplies a smooth global gauge quotient, a registered hyperbolic
existence theorem, characteristic remainder estimates for actual solutions,
or physical completeness.
