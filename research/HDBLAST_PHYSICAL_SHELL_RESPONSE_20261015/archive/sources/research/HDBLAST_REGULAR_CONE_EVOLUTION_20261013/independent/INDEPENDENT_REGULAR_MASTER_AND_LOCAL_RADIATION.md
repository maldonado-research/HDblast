# Independent regular-master derivation and local coupled radiation

Date: 13 October 2026. Independent internal mathematical derivation, not external peer review.

## Scope and assumptions

The registered background and its regular analytic branch are inputs; this note does not recertify their existence. Use five-dimensional signature `(-++++)`,

\[
g=dy^2+\rho^2\gamma,\quad \operatorname{Ric}\gamma=3\gamma,
\quad a=\rho'/\rho,\quad p=\phi',
\]

\[
a'=-\rho^{-2}-p^2/3,\quad p'=U_{,\phi}-4ap,
\quad \rho=y+O(y^3),\quad p=\kappa y+O(y^3),\quad\kappa\ne0.
\]

Assume the regular coefficients are real analytic in `y²`, as stipulated for the input branch, and the potential is the registered analytic potential. The statements below are local in ambient time and conditional on this actual background, rather than a flat test-field comparison. They neither identify a complete physical phase space nor prove global evolution after the perturbation encounters the shell.

## 1. Exact change of master variable

Define the exact conformal radius by `dr/r=dy/rho`, `r/y→1`. Analyticity gives

\[
s=r^2,\quad \rho=r\Omega(s),\quad dy=\Omega(s)dr,
\quad\phi=\Phi(s),\quad p=2r\Phi_s/\Omega,
\quad\Omega(0)=1,\quad\Phi_s(0)=\kappa/2\ne0.
\]

In ambient coordinates `X^A=(T,x^i)`, `s=η_AB X^A X^B`, the metric is `g=Ω(s)²η`; `E=X^A∂_A`. Near zero both `Ω` and `Φ_s` are nonvanishing. Their real analytic continuations to negative `s` solve the background equations there: the Cartesian equation residuals are analytic and vanish on the open positive-`s` region.

Starting from the supplied master equation, with `D=Box_gamma+4`,

\[
-(\alpha P_y)_y+qP=wDP,\quad
\alpha=(\rho^2p^2)^{-1},\quad q=2/(3\rho^2),\quad
w=(\rho^4p^2)^{-1},
\]

multiply by `−Ω²/α` and change to `r`. The result is

\[
P_{rr}-\left({4\over r}+2rL\right)P_r
 +r^{-2}(\Box_\gamma+4)P-{8\over3}r^2\Phi_s^2P=0,
\qquad L={\Omega_s\over\Omega}+2{\Phi_{ss}\over\Phi_s}.
\]

Now put `P=r⁴F`. The `r⁻²F` terms cancel exactly; using
`Box_eta=∂_r²+4∂_r/r+Box_gamma/r²` gives

\[
\boxed{\mathcal LF=\Box_\eta F-2L(s)EF
 -\left(8L(s)+{8\over3}s\Phi_s(s)^2\right)F=0.} \tag{1}
\]

This operator has smooth analytic lower-order coefficients and Minkowski principal symbol throughout `s=0`. It is normally hyperbolic there. With `V(r)=Φ_s(r²)`, the precise algebraic identity checked independently is

\[
-(\alpha P_y)_y+qP-wDP=-{1\over4\Omega^2V^2}\mathcal LF.
\]

The longitudinal reconstruction is equally regular:

\[
\Psi={sF\over\Omega^2},\quad N=-2\Psi,\quad
\boxed{h_{AB}=2sF\eta_{AB}-6F X_AX_B,\qquad
\chi=-{3\over2\Omega^2\Phi_s}(4+E)F.} \tag{2}
\]

No division by `s` remains. In particular smooth `F` produces smooth Cartesian perturbations through the cone, although the original radial gauge variables use a chart which degenerates there.

There is also a regular scalar conjugation. With `f=Phi_s` and

\[
A={1\over\Omega f^2}>0,\qquad A_s/A=-L,\qquad \psi=\sqrt A\,F,
\]

one has the exact identity

\[
\mathcal L(A^{-1/2}\psi)=A^{-1/2}(\Box_\eta-V_{\rm can})\psi,
\quad
V_{\rm can}=3L-2sL_s+sL^2+\tfrac83sf^2.
\]

This is a useful scalar PDE energy representation. By itself it does not equate this scalar action or its norm with the covariant constrained gravitational phase space.

## 2. Why these fields solve all bulk linearized equations

This is a sufficiency statement, without any claim that all perturbations have this form. On `s>0`, write `H=rho⁻²` and `D=Box_gamma+4`. The definition of `chi` and the master equation give

\[
\Psi'=-2a\Psi-p\chi/3,\qquad
\chi'=(p'/p)\chi+(3HD/p-2p)\Psi. \tag{3}
\]

For general functions of the de Sitter variables the unreduced longitudinal momentum, Hamiltonian and trace-free tangential equations are respectively

\[
D_\mu[-3(\Psi'-aN)-p\chi]=0,
\]

\[
12a(\Psi'-aN)+3HD\Psi-p(\chi'-pN)+U_{,\phi}\chi=0,
\]

\[
\left(D_\mu D_\nu-\tfrac14\gamma_{\mu\nu}\Box_\gamma\right)
(N+2\Psi)=0.
\]

They hold identically after (2)–(3). This direction never divides by a scalar gradient or a trace-free Hessian, so it does not omit exceptional harmonics from the sufficient family.

The full scalar residual is

\[
\mathcal S=\chi''+4a\chi'+H\Box_\gamma\chi-U_{,\phi\phi}\chi
 +6p\Psi'+4U_{,\phi}\Psi.
\]

Differentiate (3), use `H'=-2aH` and the background equations, and collect the independent `chi` and `Psi` coefficients. The coefficient of `chi` reduces to

\[
p''/p+4ap'/p-U_{,\phi\phi}-4H-4p^2/3
=-4a'-4H-4p^2/3=0.
\]

The coefficient of `D Psi` cancels. The remaining coefficient of `Psi` is
`−4p'−16ap+4U_,phi=0`. Thus `S=0`.

Let `E_AB=δ(G_AB−T_AB)`. The equations already established set the normal and mixed components to zero and leave only `E_mu nu=f g_mu nu`. Linearized Bianchi plus the scalar equation gives `∇^A E_AB=0`. Its normal component is `−4af=0`; because `a~1/y` and is nonzero near the cone, `f=0`. This proves the remaining tangential trace equation.

The proof extends to smooth master solutions across the cone, not merely their exterior restrictions. One precise analytic-coefficient argument is as follows. Each Cartesian Einstein/scalar residual after (2) is a linear differential operator `K` on `F`, of order at most two/three, with analytic coefficients across `s=0`. Since `T=constant` is noncharacteristic for (1), differential division by (1) eliminates all derivatives having two or more `T` derivatives, expressing

`K = B L + R`,

where `R` contains at most one `T` derivative. On any exterior open subset every master solution solves `K F=0` by the calculation above. Arbitrary Cauchy jets there imply that every coefficient of `R` vanishes. Those coefficients are analytic; they consequently vanish across the connected small Cartesian neighborhood. Thus `K=B L` there. This coefficient identity applies also to smooth, nonanalytic solutions of (1).

For the explicit construction below, an alternative proof uses only analytic local Cauchy data: Cauchy–Kowalevskaya produces an analytic germ across a chosen cone point; its full equation residuals vanish on the exterior and hence throughout the germ. Finite propagation identifies the same germ inside the solution from compact smooth cutoff data. The differential-operator argument proves the more general local smooth statement.

## 3. Smooth local Cauchy evolution and control of the remote shell

Choose `t0>0` sufficiently small for the cone sphere `T=R=t0` to lie inside the analytic extension neighborhood. Equation (1) has the form

\[
F_{TT}=\Delta_{\mathbb R^4}F-2L(s)(T F_T+x^iF_i)
 -\left(8L(s)+{8\over3}s\Phi_s^2\right)F.
\]

The usual local Cauchy theorem for smooth linear normally hyperbolic operators gives a unique smooth solution for compact smooth data `(F,F_T)` on `T=t0`, smooth dependence and finite propagation at Minkowski light speed. It applies here because the principal part is exactly Minkowski and all coefficients are smooth. Coefficients can be smoothly extended away from the required neighborhood to apply a global auxiliary Cauchy theorem; finite propagation makes the local answer independent of the chosen extension. This extension is a proof device, not a claimed physical background outside the registered neighborhood.

For a concrete support arrangement let the analytic coefficient interval contain `|s|<sigma`. Choose `t0>0` with `13 t0²<sigma` and `4 t0<r_b`, and put the initial compact support inside `t0/2<R<3t0/2`. Evolve on `−t0<T<3t0`. Finite propagation gives `R≤3t0/2+|T−t0|≤7t0/2` on the support, and consequently `|s|≤max(R²,T²)<13t0²<sigma`. Thus the entire evolving support remains inside the actual analytic neighborhood and strictly inside the shell

\[
R_b(T)=\sqrt{r_b^2+T^2}.
\]

This slab contains a neighborhood of the cone vertex `T=R=0` as well as the finite cone sphere used below. The initial angular data avoid `R=0`; subsequent evolution is performed in smooth Cartesian coordinates and remains smooth when the support reaches the origin. Every derivative in (2) has the same support property. Thus the reconstructed metric/scalar perturbations vanish on an open neighborhood of the shell throughout this slab. Set the shell displacement to zero and use equal or opposite copies in the doubled interiors. All homogeneous linearized shell conditions then hold, including the moving-shell equations, since every relevant field and derivative vanishes near the shell.

This constructs actual smooth coupled bulk-plus-shell solutions for a short time interval before boundary contact. It is not a proof of boundary evolution once the support reaches the shell. The Cauchy data of the reconstructed fields and their velocities are smooth and compactly supported; their ordinary local Sobolev/finite-energy norms are finite. This does not identify a global reduced gravitational energy norm or a completed gauge quotient.

## 4. Explicit nonzero null symplectic flux from actual Cauchy data

Use a real normalized `S³` harmonic `H_l` with

\[
-\Delta_{S^3}H_l=\lambda_l H_l,\quad
\lambda_l=l(l+2),\quad l\ge2,\quad\int_{S^3}H_l^2d\Omega_3=1.
\]

Let `eta(R)` be a compact smooth radial cutoff on the annulus above, equal to one near `R=t0`. Prescribe

\[
(F_1,F_{1T})|_{T=t_0}=(a_0\eta H_l,0),\qquad
(F_2,F_{2T})|_{T=t_0}=(0,b_0\eta H_l),\quad a_0b_0>0. \tag{4}
\]

These data are smooth in Cartesian coordinates since their support avoids `R=0`. Rotational invariance preserves their `l` sector. The constant and Hessian-degenerate de Sitter subspaces have only spatial `l=0,1`; therefore `l≥2` also avoids those exceptions throughout this evolution.

Use null coordinates `u=R−T`, `v=R+T`. At the cone sphere `u=0,v0=2t0`, put

\[
C={3\over2\Phi_s(0)}={3\over\kappa},\qquad L_0=L(0).
\]

All radial derivatives of the initial cutoff vanish there. The master equation determines the only needed second time derivatives:

\[
F_{1TT}=a_0[-\lambda_l/t_0^2-8L_0]H_l,
\qquad F_{2TT}=-2L_0t_0b_0H_l.
\]

Since `∂v=(∂T+∂R)/2` at fixed `u`, and `E=v∂v` on the cone, (2) gives

\[
\chi_1=-4Ca_0H_l,\quad\chi_2=-Ct_0b_0H_l,
\]

\[
\partial_v\chi_1=Ca_0\left({\lambda_l\over2t_0}+4t_0L_0\right)H_l,
\quad
\partial_v\chi_2=-Cb_0\left({5\over2}-t_0^2L_0\right)H_l.
\]

Consequently

\[
\boxed{\chi_1\partial_v\chi_2-\chi_2\partial_v\chi_1
=C^2a_0b_0\left(10+{\lambda_l\over2}\right)H_l^2}
\quad\text{at }v=v_0. \tag{5}
\]

All background-drift terms cancel. For `l=2` the coefficient is `14C²a0b0`. This identity comes from actual solution jets constrained by (1); it does not prescribe arbitrary candidate null traces.

The inherited action-derived current, with its stated Lee–Wald convention, has the smooth-cone limit per interior

\[
\omega|_{u=0}={v^3\over8}
 (\chi_1\partial_v\chi_2-\chi_2\partial_v\chi_1)
 \,dv\,d\Omega_3. \tag{6}
\]

Indeed smooth (2) implies the differentiated asymptotics required by that current limit, while `rho⁴ cosh³t (v/rho) dt` tends to `v³dv/8`. The scalar–metric contribution vanishes. There is no assumption of a fixed finite de Sitter-time limit here.

By smoothness the sphere integral of (5) remains strictly positive on some interval `[v0−epsilon,v0+epsilon]` inside the time slab. Equation (6) therefore gives a strictly nonzero integrated null symplectic flux. The same-parity doubled flux is twice (6); the outward-from-exterior orientation reverses its sign. Antisymmetric symplectic flux is not a positive quadratic energy; choosing `a0b0>0` only fixes the sign of this particular ordered pair.

## 5. Local gauge check and limits of the conclusion

One can also verify that the fields above are not pure diffeomorphism perturbations locally. Since `Phi_s≠0`, define the background scalar function `J` implicitly by

\[
J(\Phi(s))={4s\Phi_s^2\over\Omega^2}.
\]

The covariant scalar `K=(∇phi)²−J(phi)` vanishes identically on the background, so its first variation is gauge invariant. On the cone, `h^{AB}∂_Aphi∂_Bphi=0` for (2), and `J_,phi=4Phi_s(0)`. Hence

\[
\delta K=4\Phi_s(0)(E\chi-\chi).
\]

At the initial cone sphere the two perturbations satisfy

\[
\delta K_1=6a_0(\lambda_l+4+8t_0^2L_0)H_l,
\quad
\delta K_2=6t_0b_0(-4+2t_0^2L_0)H_l.
\]

For sufficiently small positive `t0` both are nonzero. This excludes pure gauge locally and supplements, rather than replaces, the exceptional-sector restriction. It does not define the allowed gauge group or corner charges of a global phase space. A nonzero finite-segment presymplectic flux alone would not have sufficed for that exclusion, since gauge transformations may carry endpoint charges.

The result is therefore a local existence and nonzero-flux statement for actual smooth coupled perturbations, with a short-time remote-shell realization. It does not establish global radiative scattering, global shell well-posedness, completeness of the longitudinal family, a completed physical gauge quotient, nonlinear stability, a cosmological origin, or priority over earlier work.

## Reproducibility

`check_regular_master_independently.py` independently transforms the original radial differential expression, checks regular reconstruction, verifies the unreduced constraints and full scalar residual, and checks the Cauchy-jet Wronskian and current scaling exactly with SymPy. Its commuting formal `lambda` in the scalar check stands for `Box_gamma+4`; the resulting identity is polynomial in that operator and uses no harmonic inversion.

The checker uses explicit pass/fail conditions rather than Python `assert`. `run_verification.py` recorded 22 successful checks in each of normal and optimized mode. Three deliberate mutations (master potential, scalar source, null-jet derivative) each exited with failure in both modes. All eight outcomes are recorded in `VERIFICATION_RECEIPT.json` and the per-run JSON files. These receipts verify exact algebra, not the functional-analytic Cauchy theorem or completeness.

Public capsule locations: the preceding run_verification/receipt names refer
to the preserved authoring audit. Here the original checker is under
`independent/original/`, the reviewed CLI adaptation is registered in
`independent/INDEPENDENT_ADAPTATION.json` and
`independent/ADAPTATION_VERIFICATION.json`, and the expected public result is
`results/independent.json`. The top-level runner regenerates all outcomes
outside the capsule. The 22 rows comprise 18 identities and four wrong-formula
controls; they are not 22 separate positive identities.
