# Independent local current and characteristic-limit calculation

13 October 2026. This root calculation was derived independently from the
nonlinear metric-family symplectic potential, then compared with the parallel
covariant-current derivation. It is internal mathematical review, not external
peer review. The algebra script checks the stated local identities; it does
not prove the existence or global admissibility of the hypothesized solutions.

## Conventions and direct potential contraction

Use the previously registered action with bulk density
`sqrt(-g) [R/2 - (d phi)^2/2 - U]`. For covariant metric variations,

\[
 \theta^A_g=\tfrac12\sqrt{-g}
 (g^{BC}\delta\Gamma^A_{BC}-g^{AB}\delta\Gamma^C_{BC}),
 \qquad \theta^A_\phi=-\sqrt{-g}\nabla^A\phi\,\delta\phi,
 \quad \omega=\delta_1\theta(\delta_2)-\delta_2\theta(\delta_1).
\]

Consider the commuting nonlinear coordinate family
`g_yy=exp(2 n)`, `g_yμ=0`, `g_μν=rho² exp(2 psi) gamma_μν`,
`phi=phi_background+X`. Here `a=rho'/rho`, `p=phi_background'`.
All displayed potentials below have `sqrt(-gamma)` factored out, and
`D^mu` raises with the fixed unit de Sitter metric.

The radial connection traces are

\[
g^{BC}\delta\Gamma^y_{BC}
=e^{-2n}[\delta n_y-4\delta\psi_y
-8(a+\psi_y)(\delta\psi-\delta n)],
\quad
g^{yB}\delta\Gamma^C_{BC}=e^{-2n}(\delta n_y+4\delta\psi_y).
\]

For a tangential index, one can instead use
`theta_g^mu=sqrt(-g)(nabla_B h^{mu B}-nabla^mu h)/2`.
The term `Gamma^mu_yy=-exp(2n)g^{mu nu}n_nu` supplies the
`(delta n-delta psi)D^mu n` contribution. Together these give

\[
\theta^y=-\rho^4e^{-n+4\psi}
\{4\delta\psi_y+4(a+\psi_y)(\delta\psi-\delta n)
+(p+X_y)\delta X\},
\]
\[
\theta^\mu=-\rho^2e^{n+2\psi}
\{D^\mu\delta n+3D^\mu\delta\psi
+(\delta n-\delta\psi)D^\mu n+D^\mu X\delta X\}.
\]

Restrict this family to `n=-2 psi`. Its tangent vectors at the background
are exactly the generic longitudinal perturbations used in the prior note.
The antisymmetric current is independent of this choice of commuting
second-order extension; no second-order field equation is being assumed.
The restricted potentials are

\[
\theta^y=-\rho^4e^{6\psi}
[4\delta\psi_y+12(a+\psi_y)\delta\psi+(p+X_y)\delta X],
\quad
\theta^\mu=-\rho^2[D^\mu\delta\psi
+6\delta\psi D^\mu\psi+D^\mu X\delta X].
\]

Taking their field-space exterior derivative at the background gives

\[
j^y/\sqrt{-\gamma}=\rho^4
[X_1X_{2,y}-X_2X_{1,y}
-12(\Psi_1\Psi_{2,y}-\Psi_2\Psi_{1,y})
-6p(\Psi_1X_2-\Psi_2X_1)],
\]
\[
j^\mu/\sqrt{-\gamma}=\rho^2
[X_1D^\mu X_2-X_2D^\mu X_1
+6(\Psi_1D^\mu\Psi_2-\Psi_2D^\mu\Psi_1)].
\]

Using only the momentum constraint
`Psi_y=-2 a Psi-p X/3` in the radial expression reduces its last two
terms to `-2p(Psi1 X2-Psi2 X1)`. This independently reproduces the
parallel derivation's radial and tangential currents, including their
relative normalization and the metric-scalar cross term.

## Separated components, packets and the shell

For two separated solutions with eigenvalues `Lambda_i=mu_i²+4`, the
first-order radial equations yield, with `P_i=rho² Psi_i`,

\[
 j^y/\sqrt{-\gamma}=9\alpha
 [\Lambda_1P_1P_{2,y}-\Lambda_2P_2P_{1,y}],
 \qquad \alpha=(\rho^2p^2)^{-1}.
\]

This reduces to `9 Lambda` times the ordinary Green current only for
equal eigenvalues. The full mixed-eigenvalue current obeys
`partial_y j^y + partial_mu j^mu=0` on the bulk equations. A packet must
retain mixed-eigenvalue terms or use the full local fields directly.

The even shell condition `alpha P_y=Lambda P/beta` cancels the current
for unequal as well as equal eigenvalues. For the odd condition
`alpha P_y=-c_o P`, doubling the two sides gives
`j^y(T)/sqrt(-gamma)=18 c_o(Lambda2-Lambda1)P1 P2`.
This is the tangential divergence of
`j_shell^mu/sqrt(-gamma)=18 c_o(P1 D^mu P2-P2 D^mu P1)`.
With `j^y` oriented toward increasing `y`, the integrated bulk plus shell
current is therefore conserved except for the lower-cone boundary term.
This check uses the shell current supplied by the prior moving-shell
action calculation; it is not a replacement derivation of that action.

## Conditional null-cone limit

Use global de Sitter coordinates, `gamma=-dt²+cosh²(t)dOmega_3²`, and
`u=y exp(-t)`, `v=y exp(t)`. The exterior patch has `u,v>0`,
`y²=uv`; ambient flat coordinates at leading order are
`T=(v-u)/2`, `r=(v+u)/2`. On a compact cone segment `0<v0<=v<=v1`,
assume actual on-shell perturbations have uniform differentiable expansions
`P_i=y⁴ A_i(v,Omega)+O(y⁶)` and regular background
`rho=y+O(y³)`, `p=kappa y+O(y³)`, `kappa!=0`.
The derivative remainder assumptions are essential.

Then `Psi_i=y² A_i+O(y⁴)` and
`X_i -> F_i=-3(4 A_i+v partial_v A_i)/kappa`.
The metric perturbation has the manifestly regular leading expression

\[
 h=2uv A\,g-\tfrac32 A(v\,du+u\,dv)^2+O(uv)
\]

in its tensor components with the corresponding differentiated remainders;
the exact leading identity is `h=2 Psi g-6 Psi dy²`.
This is a compatibility check near `u=0,v>0`, not a global gauge theorem
or an existence construction. It says nothing about the vertex `u=v=0`.

For per-unit-solid-angle densities and `W=F1 partial_v F2-F2 partial_v F1`,
the leading expressions are
`j^y=(u+v)^3 v W/8`, `j^t=-(u+v)^3 v W/(8 sqrt(uv))`.
A vector density transforms with the Jacobian, so

\[
 J^u={ (u/y)j^y-u j^t\over2y}
 \longrightarrow {v^3\over8}
 (F_1\partial_vF_2-F_2\partial_vF_1).
\]

The per-copy result doubles for parity sectors. The current is oriented
toward increasing `u`; outward flux at the lower boundary has the opposite
sign. This is an antisymmetric symplectic flux, not a positive energy flux.
Dropping `j^t` in the transformation loses a factor of two. Integrating
`j^y dt` along `y=epsilon` with `v` fixed produces the same limiting
measure because `dt=dv/v`.

The finite-time endpoint limit and this characteristic limit have different
integration regions. Neither the limit-point property nor vanishing Green
forms at fixed finite `t` proves the displayed null flux vanishes.
Existence of exact registered solutions with these data, their shell
matching, allowed gauge charges and complete physical state space remain
separate unresolved questions.

## Algebra verification

`check_current_independently.py` derives the current as the antisymmetric
variation of the nonlinear-family potential, checks the constraint
reduction, mixed-eigenvalue conservation, both shell balances and the
characteristic density transformation. It deliberately changes six terms
or assumptions to verify that those changes fail the identities. Normal
and optimized Python runs are recorded separately. The counts describe
algebra checks, not independent physical theorems.
