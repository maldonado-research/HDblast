# Independent source/boundary cross-review

14 October 2026. Reviewed inputs:

- `astra-bulk/FULL_QUADRATIC_TT_BULK_SOURCE.md`;
- `boundary/FORCED_TENSOR_CONTACT_IBVP.md`;
- the independently derived `sol-bulk/FULL_QUADRATIC_TT_SOURCE.md` and its
  ordinary/optimized exact-check receipts;
- the frozen regular-shell-contact manuscript, Section 5, and the published
  material second-order junction supplement, equations (2), (7), and (9).

**Judgment:** the full quadratic source, linear tensor operator, tortoise
transformation, canonical force, and even straight-shell boundary law are
consistent. There is no missing sign, factorial, warp factor, or tensor
curvature shift in those formulas. The claimed contact advance must remain
the isolated forced tensor coefficient problem in the specified first-order
gauge. The note correctly distinguishes this from the full coupled
second-order contact problem.

The pure copy-odd first-order case requires a separate graph conversion:
its quadratic bulk source is even, while its Eulerian even tensor Neumann
law can be forced. This review accepts that parity and conversion caveat,
without accepting an exact odd-graph coefficient. The published advance
selects the copy-even first-order straight-shell family.

## 1. Common convention and source sign

Both notes use

\[
g_\epsilon=g_0+\epsilon h_1+\tfrac12\epsilon^2h_2,\qquad
\phi_\epsilon=\phi_0+\epsilon\chi_1+\tfrac12\epsilon^2\chi_2,
\]
\[
\mathscr E_{AB}=R_{AB}-\phi_A\phi_B-\tfrac23U(\phi)g_{AB},\qquad
D\mathscr E[h_2,\chi_2]=-D^2\mathscr E[(h_1,\chi_1)]^2.
\]

Set \(c(t)=\cosh t\), \(a=\rho_y/\rho\), and write

\[
\Psi=u(y,t)Y,\quad\chi_1=w(y,t)Y,\quad Y=n_1^2-n_2^2,
\quad C=16u^2+\tfrac83w^2.
\]

The scalar amplitude is denoted \(w\) here to avoid confusing the bulk
source note's scalar \(v\) with the boundary note's canonical variable
\(v\). With the shared Weyl normalization,

\[
(D^2\mathscr E)^{\rm TT}_{ij}=C\mathcal T_{ij},\qquad
(D\mathscr E[h_2,\chi_2])^{\rm TT}_{ij}=-C\mathcal T_{ij}.
\]

The independent sol-bulk derivation confirms the entire source:
\(D^2R^{\rm TT}=-12\Pi_{\rm TT}(d\Psi\otimes d\Psi)\), the scalar
term is \(-2\Pi_{\rm TT}(d\chi_1\otimes d\chi_1)\), and
\(\Pi_{\rm TT}(dY\otimes dY)=-4\mathcal T/3\). The potential term
has zero spatial TT projection. The Hessian has twice the coefficient
of \(\epsilon^2\), which is already included in \(C\).

For
\((h_2)^{\rm TT}_{ij}=\rho^2c^2H\mathcal T_{ij}\), both notes give

\[
(D\mathscr E[h_2,\chi_2])^{\rm TT}_{ij}
       =-\tfrac12\rho^2c^2\mathcal L_T H\mathcal T_{ij},
\]
\[
\mathcal L_T=\partial_y^2+4a\partial_y
-\rho^{-2}(\partial_t^2+3\tanh t\,\partial_t+8c^{-2}).
\]

Thus the signs and factors imply

\[
\boxed{\mathcal L_T H=J=\frac{2C}{\rho^2c^2}
=\frac{16}{3\rho^2c^2}(6u^2+w^2).}
\tag{1}
\]

Equivalently, the convention with a positive temporal principal part is

\[
H_{tt}+3\tanh t\,H_t-\rho^2(H_{yy}+4aH_y)+8c^{-2}H
                         =-2C c^{-2}.
\tag{2}
\]

There is no disagreement between the positive source in (1) and the
negative source in (2): the entire differential equation was multiplied
by \(-\rho^2\). A coefficient-of-\(\epsilon^2\) amplitude would be
\(H/2\), with half these forces.

The rough tensor eigenvalue on \(S^3\) is \(-6\). The de Sitter
curvature shift converts it to \(-8\) in this operator; using 6 in
place of 8 would be incorrect. The spatial TT harmonic extends to a
four-dimensional covariant TT tensor \(q_{ij}=c^2H\mathcal T_{ij}\),
\(q_{0a}=0\), because its covariant temporal divergence is its spatial
trace and its spatial divergence is its sphere divergence. This validates
the boundary derivation's use of the four-dimensional tensor Ricci formula.

## 2. Independent tortoise and canonical derivation

Let \(dx=dy/\rho\), \(x_b=0\). Then

\[
H_y=\rho^{-1}H_x,\qquad
\rho^2(H_{yy}+4aH_y)=H_{xx}+3\rho_yH_x.
\]

Define

\[
S=\tfrac32\partial_x\log\rho=\tfrac32\rho_y,\quad
K=\tfrac32\tanh t,\quad
\mu=\rho^{3/2}c^{3/2},\quad V=\mu H.
\]

The two conjugations are

\[
H_{xx}+2SH_x=\mu^{-1}[V_{xx}-(S^2+S_x)V],
\]
\[
H_{tt}+2KH_t=\mu^{-1}[V_{tt}-(K^2+K_t)V],
\quad K^2+K_t=\tfrac94-\tfrac34c^{-2}.
\]

Applying them to (1) gives

\[
\boxed{V_{tt}-V_{xx}+W_2V=F,\quad
W_2=S^2+S_x-\tfrac94+\tfrac{35}{4}c^{-2},}
\]
\[
\boxed{F=-\rho^{7/2}c^{3/2}J
       =-2\rho^{3/2}c^{-1/2}C
       =-\rho^{3/2}c^{-1/2}(32u^2+\tfrac{16}{3}w^2).}
\tag{3}
\]

These powers and signs agree exactly with the boundary note. The static
potential is \(V_+=S^2+S_x\), not its partner \(S^2-S_x\).

For general Neumann data \(H_y|_b=g(t)\) and odd Dirichlet data
\(H_o|_b=d(t)\), the canonical laws follow directly:

\[
\boxed{V_x-S_bV=\rho_b^{5/2}c^{3/2}g,
\qquad V_o|_b=\rho_b^{3/2}c^{3/2}d.}
\tag{4}
\]

The outward coordinate derivative is \(+\partial_x\); after changing
to inward depth \(z=-x\), the Robin coefficient is \(\alpha=-S_b\).
Its negative boundary-energy term is the actual Robin term. The boundary
note's trace-inequality absorption and bulk shift give a valid finite-time
auxiliary estimate without changing that coefficient.

For comparison with the bulk note's ambient equation, multiplying (1)
by \(\Omega^2\), using \(\rho=r\Omega\) and \(R=rc\), gives the
force \(2C/R^2=16(6u^2+w^2)/(3R^2)\). This also agrees with its
equation (9).

## 3. Copy parity and the homogeneous even-shell law

On each interior the proper-radial coordinate increases outward. The
TT matching and Israel laws in the even longitudinal family with
\(\zeta_1=0\) are

\[
H_+=H_-,\qquad H_{+,y}+H_{-,y}=0.
\]

Thus \(H_o|_b=0\), \(H_{e,y}|_b=0\). The background contribution
\(2a\rho^2c^2H\) from the sum of extrinsic curvatures cancels the
one common tension term \(\sigma_0\rho^2c^2H/3\). This is a law
for \(H_y\), not for \(\partial_y(\rho^2c^2H)\).

For the straight-shell even first-order representative, the first
material induced metric is conformal, the tension Hessian is angular
trace, and the quadratic lapse/radial-metric terms in the extrinsic
curvature are angular trace. The first Gaussian straightening generator
has \(V^y=O(y-y_b)\) and \(V^a=O((y-y_b)^2)\). Its quadratic angular
TT correction has zero boundary value and first normal derivative.
The homogeneous Eulerian laws are therefore justified in this family.

Copy parity at order one does not prescribe arbitrary order-two Cauchy
data. To make this precise, write \(u_\pm=u_e\pm u_o\) and
\(w_\pm=w_e\pm w_o\). The exact source decomposition is

\[
C_e=16(u_e^2+u_o^2)+\tfrac83(w_e^2+w_o^2),\qquad
C_o=32u_eu_o+\tfrac{16}{3}w_ew_o.
\tag{5}
\]

Pure even and pure odd first-order fields both give \(C_o=0\).
Mixed first-order parity can force the odd tensor channel. Even when
\(C_o=0\), an independently supplied homogeneous odd second-order
initial tensor mode need not vanish. Vanishing follows only with even
second-order initial data and homogeneous odd boundary data, or from an
imposed parameter-reversal family
\(g_-(\epsilon)=g_+(-\epsilon)\), which has an even second derivative.
Assigning an odd second-order field merely because the first derivative
is odd is incorrect.

## 4. Odd first-order graphs require the material conversion

The published odd packet uses opposite normal graphs, with
\(\zeta_o=-\Psi_{o,b}/a_b\), rather than the even family's
\(\zeta_1=0\). Its first material induced metric and scalar trace vanish,
but second-order graph and normal-coordinate corrections can carry an
angular TT contribution. Products of two copy-odd first-order fields are
copy-even, including the material-pullback correction

\[
B=\Pi_{\rm TT}\frac{2\mathcal L_Vh_1+\mathcal L_V^2g_0}{\rho^2c^2}.
\]

Thus pure odd first-order data do not justify an odd second-order
Dirichlet assignment. They also do not justify setting the Eulerian
even Neumann force to zero. The boundary note's conversion
\(H_o^w|_b=-B_o|_b\), \(H_{e,y}^w|_b=-B_{e,y}|_b\) is required.
For a pure odd first-order state the quadratic \(B\) is even, so matching
has no odd quadratic forcing, while the even normal force can be nonzero.
An arbitrary first-order wave gauge requires its actual embedding and
Gaussian-collar generator, including the prescribed tangential identification.

No exact odd-graph force coefficient is accepted in this review. A separate
local calculation is retained in `ODD_GRAPH_SOURCE_CANDIDATE.md` pending
independent coefficient review. The accepted contact result selects the
copy-even first-order straight-shell family and does not depend on that
candidate.

## 5. Precisely accepted scope

The finite-slab forced tensor argument is valid on a controlled
positive-cone collar with actual smooth background coefficients, a known
smooth first-order packet, prescribed compatible second-order tensor
initial data, and the appropriate TT boundary forces. Finite de Sitter
time and compact tortoise support keep the construction away from the
cone endpoint. The characteristic boundary count, finite propagation,
Robin energy cancellation, and smooth boundary lifts concern this tensor
mode only.

Contact in this result means evolution beyond the first-order packet's
causal shell-arrival time. It does not prove a nonzero tensor trace at that
arrival front or a radiation observable. Constraint-corrected tensor
initial data are inputs, rather than freely selectable replacements for
the full gravitational initial data.

No conclusion is made here about coupled scalar/vector/gauge incoming
boundary fields, propagation of the complete constraints through the
material shell, existence of an exact nonlinear family, ambient extension
through the cone, or global stability. Those remain outside the accepted
tensor coefficient result.

## 6. Verification evidence and review method

The existing sol-bulk checker passed normally and under `python -O`
before this review; it was not rerun. It computes direct five-dimensional
Ricci and global sphere identities, independently confirming the full
source. Astra's operator and sphere check receipts and the boundary
operator/conjugation receipt were inspected. Their recorded scope is
algebraic, as it should be.

This review independently derived equations (1)--(5), including the source
factorial, both coordinate-conjugation powers, the one-tension warp
cancellation, and parity expansion. The odd-graph coefficient candidate is
separate and pending independent review. No frozen or tracked manuscript
was edited by this review.
