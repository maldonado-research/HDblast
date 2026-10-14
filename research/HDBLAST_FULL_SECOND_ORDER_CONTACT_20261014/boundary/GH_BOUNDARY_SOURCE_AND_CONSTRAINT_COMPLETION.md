# Relative wave gauge: exact shell sources and subsidiary completion

14 October 2026. This supplements the Eulerian scalar harmonic boundary note.
It supplies exact five-dimensional boundary algebra for a possible finite-mode
coefficient construction; the bulk and characteristic existence proof must be
read separately. No vacuum theorem is silently transferred to the matter shell.

## 1. Exact boundary gauge

Use the reference background connection \(\widehat\Gamma=\Gamma[g_0]\), and
\[
C^A=g^{BC}(\Gamma[g]^A{}_{BC}-\widehat\Gamma^A{}_{BC}).
\]
The reduced trace-reversed Einstein equation is
\(E_{AB}-\nabla_{(A}C_{B)}=0\), with
\(E_{AB}=R_{AB}-\phi_A\phi_B-2U(\phi)g_{AB}/3\), and scalar equation
\(S=\Box_g\phi-U'(\phi)=0\). Parentheses include the factor one half.
The shell is fixed at \(y_b\). Impose the four coordinate Dirichlet rows
\(g_{ya}|_b=0\), the physical copy-even junctions, and only the normal
harmonic condition \(C_y|_b=0\). Do not prescribe all five \(C_A\) as extra
independent physical boundary data.

At the boundary write \(g_{yy}=\mathcal N^2\), \(g_{ab}=h_{ab}\). The exact
physical and normal relative-gauge laws are
\[
h_{ab,y}=\frac{\sigma(\phi)}3\mathcal N h_{ab},\qquad
\phi_y=-\frac{\sigma'(\phi)}2\mathcal N,
\]
\[
\boxed{\mathcal N_y=\frac{2\sigma(\phi)}3\mathcal N^2
-a_b\mathcal N^3\operatorname{tr}_{h}(h_0).}\tag{1}
\]
Indeed \(\widehat\Gamma^y_{ab}=-a_bh_{0,ab}\) and
\[
C^y=\mathcal N_y/\mathcal N^3-2\sigma/(3\mathcal N)
+a_b\operatorname{tr}_h(h_0).
\]
On the background \(\mathcal N=1\), \(h=h_0\), \(\sigma_0=6a_b\), (1)
gives \(\mathcal N_y=0\), as required in the original proper-normal chart.
Absolute \(\Gamma^y=0\) would not give this condition; relative gauge matters.

## 2. Complete factorial second sources

Expand \(\mathcal N=1+\varepsilon n_1+\varepsilon^2n_2/2\),
\(h=h_0+\varepsilon h_1+\varepsilon^2h_2/2\), and
\(\phi=\phi_0+\varepsilon\psi_1+\varepsilon^2\psi_2/2\).
Here the first fields must already be the full fields in the relative wave
gauge; they generally differ from the simple longitudinal representative.
Let \(A_1=h_0^{-1}h_1\), \(\tau_1=\operatorname{tr}A_1\), and
\(\tau_2=\operatorname{tr}(h_0^{-1}h_2)\). Exact Taylor differentiation gives
\[
h_{2,ab,y}-2a_bh_{2,ab}
-\left(\frac{\sigma'}3\psi_2+2a_bn_2\right)h_{0,ab}=R_{h,ab},
\]
\[
R_{h,ab}=\left(\frac{2\sigma'}3\psi_1+4a_bn_1\right)h_{1,ab}
+\left(\frac{\sigma''}3\psi_1^2+\frac{2\sigma'}3\psi_1n_1\right)h_{0,ab},
\tag{2}
\]
\[
\psi_{2,y}+\frac{\sigma''}2\psi_2+\frac{\sigma'}2n_2=R_\phi,
\qquad R_\phi=-2\psi_1^2-\sigma''\psi_1n_1,
\tag{3}
\]
\[
n_{2,y}+4a_bn_2-a_b\tau_2-\frac{2\sigma'}3\psi_2=R_n,
\]
\[
R_n=\frac{2\sigma''}3\psi_1^2+\frac{8\sigma'}3\psi_1n_1
-16a_bn_1^2-2a_b\operatorname{tr}(A_1^2)+6a_bn_1\tau_1.
\tag{4}
\]
The scalar coefficient \(-2\psi_1^2\) uses the actual registered
\(\sigma'''=4\). There is one shell tension. The inverse-metric trace term
\(\operatorname{tr}(A_1^2)\) includes all tangential metric components.

For direct metric wave unknowns \(k=h_{2,yy}\), remember
\(k=2n_2+2n_1^2\). Thus (2)–(4) are equivalently
\[
h_{2,ab,y}-2a_bh_{2,ab}
-\left(\frac{\sigma'}3\psi_2+a_bk\right)h_{0,ab}
=R_{h,ab}-2a_bn_1^2h_{0,ab},
\]
\[
\psi_{2,y}+\frac{\sigma''}2\psi_2+\frac{\sigma'}4k
=R_\phi+\frac{\sigma'}2n_1^2,
\]
\[
k_y+4a_bk-2a_b\tau_2-\frac{4\sigma'}3\psi_2
=2R_n+4n_1n_{1,y}+8a_bn_1^2.
\tag{5}
\]
Every non-Dirichlet boundary derivative is a normal derivative of its own
unknown; the couplings in (5) are bounded zeroth-order matrices on a
controlled finite background interval. No second unknown's tangential
derivative occurs in the boundary operator. This is the useful structure
for a finite-angular characteristic trace-ODE construction.

## 3. Subsidiary boundary equations

The exact Codazzi equation and the physical shell laws give
\[
E_{na}=D_bK^b{}_a-D_aK-pD_a\phi
=-3D_a(\sigma/6)-pD_a\phi
=(-\sigma'/2-p)D_a\phi=0.
\tag{6}
\]
This cancellation requires the scalar junction; it is valid for variable
\(\sigma(\phi)\), rather than assuming constant umbilicity. If the reduced
equation holds and \(C_n=0\) at the shell, (6) gives
\[
(\nabla_nC)_a-K_a{}^bC_b=0.
\]
In boundary-orthogonal coordinates this is
\[
\boxed{\partial_yC_a-2\mathcal N K_a{}^bC_b=0.}\tag{7}
\]
It is a Robin equation for covector coordinate components, not literal
coordinate Neumann. With \(C_0=C_1=0\), its factorial order-two form is
\(C_{2,y}=0\) for the normal component and
\(\partial_yC_{2,a}-2a_bC_{2,a}=0\) for tangent components. Here
\(C_{2,y}\) denotes the component with index \(y\), rather than a radial
derivative of a normal component.

The trace-reversed Bianchi identity in five dimensions is
\[
\nabla^AE_{AB}-\tfrac12\nabla_B\operatorname{tr}_gE
=-S\,\partial_B\phi.
\tag{8}
\]
Consequently the reduced Einstein and scalar equations imply
\[
(\Box_g C)_B+R_B{}^AC_A=0.
\tag{9}
\]
If \(C_0=C_1=0\), Taylor order two has **no independent quadratic forcing**:
\((\Box_{g_0}+\operatorname{Ric}_{g_0})C_2=0\), with the homogeneous mixed
Dirichlet/Robin boundary conditions above. Initial vanishing of \(C_2\) must
be supplied by the initial gauge jets; initial vanishing of its time derivative
then follows from physical initial constraints and the reduced equations.
These initial conditions must be checked in the actual construction.

## 4. First wave-gauge conversion and mode-finite initial correction

A shell-preserving first generator satisfies \(\xi_y|_b=0\) and
\(\partial_y\xi_a-2a_b\xi_a=0\). This is the correct condition on covariant
components, since
\((\mathcal L_\xi g_0)_{ya}|_b=\partial_y\xi_a-2a_b\xi_a\).
The corresponding tangential contravariant components have ordinary Neumann
data. The generator can be chosen by
\[
(\Box_{g_0}+\operatorname{Ric}_{g_0})\xi=-C_1^{\rm long}
\]
with those boundary data and chosen initial gauge data. The resulting first
fields satisfy the relative wave gauge if this forced vector problem is
constructed. Its known sources have finite angular representation support.
The whole metric-scalar pair, its initial jets, its boundary sources (2)–(4),
and its bulk quadratic source must be transformed; the longitudinal forcing
cannot be retained unchanged.

An SO(4)-equivariant way to ensure the second physical initial correction is
mode-finite is available without assuming the old right inverse is equivariant.
Let \(L\) be the linearized physical constraint operator on the rotational
background, and let \(Lq=S\) be an existing smooth localized correction. Let
\(P\) be the compact SO(4) isotypic projection onto the finite representation
set containing \(S\), so \(PS=S\). Equivariance gives
\[
L(Pq)=P(Lq)=S.
\tag{10}
\]
The projected correction preserves radial support, because rotations preserve
radii. It therefore preserves an open shell neighborhood where the correction
vanishes. One must use the full tensor and covector group representations;
a componentwise scalar angular cutoff need not commute with \(L\).

For the selected single quadrupole the finite set can be specified exactly.
The first angular state is the SO(4) irrep \((j_L,j_R)=(1,1)\), dimension nine.
The equivariant linear wave-gauge correction stays in the same representation.
Every symmetric quadratic source lies in
\[
\operatorname{Sym}^2(1,1)
=(0,0)\oplus(1,1)\oplus(2,2)\oplus(2,0)\oplus(0,2).
\tag{11}
\]
The dimensions are \(1+9+25+5+5=45=9\cdot10/2\). The diagonal terms supply
scalar \(l=0,2,4\), and the last pair supplies both \(l=2\) tensor orientations.
Difference-one vector representations are absent. This applies to the same
fixed angular quadrupole in every first radial/time amplitude; it does not
assert the same selection for unrelated first angular states. The standalone
checker verifies the full SU(2) character identity, not only the dimension sum.
Projecting onto this complete set in (10) makes the initial coefficient fields
finite-angular while preserving their physical constraints and radial support.

The principal checker verifies the finite-dimensional boundary rank and its
subsidiary algebra. Equations (1)–(5) are checked independently by exact Taylor
differentiation in `verify_scalar_boundary.py`. These boundary calculations
remove the naive gauge-rank defect and supply the explicit inhomogeneous
normal laws. A complete finite-angular coefficient existence statement also
requires the companion bulk/gauge construction, its characteristic estimate,
and its compatible initial jets; this note alone does not assert it.
