# A forced tensor sector through shell contact

14 October 2026, America/Los_Angeles. Analytic continuation of the frozen
Einstein–scalar checkpoints. This note proves a finite-mode, linear forced
tensor initial-boundary result. It does not solve the coupled second-order
scalar, vector, gauge or constraint boundary problem, and does not prove an
exact nonlinear evolution through contact.

## 1. Fixed background and angular tensor convention

On the actual positive-cone background use

\[
 g_0=dy^2+\rho(y)^2\gamma,\qquad
 \gamma=-dt^2+\cosh^2t\,\kappa_{ij}dn^idn^j,
 \quad \phi_0=\phi_0(y),\quad a=\rho'/\rho.
\]

Here \(\kappa\) is the unit round metric on \(S^3\); \(y\) increases
outward on each interior, the shell is \(y=y_b\), and
\(a_b=\sigma_0/6\). Fix a smooth tensor harmonic \(T_{ij}\) with

\[
 D^iT_{ij}=0,\qquad \kappa^{ij}T_{ij}=0,\qquad
 -D^kD_kT_{ij}=k_T^2T_{ij},\quad
 k_T^2=l(l+2)-2,
 \quad \lambda=k_T^2+2=l(l+2),\quad l\ge2.
\]

The second metric derivative, using the convention
\(g(\varepsilon)=g_0+\varepsilon g_1+\varepsilon^2g_2/2+\cdots\), has
tensor component

\[
 (g_2)_{ij}^{\rm TT}=\rho^2\cosh^2t\,H(y,t)T_{ij}.
 \tag{1}
\]

All tensor projections in this note are the orthogonal spatial \(S^3\)
TT projection at fixed \((y,t)\). This is a compact, positive Riemannian
projection, not a putative elliptic TT projection on the whole Lorentzian
de Sitter spacetime. Integration by parts on \(S^3\) proves that it
annihilates angular trace and every angular symmetric derivative. Hence
linear scalar and vector metric fields do not mix into (1).

## 2. Derivation of the tensor operator and its source sign

Let \(q_{ab}\) denote a four-dimensional covariant transverse trace-free
metric perturbation of \(\gamma\). The radial Gauss formula is

\[
 R^{(5)}_{ab}=R[h]_{ab}-K'_{ab}+2K_{ac}K^c{}_b-KK_{ab}.
\]

For \(h_{ab}=\rho^2(\gamma_{ab}+\epsilon q_{ab})\), linearization gives

\[
 \delta R^{(5)}_{ab}
 =\delta R[\gamma]_{ab}
 -\rho^2\{(a'+4a^2)q_{ab}+2a q'_{ab}+q''_{ab}/2\}.
\]

The actual scalar background obeys
\(3-\rho^2(a'+4a^2)=2\rho^2U/3\). Its tensor Einstein residual therefore
equals its trace-reversed residual
\(\delta R^{(5)}_{ab}-2U\delta g_{ab}/3\). On unit de Sitter,
\(\delta R[\gamma]_{ab}-3q_{ab}=-(\Box_\gamma-2)q_{ab}/2\). Thus

\[
 \delta E^{\rm TT}_{ab}
 =-\tfrac12\{(\Box_\gamma-2)q_{ab}
       +\rho^2(q''_{ab}+4a q'_{ab})\}.
 \tag{2}
\]

For \(q_{ij}=\cosh^2t\,H T_{ij}\), \(q_{0a}=0\), the spacetime
divergence is zero: its temporal component is the spatial trace and its
spatial components are the angular divergence. More explicitly, if
\(c=\cosh t\) and \(h=\dot c/c\), then

\[
 \nabla_tq_{ij}=c^2H_tT_{ij},\qquad
 (\Box_\gamma q)_{ij}
 =c^2[-H_{tt}-3hH_t+(2h^2-k_T^2/c^2)H]T_{ij}.
\]

The \(2h^2\) term comes from differentiating the otherwise zero mixed
components \(q_{0i}\); omitting it gives an incorrect curvature shift.
Since \(h^2-1=-c^{-2}\),

\[
 ((\Box_\gamma-2)q)_{ij}
 =c^2[-H_{tt}-3\tanh(t)H_t-\lambda\operatorname{sech}^2(t)H]T_{ij}.
\]

Consequently

\[
 \boxed{\delta E_{ij}^{\rm TT}
 =-\tfrac12\rho^2c^2\,\mathcal L_T H\,T_{ij}},\qquad
 \mathcal L_T H=H_{yy}+4aH_y
 -\rho^{-2}[H_{tt}+3\tanh(t)H_t+\lambda\operatorname{sech}^2(t)H].
 \tag{3}
\]

If the *complete* quadratic Einstein residual in this tensor harmonic is
\(D^2E[g_1,\phi_1;g_1,\phi_1]^{\rm TT}_{ij}=C(y,t)T_{ij}\), its second
equation is \(\delta E[g_2,\phi_2]^{\rm TT}=-C T\). Therefore

\[
 \boxed{\mathcal L_T H=J,\qquad J=\frac{2C}{\rho^2c^2}}.
 \tag{4}
\]

This coefficient uses the derivative convention in (1). A Taylor
coefficient without its factorial has a different source normalization.
The matter-only Hessian cannot be substituted for the complete \(C\).
The source proved in the separate complete bulk-source calculation for the inherited scalar
reconstruction is \(C=16u^2+8v^2/3\), with
\(\Psi=u(y,t)Y(n)\), \(\chi=v(y,t)Y(n)\). Equation (3), independently of
that separate source calculation, fixes its sign and prefactor:

\[
 H_{tt}+3\tanh(t)H_t-\rho^2(H_{yy}+4aH_y)
 +8\operatorname{sech}^2(t)H
 =-(32u^2+16v^2/3)\operatorname{sech}^2(t)
\]

for the \(l=2\) tensor. The displayed \(C\) uses that separate complete
bulk-source proof and its checks; the operator normalization and boundary
construction are derived independently here.

## 3. Shell law in the inherited first-order longitudinal gauge

First take the reflection-even scalar-master family. Its inherited
first-order representative is

\[
 (g_1)_{yy}=-4\Psi,\quad (g_1)_{ya}=0,\quad
 (g_1)_{ab}=2\rho^2\Psi\gamma_{ab},\qquad \zeta_1=0.
 \tag{5}
\]

The induced first metric is conformal to the background, so the complete
quadratic material Israel source
\(2\sigma'_0\psi_1h_1/3+\sigma''_0\psi_1^2h_0/3\) has zero angular TT
projection. In a fixed-shell radial ADM chart, all products of the
first-order lapse with the first-order radial metric derivative in
\(K_{ab}=(\partial_yh_{ab}-D_a\beta_b-D_b\beta_a)/(2N)\) are also
angular trace. An arbitrary second-order shift \((g_2)_{yi}\) contributes
only an angular symmetric derivative. An arbitrary second-order normal
embedding contributes a background normal derivative and an angular
Hessian; tangential reidentification contributes a symmetric derivative.
Their angular TT projections are zero. Thus the second TT matching and
Israel equations in this **specified first-order gauge** are

\[
 H_+(y_b,t)=H_-(y_b,t),\qquad
 \partial_yH_+(y_b,t)+\partial_yH_-(y_b,t)=0.
 \tag{6}
\]

In particular, with \(H_e=(H_++H_-)/2\), \(H_o=(H_+-H_-)/2\),

\[
 \boxed{\partial_yH_e|_b=0,\qquad H_o|_b=0.}
 \tag{7}
\]

The background warp term cancels: each
\(K_{2,ij}^{\rm TT}=\rho^2c^2(aH+H_y/2)T_{ij}\), while the common
Israel tension subtracts \(2a\rho^2c^2HT_{ij}\). A boundary condition
\(\partial_y(g_2)_{ij}^{\rm TT}=0\) would lose that cancellation and is
incorrect. No exact nonlinear reflection symmetry was needed to obtain
the two conditions in (6).

For identical second-order sources and even initial tensor data, uniqueness
below makes \(H_o=0\). If instead first-order fields are copy-odd, their
quadratic bulk source is copy-even. A family obeying
\(g_-(\varepsilon)=g_+(-\varepsilon)\) has an **even second derivative**.
It must not be assigned an odd second-order Dirichlet law merely because
its first derivative is odd. Such first-order odd states generally have
\(\zeta_1\ne0\), however, so quadratic graph terms can force the Eulerian
even boundary law. The homogeneous Eulerian Neumann conclusion in (7)
uses (5); generic odd graph states need Section 4.

## 4. Conversion to an arbitrary bulk wave gauge

In a background radial collar let \(F_\varepsilon\) be the material
straightening map, with first generator \(V\), second generator derivative
\(W\), and first wave-gauge field \(g_1^w\). The material tensor derivative
is

\[
 H^{\rm mat}=H^w+B,\qquad
 B=\Pi_{\rm TT}\!\left[
 \frac{2\mathcal L_Vg_1^w+\mathcal L_V^2g_0}{\rho^2c^2}
 \right].
 \tag{8}
\]

The expression means harmonic amplitude projection, with the radial and
time factors in (1). The \(\mathcal L_Wg_0\) term has zero angular TT
projection: its normal and temporal pieces are angular trace, and its
angular piece is a symmetric derivative. Its normal derivative has the
same property. Thus \(B\) is fully determined by first-order gauge and
embedding data, rather than by unknown second-order gauge fields.

The first-order Gaussian collar generator is found from

\[
 2\partial_yV^y=-(g_1^w)_{yy},\qquad
 \rho^2\gamma_{ab}\partial_yV^b+\partial_aV^y=-(g_1^w)_{ya},
 \quad V^y|_b=\zeta_1,
\]

with its prescribed tangential shell identification as the remaining
integration data. These equations retain displacement and normal
normalization. For the pure-trace material source considered above, (6)
becomes

\[
 H_o^w|_b=-B_o|_b,\qquad
 \partial_yH_e^w|_b=-\partial_yB_e|_b.
 \tag{9}
\]

If a more general first-order material state has nonzero TT Israel source,
add its prescribed amplitude to the right side of the derivative condition.
Similarly, the wave-gauge bulk forcing is shifted by
\(-\mathcal L_T B\) when converting a material equation to \(H^w\).
That last formula assumes a material equation with its **actual material
source already computed**. The coefficient \(C=16u^2+8v^2/3\) in Section 2
is the source in the specified longitudinal first gauge; it is neither a
universal material source nor an unchanged arbitrary wave-gauge source.
Changing first-order gauge requires transforming or recalculating the
complete quadratic bulk forcing as well as the boundary jets of \(B\).
No Gaussian conversion of the interior fields is needed for the direct
longitudinal construction in Sections 3 and 5–7.
Explicitly, if \(B_{\rm long}\) straightens the longitudinal fields and
\(B_w\) straightens the wave-gauge fields using the same material chart,
then
\(J_{\rm mat}=J_{\rm long}+\mathcal L_TB_{\rm long}\) and
\(J_w=J_{\rm mat}-\mathcal L_TB_w\). These are transformed complete
sources, with \(J_{\rm long}=2C/(\rho^2c^2)\) as in (4).
At linear second order the generalized-harmonic gauge-vector derivative is
an angular symmetric derivative and has zero TT projection; products with
the first gauge residual vanish when that residual is zero. This isolates
the same tensor equation, but supplies no conditions for the other incoming
gauge or constraint fields.

For (5), \(V|_b=0\), and writing \(d=y-y_b\),
\(V^y=2\Psi_b d+O(d^2)\),
\(V^a=-\rho_b^{-2}D^a\Psi_b d^2+O(d^3)\). The angular TT part of both
quadratic terms in (8) is \(O(d^2)\): angular derivatives of \(V^y\)
occur in products, whereas every angular \(V^a\) is already \(O(d^2)\).
Hence \(B|_b=B_y|_b=0\), recovering the direct longitudinal-gauge proof.
An arbitrary wave gauge need not have these zero jets. Prescribing
homogeneous Neumann there without (8) would be an unjustified gauge change.

## 5. Tortoise coordinates and the exact forced collar equation

Let \(r\) be the inherited ambient spacelike radius, so
\(dy=\Omega(r^2)dr\), \(\rho=r\Omega(r^2)\), and

\[
 x=\log(r/b),\qquad dx=dy/\rho,\qquad x_b=0,
 \quad T=r\sinh t,\quad R=r\cosh t.
\]

The cone is \(x=-\infty\), not a freely selectable boundary condition.
The physical shell is fixed at \(x=0\). Set

\[
 S(x)=\frac32\partial_x\log\rho=\frac32\rho',\quad
 V_+(x)=S^2+S_x,\quad
 v(x,t)=\rho(x)^{3/2}\cosh(t)^{3/2}H(x,t).
\]

Direct conjugation of (4) gives

\[
 \boxed{v_{tt}-v_{xx}+W_l(x,t)v=F(x,t)},\qquad
 W_l=V_+(x)-\frac94+(\lambda+\frac34)\operatorname{sech}^2t,
 \quad F=-\rho^{7/2}\cosh(t)^{3/2}J.
 \tag{10}
\]

In particular the full-source normalization from (4) is
\(F=-2\rho^{3/2}\cosh(t)^{-1/2}C\). The inherited static tensor operator
\(V_+\) and its partner are useful consistency checks, but their spectral
positivity is not needed for finite-time well-posedness here.

For a prescribed Neumann forcing \(H_y|_b=g(t)\), and a prescribed odd
trace \(H_o|_b=d(t)\), the canonical shell equations are

\[
 v_x(0,t)-S_bv(0,t)=G(t),\quad
 G=\rho_b^{5/2}\cosh(t)^{3/2}g(t),\qquad
 v_o(0,t)=D(t)=\rho_b^{3/2}\cosh(t)^{3/2}d(t).
 \tag{11}
\]

For the even longitudinal family \(G=D=0\). Notice the negative Robin
potential \(-S_b\) in the even canonical variable; replacing it by an
arbitrary positive boundary term would destroy the exact flux identity.

## 6. Finite-slab construction and boundary count

Fix \(t_*<\infty\). Assume smooth compatible initial data have compact
support in \([-L,-a]\), \(0<a<L<\infty\), and the given smooth source
vanishes below \(x=-L-t\) for \(0\le t\le t_*\). It suffices to know the
actual background on \([-L-t_*,0]\); its smooth coefficients are bounded
there. Extend those coefficients boundedly to the left only as a device
for the characteristic construction. Finite propagation below proves
that the extension is never used on the physical solution's support.
Boundary sources may be included if smooth and compatible. For the
inherited compact collar data, the bulk source and the boundary source
are zero on a neighborhood of the initial shell corner, so all corner
compatibility identities hold there.

For completeness, in the inward variable of the next paragraph define
\(v_j=\partial_t^jv|_{t=0}\). The bulk recursion is
\(v_{j+2}=\partial_z^2v_j-\sum_{k=0}^j\binom jk
\partial_t^kW_l|_0\,v_{j-k}+\partial_t^jF|_0\).
Constant-Robin corner compatibility is
\(-\partial_zv_j|_0+\alpha v_j|_0=G^{(j)}(0)\); Dirichlet compatibility
is \(v_j|_0=D^{(j)}(0)\), for every required \(j\). This explicitly
specifies the smooth compatible class used by the construction.

Use inward depth \(z=-x\ge0\), write \(\alpha=-S_b\), and define
\(p=v_t+v_z\), \(m=v_t-v_z\). Then

\[
 p_t-p_z=F-W_lv,\qquad m_t+m_z=F-W_lv.
 \tag{12}
\]

For the outward Robin law \(-v_z+\alpha v=G\), the shell trace
\(q(t)=v(0,t)\) obeys

\[
 q'=p(0,t)-\alpha q+G,\qquad
 m(0,t)=p(0,t)-2\alpha q+2G.
 \tag{13}
\]

For Dirichlet data \(q=D\), simply
\(m(0,t)=2D'(t)-p(0,t)\). Thus the outgoing \(p\) is determined by the
bulk; the boundary supplies precisely one incoming characteristic \(m\)
in each parity. These are two scalar boundary conditions for the two-copy
tensor mode, not a count for the full Einstein–scalar system.

For example, the outgoing Volterra formula is

\[
 p(t,z)=p_0(z+t)+\int_0^t(F-W_lv)(s,z+t-s)\,ds,
\]

and the incoming formula integrates (12) back to the initial slice for
\(z\ge t\), or to \((t-z,0)\) with (13) for \(z<t\). Finally
\(v(t,z)=v_0(z)+\frac12\int_0^t(p+m)(s,z)ds\), and the trace ODE (13)
is integrated in time. Substitute the outgoing formula into the incoming
formula and into the trace ODE. On a time interval of length \(\delta\),
differences of outgoing values are bounded by
\(\delta\|W_l\|_\infty\|\Delta v\|_\infty\). Incoming differences are
bounded by \(2|\alpha|\|\Delta q\|_\infty\) plus a constant times
\(\delta\|W_l\|_\infty\|\Delta v\|_\infty\); final time integration
therefore bounds the joint map difference by
\(C\delta(\|\Delta v\|_\infty+\|\Delta q\|_\infty)\), with a finite
constant depending only on bounded coefficients. Dirichlet data have a
prescribed trace and the same argument with only \(v\) unknown. Choosing
\(C\delta<1\) gives a unique continuous solution by contraction.

Differentiating these Volterra identities gives all derivatives for smooth
sources and data; across the characteristic leaving the initial corner,
the matching conditions are precisely the usual recursive corner
compatibility identities. Data and forcing zero near that corner give
zero jets of every order. The defects in
\(v_z=(p-m)/2\) and \(v|_0=q\) have zero initial value and zero evolution
by (12)–(13), so the constructed system solves the original equation.
Successive bounded intervals cover any specified finite \(t_*\).

The same construction gives finite propagation:
\(v(t,z)=0\) for \(z>L+t\), when the initial data and the source have
the stated support. The left coefficient extension consequently cannot
affect the physical result. The cone \(z=\infty\) cannot be reached by
this compact support at finite global de Sitter time.

For the even scalar family, the full quadratic forcing is smooth and
supported within the causal support of the known first-order master
packet. Products and derivatives do not enlarge that support. Each smooth
second-order initial tensor profile localized in the same controlled collar
therefore fits this forced theorem. Initial profiles supplied by a chosen
localized constraint correction are fixed data for this theorem, not
additional unconstrained choices for the full gravitational initial state.

## 7. Auxiliary energy, uniqueness and smooth boundary sources

For homogeneous Robin data in inward \(z\), define the **actual**
Robin form, with a constant bulk shift \(K\),

\[
 E_K=\frac12\int_0^\infty
 [v_t^2+v_z^2+(W_l+K)v^2]dz+\frac\alpha2v(0,t)^2.
 \tag{14}
\]

The trace inequality
\(|v(0)|^2\le\epsilon\|v_z\|_2^2+\epsilon^{-1}\|v\|_2^2\)
first absorbs the negative boundary term; then a sufficiently large
\(K\) absorbs the bounded negative bulk potential. Thus \(E_K\) controls
\(\|v_t\|_2^2+\|v_z\|_2^2+\|v\|_2^2\). Integration by parts and the
exact Robin law give

\[
 E_K' =\int_0^\infty Fv_t\,dz
 +\frac12\int_0^\infty (W_l)_t v^2dz
 +K\int_0^\infty vv_tdz.
 \tag{15}
\]

The boundary flux cancels exactly with the derivative of
\(\alpha v(0)^2/2\). There is no uncontrolled boundary trace of \(v_t\)
in the right side. With a time-dependent Robin coefficient, its additional
term is \(\alpha'v(0)^2/2\), controlled by the same trace inequality.
For Dirichlet data omit the boundary term. Cauchy–Schwarz gives
\(E_K'\le CE_K+C\|F\|_2^2\), proving finite-slab uniqueness and
continuous dependence. This is a positive **auxiliary tensor control norm**,
not the nonlinear gravitational Hamiltonian or a stability/ghost statement.

For smooth inhomogeneous Robin data, choose a compact collar function
\(\xi\) with \(\xi(0)=0\), \(\xi_z(0)=-1\). Then
\(w=v-\xi G\) has homogeneous Robin law and forced source
\(F-\xi G''+\xi''G-W_l\xi G\). For Dirichlet use a compact
\(\eta\) with \(\eta(0)=1\), and \(w=v-\eta D\) similarly. These explicit
lifts justify the source-dependent energy estimate for smooth sources;
they require the corresponding two time derivatives, rather than claiming
a boundary \(L^2\) estimate that the basic interior energy cannot support.

## 8. Actual contact and scope of the advance

There is a simple compatible concrete choice for the initial tensor state.
At the time-symmetric initial slice \(t=T=0\), the four-metric is
\(\gamma_4=dy^2+\rho^2\kappa\), its extrinsic curvature and scalar
momentum vanish, and its scalar gradient is radial. An angular tensor
metric datum \(q_{ij}=f(y)T_{ij}\), \(q_{yA}=0\), has zero four-trace
and four-divergence: the radial part of its divergence is proportional to
the angular trace, and the angular part to \(D^iT_{ij}\). The angular
Ricci tensor of \(\gamma_4\) is proportional to \(\kappa\), so
\(\operatorname{Ric}_4\!:\!q=0\). Therefore
\(\delta R=\nabla\nabla q-\Delta\operatorname{tr}q-
\operatorname{Ric}_4\!:\!q=0\), and the metric variation of the scalar
gradient norm is also zero because \(q_{yA}=0\). This metric datum is in
the kernel of the linear Hamiltonian constraint.

Similarly an angular TT momentum density
\(\pi^{ij}=\sqrt{\gamma_4}\,f(y)T^{ij}\), with indices raised by
\(\gamma_4\), has zero covariant divergence and lies in the kernel of
the linear momentum constraint. Smooth compact radial profiles preserve
support. The second-order initial constraints have exactly these same
linear operators acting on their unknown second derivatives, with the
first-data quadratic source already fixed. Consequently subtracting the
TT metric and momentum parts of any second-order constraint correction
leaves those sourced constraints unchanged. One can choose
\(H|_{t=0}=H_t|_{t=0}=0\) without changing the other second-order
constraint fields. The first-order lapse/kinetic products in (5) have
only angular trace, so this zero canonical velocity is the same zero TT
momentum choice. These are free homogeneous tensor directions, not a
claim that arbitrary other gravitational initial data satisfy constraints.
The forced theorem now supplies a concrete zero-TT-data tensor coefficient
for the previously solved sourced initial constraint state.

In these coordinates radial null characteristics are \(t\pm x=\text{constant}\).
An initial support outer edge \(R_0<b\) reaches the shell causally at

\[
 t_c=\log(b/R_0),\qquad
 T_c=b\sinh t_c=\frac{b^2-R_0^2}{2R_0}.
 \tag{16}
\]

Choose any controlled finite interval extending past this \(t_c\) on which
the existing first-order contact solution is known. The forced tensor
construction above then exists smoothly throughout it. It does not infer
that a specific forced tensor trace is nonzero at the causal support edge;
that needs a separate source/data argument. This chart covers the positive
\(s=r^2\) region and the timelike shell, not a full ambient-time moving ball
through the cone or a global spacetime slab. Finite \(t\) bounded support
stays strictly within the actual regular background. One can combine it
with an independently established regular ambient construction if needed,
but no such extension is claimed by this tensor theorem.

The controlled first-order interval can be chosen explicitly from the
inherited conformal contact construction. If that solution is known on
\(|\tau|\le\tau_*\), its coordinate map on \(x\le0\), \(t\ge0\) is

\[
 \tau(x,t)=\arctan(e^{x+t})-\arctan(e^{x-t})
 \le\tau(0,t)=\arctan(\sinh t).
\]

For positive \(t\) the difference increases with \(e^x\le1\), as seen
by differentiating; the derivative has numerator proportional to
\((e^t-e^{-t})(1-e^{2x})\). Thus choosing
\(t_c<t_*<\operatorname{arsinh}(\tan\tau_*)\) keeps the entire tensor
support in the known first-order conformal time slab. A packet whose
outer edge is sufficiently close to \(b\) has arbitrarily small
\(t_c=\log(b/R_0)\), so such a nonempty postcontact interval exists.
This uses the actual inherited contact background and master solution;
it does not assume their extension to a generic long time interval.

The progress is a solved, source-forced **tensor coefficient** with the
correct shell law in a controlled first-order gauge, and an explicit
translation rule for another gauge. It supplies a required component of
the second-order contact problem. A TT solution alone does not establish
the full second-order Einstein–scalar equations, their constraint
propagation or the incoming scalar/vector/gauge boundary conditions.

## 9. Verification and references

The companion symbolic verifier independently checks the radial curvature
cancellation, the de Sitter tensor curvature shift, both conjugations,
source normalization, boundary warp cancellation, correct canonical Robin
sign, characteristic count, energy flux cancellation and the contact-time
identity. Wrong-formula controls deliberately change these signs and
factors. Those checks support the algebra; the Volterra construction and
energy arguments in Sections 6–7 are the analytic proof.

Model and conventions: the frozen [shell-Noether manuscript](https://github.com/maldonado-research/HDblast/blob/main/research/HDBLAST_SHELL_NOETHER_20261013/manuscript/SHELL_NOETHER.md).
Full material junction source and second-order pullback: [material supplement](https://github.com/maldonado-research/HDblast/blob/main/research/HDBLAST_MATERIAL_JUNCTIONS_20261014/MATERIAL_SECOND_ORDER_JUNCTIONS.md).
Actual first-order contact solution and regular background domain: [contact construction](https://github.com/maldonado-research/HDblast/blob/main/research/HDBLAST_REGULAR_SHELL_CONTACT_20261013/existence/FINITE_SLAB_SHELL_CONTACT_EXISTENCE.md).
Inherited tensor spectral reduction and endpoint conventions: [tensor partner note](https://github.com/maldonado-research/HDblast/blob/main/research/HDBLAST_REGULAR_SHELL_SPECTRA_20261013/tensor/TENSOR_PARTNER_GAP.md).
The positive auxiliary control norm does not import a vacuum timelike
Einstein-boundary theorem into this five-dimensional scalar-coupled model.

Related general linear methods include [Ginoux–Murro, arXiv:2007.02544](https://arxiv.org/abs/2007.02544),
Theorem 6.8 and its admissible covariant-Neumann example, and the recent
[Costeri–Dappiaggi–Juárez-Aubry, arXiv:2609.24465v1](https://arxiv.org/abs/2609.24465v1),
Proposition 3.5 and Theorem 3.6 for the scalar Klein–Gordon Robin problem
and causal fundamental solutions. The latter paper takes an inward normal;
its Robin sign must be converted to the outward convention here. These
are neighboring linear wave frameworks; the direct finite-mode proof
above neither transfers their conclusions to the coupled nonlinear
Einstein–scalar shell system nor asserts newness of the wave methods.
