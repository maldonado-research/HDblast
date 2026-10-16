# First relative-harmonic gauge vector: explicit warped and scalar-amplitude system

15 October 2026. Bounded unpublished next-step calculation. This note derives the actual forced gauge-vector operator and boundary rows. It does **not** solve the vector wave equation, prove its initial-boundary well-posedness, certify an evolved gauge vector, or compute transformed second-order forcing. It belongs to future work and is separate from the completed initial-data publication.

## 1. Background, source and sign

Retain

\[
g_0=dy^2+\rho(y)^2\bar g,\qquad
\bar g=-dt^2+\cosh^2t\,\kappa_{S^3},\qquad a=\rho_y/\rho.
\]

The longitudinal perturbation is

\[
h_{yy}=-4\Psi,\qquad h_{ya}=0,\qquad h_{ab}=2\Psi(g_0)_{ab},
\qquad \Psi=u(y,t)Y,\quad Y=n_1^2-n_2^2.
\]

Use \(\Box=\nabla^A\nabla_A\), the inherited Ricci convention, and

\[
C^A=g^{BC}(\Gamma[g]^A_{BC}-\Gamma[g_0]^A_{BC}),\qquad
DC[h]_A=\nabla^Bh_{BA}-\tfrac12\nabla_A\operatorname{tr}_0h.
\]

Here \(\operatorname{tr}_0h=4\Psi\). Direct differentiation gives

\[
C^{\rm long}_{1,y}=-6(\Psi_y+4a\Psi),\qquad
C^{\rm long}_{1,a}=0.
\]

For example \(\nabla^Bh_{By}=-4\Psi_y-24a\Psi\); subtracting half the trace derivative gives the displayed radial source. The tangential divergence equals \(2\partial_a\Psi\), canceling the trace term. Covariance gives

\[
DC[\mathcal L_\xi g_0]=(\Box_0+\operatorname{Ric}_0)\xi.
\]

Thus the desired first conversion \(h^{\rm wave}=h^{\rm long}+\mathcal L_\xi g_0\) requires

\[
P\xi:=(\Box_0+\operatorname{Ric}_0)\xi=-C^{\rm long}_1,
\qquad \xi|_{t=0}=\partial_t\xi|_{t=0}=0.                       \tag{1}
\]

The radial contravariant source is \(+6(\Psi_y+4a\Psi)\), since \((g_0)_{yy}=1\). The tangential source is zero.

## 2. Exact contravariant warped split

Write \(X=\xi^y\), \(V^a=\xi^a\), with barred operations taken in \(\bar g\). The vector rough Laplacian and Ricci contributions combine to

\[
\boxed{(P\xi)^y=X_{yy}+4aX_y+\rho^{-2}\bar\Box X
-2a\,\overline{\operatorname{div}}V-(4a_y+8a^2)X,}             \tag{2}
\]
\[
\boxed{(P\xi)^a=V^a_{yy}+6aV^a_y
+\rho^{-2}(\bar\Box V^a+3V^a)
+2a\rho^{-2}\bar\nabla^aX.}                                 \tag{3}
\]

The operator \(\bar\Box V\) here is the rough vector wave operator, not the scalar wave operator applied separately to components.

For verification, the warped connection is

\[
\Gamma^y_{ab}=-a(g_0)_{ab},\quad \Gamma^a_{yb}=a\delta^a_b,
\quad \Gamma^a_{bc}=\bar\Gamma^a_{bc},
\]

and

\[
\operatorname{Ric}^y{}_y=-4(a_y+a^2),\qquad
\operatorname{Ric}^a{}_b=(3/\rho^2-a_y-4a^2)\delta^a_b.
\]

Before adding Ricci, the radial rough-wave zeroth-order term is \(-4a^2X\). The tangential one is \((a_y+4a^2)V^a\); adding Ricci cancels this entire tangential warp potential, leaving \(3V^a/\rho^2\) in (3). This cancellation prevents an extra spurious radial-potential term in the tangential equations.

## 3. Five shell-preserving boundary rows

Keep the shell at the fixed coordinate worldtube \(y=y_b\). Its infinitesimal preservation requires

\[
\boxed{X|_b=0.}                                                \tag{4}
\]

The cross metric variation is

\[
(\mathcal L_\xi g_0)_{ya}
=\partial_aX+\rho^2\bar g_{ab}\partial_yV^b
=\partial_a\xi_y+(\partial_y-2a)\xi_a.
\]

Since the tangential derivatives of \(X|_b=0\) vanish, preserving \(h^{\rm wave}_{ya}|_b=0\) gives exactly four more rows:

\[
\boxed{\partial_yV^a|_b=0
\quad\Longleftrightarrow\quad
(\partial_y-2a_b)\xi_a|_b=0.}                                 \tag{5}
\]

There are five rows for five vector-wave components: one Dirichlet row and four Neumann rows in these contravariant coordinates. The scalar principal operator is the background wave operator, and the shell is timelike. This establishes the algebraic boundary form and count, not an energy estimate or an existence theorem.

Crucially, (5) does **not** set \(V^a|_b\) to zero. The vector may generate a tangential shell reparametrization, including a change of its time label. The physical junction residuals are natural tensors and scalars on the shell. A diffeomorphism tangent to the shell pulls those residuals back, so their first variations remain zero when the original physical laws hold. The shell's embedding remains at \(y_b\), while its tangential labels and clock parametrization may change. Normal variations and the transformed normal derivative of the scalar must be computed geometrically; they are not separately frozen.

No extra tangential harmonic rows or zero tangential-vector values have been imposed. If (1) is solved, \(C_1^{\rm wave}=0\) follows in the bulk and hence at its regular boundary by the operator identity. It is not a sixth or later independent boundary condition. Proper-clock or angularly labeled observables require the corresponding transport corrections when compared between gauges.

## 4. Closed three-amplitude system for the actual quadrupole

Let \(c(t)=\cosh t\), \(H(t)=\tanh t\), and \(\lambda=8\). Use

\[
\xi^y=r(y,t)Y,\qquad
\xi^t=b(y,t)Y,\qquad
\xi^i=d(y,t)D_\kappa^iY.
\]

The symbol \(d\) in this section is the angular gauge-vector amplitude, not a boundary defining function or the previous scalar-reconstruction coefficient. The notation \(r\) denotes an amplitude, not the physical radius \(\rho\).

The required de Sitter identities are

\[
\overline{\operatorname{div}}V=(b_t+3Hb-\lambda d)Y,
\qquad
\bar\Box(rY)=(-r_{tt}-3Hr_t-\lambda r/c^2)Y,
\]
\[
(\bar\Box+3)V^t=
[-b_{tt}-3Hb_t+(3+3H^2-\lambda/c^2)b-2\lambda Hd]Y,
\]
\[
(\bar\Box+3)V^i=
[-d_{tt}-5Hd_t+(4-\lambda)d/c^2+2Hb/c^2]D_\kappa^iY.
\]

Consequently (1) is precisely the coupled system

\[
\boxed{\begin{aligned}
r_{yy}+4ar_y+\rho^{-2}(-r_{tt}-3Hr_t-\lambda r/c^2)
-2a(b_t+3Hb-\lambda d)\qquad\\
-(4a_y+8a^2)r=6(u_y+4au),
\end{aligned}}                                                \tag{6}
\]
\[
\boxed{b_{yy}+6ab_y+\rho^{-2}
[-b_{tt}-3Hb_t+(3+3H^2-\lambda/c^2)b-2\lambda Hd-2ar_t]=0,}     \tag{7}
\]
\[
\boxed{d_{yy}+6ad_y+\rho^{-2}
[-d_{tt}-5Hd_t+(4-\lambda)d/c^2+2Hb/c^2+2ar/c^2]=0.}           \tag{8}
\]

The boundary and initial data reduce to

\[
\boxed{r|_b=0,\qquad b_y|_b=d_y|_b=0,}
\qquad
\boxed{r=b=d=r_t=b_t=d_t=0\quad\hbox{at }t=0.}                 \tag{9}
\]

The radial source and the covariant operator preserve this scalar angular representation, so no coexact vector amplitude is generated by these equations. Boundary rows (9) also preserve the representation. This closure is a property of the displayed operator and data, not a demonstrated solution construction.

## 5. Checks, initial match and remaining work

At \(t=0\), equations (6)–(8) and the zero initial vector data imply

\[
r_{tt}=-6\rho^2(u_y+4au),\qquad b_{tt}=d_{tt}=0.
\]

This agrees with the already derived initial harmonic shift acceleration \(\partial_t^2\xi^y=\dot\beta_1^y\). No sign or factor change is needed between the initial completion and the proposed future gauge-wave system.

Two independent consistency checks of the tangential formulas are useful. For any scalar \(\chi\),

\[
(\bar\Box+\operatorname{Ric})\bar\nabla\chi
=\bar\nabla(\bar\Box\chi+6\chi),
\]

since \(\operatorname{Ric}[\bar g]=3\bar g\). Substituting \(b=-\chi_t\), \(d=\chi/c^2\) in the displayed component formulas reproduces this identity. For the de Sitter boost mode \(\lambda=3\), \(b=1\), \(d=H\), both displayed tangential amplitudes vanish, as required for a Killing vector. A tangential de Sitter Killing vector independent of \(y\) also solves the full homogeneous warped system with \(X=0\).

The zero initial vector and velocity preserve the physical Cauchy data. The initial compact shell collar supplies the elementary zero corner data; higher compatibility and regularity must be checked at the level of any future chosen solution theorem. This note neither selects an inner artificial boundary nor proves continuation through contact.

Once a solution \(\xi\) is actually constructed, the complete first metric-scalar pair must be transformed. At second order, the metric transformation includes \(2\mathcal L_\xi h_1+\mathcal L_\xi^2g_0\), and the scalar has the analogous terms. The bulk and physical boundary sources must then be recomputed or covariantly transported in that same representative. The longitudinal source cannot be retained unchanged, and tangential clock transport cannot be suppressed by relabeling (5) as zero tangential-vector data.

The independent root calculation agrees with (2)–(9), including the tangential potential cancellation and the initial acceleration. This is internal algebraic cross-checking only. The gauge-wave solve and a quantitative bound for its solution remain future work.
