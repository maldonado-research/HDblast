# Explicit initial relative-harmonic completion of the physical second jet

15 October 2026. This is a local initial-gauge calculation, conditional on the retained background and the physical Cauchy jets in `EXPLICIT_COMPACT_SECOND_INITIAL_JET.md`. It gives five explicit initial gauge velocities at each of the first two parameter orders. It does not solve a gauge wave equation, evolve the second perturbation, or solve a shell initial-boundary problem. All frozen inputs remain unchanged.

## 1. Input and notation

The actual background and initial physical conventions are

\[
g_0=dy^2+\rho(y)^2[-dt^2+\cosh^2t\,\kappa_{S^3}],\quad
\gamma_0=dy^2+\rho^2\kappa,\quad
a=\rho_y/\rho,
\]
\[
K=\tfrac12\mathcal L_n\gamma,\qquad \Pi=n\phi,
\qquad K_0=\Pi_0=0.
\]

Here \(n\) is the future initial-slice unit normal, not the outward shell normal. The fields \(K,\Pi\) are geometric variables, not momentum densities. Write

\[
U=uY,\quad T=u_tY,\quad V=vY,\quad T_\phi=v_tY,
\qquad f=\sum_{l=0,2,4}f_lW_l,\quad q=\frac F3Y^2,
\]

where \(Y=n_1^2-n_2^2\), \(W_l=(1,Z_2,H_4)\), and \(f_l,F,G\) are the explicit profiles in the physical-jet input. The symbol \(q\) here is the scalar coefficient of \(K_2\), not the first master bump. Thus

\[
\begin{array}{lll}
(\gamma_1)_{yy}=-4U,&(\gamma_1)_{yi}=0,&(\gamma_1)_{ij}=2U(\gamma_0)_{ij},\\
(K_1)_{yy}=-2T/\rho,&(K_1)_{yi}=0,&(K_1)_{ij}=(T/\rho)(\gamma_0)_{ij},\\
(\gamma_2)_{yy}=(\gamma_2)_{yi}=0,
&&(\gamma_2)_{ij}=f(\gamma_0)_{ij},\\
K_2=q\gamma_0,&&\Pi_2=(G-F_y)Y^2/p,\quad \Pi_1=T_\phi/\rho.
\end{array}                                                     \tag{1}
\]

All parameter expansions use \(X_\epsilon=X_0+\epsilon X_1+\epsilon^2X_2/2\). We keep the first lapse value from the longitudinal representative and choose the second lapse and both shift values to vanish:

\[
N_\epsilon=\rho(1+\epsilon U),\qquad \beta_\epsilon^I=0
\quad\hbox{on }t=0.                                             \tag{2}
\]

Equation (2) fixes lapse and shift **values**, not their time derivatives. In particular \(N_2=0\) does not imply \((g_{tt})_2=0\). The initial spacetime metric is built by

\[
g=-N^2dt^2+\gamma_{IJ}(dx^I+\beta^Idt)(dx^J+\beta^Jdt).
\]

The relative harmonic vector and reference contraction are

\[
C^A=g^{BC}(\Gamma[g]^A_{BC}-\widehat\Gamma^A_{BC}),\qquad
\widehat\Gamma=\Gamma[g_0],\qquad R^A=g^{BC}\widehat\Gamma^A_{BC}.
\]

The contraction defining \(R\) uses the perturbed inverse metric.

## 2. Exact initial gauge prescription before Taylor expansion

The contracted-Christoffel identity gives, at zero initial shift,

\[
\Gamma[g]^t=\frac{\operatorname{tr}_\gamma K}{N}-\frac{\dot N}{N^3},
\qquad
\Gamma[g]^I=-N^{-2}\dot\beta^I+\Gamma_\gamma^I
-N^{-1}\gamma^{IJ}\partial_JN.
\]

Therefore the exact local prescription is

\[
\boxed{\dot N=N^2\operatorname{tr}_\gamma K-N^3R^t,\qquad
\dot\beta^I=N^2(\Gamma_\gamma^I-R^I)-N\gamma^{IJ}\partial_JN.}       \tag{3}
\]

It enforces \(C=0\) initially for every sufficiently small parameter. It changes no physical Cauchy datum. The metric/scalar velocities remain

\[
\dot\gamma=2NK,\qquad \dot\phi=N\Pi                              \tag{4}
\]

because the shift vanishes identically as a function of the spatial coordinates on the initial slice. In the present background \(R^t=0\) on that slice: the reference \(\widehat\Gamma^t_{ty}=a\) has zero inverse-metric contraction, while the time derivatives of \(\cosh^2t\) vanish at \(t=0\).

For an all-point expansion of (3), put

\[
A=1-4\epsilon U,\qquad B=1+2\epsilon U+\epsilon^2f/2,\qquad
P=1+\epsilon U.
\]

Then \(\gamma=A\,dy^2+\rho^2B\kappa\), \(N=\rho P\), and

\[
\operatorname{tr}_\gamma K
=\frac{-2\epsilon T/\rho+\epsilon^2q/2}{A}
+3\frac{\epsilon T/\rho+\epsilon^2q/2}{B}.
\]

The reference radial contraction is \(R^y=-a/P^2-3a/B\). Cancellation of the reference sphere connection gives the following coordinate-invariant sphere formulas, with \(D^i=\kappa^{ij}D_j\):

\[
\boxed{\dot N=\rho^2P^2\left[
\frac{-2\epsilon T/\rho+\epsilon^2q/2}{A}
+3\frac{\epsilon T/\rho+\epsilon^2q/2}{B}\right],}                 \tag{5}
\]
\[
\boxed{\dot\beta^y=\rho^2P^2\left[
\frac{A_y}{2A^2}-\frac{3B_y}{2AB}-\frac{4a}{A}
-\frac{P_y}{AP}+\frac{a}{P^2}+\frac{3a}{B}\right],}              \tag{6}
\]
\[
\boxed{\dot\beta^i=-\frac{P^2}{B}
\left(\frac{D^iA}{2A}+\frac{D^iB}{2B}+\frac{D^iP}{P}\right).}    \tag{7}
\]

Equations (5)–(7) are exact local formulas at \(t=0\), for the polynomial physical data and the lapse values (2). The physical constraints of that polynomial need hold only through second order for the present application. Exact initial \(C=0\) itself does not require exact constraints.

## 3. Explicit first and second gauge velocities

Expansion of (5)–(7), retaining the factorial convention, gives

\[
\boxed{\dot N_1=\rho T,\qquad
\dot\beta_1^y=-6\rho^2(U_y+4aU),\qquad \dot\beta_1^i=0,}         \tag{8}
\]
\[
\boxed{\dot N_2=4\rho^2q-24\rho UT,}                            \tag{9}
\]
\[
\boxed{\dot\beta_2^y=\rho^2\left[
-74UU_y-\frac32f_y-194aU^2-3af\right],\qquad
\dot\beta_2^i=22U D^iU-\frac12D^if.}                           \tag{10}
\]

For example \(\operatorname{tr}_0K_1=T/\rho\),
\(\langle\gamma_1,K_1\rangle_0=14UT/\rho\), and
\(\operatorname{tr}_0K_2=4q\). Hence

\[
(\operatorname{tr}_\gamma K)_2=4q-28UT/\rho,
\]

and multiplication by \(N^2\) adds \(4UT/\rho\) to that second-order trace contribution, yielding (9). In (6), the bracket has coefficients

\[
\epsilon[-6U_y-24aU]
+\epsilon^2[-25UU_y-\tfrac34 f_y-49aU^2-\tfrac32af].
\]

Multiplication by \(P^2\) and the factor two converting the quadratic coefficient to the second variation yield (10). In (7) the first coefficient cancels, and its bracket is \(\epsilon^2[\tfrac14D^if-11U D^iU]\).

The complete second spacetime metric values include

\[
(g_{tt})_2=-2\rho^2U^2,\quad (g_{tI})_2=0,\quad
(g_{yy})_2=(g_{yi})_2=0,\quad(g_{ij})_2=f(\gamma_0)_{ij}.          \tag{11}
\]

The corresponding second coordinate velocities are

\[
\begin{aligned}
(\dot g_{tt})_2&=-8\rho^3q+44\rho^2UT,\\
(\dot g_{ty})_2&=\dot\beta_2^y-8U\dot\beta_1^y
 =\rho^2[-26UU_y-\tfrac32f_y-2aU^2-3af],\\
(\dot g_{ti})_2&=\rho^2\kappa_{ij}\dot\beta_2^j
 =\rho^2[22U D_iU-\tfrac12D_if],\\
(\dot g_{yy})_2&=2\rho q-8UT,\qquad (\dot g_{yi})_2=0,\\
(\dot g_{ij})_2&=(2\rho q+4UT)(\gamma_0)_{ij},\\
\dot\phi_2&=\rho\Pi_2+2UT_\phi.
\end{aligned}                                                  \tag{12}
\]

The mixed conversion is especially necessary: \(\dot g_{tI}=\gamma_{IJ}\dot\beta^J\) at zero shift, so its second variation contains \(2(\gamma_1)_{IJ}\dot\beta_1^J\). Setting \((\dot g_{ty})_2=\dot\beta_2^y\) loses the term \(-8U\dot\beta_1^y\). Similarly \(\dot g_{tt}=-2N\dot N\) contributes \(-4N_1\dot N_1\) to the second variation. All these terms are included in (12).

## 4. Relation to the longitudinal first representative

The original longitudinal first perturbation has

\[
C_{1,y}^{\mathrm{long}}=-6(U_y+4aU),\qquad
C_{1,t}^{\mathrm{long}}=C_{1,i}^{\mathrm{long}}=0.                 \tag{13}
\]

The radial raised component is identical since \((g_0)_{yy}=1\). The longitudinal lapse velocity already equals (8), but its shift velocity is zero. Equation (8) changes

\[
(\dot h_{ty})_1=\dot\beta_1^y=\rho^2C_1^{y,\mathrm{long}},
\]

while leaving the initial first metric values, spatial metric velocities and scalar velocity unchanged. Thus it preserves the physical first Cauchy tangent, but not the entire first spacetime one-jet in the longitudinal gauge.

If a first generator \(\xi\) solves

\[
(\Box_0+\operatorname{Ric}_0)\xi=-C_1^{\mathrm{long}},\qquad
\xi|_{t=0}=\partial_t\xi|_{t=0}=0,
\]

then its initial acceleration is fixed locally by that equation:

\[
\partial_t^2\xi^y=\rho^2C_1^{y,\mathrm{long}}
=\dot\beta_1^y,\qquad \partial_t^2\xi^t=\partial_t^2\xi^i=0.
\]

Consequently the initial first jet constructed here agrees with the one such a gauge-wave conversion would produce. This identifies its initial acceleration, not its subsequent solution or shell boundary behavior. Smooth compactly supported local coordinate jets can also realize the initial change without solving that evolution problem.

The second relative-harmonic condition is

\[
C_2=DC[g_2]+D^2C[g_1^{\mathrm{wave}},g_1^{\mathrm{wave}}]=0.
\]

Its quadratic term uses the corrected first metric velocity. Equations (9)–(12) supply an initial second jet satisfying this full condition. If one instead obtains the wave representative by a first gauge transformation, the second spacetime transformation contains

\[
g_2^{\mathrm{wave}}=g_2^{\mathrm{seed}}+\mathcal L_{\xi_2}g_0
+2\mathcal L_\xi g_1^{\mathrm{long}}+\mathcal L_\xi^2g_0,
\]

and the analogous scalar formula. The quadratic transformation terms cannot be discarded. A transformed first evolution also changes its quadratic bulk and boundary sources. Initial gauge completion supplies no justification for attaching the old longitudinal forcing unchanged to a harmonic wave problem.

## 5. Initial normal subsidiary derivatives: the precise conditional statement

Let

\[
E_{AB}=\operatorname{Ric}_{AB}-\phi_A\phi_B-\tfrac23\mathcal U(\phi)g_{AB},\qquad
\mathscr H_{AB}=E_{AB}-\tfrac12g_{AB}\operatorname{tr}_gE.
\]

Here \(\mathcal U\) denotes the scalar potential, to distinguish it from the packet amplitude \(U=uY\). The physical Hamiltonian and momentum constraints are the normal projections of \(\mathscr H\), rather than \(E_{nn}\) alone. Suppose a reduced solution jet through these Cauchy data satisfies

\[
E_{AB}=\nabla_{(A}C_{B)}
\]

through the required parameter orders. On the initial slice \(C=0\) for the exact gauge prescription, so tangential derivatives vanish. With \(X_A=(\nabla_nC)_A\), an orthonormal \((n,e_I)\) frame gives

\[
E_{nn}=X_n,\qquad E_{nI}=\tfrac12X_I,\qquad E_{IJ}=0,\qquad
\operatorname{tr}E=-X_n,
\]
\[
\mathscr H_{nn}=\tfrac12X_n,\qquad \mathscr H_{nI}=\tfrac12X_I.
\]

The verified first and second physical constraints therefore give

\[
(\nabla_nC)_1=(\nabla_nC)_2=0,
\qquad (\partial_tC)_1=(\partial_tC)_2=0
\quad\hbox{initially}.                                        \tag{14}
\]

When differentiating this statement, lower-order normal/frame/connection variations multiply already zero lower-order constraint or subsidiary jets. They do not create an omitted second-order source. Equation (14) uses the reduced Einstein equations to determine the needed second time derivatives; those derivatives are not extra independently prescribed Cauchy data. With the scalar equation, Bianchi subsequently gives the homogeneous subsidiary equation, but no solution or boundary uniqueness result is being supplied here.

## 6. Support, derivative count and checks

All formulas are finite differentiation and multiplication by smooth fixed background factors. They preserve the same closed packet annulus and its flat endpoints, and vanish on the initial open shell and vertex collars. The first radial gauge velocity has degree two. The second velocities have scalar degrees zero, two, four and their sphere gradients, with no independent coexact vector mode. The construction applies identically to the two copies.

For profile values, the new highest radial derivative is \(f_y\); since \(f\) uses first-master derivatives through order two, this needs master order three. One physical spatial derivative of the gauge velocities adds \(f_{yy}\), requiring master order four. This matches the derivative count already needed for one spatial derivative of the physical scalar correction \(\phi_2=-3\sum f_{l,y}W_l/(2p)\). This is a derivative count, not a numerical norm enclosure; the latter additionally needs all displayed coefficient bounds.

`check_initial_harmonic_completion.py` independently forms the five-dimensional inverse metric series, all Christoffel symbols, and the contraction \(g^{BC}(\Gamma[g]-\Gamma[g_0])^A_{BC}\), with sphere normal coordinates only for evaluating point jets. Three nonspecial rational point jets verify both parameter orders and independently recover (13). Six mutation controls detect a missing lapse quadratic term, a halved radial quadratic term, a missing \(f_y\), a halved angular quadratic term, the missing mixed metric conversion, and the missing second metric lapse value. The all-point proof is the exact ADM/contracted-Christoffel calculation (3)–(7), not the finite point sample.

```sh
python check_initial_harmonic_completion.py --run-controls --output INITIAL_HARMONIC_COMPLETION.json
python -O check_initial_harmonic_completion.py --run-controls --output INITIAL_HARMONIC_COMPLETION_OPTIMIZED.json
```

Both modes pass with identical mathematical reports. This continuation closes the initial five gauge-velocity choices for the selected physical second jet. The first gauge-wave solve, its transformed evolving sources, validated second evolution and shell boundary problem remain separate tasks.

Input hashes:

- Physical-jet note: `2cca9c9b87dfa5c2c21712ae81ea0c91d01b939b902e82d54dd62a6eae6bf0ed`.
- Frozen GH boundary/source note: `5b28b3065515a3b8b1fe304862f7982780d1bcc064f15c733c038f66d281977c`.
- Prior general initial gauge/subsidiary note: `74bce09080d4015b2eccf8c183b33f47b80516216a9547837e3e0e34251d7b03`.
