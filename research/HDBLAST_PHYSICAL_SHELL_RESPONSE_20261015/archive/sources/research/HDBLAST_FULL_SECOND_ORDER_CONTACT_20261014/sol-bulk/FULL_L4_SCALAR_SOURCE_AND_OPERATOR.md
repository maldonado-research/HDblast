# Full l=4 quadratic scalar bulk residual and its general linear operator

14 October 2026. Independent derivation in the inherited first longitudinal
gauge. This is a bulk differential calculation; it is not a coupled shell
contact existence theorem, a novelty assertion, or a cosmological claim.

## 1. Conventions and angular projection

Use

\[
g_0=dy^2+\rho^2(-dt^2+c^2\kappa),\quad
c=\cosh t,\quad b=\tanh t,\quad a=\rho_y/\rho,
\quad p=\phi_{0,y},
\]
\[
\mathscr E_{AB}=\operatorname{Ric}_{AB}-\phi_{,A}\phi_{,B}
 -\tfrac23U g_{AB},\qquad
\mathscr S=\Box_g\phi-U_{,\phi}.
\]

The first fields are \(\Psi=u(y,t)Y\), \(\chi=v(y,t)Y\),
\(h_{yy}=-4\Psi\), \(h_{y\mu}=0\),
\(h_{\mu\nu}=2\rho^2\Psi\gamma_{\mu\nu}\), with
\(v=-3(\rho^2u)_y/(\rho^2p)\) for the actual inherited solution.
The source identities below do **not** require substituting that relation,
the first field equation, or the background equations.

Write \(g_\epsilon=g_0+\epsilon h+\epsilon^2 k/2+\cdots\),
\(\phi_\epsilon=\phi_0+\epsilon\chi+\epsilon^2w/2+\cdots\).
The complete quadratic residual is \(J=D^2\mathscr E[h,\chi;h,\chi]\),
\(J_S=D^2\mathscr S[h,\chi;h,\chi]\). It equals twice the
\(\epsilon^2\) coefficient along the linear field path. Thus the
second equations are \(D\mathscr E[k,w]=-J\),
\(D\mathscr S[k,w]=-J_S\).

On the unit sphere set

\[
Y=n_1^2-n_2^2,\quad F_2=n_1^2+n_2^2-\tfrac12,\quad
Q=Y^2-\tfrac12(n_1^2+n_2^2)+\tfrac1{12}.
\]

Then \(\Delta Y=-8Y\), \(\Delta Q=-24Q\), and
\(Y^2=Q+F_2/2+1/6\). Products \(Y^2\) have coefficient one
in the \(Q\) channel. Since
\(\Delta(Y^2)=-16Y^2+2|DY|^2\), the \(Q\) coefficient of
\(|DY|^2\) is \(-4\). The scalar part of
\((DY\otimes DY)^{\rm TF}\) has \(Q\)-component
\((\operatorname{Hess}Q)^{\rm TF}/3\). A direct component check is
provided in the companion checker.

## 2. Complete quadratic l=4 rows

Define

\[
\begin{split}
K={}&8u_y^2+8u u_{yy}+104a u u_y
 +48(a_y+4a^2)u^2
 +\tfrac23U_{,\phi\phi}v^2+\tfrac83U_{,\phi}uv,\\
A={}&-4u u_{tt}-12b u u_t-4u_t^2-48c^{-2}u^2-\rho^2K.
\end{split}
\]

All potential derivatives are evaluated at \(\phi_0(y)\). The three
base residual coefficients, \(J_{ab}|_{l=4}=j_{ab}Q\), are

\[
\boxed{\begin{split}
j_{yy}={}&-8u_y^2+16u u_{yy}-32a u u_y-2v_y^2
 -\tfrac23U_{,\phi\phi}v^2+\tfrac{16}3U_{,\phi}uv\\
&+\rho^{-2}(8u u_{tt}+24b u u_t-16u_t^2),\\
j_{yt}={}&12u u_{yt}-48a u u_t-2v_yv_t,\\
j_{tt}={}&24u u_{tt}+12u_t^2-2v_t^2-A.
\end{split}}
\tag{1}
\]

The mixed scalar gradients are

\[
\boxed{\begin{split}
J_{yi}|_{l=4}&=(6u u_y-24a u^2-vv_y)D_iQ,\\
J_{ti}|_{l=4}&=(18u u_t-12b u^2-vv_t)D_iQ.
\end{split}}\tag{2}
\]

Define the angular trace coefficient by
\(\kappa^{ij}J_{ij}|_{l=4}/(3c^2)=j_\Omega Q\). Then

\[
\boxed{j_\Omega=A-24b u u_t
 +c^{-2}(-80u^2+\tfrac83v^2).}\tag{3}
\]

The angular scalar trace-free coefficient is

\[
\boxed{J_{ij}^{\rm scalar,TF}|_{l=4}
 =(8u^2-\tfrac23v^2)(\operatorname{Hess}Q)^{\rm TF}_{ij}.}
\tag{4}
\]

The scalar residual is

\[
\boxed{\begin{split}
j_S={}&12u_yv_y+8u v_{yy}+32a u v_y+48p u u_y
 +32u^2(p_y+4ap)\\
&+4\rho^{-2}u(v_{tt}+3b v_t+8c^{-2}v)
 -U_{,\phi\phi\phi}v^2,
\qquad J_S|_{l=4}=j_SQ.
\end{split}}\tag{5}
\]

The use of \(p_y+4ap\) preserves an off-background identity; on the
actual background it equals \(U_{,\phi}\). No source is omitted by
leaving this substitution undone.

## 3. Exact radial derivation

For the first linear metric path set
\(N^2=1-4\epsilon\Psi\),
\(q_{\mu\nu}=\rho^2e^{2\omega}\gamma_{\mu\nu}\),
\(e^{2\omega}=1+2\epsilon\Psi\), and \(\nu=\log N\).
Put \(B=a+\omega_y\). Direct radial curvature identities are

\[
\begin{split}
R_{yy}&=-4(B_y+B^2-B\nu_y)-N\Delta_qN,\\
R_{y\mu}&=-3(\omega_{y\mu}-B\nu_{,\mu}),\\
R_{\mu\nu}&=3\gamma_{\mu\nu}
 -2\omega_{;\mu\nu}-\nu_{;\mu\nu}
 +2\omega_{,\mu}\omega_{,\nu}-\nu_{,\mu}\nu_{,\nu}\\
&\quad+\omega_{,\mu}\nu_{,\nu}+\omega_{,\nu}\nu_{,\mu}
 -\gamma_{\mu\nu}(\Box_\gamma\omega+2|d\omega|^2
                         +d\omega\cdot d\nu)\\
&\quad-\rho^2e^{2(\omega-\nu)}
 (B_y+4B^2-B\nu_y)\gamma_{\mu\nu}.
\end{split}\tag{6}
\]

Also

\[
\Box_g\phi=e^{-2\nu}[\phi_{yy}+(4B-\nu_y)\phi_y]
 +\rho^{-2}e^{-2\omega}
 [\Box_\gamma\phi+(2d\omega+d\nu)\cdot d\phi].\tag{7}
\]

Expanding \(\omega=\epsilon\Psi-\epsilon^2\Psi^2\),
\(\nu=-2\epsilon\Psi-4\epsilon^2\Psi^2\), and adding **both**
the scalar-gradient and potential terms gives (1)–(5).
The tangential non-trace quadratic residual before sphere decomposition is

\[
J_{\mu\nu}\equiv24\Psi\Psi_{;\mu\nu}
 +12\Psi_{,\mu}\Psi_{,\nu}-2\chi_{,\mu}\chi_{,\nu}
 \pmod{\gamma_{\mu\nu}}.\tag{8}
\]

Equation (8), together with the retained trace in (6), fixes the signs
in (2)–(4). Its TT coefficient is the inherited
\(16u^2+8v^2/3\), consistently separate from (4).

## 4. General scalar-harmonic linear operator on the two-dimensional base

Let \(x^a=(y,t)\), \(q_{ab}dx^adx^b=dy^2-\rho^2dt^2\), and
\(R=\rho c\). Let \(D\) be the Levi-Civita derivative of \(q\),
\(r_a=D_aR/R\), and \(P_a=D_a\phi_0=(p,0)\).
For \(l=4\), the three scalar gauge functions can set the two mixed
metric gradients and the angular scalar trace-free metric amplitude to
zero. In that second-order gauge write

\[
k_{ab}=f_{ab}Q,\qquad k_{ai}=0,\qquad
k_{ij}=2R^2H Q\kappa_{ij},\qquad \chi_2=wQ.
\tag{9}
\]

This is an order-two gauge choice. The first longitudinal gauge is
unchanged, so the quadratic source (1)–(5) is unchanged. Define

\[
f=q^{ab}f_{ab},\quad
J^a{}_{bc}=\tfrac12q^{ad}(D_bf_{cd}+D_cf_{bd}-D_df_{bc}),
\quad T=\tfrac12f+3H,\quad T_a=D_aT,
\]
\[
Z^a=r_bf^{ab}-2r^aH-D^aH,\qquad \lambda=24.
\]

The complete linear base Einstein operator is

\[
\boxed{\begin{split}
\mathcal E_{ab}={}&D_cJ^c{}_{ab}+3r_cJ^c{}_{ab}
 +\frac\lambda{2R^2}f_{ab}
 -3(r_aD_bH+r_bD_aH)-D_bT_a\\
&-P_aD_bw-P_bD_aw
 -\tfrac23(Uf_{ab}+U_{,\phi}wq_{ab}).
\end{split}}\tag{10}
\]

Its two mixed gradient coefficients are

\[
\boxed{\mathcal E_a=
 \tfrac12(D_bf^b{}_a-D_af+r^bf_{ab}+r_af)
 -2D_aH-P_aw.}\tag{11}
\]

The angular tensor equation is

\[
\boxed{D\mathscr E_{ij}=
 B(\operatorname{Hess}Q)_{ij}+R^2LQ\kappa_{ij},}
\tag{12}
\]
\[
B=-\tfrac12f-H,\qquad
L=D_aZ^a+3r_aZ^a+2r^aD_aH-r^aT_a
 +\lambda H/R^2-\tfrac43UH-\tfrac23U_{,\phi}w.
\]

Thus the normalized angular trace row is
\(\mathcal E_\Omega=L-\lambda B/(3R^2)\), and the scalar-TF
row has coefficient \(B\). The scalar operator is

\[
\boxed{\begin{split}
\mathcal S={}&D^aD_aw+3r^aD_aw-(\lambda/R^2+U_{,\phi\phi})w
 -f^{ab}D_aP_b\\
&-P^b(D_af^a{}_b+3r_af^a{}_b
                  -\tfrac12D_bf-3D_bH).
\end{split}}\tag{13}
\]

These formulas follow from the connection difference
\(C^A{}_{BC}=\frac12g_0^{AD}(\nabla_Bk_{CD}+\nabla_Ck_{BD}
-\nabla_Dk_{BC})\), with
\(\delta R_{AB}=\nabla_CC^C{}_{AB}-\nabla_BC^C{}_{AC}\).
The nonzero block amplitudes are

\[
C^a{}_{bc}=J^a{}_{bc}Q,\quad
C^i{}_{ab}=-f_{ab}D^iQ/(2R^2),\quad
C^a{}_{bi}=f^a{}_bD_iQ/2,
\]
\[
C^i{}_{aj}=D_aH Q\delta^i_j,\quad
C^a{}_{ij}=R^2Z^aQ\kappa_{ij},\quad
C^i{}_{jk}=H(\delta^i_jD_kQ+\delta^i_kD_jQ-\kappa_{jk}D^iQ).
\tag{14}
\]

The source equations are therefore
\(\mathcal E_{ab}=-j_{ab}\),
\(\mathcal E_y=-(6uu_y-24au^2-vv_y)\),
\(\mathcal E_t=-(18uu_t-12bu^2-vv_t)\),
\(\mathcal E_\Omega=-j_\Omega/\rho^2\),
\(B=-(8u^2-2v^2/3)\), and \(\mathcal S=-j_S\).

## 5. What the inherited one-master ansatz cannot close

If one uses only the inherited conformal second field, then
\(f_{yy}=-4\Psi_2\), \(f_{yt}=0\),
\(f_{tt}=-2\rho^2\Psi_2\), \(H=\Psi_2\).
It has \(f=-2\Psi_2\), hence \(B=0\) identically. Equation (4)
then cannot be balanced except under the nongeneric identity
\(v^2=12u^2\). Even if that identity holds at one point, it does
not establish any of the other sourced rows.
Adding a pure order-two diffeomorphism cannot repair this failure:
\(D\mathscr E[\mathcal L_\xi g_0,\mathcal L_\xi\phi_0]
=\mathcal L_\xi\mathscr E_0=0\) on the solved background.

For a nonzero compact displacement profile on a connected collar where
\(p\ne0\), that exceptional identity cannot hold throughout a time
slice. Substituting the actual first reconstruction would give
\(|u_y+2au|=2|p||u|/\sqrt3\), hence
\(|u_y|\le(2|a|+2|p|/\sqrt3)|u|\). Beginning outside the compact
support, the integral form and Gronwall force \(u=0\) throughout the
collar, contrary to the assumed nonzero profile. Thus (4) is nonzero
somewhere on every such nonzero displacement slice. Pure-velocity first
data have \(u=v=0\) on their initial slice and are not covered by that
claim.

Every general scalar solution instead obeys the compulsory algebraic row

\[
\boxed{f+2H=16u^2-\tfrac43v^2.}\tag{15}
\]

In component conventions \(f_{yy}=N\), \(f_{yt}=B_{yt}\),
\(f_{tt}=\rho^2T\), and \(k_{ij}=R^2KQ\kappa\), this reads
\(N-T+K=16u^2-4v^2/3\).

An algebraic particular completion of (15) is necessary, but it is not
by itself a proof that one forced inherited matter master suffices for
all remaining rows. The general scalar metric has an independent
gravitational scalar polarization. There are seven scalar metric
amplitudes before gauge fixing, three scalar gauge functions, and the
matter scalar; equivalently (9) has five invariant amplitudes before
constraints. The Hamiltonian and two scalar momentum constraints leave
two physical Cauchy pairs in the generic \(l\ge2\) sector. The
inherited reconstruction proves a sufficient one-pair family, not its
completeness. A full reduction must retain or prove consistent removal of
the second pair, as well as solve the algebraic completion and sourced
constraints. Equations (10)–(15) supply an explicit system with which
such a claim can be tested.

An explicit mixed row makes the additional completion visible. In these
component conventions (11) is

\[
\begin{split}
\mathcal E_y&=\tfrac12T_y+\tfrac32aN
 -\tfrac1{2\rho^2}(\partial_tB_{yt}+bB_{yt})-K_y-pw,\\
\mathcal E_t&=\tfrac12\partial_yB_{yt}-\tfrac12N_t
 +aB_{yt}+\tfrac12bN-bT-K_t.
\end{split}\tag{16}
\]

Set \(d=16u^2-4v^2/3\) and \(G=T+K\). Using (15), the
forced temporal mixed row becomes

\[
\boxed{(\partial_y+2a)B_{yt}-(\partial_t+b)G
 =-4u u_t-\tfrac23v v_t+8b u^2+\tfrac43b v^2.}\tag{17}
\]

The inherited conformal values have \(B_{yt}=G=0\). Even an
algebraic repair that changes only \(N\) by \(d\) keeps these
two values zero and fails (17) generically. Hence the scalar-TF trace
repair alone does not close the actual mixed system. Equation (17) is
a forced constraint relation for extra scalar gravitational amplitudes;
it is not itself a proof of a new propagating master or of its boundary
well-posedness.

On a controlled collar where \(p\ne0\), the two mixed rows give an
explicit constrained reconstruction in terms of \(K,G,B_{yt}\):

\[
N=d+G-2K,\qquad T=G-K,
\]
\[
\boxed{pw=\tfrac12(G_y+3aG)-\tfrac32(K_y+2aK)
 +\tfrac32ad-\frac{\partial_tB_{yt}+bB_{yt}}{2\rho^2}
 +(6u u_y-24au^2-vv_y).}\tag{18}
\]

Equation (17) specifies \((\partial_y+2a)B_{yt}\). For example,
prescribing the even-shell value \(B_{yt}(y_b,t)=0\) gives its exact
radial Volterra reconstruction from \(G_t+bG\) and the known source.
Equations (10), (12), and (13) must still be imposed and shown to reduce
to an admissible coupled evolution. Formula (18) retains the inherited
matter reconstruction only in its homogeneous \(G=B_{yt}=d=0\) part.
The collar hypothesis does not prove regularity of this division at the
cone where \(p=0\).

There is a precise local consistency identity, short of a boundary
propagation theorem. The exact residuals obey

\[
\nabla^A\mathscr E_{AB}
 -\tfrac12\nabla_B(g^{AC}\mathscr E_{AC})
 =-\mathscr S\,\phi_{,B}.\tag{19}
\]

On a solved background and solved first field, taking two parameter
derivatives gives
\(\nabla_0^AJ_{AB}-\frac12\nabla_{0,B}\operatorname{tr}_0J
=-J_SP_B\). Terms containing variations of the connection, inverse
metric, or \(d\phi\) multiply the vanishing zeroth or first
residuals. Thus the complete source automatically obeys the linear
compatibility identity needed by the sourced equations.

For its scalar angular piece put \(\tau=j_\Omega/\rho^2\),
\(\beta=8u^2-2v^2/3\), \(j_a\) equal to (2),
\(\mathfrak t=q^{ab}j_{ab}+3\tau\), and
\(L_J=\tau+\lambda\beta/(3R^2)\). The harmonic compatibility rows are

\[
D^aj_{ab}+3r^aj_{ab}-\lambda j_b/R^2-3r_b\tau
 -\tfrac12D_b\mathfrak t=-P_bj_S,
\]
\[
D^aj_a+3r^aj_a+L_J+(2-\lambda)\beta/R^2
 -\tfrac12\mathfrak t=0.\tag{20}
\]

The transverse tensor forcing has zero trace and divergence and does
not enter these rows. These identities prevent treating the Einstein
and scalar forcing as unrelated arbitrary functions. They do not
prescribe incoming boundary fields or prove that constraints propagate
under the shell conditions.

## 6. Coexact-vector source selection from the complete equations

The complete base rows have only scalar products. Both mixed rows (2)
are exact angular gradients. For the angular block,

\[
D^i(D_iY D_jY)^{\rm TF}
 =D_j(-4Y^2+|DY|^2/6),
\]

and
\(Y\operatorname{Hess}Y=\frac12\operatorname{Hess}(Y^2)-DY\otimes DY\).
Thus the complete angular trace-free source also has gradient divergence.
The divergence of a nonzero vector-derived tensor built from a coexact
\(S^3\) vector harmonic is a nonzero coexact harmonic, apart from
Killing vectors whose symmetric derivative is zero. Orthogonality to
gradients therefore removes every nontrivial vector-derived source.
This is a consequence of the complete mixed and angular equations,
rather than reflection parity alone. It does not forbid independent
homogeneous vector data; choosing their coefficient zero requires the
corresponding homogeneous boundary and initial conditions.

## 7. Scope and verification status

The calculation fixes the first-order longitudinal representative and
retains the full geometric, scalar-gradient and potential residuals.
Changing that first gauge requires transforming the quadratic source.
The second-order scalar gauge choice in (9) adds only the homogeneous
linear gauge direction and does not change the first source.

The bulk rows establish a genuine obstruction to an uncompleted inherited
scalar reconstruction. They supply the forced general scalar harmonic
system and an exact zero coexact-vector forcing selection. They do not
yet prove a two-master reduction, constraint propagation, an admissible
incoming boundary count, corner compatibility, or a coupled scalar shell
contact estimate. The shell supplement must use the same Eulerian gauge,
or perform its complete material conversion.

`checker/direct_metric_checker.py` is a standalone checker requiring only
the Python standard library. It constructs metric coordinate two-jets,
Christoffels, their derivatives, all Ricci components, and the scalar
Laplacian with exact rational arithmetic. It imports no inherited scientific
code. Three unconstrained local background/first-field jets at five sphere
points compare 30 quadratic coefficients. Three unconstrained background
jets and independently assigned \(f_{yy},f_{yt},f_{tt},H,w\) jets at
three sphere points compare 24 linear coefficients, including nonzero
\(f_{yt}\). A separate two-dimensional covariant implementation evaluates
(10)–(13), so these are direct five-dimensional comparisons rather than
checks of a copied operator expression against itself.
Fifteen further exact zero checks compare the explicit component mixed
rows and forced reconstruction (16)–(18) against that verified operator.

Ordinary and optimized Python runs pass; nine wrong-formula controls fail
in each mode. They reverse Ricci, remove connection products or scalar
gradient squares, alter potential normalization, lose the Taylor factorial,
corrupt angular projection, or omit base/mixed warp and angular potential
terms. The exact sample checks supplement the analytic derivation; they
do not prove equality for every possible jet solely by finite sampling.

`verify_angular_selection.py` separately constructs the unit-sphere metric
and Christoffels in the rational Hopf coordinate \(x=\sin^2\theta\).
It proves the eigenvalue identities, complete scalar/TT decomposition,
tensor trace/divergence, and gradient divergence for all points of that
chart using symbolic simplification. All 38 symbolic component checks pass
with and without Python optimization. Three independent angular corruptions
fail in both modes. Checkers use explicit exceptions, never Python
`assert`, and accept `--output`.

`VERIFICATION_RECEIPT.json` records the final artifact hashes, input hashes,
runtime, commands and outcomes. These checks support the stated formulas;
they do not prove a boundary condition or an evolution theorem.
