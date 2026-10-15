# Complete quadratic junction sources in material Gaussian normal charts

14 October 2026 (America/Los_Angeles). Separate research continuation;
the published/frozen localized capsule is unchanged. Substantial OpenAI
Codex assistance and internal independent automated review; no external
peer review or novelty certification is asserted.

## 1. Setting and the role of the embeddings

Use the established five-dimensional Einstein–scalar action with the
same smooth bulk potential on both interiors and **one** common shell
tension \(\sigma(\varphi)\). The shell is timelike, with common intrinsic
coordinates \(y^a\), \(a=0,\ldots,3\), induced Lorentzian metric \(h\),
and scalar trace \(\varphi\). Each copy has its own outward spacelike
normal. For a smooth parameterized configuration admitting normal collars,
choose one family-adapted Gaussian normal chart on each side:

\[
g_i=dr_i^2+h_i(r_i,y;\varepsilon)_{ab}\,dy^a dy^b,
\qquad W:\ r_i=0,\qquad n_i=\partial_{r_i},
\]

\[
K_{i,ab}=\tfrac12\partial_{r_i}h_{i,ab}|_0,
\qquad p_i=\partial_{r_i}\phi_i|_0.
\tag{1}
\]

The interior coordinate can be taken \(r_i<0\); thus both increasing
\(r_i\) directions are outward. The two normal derivatives are summed,
not identified by opposite signs in a single coordinate.

All fields in this note are **material pullbacks**. The chart construction
uses the actual embedding and actual unit normal for each parameter value.
It absorbs displacement, graph evaluation, tilt, and normal-normalization
effects into the coefficient fields; those effects are not discarded.
For example, after extending the shell identification to a collar
diffeomorphism \(F_{i,\varepsilon}\), let \(V_i\) be its parameter generator
and \(W_i=\partial_\varepsilon V_i|_0\). A tensor \(Q_i\) has material
derivatives

\[
\begin{aligned}
Q_i^{\mathrm{mat},1}&=Q_i^{(1)}+\mathcal L_{V_i}Q_{i,0},\\
Q_i^{\mathrm{mat},2}&=Q_i^{(2)}+2\mathcal L_{V_i}Q_i^{(1)}
+\mathcal L_{W_i}Q_{i,0}+\mathcal L_{V_i}^2Q_{i,0},
\end{aligned}\tag{2}
\]

followed by restriction to the shell. Here \(W_i\) is the derivative of
the generator, not an unqualified coordinate acceleration. Applying (2)
to \(g_i\) and \(\phi_i\) illustrates why a background-fixed Eulerian
coefficient cannot simply be substituted for a material one. The Gaussian
collar construction also fixes the transverse chart and its unit normal.

This is a local geometric gauge choice. It is not assumed identical to
the generalized harmonic gauge used for bulk hyperbolic evolution.

For a simple kinematic illustration in a background-fixed proper-normal
coordinate \(\ell\), take an embedding
\(\ell(\varepsilon,y)=\ell_0+\varepsilon\zeta(y)+\varepsilon^2\eta(y)/2\).
Scalar pullback gives

\[
\psi_1=\phi_1|_0+p_0\zeta,\qquad
\psi_2=\phi_2|_0+2\zeta\partial_\ell\phi_1|_0
+p_0\eta+(\partial_\ell^2\phi_0)|_0\zeta^2.
\tag{2a}
\]

Even when the bare Eulerian \(\phi_1\) vanishes, the material
\(\psi_1=p_0\zeta\) need not. This is a pullback calculation, not by
itself an Einstein–scalar solution. It explicitly prevents interpreting
fixed material shell coordinates as the deletion of bending sources.

## 2. Exact matching and junctions

The configuration matching is

\[
h_+(0,y;\varepsilon)=h_-(0,y;\varepsilon)=h(y;\varepsilon),
\quad\phi_+(0,y;\varepsilon)=\phi_-(0,y;\varepsilon)
=\varphi(y;\varepsilon).
\tag{3}
\]

With the established outward-normal conventions, define the exact residuals

\[
I_{ab}=\sum_i(K_{i,ab}-K_i h_{ab})+\sigma(\varphi)h_{ab},
\quad P=\sum_i p_i+\sigma'(\varphi).
\tag{4}
\]

Taking the trace of Israel in intrinsic dimension four gives
\(\operatorname{tr}I=-3\sum_iK_i+4\sigma\). Hence the invertible
trace reversal is

\[
\mathcal R_{ab}:=I_{ab}-\tfrac13(\operatorname{tr}I)h_{ab}
=\sum_iK_{i,ab}-\tfrac13\sigma h_{ab},
\qquad I_{ab}=\mathcal R_{ab}-(\operatorname{tr}\mathcal R)h_{ab}.
\tag{5}
\]

Thus the full ten metric junctions are exactly
\(\sum_iK_{i,ab}=\sigma h_{ab}/3\), and the scalar junction is
\(\sum_ip_i=-\sigma'\). No angular-sector or scalar-master restriction
has been made.

## 3. Explicit second-order material sources

Expand every field with the derivative/factorial convention

\[
h=h_0+\varepsilon h_1+\tfrac12\varepsilon^2h_2,
\quad\varphi=\varphi_0+\varepsilon\psi_1+\tfrac12\varepsilon^2\psi_2,
\quad K_i=K_{i,0}+\varepsilon K_{i,1}+\tfrac12\varepsilon^2K_{i,2},
\quad p_i=p_{i,0}+\varepsilon p_{i,1}+\tfrac12\varepsilon^2p_{i,2}.
\tag{6}
\]

In (6), \(K_{i,j}\) denotes the lower-index tensor when it carries
\(ab\); its scalar trace must always use the varying \(h\).
All \(\sigma_0\) derivatives below are evaluated at \(\varphi_0\).
The exact matching conditions give

\[
h_{+,j}|_0=h_{-,j}|_0=h_j,
\qquad\phi_{+,j}|_0=\phi_{-,j}|_0=\psi_j,\qquad j=1,2.
\tag{7}
\]

There is no additional quadratic matching term in this material chart;
the corresponding graph terms have already entered (2).
The first-order boundary equations are

\[
\sum_iK_{i,1,ab}-\frac{\sigma_0}{3}h_{1,ab}
-\frac{\sigma'_0}{3}\psi_1h_{0,ab}=0,
\qquad\sum_ip_{i,1}+\sigma''_0\psi_1=0.
\tag{8}
\]

Twice differentiating (5) and (4) gives the complete order-two equations:

\[
\boxed{\sum_iK_{i,2,ab}-\frac{\sigma_0}{3}h_{2,ab}
-\frac{\sigma'_0}{3}\psi_2h_{0,ab}
=\frac{2\sigma'_0}{3}\psi_1h_{1,ab}
+\frac{\sigma''_0}{3}\psi_1^2h_{0,ab},}
\tag{9}
\]

\[
\boxed{\sum_ip_{i,2}+\sigma''_0\psi_2
=-\sigma'''_0\psi_1^2.}\tag{10}
\]

Because (1) holds for every parameter value,
\(K_{i,2,ab}=\partial_{r_i}h_{i,2,ab}/2\) and
\(p_{i,2}=\partial_{r_i}\phi_{i,2}\). Thus (7), (9), and (10) are explicit
inhomogeneous matching/normal-derivative conditions for all second-order
metric and scalar coefficients in these charts. Their bulk second-order
equations must still be solved.

For an exactly reflection-symmetric two-copy family, each copy has the
same \(K,p\), and the exact laws become
\(K_{ab}=\sigma h_{ab}/6\), \(p=-\sigma'/2\). Equations (9)–(10) divide
by two:

\[
K_{2,ab}-\tfrac16\sigma_0h_{2,ab}-\tfrac16\sigma'_0\psi_2h_{0,ab}
=\tfrac13\sigma'_0\psi_1h_{1,ab}
+\tfrac16\sigma''_0\psi_1^2h_{0,ab},
\quad p_2+\tfrac12\sigma''_0\psi_2=-\tfrac12\sigma'''_0\psi_1^2.
\tag{11}
\]

This divides the two-sided equations; it does not add a second shell tension.

For the **actual registered tension** of the shell-Noether checkpoint,
\(z=\varphi-1\),

\[
\sigma=2\big(\tfrac13+z^2+\tfrac13z^3\big)
+\mathfrak d(1+c_*+c_*z),\quad
\sigma'=4z+2z^2+\mathfrak dc_*,\quad
\sigma''=4\varphi,\quad\sigma'''=4.
\]

Equation (10) therefore specializes to
\(\sum_i p_{i,2}+4\varphi_0\psi_2=-4\psi_1^2\).
For an exact two-copy reflection-symmetric family its single-copy form is
\(p_2+2\varphi_0\psi_2=-2\psi_1^2\). The inhomogeneous source coefficient
is independent of detuning and of the shell background value. It is a
boundary Taylor coefficient, with no new energy, force, or observational
interpretation assigned to its sign.

## 4. Normal force and tangential consistency

Let \(E_{AB}=G_{AB}-T_{AB}\) be the bulk Einstein residual. In a spacelike
normal chart the Gauss and Codazzi constraints are

\[
E_{nn,i}=\tfrac12\{K_i^2-K_{i,ab}K_i^{ab}-R[h]}
-\tfrac12p_i^2+\tfrac12|D\varphi|_h^2+U(\varphi),
\tag{12}
\]

\[
E_{na,i}=D_bK_i^b{}_a-D_aK_i-p_iD_a\varphi.
\tag{13}
\]

Write \(\Delta K_{ab}=K_{+,ab}-K_{-,ab}\),
\(\Delta K=h^{ab}\Delta K_{ab}\), and \(\Delta p=p_+-p_-\).
Using matching and the same bulk potential, subtracting (12) and using
only the definitions (4) gives the exact **off-shell** identity

\[
\boxed{\sigma\Delta K+\sigma'\Delta p
=2(E_{nn,+}-E_{nn,-})+I_{ab}\Delta K^{ab}+P\Delta p.}
\tag{14}
\]

The terms involving intrinsic curvature, scalar tangential gradients,
and potential cancel because their induced fields are common.
Summing (13) gives the exact tangential identity

\[
\boxed{D_bI^b{}_a=E_{na,+}+E_{na,-}+P D_a\varphi.}\tag{15}
\]

Consequently, if matching, the bulk constraints, and the two junction
residuals vanish through second order, the normal-force residual and
tangential conservation do too. These are consistency consequences,
not extra independently prescribed shell equations and not a proof that
constraints propagate in an unspecified hyperbolic boundary problem.

In applying (14), second-order traces must include inverse-metric terms:

\[
(\Delta K)_2=h_0^{ab}\Delta K_{2,ab}
-2h_1^{ab}\Delta K_{1,ab}
+\big[2(h_1h_0^{-1}h_1)^{ab}-h_2^{ab}\big]\Delta K_{0,ab}.
\tag{16}
\]

All raised indices on the right use \(h_0\). On the symmetric background
\(\Delta K_{0,ab}=\Delta p_0=0\). If first-order junction residuals
vanish, the second derivative of (14) then reduces to
\((\sigma\Delta K+\sigma'\Delta p)_2=2\Delta E_{nn,2}\), with the
full trace (16). Omitting its cross term changes the force residual.

## 5. What this supplies for shell contact

Equations (7), (9), and (10) supply the explicit complete quadratic
interface algebra in a geometric gauge that includes moving embeddings.
They materially formulate the next boundary calculation. They do not
construct its solution through shell contact.

The metric equation retains every tensor component. Its source can have
trace-free angular pieces for general material \(h_1\). For the specific
scalar family with \(h_1=f(y)h_0\), both displayed metric-source terms
are pure trace. No nonzero tensor shell source is inferred merely from
the previously established tensor component of the **bulk matter** source.
A tensor second-order bulk sector could instead have homogeneous tensor
junction conditions here. Likewise, pure odd-copy first-order data can
have zero common material \(h_1,\psi_1\), even though graph effects remain
inside their second-order material coefficients.

There is a useful necessary boundary simplification for the inherited
reflection-even scalar-master family. Its representative has zero
first-order even shell displacement and
\(\delta g=2F(s\eta-3X^\flat\otimes X^\flat)\). On \(s=b^2\), every
shell tangent obeys \(X\cdot dX=0\), so
\(h_1=(2b^2F_b/\Omega_b^2)h_0\). The material induced metric is unchanged
by extending that fixed intrinsic identification to a Gaussian collar.
Consequently the right side of (9) has zero trace-free part.

For an exact two-copy reflection-symmetric continuation, write the normal
background as \(h_0(r)=\rho(r)^2\gamma_{\mathrm{dS}}\) and
\(a=\rho'/\rho\). At the shell \(a_b=\sigma_0/6\). The trace-free part
of (11) gives

\[
\tfrac12\partial_r h_{2}^{\mathrm{TF}}-a_bh_{2}^{\mathrm{TF}}=0.
\]

The trace-free projector defined by \(\gamma_{\mathrm{dS}}\) is
independent of \(r\). With \(h_2^{\mathrm{TF}}=\rho^2t_2^{\mathrm{TF}}\),
this becomes \(\partial_rt_2^{\mathrm{TF}}|_b=0\). In particular each
angular transverse trace-free tensor sector has the necessary homogeneous
Neumann condition

\[
\boxed{\partial_rt_{2,\mathrm{TT}}|_b=0.}\tag{17}
\]

This condition concerns the rescaled tensor field; setting
\(\partial_rh_{2,\mathrm{TT}}=0\) instead would omit the background warp
factor. Equation (17) remains a necessary boundary relation, even if the
full second-order bulk forcing excites a tensor sector. It supplies
neither its bulk equation nor a well-posed tensor evolution theorem.

To use these equations with background-fixed wave-map unknowns, one must
compute the material/Eulerian conversion, including the embedding and
normal-coordinate jets, or establish a compatible hyperbolic formulation
using this boundary gauge. A full forced Einstein–scalar initial-boundary
problem requires bulk source equations, gauge/constraint boundary
conditions, corner compatibility, and a well-posedness argument. None is
obtained simply by writing the normal-derivative laws (9)–(10), and the
homogeneous scalar-master Wentzell theorem does not supply them.

## 6. Verification

The independent symbolic verifier imports no scientific code from earlier
capsules. It verifies the Lorentzian tensor trace reversal, all ten
independent material metric-source components and scalar source, the
two-copy versus one-tension normalization, material pullback Taylor
factors, inverse-metric trace variations, and the algebraic off-shell
normal-force identity. It also checks the tangential identity's product
rule. Wrong-formula controls include missing factors of two, double
tension, the wrong scalar third derivative or sign, omitted pullback
mixing, and an omitted inverse-metric cross term. Completed run counts
are recorded separately. These are algebra checks, not a numerical
shell-contact evolution or an independent PDE existence theorem.

Here the shell-evaluation label `b` in equation (17) identifies the inherited shell; the material proper-normal coordinate is `r=0` there. It is not the ambient radial coordinate.

## Conventions and geometric references

The action, registered tension and residual signs are those of the [shell-Noether manuscript](https://github.com/maldonado-research/HDblast/blob/main/research/HDBLAST_SHELL_NOETHER_20261013/manuscript/SHELL_NOETHER.md). Moving-hypersurface variation formulas are discussed in [Mars, gr-qc/0507005](https://arxiv.org/abs/gr-qc/0507005); the material-chart identities here follow by differentiating the displayed exact junctions and do not import a vacuum boundary existence theorem.
