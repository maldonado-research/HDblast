# Forced degree-four scalar shell boundary system

14 October 2026. Conditional local mathematics for the registered five-dimensional
Einstein–scalar background and one common shell tension. Internal AI-assisted
derivation and independent checks; no external peer review or novelty claim.

## 1. Scope, convention and fields

This note supplies the complete scalar harmonic **physical boundary rows** at
second order for the selected copy-even first perturbation. It also gives their
Regge–Wheeler reduction, a useful symmetrizable Robin pair, and the precise
principal-symbol gap between writing junction equations and proving a coupled
initial-boundary-value theorem. The bulk operators are independently derived in
the companion bulk notes; their angular trace-free algebraic row is used explicitly
only in section 6. A forced tensor construction does not supply this scalar proof.

Write

\[
g_0=dy^2+\rho(y)^2\gamma,\quad
\gamma=-dt^2+c(t)^2\kappa_{ij}d\theta^id\theta^j,
\quad c=\cosh t,\quad h=\tanh t,
\]

\[
a=\rho_y/\rho,\quad p=\phi_{0,y},\quad
a_y=-\rho^{-2}-p^2/3.
\]

Increasing \(y\) is outward on **both** copies. At the shell \(y=y_b\),
\(a=\sigma_0/6\), \(p=-\sigma'_0/2\). All tension derivatives below are
evaluated at the background scalar there. There is one common tension, not one
tension for each interior. Use derivative coefficients

\[
g=g_0+\varepsilon h_1+\tfrac12\varepsilon^2h_2,
\qquad\phi=\phi_0+\varepsilon\chi_1+\tfrac12\varepsilon^2\chi_2.
\tag{1}
\]

The first representative is longitudinal, with zero first shell displacement
and fixed intrinsic identification:

\[
(h_1)_{yy}=-4uY,\quad(h_1)_{ya}=0,\quad
(h_1)_{ab}=2uY(h_0)_{ab},\quad\chi_1=vY,
\quad Y=n_1^2-n_2^2.
\tag{2}
\]

For \(J=n_1^2+n_2^2\), define
\(F_2=J-1/2\) and \(Q=F_4=Y^2-J/2+1/12\). On the unit three-sphere,

\[
Y^2=Q+F_2/2+1/6,\qquad \Delta_\kappa Q=-24Q.
\tag{3}
\]

Thus projection of a scalar product proportional to \(Y^2\) onto this precise
\(Q\) has coefficient one. No normalization or factorial is silently changed.
The formulas below hold for any \(l\ge2\), with \(\lambda=l(l+2)\), when its
quadratic source is projected consistently. Here \(\lambda=24\).

Keep every scalar second metric amplitude before choosing its gauge:

\[
\begin{array}{lll}
h_{2,yy}=NQ,&h_{2,yt}=BQ,&h_{2,yi}=C D_iQ,\\
h_{2,tt}=\rho^2TQ,&h_{2,ti}=\rho^2J_2D_iQ,&
h_{2,ij}=\rho^2c^2(KQ\kappa_{ij}+G S_{ij}),
\end{array}
\]

\[
S_{ij}=D_iD_jQ+(\lambda/3)\kappa_{ij}Q,\quad
\chi_2=wQ,\quad\zeta_2=ZQ.
\tag{4}
\]

The function \(J_2\) is a mixed metric amplitude, distinct from the angular
polynomial \(J\). The embedding is \(y=y_b+\varepsilon^2ZQ/2\). Intrinsic
second tangential identification is retained as a gauge freedom and described
below. Fields and radial derivatives in the displayed boundary equations are
evaluated at \(y_b\).

## 2. Exact quadratic sources in the original first gauge

For the first linear parameter path, the lapse is
\(\mathcal N^2=1-4\varepsilon uY\), the shift vanishes, and the tangential
metric is \((1+2\varepsilon uY)h_0\). Direct differentiation of
\(K_{ab}=\partial_yh_{ab}/(2\mathcal N)\) and
\(p(\varepsilon)=\partial_y\phi/\mathcal N\) gives

\[
K_{1,ab}=(u_y+4au)Yh_{0,ab},\quad
D^2K[h_1,h_1]_{ab}=(4uu_y+20au^2)Y^2h_{0,ab},
\]

\[
p_1=(v_y+2pu)Y,\qquad
D^2p[h_1,\chi_1]= (4uv_y+12pu^2)Y^2.
\tag{5}
\]

These terms include the unit normal normalization. Replacing the normal by
\(\partial_y\) for the perturbed metric loses them. The actual first shell rows
are

\[
u_y+2au+pv/3=0,\qquad
v_y+2pu+\tfrac12\sigma''v=0.
\tag{6}
\]

The exact copy-even junctions are \(K_{ab}=\sigma h_{ab}/6\) and
\(p=-\sigma'/2\). In a family-adapted material Gaussian chart their factorial
second coefficients obey

\[
K_2-ah_2-\frac{\sigma'}6\psi_2h_0
=\left(\frac{2\sigma'}3uv+\frac{\sigma''}6v^2\right)Y^2h_0,
\quad p_2+\frac{\sigma''}2\psi_2=-2v^2Y^2.
\tag{7}
\]

The scalar right side uses the registered \(\sigma'''=4\). Equation (7) has
material coefficient fields. Returning to the Eulerian fields (4), moving the
quadratic normal terms (5) to the right, and applying (6), gives

\[
\boxed{S_K=\frac{\sigma''}6v^2-12au^2,\qquad
S_p=-4pu^2+2\sigma''uv-2v^2.}
\tag{8}
\]

Neither \(S_K\) nor \(S_p\) is an additional physical energy or force. These are
inhomogeneous boundary Taylor coefficients in the stated Eulerian first gauge.

## 3. General geometric and component boundary rows

For a second perturbation with arbitrary normal lapse and shift components,
the linear variation of lower-index extrinsic curvature is

\[
DK[h_2]_{ab}=\tfrac12\partial_yh_{2,ab}
-D_{(a}h_{2,b)y}-\tfrac12a h_{2,yy}h_{0,ab}.
\tag{9}
\]

A second normal graph contributes
\(Z\partial_yK_{0,ab}-D_aD_b(ZQ)\) to \(K_2\),
\(2aZQh_{0,ab}\) to the induced metric, and \(pZQ\) to the scalar trace.
Using \(a_y-\sigma'p/6=-\rho^{-2}\), the complete metric row is

\[
\boxed{\tfrac12\partial_yh_{2,ab}-ah_{2,ab}
-D_{(a}h_{2,b)y}-\tfrac12a h_{2,yy}h_{0,ab}
-\tfrac16\sigma'\chi_2h_{0,ab}
-(D_aD_b+\rho^{-2}h_{0,ab})(ZQ)
=S_KQh_{0,ab}.}
\tag{10}
\]

The scalar row is

\[
\boxed{w_y-\tfrac12pN+\tfrac12\sigma''w
+(p_y+\tfrac12\sigma''p)Z=S_p.}
\tag{11}
\]

Equations (10)–(11) include graph evaluation, bending Hessian, lapse, shift,
and one-tension source normalization. The intrinsic tangential identification
does not add a physical source: its Lie derivatives cancel between \(K_2\)
and \(ah_2\), because \(K_0=ah_0\), and the background shell scalar is constant.

With the amplitude units (4), all independent scalar rows are explicitly

\[
\begin{split}
\tfrac12T_y-\rho^{-2}B_t+\tfrac12aN+\tfrac16\sigma'w
-\rho^{-2}(Z_{tt}-Z)&=-S_K,\\
\tfrac12(J_2)_y-\frac{C_t+B-2hC}{2\rho^2}
-\frac{Z_t-hZ}{\rho^2}&=0,\\
\tfrac12G_y-\frac{C+Z}{\rho^2c^2}&=0,\\
\tfrac12K_y+\frac{\lambda(C+Z)}{3\rho^2c^2}
+\frac{h(B+Z_t)}{\rho^2}-\tfrac12aN-\tfrac16\sigma'w
-\frac{Z}{\rho^2}&=S_K,\\
w_y-\tfrac12pN+\tfrac12\sigma''w
+(p_y+\tfrac12\sigma''p)Z&=S_p.
\end{split}\tag{12}
\]

The de Sitter covariant friction terms \(hC\), \(hZ\), and \(hZ_t\) are
essential. An angular trace-free scalar Hessian is distinct from a spatial
transverse tensor; it is retained here because it fixes scalar bending.

## 4. Second gauge and Regge–Wheeler boundary reduction

Let the second-only scalar diffeomorphism generator be
\(\xi^y=\alpha Q\), \(\xi^t=\beta Q\),
\(\xi^i=\ell\kappa^{ij}D_jQ\). With the active convention
\(h_2\mapsto h_2+\mathcal L_\xi g_0\), the embedding changes
\(Z\mapsto Z-\alpha_b\), and

\[
\begin{array}{ll}
\delta N=2\alpha_y,&\delta B=\alpha_t-\rho^2\beta_y,\\
\delta C=\alpha+\rho^2c^2\ell_y,&
\delta T=-2a\alpha-2\beta_t,\\
\delta J_2=c^2\ell_t-\beta,&
\delta K=2a\alpha+2h\beta-2\lambda\ell/3,\\
\delta G=2\ell,&\delta w=p\alpha.
\end{array}\tag{13}
\]

The verifier checks that every row (12) is unchanged under (13), including
the graph terms and background equations. A second-only gauge transformation
leaves the fixed first perturbation and its bulk quadratic source unchanged;
\(DE[\mathcal L_\xi g_0]=0\) on an exact background.

For \(l\ge2\), use \(\ell=-G/2\),
\(\beta=J_2+c^2\ell_t\), and
\(\alpha=-C-\rho^2c^2\ell_y\) to set

\[
G=C=J_2=0.
\tag{14}
\]

This is the scalar Regge–Wheeler gauge. The angular trace-free row of (12)
now gives \(Z=0\), and the time-angular row gives \(B_b=0\). For smooth
traces, \(B_t|_b=0\) follows. The complete physical boundary rows reduce to

\[
\boxed{Z=0,\quad B_b=0,\quad
K_y=aN+\frac{\sigma'}3w+2S_K,\quad
T_y=-K_y,\quad
w_y-\frac p2N+\frac{\sigma''}2w=S_p.}
\tag{15}
\]

Thus the second graph is straight in this gauge as a consequence of the
trace-free physical equation, rather than a discarded embedding equation.
These are boundary conditions on the rescaled metric amplitudes. Setting
\(\partial_yh_{2,ij}=0\) instead would omit the warp derivative.

## 5. Material conversion and bulk source conversion

Straightening first fields changes the first perturbation away from the shell.
Let \(V\) generate the first collar pullback and \(W\) its generator derivative.
For any field \(q\),

\[
q_1^{\rm mat}=q_1+\mathcal L_Vq_0,\qquad
q_2^{\rm mat}=q_2+2\mathcal L_Vq_1
+\mathcal L_V^2q_0+\mathcal L_Wq_0.
\tag{16}
\]

For (2) with fixed first shell identification, \(V_b=0\),
\(V^y_y|_b=2uY\), and \(V^a_y|_b=0\). Defining
\(A_g=2\mathcal L_Vh_1+\mathcal L_V^2g_0\) and
\(A_\phi=2V\chi_1+V(V\phi_0)\), their boundary jets are

\[
(A_g)_{ab}|_b=0,\quad
\partial_y(A_g)_{ab}|_b=(8uu_y+24au^2)Y^2h_{0,ab},
\quad (A_g)_{yy}|_b=-16u^2Y^2,
\]

\[
A_\phi|_b=0,\qquad
\partial_yA_\phi|_b=(4uv_y+4pu^2)Y^2.
\tag{17}
\]

The second Gaussian condition imposes

\[
W^y_y|_b=-h_{2,yy}/2+8u^2Y^2,\quad W^y_b=ZQ,
\quad W^a_y|_b=-h_0^{ab}(h_{2,yb}+D_b(ZQ)).
\tag{18}
\]

The coordinate acceleration derivative equals
\(-h_{2,yy}/2+12u^2Y^2\), since it also includes
\(\partial_y(V^BV^y_{,B})=4u^2Y^2\). Confusing it with \(W^y_y\)
loses a quadratic term. An independent coordinate chain-rule check recovers
the complete terms (5) from (17)–(18).

If the original full coupled second equation is
\(Lq_2=J^{\rm long}=-D^2\mathcal E[q_1,q_1]\), where \(L=D\mathcal E\)
acts on metric and scalar together, then

\[
\boxed{J^{\rm mat}=J^{\rm long}+L(A_g,A_\phi).}
\tag{19}
\]

The linear pure-gauge pair generated by \(W\) has zero \(L\)-image on the
exact background. Equation (19) is the full bulk source conversion, not only
an angular projection or a boundary correction. A material boundary equation
may therefore be paired with its material bulk source, or (10)–(15) with the
original longitudinal bulk source. Interchanging them without (19) is invalid.

## 6. A concrete two-variable Robin matrix

Both independent bulk derivations obtain the angular trace-free equation

\[
N-T+K=\mathcal C,\qquad\mathcal C=16u^2-\tfrac43v^2.
\tag{20}
\]

This is used here as a cross-checked companion result; the boundary verifier
checks the join algebra rather than re-deriving its bulk geometry. Combining
\(T_y=-K_y\) with (20) gives \(N_y+2K_y=\mathcal C_y\). Equations (15)
therefore imply

\[
\boxed{
\partial_y\begin{pmatrix}N\\w\end{pmatrix}
+\begin{pmatrix}2a&-4p/3\\-p/2&\sigma''/2\end{pmatrix}
\begin{pmatrix}N\\w\end{pmatrix}
=\begin{pmatrix}E_N\\S_p\end{pmatrix},}
\tag{21}
\]

\[
E_N=\mathcal C_y-4S_K
=-16au^2-\tfrac{16}3puv+\tfrac23\sigma''v^2.
\tag{22}
\]

The positive matrix \(D=\operatorname{diag}(3,8)\) symmetrizes this Robin
matrix exactly:

\[
DM=\begin{pmatrix}6a&-4p\\-4p&4\sigma''\end{pmatrix},
\qquad\det(DM)=24a\sigma''-16p^2.
\tag{23}
\]

Positive weights do not imply \(DM\) is positive definite. More importantly,
\(N\) and \(w\) have not thereby been shown to be the two unconstrained
hyperbolic scalar master channels. The remaining bulk constraint equations
must be eliminated or propagated, with their actual principal kinetic matrix.
An energy estimate for (21) alone is not an Einstein–scalar existence theorem.

## 7. Incoming modes, gauge rows and the remaining analytic gap

The separate `PRINCIPAL_SYMBOL_REVIEW.md` gives a frozen generalized-harmonic
calculation. In five-dimensional wave reduction there are fifteen metric
incoming amplitudes and one scalar amplitude. The ten Israel rows and scalar
Robin row supply physical normal boundary information. Appending all five
harmonic constraints to the physical rows adds only one independent metric
row on decaying wave modes; it leaves four tangential metric gauge modes.
For the scalar harmonic block, seven metric amplitudes plus \(w\) remain;
four physical metric rows, the scalar row, and all scalar harmonic constraints
have total rank six of eight. Raw equation counting misses two tangential
scalar gauge freedoms.

A concrete **frozen principal candidate** supplies four boundary mixed-metric
Dirichlet coordinate rows, \(g_{na}|_b=0\), and one normal harmonic constraint
row \(C_n=0\), alongside the physical rows. In the scalar block these are two
mixed Dirichlet rows and one normal constraint row. The stable symbol is full
rank for strictly \(\operatorname{Re}s>0\). Its determinant tends to zero at
glancing, as for a Neumann wave problem; one cannot claim a strict uniform
unweighted boundary-map inverse from that determinant alone.

At the frozen principal level, Codazzi plus the reduced equations give
homogeneous Neumann data for tangential subsidiary constraints, while the
normal constraint has homogeneous Dirichlet data. The exact nonconstant
umbilic matter cancellation is

\[
D_bK^b{}_a-D_aK-pD_a\phi
=-3D_a(\sigma/6)-pD_a\phi
=(-\sigma'/2-p)D_a\phi=0.
\tag{24}
\]

This supports the candidate constraint boundary mechanism, provided the
scalar junction and matching hold; it is not a replacement for deriving the
full curved subsidiary boundary equations and their inhomogeneous sources.

Fournodavlos–Smulevici, *The initial boundary value problem in General
Relativity: the umbilic case*, arXiv:2104.08851v1, supplies a closely related
boundary-orthogonal harmonic strategy and constraint mechanism for four-
dimensional vacuum with constant umbilicity. Its dimension, matter content,
and boundary coefficient differ. Their absolute harmonic normal condition
would give \(\partial_ng_{nn}=2d\lambda g_{nn}^{3/2}\), with intrinsic
dimension \(d\); here \(d=4\), hence coefficient eight rather than six.
Our proper-normal background has \(\Gamma_{0,y}=-4a\), so one must use a
compatible relative wave gauge/source or reparameterized harmonic background;
imposing \(\Gamma_y=0\) directly in the original Gaussian background is
inconsistent. Any first-gauge change also requires the full source conversion
(19). The cited theorem has not been transferred to this model.

A coupled contact theorem must construct physical wave channels or a
constrained wave system with a proved estimate, carry the lower-order
\(\sigma(\phi)\) couplings and sources through that gauge, propagate the
curved subsidiary constraints, and handle boundary traces and corner
compatibility. For a fixed finite angular representation set, the all-angular
glancing limitation is not an obstruction: a bounded matrix Robin law can be
solved by a characteristic boundary trace ODE and Volterra construction.
The accompanying `GH_BOUNDARY_SOURCE_AND_CONSTRAINT_COMPLETION.md` supplies
the exact curved boundary laws, subsidiary Robin equations, and a finite
SO(4) projection of localized initial constraint corrections for this route.
Its bulk/characteristic existence proof is a separate companion result.
The scalar Robin symmetrizer (23) alone does not complete those analytic steps.

## 8. Verification and scientific claim

`verify_scalar_boundary.py` independently differentiates the exact lapse
normalization, checks the one-shell factorial sources and material coordinate
jets, verifies the degree-four harmonic decomposition, tests every component
row under arbitrary second scalar gauge functions, and checks the Robin
join/symmetrizer. It imports no prior scientific code, uses explicit exceptions,
and has an actual `--inject-error` mode that changes the extrinsic quadratic
coefficient twenty to eighteen. Normal and optimized runs and deliberate
corruption runs are saved separately. The principal checker verifies frozen
rank and constraint boundary algebra separately.

This advances the explicit second-order scalar boundary calculation and
identifies a viable principal gauge structure. This note supplies boundary
data rather than a complete coupled evolution proof; the status and domain
of a finite-angular coefficient construction are stated in its companion
theorem note. A second-order coefficient construction does not establish an
exact nonlinear family through contact. This is conditional mathematics,
not evidence for a higher-dimensional origin of the Big Bang, nor an
observational or publication-readiness claim.
