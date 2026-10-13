---
title: Covariant flux and physical-domain limits of regular de Sitter shells
author: Ricardo Maldonado
affiliation: Independent Researcher
orcid: 0009-0009-3937-6527
email: ricardomaldonadoandrade@gmail.com
date: 13 October 2026 — technical draft for author review
abstract: |
  A positive self-adjoint radial operator does not by itself identify the
  physical perturbation space at a characteristic cone. For a registered
  doubled Einstein–scalar de Sitter shell, we derive the covariant
  presymplectic current on its generic longitudinal scalar solution slice.
  Its radial component is proportional to the reduced Green current at a
  common spectral value, but has a different mixed-spectral expression for
  general pairs and packets. We verify the even shell cancellation and the
  completed moving-shell current in the odd sector. Under explicit uniform
  characteristic asymptotics for actual solutions, the cone limit is a
  generally nonzero scalar null symplectic bilinear. This conditional
  statement is distinguished from an exact smooth finite-energy test-scalar
  example in a flat cone. The known singularity of a normalizable tensor
  zero mode is recovered componentwise, while cutoff estimates retain that
  mode in the energy completion. Smooth residual translations have no
  shrinking-vertex charge, without implying zero charge at finite null cuts.
  These calculations refine the scope of the existing radial theorem; they
  neither construct global radiative solutions of the registered shell
  problem nor establish physical completeness or nonlinear stability.
---

**Status.** Technical continuation draft; not an APS submission. The earlier
spectral capsule [1] is an immutable input, not a new edition of this note.

# 1. Problem, model and inherited result

The regular-shell checkpoint [1] proves positive scalar spectral bounds and
the tensor mass threshold for specified radial realizations. Its next
physical question is which covariant perturbations, boundary data and gauge
quotient those realizations describe. Three distinctions are essential:
smooth configurations versus their finite-energy completion; a shrinking
vertex limit versus a finite segment of the null cone; and an antisymmetric
symplectic flux versus a positive quadratic energy. We keep all three
explicit. Our main result is the local current and its shell balance for
the generic scalar sector of the registered model.

Take signature $(-++++)$, $8\pi G_5=1$, and two identical interiors with one
common shell. The action, including its variational boundary completion, is
\[
 S=\sum_{i=\pm}\left\{\int_{M_i}\!\epsilon_g
 \left[\frac R2-\frac12(\nabla\phi)^2-U(\phi)\right]
 +\int_\Sigma\!\epsilon_h K_i\right\}
 -\int_\Sigma\!\epsilon_h\sigma(\phi).                 \tag{1}
\]
Both shell normals point toward increasing interior radius. With
$\eta=\phi-1$, the registered functions are
\[
 W=\frac13+\eta^2+\frac{\eta^3}{3},\qquad
 U=\frac12W_{,\eta}^2-\frac23W^2,\qquad
 \sigma=2W+\delta(1+c+c\eta),\quad
 c=\frac2{1.0357712571566784}-\frac43 .                 \tag{2}
\]
The decimal denotes an exact rational. The existing regular branch covers
$0<\delta\le0.002$ and $|\eta|\le1/100$. We do not recertify that branch.
Write
\[
 g=dy^2+\rho(y)^2\gamma_{\mu\nu}dx^\mu dx^\nu,
 \quad \operatorname{Ric}(\gamma)=3\gamma,
 \quad 0<y\le T_b,\quad a=\rho'/\rho,\quad p=\phi'.
\]
The equations and background junctions are
\[
 p'=U_{,\phi}-4ap,\quad a'=-\rho^{-2}-p^2/3,\quad
 a^2=\rho^{-2}+p^2/12-U/6,\qquad
 a_b=\sigma_b/6,\quad p_b=-\sigma_{,\phi,b}/2.          \tag{3}
\]
At a fixed background, $p<0$ on $(0,T_b]$, $a_b>0$, and
$B_b=(p'/p+\sigma_{,\phi\phi}/2)_b>0$. The regular expansions have
$\rho=y+O(y^3)$ and $p=\kappa y+O(y^3)$, $\kappa\ne0$, with the
differentiated remainders used below. Smooth geometric extension statements
use coefficient functions smooth in $y^2$, as in the analytic regular branch.

For orientation, the inherited generic master coefficients are
\[
 P=\rho^2\Psi,\quad \alpha=(\rho^2p^2)^{-1},\quad
 q=\frac2{3\rho^2},\quad w=(\rho^4p^2)^{-1},\quad
 -(\alpha P')'+qP=\Lambda wP,\quad \mu^2=\Lambda-4.
                                                               \tag{4}
\]
The even endpoint condition is $\alpha_bP'_b=\Lambda P_b/\beta$ with
$\beta=\rho_b^4p_b^2B_b>0$; the odd one is
$\alpha_bP'_b=-c_oP_b$, $c_o=(3a_b\rho_b^2)^{-1}>0$.
The even Hilbert norm includes $|P_b|^2/\beta$, whereas the odd norm
contains only $\int w|P|^2dy$. The checkpoint establishes a strictly
positive scalar $\mu^2$ lower bound at each fixed detuning, essential
threshold $9/4$, and the even tensor spectral set
$\{0\}\cup[9/4,\infty)$, with $[9/4,\infty)$ for odd tensors.
These remain statements about their specified operators. No uniform
positive physical mass floor follows as $H_b=1/\rho_b\to0$.

# 2. Action-derived current on the generic scalar slice

For $h_{AB}=\delta g_{AB}$, $h=g^{AB}h_{AB}$ and $\chi=\delta\phi$,
fix the bulk Lee–Wald potential by
\[
 \theta^A=\sqrt{-g}\left[\frac12(\nabla_Bh^{AB}-\nabla^Ah)
                         -\nabla^A\phi\,\chi\right],\qquad
 \omega^A(1,2)=\delta_1\theta^A(2)-\delta_2\theta^A(1). \tag{5}
\]
Here $h^{AB}$ raises the covariant variation, rather than denoting
$\delta g^{AB}$. Variations commute. The displayed currents are coordinate
vector densities. In particular $\omega^y$ points toward increasing $y$;
the outward flux from $\epsilon<y<T_b$ at its lower boundary is
$-\omega^y$. Complex solutions may be paired by conjugating one argument
after complexification.

The formalism and its boundary ambiguities are established [2–5]. Adding a
boundary Lagrangian changes the potential by a field variation; a spacetime
exact improvement changes the current by a corner term. Equation (5), the
two GHY terms, and the one shell action specify our convention. We add no
cone counterterm and do not infer a null-boundary phase space from (5).

An explicit derivation uses the commuting nonlinear family
$g_{yy}=e^{2n}$, $g_{y\mu}=0$, $g_{\mu\nu}=\rho^2e^{2\psi}\gamma_{\mu\nu}$,
$\phi=\phi_0+X$. This family evaluates the second variation; it is not
assumed to satisfy nonlinear equations. With $K=4(a+\psi')e^{-n}$,
direct contraction in (5) gives
\[
 \theta^y+\delta(\sqrt{-h}K)
 =3\sqrt{-h}K\,\delta\psi-\sqrt{-h}e^{-n}\phi'\delta\phi
 =\pi^{\mu\nu}\delta h_{\mu\nu}+\pi_\phi\delta\phi,
                                                               \tag{6}
\]
where $\pi^{\mu\nu}=\sqrt{-h}(Kh^{\mu\nu}-K^{\mu\nu})/2$ and
$\pi_\phi=-\sqrt{-h}e^{-n}\phi'$. The field-space curl of the GHY
variation vanishes. Thus the radial canonical and covariant currents agree
pointwise on this zero-shift slice. General variations with shift still
require their ordinary tangential corner term.

In the generic scalar subspace the longitudinal constraints give
\[
 N=-2\Psi,\qquad
 \Psi'=-2a\Psi-pX/3,\qquad
 pX'=p'X+3\rho^{-2}(\Box_\gamma+4)\Psi-2p^2\Psi.        \tag{7}
\]
Here $N$ is the linear lapse perturbation. These implications require
independent scalar gradient and trace-free Hessian tensors. The constant
and Hessian-degenerate harmonic subspaces require separate treatment;
entire mass values are not excluded.

For a direct check, restrict the nonlinear family to $n=-2\psi$. Factoring
out $\sqrt{-\gamma}$, its potentials are
\[
 \theta^y=-\rho^4e^{6\psi}
 [4\delta\psi'+12(a+\psi')\delta\psi+(p+X')\delta X],
 \quad
 \theta^\mu=-\rho^2[D^\mu\delta\psi
                  +6\delta\psi D^\mu\psi+D^\mu X\delta X].
\]
Their antisymmetric variations give a radial gravity term
$-12\rho^4(\Psi_1\Psi'_2-\Psi_2\Psi'_1)$ and a scalar–metric term
$-6\rho^4p(\Psi_1X_2-\Psi_2X_1)$. Applying the momentum constraint in
(7) combines them into the following result, per interior:
\[
 \boxed{\omega^y=\sqrt{-\gamma}\rho^4
 [X_1X'_2-X_2X'_1-2p(\Psi_1X_2-\Psi_2X_1)]},           \tag{8}
\]
\[
 \boxed{\omega^\mu=\sqrt{-\gamma}\rho^2
 [X_1D^\mu X_2-X_2D^\mu X_1
  +6(\Psi_1D^\mu\Psi_2-\Psi_2D^\mu\Psi_1)]}.          \tag{9}
\]
The $-2p$ coefficient needs both gravity and scalar variations. Equation
(9) also agrees with the local constrained kinetic action. Neither formula
is guessed from a master-variable norm.

Put $\mathcal D=\Box_\gamma+4$. The reconstruction
$X=-3P'/(\rho^2p)$ and (7) convert (8) to
\[
 \boxed{\omega^y=9\sqrt{-\gamma}\alpha
 [(\mathcal DP_1)P'_2-(\mathcal DP_2)P'_1]}.             \tag{10}
\]
For $P_i(y)Y_i(x)$ with $\Box Y_i=(\Lambda_i-4)Y_i$, this is
\[
 \begin{split}
 \omega^y&=9\sqrt{-\gamma}\alpha
   [\Lambda_1P_1P'_2-\Lambda_2P'_1P_2]Y_1Y_2,\\
 \omega^\mu&=9\sqrt{-\gamma}
   [\alpha P'_1P'_2+qP_1P_2]
   (Y_1D^\mu Y_2-Y_2D^\mu Y_1).
 \end{split}                                         \tag{11}
\]
Consequently the multiplier is $9\Lambda$ times the reduced Green current
only at a common spectral value. A packet requires (10) or the mixed terms
in (11). The independent identity
\[
 \partial_y\{\alpha[\Lambda_1P_1P'_2-\Lambda_2P'_1P_2]\}
 = (\Lambda_1-\Lambda_2)(\alpha P'_1P'_2+qP_1P_2)       \tag{12}
\]
follows from (4) and cancels the tangential divergence in (11). A
common-$\Lambda$ replacement for unequal profiles fails this check.

# 3. Doubled-shell balance

Varying the common induced fields in (1) gives
\[
 \sum_i(K^{(i)}_{\mu\nu}-K_i h_{\mu\nu})+\sigma h_{\mu\nu}=0,
 \qquad \sum_i n_i\!\cdot\!\nabla\phi+\sigma_{,\phi}=0. \tag{13}
\]
Two pure-even or two pure-odd perturbations with the same one-side amplitude
double (8)–(11); cross-parity bulk pairings cancel. In the generic even
sector the shell is straight, and its condition in (4) cancels
$\omega^y_+(T_b)+\omega^y_-(T_b)$ even for unequal $\Lambda_i$.
There is no extra tangential kinetic shell current on this slice. The
endpoint mass in the even reduced norm is not an independent physical shell
kinetic field in this gauge.

The odd shell moves in longitudinal gauge. The complete moving-boundary
action in [1] retains graph GHY minus slice GHY, the Einstein tangential
divergence, the shifted bulk integration domain, and one common shell
volume. In a graph $y=T_b+\zeta(x)$, the first two terms contribute
$+3a_b\rho_b^2(D\zeta)^2$ per copy, while the shifted domain contributes
$+3\rho_b^2D\zeta\cdot D\psi_b$ per copy. The single shell contributes
$-\sigma_b\rho_b^2(D\zeta)^2/2$. Using odd continuity
$\psi_b=-a_b\zeta$, $X_b=-p_b\zeta$ and (3), their sum is
$-3a_b\rho_b^2(D\zeta)^2$.

With $\zeta=-P_b/(a_b\rho_b^2)$, this fixes the completed kinetic action
and its current:
\[
 S^{(2)}_{\Sigma,\mathrm{kin}}=-9c_o\int\sqrt{-\gamma}(DP_b)^2,
 \qquad
 \omega^\mu_\Sigma=18c_o\sqrt{-\gamma}
 (P_{1b}D^\mu P_{2b}-P_{2b}D^\mu P_{1b}).              \tag{14}
\]
It includes both gravitational copies and one shell, rather than two shell
tensions. The odd condition in (4) gives the exact balance
\[
 \partial_\mu\omega^\mu_\Sigma
 =\omega^y_+(T_b)+\omega^y_-(T_b).                     \tag{15}
\]
Indeed both sides equal
$18c_o\sqrt{-\gamma}(\Lambda_2-\Lambda_1)P_{1b}P_{2b}Y_1Y_2$
for separated pairs. Thus
\[
 \partial_\mu\left[\int_\epsilon^{T_b}
 (\omega^\mu_++\omega^\mu_-)dy+\omega^\mu_\Sigma\right]
 =\omega^y_+(\epsilon)+\omega^y_-(\epsilon).            \tag{16}
\]
This is the conserved bulk-plus-shell current, subject to the lower-boundary
flux. It is proved on the constrained generic class, not for arbitrary
unreduced embeddings or additional boundary degrees of freedom.

At a common spectral root the integrated kinetic coefficients reproduce
$K_e=18\Lambda N_e$ and $K_o=18\Lambda N_o$ from [1]. Their healthy sign
does not make (8), (14) or (16) positive bilinears: presymplectic currents
are antisymmetric and vanish when their two arguments coincide.

# 4. Vertex geometry and the radial endpoint

Use global $dS_4$ coordinates,
$\gamma=-dt^2+\cosh^2t\,d\Omega_3^2$. The exact conformal radius obeys
\[
 \frac{dr}{r}=\frac{dy}{\rho},\qquad r/y\to1,
 \qquad \rho=r\Omega(r^2),\qquad dy=\Omega(r^2)dr.
\]
With ambient coordinates
$\mathsf T=r\sinh t$, $\mathsf R=r\cosh t$ and
$u=re^{-t}$, $v=re^t$, the metric is exactly
\[
 g=\Omega(uv)^2\left[du\,dv+\frac{(u+v)^2}{4}d\Omega_3^2\right],
 \quad \mathsf T=(v-u)/2,\quad \mathsf R=(v+u)/2.        \tag{17}
\]
Here $\Omega(0)=1$. Smooth even regular expansions give local smooth
background extensions; they do not prescribe a global continuation.
The radial chart covers the spacelike exterior of the cone. At fixed finite
$t$, $y\to0$ reaches the vertex. Reaching $u=0$, $v>0$ requires
$t=\log(v/r)\to+\infty$. These are different limits, as already emphasized
in the de Sitter-brane geometry of [6].

For smooth perturbations with bounded Cartesian first derivatives, the
geometric current (5) is bounded near the vertex. The proper volume of
$y=\epsilon$ over a fixed compact time interval is $O(\epsilon^4)$.
Its flux therefore vanishes. This sufficient statement is independent of
the divided master variable and does not control a time interval translated
to infinity.

For the radial operator set
\[
 x=\int_y^{T_b}\frac{ds}{\rho(s)},\quad
 F=\frac{P}{\rho^{3/2}|p|},\quad
 s_L=\rho\left(\frac32a+\frac{p'}p\right).
\]
The potential is $V_L=s_{L,x}+s_L^2+2\rho^2p^2/3
=25/4+O(e^{-2x})$, and
\[
 \alpha(P_1P'_2-P'_1P_2)=-W_x(F_1,F_2).               \tag{18}
\]
For the operator domain, $F\in H^2$ near infinity, so both $F$ and $F_x$
tend to zero. This proves the radial Green-form limit. The form domain
$H^1$ alone gives no pointwise endpoint trace for $F_x$. Smooth disjoint
bumps can give two $H^1$ functions whose Wronskian is one at an unbounded
sequence and zero between supports; the closed variational identity does
not assert a pointwise limit for all form vectors.

For the physical packet current (10), let
$G=(\mathcal DP)/(\rho^{3/2}|p|)$. Then
\[
 \frac{\omega^y}{9\sqrt{-\gamma}}
 =-G_1F_{2x}+G_2F_{1x}+s_L(G_1F_2-G_2F_1).            \tag{19}
\]
A sufficient strong finite-time condition is continuity on a compact
time interval into $H_x^2(L^2(S^3))$ for $F_i$ and into
$H_x^1(L^2(S^3))$ for $G_i$, with enough angular regularity to define
the current. Compactness of these time-indexed sets yields uniform tail
control and zero integrated limit. No bound here is uniform over time
intervals tending to infinity. Weak energy-space flux requires a separate
trace or approximation construction.

\begin{figure}[tb]
\centering
\IfFileExists{../figures/cone-limits.pdf}{%
\includegraphics[width=\linewidth]{../figures/cone-limits.pdf}}{%
\fbox{\parbox{0.9\linewidth}{Optional figure pending: exact conformal cone
geometry and the exact flat test-scalar comparison.}}}
\caption{Distinct cone limits. The geometric panel uses conformal radius
$r$ and ambient $(\mathsf T,\mathsf R)$. Any numerical flux panel refers
only to the exact flat test-scalar example, not a registered gravitational
wave packet.}
\end{figure}

# 5. Characteristic flux: a conditional result and an exact comparator

**Conditional registered-model statement.** Let two actual generic
linearized solutions admit, uniformly on a compact segment
$0<v_1\le v\le v_2$, the differentiated asymptotics
\[
 P_i=y^4 A_i(v,\vartheta)+O(y^6),\qquad
 \Psi_i=y^2 A_i+O(y^4),\qquad X_i=X_{0i}+O(y^2),
 \quad X_{0i}=-\frac3\kappa(4A_i+v\partial_vA_i).       \tag{20}
\]
Here $v=re^t$ and $\vartheta\in S^3$; remainders are smooth in
$(u,v,\vartheta)$ with the uniform derivatives entering the current and
equations. A $C^1$ bound suffices for current limits, and the invoked
field equations require corresponding $C^2$ control. The last relation
follows from radial momentum reconstruction. Replacing $r$ by $y$ changes
only controlled subleading orders on this segment.

These leading fields are compatible with smooth metric components at the
finite null surface: since $\delta g=2\Psi g-6\Psi dy^2$, their leading
expression is
\[
 \delta g=2uvA\,g-\frac32A(v\,du+u\,dv)^2
          +\text{smooth subleading tensor coefficients}.        \tag{21}
\]
This does not prove existence, a smooth global longitudinal representative,
or smoothness at $u=v=0$. An $O(u^2)$ remainder in $\Psi$ can contribute
$O(u)$ to the $du^2$ coefficient and must not be mislabeled $O(u^2)$
as a tensor remainder.

Directly applying (8) gives, per interior,
\[
 \omega^y\,dt\longrightarrow
 \frac{v^3}{8}\sqrt{\gamma_{S^3}}
 (X_{01}\partial_vX_{02}-X_{02}\partial_vX_{01})\,dv.   \tag{22}
\]
Indeed $\rho^4\cosh^3t=(v^3/8)y+O(y^3)$ and
$X_i'|_t=(v/y)\partial_vX_{0i}+O(y)$; the scalar–metric term vanishes
in this limit. The density transformation verifies that this is the actual
null pullback. Its exact Jacobian in (17) gives
\[
 \omega^u=\frac{\Omega}{2r}
 \left(\frac{u}{\rho}\omega^y-u\omega^t\right)
 \longrightarrow\frac{v^3}{8}\sqrt{\gamma_{S^3}}
 (X_{01}\partial_vX_{02}-X_{02}\partial_vX_{01}).        \tag{23}
\]
Both components are necessary: omitting $\omega^t$ loses a factor of two.
Equivalently $\sqrt{-g}g^{uv}\to v^3\sqrt{\gamma_{S^3}}/8$.
The doubled parity flux in the $+y$ convention is consequently
\[
 \boxed{\lim_{y\to0}\int\omega^y_{\mathrm{doubled}}dt\,d\vartheta
 =\frac14\int_{v_1}^{v_2}v^3dv\int_{S^3}\!d\Omega_3
 (X_{01}\partial_vX_{02}-X_{02}\partial_vX_{01})}.       \tag{24}
\]
The outward flux from the exterior has the opposite sign. This is an
antisymmetric symplectic flux, not a positive energy flux.

The leading bilinear is not identically zero: choosing
$X_{01}=f(v)H(\vartheta)$, $X_{02}=vf(v)H(\vartheta)$ gives Wronskian
$f^2H^2$. Equation (20) can be solved for the corresponding leading $A_i$
on the positive segment. These are candidate leading data only. No pair of
global registered solutions with these data, the distant shell conditions
and the required gauge regularity has been constructed here.

**Unconditional flat comparator.** A massless test scalar in five-dimensional
Minkowski space has radial operator $-y^{-2}(y^4 f')'$, Hilbert weight
$y^2$, and form $\int y^4|f'|^2dy$. At zero mass its radial solutions
are $1$ and $y^{-3}$; only the first is square integrable at zero. The
endpoint is limit-point, with Liouville potential $9/4$. Consider
\[
 \varphi_1=1,\qquad \varphi_2=\mathsf T=y\sinh t.       \tag{25}
\]
These are smooth local exact solutions. To realize them with finite ambient
energy, choose compact smooth initial data at $\mathsf T=0$ equal to
$(\varphi_1,\dot\varphi_1)=(1,0)$ and
$(\varphi_2,\dot\varphi_2)=(0,1)$ on a ball of radius $R_0$.
Finite propagation speed gives (25) wherever $|\mathsf T|+\mathsf R<R_0$.
Thus any cone segment $0<a\le v\le b<R_0$ has a neighborhood where
the finite-energy solutions equal (25). A remote boundary outside that
causal diamond does not affect the calculation.

At each finite $t$ the Liouville profiles are $y^{3/2}$ and
$y^{5/2}\sinh t$, both $H^2$ near the endpoint. Their current per solid
angle on $y=\epsilon$ obeys
\[
 \int_{t_1}^{t_2}\!j^y dt
 =\frac{\epsilon^4}{4}[\cosh^4t_2-\cosh^4t_1]\longrightarrow0,
                                                               \tag{26}
\]
whereas on $t\in[\log(a/\epsilon),\log(b/\epsilon)]$,
\[
 j^y dt=\frac{(v+\epsilon^2/v)^3(v-\epsilon^2/v)}{16v}dv
 \longrightarrow\frac{v^3}{16}dv,\qquad
 \int j^y dt\longrightarrow\frac{b^4-a^4}{64}\ne0.     \tag{27}
\]
The convergence is uniform since $a>0$. The second field separately has
nonzero null derivative and a nonzero quadratic stress-energy component;
this is distinct from the antisymmetric pairing in (27).

The example disproves an inference from radial limit-point closure, smooth
extension and finite energy to zero characteristic flux. It is a test-field
comparator, not a solution of the registered coupled Einstein–scalar shell
problem. It contradicts no self-adjoint radial theorem. Free fields on
de Sitter leaves already admit established ambient spectral decompositions
under additional state and convergence hypotheses [7].

**Energy-completion comparator.** The same exact flat test field also
distinguishes two energy topologies. On $\mathsf T=t=0$, take zero velocity
and $f(y)=\chi_{\rm out}(y)y^{-\nu}$, where the outer cutoff is smooth
and equals one near zero. For $1\le\nu<3/2$,
\[
 N[f]\sim\int_0 y^{2-2\nu}dy<\infty,\qquad
 Q[f]\sim\nu^2\int_0 y^{2-2\nu}dy<\infty,\qquad
 E_{\mathsf T,\mathrm{grad}}[f]\sim
 \nu^2\int_0 y^{1-2\nu}dy=\infty.                    \tag{27a}
\]
For $\nu=5/4$ the powers are $-1/2$ and $-3/2$, respectively.
Inner smooth cutoffs converge in $N+Q$, with cost
$O(\epsilon^{3-2\nu})$, while their inertial gradient energies diverge.
The different weights reflect the lapse $y$ of the de Sitter time clock:
its radial gradient contribution is $\int y^4|f'|^2dy$, whereas an
inertial slice has $\int y^3|f'|^2dy$. This example uses a singular
energy-limit datum, in contrast to the smooth finite-energy solutions
in (25)–(27). It disproves a blanket identification of these two test-field
completions. It establishes no equivalence or inequivalence of registered
gravitational phase spaces. Specifying which energy is completed is
indispensable.

# 6. Tensor regularity and energy completion

The even tensor zero profile $h_T(y)=1$ has norm $\int_0^{T_b}\rho^2dy$.
Let $q_{\mu\nu}$ be a slice TT Einstein graviton,
$\delta R_{\mu\nu}[\gamma;q]=3q_{\mu\nu}$. Its lift is
$\delta g_{\mu\nu}=\rho^2q_{\mu\nu}$, $\delta\phi=0$.
For any four-dimensional Einstein slice, warped curvature gives
\[
 C^{(5)}_{\mu\nu\alpha\beta}=\rho^2C^{(4)}_{\mu\nu\alpha\beta},
 \qquad C^{(5)}_{y\mu y\nu}=C^{(5)}_{y\mu\nu\alpha}=0.
\]
Linearization and an orthonormal frame therefore give
\[
 \delta C^{(5)}_{\hat a\hat b\hat c\hat d}
 =\rho^{-2}\delta C^{(4)}_{\hat a\hat b\hat c\hat d}.   \tag{28}
\]
The background Weyl tensor vanishes, so this perturbation is gauge invariant.
Along a finite-$t$ ray the frame has a bounded Cartesian limit. Nonzero
slice Weyl implies a $y^{-2}$ divergence and excludes a $C^2$ metric
extension through the vertex. Tensor components, rather than the Lorentzian
contraction $C^2$, are needed because nonzero radiative Weyl tensors can
have vanishing scalar contractions. Directions with zero slice Weyl are
not excluded by this argument.

This diagnostic is established prior art: Garriga and Sasaki [6, section
III, p. 9] display $C_5^2=C_4^2/a^4$ and explain the singular normalizable
zero mode as part of a mode expansion. Equation (28) applies the mechanism
componentwise to the nonconstant scalar warp; it is not a novel discovery.

The smoothness obstruction does not eliminate the mode from an energy
completion. For a smooth cutoff $h_\epsilon=\chi(y/\epsilon)$, zero below
$\epsilon$ and one above $2\epsilon$, the Neumann shell condition is
unchanged and
\[
 \int\rho^2|h_\epsilon-1|^2dy=O(\epsilon^3),\qquad
 \int\rho^4|h_\epsilon'|^2dy=O(\epsilon^3).            \tag{29}
\]
These are initial-data approximants, not separated zero-mass eigenfunctions.
For global $S^3$ tensor harmonics, data on the regular slice
$\mathsf T=0$ have bounded Cartesian metric components, with
derivatives and velocities $O(r^{-1})$. Their local $H^1\times L^2$
energy tail is $\int_0^\epsilon r^3r^{-2}dr=O(\epsilon^2)$; the cutoff
derivative costs the same order. Radial cutoffs preserve angular TT metric
and extrinsic-curvature data. Their full spatial divergence remains zero
because its radial term is proportional to the angular trace. On the
time-symmetric background slice the angular Ricci tensor is isotropic,
the scalar gradient is radial, and $\delta\phi=\delta\pi_\phi=0$;
the linearized Hamiltonian and momentum constraints are preserved.

Thus pointwise smooth configurations and finite-energy limits differ here.
No deletion of massless gravity, inequivalence of completed physical spaces,
or global evolution theorem follows from (28)–(29).

# 7. Gauge charges and exceptional subspaces

For a field-independent vector $\xi$, the on-shell covariant identity is
\[
 \boldsymbol\omega(\delta,\mathcal L_\xi)
 =d\boldsymbol k_\xi(\delta),\quad
 \boldsymbol k_\xi=\delta\boldsymbol Q_\xi-i_\xi\boldsymbol\theta,
 \quad \boldsymbol Q_\xi=-\tfrac12\star d\xi^\flat.     \tag{30}
\]
Field-dependent vectors require the $-Q_{\delta\xi}$ correction. Smooth
Cartesian fields, perturbations and vector generators with bounded required
derivatives have bounded $k_\xi$. A cut $y=\epsilon,t=t_0$ has area
$O(\epsilon^3)$ at fixed finite $t_0$, hence
\[
 \int_{S^3_{\epsilon,t_0}}k_\xi=O(\epsilon^3)\to0.      \tag{31}
\]
This sufficient vertex-charge statement gives no result on a finite-radius
null cut. At the shell the embedding must transform with the bulk fields.
Generators supported away from physical boundaries are degeneracies; other
generators require the total boundary pairing in the specified phase space.

For precision, the Hessian-degenerate subspace is
$D_\mu D_\nu Y=-\gamma_{\mu\nu}Y$, locally the five embedding harmonics
$Y=c_An^A$ with $X^A=rn^A$. Its Gaussian-normal residual is
$\xi^y=e_0Y$, $\xi^\mu=(\ell_0-e_0F)D^\mu Y$, $F'=\rho^{-2}$.
Fix the additive constant by the regular odd Laurent expansion. With
$s=r^2$, $A=-rF$ and $B=(\Omega^{-1}-A)/s$, its translation part is
\[
 \xi^A=e_0[A(s)c^A+B(s)(c_BX^B)X^A],\qquad
 A-2sA_s=\Omega^{-1},\quad A(0)=1.                    \tag{32}
\]
Both coefficients are smooth. For example
$F=-1/r+\int_0^r[\Omega(v^2)^{-1}-1]v^{-2}dv$ exhibits this explicitly.
The apparent $1/y$ tangential coefficient represents a regular Cartesian
translation. In contrast, the remaining vector is
$\ell_0[rc^A-(c_BX^B)X^A/r]$, generically not continuous in the full
ambient approach to the vertex, and hence not $C^1$. For example, with
$c=(0,1,0,\ldots)$, $X^1=\epsilon$ and
$X^0=\sqrt{\epsilon^2-\epsilon^6}$, one has $r=\epsilon^3$ and
$\xi^1=\ell_0(\epsilon^3-\epsilon^{-1})$. Bounded finite-$t$ rays do
not detect this divergence.
A smooth gauge convention fixes that constant; no charged physical mode
is inferred from its exclusion. A weaker completed gauge group needs its
own pairing analysis.

The background equation also gives
$a+F=-\frac13\int_0^yp(v)^2dv=O(y^3)$. The earlier relative-displacement
exclusion uses $p_bB_b\ne0$ and simultaneous shell matching. The constant
harmonic instead uses the static Jacobi determinant in [1]. Neither
argument permits dropping all harmonics at $\mu^2=0$ or $\mu^2=-4$;
other harmonics at these values remain generic. We prove no global
exceptional-charge classification or global longitudinal-gauge coverage.

# 8. Interpretation, prior art and reproducibility

Equations (8)–(16) supply an explicit action-derived current and shell
balance for this constrained solution class, including unequal spectral
values. Equation (24) is a conditional identity for actual solutions with
the stated uniform asymptotics. Equation (27) is an unconditional
counterexample to a general zero-characteristic-flux inference. Equations
(28)–(32) distinguish smoothness, completion and local gauge charge.
All leave the earlier radial positivity theorem intact.

The underlying methods and several diagnostics predate this calculation.
Iyer–Wald [2] describes the potential and charge ambiguities; Wald–Zoupas [3]
requires a boundary field space and suitable extension of the full current
for radiative charges. Harlow–Wu [4] treats variational and corner terms,
explicitly excluding null boundaries from its main construction. Null
charges in [5] concern a specified four-dimensional vacuum field space,
not this expanding five-dimensional matter cone and shell. Garriga–Sasaki
[6] supplies the direct horizon and zero-mode precedent. Bertola et al. [7]
give an ambient spectral dictionary under additional hypotheses. Gen–Sasaki
[8] treats displacement and moving-boundary terms, while
Kofman–Martin–Peloso [9] derives physical normalization including boundary
terms.

The recent null-hypersurface analysis of Ciambelli–Freidel–Leigh [10] was
also inspected. Its selected sections explicitly retain boundary-kernel,
corner and function-space choices; its four-dimensional vacuum kinematical
reduction is not a completed constrained phase space for this model.
The accompanying literature report records bounded queries and primary
passages read. It establishes neither exhaustive coverage nor priority for
this exact application. The claim is the displayed model-specific
calculation, not a new phase-space method, threshold, or physical mechanism.

The four accompanying derivations and their exact checkers verify the
potential curl, mixed-spectral conservation, shell signs, Liouville
conversion, density Jacobian, comparator flux and tensor/gauge identities.
Deliberately wrong formulas and optimization-mode runs test the checkers.
Independent internal calculations agree on the current. These receipts
support algebra and exact scaling calculations; they are not formal proofs
of functional analysis, hyperbolic existence or completeness.

The next task is to specify a regular Cauchy-data completion or a
characteristic boundary-data space, prove existence and shell compatibility,
control the reconstructed current under that completion, and identify its
zero-charge gauge quotient. Compactly supported tangential test variations
used in local action integrations are not arbitrary nonzero on-shell
solutions compactly supported in time. Finite-slab evolution keeps its
Cauchy terms, and the global $S^3$ sections avoid unmentioned spatial flux.
No full registered radiative packet, global stability theorem, branch
transition, reheating process, or evidence for a higher-dimensional
cosmological origin is established here.

**Source and data availability.** The prior checkpoint is
DOI [10.5281/zenodo.23335485](https://doi.org/10.5281/zenodo.23335485),
*HDBLAST: Conditional Scalar Spectral Positivity and Tensor Mass Thresholds
on Regular de Sitter Shells*, with immutable source
[401efdfaaef657bcbdaec588b21b82789daddd0e](https://github.com/maldonado-research/HDblast/tree/401efdfaaef657bcbdaec588b21b82789daddd0e/research/HDBLAST_REGULAR_SHELL_SPECTRA_20261013).
This continuation and its technical notes, checkers and review records are
distributed in the accompanying
[GitHub capsule](https://github.com/maldonado-research/HDblast/tree/main/research/HDBLAST_CONE_PHASE_SPACE_20261013).
There is no separate DOI for this continuation. DOI
10.5281/zenodo.23335485 identifies the prior spectral checkpoint only;
that edition and its immutable source remain unchanged.

**Declarations.** No funding was received. The author declares no known
competing interests. Substantial OpenAI Codex assistance was used for
derivations, symbolic-check implementation, literature organization and
drafting. Separate agent processes supplied internal review; this is not
external peer review. The author retains responsibility for the work and
this technical draft remains subject to author review.

# References

1. R. Maldonado, *HDBLAST: Conditional Scalar Spectral Positivity and Tensor
   Mass Thresholds on Regular de Sitter Shells* (2026),
   [doi:10.5281/zenodo.23335485](https://doi.org/10.5281/zenodo.23335485).
2. V. Iyer and R. M. Wald, *Some Properties of Noether Charge and a Proposal
   for Dynamical Black Hole Entropy*, Phys. Rev. D **50**, 846 (1994),
   [doi:10.1103/PhysRevD.50.846](https://doi.org/10.1103/PhysRevD.50.846),
   [gr-qc/9403028v1](https://arxiv.org/abs/gr-qc/9403028v1).
3. R. M. Wald and A. Zoupas, *A General Definition of “Conserved Quantities”
   in General Relativity and Other Theories of Gravity*, Phys. Rev. D **61**,
   084027 (2000), [doi:10.1103/PhysRevD.61.084027](https://doi.org/10.1103/PhysRevD.61.084027),
   [gr-qc/9911095v1](https://arxiv.org/abs/gr-qc/9911095v1).
4. D. Harlow and J.-q. Wu, *Covariant phase space with boundaries*, JHEP
   **10** (2020) 146, [DOI](https://doi.org/10.1007/JHEP10(2020)146),
   [1906.08616v4](https://arxiv.org/abs/1906.08616v4).
5. V. Chandrasekaran, É. E. Flanagan and K. Prabhu, *Symmetries and charges
   of general relativity at null boundaries*, JHEP **11** (2018) 125,
   [doi:10.1007/JHEP11(2018)125](https://doi.org/10.1007/JHEP11(2018)125),
   [1807.11499v3](https://arxiv.org/abs/1807.11499v3).
6. J. Garriga and M. Sasaki, *Brane-world creation and black holes*, Phys.
   Rev. D **62**, 043523 (2000),
   [doi:10.1103/PhysRevD.62.043523](https://doi.org/10.1103/PhysRevD.62.043523),
   [hep-th/9912118v2](https://arxiv.org/abs/hep-th/9912118v2).
7. M. Bertola, J. Bros, V. Gorini, U. Moschella and R. Schaeffer,
   *Decomposing Quantum Fields on Branes*, Nucl. Phys. B **581**, 575–603 (2000),
   [doi:10.1016/S0550-3213(00)00280-7](https://doi.org/10.1016/S0550-3213(00)00280-7),
   [hep-th/0003098v1](https://arxiv.org/abs/hep-th/0003098v1).
8. U. Gen and M. Sasaki, *Radion on the de Sitter brane*, Prog. Theor. Phys.
   **105**, 591 (2001), [doi:10.1143/PTP.105.591](https://doi.org/10.1143/PTP.105.591),
   [gr-qc/0011078v2](https://arxiv.org/abs/gr-qc/0011078v2).
9. L. Kofman, J. Martin and M. Peloso, *Exact identification of the radion
   and its coupling to the observable sector*, Phys. Rev. D **70**, 085015
   (2004), [doi:10.1103/PhysRevD.70.085015](https://doi.org/10.1103/PhysRevD.70.085015),
   [hep-ph/0401189v1](https://arxiv.org/abs/hep-ph/0401189v1).
10. L. Ciambelli, L. Freidel and R. G. Leigh, *The Phase Space of Gravity on
    Null Hypersurfaces* (2026),
    [arXiv:2608.14449v1](https://arxiv.org/abs/2608.14449v1), submitted
    14 August 2026. No DOI is assigned here.
