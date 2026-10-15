---
title: Smooth Einstein–scalar perturbations through timelike shell contact
author: Ricardo Maldonado
affiliation: Independent Researcher
orcid: 0009-0009-3937-6527
email: ricardomaldonadoandrade@gmail.com
date: 13 October 2026 (America/Los_Angeles) — technical draft for author review
abstract: |
  We construct smooth linearized Einstein–scalar perturbations that evolve
  through contact with a timelike shell in a fixed five-dimensional model.
  A conformal map sends the hyperboloidal shell to the equator of a round
  hemisphere and converts its exact even boundary law into a dynamic
  boundary equation with positive auxiliary inertia; the odd law becomes
  Robin. For finitely many angular harmonics and arbitrary smooth radial
  collar data, a characteristic integral system coupled to a boundary
  ordinary differential equation proves existence, uniqueness and smoothness
  on an explicit finite slab. Time reversal and one-dimensional sideways
  uniqueness prove actual boundary activation for pure displacement or
  pure velocity data. Direct graph geometry then
  verifies all metric and scalar junctions, including the odd shell
  displacement. The action-derived covariant current requires a moving-shell
  cut term in odd parity. A radial Poincaré bound gives a nondegenerate
  physical Cauchy pairing on a sufficiently thin compact annulus and a
  nonzero total pairing conserved through contact. No nonzero local shell
  current, general completed phase space, nonlinear stability or cosmological
  origin claim is inferred.
---

**Status.** Technical continuation for author review. This is not an APS
submission or external peer review. Detailed derivations and executable
identity checks accompany the draft; verification scope and declarations
appear in Section 8.

# 1. Fixed model and the statement proved

Use signature $(-++++)$, $8\pi G_5=1$, and two interior copies joined at a
common timelike shell:
\[
 S=\sum_{i=\pm}\left\{\int_{M_i}\epsilon_g
 [R/2-(\nabla\phi)^2/2-U(\phi)]+\int_\Sigma\epsilon_h K_i\right\}
 -\int_\Sigma\epsilon_h\sigma(\phi).                       \tag{1}
\]
Both outward normals point toward increasing interior radius. The fixed
registered functions are
\[
 \begin{aligned}
 &\eta=\phi-1,\quad W=\tfrac13+\eta^2+\tfrac13\eta^3,\quad
 U=\tfrac12W_{,\phi}^2-\tfrac23W^2,\\
 &
 \sigma=2W+\delta(1+c_*+c_*\eta),
 \qquad c_*=\frac2{1.0357712571566784}-\frac43,\quad
 0<\delta\le0.002.
 \end{aligned}                                         \tag{2}
\]
The decimal is exact. Fix one background on the previously registered
nonconstant regular branch [1,2]. In its exterior-cone chart,
\[
 \begin{aligned}
 &g=dy^2+\rho^2\gamma,\quad \operatorname{Ric}\gamma=3\gamma,
 \quad a=\rho'/\rho,\quad p=\phi',\\
 &
 a'=-\rho^{-2}-p^2/3,\quad p'=U_{,\phi}-4ap,
 \quad a_b=\sigma_b/6>0,\quad p_b=-\sigma'_b/2\ne0.
 \end{aligned}           \tag{3}
\]
We assume the branch property
$B_b=p'_b/p_b+\sigma''_b/2>0$. This paper does not recertify the background
branch or take a uniform zero-detuning limit.

Let $dr/r=dy/\rho$, with $r/y\to1$ at the regular endpoint, and write
\[
 g=\Omega(s)^2\eta_{AB}dX^AdX^B,\quad
 \phi=\Phi(s),\quad s=R^2-T^2,\quad f=\Phi_s,
 \quad b=r_b,\quad\rho_b=b\Omega_b.                         \tag{4}
\]
The **actual Einstein–scalar background** is assumed analytically continued
on $-s_-<s\le b^2$, with $\Omega>0$ and $f\ne0$. The prior regular-cone
calculation explains this local continuation and its operator identities
[1]. An arbitrary extension of coefficient functions would not establish
the background equations.

The result is a sufficient family of coupled linearized solutions. For a
fixed finite set of $S^3$ harmonics, arbitrary smooth radial data supported
in a sufficiently near-shell annulus evolve on an explicit slab. Selected
pure displacement and pure velocity states make actual shell contact within
that slab. For $l\ge2$ the resulting nonzero states are not smooth pure
gauge. On a sufficiently thin fixed annulus the actual physical pairing is
weakly nondegenerate; selected displacement–velocity pairs both contact the
shell and retain a nonzero total pairing. Scalar uniqueness is proved for
this family. Completeness of gravitational perturbations and uniqueness in
a completed gauge quotient are separate questions.

# 2. Exact regular master and shell laws

Define $E=X^A\partial_A$ and $L=\Omega_s/\Omega+2f_s/f$. The regular master
and reconstruction from [1] are
\[
 \mathcal LF=\Box_\eta F-2LEF-[8L+(8/3)sf^2]F=0,
 \qquad P=r^4F,
\]
\[
 h_{AB}=2F(s\eta_{AB}-3X_AX_B),\qquad
 \chi=-\frac3{2\Omega^2f}(E+4)F.                           \tag{5}
\]
For arbitrary smooth $F$, the complete bulk residuals factor as
\[
 \delta\mathscr E_{AB}=\Omega^{-2}(3X_AX_B-s\eta_{AB})\mathcal LF,
 \qquad
 \delta\mathscr S=-\frac3{2\Omega^4f}(E+6)\mathcal LF.       \tag{6}
\]
These differential identities prove sufficiency without dividing by a
tangential harmonic eigenvalue and remain smooth at the cone and vertex.

With $A=(\Omega f^2)^{-1}$ and $\psi=\sqrt A F$, the canonical equation is
\[
 (\Box_\eta-V_c)\psi=0,\qquad
 V_c=3L-2sL_s+sL^2+(8/3)sf^2.                              \tag{7}
\]
The shell is $R_b(T)=\sqrt{b^2+T^2}$. Its exact master laws are
\[
 n_\eta\psi-c\Box_{\Sigma,\eta}\psi+d_e\psi=0,
 \qquad n_\eta\psi+d_o\psi=0,\qquad c=(\Omega_b B_b)^{-1}>0,
\]
\[
 d_j=\Omega_b\gamma_j+bL_b,\quad
 \gamma_e=4/\rho_b-4/(\rho_b^2B_b),\quad
 \gamma_o=4/\rho_b+p_b^2/(3a_b).                            \tag{8}
\]
The first law is even and the second odd under exchange of the interiors.
The shifted coefficients $d_j$ need not be positive. Positive $c$ is an
auxiliary scalar boundary inertia, not an independently added physical shell
field or a conclusion about the gravitational kinetic form.

# 3. Fixed hemisphere and controlled physical domain

Introduce $0\le\vartheta\le\pi/2$, $|\tau|<\pi/2$, and
\[
 Q=\frac b{\cos\tau+\cos\vartheta},\qquad
 T=Q\sin\tau,\quad R=Q\sin\vartheta,
\]
\[
 \eta=Q^2[-d\tau^2+d\vartheta^2+\sin^2\vartheta\,d\Omega_3^2],
 \quad
 s=b^2\frac{\cos\tau-\cos\vartheta}{\cos\tau+\cos\vartheta}. \tag{9}
\]
The shell is the equator of $M=S^4_+$, with $Q_b=b\sec\tau$ and boundary
$\Gamma=S^3$. Choose
\[
 0<\tau_*<\min\{\pi/4,\,2\arctan(\sqrt{s_-}/b)\},\qquad
 0<a_0<d_0<\tau_*/2.                                     \tag{10}
\]
The whole rectangle lies in the actual background because
$-b^2\tan^2(\tau_*/2)\le s\le b^2$. It contains every point of the full
ambient moving-ball slab
\[
 |T|\le T_0=b\tan(\tau_*/2),\qquad R\le R_b(T).            \tag{11}
\]
Indeed $\tau=\arctan((T+R)/b)+\arctan((T-R)/b)$ is maximal at $R=0$ for
fixed positive $T$. The shell itself reaches $T=b\tan\tau_*$. These are
different coverage bounds; a hemisphere rectangle is not an arbitrary
global ambient-time slab.

Set $\psi=Q^{-3/2}u$. Direct differentiation gives
\[
 u_{\tau\tau}-\Delta_Mu+\mathcal V u=0,
 \qquad\mathcal V=9/4+Q^2V_c(s).                         \tag{12}
\]
The curvature term is the five-dimensional conformal coupling $3/16$ times
the cylinder scalar curvature $12$. Equation (8) becomes
\[
 \partial_\nu u+\partial_\tau(\mu u_\tau)
       -\mu\Delta_\Gamma u+\beta u=0,
 \quad\mu=\frac{c\cos\tau}{b}>0,
\]
\[
 \beta=Q_b(d_e-3/(2b))
       -\mu[(3/2)\sec^2\tau+(3/4)\tan^2\tau],             \tag{13}
\]
or $\partial_\nu u+\alpha_o u=0$, with
$\alpha_o=Q_b(d_o-3/(2b))$. Here $\nu=\partial_\vartheta$.
All undifferentiated coefficients are time even; $\mu'$ is odd.

For reference, positive control of smooth even solutions follows from
\[
 \mathcal E=\tfrac12\int_M(u_\tau^2+|\nabla u|^2+m^2u^2)
  +\tfrac12\int_\Gamma[\mu(z_\tau^2+|\nabla z|^2)+kz^2],
 \qquad z=u|_\Gamma,
\]
\[
 \mathcal E'=\int_M(m^2-\mathcal V)uu_\tau
  +\tfrac12\int_\Gamma\mu'(|\nabla z|^2-z_\tau^2)
  +\int_\Gamma(k-\beta)zz_\tau.                           \tag{14}
\]
For $m,k>0$, bounded coefficients and $\mu_{\min}>0$ give
$|\mathcal E'|\le C\mathcal E$. In odd parity retain the actual boundary
potential $\frac12\int_\Gamma\alpha_o z^2$ and choose a sufficiently large
bulk mass. The epsilon trace inequality then gives coercivity and controls
$\alpha'_o z^2$. An independently shifted boundary potential cannot be
estimated by assuming an $L^2$ trace of the bulk velocity. Neither energy
is asserted to be the physical Einstein–scalar Hamiltonian.

# 4. Smooth construction and actual contact

Fix $-\Delta_{S^3}H_l=\lambda_l H_l$, $\lambda_l=l(l+2)$, and put
$x=\pi/2-\vartheta$. Take smooth initial profiles supported in
$[a_0,d_0]$. Such data vanish on an open shell neighborhood and therefore
satisfy every corner compatibility condition: recursively differentiate
(12) to obtain the bulk time jets and insert them into each time derivative
of (13) or the Robin law. All these jets vanish at the boundary. The full
hierarchy is recorded in the companion existence proof.

Writing $u=u_l(\tau,x)H_l$ and $v=\cos^{3/2}x\,u_l$ gives
\[
 v_{\tau\tau}-v_{xx}+W_l v=0,
 \qquad W_l=Q^2V_c+(\lambda_l+3/4)\sec^2x.                \tag{15}
\]
Let $p=v_\tau+v_x$, $m=v_\tau-v_x$, and $z=v(\tau,0)$.
For initial data $(v_0,v_1)$, set $p_0=v_1+\partial_xv_0$ and
$m_0=v_1-\partial_xv_0$. Collar support gives $z(0)=z'(0)=0$.
Then $p_\tau-p_x=-W_lv$ and $m_\tau+m_x=-W_lv$. The exact incoming
boundary prescription is
\[
 \mu z''+(\mu'+1)z'+(\mu\lambda_l+\beta)z=p(\tau,0),
 \qquad m(\tau,0)=2z'-p(\tau,0),                         \tag{16}
\]
or, in odd parity, $z'=p-\alpha_o z$ and $m=p-2\alpha_o z$.
Thus the dynamic boundary is a positive-mass ODE driven by the outgoing
characteristic, not two freely prescribed incoming data.

For $t=\tau\ge0$, the construction starts from
\[
 p(t,x)=p_0(x+t)-\int_0^tW_l(s,x+t-s)v(s,x+t-s)\,ds.       \tag{17}
\]
Integrate the other characteristic from the initial slice when $x\ge t$
and from (16) when $x<t$, then recover
$v(t,x)=v_0(x)+\frac12\int_0^t(p+m)(s,x)\,ds$.
Together with the integrated boundary ODE these are Volterra equations.
After substituting (17), their map on $(v,z,z')$ has Lipschitz constant
$C_l\Delta t$ in the uniform norm on a short interval; $C_l$ depends only
on bounded coefficients, $\mu_{\min}^{-1}$ and fixed $\lambda_l$.
The odd case uses $(v,z)$. Contraction gives uniqueness and existence.
Successive differentiation gives smoothness; the zero initial corner
neighborhood removes all matching jumps along its outgoing characteristic.
The characteristic constraints $v_x=(p-m)/2$ and $v(t,0)=z(t)$ propagate
from the initial data.

The interval length depends on coefficient bounds, not solution amplitude,
so iteration covers $|\tau|\le\tau_*$. The characteristic equations imply
\[
 \operatorname{supp}v(\tau,\cdot)
 \subset[\max(0,a_0-|\tau|),\,d_0+|\tau|].                \tag{18}
\]
Since $d_0+\tau_*<\pi/2$, no signal reaches the polar origin. Smoothly
extending the coefficient outside this causal collar is only a construction
device; the solution is zero before that extension can affect it. Extension
by zero near the origin yields a smooth full-hemisphere solution. Finite
linear combinations of harmonics follow by superposition. No uniform
$l\to\infty$ estimate or general smoothness theorem from Galerkin
compactness is presumed.

Finite propagation alone proves only that the shell trace is zero for
$|\tau|<a_0$. Actual activation needs another argument. For nonzero
displacement with support infimum exactly $a_0$ and zero velocity, time
reversal makes $v$ even. If its trace
were zero on $(0,\tau_1)$, $\tau_1>a_0$, it would vanish on
$(-\tau_1,\tau_1)$ and the boundary law would give $v_x=0$ there.
Regarding (15) as the ordinary one-dimensional wave equation
$v_{xx}-v_{\tau\tau}=W_lv$ with $x$ as time, zero Cauchy data at $x=0$
imply $v=0$ in $x+|\tau|<\tau_1$. This contradicts an initial profile
with support infimum $a_0$. For zero displacement and nonzero velocity,
the solution is time odd and the same argument forces
$v_\tau(0,x)=0$, giving the same contradiction. This sideways argument
is specifically one dimensional; it is not a general higher-dimensional
timelike Cauchy theorem.

Thus either choice has nonzero shell trace somewhere in every interval
immediately after $a_0$. No nonzero value or jet at exactly the smooth
support edge is claimed. In ambient time,
\[
 T_c=b\tan a_0=\frac{b^2-R_0^2}{2R_0}<T_0,
 \qquad R_0=b\tan(\pi/4-a_0/2).                          \tag{19}
\]
For prescribed master data, set $w=(\Omega f^2)^{1/2}$ and
$Q_0=b/(1+\cos\vartheta)$. On $T=\tau=0$ the exact conversion is
\[
 u_0=Q_0^{3/2}w^{-1}F_0,\qquad
 u_1=Q_0^{5/2}w^{-1}(F_T)_0.                             \tag{20}
\]
The weights are smooth and positive. Therefore both master states
$(q_0H_l,0)$ and $(0,q_0H_l)$ have the same contact infimum, where
$a_0$ is the infimum of the nonzero profile's support in $x$.

![Possible causal support for illustrative collar data, not an evolved waveform. The trace activates arbitrarily soon after the arrival infimum; smooth data remain flat at the exact front. The dashed curve bounds the fully covered ambient-time slab. Its scale is illustrative and does not certify a background continuation.](../figures/contact-support.pdf){width=80%}

Optionally imposing $d_0+2\tau_*<\pi/2$ keeps the whole causal support
in $\vartheta>|\tau|$, hence $s>0$. Contact can therefore occur entirely
within the original regular positive-$s$ collar. These data are a new
family, not a continuation of the earlier near-cone nonzero-flux pair.

# 5. Full junction sufficiency

Scalar existence becomes a coupled result only after checking the complete
shell system. In the exterior chart put $\Psi=P/\rho^2$; reconstruction
implies, for arbitrary smooth packets,
\[
 \Psi'+2a\Psi=-p\chi/3,\qquad
 p\chi'=p'\chi+3\rho^{-2}(\Box_\gamma+4)\Psi-2p^2\Psi.    \tag{21}
\]
For normal graphs $y_i=y_b+\zeta_i$ in a common shell chart, direct
linearized graph geometry gives
\[
 \delta h^{\rm ind}_{\mu\nu}=2\rho_b^2(\Psi+a\zeta)\gamma_{\mu\nu},
 \quad\delta\phi^{\rm ind}=\chi+p\zeta,
\]
\[
 \delta K^\mu{}_{\nu}=(\Psi'+2a\Psi+a'\zeta)\delta^\mu{}_{\nu}
       -\rho^{-2}D^\mu D_\nu\zeta,
 \quad\delta p_n=\chi'+2p\Psi+p'\zeta.                    \tag{22}
\]
The graph evaluation and lapse terms are essential. Write $P_\pm=P_e\pm P_o$
and choose
\[
 \zeta_e=0,\qquad\zeta_o=-P_{ob}/(a_b\rho_b^2),\qquad
 \zeta_\pm=\pm\zeta_o.                                  \tag{23}
\]
Induced metric continuity is immediate; scalar continuity is precisely
$P_o'=-p_b^2P_o/(3a_b)$, the odd law. The full trace-reversed Israel
residual and summed scalar residual reduce to
\[
 \delta\mathcal R^\mu{}_{\nu}
       =-2\rho_b^{-2}(D^\mu D_\nu+\delta^\mu{}_{\nu})\zeta_e,
\]
\[
 \delta\mathcal S_\Sigma
 =-\frac6{\rho_b^2p_b}
       [B_bP_e'-\rho_b^{-2}(\Box_\gamma+4)P_e]
       +2p_bB_b\zeta_e.                                  \tag{24}
\]
Both vanish by (23) and the even law. This checks all tensor components,
not merely the trace-free part. The difference of the bulk Hamiltonian
constraints supplies normal-force consistency:
$3a_b\delta K_o=p_b\delta p_{n,o}$, equivalently
$\sigma_b\delta K_o+\sigma'_b\delta p_{n,o}=0$.
Codazzi gives tangential conservation. There is no extra odd wave equation
to impose.

These are sufficiency statements without harmonic division. Necessity for
all embeddings is different: Israel permits
$(D_\mu D_\nu+\gamma_{\mu\nu})\zeta_e=0$, with the five de Sitter
ambient-coordinate solutions. The scalar condition then still contains
$\zeta_e$. Setting it to zero selects a sufficient family and does not
classify these exceptional representatives.

# 6. Physical current and nonzero conserved pairing

Fix the action-derived potential convention
\[
 \theta^A=\sqrt{-g}\{\tfrac12(\nabla_Bh^{AB}-\nabla^Ah)
                        -\nabla^A\phi\,\chi\},\qquad
 \omega^A=\delta_1\theta^A(2)-\delta_2\theta^A(1).
\]
For one interior the exact reconstructed currents are
\[
 \omega^y=9\sqrt{-\gamma}\,\mathfrak a
 [((\Box_\gamma+4)P_1)P'_2-((\Box_\gamma+4)P_2)P'_1],
 \qquad\mathfrak a=(\rho^2p^2)^{-1},
\]
\[
 \omega^\mu=\sqrt{-\gamma}\rho^2
 [\chi_1D^\mu\chi_2-\chi_2D^\mu\chi_1
       +6(\Psi_1D^\mu\Psi_2-\Psi_2D^\mu\Psi_1)].          \tag{25}
\]
For two even states the boundary law makes $\omega^y_b=0$ pointwise;
there is no additional straight-even-shell current in this slice. For two
odd states the completed graph and shell derivative action is
\[
 S^{(2)}_{\Sigma,\rm der}
   =-3a_b\rho_b^2\int\sqrt{-\gamma}(D\zeta_o)^2
   =-9c_o\int\sqrt{-\gamma}(DP_{ob})^2,
 \qquad c_o=(3a_b\rho_b^2)^{-1}.                          \tag{26}
\]
It includes both GHY contributions, the Einstein tangential-divergence
corner, the shifted bulk integration domains, and one common shell tension.
Omitting the moving-domain term reverses the resulting coefficient. The
companion junction calculation derives these terms independently.

Polarization gives the physical shell current and balance
\[
 \omega_\Sigma^\mu=18c_o\sqrt{-\gamma}
 (P_{1b}D^\mu P_{2b}-P_{2b}D^\mu P_{1b}),\qquad
 \partial_\mu\omega_\Sigma^\mu=\omega_+^y(y_b)+\omega_-^y(y_b). \tag{27}
\]
Mixed reflection parities cancel. In ambient polar coordinates the normal
moving-boundary flux is
$\omega^R-R_b'\omega^T=\omega^y/R_b$. Thus, for a full interior with
no additional boundary, the future-oriented total pairing
\[
 \sigma_T=-\left[\sum_i\int_{R\le R_b(T)}\omega_i^T\,dR\,d\Omega
                     +\int_{S^3}\omega_\Sigma^T\,d\Omega\right]
 \quad\hbox{satisfies}\quad\frac{d\sigma_T}{dT}=0.          \tag{28}
\]
The overall minus sign specifies the future-Cauchy convention for the
displayed vector densities; reversing the common orientation changes no
conservation statement. The cone is internal to the smooth domain. A region
artificially cut there must retain its cone flux separately.

At $T=0$, let $(q,v)=(F,F_T)$ have compact support in an annulus inside the
shell. The shell cut term is initially zero. Per interior, direct pullback
of (25) yields
\[
 \sigma((q_1,v_1),(q_2,v_2))=K(q_1,v_2)-K(q_2,v_1),
\]
\[
 K(q,v)=\int dr\,d\Omega_3
 \left[m_r(D_r+4)q(D_r+5)v+\frac{6r^7}{\Omega}qv\right],
 \quad D_r=r\partial_r,\quad m_r=\frac{9r^3}{4\Omega f^2}>0. \tag{29}
\]
The unit difference between the two first-order factors follows from
$[\partial_T,E]=\partial_T$. Compact radial integration by parts gives
\[
 K(q,q)=\int dr\,d\Omega_3
 \left[m_r\{(D_rq)^2+(2+9sL)q^2\}+\frac{6r^7}{\Omega}q^2\right]. \tag{30}
\]
For $I=(r_1,r_2)$, write $w_I=r_2-r_1$,
$m_* =\min_{\bar I}(m_rr^2)>0$ and
$M_* =\max_{\bar I}m_r(-2-9sL)_+$. Dirichlet Poincaré gives
\[
 K(q,q)\ge\left(\frac{\pi^2m_*}{w_I^2}-M_*\right)
                  \int_I dr\,d\Omega_3\,q^2.             \tag{31}
\]
Every regular near-shell collar contains an annulus thin enough for the
bracket to be positive while satisfying (10). Hence two same-reflection-
parity states $(q,0)$ and $(0,q)$ have doubled physical pairing
$2K(q,q)>0$. Section 4 proves actual contact for both; (28) conserves their
nonzero **total** pairing afterward, including the odd cut term. Nonzero
boundary traces do not establish a nonzero local shell-current Wronskian,
and no such claim is made.

On all smooth data supported in this fixed thin annulus with a fixed finite
harmonic set, the pairing is weakly nondegenerate. A kernel element $(q,v)$
tested against $(0,q)$ has $q=0$ by (31); testing next against $(v,0)$ gives
$v=0$. Smooth uniqueness and (28) transport this statement to the solution
image. $K$ need not be symmetric. Positivity of its diagonal is not a
physical energy theorem or a completed-space identification.

# 7. Gauge obstruction and scope

For $l\ge2$, non-gauge character also follows directly on a polar annulus.
Write $g=g_{ab}(z)dz^adz^b+w(z)^2q_{ij}d\omega^id\omega^j$ and the scalar
part of a gauge covector as $\xi_a=f_aH_l$, $\xi_i=JD_iH_l$. Then
\[
 (\mathcal L_\xi g)_{ij}^{\rm TF}=2J(D_iD_jH_l)^{\rm TF},\quad
 (\mathcal L_\xi g)_{ai}=[\partial_aJ-2(w_a/w)J+f_a]D_iH_l,
\]
\[
 \int|(\operatorname{Hess}H_l)^{\rm TF}|^2
    =\tfrac23\lambda_l(\lambda_l-3)\int|H_l|^2>0.          \tag{32}
\]
The reconstruction has angular trace only and zero mixed component, so
$J=0$ and then $f_a=0$. Orthogonal angular sectors cannot cancel this
scalar projection. A nonzero reconstructed $l\ge2$ state is therefore not
smooth pure gauge. This does not classify $l=0,1$, finite-cut charges or
weaker gauge groups.

The result combines four separate ingredients: actual-background bulk
reconstruction, smooth finite-harmonic boundary evolution, complete junction
sufficiency, and the physical current including shell motion. A scalar
boundary theorem alone would not imply the other three. General angular
Sobolev completion, unrestricted initial-boundary wellposedness, all-sector
completeness, nonlinear integrability and stability remain open here. The
construction gives no reheating mechanism, observational discrimination or
evidence for a higher-dimensional origin of the Big Bang.

# 8. Prior art, verification and declarations

Positive dynamic-boundary inertia, scalar compatibility and energy methods
are established subjects. Zahn [3] treats static flat generalized Wentzell
boundaries; Vitillaro [4] proves nonlinear energy and stronger-regularity
results on fixed bounded domains. Marta [5] studies time-dependent compact
Riemannian metrics. These are precedents, not already verified theorems for
the present moving shell and coupled reconstruction. We use our own finite-
coefficient Gronwall bound and explicit characteristic proof.

The recent Robin treatment of Costeri, Dappiaggi and Juárez-Aubry [6]
provides smooth interior-compact Cauchy existence and causal support on
spacetimes with timelike boundary. Its inward-normal convention maps to
$\kappa=-d_o$. Extending its stated curvature-mass family to $V_c$ requires
an explicit smooth lower-order adaptation. Its later microlocal reflection
hypotheses are not silently added to this finite-slab smooth construction.
The accompanying literature register records primary-source reading and
failed searches. It supports a bounded comparison, not external or worldwide
originality. No new general scalar boundary theory is claimed.

The package contains detailed existence, junction/current and independent
audit notes, plus original symbolic programs and strict wrong-formula
controls. They check coordinate and conformal transformations, full graph
geometry, arbitrary-packet residuals, moving-shell coefficients, and the
annular Cauchy identity. The explicit contracts contain **84 positive
algebraic identities and 14 wrong-formula controls per execution mode**:
21/4 for the coordinate checker and 63/10 for junctions and currents.
Normal and optimized individual outputs agree byte for byte. The portable
replay additionally checks relocation, input integrity, receipt equality and
injected failures; its receipts record the actual execution status.
Algebraic check counts are not a count of independent physical results or
a machine-certified PDE theorem.

**Author and support.** Ricardo Maldonado, Independent Researcher,
ORCID 0009-0009-3937-6527; correspondence:
ricardomaldonadoandrade@gmail.com. No external funding is reported. The
author reports no known conflicts of interest.

**AI assistance and review.** Substantial OpenAI Codex assistance was used
for derivation, symbolic programming, literature organization, internal
checking and drafting. Separate agent processes supplied internal reviews;
these are not external human peer review. The author retains responsibility
for verification, interpretation and any submission. APS submission and
external novelty assessment remain open; neither acceptance nor readiness
for acceptance is asserted.

# References

1. R. Maldonado, *Regular cone evolution and local symplectic flux in an
   Einstein–scalar shell model*, technical checkpoint, 13 October 2026,
   frozen `HDBLAST_REGULAR_CONE_EVOLUTION_20261013` capsule, especially its
   reconstruction and boundary notes.
2. R. Maldonado, *From scalar-boundary Jacobi response to spectral positivity
   of regular de Sitter shells*, registered spectral checkpoint (2026),
   DOI [10.5281/zenodo.23335485](https://doi.org/10.5281/zenodo.23335485).
3. J. Zahn, *Generalized Wentzell boundary conditions and quantum field
   theory*, Ann. Henri Poincaré **19**, 163–187 (2018),
   [arXiv:1512.05512v3](https://arxiv.org/abs/1512.05512v3).
4. E. Vitillaro, *On the wave equation with hyperbolic dynamical boundary
   conditions, interior and boundary damping and source*, Arch. Ration.
   Mech. Anal. **223**, 1183–1247 (2017),
   [arXiv:1506.00910v4](https://arxiv.org/abs/1506.00910v4).
5. A. Marta, *The nonlinear wave equation with nonlinear Wentzell boundary
   conditions on time-dependent compact Riemannian manifolds* (2024),
   [arXiv:2404.05855v2](https://arxiv.org/abs/2404.05855v2).
6. B. Costeri, C. Dappiaggi and B. A. Juárez-Aubry, *The Hadamard parametrix
   on globally hyperbolic spacetimes with Robin boundary conditions:
   Fundamental solutions and Hadamard states* (2026),
   [arXiv:2609.24465v1](https://arxiv.org/abs/2609.24465v1).

**Companion proof files.**
[Constructive existence](../existence/FINITE_SLAB_SHELL_CONTACT_EXISTENCE.md);
[junctions and physical current](../junctions/FULL_PACKET_JUNCTION_SUFFICIENCY_AND_CURRENT.md);
[independent energy audit](../existence/independent/ENERGY_EXISTENCE_AUDIT.md);
[bounded primary-source comparison](../literature/SHELL_CONTACT_PRIOR_ART_AND_ASSUMPTIONS.md).
