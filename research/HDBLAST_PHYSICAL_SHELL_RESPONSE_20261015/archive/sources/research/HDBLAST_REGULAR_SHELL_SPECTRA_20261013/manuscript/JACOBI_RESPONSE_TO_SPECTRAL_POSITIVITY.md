---
title: "From scalar-boundary Jacobi response to spectral positivity of regular de Sitter shells"
author: "Ricardo Maldonado"
date: "13 October 2026 — technical continuation draft"
fontsize: 11pt
geometry: margin=1in
colorlinks: true
header-includes:
  - \usepackage{mathrsfs}
---

**Ricardo Maldonado**  
Independent Researcher  
ORCID: [0009-0009-3937-6527](https://orcid.org/0009-0009-3937-6527)  
Email: ricardomaldonadoandrade@gmail.com

**Draft status.** This is a separate technical continuation, not a submitted
manuscript. The local even and completed odd kinetic calculations have
received independent internal review. Global gauge, cone-flux and physical
state-space assumptions remain explicit.

## Abstract

We connect the scalar-boundary Jacobi response of a certified family of
regular de Sitter shells to positivity of its radial scalar operators. A
specific combination of cone Jacobi fields gives a nodeless zero-mass
reference whose shell residual is strictly positive throughout the
nonconstant branch \(0<\delta\le 0.002\). For the reflection-even scalar
problem, we construct the self-adjoint realization with a positive endpoint
mass, prove the ground-state identity on its full closed form domain, and
derive an explicit positive gap as a function of each fixed background.
The reflection-odd matching instead gives an ordinary Robin realization
with an analogous bound. The scalar essential threshold is
\(\mu^2=9/4\). A tensor partner comparison excludes every nonzero tensor
spectral value below \(9/4\) on the same backgrounds. These are applications
of established spectral methods. Physical interpretation also requires
the separately stated gauge, cone-domain, exceptional-harmonic and
kinetic-action assumptions; no complete nonlinear stability or uniform
positive physical mass gap is asserted.

## 1. Problem and relation to established methods

Stabilizing a static shell does not by itself exclude a scalar instability.
The coupled Einstein–scalar equations on de Sitter slices contain a
curvature shift, so positivity of the shifted radial eigenvalue is weaker
than positivity of the four-dimensional mass [1,2]. Soft scalar junctions
already lead to eigenparameter-dependent boundary conditions and boundary
terms in the Rayleigh denominator [2–4]. Enlarged Hilbert spaces for such
conditions [5], ground-state transformations [6], and factorized tensor
operators [7,8] are established tools.

The contribution considered here is a particular exact application. The
certified static branch supplies pointwise Jacobi signs. An explicit
combination of those Jacobi fields solves the generic scalar perturbation
equations at zero mass and has a strict shell residual. This additional
structure, together with a closed-domain argument, upgrades a static
response calculation to a statement about the entire spectrum of specified
radial operators. The decisive estimate does not divide an approximation
error by the small central scalar slope and does not require an eigenvalue
scan.

The accompanying bounded primary-source review compares these identities
with prior scalar junction, radion, and tensor analyses. It does not
establish exhaustive priority or global novelty of the mathematics.
In particular, a positive radial Hilbert norm must be distinguished from a
physical kinetic coefficient [9]. A doubled background also allows odd
perturbations unless a reflection quotient is imposed; a calculation of
the even sector alone does not settle that larger system [10,11].

The single-brane instability result of Lesgourgues and Sorbo [3] warrants
a specific assumption comparison. Their setup has flat Minkowski slices,
a compact orbifold and a smooth second fixed point with a nonzero warp
factor and vanishing first scalar and warp derivatives. Their Eqs. (16)–(17)
already contain the spectral parameter in the scalar junction. Here the
other endpoint is a regular de Sitter cone with \(\rho(0)=0\),
\(a\sim1/y\), and an infinite conformal throat. Their theorem therefore
does not directly decide this boundary-value problem. The relevant
distinction is the slicing and endpoint hypotheses, not the label
“one brane”; the present application is not a counterexample to their
theorem.

## 2. Exact background input and conventions

Use unit-curvature \(dS_4\), with signature \((-+++)\), and two identical
interiors joined at \(y=T\):
\[
 ds^2=dy^2+\rho(y)^2\gamma_{\mu\nu}dx^\mu dx^\nu,\qquad
 p=\eta',\quad a=\rho'/\rho,\quad {\cal H}=\rho^{-2}.
 \tag{1}
\]
The action includes the boundary completion,
\[
 S=\sum_{i=\pm}\left[
 \int_{M_i}\sqrt{-g}\left(\frac R2-\frac12(\nabla\eta)^2-U\right)
 +\int_\Sigma\sqrt{-h}\,K_i\right]
 -\int_\Sigma\sqrt{-h}\,\sigma(\eta).
 \tag{2}
\]
Both shell normals point toward increasing radius in their own interiors.
The background equations and identical-copy junctions are
\[
 \begin{aligned}
 p'&=U'-4ap,&
 a'&=-{\cal H}-p^2/3,&
 a^2&={\cal H}+p^2/12-U/6,\\
 a(T)&=\sigma/6,&p(T)&=-\sigma'/2.
 \end{aligned}                                                   \tag{3}
\]
Primes on potentials denote field derivatives.

The registered model is
\[
 \begin{aligned}
 W(\eta)&=\frac13+\eta^2+\frac{\eta^3}{3},&
 U&=\frac12(W')^2-\frac23W^2,\\
 \sigma&=2W+\delta(1+c+c\eta),&
 c&=\frac{2}{1.0357712571566784}-\frac43
   =\frac{96714321425692}{161839258930731}.
 \end{aligned}                                                   \tag{4}
\]
The decimal in (4) is interpreted exactly as a rational number.
The October 11 static theorem supplies, for every
\(0<\delta\le1/500\), one regular solution in its specified small Dirichlet
tube, with \(T\ge9\), \(|\eta|\le1/100\), \(a>0\), and \(p<0\) on
\((0,T]\). Uniqueness is restricted to that tube. The theorem is a
mathematical branch statement; it does not certify historical floating-point
trajectories individually.

The field shift \(\eta=\phi-1\) rewrites the earlier
\(W_{\rm old}(\phi)=1-\phi+\phi^3/3\) exactly as (4).
The previously investigated unstable hilltop configuration and the
regular branch treated here are distinct backgrounds of the same action:
the former has a cone near \(\phi=-1\), whereas this certified tube is
near \(\phi=+1\). Tube-local uniqueness does not exclude the former
background. No dynamical evolution connecting the branches is established.
By contrast, replacing the affine detuning by a scalar-independent
detuning changes the boundary theory; that earlier curvature comparator
is not another branch of the same shell action.

Write \(H_b^2=\rho(T)^{-2}\). Our harmonic convention is
\[
 \Box_\gamma Y=\mu^2Y,\qquad m_4^2=H_b^2\mu^2,\qquad
 \Lambda=\mu^2+4.                                      \tag{5}
\]
Thus exclusion of nonpositive scalar masses requires exclusion of
\(\Lambda\le4\), not merely \(\Lambda\le0\). The constant-scalar
\(\delta=0\) limit is outside the divided scalar formulation below.

## 3. The Jacobi reference and its strict shell residual

Let \(u=\eta(0)\) label the regular cone family at fixed radial coordinate,
and define
\[
 h=\eta_u,\quad q=p_u=h',\quad
 v=(\log\rho)_u,\quad b=a_u=v',\quad
 j=b+ph/3,\quad {\cal D}=hp'-pq.
 \tag{6}
\]
The static proof establishes \(h\ge1\), \(q>0\), \(pj>0\), and
\(p{\cal D}>0\) on a nonconstant interior. It also gives
\[
 J=\frac{q(T)}{h(T)}+\frac{\sigma''}{2}\ge1,\qquad
 B_\partial=\frac{p'(T)}{p(T)}+\frac{\sigma''}{2}>J,
 \quad B_\partial-J=\frac{{\cal D}(T)}{p(T)h(T)}.
 \tag{7}
\]
The relevant constrained Jacobi identities are
\[
 j'+4aj=-2{\cal H}v,\qquad
 {\cal D}'+4a{\cal D}=4p({\cal H}h+pj),\qquad
 pq-p'h=12aj+12{\cal H}v.
 \tag{8}
\]

For a generic scalar harmonic, the longitudinal gauge-invariant fields
\((\Psi,X)\) obey
\[
 \Psi'=-2a\Psi-pX/3,\qquad
 pX'=p'X+(3\Lambda{\cal H}-2p^2)\Psi.             \tag{9}
\]
Set
\[
 K=\rho^2j,\qquad
 X_0=h+pK,\qquad
 \Psi_0=v+aK=-\rho^2{\cal D}/12.
 \tag{10}
\]
Equation (8) gives \(K'=-2\Psi_0\); substitution proves both equations
(9) at \(\Lambda=4\). Thus (10) is an actual zero-mass bulk solution,
not the uncorrected static variation \((h,v)\). Moreover
\(X_0>h\ge1\) and \(p\Psi_0<0\). On the registered branch \(\Psi_0>0\),
so \(P_0=\rho^2\Psi_0\) is positive throughout the interior.

The full even scalar shell mismatch is
\[
 M=X'+\frac{\sigma''}{2}X+2p\Psi.
 \tag{11}
\]
Retaining its gravitational term, (10) gives
\[
 M_0=hJ+B_\partial pK,\qquad
 \frac{M_0}{X_0}
 =\frac{J+B_\partial\alpha_K}{1+\alpha_K}>J\ge1,
 \qquad \alpha_K=\frac{pK}{h}>0.                       \tag{12}
\]
This is a normalization-independent strict margin. It requires no lower
bound on the small central slope. The cone expansions are
\[
 p=\frac{U'(u)}5y+O(y^3),\quad
 X_0=1+O(y^2),\quad
 \Psi_0=-\frac{U'(u)}{60}y^2+O(y^4),\quad
 P_0=-\frac{U'(u)}{60}y^4+O(y^6).
 \tag{13}
\]
The static junction determinant is also positive:
\(\det D_{(u,T)}(a-\sigma/6,p+\sigma'/2)={\cal H}_bM_0>0\).
Nodelessness follows from the pointwise signs in (10), separately from
that determinant.

## 4. A closed scalar realization and a positive gap

### 4.1 Endpoint mass and form domain

For \(P=\rho^2\Psi\), define
\[
 \alpha=\frac1{\rho^2p^2},\qquad
 q_{SL}=\frac2{3\rho^2},\qquad
 w_{SL}=\frac1{\rho^4p^2},\qquad
 \beta_{SL}=\rho(T)^4p(T)^2B_\partial>0.
 \tag{14}
\]
Here \(\alpha\) is the Sturm–Liouville leading coefficient, not the
background momentum. Equations (9) and (11) give
\[
 -(\alpha P')'+q_{SL}P=\Lambda w_{SL}P,\qquad
 \alpha(T)P'(T)=\frac{\Lambda}{\beta_{SL}}P(T).
 \tag{15}
\]
In particular \(P'=-\rho^2pX/3\). The reference residual is
\[
 r_0=\frac{\alpha(T)P_0'(T)}{P_0(T)}-\frac4{\beta_{SL}}
 =-\frac{M_0}{3\rho(T)^2p(T)B_\partial\Psi_0(T)}>0.
 \tag{16}
\]

The even Hilbert space is
\[
 {\mathscr H}_e=L^2(w_{SL}dy)\oplus\mathbb C,\qquad
 \|(P,z)\|_e^2=\int_0^T w_{SL}|P|^2dy+\frac{|z|^2}{\beta_{SL}}.
 \tag{17}
\]
Its components are independent. Trace identification \(z=P(T)\) belongs
to the form domain
\[
 V_e=\left\{(P,P(T)):
 P\in AC_{\rm loc}((0,T]),\
 \int_0^T(w_{SL}|P|^2+\alpha|P'|^2)dy<\infty\right\}.
 \tag{18}
\]
The energy is
\({\mathfrak e}[P]=\int(\alpha|P'|^2+q_{SL}|P|^2)\).
Since \(q_{SL}/w_{SL}\) is bounded, it introduces no additional
integrability condition. Local \(H^1\) convergence proves closedness
and identifies the shell trace. Narrow shell bumps establish density in
the enlarged Hilbert space.

Smooth functions vanishing near the cone, with unrestricted shell trace,
form a core. For a cutoff over \([\epsilon,2\epsilon]\), its extra
derivative energy is bounded by a constant times the weighted \(L^2\)
tail, because \(\alpha/w_{SL}=\rho^2=O(y^2)\).
The representation theorem therefore gives the self-adjoint Friedrichs
operator
\[
 A_e(P,P(T))=
 \left(w_{SL}^{-1}[-(\alpha P')'+q_{SL}P],\
       \beta_{SL}\alpha(T)P'(T)\right).                 \tag{19}
\]
Its domain is (18) with \(\alpha P'\) locally absolutely continuous and
the first component in \(L^2(w_{SL}dy)\). Equation (15) is precisely its
eigenvalue equation.

At \(\Lambda=4\), reduction of order gives a second solution
\(P_0\int(\alpha P_0^2)^{-1}dy\sim y\); \(P_0\sim y^4\).
Only the latter is in \(L^2(w_{SL}dy)\), since \(w_{SL}\sim
A_c^{-2}y^{-6}\), \(A_c=U'(u)/5\). The cone is limit-point.
No independent cone boundary parameter is introduced.

### 4.2 Ground-state identity and quantitative coercivity

Let
\[
 a_g=\alpha P_0^2,\quad m_g=w_{SL}P_0^2,\quad
 b_g=P_0(T)^2/\beta_{SL},\quad k_g=r_0P_0(T)^2>0.
 \tag{20}
\]
For \(P=P_0f\), the closed ground-state identity is
\[
 {\mathfrak e}[P]-4\|(P,P(T))\|_e^2
 =\int_0^T a_g|f'|^2dy+k_g|f(T)|^2=:D[f].
 \tag{21}
\]
It holds first on the cone core by integration by parts. On that core,
the original form norm equals \(D+5\|P\|_e^2\); the same cutoff argument,
now using \(a_g/m_g=\rho^2\), proves the identity on the full form domain.
This avoids assigning Frobenius expansions to arbitrary form vectors.

Strict positivity of (21) for each nonzero vector would not alone exclude
zero from continuous spectrum. A weighted Hardy estimate supplies the
stronger conclusion. Define
\[
 \alpha_*=\inf_{0<y\le T}\frac{a_g(y)}{y^4}>0,\quad
 m_*=\sup_{0<y\le T}\frac{m_g(y)}{y^2}<\infty,\quad
 M_g=\int_0^T m_gdy.
 \tag{22}
\]
The limits \(a_g/y^4,m_g/y^2\to1/144\), together with positivity on
compact subintervals, establish these bounds at each fixed background.
For \(g=f-f(T)\), integration by parts on \([\epsilon,T]\) gives
\[
 3\int_\epsilon^T y^2|g|^2
 =-\epsilon^3|g(\epsilon)|^2
 -2\Re\int_\epsilon^T y^3\bar g g'
 \le2\left(\int_\epsilon^T y^2|g|^2
             \int_\epsilon^T y^4|g'|^2\right)^{1/2}.
 \tag{23}
\]
Hence \(\int m_g|f-f(T)|^2\le
4m_*/(9\alpha_*)\int a_g|f'|^2\). The Hilbert-space triangle and
Cauchy–Schwarz inequalities then yield
\[
 \boxed{\quad A_e-4I\ge\varepsilon_e I,\qquad
 \varepsilon_e=
 \left[\frac{4m_*}{9\alpha_*}+\frac{M_g+b_g}{k_g}\right]^{-1}>0.\quad}
 \tag{24}
\]
Thus the entire shifted spectrum, including its continuous part, lies
above a positive number. Self-adjointness also excludes complex spectral
values. Formula (24) is an explicit function of a background; no numerical
value or lower bound uniform in \(\delta\) has been evaluated here.
Its value is unchanged by rescaling \(P_0\).

### 4.3 Essential spectrum and odd Robin realization

Assume the regular-cone expansions through two logarithmic derivatives.
The Liouville map
\[
 x=\int_y^T\frac{dt}{\rho(t)},\qquad
 F=\frac{P}{\rho^{3/2}|p|},\qquad
 s_L=\rho\left(\frac32a+\frac{p'}p\right)
 \tag{25}
\]
puts the cone at \(x=\infty\) and the shell at \(x=0\). The transformed
potential is
\[
 V_L=\frac23\rho^2p^2+(s_L)_x+s_L^2
       =\frac{25}{4}+O(e^{-2x}),
 \quad
 -F_x(0)+s_L(0)F(0)=
       \frac{\Lambda}{\rho(T)B_\partial}F(0).
 \tag{26}
\]
The endpoint mass remains positive. A decaying bounded potential
perturbation is relatively compact on the half-line. Comparing the
endpoint realization with a Dirichlet operator having the same actual
potential changes its resolvent by finite rank. Consequently
\[
 \sigma_{\rm ess}(A_e-4I)=[9/4,\infty).
 \tag{27}
\]
Any spectrum below that threshold is discrete. The best actual gap
cannot exceed \(9/4\).

If both interiors fluctuate independently, induced scalar and metric
continuity give \(X_o+p\xi_o=0\) and \(\Psi_o+a\xi_o=0\).
The generic summed trace-free junction fixes \(\xi_e=0\).
Thus the even mismatch is (11), whereas odd continuity implies
\[
 \alpha(T)P_o'(T)+c_oP_o(T)=0,\qquad
 c_o=\frac1{3a(T)\rho(T)^2}>0.                         \tag{28}
\]
The odd Hilbert norm has no endpoint mass:
\(\|P\|_o^2=\int w_{SL}|P|^2\).
Its form includes \(+c_o|P(T)|^2\).
With \(k_o=[\alpha P_0'/P_0+c_o]_T P_0(T)^2>0\), the same proof gives
\[
 A_o-4I\ge\varepsilon_o I,\qquad
 \varepsilon_o=
 \left[\frac{4m_*}{9\alpha_*}+\frac{M_g}{k_o}\right]^{-1}>0,
 \quad \sigma_{\rm ess}(A_o-4I)=[9/4,\infty).
 \tag{29}
\]
This operator statement assumes the displayed two-copy matching and its
finite-energy realization. It does not identify the odd radial norm with a
physical kinetic norm by definition.

## 5. Exceptional harmonic subspaces and physical kinetic information

### 5.1 Covariant exceptional classification

The generic reduction uses both \(D_\mu Y\) and the trace-free Hessian
\(T_{\mu\nu}[Y]=D_\mu D_\nu Y-\gamma_{\mu\nu}\Box_\gamma Y/4\).
Their degeneracies require separate analysis [12]. On unit \(dS_4\),
\[
 D^\mu T_{\mu\nu}[Y]=\frac34(\mu^2+4)D_\nu Y.
 \tag{30}
\]
Thus \(DY=0\) identifies the constant subspace at \(\mu^2=0\).
If \(T[Y]=0\) but \(DY\ne0\), then \(\mu^2=-4\) and
\(D_\mu D_\nu Y=-\gamma_{\mu\nu}Y\).
Other harmonics at these same mass values remain generic. For example,
the homogeneous de Sitter exponents \(r=0,1\) are exceptional, while
\(r=-3,-4\) have the same respective masses and are generic.

For the constant sector, assume the regular static Jacobi family spans
the admitted perturbations modulo radial gauge. The even junction
determinant is \({\cal H}_bM_0>0\). The odd induced-field matrix has
determinant
\[
 \det\begin{pmatrix}h&p\\v&a\end{pmatrix}
 =ah-pv=aX_0-p\Psi_0>0.                               \tag{31}
\]
There is therefore no constant modulus within that complete regular family.
For the Hessian-degenerate sector, the local bulk spanning calculation is
pure gauge. If the required transformations are globally admissible,
gauge the bulk fields to zero and retain the two embedding perturbations
\(\Delta_\pm\). Scalar continuity and the summed scalar junction give
\[
 p(\Delta_+-\Delta_-)=0,\qquad
 pB_\partial(\Delta_++\Delta_-)=0.
 \tag{32}
\]
Their determinant has magnitude \(2p^2B_\partial>0\), excluding the
relative embeddings. These conclusions depend on the stated completeness
and global gauge assumptions; they are not obtained by substituting
exceptional mass values into a divided generic equation.

### 5.2 Local even and completed odd kinetic coefficients

The quadratic action, including Einstein–Hilbert, GHY and shell terms,
must decide the kinetic sign. In a zero-shift radial ADM decomposition,
the tangential two-derivative density before the generic constraint is
\(\rho^2[3DN\!\cdot\!D\psi+3(D\psi)^2-(DX)^2/2]\).
On the longitudinal constraint surface \(N=-2\Psi\), the doubled even
kinetic coefficient is
\[
 K_e=2\int_0^T\rho^2(|X|^2+6|\Psi|^2)dy
     =18{\mathfrak e}[P]
     =18\Lambda\|(P,P(T))\|_e^2.
 \tag{33}
\]
The last equality holds at spectral roots. The first equality is a local
action calculation with a straight shell; the middle equality follows
from \(P'=-\rho^2pX/3\). The factor \(\Lambda\) is essential:
the radial norm alone is not the physical kinetic coefficient.
This gives the healthy kinetic sign for finite-energy nonzero generic
even profiles with the admitted longitudinal representative.

An odd profile has a moving shell in longitudinal gauge. The bulk
calculation alone gives
\[
 K_{o,\mathrm{bulk}}=
 18\Lambda\int w_{SL}|P_o|^2dy-18c_o|P_o(T)|^2.
 \tag{34}
\]
The completed graph-boundary calculation includes the terms missing from
that bulk expression. Write the shell as \(y=T+\epsilon\zeta Q(x)\), with
odd continuity \(\psi_b=-a_b\zeta\), \(X_b=-p_b\zeta\).
The tangential Einstein divergence on a moving radial domain supplies a
shell flux. Together with graph GHY minus coordinate-slice GHY, its
one-copy quadratic action coefficient is
\(+3a_b\rho_b^2\zeta^2(DQ)^2\).
The shifted upper limit of the radial ADM intrinsic-curvature term
supplies \(+3\rho_b^2\zeta\psi_b(DQ)^2\) per copy.
The one common shell volume contributes
\(-\sigma_b\rho_b^2\zeta^2(DQ)^2/2\).

Adding both gravitational copies and one shell, then using
\(\sigma_b=6a_b\) and continuity, gives the boundary action coefficient
\(-3a_b\rho_b^2\zeta^2(DQ)^2\). Thus
\[
 K_{o,\partial}=6a_b\rho_b^2\zeta^2
               =18c_o|P_o(T)|^2,\qquad
 \boxed{K_o=18[{\mathfrak e}[P_o]+c_o|P_o(T)|^2]
             =18\Lambda\int_0^T w_{SL}|P_o|^2dy.}
 \tag{35}
\]
The last equality again holds at spectral roots. The two gravitational
edge terms cancel after odd continuity is imposed; this explains why
the final contribution equals the shell-volume value. Assuming that
equality before calculating the moving Einstein and GHY terms would
miss essential sign checks.

Compactly supported tangential test variations justify the integrations
by parts. They vanish near auxiliary initial and final time surfaces,
where the graph therefore has no perturbed corner angle. This determines
the local kinetic bilinear without assuming a moving graph is orthogonal
to arbitrary time boundaries. An independent covariant graph derivation
agrees with the direct normal calculation.

Equations (33) and (35) supply healthy local kinetic signs for generic
admitted finite-energy profiles. They do not prove that every physical
state has a regular longitudinal representative or that the global
covariant symplectic flux has the prescribed cone limit. Those state-space
and gauge questions remain distinct from the closed radial form theorem.

## 6. Tensor partner bound throughout the same branch

The pure transverse-traceless reduction is
\(h''+4ah'+\mu_T^2\rho^{-2}h=0\).
Two-copy tensor continuity and the sum junction give even Neumann and
odd Dirichlet conditions. With
\[
 dz=dy/\rho,\quad z(T)=0,\quad
 \psi=\rho^{3/2}h,\quad S=\frac32\rho',\quad
 D=\partial_z-S,\quad D^*=-\partial_z-S,
 \tag{36}
\]
the original tensor operator is \(D^*D\). The tensor norm is
\(\int\rho^2|h|^2dy=\int|\psi|^2dz\).
Factorization and this massless profile are standard [7]; a de Sitter
continuum threshold \(3H_b/2\) is also established in the scalar-free
setting [8].

The useful comparison on the nonconstant branch is
\[
 V_-=S^2-S_z
 =\frac94+\rho^2\left(\frac9{16}p^2-\frac18U\right)
 \ge\frac94\qquad(U\le0).                              \tag{37}
\]
For even fields, take \(\operatorname{Dom}D=H^1(-\infty,0]\);
the adjoint has zero shell trace. The original condition is
\(D\psi(0)=0\), and the partner \(DD^*\) is Dirichlet.
For odd fields, reverse these first-order domains. The partner condition
is \(\chi'(0)+S_b\chi(0)=0\), whose form has the positive endpoint term
\(S_b|\chi(0)|^2\), since \(S_b>0\).
Both partner forms are bounded below by \(9/4\).

The polar decomposition of the closed first-order operator compares the
entire nonzero spectra of \(D^*D\) and \(DD^*\), with their correct domains
[13]. The even kernel is the normalizable \(\rho^{3/2}\); the odd
Dirichlet condition removes it. Therefore the even tensor spectrum is
\(\{0\}\cup[9/4,\infty)\), and the odd spectrum is \([9/4,\infty)\)
as spectral sets. Both have essential spectrum \([9/4,\infty)\);
this statement does not need an additional claim about embedded
eigenvalues or threshold spectral type.

The model gate in (37) is exact. Expanding (4),
\[
 U=-\frac2{27}+\frac{14}{9}\eta^2+\frac{50}{27}\eta^3
 -\frac16\eta^4-\frac49\eta^5-\frac2{27}\eta^6.
 \tag{38}
\]
On \(|\eta|\le1/100\), its coefficient upper bound is
\[
 U\le U(0)+\sum_{j\ge1}|U_j|(1/100)^j
 =-\frac{997874976899}{13500000000000}<0.
 \tag{39}
\]
The comparison thus applies to every background in the certified interval.
Nonzero physical tensor masses satisfy \(m_T\ge3H_b/2\).
This is a reduced tensor conclusion, not a complete ghost analysis of
all gravitational and shell degrees of freedom.

## 7. Evidence, limitations and reproducibility

The mathematical conclusions are full-spectrum statements for the specified
finite-energy radial realizations. The physical conclusions additionally
require the matching conditions, mode reductions, and admissible global
gauge/cone class described above. In particular, neither normalizability
in the chosen radial norm nor a local pure-gauge calculation disposes of
externally imposed cone data or every global extension beyond the cone.

The scalar gap formulas (24) and (29) are positive for each fixed
\(\delta>0\). They have not been numerically evaluated or bounded uniformly
in dimensionless units. In physical units the essential threshold is
\((9/4)H_b^2\), and \(H_b^2\to0\) on the branch as \(\delta\to0\).
A strictly positive physical mass floor uniform over that limit is
therefore incompatible with this spectral family. The tensor
dimensionless threshold is exactly controlled, while its physical scale
also follows \(H_b\).

The supplementary domain proof gives the full trace, closure, limit-point,
Hardy and essential-spectrum arguments. Its portable algebra checker
retains ten original identities and adds full spectral-shift and boundary
checks, for thirteen positive checks and four deliberately wrong formulas
with exact nonzero witnesses. Normal and optimized Python produce
byte-identical receipts; injected failures remain failures in both modes.
The tensor checker passes eight identities/sign checks and rejects four
wrong formulas. The completed odd action has eleven exact checks,
including two wrong-formula controls, and a separate covariant
reconstruction. These computations support algebra and exact polynomial
arithmetic. They do not constitute formal verification of functional
analysis or physical completeness. Independent internal reviews inspect
those arguments separately; they are not external peer review.

No broad floating-point spectral scan is used to prove the new inequalities.
No claim is made about vector-sector completeness, nonlinear attraction,
vacuum decay, reheating, or cosmological observables.

### Source and data availability

The October 11 static theorem and exact arithmetic certificate are the
dated background basis. The October 11 Jacobi-to-zero-mode bridge is
identified separately as a dated prospective calculation. The October 13
domain proof, completed operator audit, tensor comparison, literature
review, exact checkers and receipts form the new dated supplementary
record.

The prior static basis is available at the immutable
[October 11 source commit](https://github.com/maldonado-research/HDblast/tree/ca001b58b73ca45117e3849f74ac241489583c66/research/HDBLAST_SCALAR_RESPONSE_MONOTONICITY_20261011).
The capsule contains an unchanged [static basis](../basis/static/README.md),
the dated [prospective bridge](../basis/bridge/PROSPECTIVE_ZERO_MODE_STRATEGY.md),
and [source identities](../basis/SOURCE_IDENTITIES.json).
The continuation's relative [manifest](../MANIFEST.json) identifies its
contents; external release metadata binds the final source commit.

The [domain manifest](../domain/PORTABLE_MANIFEST.json) pins the unchanged
proof and distinguishes original-authoring hashes from adapted CLI hashes.
Historical public and private source files, including the previously
approved manuscript packet, are unchanged.

### Declarations

**Funding:** No funding was received for this work.  
**Competing interests:** The author declares no known competing interests.  
**AI assistance:** Substantial OpenAI Codex assistance was used for
mathematical exploration and drafting, symbolic-check implementation,
literature organization, and internal reviews by separate agent processes.
These reviews do not substitute for human scientific judgment or external
peer review. The author retains responsibility for the scientific
content; this draft remains subject to author review.

## References

1. A. V. Frolov and L. Kofman, “Can Inflating Braneworlds be Stabilized?”
   *Physical Review D* **69**, 044021 (2004).
   [doi:10.1103/PhysRevD.69.044021](https://doi.org/10.1103/PhysRevD.69.044021).
2. D. Konikowska, M. Olechowski and M. G. Schmidt, “Radion stabilization
   with(out) Gauss–Bonnet interactions and inflation,” *Physical Review D*
   **73**, 105018 (2006).
   [doi:10.1103/PhysRevD.73.105018](https://doi.org/10.1103/PhysRevD.73.105018).
3. J. Lesgourgues and L. Sorbo, “Goldberger–Wise variations: stabilizing
   brane models with a bulk scalar,” *Physical Review D* **69**, 084010
   (2004).
   [doi:10.1103/PhysRevD.69.084010](https://doi.org/10.1103/PhysRevD.69.084010).
4. M. Olechowski, “k-stabilization in brane models,” *Physical Review D*
   **78**, 084036 (2008).
   [doi:10.1103/PhysRevD.78.084036](https://doi.org/10.1103/PhysRevD.78.084036).
5. M. Shahriari, A. Jodayree Akbarfam and G. Teschl, “Uniqueness for
   Inverse Sturm–Liouville Problems with a Finite Number of Transmission
   Conditions,” *Journal of Mathematical Analysis and Applications* **395**,
   19–29 (2012).
   [doi:10.1016/j.jmaa.2012.04.048](https://doi.org/10.1016/j.jmaa.2012.04.048).
6. Y. Pinchover and K. Tintarev, “A ground state alternative for singular
   Schrödinger operators,”
   [arXiv:math/0411658v1](https://arxiv.org/abs/math/0411658v1),
   Lemma 2.4, Eq. (2.6). The verified primary record used here is the
   arXiv version; no unverified DOI is supplied.
7. O. DeWolfe, D. Z. Freedman, S. S. Gubser and A. Karch, “Modeling the
   fifth dimension with scalars and gravity,” *Physical Review D* **62**,
   046008 (2000).
   [doi:10.1103/PhysRevD.62.046008](https://doi.org/10.1103/PhysRevD.62.046008).
8. J. Garriga and M. Sasaki, “Brane-world creation and black holes,”
   *Physical Review D* **62**, 043523 (2000).
   [doi:10.1103/PhysRevD.62.043523](https://doi.org/10.1103/PhysRevD.62.043523).
9. L. Kofman, J. Martin and M. Peloso, “Exact identification of the radion
   and its coupling to the observable sector,” *Physical Review D* **70**,
   085015 (2004).
   [doi:10.1103/PhysRevD.70.085015](https://doi.org/10.1103/PhysRevD.70.085015).
10. U. Gen and M. Sasaki, “Radion on the de Sitter brane,”
    *Progress of Theoretical Physics* **105**, 591 (2001).
    [doi:10.1143/PTP.105.591](https://doi.org/10.1143/PTP.105.591).
11. M. Olechowski, “Stability of multibrane models,” *Nuclear Physics B*
    **1011**, 116807 (2025).
    [doi:10.1016/j.nuclphysb.2025.116807](https://doi.org/10.1016/j.nuclphysb.2025.116807).
12. H. Kodama and A. Ishibashi, “A master equation for gravitational
    perturbations of maximally symmetric black holes in higher dimensions,”
    *Progress of Theoretical Physics* **110**, 701 (2003).
    [doi:10.1143/PTP.110.701](https://doi.org/10.1143/PTP.110.701).
13. A. Kostenko, A. Sakhnovich and G. Teschl, “Commutation Methods for
    Schrödinger Operators with Strongly Singular Potentials,”
    *Mathematische Nachrichten* **285**, 392–410 (2012).
    [doi:10.1002/mana.201000108](https://doi.org/10.1002/mana.201000108).
