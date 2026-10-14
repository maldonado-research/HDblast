---
title: A rotational obstruction to nonlinear continuation of Einstein–scalar shell perturbations
author: Ricardo Maldonado
affiliation: Independent Researcher
orcid: 0009-0009-3937-6527
email: ricardomaldonadoandrade@gmail.com
date: 13 October 2026 (America/Los_Angeles) — technical draft for author review
abstract: |
  We derive a quadratic obstruction to nonlinear continuation of actual
  smooth linearized solutions in a five-dimensional Einstein–scalar model
  with two interiors and one dynamical timelike shell. An unrestricted
  boundary variation retains both outward Gibbons–Hawking–York terms and
  the common tension. Its completed Noether cut potential is explicitly
  proportional to the Israel residual, proving cancellation for common
  shell-tangent diffeomorphisms. Moving embeddings are handled by exact
  pullbacks. Differentiation at the background gives ten necessary Lorentz
  quadratic conditions for parameter-C² exact families with closed matched
  Cauchy sections. Independently, an ADM momentum-plus-Israel calculation
  retains arbitrary shell tilt and proves that no second-order shell
  correction can cancel the rotational constraint Hessian. Compact
  degree-two data yield a nonzero Hessian and hence a nonintegrable smooth
  non-gauge linear state, which can also undergo shell contact. Within an
  eighteen-dimensional quadrupole data space with fixed radial profile,
  time scale and reflection parity, all ten necessary tests vanish exactly
  when the displacement and velocity matrices commute. Passing these tests
  does not prove second-order solvability or nonlinear integrability.
---

**Status.** Technical continuation for author review, with companion proofs
and executable algebra checks. This is not an APS submission or external
peer review. The previous linear shell-contact capsule [1] is unchanged.

# 1. Model, hypotheses and result

The first-order construction [1] supplies smooth finite-angular solutions
of the complete Einstein–scalar junction problem, including states that
reach the shell. It does not establish that every such solution is tangent
to an exact nonlinear family. We prove an explicit obstruction to that
extension under a specified classical regularity assumption.

Use signature $(-++++)$, $8\pi G_5=1$, outward spacelike normals $n_i$, and
$K_{i,ab}=h_a{}^Ah_b{}^B\nabla_A n_{i,B}$. The action is
\[
 S=\sum_{i=\pm}\int_{M_i}\epsilon_i
 \left[\frac{R_i}{2}-\frac{(\nabla\phi_i)^2}{2}-U(\phi_i)\right]
 +\int_W\epsilon_W\,[K_++K_--\sigma(\varphi)].             \tag{1}
\]
There are two gravitational boundary completions and **one** shell tension.
The common induced fields obey the exact configuration requirements
\[
 X_+^*g_+=X_-^*g_-=h,\qquad
 X_+^*\phi_+=X_-^*\phi_-=\varphi.                         \tag{2}
\]
The fixed functions are inherited from [1]: with $z=\phi-1$,
\[
 \begin{gathered}
 W=\tfrac13+z^2+\tfrac13z^3,\qquad
 U=\tfrac12W_{,\phi}^2-\tfrac23W^2,\\
 \sigma=2W+\mathfrak d(1+c_*+c_*z),\quad
 c_*=\frac2{1.0357712571566784}-\frac43,\quad
 0<\mathfrak d\le0.002 .
 \end{gathered}                                           \tag{3}
\]
The decimal is exact; fix one background on the registered regular branch,
without a new branch certification or a uniform zero-detuning limit. It has
\[
 g_0=\Omega(s)^2\eta,\quad\phi_0=\Phi(s),\quad
 s=R^2-T^2,\quad f=\Phi_s\ne0,\quad\Omega>0,
 \qquad W:\ s=b^2.                                      \tag{4}
\]
All coefficients are the actual background functions on the controlled
interval of [1]. Its shell hypotheses include $a_b=\sigma_b/6>0$,
$p_b=-\sigma'_b/2\ne0$ and $B_b=p'_b/p_b+\sigma''_b/2>0$,
where $a=\rho'/\rho$, $p=\phi'$ in the radial proper-distance chart.

The complete Cauchy section comprises two compact four-balls glued at their
common closed $S^3$ cut. No exterior boundary or cone excision is present;
the smooth cone and its vertex remain internal. A hypothetical exact family
is assumed twice continuously differentiable in its parameter in a
piecewise $C^k$, $k\ge4$, norm on a compact reference slab, with smooth
timelike embeddings admitting collars, uniformly nondegenerate metrics and
uniformly spacelike Cauchy sections. Pullbacks preserve (2) throughout the
family. All bulk equations and full metric/scalar junctions hold. These are
sufficient classical assumptions for the products, traces and derivatives
below; no solution manifold or weaker-space extension is presumed.

Under these hypotheses, common shell-tangent diffeomorphisms are null
directions of the completed presymplectic form. The ten background Lorentz
symmetries consequently impose necessary quadratic conditions. We exhibit
a smooth, non-gauge linear solution violating a rotational condition,
allowing arbitrary second-order corrections in every bulk and shell sector.
Within the restricted quadrupole data space of Section 6, all ten conditions
vanish if and only if two real symmetric matrices commute.

# 2. Unrestricted boundary variation and cut cancellation

Fix the bulk Lee–Wald potential of [1]. For $q_{AB}=\delta g_{AB}$,
$q=g^{AB}q_{AB}$ and $q^{AB}=g^{AC}g^{BD}q_{CD}$,
\[
 \theta^A=\sqrt{-g}\left[
 \tfrac12(\nabla_Bq^{AB}-\nabla^Aq)-\nabla^A\phi\,\delta\phi\right],
 \qquad \omega=\delta\theta .                            \tag{5}
\]
Each outward boundary form is represented in a common co-oriented density
convention on $W$, denoted $\epsilon_W$. The conversion is applied to
$\theta_i,Q_i,\beta_i$ as well as GHY. Using raw oppositely oriented
pullbacks together with unsigned Israel sums would be inconsistent.

Put $p_i=n_i\phi_i$ and
$v_i^a=h^{ab}n_i^Aq_{i,Ab}$. The unrestricted GHY identity is
\[
 \theta_i|_W+\delta(K_i\epsilon_W)
 =-\epsilon_W\left[\tfrac12(K_i^{ab}-K_ih^{ab})\delta h_{ab}
                 +p_i\delta\varphi\right]+d_W\beta_i,
 \quad \beta_i=-\tfrac12 i_{v_i}\epsilon_W.               \tag{6}
\]
Normal and mixed metric variations have not been gauge-fixed away. Define
\[
 \begin{gathered}
 I_{ab}=\sum_i(K_{i,ab}-K_i h_{ab})+\sigma h_{ab},\qquad
 P=p_++p_-+\sigma',\\
 \ell=(K_++K_--\sigma)\epsilon_W,\quad\beta=\sum_i\beta_i,
 \quad\alpha=-\epsilon_W[\tfrac12I^{ab}\delta h_{ab}+P\delta\varphi].
 \end{gathered}                                           \tag{7}
\]
Both sides share the variations in (7) because of (2). Summing (6) and
varying the single tension gives $\sum_i\theta_i|_W+\delta\ell
=\alpha+d_W\beta$.

Fix reference bulk vectors $\xi_i$ tangent to $W$ and inducing the same
vector $\xi$ there. Their normal derivatives are arbitrary. Let
$G_\xi\Phi_i=\mathcal L_{\xi_i}\Phi_i$ and
$Q_{\xi,i}=-\tfrac12*_i d\xi_i^\flat$. Minimal scalar coupling adds no
separate Komar potential. The essential local cancellation is
\[
 B_\xi:=\sum_iQ_{\xi,i}|_W+i_\xi\ell-\beta(G_\xi)
       =-i_{I\cdot\xi}\epsilon_W .                       \tag{8}
\]
Indeed, if $A_i^a=h^a{}_B n_i^A\nabla_A\xi_i^B$, then
$Q_{\xi,i}|_W=-\tfrac12i_{A_i+K_i\cdot\xi}\epsilon_W$ and
$\beta_i(G_\xi)=-\tfrac12i_{A_i-K_i\cdot\xi}\epsilon_W$.
Subtracting cancels the normal extension $A_i$; the remaining GHY and
tension terms give precisely $I$. Thus (8) is a residual identity, not an
assumption that the boundary charge vanishes.

For $k_\xi=\sum_i(\delta Q_{\xi,i}-i_{\xi_i}\theta_i)|_W$,
covariance, fixed parameters and Cartan's identity yield
\[
 k_\xi-(\delta\beta)(\delta,G_\xi)
 =\delta B_\xi-i_\xi\alpha+d_W(i_\xi\beta(\delta)).       \tag{9}
\]
The final term integrates to zero on the closed cut $S$. The timelike
boundary's future-end orientation is opposite to $\partial C_i$, giving
\[
 \Theta_C^{\rm rel}=\sum_i\int_{C_i}\theta_i-\int_S\beta,
 \qquad\Omega^{\rm rel}=\delta\Theta_C^{\rm rel}.         \tag{10}
\]
The corner is included once. Write $E_{i,AB}=G_{i,AB}-T_{i,AB}$,
$C_{\xi,i}=E_i{}^A{}_B\xi_i^B\epsilon_{i,A}$ and
$\boldsymbol E_i\delta\Phi_i=\epsilon_i[-E_i^{AB}\delta g_{i,AB}/2
+(\Box\phi_i-U')\delta\phi_i]$. The bulk identity and (9) give
\[
 \begin{gathered}
 H_\xi=\sum_i\int_{C_i}C_{\xi,i}+\int_S B_\xi,\\
 \Omega^{\rm rel}(\delta,G_\xi)
 =\delta H_\xi+\sum_i\int_{C_i}i_{\xi_i}
           (\boldsymbol E_i\delta\Phi_i)-\int_S i_\xi\alpha .
 \end{gathered}                                           \tag{11}
\]
On a full solution and a full tangent solution, all residuals and their
linearizations vanish. Hence $\Omega^{\rm rel}(\delta,G_\xi)=0$.
Curling the boundary variation similarly gives conservation of (10).

Embedding consistency is also nonlinear. Gauss–Codazzi and common induced
fields imply, with $\Delta K_{ab}=K_{+,ab}-K_{-,ab}$ and
$\Delta p=p_+-p_-$,
\[
 \begin{gathered}
 \sum_iE_{i,n_i a}=D_bI^b{}_a-PD_a\varphi,\\
 \sigma\Delta K+\sigma'\Delta p
 =2(E_{+,n_+n_+}-E_{-,n_-n_-})+I_{ab}\Delta K^{ab}+P\Delta p .
 \end{gathered}                                           \tag{12}
\]
Thus tangential conservation and normal force are retained residual
relations, rather than extra equations dropped using a first-order identity.

# 3. Moving embeddings and the ten necessary conditions

Straighten each moving region by $Y_i:\bar M_i\to M_i$, with
$Y_i|_W=X_i$, and set $\chi_i=\delta Y_i\circ Y_i^{-1}$.
The exact action equals its barred version, and
\[
 \delta\bar\Phi_i=Y_i^*(\delta\Phi_i+\mathcal L_{\chi_i}\Phi_i),
 \quad
 \Theta_C^{\rm rel}
 =\sum_i\int_{Y_i(\bar C_i)}(\theta_i+i_{\chi_i}L_i+C_{\chi_i,i})
   +\int_S\sum_iQ_{\chi_i,i}-\int_{\bar S}\bar\beta .     \tag{13}
\]
Pullbacks and the common orientation conversion are understood in the cut
terms. The $i_\chi L$ term records moving domains; it generally survives
with matter. The $Q_\chi$ and $\bar\beta$ terms retain the gravitational
cut. $C_\chi$ is removed only on shell, and from the curl only with the
linearized constraints. Definition (10) in barred variables fixes that curl
without omitting moving-location or normal-normalization variations.
Twice differentiating (2) includes bending and evaluation terms: for example,
$\partial_\varepsilon^2(Y_\varepsilon^*g_\varepsilon)
=Y_\varepsilon^*(\ddot g+2\mathcal L_\chi\dot g+
\mathcal L_{\dot\chi}g+\mathcal L_\chi^2g)$.

Since $f\ne0$, a common metric–scalar Killing field satisfies $Ks=0$ and
then $\mathcal L_K\eta=0$. Flat translations fail $Ks=0$ on an open set;
the remaining algebra is exactly the six rotations and four boosts of
$\mathfrak{so}(1,4)$. The shell requires the same action on both copies.
This counts background symmetries, not every possible adjoint kernel in an
unspecified constraint domain.

Use the fixed compensator $k_0=K|_W$:
\[
 G_{K,k_0}X_i=-K\circ X_i+(X_i)_*k_0,
 \qquad G_{K,k_0}(\Phi_0)=0 .                            \tag{14}
\]
Extend $k_0$ by fixed shell-tangent reference vectors $\xi_i$, equal to the
background $K$ when $Y_{i0}=\mathrm{id}$. The compensated transformation
then gives $G_{K,k_0}\bar\Phi_i=\mathcal L_{\xi_i}\bar\Phi_i$, so (11)
applies to the moving system through its barred action.
Its linear action $\mathcal K=DG_{K,k_0}|_0$ is the Lie derivative on the
reconstructed fields and their normal graphs. For an exact family
$\Phi_\varepsilon=\Phi_0+\varepsilon v+\varepsilon^2w/2+o(\varepsilon^2)$,
differentiate the proved degeneracy along that family:
\[
 \Omega^{\rm rel}_0(v,\mathcal K v)=0,
 \qquad D^2H_K|_0[v,v]=\Omega^{\rm rel}_0(v,\mathcal K v). \tag{15}
\]
The Hessian uses a chart on the exact matching space. Equation (11) gives
$DH_K|_0[w]=0$ for arbitrary matching variations; derivatives of its
residual terms vanish for $v$. Thus (15) also follows by differentiating
the exact constraint. A Taub functional defined as the Taylor coefficient
is half this Hessian.

Our relative convention is canonical: $\Omega_{\rm future}=-\Omega^{\rm rel}$.
For the compact initial-support states below, choose the first variation of
the straightening map zero near the **entire initial Cauchy section**,
including its required time jet. All extra terms in (13) then drop from
the initial quadratic pairing, and $\Omega_{\rm future}=\Omega_{\rm phys}$
of [1]. Vanishing merely at the cut would not remove bulk $i_\chi L$ terms.
No new component-by-component identification of arbitrary after-contact cut
currents is claimed; the established total-current conservation of [1]
already transports its initial pairing through contact.

# 4. Independent ADM derivation, with arbitrary shell tilt

The rotational obstruction admits a separate constraint proof. Let $\gamma$
be the spatial metric, $u$ its future unit normal and $k$ its extrinsic
curvature. The normalization in (1) gives
\[
 \pi^{ab}=\frac{\sqrt\gamma}{2}(k^{ab}-k\gamma^{ab}),\quad
 \varpi=\sqrt\gamma\,u\phi,\quad
 \mathcal H_a=-2D_b\pi^b{}_a+\varpi D_a\phi
             =-\sqrt\gamma\,E_{ua}.                     \tag{16}
\]
Let $s_i$ be the outward spatial normal to the cut, $\tau$ the common
future normal to it within the shell, and $c$ its common metric. Define
the boost angle $\eta_i$ by
\[
 n_i=\cosh\eta_i\,s_i+\sinh\eta_i\,u_i,
 \qquad\tau=\cosh\eta_i\,u_i+\sinh\eta_i\,s_i .         \tag{17}
\]
For a fixed common cut-tangent rotation $J$, differentiating this normal
frame along $J$ gives $k_i(s_i,J)=-K_i(\tau,J)-J\eta_i$.
Also $I(\tau,J)=\sum_iK_i(\tau,J)$: the trace and pure-tension terms have
zero mixed projection. Integration by parts therefore derives the exact
functional
\[
 \begin{aligned}
 \mathcal C_J
 &:=\sum_i\int_{C_i}\mathcal H_{i,a}J^a
        -\int_S\sqrt c\,I(\tau,J)\\
 &=M_J+\int_S\sqrt c\,J(\eta_++\eta_-),\qquad
 M_J=\sum_i\int_{C_i}
       (\pi_i^{ab}\mathcal L_J\gamma_{i,ab}+\varpi_iJ\phi_i).
 \end{aligned}                                           \tag{18}
\]
Every exact constraint-and-junction datum has $\mathcal C_J=0$. The angle
term is essential for a tilted cut. Scalar-dependent tension has no
independent intrinsic kinetic momentum or extra mixed angular stress.

At $T=0$, $\pi_0=\varpi_0=\eta_{i0}=0$ and
$\mathcal L_J\gamma_0=J\phi_0=J\sqrt{c_0}=0$. For any second-order
correction $w$, even one nonzero on the cut,
\[
 D\mathcal C_J|_0[w]
 =\int_S\sqrt{c_0}\,J\eta_{\rm sum}[w]=0 .              \tag{19}
\]
The equality uses the closed cut and its invariant background density.
For a first-order state vanishing on an open shell neighborhood, the corner
Hessian $2\int_S\delta\sqrt c\,J\delta\eta_{\rm sum}$ is zero. Thus
\[
 0=D^2\mathcal C_J[v,v]
 =2\sum_i\int_{C_i}
  [\delta\pi_i^{ab}\mathcal L_J\delta\gamma_{i,ab}
                       +\delta\varpi_iJ\delta\phi_i].   \tag{20}
\]
The support restriction removes this *quadratic* cut term; it does not
remove $w$. Equation (19) is the separate cancellation of its arbitrary
second-order shell contribution. No scalar-sector restriction is imposed
on $w$.

For master initial data $q=F|_0$, $p=F_T|_0$, the inherited reconstruction is
\[
 \delta g_{AB}=2F(s\eta_{AB}-3X_AX_B),\qquad
 \delta\phi=-\frac{3(E+4)F}{2\Omega^2f},\quad E=X^A\partial_A.
                                                               \tag{21}
\]
Direct ADM variation at $T=0$, with $D=r\partial_r$, gives
\[
 \begin{gathered}
 \delta\gamma_{ij}=2q(r^2\delta_{ij}-3x_ix_j),\qquad
 \delta\pi^{ij}=-\frac{3p}{2\Omega}x_ix_j,\\
 \delta\phi=-\frac{3(D+4)q}{2\Omega^2f},\qquad
 \delta\varpi=-\frac{3\Omega(D+5)p}{2f}.
 \end{gathered}                                           \tag{22}
\]
The $D+5$ follows from $[\partial_T,E]=\partial_T$. Substitution in (20)
recovers the actual Cauchy bilinear
\[
 \mathcal B(q,p)=\int dr\,d\Omega_3
 [m_r(D+4)q(D+5)p+n_rqp],\quad
 m_r=\frac{9r^3}{4\Omega f^2},\quad n_r=\frac{6r^7}{\Omega}.\tag{23}
\]
Per copy, [1] has
$\Omega_1(v_1,v_2)=\mathcal B(q_1,p_2)-\mathcal B(q_2,p_1)$.
The canonical $\delta\pi\wedge\delta\gamma+\delta\varpi\wedge\delta\phi$
orientation is its negative. For two copies of fixed reflection parity,
\[
 D^2\mathcal C_J[v,v]=4\mathcal B(Jq,p)
                    =-\Omega_{\rm phys}(v,Jv).          \tag{24}
\]
This independently fixes both the sign and the Hessian factor two.

# 5. A smooth non-gauge state with nonzero obstruction

Choose a real harmonic $H$ with $JH\ne0$, a nonzero real
$Q\in C_c^\infty(I)$, and a fixed real $\omega\ne0$ with inverse-time
units. Take
\[
 q=QH,\qquad p=\omega QJH .                              \tag{25}
\]
Angular integration by parts gives, without assuming symmetry of the radial
bilinear,
\[
 D^2\mathcal C_J[v,v]
 =4\omega\mathcal B_{\rm rad}(Q,Q)\int_{S^3}(JH)^2 .    \tag{26}
\]
Here $\mathcal B_{\rm rad}$ is (23) without angular integration. With
$L=\Omega_s/\Omega+2f_s/f$, compact radial integration gives
\[
 \mathcal B_{\rm rad}(Q,Q)
 =\int_I[m_r\{(DQ)^2+(2+9r^2L)Q^2\}+n_rQ^2]\,dr
 \ge\left(\frac{\pi^2m_*}{w_I^2}-M_*\right)\int_IQ^2dr>0,\tag{27}
\]
where $m_* =\min_{\bar I}(m_rr^2)>0$,
$M_* =\max_{\bar I}m_r(-2-9r^2L)_+$ and $w_I$ is the annulus width.
Choose $I$ sufficiently thin within the admissible collar of [1]. This is
positivity of one radial diagonal, not a positive gravitational energy
statement.

An explicit quadrupole uses $n=x/r$ and
\[
 \begin{gathered}
 H=n_1^2-n_2^2,\quad J=x_1\partial_2-x_2\partial_1,
 \quad JH=-4n_1n_2,\quad J^2H=-4H,\\
 \int H^2=\pi^2/3,\quad\int(JH)^2=4\pi^2/3,\qquad
 \boxed{D^2\mathcal C_J[v,v]
       =\frac{16\pi^2}{3}\omega\mathcal B_{\rm rad}(Q,Q)\ne0.}
 \end{gathered}                                           \tag{28}
\]
Either reflection parity doubles the same quadratic contribution. Equations
(19) and (28) exclude a parameter-$C^2$ exact coupled family with this
first derivative, including families with unrestricted second-order sectors.

These are actual solutions of the inherited finite-angular linear problem,
not formal unconstrained tensors. The compact support satisfies every
initial shell compatibility condition; [1] proves smooth evolution and full
junction sufficiency. In its hemisphere collar coordinate $x$, choose
support $[a_0,d_0]$ and a slab with
$0<a_0<d_0<\tau_*$, $d_0+\tau_*<\pi/2$, and coefficients within the
controlled background interval. The support infimum is $a_0$. The
displacement harmonic $H$ and velocity harmonic $JH$ are orthogonal since
$\int HJH=0$. Their separate boundary activations proved in [1] therefore
cannot cancel as a full angular trace. The obstructed state can be chosen
among smooth contacting solutions; no nonzero jet at the exact support
front is asserted.

For $l=2$, non-gauge character follows from the angular trace-free metric
and mixed components in (21). A scalar gauge covector with angular part
$Z D_aH_l$ would produce $2Z(\operatorname{Hess}H_l)^{\rm TF}$, whose
squared norm is $\tfrac23\lambda_l(\lambda_l-3)\int H_l^2>0$.
The absent trace-free part forces $Z=0$, and the absent mixed part then
forces the remaining scalar gauge coefficients to zero. Other angular
sectors cannot cancel this projection.

# 6. Quadrupole classification of the ten tests

Fix the same $Q,\omega$ and reflection parity. Every real degree-two
harmonic is $H_A(n)=n^TAn$ for a symmetric traceless $4\times4$ matrix.
Consider the eighteen-dimensional data space
\[
 q=QH_A,\qquad p=\omega QH_B,
 \quad A=A^T,\ B=B^T,\ \operatorname{tr}A=\operatorname{tr}B=0.
                                                               \tag{29}
\]
For antisymmetric $M$, put $J_M=(Mn)\cdot\partial_n$ and
$V_3=2\pi^2$. Exact sphere moments give
\[
 \int H_AH_B=\frac{V_3}{12}\operatorname{tr}(AB),\qquad
 J_MH_A=H_{[A,M]},\qquad
 D^2\mathcal C_{J_M}[v,v]
 =\frac{V_3}{3}\omega\mathcal B_{\rm rad}(Q,Q)
                        \operatorname{tr}(B[A,M]).      \tag{30}
\]
For $J_{ij}=X_i\partial_j-X_j\partial_i$, this is
$-(2V_3/3)\omega\mathcal B_{\rm rad}(Q,Q)[A,B]_{ij}$.
Since $[A,B]$ is antisymmetric and the prefactor is nonzero, all six
rotational tests vanish exactly when $[A,B]=0$.

The four boosts vanish for every pair (29). The master equation gives
$F_{TT}|_0=\mathcal A q$, where
$\mathcal A=\Delta_{\mathbb R^4}-2L D-[8L+(8/3)r^2f^2]$.
A boost acts on initial data by
\[
 (q,p)\longmapsto (X_i p,\ \partial_iq+X_i\mathcal A q). \tag{31}
\]
Both entries in (29) have even antipodal parity; (31) has odd parity.
The radial operators and weights in (23) preserve it, so each full-sphere
pairing is zero. The transformed state lies in the enlarged finite-angular
family; (29) itself is not asserted to be Lorentz invariant.

Thus **all ten necessary Lorentz tests vanish within (29) if and only if
$[A,B]=0$**. Commuting real symmetric matrices are simultaneously
orthogonally diagonalizable. The degree-two specialization of the previous
same-real-harmonic span is $A=aC$, $B=bC$, including pure displacement and
pure velocity. It passes
these tests, with nonlinear extension still unproved. Neither this
commuting criterion nor ten zero charges characterizes the full nonlinear
constraint range or guarantees higher-order compatibility.

# 7. Prior work, verification and declarations

The general Noether and linearization-stability mechanisms are established.
Iyer–Wald [2] supplies the bulk Hamiltonian identity; Harlow–Wu [3] derives
the completed boundary potential, but its GR example fixes an external
induced metric. The physical interface here instead varies common induced
fields and cancels the two sides against the single tension. Speranza [4]
provides exact embedding pullback and moving-cut formulas; auxiliary
edge-mode degeneracy alone is not used to prove (8)–(11).

Hájíček–Kijowski [5] is direct physical-shell precedent for induced matching,
the area/boost-angle cut pair and interface constraints. Its erratum is
essential: items 13–14 restore density factors and place the last integral
of its Eqs.(99)–(100) on the external boundary, rather than the internal
shell. Its fluid system and differentiability discussion are not silently
imported as an Einstein–scalar existence theorem. Altas–Tekin [6] explicitly
exhibits the smooth compact-Cauchy Taub obstruction. Classical
Fischer–Marsden, Moncrief and Arms references, distinguished from originals
actually read here, appear in the companion literature register. The recent
field-dependent-generator analysis [7] reinforces the need for additional
parameter variations if the compensator in (14) were field dependent; we
keep it fixed. The contribution here is the explicit application to this
matched action, its independent rotational normalization, and the evaluated
finite-angular obstruction, without a claim of general or worldwide novelty.

The recorded normal and optimized-Python runs each contain **81 positive
algebra checks and 21 wrong-formula controls**: 6/8 for nonlinear residuals
and cut coefficients, 46/9 for the independent rotational calculation, and
29/4 for quadrupole algebra. The replay records identify both execution
modes and their results. These are finite algebra checks, not machine-certified PDE,
functional-analytic or complete differential-form proofs. The linear
existence theorem is inherited from [1]; exact-family differentiability is
an explicit hypothesis. No dynamical growth, positive physical energy,
stability or cosmological-origin inference follows from the obstruction.

**Code and data availability.** Companion material is provided at relative
capsule paths: [nonlinear identity](../derivation/FULL_NONLINEAR_SHELL_NOETHER_IDENTITY.md),
[independent ADM proof](../independent/ROTATIONAL_SECOND_ORDER_SHELL_OBSTRUCTION.md),
[quadrupole proof](../root-audit/QUADRUPOLE_ADMISSIBILITY.md), and
[primary-source comparison](../literature/NONLINEAR_SHELL_NOETHER_PRIOR_ART.md).
The original checkers are
[Noether algebra](../derivation/check_nonlinear_shell_noether.py),
[rotational constraints](../independent/check_rotational_constraints.py), and
[quadrupoles](../root-audit/check_quadrupole_admissibility.py).
All paths are relative to the `HDBLAST_SHELL_NOETHER_20261013` capsule,
continuation `2026.10.13-shell-noether`; its integrity and replay records
identify the packaged inputs and executed checks.

**Author and support.** Ricardo Maldonado, Independent Researcher,
ORCID 0009-0009-3937-6527; correspondence:
ricardomaldonadoandrade@gmail.com. No external funding is reported. The
author reports no known conflicts of interest.

**AI assistance and review.** Substantial OpenAI Codex assistance was used
for derivation, code, literature work, internal review and drafting. Separate
agent processes supplied internal checks; these are not external human peer
review. The author retains responsibility for verification, interpretation
and any submission. No APS submission, external peer review or novelty
certification is claimed.

# References

1. R. Maldonado, *Smooth Einstein–scalar perturbations through timelike
   shell contact*, technical checkpoint, 13 October 2026, frozen
   `HDBLAST_REGULAR_SHELL_CONTACT_20261013` capsule, especially the
   existence, full-junction and physical-current proofs.
2. V. Iyer and R. M. Wald, *Some Properties of Noether Charge and a Proposal
   for Dynamical Black Hole Entropy*, Phys. Rev. D **50**, 846–864 (1994),
   [DOI](https://doi.org/10.1103/PhysRevD.50.846),
   [arXiv:gr-qc/9403028v1](https://arxiv.org/abs/gr-qc/9403028v1).
3. D. Harlow and J.-q. Wu, *Covariant phase space with boundaries*, JHEP
   **10** (2020) 146,
   [DOI](https://doi.org/10.1007/JHEP10(2020)146),
   [arXiv:1906.08616v4](https://arxiv.org/abs/1906.08616v4).
4. A. J. Speranza, *Local phase space and edge modes for
   diffeomorphism-invariant theories*, JHEP **02** (2018) 021,
   [DOI](https://doi.org/10.1007/JHEP02(2018)021),
   [arXiv:1706.05061v2](https://arxiv.org/abs/1706.05061v2).
5. P. Hájíček and J. Kijowski, *Lagrangian and Hamiltonian formalism for
   discontinuous fluid and gravitational field*, Phys. Rev. D **57**,
   914–935 (1998),
   [DOI](https://doi.org/10.1103/PhysRevD.57.914),
   [arXiv:gr-qc/9707020v1](https://arxiv.org/abs/gr-qc/9707020v1);
   erratum, Phys. Rev. D **61**, 129901 (2000),
   [DOI](https://doi.org/10.1103/PhysRevD.61.129901).
6. E. Altas and B. Tekin, *Second Order Perturbation Theory in General
   Relativity: Taub Charges as Integral Constraints*, Phys. Rev. D **99**,
   104078 (2019),
   [DOI](https://doi.org/10.1103/PhysRevD.99.104078),
   [arXiv:1903.11982v2](https://arxiv.org/abs/1903.11982v2).
7. H.-F. Ding, *Conserved charges and the first law of black holes for field
   dependent symmetry generators in the covariant phase space formalism*
   (2026), [arXiv:2609.27454v2](https://arxiv.org/abs/2609.27454v2).
