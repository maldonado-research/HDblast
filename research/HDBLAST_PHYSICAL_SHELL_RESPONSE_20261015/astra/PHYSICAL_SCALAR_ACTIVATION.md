# Physical shell scalar activation for a selected annular quadrupole

15 October 2026. Independent analytic derivation. This note retains the
registered model and the hypotheses of the conditional finite-regularity
shell-evolution theorem. It does not infer response from existence alone.
The selected datum below is a **new pure-displacement member** of the
admissible reflection-invariant class. It is not the previously selected
state with both displacement and velocity nonzero.

## 1. Coupled equations, shell, and the selected tangent

Fix one actual registered background

\[
g_0=dy^2+\rho(y)^2\gamma,\qquad
\gamma=-dt^2+\cosh^2t\,\kappa_{S^3},\qquad
a=\rho_y/\rho,\quad p=\phi_{0,y},
\]

with \(\rho>0,p\ne0\) on the closed controlled positive-cone collar.
The field equations and identical-copy outward-normal shell laws are

\[
R_{AB}-\phi_A\phi_B-\tfrac23U(\phi)g_{AB}=0,\qquad
\Box_g\phi-U'(\phi)=0,
\quad K_{ab}=\tfrac16\sigma(\phi)h_{ab},\quad
n\phi=-\tfrac12\sigma'(\phi).
\tag{1}
\]

Here \(W=1-\phi+\phi^3/3\), \(U=(W')^2/2-2W^2/3\),
\(\sigma=2W+\delta(1+c\phi)\),
\(c=96714321425692/161839258930731\), and \(0<\delta\le1/500\).
The background gives \(a_b=\sigma_b/6\),
\(p_b=-\sigma'_b/2\), and
\(B_b=p_{y,b}/p_b+\sigma''_b/2>0\).
No uniform limit as \(\delta\to0\) is used.

Use \(dr/r=dy/\rho\), \(r_b=b\), and inward depth
\(z=\log(b/r)\). Then \(\partial_y=-\rho^{-1}\partial_z\).
Let \(Y(n)=n_1^2-n_2^2\), so \(\Delta_{S^3}Y=-8Y\).
Choose numbers \(0<a_*<d_*<T\) and a closed collar
\(0\le z\le Z\) with \(Z>d_*+T\). The weighted initial-data
annulus must contain the chosen support strictly; the nonlinear causal
margin must contain the support of its entire constraint correction,
as in the retained exact theorem. Set

\[
h(z)=\begin{cases}
\exp[-1/((z-a_*)(d_*-z))],&a_*<z<d_*,\\
0,&\text{otherwise}.
\end{cases}
\tag{2}
\]

In ambient coordinates \(X=(T_{\rm amb},\mathbf X)\), the actual
master Cauchy data at \(T_{\rm amb}=0\) are

\[
F_0=G(s)(X_1^2-X_2^2),\qquad F_{T_{\rm amb},0}=0,\qquad
G(s)=s^{-3}h\!\left(\log\frac b{\sqrt s}\right),
\quad s=|\mathbf X|^2.
\tag{3}
\]

The expression is extended by zero; its annular support makes it smooth
at the origin. The displacement matrix is
\(A=\operatorname{diag}(1,-1,0,0)\), and the velocity datum is zero.
One could formally set \(B=A\) with profile \(K=0\), but the effective
physical pair is \((A,0)\), not the requested state with both aligned
sectors active. This is explicitly a complementary selected tangent. The
inherited coupled ADM reconstruction solves the linearized constraints
and is invariant under each coordinate reflection. The exact localized
initial-data theorem applies to this tangent, including its prescribed
scalar sectors; the exact shell-evolution proof applies with the same
background, regularity, support, and finite-interval hypotheses. Thus it
constructs a new exact family rather than changing the tangent of an
already fixed family.

## 2. Full physical linear solution and its boundary observation

The retained reconstruction, which is a proved sufficient family of the
complete linearized Einstein--scalar equations, is

\[
P=r^4F,\quad \Psi=P/\rho^2,\quad
\chi=-3P_y/(\rho^2p),\qquad
k_{yy}=-4\Psi,\quad k_{ya}=0,\quad
k_{ab}=2\rho^2\Psi\gamma_{ab}.
\tag{4}
\]

It obeys all bulk equations whenever

\[
P_{yy}-(2a+2p_y/p)P_y+
\rho^{-2}(\Box_\gamma+4)P-\tfrac23p^2P=0.
\tag{5}
\]

This is a coupled gravitational--scalar construction, not a scalar test
field on a frozen metric. The unreduced linear equations are explicitly

\[
\delta R_{AB}-\phi_{0,A}\partial_B\chi-\phi_{0,B}\partial_A\chi
-\tfrac23(U'_0\chi g_{0,AB}+U_0k_{AB})=0,
\]
\[
\Box_0\chi-U''_0\chi-k^{AB}\nabla_A\nabla_B\phi_0
-(\nabla^Ak_A{}^B-\tfrac12\nabla^B\operatorname{tr}k)\partial_B\phi_0=0.
\tag{6}
\]

For identical copies choose the sufficient even embedding
\(\zeta_+=\zeta_-=0\). The retained complete graph calculation gives
the tensor Israel residual
\(-2\rho_b^{-2}(D^aD_b+\delta^a_b)\zeta_e\), and the summed scalar
residual

\[
-\frac6{\rho_b^2p_b}
\left[B_bP_y-\rho_b^{-2}(\Box_\gamma+4)P\right]
+2p_bB_b\zeta_e.
\tag{7}
\]

Consequently every matching and physical junction row is satisfied by
(4)--(5) and

\[
B_bP_y=\rho_b^{-2}(\Box_\gamma+4)P\quad\text{at the shell}.
\tag{8}
\]

The physical first variation of the shell scalar is the pullback

\[
\Phi_1=\chi_b+p_b\zeta_e=\chi_b.
\tag{9}
\]

The zero displacement is a property of this sufficient even representative;
it cannot be imposed in an arbitrary gauge while dropping the displacement
term in (9). The quadrupole avoids the finite-dimensional de Sitter
embedding kernel, but no harmonic division is used in (7).

Write \(P(t,z,n)=V(t,z)Y(n)\). Direct transformation of (5) and (8)
gives the complete radial system

\[
V_{tt}-V_{zz}+3\tanh t\,V_t-C(z)V_z+
\left[8\operatorname{sech}^2t-4+\tfrac23\rho^2p^2\right]V=0,
\quad C=\rho(3a+2p_y/p),
\tag{10}
\]
\[
V_{tt}+3\tanh t\,V_t+
(8\operatorname{sech}^2t-4)V-\rho_bB_bV_z=0\quad(z=0),
\tag{11}
\]
\[
V(0,z)=h(z),\qquad V_t(0,z)=0,\qquad
\Phi_1(t,n)=S(t)Y(n),\quad
S(t)=\frac{3V_z(t,0)}{\rho_b^3p_b}.
\tag{12}
\]

The Cauchy conversion is exact: on the ambient zero-time slice,
\(r=|\mathbf X|\), \(\partial_t=r\partial_{T_{\rm amb}}\), so
\(P_0=r^6G(r^2)Y=hY\) and \(P_t=0\).

For clarity, the finite-harmonic smooth solver used here can be obtained
on this very collar. Put

\[
V=\alpha(z)\cosh^{-3/2}t\,w,\qquad
\alpha(0)=1,\quad\alpha'/\alpha=-C/2.
\]

Then

\[
w_{tt}-w_{zz}+\mathcal W(t,z)w=0,\quad
\mathcal W=\tfrac{35}4\operatorname{sech}^2t-\tfrac{25}4
+\tfrac12C'+\tfrac14C^2+\tfrac23\rho^2p^2,
\tag{13}
\]

and the boundary is

\[
-w_z+\mu w_{tt}+\beta(t)w=0,\qquad
\mu=(\rho_bB_b)^{-1}>0,\quad
\beta=\mu(\tfrac{35}4\operatorname{sech}^2t-\tfrac{25}4)+C_b/2.
\tag{14}
\]

Let \(p_c=w_t+w_z\), \(m_c=w_t-w_z\), \(q=w(t,0)\),
\(v=q'\). Their equations are
\((\partial_t-\partial_z)p_c=-\mathcal Ww\),
\((\partial_t+\partial_z)m_c=-\mathcal Ww\), and
\(\mu v'+v+\beta q=p_c(t,0)\), \(q'=v\),
\(m_c(t,0)=2v-p_c(t,0)\). Integrating along characteristics and
integrating the boundary ODE gives the same linear Volterra contraction
as the retained contact proof. Its step size depends on bounded
coefficients and \(\mu^{-1}\), not amplitude. Compatible smooth data
therefore give a smooth solution on the whole prescribed finite slab by
finitely many steps, with support in
\([\max(0,a_*-|t|),d_*+|t|]\). All corner jets vanish because (2)
vanishes on an open shell collar. The strict depth margin makes any
coefficient extension beyond the collar irrelevant. This argument requires
no continuation of a new perturbation through the cone and no global
spectral expansion.

## 3. Nonvanishing is a separate uniqueness argument

**Activation theorem.** For every \(\eta>0\) with
\(a_*+\eta\le T\), there is \(t_*\in(a_*,a_*+\eta)\) with
\(S(t_*)\ne0\).

The finite propagation and homogeneous boundary ODE give \(S=0\)
for \(0\le t<a_*\). Suppose, in contradiction to the theorem, that
\(S=0\) also throughout \((a_*,a_*+\eta)\). Equation (12) gives
\(V_z(t,0)=0\) for \(0\le t<a_*+\eta\). Equation (11) then
says that \(g(t)=V(t,0)\) solves

\[
g''+3\tanh t\,g'+(8\operatorname{sech}^2t-4)g=0,
\qquad g(0)=g'(0)=0.
\tag{15}
\]

Ordinary ODE uniqueness implies \(g=0\) on this entire interval.
The coefficients of (10)--(11) are time even, except for the time-odd
coefficient multiplying one time derivative. The initial velocity is zero.
Linear uniqueness therefore gives the even extension
\(V(-t,z)=V(t,z)\). Both \(V\) and \(V_z\) vanish at \(z=0\)
for \(|t|<a_*+\eta\).

For this one fixed harmonic, regard (10) as a wave equation with \(z\)
as time and \(t\) as its one spatial coordinate. Its principal part is
\(V_{zz}-V_{tt}\); the remaining coefficients are bounded and smooth
on each closed subtriangle. Zero Cauchy data at \(z=0\) imply the zero
solution in \(z+|t|<a_*+\eta\), by the elementary one-dimensional
hyperbolic energy estimate or its characteristic integral formula. At
\(t=0\) this would give \(h(z)=0\) for
\(0<z<a_*+\eta\), contradicting (2). This proves activation.

This is a sideways uniqueness proof for one fixed radial equation. It
does not assert that an arbitrary multidimensional timelike Cauchy problem
is well posed. Nor does a nonzero master trace alone imply a nonzero
scalar normal derivative: the boundary ODE and zero initial jets are the
essential extra argument. Smooth flat support edges permit every jet at
\(t=a_*\) to vanish, so no nonzero jet or value at that exact instant is
claimed. Cancellation at some later times is allowed.

## 4. Identification with the exact family's first variation

The longitudinal reconstruction need not be relative harmonic. Let
\(C^{(1)}[k]_A=\nabla^Bk_{AB}-\tfrac12\nabla_A\operatorname{tr}k\).
Solve the retained vector wave problem

\[
(\Box_0+\operatorname{Ric}_0)\xi=-C^{(1)}[k],\qquad
\xi_y|_b=0,\quad (\partial_y\xi_a-2a_b\xi_a)|_b=0,
\qquad \xi|_0=\xi_t|_0=0.
\tag{16}
\]

The established linear mixed estimate applies; the source is smooth and
vanishes on the initial shell collar. Then
\(k^w=k+\mathcal L_\xi g_0\),
\(\chi^w=\chi+\xi^y p\) obey the full linearized reduced equations,
\(C^{(1)}=0\), and the fixed-shell coordinate rows. The boundary vector
is tangential and preserves the shell. Its zero initial first jet preserves
the initial canonical geometry. The exact initial harmonic completion can
be chosen to have precisely these full first reduced Cauchy data, as in
the retained theorem. Differentiating the exact family and invoking
linear uniqueness identifies its first variation with \((k^w,\chi^w)\).
At the shell \(\xi^y=0\); hence its physical scalar variation is still
(12). This supplies the bridge between an actual coupled master solution
and the tangent of the exact nonlinear family.

## 5. Intrinsic clock and a positive exact nonlinear observable

Let \(\Phi_\varepsilon=\iota_\varepsilon^*\phi_\varepsilon\)
be the actual common shell scalar. Choose a regular intrinsic clock
\(\tau_\varepsilon\) by proper time along the shell's timelike geodesics
normal to the initial shell section. On any fixed finite interval where
this congruence remains regular, let \(j_{\varepsilon,\tau}:S^3\to\Sigma\)
identify its clock cut with the initial section. At the background this
clock is \(\tau=\rho_b t\) and its cut metric is
\(\rho_b^2\cosh^2t\,\kappa\). The background normal congruence is
regular on a finite interval; after shrinking the exact family parameter,
regularity persists by finite-interval ODE dependence and the nonzero
Jacobian of the background flow. Alternatively regularity of any other
specified intrinsic clock may be stated explicitly as a hypothesis.

Define

\[
Q_\varepsilon(\tau)=
\int_{S^3}
\left[j_{\varepsilon,\tau}^*\Phi_\varepsilon-\phi_b\right]^2
d\mu_{\varepsilon,\tau}.
\tag{17}
\]

The retained \(C^2\) curve into \(\mathcal X_{10}\), with its lower
piecewise \(C^4\) control, is sufficient for the pullback, clock, induced
measure, and twice differentiated integral. The background scalar
\(\phi_b\) is constant on the entire shell, so the first derivative of
the cut map contributes \(\mathcal L_{\dot j}\phi_b=0\).
The embedding contribution \(p_b\zeta\) is already included in (9).
Thus, independently of the regular intrinsic clock's first variation,

\[
Q_0=Q'_0=0,\qquad
Q''_0(\rho_b t)=2\int_{S^3}\Phi_1(t,n)^2\,d\mu_{0,t}.
\tag{18}
\]

Terms involving the first or second variation of the measure or the
second variation of the scalar vanish here because the unperturbed
squared residual and its first derivative are zero. Since
\(\int_{S^3}Y^2d\Omega_3=\pi^2/3\), equation (18) is exactly

\[
\boxed{
Q''_0(\rho_b t)=
\frac{2\pi^2}{3}\rho_b^3\cosh^3t\,S(t)^2
=\frac{6\pi^2\cosh^3t}{\rho_b^3p_b^2}\,V_z(t,0)^2.}
\tag{19}
\]

At every activation witness \(t_*\) this is strictly positive, and

\[
Q_\varepsilon(\rho_b t_*)=
\frac{\varepsilon^2}{2}Q''_0(\rho_b t_*)+o(\varepsilon^2)>0
\]

for all sufficiently small nonzero \(\varepsilon\). This is a concrete
gauge-invariant exact nonlinear shell response under the retained
existence hypotheses. It is the squared total deviation from the fixed
background shell value, not a mean-subtracted variance. There is one
physical common shell, so no additional factor of two from the identical
bulk copies occurs. The unknown numerical response size is not estimated.
The actual exact family may have an earlier support arrival from nonlinear
constraint corrections; no exact finite-amplitude contact time is asserted.

## 6. Provenance and limits of computation

The source equations and full junction sufficiency are retained in
`HDBLAST_SHELL_NOETHER_20261013/basis/prior-checkpoint/junctions/FULL_PACKET_JUNCTION_SUFFICIENCY_AND_CURRENT.md`.
The Volterra construction and pure time-parity activation method are
retained in the adjacent `existence/FINITE_SLAB_SHELL_CONTACT_EXISTENCE.md`.
The boundary-preserving harmonic bridge is
`HDBLAST_FULL_SECOND_ORDER_CONTACT_20261014/FULL_FINITE_ANGULAR_SECOND_ORDER_CONTACT.md`,
section 3. Exact initial data and exact evolution/parameter dependence
are the frozen `HDBLAST_QUASILINEAR_SHELL_EVOLUTION_20261015` capsule.
No prior source is modified. The new step is the physical scalar trace
argument (15) and its connection to the exact observable (19).

The companion symbolic checker independently verifies coordinate
conversion, the master conjugation and boundary ODE, the scalar
observation factor, sphere normalization, and the second variation of
the observable. Deliberately wrong formulas must fail. These finite
checks support specific algebra; they do not establish analytic existence,
uniqueness, a numerical lower bound, or novelty. All scientific claims
remain conditional on the actual registered background and the retained
exact initial-data/evolution hypotheses.
