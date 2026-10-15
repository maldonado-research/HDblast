# A nonzero physical scalar response in the conditional exact shell family

15 October 2026 UTC. Conditional analytic research checkpoint.
The background, exact constraint and finite-time evolution hypotheses below
are retained inputs. No numerical response amplitude, external peer review,
certified novelty or cosmological-origin inference is asserted.

## 1. Result and a specified initial state

Retain the registered five-dimensional Einstein–scalar model, identical
interior copies, one common shell tension, and the background

\[
g_0=dy^2+\rho(y)^2[-dt^2+\cosh^2t\,\kappa_{S^3}],\quad
a=\rho_y/\rho,\quad p=\phi_{0,y}\ne0,\quad
z=\log(b/r),\quad dy=-\rho\,dz.
\tag{1}
\]

The shell is at \(z=0\). Both copies use their outward normals toward it.
The physical laws are \(K_{ab}=\sigma h_{ab}/6\) and
\(n\phi=-\sigma'/2\). In particular
\(a_b=\sigma_b/6\), \(p_b=-\sigma'_b/2\), and the registered regular
branch has \(B_b=p_{y,b}/p_b+\sigma''_b/2>0\).

Fix a smooth actual positive-cone collar \([0,Z]\times S^3\), a finite
time interval \([0,T]\), and the strict causal margin required by the
retained exact-evolution theorem. Choose \(\ell>0\) with \(\ell<T\)
and an enclosing correction annulus whose inward depth endpoints satisfy
\(0<A_0<\ell<2\ell<L\) and \(L+c_*T<Z\), where \(c_*\)
bounds the nearby actual radial characteristic speeds. The entire exact
constraint correction, including its larger support, must satisfy this
margin. Define

\[
q_\ell(z)=\begin{cases}
\exp\!\left(4-\dfrac{\ell^2}{(z-\ell)(2\ell-z)}\right),
&\ell<z<2\ell,\\
0,&\text{otherwise},
\end{cases}\qquad
\alpha(z)=\left(\frac{\rho(z)}{\rho_b}\right)^{3/2}
\frac{p(z)}{p_b}>0.
\tag{2}
\]

This smooth bump has maximum one, support infimum \(\ell\), and is
strictly increasing on \((\ell,3\ell/2)\). Let
\(Y=n_1^2-n_2^2\), so \(\Delta_{S^3}Y=-8Y\). Register the master data

\[
P(0,z,\Omega)=\alpha(z)q_\ell(z)Y,\qquad
P_t(0,z,\Omega)=e^{-z}\alpha(z)q_\ell(z)Y.
\tag{3}
\]

These are genuine coupled, reflection-invariant aligned quadrupole data.
Indeed, in the regular ambient master \(P=r^4F\), at ambient time zero
\(\partial_t=r\partial_{T_{\rm amb}}\). With \(s=r^2\), they are exactly

\[
F|_0=G(s)(x_1^2-x_2^2),\quad
F_{T_{\rm amb}}|_0=K(s)(x_1^2-x_2^2),\quad
G(s)=s^{-3}\alpha(\log(b/\sqrt s))q_\ell(\log(b/\sqrt s)),
\quad K=G/b.
\tag{4}
\]

Both profiles are nonzero and smooth with compact annular support. Their
angular matrices are \(A=B=\operatorname{diag}(1,-1,0,0)\). The full
linearized constraints hold by the retained metric–scalar reconstruction;
they are not imposed on a test scalar alone. The general exact localized
initial-data theorem applies to this specified tangent under its stated
reflection, adjoint-kernel, weighted-inverse and support hypotheses. Its
metric and momentum corrections may occupy a larger enclosing annulus.

**Conditional response theorem.** Apply the retained exact initial-data
and finite-time evolution construction to (3)–(4), with matching completed
initial gauge data. For each fixed \(T_1\) with \(\ell<T_1<T\), a
physical-shell response functional defined in section 4 has

\[
\mathcal R_\varepsilon(\rho_bT_1)
=\varepsilon^2\mathcal C_{T_1}+o(\varepsilon^2),\qquad
\mathcal C_{T_1}>0.
\tag{5}
\]

Consequently \(\mathcal R_\varepsilon\ge
\mathcal C_{T_1}\varepsilon^2/2>0\) for all sufficiently small nonzero
\(\varepsilon\). The coefficient is independent of the second initial
correction among exact curves with this same first tangent. No numerical
value of \(\mathcal C_{T_1}\) or effective amplitude radius is supplied.
The clock construction uses the strict time margin \(T_1<T\).

## 2. The full coupled linear problem and the physical trace

Write \(P(t,z,\Omega)=V(t,z)Y(\Omega)\). The inherited longitudinal
representative is

\[
h_{yy}=-4\Psi,\quad h_{ya}=0,\quad h_{ab}=2\Psi h_{0,ab},\quad
\Psi=P/\rho^2,\quad \chi=-3P_y/(\rho^2p).
\tag{6}
\]

The retained differential identities verify all linearized Einstein and
scalar equations from

\[
P_{yy}-(2a+2p_y/p)P_y+
\rho^{-2}(\Box_\gamma+4)P-\tfrac23p^2P=0.
\tag{7}
\]

The copy-even \(l=2\) trace-free Israel equation sets the normal shell
displacement to zero in this representative. The remaining full shell
rows give

\[
B_bP_y=\rho_b^{-2}(\Box_\gamma+4)P_b.
\tag{8}
\]

For an arbitrary representative the physical scalar variation is
\(\Phi_1=\chi_b+p_b\zeta\), including the embedding variation.
Thus in the representative (6),

\[
\Phi_1(t,\Omega)=S(t)Y(\Omega),\qquad
S(t)=\frac{3V_z(t,0)}{\rho_b^3p_b}.
\tag{9}
\]

In the actual positive-cone collar, (7)–(8) become

\[
V_{tt}-V_{zz}+3\tanh t\,V_t-C(z)V_z+
[8\operatorname{sech}^2t-4+\tfrac23\rho^2p^2]V=0,
\quad C=\rho(3a+2p_y/p),
\tag{10}
\]
\[
V_{tt}(t,0)+3\tanh t\,V_t(t,0)+
[8\operatorname{sech}^2t-4]V(t,0)=\rho_bB_bV_z(t,0).
\tag{11}
\]

All initial shell jets vanish on an open collar. If \(S\) were zero on
\((0,T_1)\), (9) would give \(V_z(t,0)=0\). Equation (11) and its zero
initial values would then give \(V(t,0)=0\) throughout that interval.
Section 3 proves that simultaneous vanishing is impossible for (3).

The smooth master solution can be constructed on this actual positive-cone
collar. The conjugation in section 3 gives the wave equation (13) and
\(-w_z+\mu w_{tt}+\beta(t)w=0\), where
\(\mu=(\rho_bB_b)^{-1}>0\) and
\(\beta=\mu(35\operatorname{sech}^2t/4-25/4)+C_b/2\).
Put \(j_+=w_t+w_z\), \(j_-=w_t-w_z\), \(q=w(t,0)\), and
\(v=q'\). Their characteristic and boundary equations are

\[
(\partial_t-\partial_z)j_+=-\mathcal Ww,\quad
(\partial_t+\partial_z)j_-=-\mathcal Ww,\quad
\mu v'+v+\beta q=j_+(t,0),\quad q'=v,\quad
j_-(t,0)=2v-j_+(t,0).
\tag{11a}
\]

Integrate along characteristics, integrate the boundary ODE with its zero
initial jets, and impose \(w=w_0+\frac12\int_0^t(j_++j_-)ds\).
Eliminating \(j_\pm,q,v\) gives a short-slab Volterra contraction for
\(w\): a change in \(w\) changes \(j_+\) by at most
\(M\Delta t\|\delta w\|_\infty\), the boundary ODE introduces a
further time integral, and the final integral defining \(w\) supplies
the contraction factor. Its length depends on the bounded coefficients
and \(\mu^{-1}\). Finite stepping supplies the solution on the whole
fixed slab; compatibility and differentiated characteristic equations
give smoothness. The identities defining the derivatives propagate from
their compatible initial values. Domain of dependence and the strict
support margin remove any auxiliary inner coefficient extension. This
constructs the needed finite harmonic solution without an actual
negative-cone background or a multidimensional timelike Cauchy theorem.

To identify this physical solution with the first derivative of the exact
family, use the retained explicit harmonic gauge bridge. For the vector
covector field \(\xi\), solve

\[
(\Box_{g_0}+\operatorname{Ric}_{g_0})\xi=-C_{\rm GH}[h],\quad
\xi_y|_b=0,\quad(\partial_y-2a_b)\xi_a|_b=0,\quad
\xi|_{t=0}=\partial_t\xi|_{t=0}=0.
\tag{12}
\]

The retained variable mixed linear construction supplies this compatible
finite-regularity solve. Then \(h^{\rm GH}=h+\mathcal L_\xi g_0\) is
relative harmonic and has \(h^{\rm GH}_{ya}=0\) at the shell. The gauge
vector is tangent there, so \(\chi_b^{\rm GH}=\chi_b\). Initial canonical
data and lapse/shift values are unchanged; their time velocities are the
ones fixed by harmonic completion. Uniqueness of the completed linear
harmonic problem identifies this perturbation with the exact-family first
variation. This uses the full fixed-gauge problem, not master uniqueness
as a substitute for metric uniqueness.

## 3. Activation with simultaneous displacement and velocity

One checks \(\alpha_z/\alpha=-C/2\). Setting
\(V=\alpha(z)\cosh^{-3/2}t\,w\) converts (10) to

\[
w_{tt}-w_{zz}+\mathcal W(t,z)w=0,\quad
\mathcal W=\tfrac{35}4\operatorname{sech}^2t-\tfrac{25}4
+\tfrac12C_z+\tfrac14C^2+\tfrac23\rho^2p^2,
\quad w|_0=q_\ell,\quad w_t|_0=e^{-z}q_\ell.
\tag{13}
\]

Suppose \(S\equiv0\) on \((0,T_1)\). The preceding argument gives
\(w=w_z=0\) at \(z=0\). Extend \(w\) by zero to \(z<0\), and extend
the smooth potential boundedly across that auxiliary region. There are no
distributional interface sources because both the value and normal
derivative vanish. Local uniqueness identifies the extension with the
whole-line Cauchy solution on the relevant backward cone. No actual
background outside the positive-cone collar is being asserted.

Take \(t=\ell+\eta<T_1\), \(0<\eta<\ell/2\), and put
\(Q=q_\ell(t)>0\). On the backward diamond of \((t,0)\), initial data
sample only \([\ell,t]\), where \(q_\ell\le Q\) and
\(0\le e^{-z}q_\ell\le Q\). Finite propagation confines the potentially
nonzero part of each source slice to width at most \(\eta\). Let
\(M\) bound \(|\mathcal W|\) on the physical positive-depth part of
the fixed compact cone and \(U_*\) be the supremum of
\(|w|\) on the diamond. The one-dimensional Duhamel formula gives

\[
U_*\le Q(1+\eta/2)+\beta U_*,\qquad
\beta=Mt\eta/2.
\tag{14}
\]

The free observed value is at least \(Q/2\); the source correction has
absolute value at most \(\beta U_*\). Hence, for \(\beta<1\),

\[
w(t,0)\ge
\frac{Q[1-\beta(3+\eta)]}{2(1-\beta)}>0
\quad\hbox{if }Mt\eta(3+\eta)<2.
\tag{15}
\]

Such \(\eta\) exists in every proposed interval beyond \(\ell\). For
example, first bound the potential on the fixed compact cone and choose

\[
0<\eta\le\min\left\{\ell/4,(T_1-\ell)/2,1,
\frac1{8(1+M)(\ell+1)}\right\}.
\tag{16}
\]

This contradicts the assumed zero extension. Therefore \(S\) is nonzero
somewhere in every \((\ell,\ell+\delta)\subset(\ell,T)\). It vanishes
for \(0\le t\le\ell\) by finite propagation and initially zero shell
data. Smooth support edges have all jets zero at \(\ell\); activation
does not mean a nonzero derivative of finite order at the front.

The positive lower bound (15) belongs to the auxiliary whole-line solution
under the contradictory zero-trace hypothesis. It is not a numerical or
pointwise lower bound for the actual Wentzell shell response. Independent
pure-displacement activation via time reflection and sideways uniqueness
is recorded separately; it is a different legitimate tangent.

## 4. An exact intrinsic-clock observable

Let \(\iota_\varepsilon:\Sigma_\varepsilon\to M_\varepsilon\) be the
physical shell and \(h_\varepsilon=\iota_\varepsilon^*g_\varepsilon\).
From its specified initial spatial section, follow future-directed unit
timelike geodesics of \(h_\varepsilon\) initially normal to that section.
Write their proper time as \(\tau\), and their initial labels as
\(\Omega\in S^3\). Denote this intrinsic flow by \(G_\varepsilon\).

At the background, \(G_0(\tau,\Omega)=(\tau/\rho_b,\Omega)\).
Its compact-section Jacobi factors are \(\cosh(\tau/\rho_b)>0\), and
the sections are spacelike. On \(0\le\tau\le\rho_bT_1\), the strict
time margin and small \(C^2\) perturbation preserve existence, local
invertibility and spacelike sections. The maps remain embeddings on the
compact slab by stability of a smooth compact embedding; equivalently one
can use closeness to the background flow and degree one for the spatial
maps. Shrink the parameter interval accordingly. The geodesic ODE and the
retained \(C^2\) parameter family into a piecewise \(C^4\) norm give the
parameter differentiability needed for the following pullbacks and
integrals. No arbitrary global clock or caustic-free infinite evolution
is asserted.

Define \(\widehat\Phi_\varepsilon=G_\varepsilon^*
\iota_\varepsilon^*\phi_\varepsilon\), and let
\(c_\varepsilon(\tau)\) be the positive induced metric of each clock
section. The background value is the constant \(\phi_b\). Set

\[
Q_\varepsilon(\tau)=\int_{S^3}
(\widehat\Phi_\varepsilon(\tau,\Omega)-\phi_b)^2
\,d\mu_{c_\varepsilon(\tau)},\qquad
\mathcal R_\varepsilon(\tau_1)=\int_0^{\tau_1}Q_\varepsilon(\tau)\,d\tau.
\tag{17}
\]

These definitions use the physical induced fields, fixed initial section
and intrinsic proper time. They are unchanged by bulk gauge transformations
with the corresponding embedding change, and by shell relabelings preserving
the specified initial section. Their normalization measures a scalar response,
not gravitational energy.

Because the background shell scalar is constant, every first clock or
label shift multiplies a zero tangential derivative of that background.
Consequently

\[
\partial_\varepsilon\widehat\Phi_\varepsilon|_0
=S(\tau/\rho_b)Y(\Omega).
\tag{18}
\]

For a general bulk representative, the normal embedding term
\(p_b\zeta\) in section 2 remains essential. It is absent from (18)
only after use of the straight-shell representative or its tangent gauge
bridge.

Uniform parameter differentiation on the compact slab gives

\[
Q_0=Q'_0=0,\quad
Q''_0(\tau)=\frac{2\pi^2}{3}\rho_b^3
\cosh^3(\tau/\rho_b)S(\tau/\rho_b)^2,
\tag{19}
\]

since \(\int_{S^3}Y^2d\Omega=\pi^2/3\). In particular

\[
\boxed{\mathcal C_{T_1}=\frac{\pi^2}{3}\rho_b^4
\int_0^{T_1}\cosh^3t\,S(t)^2dt>0.}
\tag{20}
\]

Continuity of \(S\) and section 3 prove the strict inequality. Taylor's
formula for the \(C^2\) parameter curve yields (5) and its sufficiently
small-amplitude lower bound. Neither the second scalar variation nor the
first measure/clock variation enters the quadratic coefficient, because
the zeroth scalar difference in (17) is zero. Thus second-order initial
completion affects higher response terms but not (20).

## 5. What remains open

The argument proves qualitative nonvanishing for a specified family, under
the retained exact-data/evolution hypotheses. It does not produce a numerical
trace \(S(t)\), compute \(\mathcal C_{T_1}\), or give a numerical allowed
amplitude interval. A validated time trace and effective bounds remain useful
next calculations. The exact nonlinear first-arrival time is also open:
perturbed cones vary, and the exact initial correction may create
order-\(\varepsilon^2\) data in a larger annulus, potentially giving an
order-\(\varepsilon^4\) response diagnostic before the linear front.

The conclusion concerns the controlled region and finite time only. It
does not establish global stability, observational cosmology, a higher-
dimensional origin of the Big Bang, or publication readiness for an external
journal. Finite symbolic checks below address specified identities and
control formulas; they do not computer-certify the analytic proof.

## 6. Derivations and retained inputs

The independently derived [aligned positive-velocity response](sol/ALIGNED_POSITIVE_VELOCITY_SHELL_RESPONSE.md)
gives the complete characteristic transfer and the near-front argument.
An independent [pure-displacement construction](astra/PHYSICAL_SCALAR_ACTIVATION.md)
uses time reflection and sideways uniqueness for a different allowed state.
The [critical observable audit](referee/INDEPENDENT_OBSERVABLE_AUDIT.md)
records the physical trace, gauge, clock and completion checks. Pinned
acceptances accompany the frozen source manifest.

The retained [exact finite-regularity evolution theorem](inputs/quasilinear/EXACT_FINITE_REGULARITY_SHELL_EVOLUTION.md),
its [exact initial-data input](inputs/quasilinear/inputs/EXACT_LOCALIZED_INITIAL_DATA_AND_PRECONTACT.md),
and its [gauge/parameter completion](inputs/quasilinear/gauge-parameter/INITIAL_GAUGE_SUBSIDIARY_AND_PARAMETER.md)
state the conditional analytic hypotheses used here. The new profile is an
instance of the general admissible tangent class; it is not a claim that an
earlier unspecified radial bump has the same response function.

The verbatim [regular-shell contact proof](inputs/REGULAR_SHELL_CONTACT.md)
supplies the full coupled reconstruction and physical graph equations. The
[full finite-angular contact proof](inputs/FULL_FINITE_ANGULAR_SECOND_ORDER_CONTACT.md)
supplies the explicit harmonic gauge bridge. Their original Git source
identities and the accessible historical archive reads are recorded in
[the source map](archive/SHELL_RESPONSE_SOURCE_MAP.md) and its provenance
ledger. The new direct collar construction above avoids importing the older
contact proof's additional negative-cone background assumption.

The [bounded primary-literature review](literature/SHELL_RESPONSE_PRIMARY_LITERATURE.md)
distinguishes pullback identities, observability restrictions and other
brane models from the proof for this particular coupled model.
