# Complete geometric quadratic bulk-source adapter

This adapter completes the geometric Ricci contribution missing from the
numerical-response parent. It evaluates a local covariant two-jet of the
registered first metric, adds every natural Einstein/scalar matter term,
and applies the registered relative-harmonic reduction. The output is an
explicit indexed algebraic recipe, with exact coefficient maps exported and
verified by the sibling `projection` package. It is not an evolved second
waveform or an exact nonlinear cosmology.

## Conventions and inputs

The reference metric and frame are

\[
g_0=dy^2+\rho(y)^2[-dt^2+c(t)^2\kappa],\qquad
e_y=\partial_y,\quad e_0=\rho^{-1}\partial_t,\quad
e_a=(\rho c)^{-1}E_a,
\]

where \(c=\cosh t\), \(\theta=\tanh t\), \(a=\rho_y/\rho\),
\(p=\phi_{0,y}\). The parameter expansion is
\(g=g_0+\epsilon q+\epsilon^2k/2\),
\(\phi=\phi_0+\epsilon\eta+\epsilon^2w/2\). Every quadratic
source is a second parameter derivative, not an unscaled coefficient of
\(\epsilon^2\).

The seven metric amplitudes supplied to `evaluate_sources` mean exactly

\[
q_{yy}=AY,\quad q_{y0}=BY,\quad q_{ya}=C D_aY,\quad
q_{00}=DY,\quad q_{0a}=E D_aY,\quad
q_{ab}=LY\kappa_{ab}+M D_aD_bY,\qquad \eta=sY,
\]

with \(Y=n_1^2-n_2^2\). They are the actual relative-harmonic frame
amplitudes, not the raw gauge-vector amplitudes or the longitudinal master
alone. From the parent's first-field transformation,

\[
\begin{aligned}
A&=-4u+2f_y,& B&=f_t/\rho-\rho b_y,&
C&=f/(\rho c)+\rho c d_y,\\
D&=-2(u+af+b_t),&E&=c d_t-b/c,&
L&=2(u+af+\theta b),&M&=2d,&s&=v+pf.
\end{aligned}
\]

For each of `A,B,C,D,E,L,M`, input dictionary `m` contains the suffixes
`'',y,t,yy,yt,tt`. Dictionary `s` contains `s,sy,st,syy,syt,stt`.
These are **ordinary derivatives in physical y,t**, including derivatives
of the frame factors in the displayed amplitudes. They are not covariant
derivatives and not scaled x,r derivatives. The common-gauge master/gauge
three-jet suffices to construct these two-jets.

Background dictionary `b` contains
`rho,c,theta,a,ay,p,py,Vp,Vpp,Vppp,ell`. The adapter uses
\(\rho_y=a\rho\), \(\rho_{yy}=(a_y+a^2)\rho\),
\(c_t=c\theta\), \(c_{tt}=c\), and
\(\theta_t=1-\theta^2\). Its mathematical domain requires nonzero
\(\rho,c\); the registered geometry has both positive. Inputs are
assumed to be coherent jets. The routine does not silently project arbitrary
supplied jets onto a field equation or a harmonic constraint.

## Every frame and derivative slot

The code uses six storage slots: index 0 is y, index 1 is the timelike frame
vector, and indices 2 through 5 are four redundant ambient tangent slots
for the three-sphere. This is a five-dimensional geometry, not a
six-dimensional one. The tangent projector is
\(\Pi_{ij}=\delta_{ij}-n_i n_j\), with \(\operatorname{tr}\Pi=3\).
All sphere slots are tangent. Ambient Euclidean contraction on those slots
therefore equals intrinsic contraction. The inverse-frame contraction signs
are `(1,-1,1,1,1,1)`; the background covariant metric's sphere block is
\(\Pi\), not the four-dimensional identity.

The angular ring reduces polynomial expressions by
\(n_4^2=1-n_1^2-n_2^2-n_3^2\). For any tangent covariant tensor
\(T_J\), its intrinsic derivative is implemented by the exact identity

\[
D_i T_J=\Pi_{ik}\partial_k T_J+
\sum_{s:\,J_s\text{ angular}}n_{J_s}\,T_{J[J_s\mapsto i]}.
\]

There is one connection contribution for **each** angular slot. When
\(T=\nabla q\), the first slot is the derivative slot and is included.
This implements the sphere curvature commutators; ordinary component
partial derivatives would not suffice.

Writing \(h=\theta/\rho\), the remaining connection is

\[
\nabla_y e_A=0,\quad
\nabla_0e_y=ae_0,\quad\nabla_0e_0=ae_y,\quad
\nabla_i e_y=ae_i,\quad\nabla_i e_0=he_i,
\]
\[
\nabla_i e_j=(\rho c)^{-1}\omega^k{}_{ij}e_k
                 +h\Pi_{ij}e_0-a\Pi_{ij}e_y.
\]

Consequently the derivative on a lower slot contributes:

| Direction | Slot y | Slot 0 | Slot j on sphere |
|---|---|---|---|
| y | 0 | 0 | 0 |
| 0 | \(-aT_{y\mapsto0}\) | \(-aT_{0\mapsto y}\) | 0 |
| i on sphere | \(-aT_{y\mapsto i}\) | \(-hT_{0\mapsto i}\) | \(+a\Pi_{ij}T_{j\mapsto y}-h\Pi_{ij}T_{j\mapsto0}\) |

These act in addition to the component derivative
\(\partial_y\), \(\rho^{-1}\partial_t\), or
\((\rho c)^{-1}D_i\). The table acts on **every** lower index,
including every previous derivative index. `frame_two_jet` builds
`dq[d,a,b]=nabla_d q_ab` and
`ddq[e,d,a,b]=nabla_e nabla_d q_ab` this way. A small ordinary two-variable
jet algebra differentiates amplitudes and frame coefficients using exact
Leibniz rules. Only first derivatives of the connection coefficients are
needed; unused second derivatives of a and theta are never consumed.

## Indexed Ricci Hessian

The reusable `indexed_geometry` function consumes arbitrary covariant
metric two-jets and a diagonal inverse-frame signature. It does not depend
on the quadrupole, sphere polynomial package, or first-field equations.
It forms the entirely lower-index tensors

\[
P_{ABC}=\tfrac12(\nabla_Bq_{CA}+\nabla_Cq_{BA}-\nabla_Aq_{BC}),
\qquad J_{ABC}=-2q_A{}^D P_{DBC}.
\]

Metric compatibility gives

\[
\nabla_E P_{ABC}=\tfrac12(
\nabla_E\nabla_Bq_{CA}+\nabla_E\nabla_Cq_{BA}
-\nabla_E\nabla_Aq_{BC}),
\]
\[
\nabla_EJ_{ABC}=-2[(\nabla_Eq_A{}^D)P_{DBC}
                              +q_A{}^D\nabla_EP_{DBC}].
\]

Along the affine metric path, the connection difference is
\(\Gamma(g_\epsilon)-\Gamma(g_0)=\epsilon P+
\epsilon^2J/2+O(\epsilon^3)\). The exact curvature difference identity
therefore yields

\[
Q^{\rm Ric}_{AB}=\nabla_CJ^C{}_{AB}-\nabla_BJ^C{}_{AC}
 +2(P^C{}_{CD}P^D{}_{AB}-P^C{}_{BD}P^D{}_{AC}).
\]

No extra background-curvature term is omitted: the background connection
is already present in these tensorial covariant derivatives. Commuting
second derivatives into a different displayed order would introduce
curvature terms; this evaluator retains their stated order.

## Natural sources, lowering and reduction

For the residuals
\(\mathscr E_{AB}=R_{AB}-\phi_A\phi_B-\frac23Vg_{AB}\),
\(\mathscr S=\Box\phi-V'\), the complete natural Einstein Hessian is

\[
Q^E_{AB}=Q^{\rm Ric}_{AB}-2\eta_A\eta_B
-\tfrac23(V''\eta^2g_{0,AB}+2V'\eta q_{AB}).
\]

The scalar Hessian in the actual harmonic first field is

\[
Q^\phi=2q^A{}_Dq^{DB}(\nabla^2\phi_0)_{AB}
-2q^{AB}(\nabla^2\eta)_{AB}
+2q^{AB}P^C{}_{AB}(\phi_0)_C-V'''\eta^2.
\]

Here \(\nabla^2\phi_0\) has yy entry \(p_y\), 00 entry
\(-ap\), and sphere block \(ap\Pi\). `Qphi_general` retains the
additional terms \(2q^C{}_D C_1^D(\phi_0)_C-2C_1^C\eta_C\).

Let \(C_1^A=g_0^{BC}P^A{}_{BC}\). Three related second derivatives
must be distinguished off harmonicity:

\[
Q_{C,\mathrm{up}}^A=(C^A(g_\epsilon))''
=-2q^A{}_D C_1^D-2q^{BC}P^A{}_{BC},
\]
\[
Q_{C,\mathrm{background},A}=g_{0,AB}Q_{C,\mathrm{up}}^B,
\qquad
Q_{C,\mathrm{actual},A}=(g_{\epsilon,AB}C^B(g_\epsilon))''
=-2q^{BC}P_{ABC}.
\]

The \(qC_1\) terms cancel in the last expression because differentiating
the lowering metric adds \(2q_{AB}C_1^B\). Outputs are named `QCup`,
`QC_background`, and `QC_actual` respectively. `QC` aliases `QC_actual`.
Raising a time component changes its sign. The output `C1` is the
**lower** first gauge covector.

For the full variable-metric reduction
\(\mathscr E^{\rm red}=\mathscr E-\nabla^g_{(A}C_{B)}\),
the general affine Hessian is

\[
Q^{E,\mathrm{red}}_{AB}=Q^E_{AB}
-\nabla^0_{(A}Q_{C,\mathrm{actual},B)}
+2P^D{}_{AB}C_{1,D}.
\]

Parentheses contain one half. On the registered actual first harmonic
field \(C_1=0\), this reduces to the requested
\(Q^E-\nabla^0_{(A}Q_{C,B)}\). Default `harmonic=True` returns that
registered specialization; it is an assumption justified by the prior
first gauge solve, not a fact inferred from arbitrary supplied interval
jets. `QEred_general` always returns the unspecialized expression, and
`harmonic=False` selects it for `QEred` and the corresponding forcing.
`gauge_gradient` is the symmetrized derivative with one half;
`PC1` is the whole correction \(2P^D{}_{AB}C_{1,D}\).

## Unit-principal signs

The parent's second equations have
\(\Box k=\text{linear couplings}+2Q^{E,\rm red}\) and
\(\Box w=\text{linear couplings}-Q^\phi\). The principal part of
\(\rho^2\Box\) is \(-\partial_t^2+\partial_z^2\).
For \(x=t/\ell,r=z/\ell\), multiplication by
\(-\ell^2\rho^2\) gives the registered principal operator
\(\partial_x^2-\partial_r^2\). Thus the output forcing is exactly

\[
F_g=-2\ell^2\rho^2 Q^{E,\rm red},\qquad
F_\phi=+\ell^2\rho^2 Q^\phi.
\]

The source itself does not include the second unknown k,w. No scalar-only
second metric ansatz, trace-only truncation, or TT deletion occurs.

## Verification and limits

The sibling projection checker expands all 42 independent ordinary metric
jet symbols and all six independent scalar jets. It proves exact polynomial
reconstruction in the sphere quotient ring for every output sector and
exports the registered 23 coefficients. Coefficient comparison to the
parent is universal for QC, scalar Qphi, and Einstein matter terms. This
universal symbolic proof, rather than rational sample checks, establishes
the coefficient identities for arbitrary admissible input values.

The geometry checker separately uses an independent six-sector scalar
Hessian formula and the five-dimensional conformal transformation law.
For \(q=f g_0\), it checks all first and second covariant jets and

\[
Q^{\rm Ric}=3f\nabla^2f+\tfrac92df\otimes df
+g_0(f\Box f-\tfrac12|df|^2).
\]

These are exact all-sphere polynomial comparisons at populated generic
rational y,t jets. Independent review also compares the entire frame
adapter to direct inverse-metric/Christoffel epsilon-series calculations
in curved, nonnormal coordinate families. General off-harmonic inputs
exercise both lowering conventions and the surviving reduction term.
Normal and optimized Python runs must agree; named mutations fail for
specific mathematical diagnostics rather than disabled assertions.

The arithmetic core supports rational, symbolic, and compatible interval
coefficients. Integer numerical inputs remain exact rationals; raw float
inputs are rejected explicitly. Direct polynomial interval evaluation can be wide because it
loses correlations before angular cancellation. For actual numerical
enclosures, evaluate the exported exact coefficient expressions on
authenticated first-field jet intervals. This geometric algebra does not
itself certify that numerical input intervals enclose the actual master
or gauge solution, solve a second evolution, prove a nonlinear remainder,
or establish an invariant force.

The frame adapter specifically implements the displayed warped background
and fixed first quadrupole. Arbitrary first angular profiles, additional
coexact first modes, or a general nondiagonal reference-frame metric are
outside that adapter's API. The lower-level indexed evaluator can be reused
for other supplied covariant jets in a diagonal inverse frame; those jets
and their consistency are then the caller's responsibility.
