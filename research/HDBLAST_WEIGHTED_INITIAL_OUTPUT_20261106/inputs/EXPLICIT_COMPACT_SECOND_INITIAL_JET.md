# An explicit compact second initial jet for the actual aligned packet

15 October 2026. Separate follow-on to the frozen interior-master/source
checkpoint. Conditional on its retained smooth regular background and first
reconstruction. This calculation solves the five second differentiated
**initial constraints** for a specified choice of physical Cauchy variables.
It does not solve second-order evolution, complete generalized harmonic gauge,
or construct an exact nonlinear family.

The previous analytic initial-data construction chose the second scalar and
scalar momentum-density data to vanish and corrected metric and metric
momentum through a weighted inverse. Those choices are optional when selecting
a second jet. Allowing second scalar and unit-normal scalar-velocity data gives
a simpler local correction for this packet. It remains in the original packet
annulus and requires no numerical radial boundary-value solve.

## 1. Actual data and physical Cauchy conventions

At \(t=0\), use the actual time-symmetric spatial background

\[
\gamma_0=dy^2+\rho(y)^2\kappa_{S^3},\quad
K_0=0,\quad\phi_0=\phi_0(y),\quad\Pi_0=0,
\quad a=\rho_y/\rho,\quad p=\phi_{0,y}\ne0.
\]

Our convention is \(K=\tfrac12\mathcal L_n\gamma\) for the future
unit normal and \(\Pi=n\phi\). The physical constraints are

\[
\begin{aligned}
\mathcal H&=R[\gamma]+(\operatorname{tr}_\gamma K)^2-|K|_\gamma^2
             -\Pi^2-|d\phi|_\gamma^2-2U(\phi)=0,\\
\mathcal M_I&=D^J K_{JI}-D_I\operatorname{tr}_\gamma K-\Pi\partial_I\phi=0.
\end{aligned}                                                     \tag{1}
\]

Thus \(\mathcal H\) is twice the normal Einstein constraint, not the
trace-reversed Einstein \(E_{nn}\) alone. \(K\) here is a covariant
tensor and \(\Pi\) an ordinary scalar velocity. They are **not** canonical
momentum densities. Converting to density variables is a nonlinear change of
Cauchy variables and changes the quadratic source by its induced second jet.

The registered model has \(\delta=1/1000\),
\(c=96714321425692/161839258930731\), \(\eta=\phi-1\), and

\[
U(1+\eta)=-\frac2{27}+\frac{14}9\eta^2+\frac{50}{27}\eta^3
          -\frac16\eta^4-\frac49\eta^5-\frac2{27}\eta^6.
\]

The common shell tension is the same retained function and parameter choice
(in particular \(\sigma''=4\phi\), \(\sigma'''=4\)). It is not altered
by this initial correction; the data are background on an open shell collar.

The actual regular branch, its shell values and coefficient enclosures are
the retained background input, authenticated in
`inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json`. In particular
\(129\le\rho\le131\), \(-7/50000\le p\le-3/25000\), and
the background scalar equation is

\[
p_y+4ap=U'.                                                     \tag{2}
\]

The initial packet is specified completely in the inherited radial coordinate
\(dy=-\rho\,dz\) by \(\ell=10^{-6}\),

\[
q_\ell(z)=\begin{cases}
\exp\!\left(4-\dfrac{\ell^2}{(z-\ell)(2\ell-z)}\right),&\ell<z<2\ell,\\
0,&\text{otherwise},
\end{cases}\qquad Y=n_1^2-n_2^2.
\]

With \(\alpha=(\rho/\rho_b)^{3/2}p/p_b\), set

\[
C=\alpha/\rho^2,\quad B=3\alpha/(\rho^3p),\quad
k=\alpha_z/\alpha=\rho(5a/2-U'/p).
\]

At \(t=0\) all amplitude functions used below are fixed by

\[
\begin{aligned}
u&=Cq_\ell,&u_t&=Ce^{-z}q_\ell,\\
v&=B(q_{\ell,z}+kq_\ell),&
v_t&=Be^{-z}\big(q_{\ell,z}+(k-1)q_\ell\big),
\end{aligned}                                                     \tag{3}
\]

and radial derivatives mean \(\partial_y=-\rho^{-1}\partial_z\).
These are the actual first solution's initial data, so no approximation of
its evolved boundary trace is needed. The first physical Cauchy tangent is

\[
\begin{array}{lll}
(\gamma_1)_{yy}=-4uY,&(\gamma_1)_{yi}=0,
 &(\gamma_1)_{ij}=2uY(\gamma_0)_{ij},\\
(K_1)_{yy}=-2u_tY/\rho,&(K_1)_{yi}=0,
 &(K_1)_{ij}=u_tY(\gamma_0)_{ij}/\rho,\\
\phi_1=vY,&&\Pi_1=v_tY/\rho.
\end{array}                                                       \tag{4}
\]

It solves all linearized constraints by the retained coupled reconstruction.
For example, its radial momentum equation is precisely
\(pv_t=-3(u_{yt}+2au_t)\). The letter \(K_1\) denotes the spatial
extrinsic-curvature tangent, not the shell Israel tensor.

We expand physical Cauchy variables as
\(X_\epsilon=X_0+\epsilon X_1+\epsilon^2X_2/2\), where
\(X=(\gamma,K,\phi,\Pi)\). Let
\(Q_H=D^2\mathcal H[X_1,X_1]\),
\(Q_M=D^2\mathcal M[X_1,X_1]\), computed along the affine path in
these **physical Cauchy variables**. The second constraint equations are

\[
D\mathcal H[X_2]=-Q_H,\qquad D\mathcal M[X_2]=-Q_M.               \tag{5}
\]

## 2. Complete source formulas and a nonzero actual witness

Put \(G_Y=|\nabla_\kappa Y|^2\). Direct expansion gives

\[
Q_H=A\,Y^2+D\,G_Y,
\]

where

\[
\begin{aligned}
A={}&-24u_y^2-288a u u_y-24u u_{yy}
 -(192a_y+384a^2+32p^2)u^2\\
&-\rho^{-2}(336u^2+12u_t^2+2v_t^2)
 -2v_y^2-16p u v_y-2U''v^2,\\
D={}&\rho^{-2}(36u^2-2v^2).                                    \tag{6}
\end{aligned}
\]

No \(u_{tt}\) or \(v_{tt}\) remains in this initial constraint.
One independent derivation uses the radial lapse
\(L^2=1-4\epsilon uY\) and angular metric
\(h=\rho^2(1+2\epsilon uY)\kappa\). With
\(j=(a+\partial_y\log(1+2\epsilon uY)/2)/L\),

\[
R[\gamma_\epsilon]=R[h]-12j^2-6L^{-1}j_y-2L^{-1}\Delta_hL.
\]

Its Hessian is

\[
\begin{aligned}
R^{[2]}={}&\left[-24u_y^2-288a u u_y-24u u_{yy}
 -(192a_y+384a^2)u^2-336u^2/\rho^2\right]Y^2
 +36u^2G_Y/\rho^2.
\end{aligned}
\]

The other terms in (6) follow by expanding the scalar gradient, potential,
and the squares of (4). The checker independently evaluates the full
four-dimensional coordinate Ricci expansion, with the unit sphere's normal
coordinate metric two-jet, to verify the complete formula (6).

The momentum source has a radial scalar and an angular gradient only:

\[
\boxed{(Q_M)_y=G\,Y^2,\quad (Q_M)_i=F\,\partial_i(Y^2),}
\]
\[
\boxed{F=\frac{18u u_t-vv_t}{\rho},\qquad
G=\frac{12u u_{yt}-48a u u_t-6u_yu_t-2v_tv_y}{\rho}.}            \tag{7}
\]

In particular no independent coexact angular momentum source needs solving
for this aligned state. To check the factorials, write
\(P^J{}_I=K^J{}_I-\delta^J{}_I\operatorname{tr}K\) and
\(s=u_t/\rho\). Then
\(P_1=\operatorname{diag}(-3s,0,0,0)Y\) and
\(D^2P=\operatorname{diag}(12us,24us,24us,24us)Y^2\).
Differentiating the connection as well as \(P\) yields (7).

The actual packet cannot use zero second initial data. At

\[
t=0,\quad z_*=3\ell/2,\quad n=(1,0,0,0),
\]

one has \(q=1\), \(q_z=0\), \(q_{zz}=-32/\ell^2\),
\(Y=1\), and \(G_Y=0\). In particular

\[
v_y=\frac{32B}{\rho\ell^2}-\frac{B_zk+Bk_z}{\rho}.
\]

The authenticated coefficient boxes imply

\[
\frac{19}{2000}<|B|<\frac3{250},\quad |B_z|<\frac3{10},
\quad |k|<270,\quad |k_z|<24000,
\]

so rational triangle inequalities give

\[
2.3\times10^9<|v_y(z_*)|<3\times10^9.
\]

The remaining terms of \(A\) have total absolute value less than
\(5\times10^6\): coarse sufficient bounds are
\(u<10^{-4}\), \(|u_y|<1\), \(|u_{yy}|<10^6\),
\(|u_t|<1\), \(|v|,|v_t|<4\),
\(|a|,|a_y|,|p|<1\), \(|U''|<4\). Each follows directly from
(3) and the authenticated background boxes at this point. Consequently

\[
\boxed{Q_H(z_*,n)<-2(2.3\times10^9)^2+5\times10^6<-10^{19}.}    \tag{8}
\]

If all four second physical Cauchy data were zero, the left side of the
second Hamiltonian equation (5) would be zero. Equation (8) excludes that
choice. Setting all second spacetime initial metric/scalar values and
coordinate velocities to zero also fails this Hamiltonian constraint; the
Hamiltonian Hessian is unchanged by the associated second lapse conversions
at \(K_0=\Pi_0=0\).

## 3. A local Hamiltonian right inverse in the three actual sectors

Let
\(Z_2=n_1^2+n_2^2-1/2\) and
\(H_4=Y^2-(n_1^2+n_2^2)/2+1/12\). The angular identities are

\[
Y^2=\frac16+\frac12Z_2+H_4,\qquad
G_Y=\frac43+2Z_2-4H_4.
\]

Thus \(Q_H=S_0+S_2Z_2+S_4H_4\), where

\[
S_0=A/6+4D/3,\qquad S_2=A/2+2D,\qquad S_4=A-4D.               \tag{9}
\]

For \(W_l=(1,Z_2,H_4)\), \(\lambda_l=(0,8,24)\), choose the
second spatial metric and scalar amplitudes

\[
\begin{aligned}
(\gamma_2)_{yy}&=(\gamma_2)_{yi}=0,\qquad
(\gamma_2)_{ij}=\sum_{l=0,2,4}f_l(y)W_l(\gamma_0)_{ij},\\
\phi_2&=\sum_{l=0,2,4}\zeta_l(y)W_l,\qquad
\boxed{f_l=-\frac{\rho^2S_l}{2(\lambda_l-3)},\quad
\zeta_l=-\frac{3f_{l,y}}{2p}.}                                 \tag{10}
\end{aligned}
\]

To verify it, the linear curvature variation for angular metric amplitude
\(f_lW_l\) is

\[
\delta R=-3f_{l,yy}-12af_{l,y}
                  +2(\lambda_l-3)f_l/\rho^2.
\]

Its scalar correction contributes

\[
-2p\zeta_{l,y}-2U'\zeta_l
 =3f_{l,yy}+3(U'/p-p_y/p)f_{l,y}
 =3f_{l,yy}+12af_{l,y},
\]

using (2). The radial derivatives cancel exactly and

\[
D\mathcal H[(\gamma_2,0,\phi_2,0)]
 =\sum_l\frac{2(\lambda_l-3)}{\rho^2}f_lW_l=-Q_H.              \tag{11}
\]

The three denominators are **−6, 10, 42**. The excluded eigenvalue three
would be the dipole; none occurs in the actual quadratic Hamiltonian source.
The monopole is included explicitly, rather than extrapolated from a
nondegenerate angular trace-free equation.

## 4. A local momentum right inverse and simultaneous completion

Choose the remaining second physical Cauchy data as

\[
\boxed{K_2=q\gamma_0,\quad q=\frac{F}{3}Y^2,\qquad
\Pi_2=\frac{G-F_y}{p}Y^2.}                                    \tag{12}
\]

Since there are four spatial dimensions, \(K_2-\gamma_0\operatorname{tr}K_2
=-3q\gamma_0\), hence

\[
D\mathcal M[(0,K_2,0,\Pi_2)]
 =-3dq-p\Pi_2dy
 =-F\,d_{S^3}(Y^2)-G Y^2dy=-Q_M.                              \tag{13}
\]

Time symmetry makes the two constructions independent:
\(D\mathcal H\) has no \(K_2\) or \(\Pi_2\) term at
\(K_0=\Pi_0=0\); \(D\mathcal M\) has no \(\gamma_2\) or
\(\phi_2\) term there. Thus (10) and (12) together solve all five
second differentiated constraints (5). They introduce no second angular
transverse-traceless harmonic data. They select a different particular second
jet from the previous choice prescribing zero second scalar/density data.

All profiles are explicit functions of the exact initial packet (3) and
actual background. In particular there is no undetermined integration
constant, boundary shooting condition, homogeneous radial tail, or adjoint
solvability condition for these three scalar sectors.

## 5. Coordinate conversion, support and remaining evolution boundary work

There are two different affine paths. For an affine **spacetime** metric
\(g_0+\epsilon h\) with zero second coordinate velocities, the lapse is
\(N_\epsilon=\rho\sqrt{1+2\epsilon uY}\). It induces

\[
K_2^{\rm induced}=-2uY K_1,\qquad
\Pi_2^{\rm induced}=-2uY\Pi_1.
\]

Their linear momentum contribution, using the first momentum equation,
is \(6u_yu_tY^2/\rho\). Therefore the spacetime-affine radial momentum
Hessian has the same \(G\) as (7) with the term \(-6u_yu_t/\rho\)
removed. Mixing that source with the physical-Cauchy-affine completion
(12) would leave a nonzero constraint. The angular source \(F\) and
Hamiltonian source (6) agree for the two paths.

To convert our completed physical data into coordinate velocities, take the
retained first lapse \(N_1=\rho uY\), first shift zero, and choose the
second lapse and shift values zero. Then

\[
\boxed{\dot\gamma_2=2\rho(K_2+2uY K_1),\qquad
\dot\phi_2=\rho(\Pi_2+2uY\Pi_1).}                             \tag{14}
\]

Zero second lapse value means \((g_{tt})_2=-2\rho^2u^2Y^2\),
since \(g_{tt}=-N^2\) at zero shift. These are specified initial
coordinate data, not a relative-harmonic gauge completion. The five initial
gauge-velocity equations and the first gauge-vector solve still have to be
performed in a common representative, with the quadratic transformation
terms retained as in the earlier full-contact note.

The bump is smooth and flat at both \(z=\ell\) and \(z=2\ell\).
Every operation in (3), (6), (7), (9), (10), and (12) is finite
differentiation or multiplication by a smooth background factor, with
\(p\) bounded away from zero. All second physical profiles therefore have
support in the **same** closed annulus \([\ell,2\ell]\), are flat at
its endpoints, and extend by zero smoothly. They vanish on open
neighborhoods of the initial physical shell and the vertex. They are
reflection invariant and contain only scalar degrees zero, two, four.
Both identical copies can use the same construction. No sum-of-normal
convention changes an initial bulk constraint.

The earlier enlarged-annulus weighted inverse remains a valid way to impose
its different scalar prescriptions; the present local choice does not compute
that particular inverse. The new correction has the same initial support as
the packet, so its background causal enlargement over \(T=3\ell\) lies
within \(z\le5\ell\), inside the existing \(8\ell\) collar.
Perturbed nonlinear characteristic speeds and any nonlinear family's
corrections are outside this linear second-coefficient support statement.

Because all first and second physical data vanish on an open initial shell
collar, their physical shell corner conditions are the background ones.
After a **local** compatible gauge completion, its initial jets also vanish
there, and the reduced equation's corner recursion is homogeneous. This
does not construct the subsequent transformed contact forcing or certify
the entire second wave evolution.

The next quantitative tasks are to enclose the values/norms of (10) and
(12), complete the first gauge transformation and transformed second data,
and solve and validate the finite coupled initial-boundary system with all
physical/gauge/subsidiary rows. The known C3 first-master bounds are sufficient
in derivative count for the *values* of this second initial jet: \(f_l\)
uses first master derivatives through two, while \(f_{l,y}\) and \(F_y\)
use them through three, together with the needed background coefficient jets.
No such second-initial norm enclosure is supplied here.

## 6. Reproduction and limits of the exact checks

`check_initial_constraint_completion.py` uses only the Python standard library.
It authenticates the actual background report by its fixed SHA256. At three
deterministic nonspecial rational jets it independently expands inverse metric,
Christoffel symbols, Ricci tensor and scalar curvature in four dimensions,
including the unit sphere's normal-coordinate curvature; it verifies (6).
A separate direct coordinate expansion of the momentum divergence verifies
both rows of (7), including connection variations. Exact rational polynomial
checks verify the three projections, the three Hamiltonian inverse rows,
momentum completion and affine-path conversion. The packet-maximum consequence
(8) is an exact inequality using actual authenticated coefficient boxes.

```sh
python check_initial_constraint_completion.py --run-controls --output INITIAL_CONSTRAINT_COMPLETION.json
python -O check_initial_constraint_completion.py --run-controls --output INITIAL_CONSTRAINT_COMPLETION_OPTIMIZED.json
```

The seven wrong-formula controls omit the radial momentum cross term, halve
the angular momentum source, flip the scalar compensator, change a sector
denominator, halve the Hamiltonian Hessian, omit its angular-gradient term,
or remove the packet's second derivative at the witness. Each must fail
for its mutation-specific reason and leave no output report. Normal and
optimized runs retain the same mathematical results.

Generic point-jet tests supplement the analytic derivation; they do not prove
the full field identities at all points. The explicit formulas and identities
(2), (11), and (13) give the all-point constraint proof. The certificate (8)
is an actual pointwise obstruction, rather than a generic-jet example.
This follow-on supplies no external peer review, second evolution, numerical
U2 amplitude, invariant force, exact nonlinear family or remainder estimate.
