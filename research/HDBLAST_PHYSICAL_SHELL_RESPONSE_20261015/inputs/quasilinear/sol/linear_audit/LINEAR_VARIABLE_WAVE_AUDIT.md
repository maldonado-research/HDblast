# Independent variable-coefficient linear-wave audit

The angular-only construction below proves the smooth linear solver and the high-order estimate without commuting an angular projection through a coefficient. It suffices directly for data and forcing whose initial shell corner jets vanish on an open radial collar. A jet lift extends the construction to arbitrary smooth compatible data.

## Precise hypotheses

Write the scalar-principal bundle system, after multiplication by a positive scalar and placing spatial principal terms in divergence form, as

\[
 Lu=a u_{tt}+2b^iX_i u_t-X_i(c^{ij}X_ju)
       +B^t u_t+B^iX_i u+Cu=f. \tag{1}
\]

Here \(X_0=\partial_z\), while \(X_1,X_2,X_3\) are a fixed left-invariant round-\(S^3\) frame. Use a fixed split-preserving global bundle frame; changing from a geometric connection adds only the displayed lower-order terms. The three angular fields have zero round divergence, constant commutators, and commute with the scalar round Laplacian. All principal coefficients in (1) are real scalars, and \(c^{ij}=c^{ji}\). On the whole time slab suppose

\[
a\ge a_*>0,\qquad c^{ij}\xi_i\xi_j\ge c_*|\xi|^2,
\qquad b^0|_{z=0}=0,\quad c^{0A}|_{z=0}=0.\tag{2}
\]

The last conditions are precisely the time-normal and angular-normal parts of \(g^{za}=0\). At the shell impose the fixed orthogonal split

\[
u_D=0,\qquad \partial_z u_N=0.\tag{3}
\]

At the artificial far wall one can impose full Dirichlet conditions, with data supported behind a causal margin. Then no additional normal-cross-coefficient condition is needed there. Alternatively use a second D/N wall satisfying (2). The proof requires ordinary fixed positive fiber norm, and no symmetry of the lower-order matrix coefficients.

## Basic energy

For

\[
E[u]=\tfrac12\int(a|u_t|^2+c^{ij}\langle X_i u,X_j u\rangle+|u|^2),
\]

integration by parts gives

\[
E'[u]\le C E[u]+C\sqrt{E[u]}\|f\|_2.\tag{4}
\]

The time-space mixed term has boundary flux \(b^0|u_t|^2\). The spatial term has flux \(\langle u_t,c^{0j}X_j u\rangle\). Both vanish at the shell: (2), \(u_{D,t}=0\), \(u_{N,z}=0\), and orthogonality of the fixed split give this directly. At a full Dirichlet far wall \(u_t=0\). Only bounded first derivatives of the principal coefficients and bounded lower coefficients enter (4).

## Uniform high-order estimate

Let \(\mathscr T=\{\partial_t,X_1,X_2,X_3\}\) and

\[
E_m^{\rm tan}=\sum_{|I|\le m-1}E[\mathscr T^Iu],\qquad
H_m[u](t)=\sum_{j=0}^m\|\partial_t^ju(t)\|_{H^{m-j}}.
\]

Every \(\mathscr T^Iu\) satisfies (3). Its equation has forcing
\(\mathscr T^If+[L,\mathscr T^I]u\). This commutator has total order at most \(m\), including the terms in which a derivative strikes a principal coefficient. It is not legitimate to omit those terms.

Recover normal derivatives from (1), rather than commuting the boundary energy with \(\partial_z\). Solve (1) for \(c^{00}\partial_z^2u\). Apply \(\partial_z^{r-2}\mathscr T^I\), where \(r+|I|\le m\), and induct first on total order and then on the number \(r\) of normal derivatives. Every other top-order principal term has fewer normal derivatives; differentiated coefficients multiply strictly lower total derivatives. Consequently

\[
H_m[u](t)\le C_m\left(\sqrt{E_m^{\rm tan}(t)}+
 K_m[f](t)\right),\qquad
K_m[f]=\sum_{j=0}^{m-2}\|\partial_t^jf\|_{H^{m-2-j}}.\tag{5}
\]

The sum is empty when \(m=1\). This includes all time derivatives, not just spatial ones. In particular a commutator term containing \(\partial_z^2\mathscr T^{m-2}u\) has exactly the regularity supplied by (5).

Insert (5) in the sum of commuted versions of (4), then use Gronwall. The pointwise forcing terms in (5) obey

\[
\sup_t\|\partial_t^jf(t)\|_{H^{m-2-j}}
\le \|\partial_t^jf(0)\|_{H^{m-2-j}}
 +\|\partial_t^{j+1}f\|_{L^1_tH^{m-2-j}}.
\]

The initial tangential energies are controlled by the recursively computed time jets. Thus the following sufficient, non-sharp estimate holds:

\[
\sup_{0\le t\le T}H_m[u](t)\le C_m\Big[
\|u_0\|_{H^m}+\|u_1\|_{H^{m-1}}
+\sum_{j=0}^{m-1}\|\partial_t^jf\|_{L^1_tH^{m-1-j}}
+\sum_{j=0}^{m-2}\|\partial_t^jf(0)\|_{H^{m-2-j}}
\Big].\tag{6}
\]

The constants depend on positivity, the time slab, and finitely many coefficient norms, but not an angular cutoff. For nonlinear iteration, a tame version must retain the Moser product structure of the commutators and normal recovery. A bare unspecified \(C_m\) depending on all high coefficient norms does not by itself establish the iteration estimates.

## Genuine construction by angular-only Galerkin

Let \(P_n\) be the componentwise scalar round-harmonic orthogonal projection in the fixed global left-invariant split frame, applied to the same finite harmonic set in every component. It commutes with \(\partial_t,\partial_z,X_A\), preserves the D/N split and radial support, and is an \(L^2\) contraction. This specially chosen frame and these globally defined angular generators are essential: arbitrary coordinate-component cutoffs need not have these commutation properties. Unlike the prior background-equivariant construction, no commutation of the full operator with \(P_n\) is claimed. In particular \(P_n\) generally does **not** commute with any variable coefficient. Solve

\[
P_n Lu_n=P_nf,\qquad u_n=P_nu_n,\tag{7}
\]

with the homogeneous boundary conditions and projected initial data.

For fixed \(n\), angular derivatives become finite matrices. Introduce \(v=u_{n,t}\), \(w=u_{n,z}\). The radial principal part is symmetric hyperbolic, with

\[
A_n=P_n aP_n>0,\quad G_n=P_nc^{00}P_n>0,\quad
\mathcal A_0=\begin{pmatrix}A_n&0\\0&G_n\end{pmatrix},\quad
\mathcal A_1=\begin{pmatrix}2P_nb^0P_n&-G_n\\-G_n&0\end{pmatrix}.
\]

The remaining angular terms are bounded zero-order matrices in \((u_n,v,w)\) for that fixed \(n\). At the shell the boundary space \(v_D=0,w_N=0\) is a maximal zero-flux subspace; \(G_n\) preserves the split and \(b^0=0\) there. The following radial variational construction spells out existence without assuming any unproved normal regularity.

First perform the jet lift described below, so initial data vanish and the residual forcing is flat in time at zero. At fixed \(n\), the projected radial equation has the form

\[
A_nu_{tt}+2H_nu_{tz}-\partial_z(G_nu_z)
+L_tu_t+L_zu_z+L_0u=f_n.\tag{8}
\]

All matrices are smooth, \(A_n,G_n\) are symmetric positive, \(H_n\) is symmetric, and \(H_n(0)=0\). Angular-normal divergence terms have been expanded into the lower-order matrices; their boundary terms vanish because \(c^{0A}|_b=0\). Use any fixed dense smooth Galerkin basis in

\[
V=\{q\in H^1([0,Z]):q_D(0)=0,\ q(Z)=0\}.
\]

Use the variational spatial form \(\int\langle G_nu_z,q_z\rangle\); Neumann is a natural condition. Coercivity of the finite mass matrix gives each finite radial ODE solution. All of its initial time jets vanish, since the residual source is flat.

Here is the essential uniform-in-radial-dimension regularity argument. Put \(v_p=\partial_t^pu\), with its usual energy \(E_p\). Differentiating the Galerkin equation \(p\) times is legitimate. The only apparently non-\(L^2\) commutator is the differentiated spatial form, whose energy contribution is

\[
-\sum_{k=1}^p\binom pk\int\left\langle
  (\partial_t^kG_n)\partial_z\partial_t^{p-k}u,
  \partial_z\partial_tv_p\right\rangle.
\]

Add the correction

\[
C_p=\sum_{k=1}^p\binom pk\int\left\langle
 (\partial_t^kG_n)\partial_z\partial_t^{p-k}u,
 \partial_zv_p\right\rangle
\]

to \(E_p\). Its derivative cancels the displayed term. Every surviving term involves only energies \(E_q\), \(q\le p\), because \(p-k+1\le p\). The time-mass commutator contains at worst \(\partial_t^{p+1}u\), and the mixed-principal commutator contains at worst \(\partial_z\partial_t^pu\), also controlled by \(E_p\). Moreover

\[
|C_p|\le\varepsilon E_p+C_{p,\varepsilon}\sum_{q<p}E_q.
\]

Adding recursively sufficiently large multiples of the lower energies gives a positive energy equivalent to \(\sum_{q\le p}E_q\) and a Gronwall estimate with \(\sum_{q\le p}\|\partial_t^qf_n\|_{L^1_tL^2_z}\). Constants may depend on \(n,p\), but not on radial Galerkin dimension. No \(u_{zz}\) estimate was used.

Take the radial weak limit. All time derivatives now belong to \(L^\infty_tH^1_z\), with the next time derivative in \(L^\infty_tL^2_z\). Equation (8) therefore gives \(u_{zz}\in L^\infty_tL^2_z\) by inversion of \(G_n\). Differentiating in time and inducting gives the same statement for every time derivative. Smooth one-sided coefficients and the equation then recover arbitrarily many normal derivatives. The weak natural condition is \((G_nu_z)_N=0\), which is equivalent to \(u_{N,z}=0\); the Dirichlet traces are preserved by \(V\). Sobolev embedding gives a smooth fixed-\(n\) solution up to the corners. This completes the finite-angular solver used in (7).

The full geometric energy and high-order estimates (4)–(6), applied to the constructed solutions, hold uniformly in \(n\). When testing (7), orthogonality removes \(P_n\) from the pairing. When commuting, retain
\(P_n[L,\mathscr T^I]u_n\); no coefficient commutation is used. For normal recovery replace division by \(c^{00}\) with inversion of

\[
G_n(t,z):\operatorname{ran}P_n\to\operatorname{ran}P_n,
\qquad \langle G_nv,v\rangle\ge c_*\|v\|_2^2.
\]

Hence \(\|G_n^{-1}\|_{2\to2}\le c_*^{-1}\), uniformly. Normal derivatives commute with \(P_n\); differentiating the equation and using the same induction proves the projected analogue of (5). No uniform Sobolev norm of \(G_n^{-1}\) is assumed.

If the initial data and all initial forcing jets vanish on an open shell collar, their angular projections retain that property. The finite projected systems then satisfy all corner identities. Uniform bounds at every order, compactness on the compact spatial cylinder, and a diagonal subsequence give a smooth solution. To pass the equation, first test against finite angular sums and then use their density. The same bounds preserve the traces. Equation (4) gives uniqueness, so subsequence choices do not matter. Finite propagation follows by the same localized energy argument for the limiting local equation; it need not hold for the angularly projected approximants.

## General compatible corners and specific traps

For arbitrary smooth compatible data, let \(J_j\) be their exact time jets from (1). Each jet obeys \((J_j)_D=0\), \(\partial_z(J_j)_N=0\). Choose a smooth Borel lift \(W(t)\) with these jets, using scalar time cutoffs applied to the \(J_j\). The lift satisfies the boundary conditions at every time. Then \(r=u-W\) has zero initial data and \(f-LW\) is flat in time at zero. Apply (7) to this remainder. Every projected problem has zero corner jets. For smooth input and a fixed finite regularity target, a sufficiently long finite Taylor lift can replace the Borel lift. This is not a construction of a bounded lift from arbitrary finite-Sobolev jets: \(J_j\in H^{m-j}\) alone does not make \(t^jJ_j\) an \(H^m\) lift. Estimate (6) follows directly for the resulting exact smooth solution, so no quantitative estimate for the Borel lift is required in (6).

The following shortcuts are invalid without further work:

* Projecting general compatible initial data and declaring the projected variable-coefficient corner recursion compatible. Projection and coefficient multiplication do not commute.
* Using full spatial mixed-Laplacian eigenprojections as uniformly bounded high-Sobolev approximations of arbitrary smooth mixed-boundary data. Powers of the auxiliary Laplacian impose additional boundary identities. Such projections also destroy radial corner support.
* Commuting normal derivatives through the energy and claiming homogeneous boundary conditions for them. Even \(\partial_zu_D\) has no homogeneous Dirichlet condition, and \(\partial_z^2u_N\) has no homogeneous Neumann condition.
* Treating \(P_nc^{00}P_n\) as multiplication by \(c^{00}\). Only its uniform coercive inverse is needed.
* Dropping \(g^{zt}|_b=0\): its boundary flux is separate from the spatial conormal flux.
* Claiming arbitrary \(L^2\) Neumann boundary control from the homogeneous energy. The boundary trace of \(u_t\) has not been estimated.

These observations establish the frozen linear step. They do not by themselves establish tame quasilinear bounds, nonlinear compatibility propagation, or a lifespan through a prescribed contact time.
