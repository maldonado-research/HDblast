# Second-order spacetime evolution before shell contact

14 October 2026 (America/Los_Angeles). A bounded continuation of the localized
second-order initial-data checkpoint. Substantial AI assistance (OpenAI
Codex and collaborating agents); internal symbolic checks and automated
reviews are not external peer review. No cosmological conclusion follows.

## 1. What is added

The localized checkpoint constructs compactly supported second-order
**constraint data**. This note extends those data to second-order
**spacetime coefficients on a nonzero time slab** during which the
perturbation remains separated from the shell. Every shell junction then
holds through second order throughout that slab, with zero shell motion.
The result does not solve the quadratic junction problem at shell contact.

For the velocity-only quadrupole, we also give the exact initial spacetime
accelerations in Gaussian proper time. A separate angular calculation
identifies a tensor-harmonic component of the quadratic scalar stress that
a scalar-only continuation would need to account for. Neither result is
a claim of an exact nonlinear solution or of scientific novelty.

## 2. A finite pre-contact existence statement

Let the actual background be

\[
g_0=\Omega(s)^2\eta,\quad\phi_0=\Phi(s),\quad
s=R^2-T^2,\quad W_0:\ R^2-T^2=b^2,
\]

with the same bulk action, potential, and two-copy shell model as the
shell-Noether and localized checkpoints. Assume it is smooth and
nondegenerate on the closed region used below; \(\Omega>0\). Smoothness
through the cone is part of the actual inherited background assumption,
not an arbitrary change of the potential. If only the positive-\(s\)
background interval is being used, impose the additional bound given below.

Write the initial canonical data as

\[
x(\varepsilon)=x_0+\varepsilon v+\tfrac12\varepsilon^2w,
\qquad DC[v]=0,\qquad DC[w]+D^2C[v,v]=0.
\tag{1}
\]

Both \(v\) and \(w\) are smooth and supported in the **same containing
annulus** \(R_-\le R\le R_+\), with \(0<R_-<R_+<b\), in each copy. It is
the union of their supports that matters: the monopole correction's chosen
bump may be larger than the first-order packet. Initial shell displacement,
tilt, and their first time derivatives vanish at both orders.

**Proposition 1 (formal second-order evolution).** On every compact smooth
background slab with

\[
0<T_*<T_c:=\frac{b^2-R_+^2}{2R_+},
\tag{2}
\]

there are smooth bulk coefficients \((h_1,\psi_1)\) and
\((h_2,\psi_2)\), unique after a fixed hyperbolic gauge and its Cauchy gauge
data are specified, that induce \(v,w\), solve all Einstein–scalar field
equations through order \(\varepsilon^2\), and vanish on an open shell
collar for \(|T|\le T_*\). All matching, Israel, scalar, and normal-force
junctions therefore hold through this order with zero first- and
second-order shell embeddings. Choose the initial lapse, shift, and their
parameter and time jets to agree with the background outside the containing
annulus. If the available background control is restricted to \(s>0\),
take also \(T_*<R_-/2\). More generally the whole possible support interval
\([\max(0,R_--T_*)^2-T_*^2,R_+^2+2R_+T_*]\) must lie in the inherited controlled
background interval. Positivity of its lower endpoint alone is insufficient.
The lower endpoint is \(R_-^2-2R_-T_*\) while \(T_*\le R_-\), and
\(-T_*^2\) after the support cone can reach the spatial origin.

*Proof.* Use a background-adapted generalized harmonic (wave-map) reduction.
Its gauge vector is

\[
C^A(g)=g^{BC}\big(\Gamma(g)^A{}_{BC}-\Gamma(g_0)^A{}_{BC}\big).
\]

The reduced Einstein–scalar equations have a diagonal normally hyperbolic
principal part, with background null cones at each perturbative order.
At order one solve the homogeneous linearized reduced system. At order two
solve the same linear operator with the known smooth quadratic source:

\[
Lz_1=0,\qquad Lz_2=-D^2\mathcal E^{\mathrm{red}}_0[z_1,z_1].
\tag{3}
\]

Here \(z_j=(h_j,\psi_j)\) and \(L=D\mathcal E^{\mathrm{red}}_0\); equation
(3) is a definition of the full coupled source, not a claim that its
components reduce to the earlier scalar master. Choose initial lapse and
shift parameter jets in the usual reduced Cauchy construction so that
\(C^{(1)}=C^{(2)}=0\) initially. The differentiated Hamiltonian and momentum
constraints in (1) give the corresponding normal gauge derivatives.
The scalar equation and the differentiated contracted Bianchi identity
then propagate zero gauge through these orders. Consequently the
coefficients solve the unreduced field equations through order two.
This uses the standard Cauchy theorem for a normally hyperbolic operator
on a smooth vector bundle over a globally hyperbolic background, with a
prescribed smooth source and smooth Cauchy data; see Bär–Ginoux–Pfäffle,
*Wave Equations on Lorentzian Manifolds and Quantization*, Theorem 3.2.11
(2007, [arXiv:0806.1036](https://arxiv.org/abs/0806.1036)). The general
Einstein hyperbolic reduction and constraint propagation are described in
Choquet-Bruhat, *General Relativity and the Einstein Equations* (2009),
Chapter VI. At finite order the scalar coupling only contributes smooth
lower-order coefficients and source terms to this vector-bundle wave
system. The argument never assumes that (1) extends to exact constraint data.

The source in (3) is local and quadratic in \(z_1\) and its derivatives.
Finite propagation for the normally hyperbolic system therefore gives,
for \(t=|T|\),

\[
\operatorname{supp}z_j(T,\cdot)\subset
\{\max(0,R_--t)\le R\le R_++t\},\qquad j=1,2.
\tag{4}
\]

One can solve on an auxiliary extension of the bulk chart; the portion
of that extension outside the controlled causal region is never used.
Equivalently one uses the Cauchy domain of dependence. The shell remains
strictly outside the upper support bound since

\[
R_b(T)^2-(R_++t)^2=b^2-R_+^2-2R_+t=2R_+(T_c-t)>0.
\tag{5}
\]

The positive-\(s\) restriction, when needed, follows from
\(s\ge(R_--t)^2-t^2=R_-(R_--2t)>0\). Strict compact-slab inequalities
give an open collar where both coefficients and all their derivatives
vanish. Pullback and geometric junction residuals there agree with the
background through order two. This proves every claimed shell condition.
\(\square\)

**Corollary 1.1 (every finite order in the reflection-invariant sector).**
Take the common-support reflection-invariant formal initial-data sequence
constructed in [Reflection symmetry and formal initial data](REFLECTION_FORMAL_INITIAL_DATA.md).
Choose \(R_-,R_+\) to enclose all its coefficient supports, as permitted
by the nested-annulus construction. On any one fixed strict slab satisfying
Proposition 1's background and causal bounds, for every finite \(N\) there
are smooth spacetime coefficients through order \(N\) inducing the first
\(N\) initial coefficients and satisfying all bulk and shell equations
through that order. The shell coefficients can be zero at every order.

Indeed, after \(z_1,\ldots,z_{n-1}\) are known, the \(n\)th differentiated
reduced equation has the form
\(Lz_n=-F_n^{\mathrm{red}}(z_1,\ldots,z_{n-1})\), with the same normally
hyperbolic operator \(L\). The source is a known smooth local differential
expression supported in the same background causal development. Apply
the Cauchy theorem, choose localized gauge jets, and use the differentiated
Bianchi identity and the \(n\)th initial constraint equation to propagate
zero gauge. Induction gives (4) at each order, and hence the same shell
collar. There is no estimate uniform in \(n\) and no convergence assertion.

This constructs a formal order-two solution (and, in Corollary 1.1's
sector, a formal solution at every finite order): the parameter Taylor
coefficients of the equations through that order vanish. The polynomial
\(g_0+\varepsilon h_1+\varepsilon^2h_2/2\) is generally not an exact
solution. There is no estimate here on a nonlinear remainder, convergence,
higher-order shell contact, or an all-angular shell solution space.

## 3. Exact initial acceleration for a velocity quadrupole

All formulas in this section use **Gaussian proper time** \(\tau\), unit
lapse and zero shift, with \(\tau=0\) the initial surface. It differs from
the inherited conformal ambient time \(T\), whose initial lapse is
\(\Omega\). The causal bounds (2)–(5) use \(T\), not \(\tau\).

Let

\[
K_{ij}=\tfrac12\partial_\tau\gamma_{ij},\quad
\pi^{ij}=\frac{\sqrt\gamma}{2}(K^{ij}-\gamma^{ij}K),\quad
u=\partial_\tau\phi=\frac{\varpi}{\sqrt\gamma}.
\]

These reproduce the checkpoint's Hamiltonian factors. The exact bulk
equations in this chart are

\[
\begin{aligned}
\partial_\tau K_{ij}
&=-R^\gamma_{ij}+2K_i{}^aK_{aj}-KK_{ij}
  +\partial_i\phi\partial_j\phi+\tfrac23U\gamma_{ij},\\
\partial_\tau^2\phi&=\Delta_\gamma\phi-Ku-U'.
\end{aligned}\tag{6}
\]

Choose the admissible velocity-only packet

\[
q=0,\quad p=\mathcal K(s)H,\quad H=x_1^2-x_2^2,
\quad s=|x|^2.
\]

Its ten quadratic tests vanish by the same commuting/reflection argument
as in the localized checkpoint. That checkpoint gives compactly supported
second-order fields \((k_{ij},\chi)\) with vanishing second-order momenta.
At \(\tau=0\) write

\[
\gamma=\gamma_0+\tfrac12\varepsilon^2k,\quad
\phi=\Phi+\tfrac12\varepsilon^2\chi,\quad
\pi=\varepsilon\pi_1,\quad\varpi=\varepsilon\varpi_1.
\tag{7}
\]

Then

\[
\begin{aligned}
K_{1,ij}&=\frac{p}{\Omega}(s\delta_{ij}-3x_ix_j),&
\operatorname{tr}_{\gamma_0}K_1&=\frac{sp}{\Omega^3},\\
u_1&=-\frac{3(2s\mathcal K'+7\mathcal K)H}{2\Omega^3\Phi'},&
\partial_\tau k|_0&=0,\qquad\partial_\tau\chi|_0=0.
\end{aligned}\tag{8}
\]

The prime in (8) is \(d/ds\). Twice differentiating the exact equations
(6) in the perturbation parameter gives the initial accelerations

\[
\begin{aligned}
\partial_\tau^2 k_{ij}|_0={}&
-2(D\operatorname{Ric}_{\gamma_0}[k])_{ij}
+2(\Phi_i\chi_j+\chi_i\Phi_j)\\
&+\tfrac43\big(Uk_{ij}+U'\chi\gamma_{0,ij}\big)
+\boxed{\frac{4sp^2}{\Omega^4}(s\delta_{ij}+9x_ix_j)},\\
\partial_\tau^2\chi|_0={}&
\Delta_{\gamma_0}\chi+(D\Delta_{\gamma_0}[k])\Phi-U''\chi
+\boxed{\frac{3s\mathcal K(2s\mathcal K'+7\mathcal K)H^2}
                  {\Omega^6\Phi'}}.
\end{aligned}\tag{9}
\]

All background potentials in (9) are evaluated at \(\Phi(s)\). The
metric-dependent scalar term is explicitly

\[
(D\Delta[k])\Phi
=-k^{ij}D_iD_j\Phi-
\left(D_i k^{ij}-\tfrac12D^j\operatorname{tr}_{\gamma_0}k\right)D_j\Phi.
\]

The boxes are the closed-form quadratic **kinetic** forcing; the other
terms are the required linear response of the second-order correction.
They cannot be discarded when computing the full acceleration.

For general displacement–velocity data, converting second-order
densitized momenta into time derivatives also generates \(h_1\pi_1\)
and \(h_1\varpi_1\) cross terms. They vanish in (7). Nor should one infer
the time evolution of the inherited reconstruction by simply putting
\(p=\partial_Tq\) into its spatial fields: reconstruction and Gaussian
gauge must be related consistently, including lapse/shift variations.

## 4. An explicit tensor sector of the quadratic scalar stress

Scalar quadrupoles on \(S^3\) form the \((1,1)\) representation of
\(SO(4)\). Its symmetric square is

\[
\operatorname{Sym}^2(1,1)
=(0,0)\oplus(1,1)\oplus(2,2)\oplus(2,0)\oplus(0,2),
\tag{10}
\]

with dimensions \(1+9+25+5+5=45\). The first three are scalar angular
sectors \(l=0,2,4\); the final two are the two transverse, trace-free
tensor sectors at \(l=2\). The scalar identity

\[
H^2=\underbrace{H^2-\tfrac12s(x_1^2+x_2^2)+\tfrac1{12}s^2}_{l=4}
+\tfrac12s\underbrace{(x_1^2+x_2^2-s/2)}_{l=2}+s^2/6
\tag{11}
\]

does not account for every tensor-valued quadratic source.

Here is a concrete exact witness. On the unit sphere let
\(B=\operatorname{diag}(1,-1,0,0)\), \(Y=n^TBn\). Form the algebraic
curvature tensor
\(R_{ikjl}=B_{ij}B_{kl}-B_{il}B_{kj}\) and its Weyl part \(C\) in four
Euclidean dimensions. This fixes the normalization of the curvature
product. Its nonzero sectional entries are

\[
C_{1212}=C_{3434}=-1/3,\quad
C_{1313}=C_{1414}=C_{2323}=C_{2424}=1/6.
\]

Set \(\mathcal T_{ij}(n)=C_{ikjl}n^kn^l\). Weyl symmetries imply
\(\mathcal Tn=0\), \(\operatorname{tr}\mathcal T=0\), and
\(\operatorname{div}_{S^3}\mathcal T=0\). It is a nonzero tensor harmonic
in the two final sectors of (10). In fact, let

\[
F_2=n_1^2+n_2^2-1/2,\qquad
F_4=Y^2-(n_1^2+n_2^2)/2+1/12.
\]

The following is an exact tensor identity everywhere on the unit sphere:

\[
\boxed{(DY\otimes DY)^{\mathrm{TF}}
=\tfrac12(\operatorname{Hess}F_2)^{\mathrm{TF}}
+\tfrac13(\operatorname{Hess}F_4)^{\mathrm{TF}}
-\tfrac43\mathcal T.}\tag{12a}
\]

For a direct calculation, put \(\Pi=\delta-nn^T\), \(A=B^2\), and
\(v=\Pi Bn\). Extend \(F_2,F_4\) to their homogeneous ambient harmonics
of degrees two and four as in (11). The trace-free part of the sphere
Hessian is \((\Pi\operatorname{Hess}_{\mathbb R^4}F_l\Pi)^{\mathrm{TF}}\);
the normal derivative term in the sphere Hessian is pure trace. The Weyl
subtraction and those two elementary Hessians give

\[
\begin{aligned}
\mathcal T&=[Y\Pi B\Pi-v\otimes v+\tfrac12\Pi A\Pi]^{\mathrm{TF}},\\
(\operatorname{Hess}F_2)^{\mathrm{TF}}&=2[\Pi A\Pi]^{\mathrm{TF}},\\
(\operatorname{Hess}F_4)^{\mathrm{TF}}
&=[8v\otimes v+4Y\Pi B\Pi-\Pi A\Pi]^{\mathrm{TF}}.
\end{aligned}
\]

Substitution leaves \(4[v\otimes v]^{\mathrm{TF}}\), precisely the
left-hand side of (12a). This is a direct tensor proof, also checked by
the independent exact verifier. Transverse
trace-free tensors are orthogonal to scalar-derived Hessians by integration
by parts on the closed sphere. Thus the TT projection of this particular
quadratic scalar angular stress is exactly \(-4\mathcal T/3\).

As a separate check, direct contraction gives

\[
\mathcal T_{ij}D^iY D^jY
=4\left[-\tfrac43n_1^2n_2^2
+\tfrac16(n_1^2+n_2^2)(n_3^2+n_4^2)\right].
\]

Since the normalized sphere moment of \(n_i^2n_j^2\), \(i\ne j\), is
\(1/24\),

\[
\frac1{\operatorname{Vol}(S^3)}
\int_{S^3}\mathcal T_{ij}D^iY D^jY=-1/9\ne0.
\tag{12}
\]

The same sphere-moment calculation gives
\(\operatorname{Vol}(S^3)^{-1}\int|\mathcal T|^2=1/12\), consistent
with the coefficient \(-4/3\) in (12a).

Thus the quadratic scalar angular stress \(D\psi_1\otimes D\psi_1\)
has a nonzero tensor-sector projection wherever its \(Y\) coefficient is
nonzero. This is **one component** of the full quadratic Einstein source;
gravitational terms may contribute to the same sector, and their total
projection has not been computed here. Consequently (12) does not prove
that the complete forcing is nonzero. It does show why importing the
homogeneous scalar-master Wentzell theorem alone is insufficient: one
must either calculate a cancellation or include the sourced tensor
sector with its junction conditions.

## 5. Verification and the remaining target

`verify_precontact_dynamics.py` imports none of the inherited scientific
code. It computes the five-dimensional Christoffels and Ricci tensor for
an independent time-dependent diagonal metric to check the universal ADM
kinetic coefficients and the parameter Taylor factors. It separately
checks the canonical density conversion, the two explicit quadrupole
forces, sphere tangency/trace/divergence, the full scalar-stress tensor
decomposition, the exact tensor-stress moment and tensor norm,
and the moving-shell support gaps. Its wrong-formula controls include
incorrect Taylor factors, an incorrect trace sign, omission of scalar
damping, an incorrect conformal density power, deletion of the tensor
sector, and a static-shell contact formula. Counts and completed normal,
optimized, and injected-failure runs belong in their receipts.

These algebra checks do not replace a hyperbolic PDE proof and do not
simulate shell contact. Proposition 1 uses the standard normally
hyperbolic Cauchy theorem and constraint propagation; (9) supplies the
complete initial-time acceleration once the inherited correction is
inserted. The next precise calculation is the full projection of
\(-D^2\mathcal E[z_1,z_1]\) into the finite scalar and tensor angular
sectors, together with the quadratic pullback, Israel, scalar, and
embedding sources at the moving shell. A full solution through contact
requires a compatible forced hyperbolic interface problem in every
nonzero sector. Nonlinear existence and convergence remain separate.
