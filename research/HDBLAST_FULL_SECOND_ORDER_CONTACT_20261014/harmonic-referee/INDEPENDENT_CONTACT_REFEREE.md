# Independent finite-angular wave-gauge contact review

14 October 2026. Independent internal mathematical review of the proposed
second-order coefficient construction. This is not external peer review.

**Verdict.** The proposed construction closes for a selected finite-angular,
smooth second-order coefficient on a controlled finite global-de-Sitter-time
slab, provided the inherited localized second-order initial-constraint lemma
and inherited smooth first-order contact solution are used with their actual
hypotheses. No fatal boundary, constraint, or analytic obstruction was found.
The argument below supplies the missing linear boundary construction; it
does not invoke the four-dimensional vacuum umbilic theorem as a theorem
about this model. It proves neither an exact nonlinear contact family nor a
uniform estimate over arbitrarily high angular frequencies.

## 1. Precise setting and claim

Use derivative Taylor coefficients, so
\(q_\epsilon=q_0+\epsilon q_1+\epsilon^2q_2/2\), with
\(q=(g,\phi)\). The five-dimensional residuals are
\[
E_{AB}=R_{AB}-\phi_{,A}\phi_{,B}-\tfrac23U(\phi)g_{AB},
\qquad S=\Box_g\phi-U'(\phi).
\]
The exact background is
\[
g_0=dy^2+\rho(y)^2\gamma,\qquad
\gamma=-dt^2+\cosh^2t\,\kappa_{S^3},\qquad \phi_0=\phi_0(y).
\]
At its fixed shell \(y=y_b\), increasing \(y\) is outward in both
identical copies. Set \(a_b=\rho_y/\rho=\sigma_0/6\) and
\(p_b=\phi_{0,y}=-\sigma'_0/2\). The identical-copy physical laws are
\[
K_{ab}=\sigma(\phi)h_{ab}/6,\qquad n\phi=-\sigma'(\phi)/2.
\tag{1}
\]
All background coefficients and their needed derivatives are smooth on the
closed causal collar in question, with \(\rho>0\).

Assume the selected aligned, copy-even \(l=2\) longitudinal first solution
is smooth on \(0\le t\le T\), solves every first bulk and physical shell
equation, and has compact initial support away from shell and inner endpoint.
Assume the inherited localized theorem gives a smooth second canonical
correction supported in a radial annulus away from both endpoints and
solving the full sourced Hamiltonian and momentum constraints. These are
actual inputs; neither follows merely from charge counting.

Then one can choose second coordinate data and construct a smooth, identical
second coefficient on this finite slab which solves all second bulk
equations and (1), with a fixed shell in a boundary-orthogonal wave chart.
The construction includes the scalar \(l=0,2,4\) and tensor \(l=2\)
representations needed by the selected first field. It makes no assertion
that the inherited single scalar master spans the second system.

## 2. Exact boundary reduction in the relative wave gauge

Define
\[
C^A=g^{BC}(\Gamma[g]^A{}_{BC}-\Gamma[g_0]^A{}_{BC}),\qquad
C_A=g_{AB}C^B,
\]
and solve \(E_{AB}-\nabla_{(A}C_{B)}=0\), \(S=0\). Parentheses include
one half. Impose four boundary coordinate rows \(g_{ya}=0\). At the shell
write \(g_{yy}=N^2\), \(g_{ab}=h_{ab}\), with \(N>0\). Tangential
derivatives of the zero mixed trace vanish. Therefore
\[
K_{ab}=\frac{1}{2N}\partial_yh_{ab},\qquad n\phi=N^{-1}\phi_y,
\]
and (1) becomes
\[
\partial_yh_{ab}=\frac{\sigma}{3}Nh_{ab},\qquad
\phi_y=-\frac{\sigma'}2N.
\tag{2}
\]
The reference connection satisfies
\(\widehat\Gamma^y{}_{ab}=-a_b(h_0)_{ab}\) and
\(\widehat\Gamma^y{}_{yy}=\widehat\Gamma^y{}_{ya}=0\). Direct contraction,
without assuming normal derivatives of the mixed metric vanish, gives
\[
C^y=\frac{N_y}{N^3}-\frac{1}{2N^2}h^{ab}h_{ab,y}
       +a_b\operatorname{tr}_h h_0.
\]
The single additional gauge row \(C_y=0\) is consequently exactly
\[
N_y=\frac{2\sigma}{3}N^2-a_bN^3\operatorname{tr}_h h_0.
\tag{3}
\]
It is satisfied by \(N=1,h=h_0,\sigma=6a_b\). Replacing relative gauge by
an absolute harmonic condition in this chart would fail that background
test.

Twice differentiating (2)–(3) and \(g_{ya}=0\) about the background gives
four Dirichlet and twelve ordinary normal-derivative/Robin conditions on
the fifteen metric components and the scalar. The derivative coefficient
on each of the twelve Robin variables is nonzero. After solving these
rows for normal derivatives, the remaining matrix is order zero. It need
not be symmetric, positive, or diagonally symmetrizable. Tangential
derivatives of the Dirichlet traces are known boundary sources and add no
incoming unknown. In particular, imposing all five components of \(C\)
in addition would be an invalid extra prescription.

## 3. A first-order gauge change with controlled boundary and source

For \(h^L_{1,yy}=-4\Psi\), \(h^L_{1,ya}=0\),
\(h^L_{1,ab}=2\Psi(h_0)_{ab}\), direct covariant differentiation gives
\[
C^{L}_{1,y}=-6(\Psi_y+4a\Psi),\qquad C^L_{1,a}=0.
\tag{4}
\]
Let \(P\xi=\Box_0\xi+\operatorname{Ric}_0(\xi)\), understood with
consistent vector or covector indices. Solve
\[
P\xi=-C^L_1,\qquad \xi_y|_b=0,\qquad
(\partial_y\xi_a-2a_b\xi_a)|_b=0,
\qquad \xi|_{t=0}=\partial_t\xi|_{t=0}=0.
\tag{5}
\]
The linear lemma in section 6 proves this problem as well. The first gauge
variation obeys \(DC[\mathcal L_\xi g_0]=P\xi\); thus
\[
q_1^W=q_1^L+\mathcal L_\xi q_0,\qquad C_1^W=0.
\]
The first condition in (5) makes \(\xi\) tangent to the shell. The second
makes \((\mathcal L_\xi g_0)_{ya}=0\) there. Hence all first physical
shell equations and the fixed-shell representation are preserved. The
zero initial value and first derivative of \(\xi\) preserve the canonical
initial data; its higher time derivatives change coordinate jets.

All second expressions must use \(q_1^W\). With a fixed first generator,
the conversion is
\[
B_g=2\mathcal L_\xi h_1^L+\mathcal L_\xi^2g_0,
\qquad B_\phi=2\mathcal L_\xi\chi_1^L+\mathcal L_\xi^2\phi_0.
\tag{6}
\]
Thus \(q_2^W=q_2^L+B\) up to a separately specified second pure gauge.
If \(J^L=-D^2(E,S)[q_1^L,q_1^L]\), the transformed physical source is
\(J^W=J^L+D(E,S)[B]\). Equivalently, recompute the defining Hessian at
\(q_1^W\); these procedures agree by covariance. Boundary and initial
expressions must be converted by the same pullback. A source from the
longitudinal bulk note cannot be attached unchanged to this wave problem.

Let \(Q_C=D^2C[h_1^W,h_1^W]\). Then
\[
C_2=DC[h_2]+Q_C,
\]
and the actual second reduced metric equation is
\[
DE[h_2,\chi_2]-\nabla^0_{(A}DC[h_2]_{B)}
=-D^2E[q_1^W,q_1^W]+\nabla^0_{(A}(Q_C)_{B)}.
\tag{7}
\]
Terms involving variations of the connection times \(C_1\) vanish.
The boundary row is the full \((C_2)_y=0\), including \((Q_C)_y\),
as is already guaranteed by differentiating (3).

## 4. Initial constraints, gauge jets, and angular closure

The radial background and its linear initial-constraint operator are
SO(4)-equivariant. Let \(V\) be the \(l=2\) scalar representation,
with spin labels \((1,1)\). Every selected quadratic source is the
image of the symmetric equivariant bilinear map on the same angular
input. Its representation lies in
\[
\operatorname{Sym}^2V=(2,2)\oplus(2,0)\oplus(0,2)
                       \oplus(1,1)\oplus(0,0).
\tag{8}
\]
These are the scalar \(l=4,2,0\) and the two tensor \(l=2\) chiral
representations. The character identity and dimensions give
\(45=25+5+5+9+1\). There is no coexact-vector representation here.
The wave-gauge transformation (5) is an equivariant linear map of the
same first \(Y\), so (6), all wave-gauge sources, and the boundary
Hessians obey the same statement. Independent homogeneous vector data
are excluded by selecting zero data in those representations.

Apply the compact-group isotypic projector \(\Pi\) for (8) to an inherited
localized second constraint correction. It commutes with the linear
constraint operator, leaves its quadratic source fixed, and preserves
the radial annular support. Its image in each relevant finite-rank tensor
bundle on \(S^3\) has finite angular multiplicity. This produces admissible
finite-angular initial data without assuming the inherited localized
solver itself preserves a componentwise harmonic truncation. No inverse
of an exceptional \(l=0\) harmonic eigenvalue is used.

Choose initial second lapse and shift values smoothly with localized
support. The physical metric time derivatives are fixed by the prescribed
canonical metric/momentum data. The five remaining time derivatives of
lapse and shift set \(C_2=0\) algebraically. For example, in an orthonormal
Cauchy frame at a point, with \(0\) timelike and \(i=1,\ldots,4\),
\[
\begin{split}
DC[k]_0&=-\tfrac12\partial_0k_{00}
 -\tfrac12\partial_0\sum_i k_{ii}+\partial_i k_{0i}+\text{lower},\\
DC[k]_i&=-\partial_0k_{0i}+\partial_jk_{ij}
 +\tfrac12\partial_i k_{00}-\tfrac12\partial_i\sum_j k_{jj}
 +\text{lower}.
\end{split}\tag{9}
\]
The coefficient matrix for the five unknown temporal derivatives is
invertible. Include \(Q_C\) and every known lower term in the right side.
The operation is equivariant and local and preserves initial collar
vanishing. There is no change of physical canonical data.

The physical constraints concern
\(H_{AB}=E_{AB}-g_{AB}\operatorname{tr}_gE/2\), not \(E_{00}\) alone.
On the initial slice, \(C_2=0\) implies every spatial derivative of
\(C_2\) is zero. The second reduced equation consequently gives
\[
(H_2)_{00}=\tfrac12\partial_0(C_2)_0,\qquad
(H_2)_{0i}=\tfrac12\partial_0(C_2)_i.
\tag{10}
\]
The sourced canonical constraints therefore imply the zero normal
derivative of every \(C_2\) component. Since all perturbative data vanish
on an initial shell collar, recursively using the reduced equations gives
zero corner jets of every order, both for (5) and for (7).

## 5. Exact subsidiary boundary equations

For the convention in (1), Codazzi gives
\[
E_{na}=D_bK^b{}_a-D_aK-(n\phi)D_a\phi
       =-3D_a(\sigma/6)-(n\phi)D_a\phi=0.
\tag{11}
\]
The cancellation uses the scalar junction. It does not assert
\(R_{na}=0\) when \(\sigma\) varies.

For the covariant reduced equation, (11) gives
\((\nabla_nC)_a+(\nabla_aC)_n=0\). Since \(C_n=0\) as a boundary
trace, \((\nabla_aC)_n=-K_a{}^bC_b\). In the boundary-orthogonal chart,
\[
C_y=0,\qquad \partial_y C_a-2NK_a{}^bC_b=0.
\tag{12}
\]
The factor two comes from the two connection terms; a literal coordinate
Neumann condition would be wrong for this definition of \(C\).

The exact differential identity is
\[
\nabla^AE_{AB}-\tfrac12\nabla_B\operatorname{tr}_gE
   =-S\phi_{,B}.
\]
Together with the reduced system it yields
\(\Box_g C_B+R_B{}^A C_A=0\). Taking second parameter derivatives,
using solved zeroth/first equations and \(C_0=C_1=0\), gives the
homogeneous linear equation
\[
(\Box_0+\operatorname{Ric}_0)C_2=0,
\qquad (C_2)_y|_b=0,\qquad
(\partial_y(C_2)_a-2a_b(C_2)_a)|_b=0.
\tag{13}
\]
All apparently inhomogeneous variations of (12) multiply \(C_0\) or
\(C_1\) and vanish. Equation (10) gives zero Cauchy data. Uniqueness
from the next section proves \(C_2=0\), hence the unreduced equations.
The inherited exact shell normal-force identity may then be differentiated:
its second coefficient vanishes from the bulk, matching, Israel, and
scalar residuals. It supplies a consistency consequence, not an additional
incoming boundary row.

## 6. Linear characteristic lemma, including a nonsymmetric Robin matrix

Use \(x=\log(r/b)\le0\), where \(dy=\rho\,dx\). On each finite angular
space all three systems used above have scalar principal part. After
multiplication by a nonzero scalar and smooth component rescaling they are
\[
u_{tt}-u_{xx}=A(t,x)u_t+B(t,x)u_x+D(t,x)u+f(t,x).
\tag{14}
\]
The matrices are smooth and bounded on the finite causal collar. Split
\(u=(d,r)\) according to Dirichlet and Robin components. A smooth lift
removes a prescribed Dirichlet trace. The remaining boundary conditions
can be written
\[
d(t,0)=0,\qquad r_x(t,0)=M(t)r(t,0)+g(t).
\tag{15}
\]
Any coupling to Dirichlet values has been absorbed in \(g\). No symmetry
or sign is imposed on \(M\).

For completeness, a direct construction follows. Put
\(p=u_t+u_x\), \(q=u_t-u_x\), and
\(F=A(p+q)/2+B(p-q)/2+Du+f\). Then
\[
(\partial_t-\partial_x)p=F,\qquad
(\partial_t+\partial_x)q=F.
\tag{16}
\]
The \(q\) characteristics, followed backwards, always meet the initial
surface in \(x\le0\); their trace \(q(t,0)\) is therefore given by a
Volterra integral of \(F\). The \(p\) characteristics either meet the
initial surface or the boundary. If \(z(t)=r(t,0)\), (15) is exactly
\[
z'=Mz+q_R(t,0)+g,\qquad
p_R(t,0)=q_R(t,0)+2Mz+2g,\qquad
p_D(t,0)=-q_D(t,0).
\tag{17}
\]
Solve the finite-dimensional ODE in (17) for its initial trace **before**
forming a contraction. Its fundamental matrix obeys
\(\|\Phi(t,s)\|\le\exp(\int_s^t\|M(v)\|dv)\), including for
nonsymmetric \(M\).

Specifically, for a candidate continuous triple \((u,p,q)\), first
compute \(F\), then compute the new outgoing \(q_*\) by (16), then
solve (17) with \(q_*\), then compute the incoming \(p_*\) by (16),
and finally put \(u_*=u_0+\tfrac12\int_0^t(p_*+q_*)ds\). For two
candidate triples with the same initial data, on a slab of length
\(\delta\),
\[
\|q_*^{(1)}-q_*^{(2)}\|_\infty\le c\delta\|X^{(1)}-X^{(2)}\|_\infty,
\]
\[
\|z_*^{(1)}-z_*^{(2)}\|_\infty
 \le\delta e^{\delta\|M\|_\infty}\|q_*^{(1)}-q_*^{(2)}\|_\infty,
\]
and hence
\(\|X_*^{(1)}-X_*^{(2)}\|_\infty\le c'\delta
 e^{c'\delta}\|X^{(1)}-X^{(2)}\|_\infty\).
The reflection coefficient one is harmless because the new outgoing
field, already containing the factor \(\delta\), is used in (17).
Using the candidate outgoing field directly at this point would fail
to prove contraction in the stated norm.

For small \(\delta\), Banach's theorem produces the continuous fixed
point. It solves the first-order equations distributionally. The identity
\(\partial_t[u_x-(p-q)/2]=0\), its zero initial value, and
\(u_t=(p+q)/2\) imply \(u\in C^1\) with the expected derivatives.
The boundary ODE and initial boundary compatibility ensure (15).
Smoothness follows by time differentiation: each differentiated problem
has the same principal/Robin operators and known lower derivative sources.
Difference quotients and uniqueness justify those derivatives successively;
normal derivatives are recovered from (16), or from (14) after two time
derivatives. Smooth compatible data consequently give a smooth solution.

The same estimates give uniqueness and a finite-slab bound in a
\(C^1\)-type characteristic norm; corresponding finite derivative norms
follow after differentiation. They are not an unproved coercive physical
energy statement. Iterating the local construction covers any prescribed
finite \(T\) on which the coefficients are bounded. A nonsymmetric or
negative Robin coefficient may cause finite exponential growth but cannot
obstruct this finite-time linear construction.

The characteristic formulas also prove finite propagation. Extend
coefficients smoothly outside a compact causal collar if desired; no value
of that extension can enter the solution in the collar's domain of
dependence. Initial radial support in \(x\ge x_-\) and local forcing
generated by the first solution give support in \(x\ge x_- -t\).
Thus on finite \(0\le t\le T\) the perturbations remain in
\(r\ge b e^{x_- -T}>0\). This does not give a result at unbounded
de Sitter time or across every nonvertex point of the null cone.

## 7. Application and limits of the verdict

Apply the lemma first to (5), then to the finite-angular second reduced
system with boundary (2)–(3), then to the homogeneous subsidiary system
(13). The first two produce the fields; the third proves that the second
field is a solution of the full physical equations. Equal data and equal
boundary laws give identical copies. The common tension is used throughout.
The construction can be continued to any controlled
\(t_c<T<\operatorname{arsinh}(\tan\tau_*)\) supplied by the inherited
first contact solution. Smooth compact support gives zero corner jets;
contact at the support edge need not have a nonzero shell trace.

This is a complete selected second Taylor coefficient through contact,
conditional on the two inherited input theorems stated in section 1. It
does not prove nonlinear convergence, exact nonlinear continuation,
uniform angular Sobolev estimates, a two-master diagonalization, late-time
stability, or a cosmological interpretation. No such result is needed for
this finite-angular coefficient theorem.

The companion independent checker verifies the relative-gauge contraction,
subsidiary connection factor, longitudinal first gauge defect, initial
constraint factors, angular character decomposition, and a nontrivial
nonsymmetric Robin trace ODE. Its explicit corruption switches change the
tested candidate expressions. Passing those algebraic checks is evidence
for the identities; the existence proof is the argument in section 6.
