# Primary check: linear homogeneous Neumann support

Prepared 14 October 2026 (America/Los_Angeles), retrieved 15 October UTC.
The issue is existence of smooth frozen linear solves after a proposed
Robin-to-Neumann field transformation. This note does not import an
unchecked quasilinear theorem.

## Shibata 1988

Yoshihiro Shibata, “On the Neumann problem for some linear hyperbolic
systems of second order,” Tsukuba Journal of Mathematics 12 (1988),
149–209. [Full primary PDF](https://core.ac.uk/download/pdf/56641220.pdf).
The full PDF was parsed. Mathematical extraction includes occasional
typographical/OCR artifacts; no potentially ambiguous printed trace norm
in the later energy formulas is used as a quantitative assertion here.

Section 1 problem (N), (1.1a)–(1.1c), has the form
\[
 P u=u_{tt}+2H^j\partial_j u_t+H^{n+1}u_t
       -\partial_i(A^{ij}\partial_j u)+A^j\partial_j u+A^{n+1}u=f,
\]
\[
 Q u=\nu_i A^{ij}\partial_j u+B^i\partial_i u+B^{n+1}u+H^0u_t=g.
\]
The fields are vectors of fixed finite dimension.

Section 2 assumptions, summarized without removing their substance:

1. (A.1): bounded C2 principal/drift/boundary coefficients and the stated
   lower-order regularity.
2. (A.2): a global-in-time coefficient extension to constant Euclidean
   principal part outside a finite time interval. Lemma 2.4 explains the
   extension of finite-time coefficients preserving the conditions.
3. (A.3): transpose symmetry A^(ij) transpose=A^(ji),
   H^j and H^0 symmetric, and B^i skew-symmetric.
4. (A.4): a coercive spatial bulk plus boundary bilinear form, up to
   a lower L2 term.
5. (A.5): nu_i B^i=0.
6. (A.6): the quadratic form nu_i H^i+H^0 is nonnegative.
   This is non-strict and includes the zero case.

Theorem 2.2(1) supplies unique E2 existence, with
\[
 E^L([0,T])=\bigcap_{j=0}^L C^j([0,T];H^{L-j}(\Omega)),
\]
for H2/H1 initial data, the stated strong time-dependent bulk/boundary
source classes and initial compatibility. In particular this is not
an existence theorem with arbitrary bare L2 boundary input in the H1
energy space.

Theorem 2.3 upgrades to E^L for L>=3 when the initial data are
H^L/H^(L-1), principal coefficients have the stated B^L regularity,
lower-order coefficients the stated B^(L-1) regularity, sources have
the displayed mixed space/time regularity, and the compatibility
conditions of Definition 2.1 hold through order L-2.

For a scalar-principal wave equation normalized by its negative
time-time coefficient, A^(ij)=a^(ij)I is symmetric and positive near
the registered background. With pure homogeneous conormal Neumann,
B^i=B^(n+1)=H^0=0. If the normal time drift vanishes at the shell,
(A.6) is an equality. The vanishing normal–tangential metric block
makes a scalar normal derivative equivalent to this conormal row.
First-order bulk matrix couplings may be retained as lower-order terms.

This provides direct support for smooth linear all-Neumann solves and
their finite-order regularity. A mixed set of Dirichlet and Neumann
field components with coupled first-order terms requires its own
combination argument, because Shibata's displayed problem (N) places
Neumann-type rows on every field. The theorem's coefficient norms
must also be matched to an actual tame Sobolev iteration before
calling it a quasilinear existence theorem.

## Why Ikawa 1969 was not substituted

Mitsuru Ikawa, “A mixed problem for hyperbolic equations of second order
with non-homogeneous Neumann type boundary condition,” Osaka Journal
of Mathematics 6 (1969), 339–374.
[Full primary PDF](https://projecteuclid.org/journals/osaka-journal-of-mathematics/volume-6/issue-2/A-mixed-problem-for-hyperbolic-equations-of-second-order-with/ojm/1200692528.pdf).

Section 1 (1.3) requires
\[
 \sigma_1 \le h\cdot\nu-\epsilon_0,\qquad\epsilon_0>0,
\]
for the stated boundary operator (1.4). Theorems 1–2 then give
existence and higher regularity for the scalar problem, including
specified boundary-source norms and compatibility. Pure undamped
Neumann with zero normal drift sits at equality, so it is outside
this theorem's displayed strict hypothesis. The introduction itself
distinguishes the earlier homogeneous equality-case paper.

## Flux check in the later wave literature

Metcalfe–Sogge–Stewart, arXiv:math/0411513v2, Section 5 (5.2),
requires the contraction of the principal perturbation between a
normal covector and a tangent covector to vanish at the Neumann
boundary. Under that condition the energy estimate (4.3) continues
to hold for (5.1). Hou–Yin–Yuan arXiv:2602.14019v1 (1.3) retains
the same structural requirement in the nonlinear coefficients.

Thus checking g^(za)=0 on the shell is essential. A Robin matrix
similarity or a scalar Neumann symbol calculation alone does not
replace that bulk/boundary structural check.
