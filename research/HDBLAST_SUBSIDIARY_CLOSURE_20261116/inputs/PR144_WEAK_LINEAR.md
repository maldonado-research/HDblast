# Weak linear energy solutions on a closed manifold are strongly continuous

14 November 2026 UTC. This is a qualitative regularity lemma for the
additional top-time-continuity argument. It does not change the evaluated
amplitude, Picard, difference or factorial-tail constants. Its boundary
application requires a separate justified doubling and commuted weak-row
argument.

## 1. Precise closed-manifold statement

Let M be a compact smooth manifold without boundary, with fixed smooth
positive measure and a finite-rank real metric bundle F with fixed smooth
metric connection. On a compact time interval I, suppose the distributional
system is

```
Y_s=M^i(s,x) nabla_i Y+L(s,x)Y+f,
Y in L^infinity(I;L2(F)),               f in L1(I;L2(F)).          (1)
```

The principal symbol M is W^(1,infinity) in time and space. It has a
symmetrizer S with

```
S=S*,      m Id<=S<=M0 Id,       m>0,
S in W^(1,infinity)(I x M),
S M^i=(S M^i)*,
L in L^infinity(I x M).                                        (2)
```

The repeated coordinate notation denotes the tensorial contraction of
the principal symbol with the fixed connection. Local connection and
volume terms are included in L or in covariant divergence. All bounds
are uniform on this compact interval. The conclusion is

```
Y has a unique representative in C(I;L2(F)).                    (3)
```

For that representative the ordinary energy identity holds for every
pair of times s0,s1:

```
E(s)=1/2 integral <S(s)Y(s),Y(s)>,
E(s1)-E(s0)=integral_(s0)^(s1) Q(r)dr,
Q=integral [1/2<S_sY,Y>-1/2<div(S M)Y,Y>
                         +<S Y,L Y>+<S Y,f>].                  (4)
```

No boundary term is present. If a non-coordinate fixed frame is used,
its smooth divergence/connection terms are retained in the same
covariant formula. In particular Q is integrable. No classical time
derivative of the original L2 solution is assumed.

## 2. First obtain the weakly continuous trace

Equation (1), integration by parts against a smooth spatial test, and
the spatial Lipschitz bounds imply

```
Y_s in L1(I;H^-1(F)).                                         (5)
```

For example `M^i nabla_i Y=div(MY)-(div M)Y`, with each zeroth-order
term bounded in L2 and the divergence bounded in H^-1. Integrating
this H^-1-valued derivative gives a continuous, absolutely continuous
H^-1 representative. The uniform almost-everywhere L2 bound extends
it to an L2 vector at every time: approach any time through regular
times, take a weak L2 subsequence, and identify its limit by the
continuous H^-1 value. Every subsequence gives the same value.

Pairing first with smooth tests and then using their density in L2
proves weak continuity in L2 at every time. Its L2 norm is no larger
than the essential bound. This supplies the initial and terminal weak
traces used below; no trace is recovered by asserting strong
continuity in advance.

## 3. The spatial Friedrichs commutator

Choose fixed, time-independent spatial regularizers J_epsilon on F,
uniformly bounded on L2, converging strongly to the identity there,
and smoothing H^-1 into smooth spatial sections for fixed epsilon.
They may be constructed in a finite smooth atlas with orthonormal
bundle trivializations, ordinary compactly supported convolution and
a quadratic partition of unity `sum chi_a^2=1`.

The required elementary commutator fact is

```
r_epsilon=J_epsilon(M^i nabla_i Y)
                      -M^i nabla_i(J_epsilon Y),
r_epsilon ->0 in L1(I;L2),
||r_epsilon(s)||2<=C ||M(s)||W^(1,infinity) ||Y(s)||2.           (6)
```

Here C is a finite fixed chart/mollifier constant. Its value is not
used in any numerical amplitude or evolution estimate: only uniform
boundedness and convergence to zero enter this lemma.

For completeness, in one Euclidean trivialization with convolution
kernel rho_epsilon, one coefficient entry a and one derivative,

```
[J_epsilon,a partial_j]h(x)
 =integral (partial_j rho_epsilon)(x-y)[a(y)-a(x)]h(y)dy
                          -J_epsilon((partial_j a)h)(x).
```

The first term is bounded on L2 by
`Lip(a) ||h||2 integral |z| |partial_j rho(z)|dz`, and the second
by `Lip(a)||h||2 ||rho||1`. These bounds are independent of epsilon.
For smooth h, the original commutator tends to zero in L2 because
convolution converges on the two L2 products. Approximate an arbitrary
h in L2 by smooth h and use the uniform bound to conclude the same
convergence for h. Matrices and the finitely many derivative directions
are finite sums of this calculation.

In a localized global construction, commutation with the smooth
cutoffs produces only zeroth-order terms. Their limiting sum vanishes
because `sum chi_a partial_i chi_a=0`; smooth frame/volume changes
are ordinary smooth lower-order multipliers. This proves the global
version. At almost every time it applies to Y(s). The uniform estimate
in (6), the uniform coefficient bounds and the L^infinity L2 bound
permit dominated convergence in time, proving its L1 L2 assertion.

Likewise, with merely bounded L,

```
J_epsilon(LY)-L J_epsilon Y ->0 in L1(I;L2),
J_epsilon f -> f in L1(I;L2).                                (7)
```

Both follow from strong L2 approximation at almost every time and
dominated convergence, using f in L1 L2 for its corresponding bound.

## 4. Regularized identity and its limit

Put `Y_epsilon=J_epsilon Y`. Equations (6)--(7) give

```
(Y_epsilon)_s=M^i nabla_i Y_epsilon+L Y_epsilon+f_epsilon,
f_epsilon ->f in L1(I;L2).                                   (8)
```

For each fixed epsilon the spatial smoothing of (5) makes
Y_epsilon absolutely continuous in time with values in L2. The
usual product rule with the Lipschitz-in-time S and one spatial
integration by parts therefore prove (4) for Y_epsilon. Its time
trace is precisely J_epsilon applied to the weak trace constructed
in Section2, since J_epsilon is continuous from H^-1 to L2 for fixed
epsilon.

For every time s the weak representative Y(s) is in L2, so
`J_epsilon Y(s)->Y(s)` strongly in L2. The regularized sections are
uniformly bounded in L^infinity L2. Thus each endpoint energy tends
to the corresponding energy in (4). The bounded quadratic production
terms converge in L1 by dominated convergence. For the forcing term,
split its difference into the L1 L2 error f_epsilon-f and the
pointwise strong approximation of Y, dominated by the integrable
norm of f. These observations pass the regularized identity to
(4) for every s0,s1, not only almost every endpoint.

It follows that E is absolutely continuous. We already have weak
L2 continuity. Since S is strongly continuous in the L^infinity
operator norm and uniformly positive, energy continuity gives strong
L2 continuity: fix s0, expand
`<S(s0)(Y(s)-Y(s0)),Y(s)-Y(s0)>`, use weak continuity for the cross
term, and use `||S(s)-S(s0)||infinity->0` for the change in the
quadratic term. The result tends to zero; its positive lower bound
is `m||Y(s)-Y(s0)||2^2`. This proves (3), including one-sided
continuity at the interval endpoints.

## 5. Consequence for a scalar-principal bundle wave

On the same closed product geometry, suppose

```
w_ss-A^i D_i w_s-B^ij H_ij w=g,
w in L^infinity H1,         w_s in L^infinity L2,
g in L1 L2,
A,B in W^(1,infinity),      B=B*,      b Id<=B<=b0 Id, b>0.
```

Introduce `v=w_s`, `V_i=D_i w`, and `Y=(v,V,w)`. The fixed product
connection is independent of time, so the gradient equation is
`(V_i)_s=D_i v` in distributions. The principal symbol has blocks

```
M^i = [[A^i, B^ij, 0], [delta_j^i, 0, 0], [0, 0, 0]],
S = diag(1,B,1).
```

The product S M^i is symmetric because B is symmetric. The full
Hessian's derivative-index connection enters the bounded lower-order
matrix; it is not omitted. The equation `w_s=v` supplies its lower
zero-order block, and the forcing is `(g,0,0)`. Hence (1)--(2) apply.
Strong continuity of all three blocks means

```
w in C(I;H1),                  w_s in C(I;L2).                 (9)
```

This proves energy regularity of a single weak wave solution under
the stated hypotheses. It is not a claim that a particular sequence
of nonlinear approximations converges strongly in top energy, nor a
claim of continuous dependence on top-order initial data. Mixed
boundary conditions require the separately proved extension and weak
boundary-row interface before this closed-manifold lemma can be used.
