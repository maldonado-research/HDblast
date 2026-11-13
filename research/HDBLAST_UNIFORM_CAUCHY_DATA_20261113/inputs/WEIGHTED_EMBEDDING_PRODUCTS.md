# All-mode weighted pointwise and product estimates

13 November 2026 UTC. Proof for the selected weighted initial
constraint map. The actual smooth background and its authenticated
coefficient box are retained premises. The constants below concern
arbitrary tensors and arbitrary angular frequencies. They do not assume
the finite angular support of the registered second source, and do not
admit a nonlinear amplitude.

## 1. Exact geometry, norms, and quantitative conclusions

Use the fixed initial metric and specified member

```
gamma0=dy^2+rho(y)^2 kappa_S3, dy=-rho dz, rho_y=a rho,
129<=rho<=131, 0<a<=112/1000,
ell=10^-6, kappa_w=10^-8,
r=(z-ell/2)(5ell/2-z)/ell^2,
q=kappa_w^2 r^2, w=exp(-1/r),
ds=dy/q, D=partial_s, alpha=q/rho, beta=qa, d=q_y.
```

All covariant derivatives below use gamma0, acting also on earlier
derivative slots. Norms are full physical Hilbert/Frobenius norms, and
the volume is the full `dmu_gamma0=rho^3 |dy| dOmega`, with
`vol(S3)=2pi^2`. Contravariant slots may first be lowered by gamma0;
this is parallel and isometric. The notation for a general positive
weight omega is

```
W_m,omega(T)^2=sum_(j=0)^m integral omega^2 q^(2j)
                                  |nabla^j T|^2 dmu_gamma0.
```

The dot, when used, means compact-smooth closure in this norm.
The proof applies to every original tensor rank at most 32. It counts
the additional three derivative slots explicitly, rather than asserting
a rank-independent tensor embedding.

The evaluated pointwise constant is

```
|T(y,n)| <=10^7 [omega(y)q(y)^2]^-1 W_3,omega(T),       (1)
```

for omega=w and for each of

```
omega_d=w^-1 q^(d-2), d=0,1,2.
```

The letter d in omega_d is an integer grade, distinct from the radial
coefficient d=q_y when that coefficient is written explicitly below.
For clarity in the proofs, the latter is denoted `d_q` from section 2
onwards.

Define the graded physical norms, as literal sums of L2 norms,

```
E_d(T)=sum_(j=0)^(14-d)
          ||w^-1 q^(j+d-2) nabla^j T||_L2(gamma0),
d=0,1,2.                                               (2)
```

The completed interfaces are

```
|nabla^j T| <=10^12 w q^(-j-d) E_d(T),
                j+3<=14-d, rank(T)+j<=32,              (3)

E_(d1+d2)(U tensor V) <=10^17 E_d1(U) E_d2(V),
                d1,d2>=0, d1+d2<=2,
                rank(U),rank(V)<=25,                  (4)

||w^-2 Q||_H12(q,w) <=10^18 E_2(Q).                    (5)
```

For applications, original tensors are scalars, one-forms or rank-two
metric/momentum tensors, so all ranks in (3)--(5) fit well inside the
stated rank bound. The input rank bound 25 in (4) allows all seven
possible low-factor derivative slots within the embedding bound 32.
A fixed parallel contraction multiplies (4) by
its actual tensor operator norm. Cross-contractions such as matrix
multiplication have norm at most one; a trace on a single factor in
dimension four has norm two. No contraction factor is silently omitted.

## 2. Background coefficients needed by the proof

Set

```
d_q=q_y, c_w=D log(w sqrt(qrho^3)),
l_w=D log w=kappa_w^2 r_y,
h_alpha=D log alpha=d_q-beta.
```

The explicit polynomial coefficient calculation retained in
`HDBLAST_WEIGHTED_REGULARITY_20261106/SCALAR_REGULARITY.md` and proved
again in the accepted all-mode `WEIGHTED_NORM_BRIDGE.md` uses
x=(z-3ell/2)/ell and b=rho^-1. It gives

```
r=1-x^2, r_y=2xb/ell,
d_q=2kappa_w^2 r r_y,
beta=kappa_w^2 r^2 a,
c_w=kappa_w^2(1+r)r_y+(3/2)beta.
```

Their polynomial coefficient sums obey

```
|d_q|_coeff<=8*10^-10,
|beta|_coeff<=4*10^-16,
|l_w|_coeff<=2*10^-10,
|c_w|_coeff<=6.000006*10^-10.
```

The actual polynomial D-vector field has coefficient sum at most
10^-8 and degree at most nine. For a polynomial initially of degree
D0 and coefficient sum C0, its kth D derivative has coefficient sum
at most
`C0 product_(i=0)^(k-1)[10^-8(D0+8i)]`.
All factors needed here, through order sixteen, are below one.
This proves, on the whole annulus including both tails,

```
max_(k<=16)(|D^k d_q|,|D^k beta|,|D^k l_w|,|D^k c_w|)<=1,
|h_alpha|<10^-8, 0<alpha<1/2, |alpha|,|beta|<=1.        (6)
```

For a general weight used in (1), write
`M_omega=omega sqrt(qrho^3)` and `c_omega=D log M_omega`.
For omega=w, c_omega=c_w. For omega_d in (2),

```
c_omega_d=-l_w+(d-3/2)d_q+(3/2)beta,
|c_omega_d|<10^-8<1.                                   (7)
```

Only this zeroth bound on c_omega is needed for the basic embedding.
The higher multiplier estimates use the derivative version already
supplied by the same polynomial calculation.

The complete s coordinate permits an interval
`I=[s0-1,s0+1]` around every interior point, with no endpoint
extension. Equation (6) gives, for 0<=b<=2,

```
alpha(s0)^b <=exp(2*10^-8) alpha(s)^b<2alpha(s)^b,
s in I.                                                (8)
```

For example `exp(t)<=1/(1-t)` for 0<=t<1 proves the last comparison.
The analogous q,w and omega_d ratios on I are also below two, from
their displayed logarithmic derivatives. Their ratios are not required
in the proof: M_omega transports their measures exactly. No cutoff or
scaled-cover overlap constant is hidden here; a direct interval trace
inequality is used instead.

## 3. Scalar sphere kernel with the full volume normalization

Let Y_lm be full-L2-unit scalar spherical harmonics on the unit S3.
Their eigenvalue is lambda_l=l(l+2) and their multiplicity is
(l+1)^2. Rotation invariance and integration of the diagonal projection
kernel give the addition theorem

```
sum_m |Y_lm(n)|^2=(l+1)^2/(2pi^2).
```

For 0<alpha0<=1/2, the evaluation kernel associated with
`1-alpha0^2 Delta_S3` is

```
K(alpha0)=1/(2pi^2) sum_(l>=0)
             (l+1)^2/[1+alpha0^2 l(l+2)]^2.
```

Put n=l+1 and N=floor(alpha0^-1)>=2. Since
`1+alpha0^2(n^2-1)>=(1/2)(1+alpha0^2n^2)`,

```
sum_(n<=N) n^2/(1+alpha0^2n^2)^2 <=N^3<=alpha0^-3,
sum_(n>N) n^2/(1+alpha0^2n^2)^2
 <=alpha0^-4 sum_(n>N)n^-2
 <=alpha0^-4/N<=2alpha0^-3.
```

The integral comparison for the decreasing n^-2 tail gives the middle
inequality. Therefore

```
K(alpha0)<=12/(2pi^2) alpha0^-3<alpha0^-3.              (9)
```

The estimate uses every scalar harmonic, with no invariant or angular
cutoff. For a function v taking values in any finite-dimensional
Hilbert space, apply weighted Cauchy--Schwarz to each component and
sum. This proves

```
|v(n)|^2<=alpha0^-3
             ||(1-alpha0^2 Delta_S3)v||_L2(S3)^2.       (10)
```

There is no multiplication by the Hilbert-space dimension: that
dimension is already summed in both squared norms. This fact will
permit tensor components, while their frame connection costs still
depend explicitly on tensor rank.

## 4. Tensor frames and complete derivative-slot accounting

Identify the unit S3 with unit quaternions. Its global orthonormal
left-invariant frame E1,E2,E3 has

```
[E_i,E_j]=2epsilon_ijk E_k,
nabla^S_(E_i)E_j=epsilon_ijk E_k.
```

Use the physical orthonormal frame e0=partial_y, e_i=rho^-1 E_i.
It is radially parallel. The physical connection is

```
nabla_e0 e_A=0,
nabla_ei e0=a e_i,
nabla_ei e_j=rho^-1 epsilon_ijk e_k-a delta_ij e0.
```

Let a rank-R tensor be represented by its full array of 4^R frame
components. That identification is an isometry even for tensors with
symmetries, since the full Frobenius norm counts all components.
For angular differentiation its component connection matrices are
`alpha L_i^(S,R)+beta L_i^(rad,R)`. Each one-slot rotation has
operator norm at most one. Acting on every slot, each matrix has
norm at most R, and the aggregate over the three directions obeys

```
||(alpha L_i^(S,R)+beta L_i^(rad,R))_i||op
 <=sqrt(3)R(|alpha|+|beta|)<=4R.                        (11)
```

All matrices and index projections are constant in this global
component frame. The scalar E_i derivatives are kept in their
original order; their noncommutativity is not discarded.

For the original rank-r tensor T define component arrays

```
U_j=M_omega q^j (nabla^j T)_components, j=0,1,2,3,
v=U0.
```

The rank of U_j is r+j, including all derivative indices. With
Pi_A selecting the newly added derivative index, the exact recursions
are

```
D U_j=Pi_0 U_(j+1)+(c_omega+j d_q)U_j,
alpha E_i U_j=Pi_i U_(j+1)-A_i^(r+j)U_j,
A_i^R=alpha L_i^(S,R)+beta L_i^(rad,R).                 (12)
```

Angular derivatives of c_omega, alpha and beta vanish. In particular,
using |c_omega|<=1 and putting a_R=4R, (11)--(12) give the following
aggregate Hilbert norm bounds at each point:

```
|v| <=|U0|,
|Dv| <=|U1|+|U0|,
|alpha E v| <=|U1|+a_r|U0|,
|alpha E Dv| <=|U2|+(a_(r+1)+1)|U1|+a_r|U0|,
|alpha^2 E E v|
 <=|U2|+(a_(r+1)+a_r)|U1|+a_r^2|U0|,
|alpha^2 E E Dv|
 <=|U3|+(a_(r+2)+a_(r+1)+1)|U2|
      +[a_(r+1)^2+a_(r+1)+a_r]|U1|+a_r^2|U0|.         (13)
```

Every E or E E array includes the sum over all three or nine ordered
angular directions. The aggregate bound (11), rather than a bound
on only one direction, is used at every step. For r<=32, the sums
of coefficients in these six lines are at most

```
1, 2, 129, 262, 16645, 34338,
their total =51377<10^5.                               (14)
```

Unitarity of the factor M_omega is exact:
`M_omega^2 ds dOmega=omega^2 dmu_gamma0`.
Therefore the sum of the L2 norms of all six arrays in (13), on any
subset of the complete cylinder, is at most `10^5 W_3,omega(T)`.

## 5. Radial trace and the weighted tensor embedding

For a Hilbert-valued smooth function F on I=[s0-1,s0+1], the fundamental
theorem and Cauchy--Schwarz give, for every t in I,

```
|F(s0)|^2<=2|F(t)|^2+2|t-s0| integral_I |F'|^2.
```

Average over the interval of length two. Since
`integral_I |t-s0|dt=1`,

```
|F(s0)|^2<=integral_I (|F|^2+|F'|^2).                  (15)
```

Apply this in the Hilbert space of full-sphere L2 component arrays,
with `F=(1-alpha0^2 Delta_S3)v` and alpha0=alpha(s0).
For a scalar component,

```
(nabla_S^2 v)_ij=E_i E_j v-epsilon_ijk E_k v,
|nabla_S^2 v|<=|E E v|+2|E v|,
|Delta_S3 v|<=sqrt(3)|nabla_S^2 v|.
```

The same statements hold after square summing all tensor components.
Equations (8), (13)--(15) show that the square root of the right
side of (15) is bounded by eight times the sum of the six array L2
norms in (13), hence by `8*10^5 W_3,omega(T)`. Here the factors
are explicit: sqrt(3)<2, alpha0<=1, the angular-Hessian correction
costs two, and each alpha0^b/alpha(s)^b comparison costs at most two.

Combine with (10):

```
|v(s0,n)|<=8*10^5 alpha0^-3/2 W_3,omega(T).
```

Finally the dimensional cancellation is exact,

```
alpha^-3/2/M_omega
 =(q/rho)^-3/2/[omega sqrt(qrho^3)]
 =1/(omega q^2).
```

Thus (1) follows with the rounded constant 10^7. The argument
counts every tensor slot through (11)--(14), includes all angular
components through (9), and never imports an unspecified local
Sobolev or cover constant.

## 6. Covariant weight multipliers and higher pointwise jets

A useful finite multiplier ledger can be proved in the same frame.
Let f(y)>0 be radial and define

```
R_k=q^k f^-1 nabla^k f, ell_f=D log f.
```

Its frame components are rotation-invariant tensors built from the
radial covector and the angular metric. Their scalar E_i derivatives
vanish and their pure sphere-rotation connection action is zero.
The exact radial/angular covariant recurrence is

```
(R_(k+1))_0=D R_k+(ell_f-k d_q)R_k,
(R_(k+1))_i=beta L_i^(rad,k) R_k.                       (16)
```

The signs of the slot rotations do not affect the norm ledger.
If all needed D derivatives of ell_f,d_q,beta have absolute value
at most one, expand R_k as monomials in those profiles and their
derivatives. A term has at most k profile factors. Differentiating
its coefficient costs at most k, multiplication by ell_f-k d_q
costs at most 1+k, and the aggregate rank-k connection action costs
at most 3k. Thus the coefficient sum grows by at most
`1+5k<=5(k+1)` at the next step. Starting at R0=1 gives

```
|q^k f^-1 nabla^k f|<=5^k k!.                          (17)
```

For f=q^j, ell_f=j d_q. For j<=16 and k<=3 the same calculation
has next-step cost at most `j+5k<=26<32`, since the last step has
k=2. Therefore

```
|q^k q^-j nabla^k(q^j)|<=32^k,
0<=k<=3, 0<=j<=16.                                    (18)
```

Apply the full tensor Leibniz formula to `S_j=q^j nabla^j T`.
All slot permutations are isometries, and no covariant derivatives
within a factor are commuted. Equation (18) gives

```
W_3,omega(S_j)
 <=[sum_(k=0)^3 sum_(l=0)^k binom(k,l)32^l]
                                W_(j+3),omega(T)
 =37060 W_(j+3),omega(T)<10^5 W_(j+3),omega(T).          (19)
```

For E_d, `W_(j+3),omega_d(T)<=E_d(T)` whenever j+3<=14-d.
Use (1), (19) and `omega_d q^2=w^-1 q^d` to prove (3).
The allowable low jets are, in particular, j<=11 for grade zero
and j<=10 for grade one. Their derivative-slot ranks are included.

For the original potentials, the same calculation with omega=w gives

```
|nabla^j N|<=10^12/(w q^(j+2)) ||N||_H16(q,w), j<=13,
|nabla^j Z|<=10^12/(w q^(j+2)) ||Z||_H14(q,w), j<=11.    (20)
```

These estimates cover arbitrary scalar and one-form inputs, including
every Gamma-invariant angular mode, rather than only the selected
finite-channel second derivative.

## 7. Low/high products with the complete weights

Take nonnegative grades d1,d2 with d=d1+d2<=2. At an external
derivative order m<=14-d, a Leibniz allocation has i derivatives
on U and k=m-i on V. The total intrinsic-plus-external order is

```
(i+d1)+(k+d2)=m+d<=14.
```

At least one of these two orders is at most seven. Put that factor
in (3) and the other in its actual weighted L2 norm. The low
embedding is available because its order plus three is at most
ten, below its retained 14-derivative graded budget, and its tensor
rank including derivative slots is at most 25+7=32. If U is low,
the exact weight product is

```
[w^-1 q^(m+d-2)]
 *[w q^(-i-d1)]
 *[w q^(2-k-d2)] = w.
```

The same identity holds when V is low. Since 0<w<=1, every
allocation is bounded by `10^12 E_d1(U)E_d2(V)`. This calculation
includes all lower external orders m, not just the top derivative.
Summing all binomial allocations and all exterior derivative norms,

```
E_d(U tensor V)
 <=10^12[sum_(m=0)^(14-d) 2^m] E_d1(U)E_d2(V)
 <=10^12(2^15-1) E_d1(U)E_d2(V)
 <=10^17 E_d1(U)E_d2(V),                               (21)
```

which is (4). The norms in (2) are literal sums, so no omitted
square-sum-to-sum factor is present. Fixed parallel contractions
commute with all covariant derivatives and multiply the result by
their explicit operator norm.

The derivative maps themselves have norm at most one:

```
E_(d+1)(nabla T)<=E_d(T), d=0,1,
E_2(nabla^2 T)<=E_0(T).                                (22)
```

This follows by reindexing the literal sums; it does not commute
covariant derivatives. Consequently the required constraint products
have the following direct interfaces, before any fixed contractions:

```
E2(h nabla^2 k) <=10^17 E0(h)E0(k),
E2(nabla h nabla k) <=10^17 E0(h)E0(k),
E2(tau sigma) <=10^17 E1(tau)E1(sigma),
E2(h nabla tau) <=10^17 E0(h)E1(tau),
E2(tau nabla h) <=10^17 E1(tau)E0(h).                   (23)
```

The same estimates hold for independent state and direction tensors.
They are bilinear bounds, not estimates for one selected second source.
Repeated products are obtained by explicit iteration of (21), retaining
one factor 10^17 per multiplication and the total grade, provided each
input rank at every iteration remains at most 25. The actual contracted
scalar, matrix and one-form products keep their ranks below this limit;
no assertion for unbounded-rank uncontracted tensor powers is needed.

For a background coefficient B that need not itself belong to E0,
if `max_(j<=14-d) sup |q^j nabla^j B|<=C_B`, then the same Leibniz
calculation gives

```
E_d(B tensor T)<=(2^(15-d)-1) C_B E_d(T)
                <=10^5 C_B E_d(T).                    (24)
```

Again a contraction has its own operator norm. The constant identity
tensor and inverse-metric expansions can be treated in a unital
algebra; this proof does not incorrectly place a nondecaying background
identity in the exponentially weighted E0 space.

## 8. All derivatives of the exact source normalization

Let Q be a scalar or one-form physical constraint expression of grade
two. To estimate the literal normalized source w^-2 Q, take f=w^-2
in (16)--(17). Here

```
ell_f=D log(w^-2)=-2l_w,
max_(k<=12)|D^k ell_f|<=1
```

by the small polynomial coefficients of section 2. Thus

```
|w^2 q^k nabla^k(w^-2)|<=5^k k!, k<=12.                (25)
```

After m exterior derivatives, a term with k derivatives on w^-2
is bounded exactly by

```
|w^2 q^k nabla^k(w^-2)|
       *|w^-1 q^(m-k) nabla^(m-k)Q|.
```

This uses the actual source weight w q^m, so it retains all
derivatives of the normalization rather than treating w^-2 as a
constant. The H12 square-sum norm is bounded by the larger literal
sum. Consequently

```
||w^-2 Q||_H12(q,w)
 <=[sum_(m=0)^12 sum_(k=0)^m binom(m,k)5^k k!] E2(Q),
sum_(m=0)^12 sum_(k=0)^m binom(m,k)5^k k!
 =145260177974605953<10^18.                             (26)
```

This proves (5). Combining (5) and (23) supplies, for example,
`||w^-2 h nabla^2 k||_H12(q,w)<=10^35 E0(h)E0(k)`
before fixed contractions and coefficient factors. The other four
products have the same 10^35 basic normalized bilinear ceiling.

## 9. Completion, scope, and use by the nonlinear certificate

All proofs begin on smooth compact tensors and preserve their actual
physical derivative ordering. Estimate (1) gives locally uniform
convergence on the open annulus for a weighted H3 Cauchy sequence,
because omega q^2 is bounded below on each compact subset. It therefore
passes to the continuous representative in the weighted completion.
The higher-jet statement follows in the same way for S_j.

The bounded bilinear products extend by density to the corresponding
graded completions. Compact tests normalize to compact tests under
w^-2, and (26) then places the limiting source in the actual dotted
H12 range. No endpoint trace of a weighted potential is assumed.
Full sphere and reflected-copy square sums introduce no additional
operator factor when used consistently on both sides.

These are the pointwise, product, coefficient-multiplier and source-
normalization interfaces needed for a uniform second derivative of
the actual constraint map. The general correction operators and the
compact affine seed must still be inserted with their own evaluated
E0/E1 bounds. Inverse-metric and determinant derivatives must still
be controlled on the entire chosen positive ball, and the fixed
canonical density rho0 must remain fixed in that assembly.

No constant here alone is a numerical amplitude radius or an evolved
observer/PDE tube. Frozen prior notes, the public worktree, and external
services are not modified by this private proof.
