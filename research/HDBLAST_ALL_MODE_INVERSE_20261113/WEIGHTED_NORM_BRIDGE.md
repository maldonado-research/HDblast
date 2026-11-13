# An explicit all-mode bridge for the actual weighted norms

13 November 2026 UTC. Reviewed norm-comparison proof. This is a norm-comparison lemma,
not an inverse estimate or an amplitude-radius theorem. It applies to the
specified current annular weight member and its background. All original
tensor slots and all covariant derivative slots use the physical positive
metric and its full Hilbert/Frobenius norms.

## 1. Statement and conventions

Use the actual initial metric and weight member

```
gamma0=dy^2+rho(y)^2 kappa_S3,  dy=-rho dz,  129<=rho<=131,
ell=10^-6, x=(z-3ell/2)/ell, r=1-x^2,
kappa=10^-8, q=kappa^2 r^2, w=exp(-1/r),
ds=dy/q, M=w sqrt(q rho^3), alpha=q/rho, beta=q a,
a=rho_y/rho, D=partial_s.
```

The letter kappa without a subscript in `kappa=10^-8` is the weight
parameter; `kappa_S3` is the round unit metric. All integrals use the full
round sphere measure. Angular harmonics and one-forms below are L2-unit
in that full measure, not in probability measure. Since q vanishes
quadratically at both endpoints, s ranges over the complete real line,
up to orientation and translation.

For an original scalar or one-form T, let

```
W_m(T)^2=sum_(j=0)^m integral_A w^2 q^(2j)|nabla_gamma0^j T|^2 dvol_gamma0.
```

The dotted space is the closure of compact smooth tests in this norm.
Let I express tensors in the radially parallel orthonormal splitting
`R dy + T*S3`, sending the physical coframe `(dy,rho theta)` to
`(dy,theta)`, on **every** original and derivative slot. This is a fiber
isometry; it is not the coordinate substitution `dy=-rho dz`. Put
`v=M I(T)` and decompose its original slots into:

* scalar harmonics E with `-Delta_S3 E=lambda E`;
* the radial one-form sector `E dy` and exact unit one-form
  `omega=dE/sqrt(lambda)`, using the same two-vector coefficient v;
* unit coexact one-forms with `Delta_H omega=mu omega`.

For each nonconstant scalar/poloidal block set
`Lambda=lambda>=8`, and for each retained coexact block set
`Lambda=mu>=9`. The scalar and radial monopoles have `Lambda=0`.
These are exactly the allowed low-end bounds after the stated reflection
restriction: degree-one scalar modes and the Killing one-form sector have
no invariant part. Write `zeta=alpha sqrt(Lambda)` and define

```
J_m(v)^2=sum_blocks sum_(a+b<=m) ||zeta^b D^a v_block||_L2(ds)^2.       (1)
```

For the monopole, terms with b>0 are zero. Vector coefficients in (1)
use their Euclidean Hilbert norm. The same definitions apply to a source
scalar or one-form.

For every integer `0<=m<=16`, the two explicit bounds are

```
W_m(T) <= sqrt(m+1) 6^m m! J_m(v),
J_m(v) <= sqrt((m+1)(m+2)/2) 8^m m! W_m(T).               (2)
```

Consequently a single evaluated constant works in both directions at
every order through sixteen:

```
Cnorm=13*8^16*16!
     =76560143353807735162404864000 <10^29,
Cnorm^-1 J_m(v) <= W_m(T) <= Cnorm J_m(v).                (3)
```

No cutoff on Lambda occurs. There is no finite-angular-polynomial
derivative estimate in this proof. The common constant also works for
Hilbert direct sums of scalar, poloidal and coexact modes. If a later
interface replaces a two-summand Hilbert norm by the sum of its two norms,
the usual additional factor at most sqrt(2) must be charged there.

## 2. Explicit radial coefficient jets

Set `d=q_y`, `c=D log M`, and `h=D log alpha=d-beta`. Direct
differentiation, with `r_y=2xb/ell` and `b=rho^-1`, gives

```
d=2kappa^2 r r_y,
beta=kappa^2 r^2 a,
c=kappa^2(1+r)r_y+(3/2)beta.                              (4)
```

The retained `HDBLAST_WEIGHTED_REGULARITY_20261106/SCALAR_REGULARITY.md`,
“Explicit coefficient bounds,” gives the actual polynomial vector field
in variables x,b,a,p,eta, each of absolute value at most one:

```
D x=-10^-10(1-x^2)^2 b,
D b=-q a b,
D a=q(-b^2-p^2/3),
D p=q[U'(eta)-4ap],
D eta=q p.
```

Every vector-field component has polynomial coefficient l1 norm at most
`10^-8` and degree at most nine. For a polynomial of coefficient norm
A and degree d0, repeated differentiation therefore has coefficient norm
at most

```
A product_(i=0)^(k-1) [10^-8(d0+8i)].                    (5)
```

This follows because the sum of coefficient norms of all first partial
derivatives is at most d0 A, and one application raises degree by at
most eight. It controls all centers on the actual annulus, including
arbitrarily close to either endpoint.

The three polynomials (4) have degree at most five and respective
coefficient norms bounded by

```
|c|_coeff <=6.000006*10^-10,
|d|_coeff <=8*10^-10,
|beta|_coeff <=4*10^-16.
```

For example the first c term has norm `6*10^-10`, and the second
has norm `6*10^-16`. In (5), even through k=64 every factor is at
most `10^-8(5+8*63)<1`. In particular the exact needed bounds are

```
max_(k<=16)(|D^k c|,|D^k d|,|D^k beta|) <=1,
max_(k<=16)|D^k h| <=2.                                  (6)
```

Thus no unevaluated coefficient-jet envelope is used in (2). The smooth
background and the authenticated background box are retained inputs, as
in the polynomial certificate; no separate numerical trajectory is used.

## 3. All-angular derivative comparison

Let `nabla_S` be the round covariant derivative acting on all original
and derivative slots, and introduce the normal-ordered intermediate norm

```
K_m(v)^2=sum_(a+b<=m) ||alpha^b nabla_S^b D^a v||_L2(ds dOmega)^2.      (7)
```

The original radial line is trivial for this connection. We prove

```
K_m(v) <= J_m(v) <=2^m K_m(v).                           (8)
```

Here is a rank-independent proof of the required upper angular estimate.
Use the six ambient rotation Killing fields
`K_ij=n_i partial_(n_j)-n_j partial_(n_i)` and their natural Lie
derivatives R_ij. At `n=e4`, the three K_i4 are tangent orthonormal
directions and have zero tangential derivative. The other three Killing
fields vanish at that point and their derivatives are the three elementary
tangent rotation matrices. Consequently, for every geometric tensor A,

```
sum_(i<j)|R_ij A|^2
 =|nabla_S A|^2+sum_(i<j)|spin_ij A|^2.                  (9)
```

There is no cross term at that point, and rotation invariance proves
(9) everywhere. The spin term is nonnegative at every tensor rank.
Natural Lie derivatives commute with the invariant covariant derivative,
including its added derivative indices.

For a scalar, the spin sum is zero, and the natural Casimir is the
positive scalar Laplacian. For a one-form omega, the three nonzero
tangent spin generators give exactly `sum|spin omega|^2=2|omega|^2`.
After integration and polarization, its natural Casimir is

```
-sum R_ij^2=nabla_S* nabla_S+2=Delta_H.                   (10)
```

The last equality is the unit-S3 Weitzenbock identity, whose Ricci
endomorphism on one-forms is `2 Id`. Thus the scalar Casimir eigenvalue
is Lambda=lambda and the one-form eigenvalue is Lambda=mu (lambda in
the exact one-form sector). Isometries commute with these operators, so
every natural rotation word remains in the same eigenspace. Iterating
(9), using derivative-slot naturality, gives

```
||nabla_S^b A||^2
 <=sum_(|I|=b)||R^I A||^2
 =Lambda^b ||A||^2.                                     (11)
```

The last equality follows inductively by summing the first Casimir
identity over all ordered words. It uses no assumption that those
rotations commute. This proof works also for a one-form with sixteen
derivatives, whose final tensor rank is seventeen.

For completeness, (11) implies the requested loose Bochner-type envelope
with an explicit constant: if the rough eigenvalue is nu, then

```
||nabla_S^b A|| <=(nu+2)^(b/2)||A||
 <=[product_(j=1)^b (nu+2j^2)]^(1/2)||A||.                (12)
```

The stronger first expression follows from the exact Casimir calculation,
so no unspecified curvature-commutator constant enters (12).

For the lower bound, let the rough Laplacian eigenvalue be nu. Successive
adjacent-pair traces of `nabla_S^(2k) A` equal `(-nu)^k A`.
For odd order, retain the first derivative and trace the remaining
adjacent pairs to obtain `(-nu)^k nabla_S A`. These statements use
metric-parallel contractions in the indicated derivative order; no two
covariant derivatives are commuted. Each trace has Hilbert operator
norm sqrt(3). Since `||nabla_S A||^2=nu||A||^2`,

```
||nabla_S^b A||
 >=nu^(b/2) 3^(-floor(b/2)/2)||A||.                      (13)
```

For scalars, nu=Lambda. For the exact unit one-forms,
`nu=lambda-2>=3lambda/4` because lambda>=8. For the retained
coexact unit one-forms, `nu=mu-2>=3mu/4` because mu>=9.
Hence in every nonconstant block

```
2^-b Lambda^(b/2)||A|| <=||nabla_S^b A||
                         <=Lambda^(b/2)||A||.            (14)
```

Indeed `(3/4)^(b/2)3^(-floor(b/2)/2)>=2^-b`.
For the scalar and radial one-form monopoles all reference angular
derivatives vanish, in agreement with zeta=0.

The radial and tangential original one-form slots are orthogonal and
remain so under nabla_S. Within any eigenspace the preceding inequalities
apply to an arbitrary linear combination, not only to an individual basis
vector. Distinct natural Casimir eigenspaces remain orthogonal after any
number of covariant derivatives, since naturality intertwines the
self-adjoint Casimir on the derivative tensor bundles. This handles
degeneracies and both coexact families without adding a mode-count factor.
Multiply (14) by alpha^b, apply it to D^a v, integrate in s and sum.
This proves (8) for arbitrary finite harmonic sums, hence by completion
for all modes.

## 4. Exact scaled physical-jet recurrences

Let the original rank be rT=0 or 1 and set

```
S_j=M q^j I(nabla_gamma0^j T),  S_0=v.
```

Their L2 square sum is exactly W_m(T)^2. In the radially parallel
physical frame, the difference between the physical connection and the
reference splitting connection is `a K_A` in each angular direction.
For one covariant slot, K_A is an elementary skew rotation mixing the
radial slot and the A angular slot, so its operator norm is one. On a
rank-R tensor it acts on each slot and has norm at most R. Collecting
the three angular directions gives

```
||K_(R)||op <=sqrt(3)R <=2R.                             (15)
```

This charges the j covariant derivative slots as well as the original
slots: here R=rT+j. These slot maps, radial/angular projections and
permutations are round-parallel tensor maps; this assertion uses the
induced round connection on their direction indices as well.

Splitting the new derivative index into its radial and angular parts
gives the exact recurrences

```
(S_(j+1))_0=(D-c-jd) S_j,
(S_(j+1))_S=alpha nabla_S S_j+beta K_(rT+j) S_j.           (16)
```

The term jd differentiates q^j. The term c differentiates the unitary
weight M. The second row includes the original and derivative-slot
connection action. Conversely,

```
D S_j=(S_(j+1))_0+(c+jd)S_j,
alpha nabla_S S_j=(S_(j+1))_S-beta K_(rT+j)S_j.            (17)
```

The new-index projections in these identities have norm at most one.
Since `(D alpha)/alpha=h`, the normal-ordered tensors

```
U_(a,b)=alpha^b nabla_S^b D^a v
```

obey

```
D U_(a,b)=U_(a+1,b)+b h U_(a,b),
alpha nabla_S U_(a,b)=U_(a,b+1).                         (18)
```

These identities explicitly charge every interleaving of radial and
angular differentiation. In particular `[D,zeta]=h zeta`, with the
relative coefficient bounded by (6); commuting them as constants is
not permitted.

## 5. Finite coefficient accounting in both directions

Here is an explicit induction bound, rather than an unexplained large
ceiling. A coefficient term is a product of derivatives of c,d,beta,h
and a fixed round-parallel slot map. Give each c,d,beta factor cost one,
each h factor cost two, and multiply by the map operator norm and the
absolute numerical coefficient. The total cost is the sum over terms.
Equation (6) bounds the actual multiplier norm by that formal cost.
At stage j every monomial has at most j scalar factors. Differentiating
its coefficient raises total cost by at most a factor j: each factor is
replaced by the next derivative, with the same assigned bound. Required
derivatives have order at most sixteen.

For the forward expansion of S_j in the U_(a,b), take an expansion with
total cost C_j. A radial step in (16) has the following costs:

* coefficient differentiation: at most j;
* the direct `U_(a+1,b)` term: one;
* the `b h` term in (18): at most 2j;
* multiplication by `c+jd`: at most 1+j.

Thus the radial contribution costs at most `4j+2`. The angular direct
term costs one and its spin term costs at most `2(rT+j)` by (15).
The full next tensor is the orthogonal radial/angular sum; bounding its
norm by the sum of these two contributions is valid. Therefore

```
C_0=1,
C_(j+1) <=(6j+2rT+3) C_j <=6(j+1) C_j,
C_j <=6^j j!,                         rT<=1.             (19)
```

Every U in this expansion has a+b<=j. It follows that
`||S_j||<=6^j j! K_m` for j<=m, and hence

```
W_m <=sqrt(m+1)6^m m! K_m.                              (20)
```

For the reverse expansion, express each U_(a,b) with a+b=n in S_j,
j<=n, using (17) and (18). If a radial step is used, coefficient
differentiation costs at most n, the direct new S row costs one,
`c+jd` costs at most 1+n, and subtracting `b h` costs at most 2n.
The total is at most `4n+2`. An angular step costs at most
`1+2(rT+n)<=2n+3`. Either is at most `4(n+1)` for rT<=1.
Thus every U_(a,b) of order n has an expansion of cost at most
`4^n n!`, and

```
||U_(a,b)|| <=4^n n! W_m,                  n<=m.
```

There are exactly `(m+1)(m+2)/2` pairs a,b with a+b<=m. Consequently

```
K_m <=sqrt((m+1)(m+2)/2)4^m m! W_m.                     (21)
```

Equations (8), (20) and (21) prove (2). For m<=16 the pair count is at
most 153<13^2. Its reverse constant dominates the forward one, so the
integer in (3) follows. The factor counts include the connection action
on all derivative indices; there is no independent per-component or
per-harmonic multiplicity reserve hidden in the result.

## 6. Source normalization, zero-order identity and completed domains

At order zero there is exact unitarity:

```
||T||_L2(w^2 dvol_gamma0)^2
 =integral |M I(T)|^2 ds dOmega
 =J_0(v)^2.                                              (22)
```

For the normalized equation `L U=f`, with f in the declared weighted
range Y, its transformed source is exactly `g=M I(f)`. Applying (2)
at m=12 gives the same two-sided source-norm bridge. If instead a symbol
`f_phys` denotes the unnormalized physical constraint source before
the w^-2 factor, then

```
f=w^-2 f_phys,
g=M I(f)=sqrt(q rho^3) w^-1 I(f_phys).                    (23)
```

This distinction agrees with the retained scalar regularity source
definition. No w^-2 factor may be dropped when passing between these
two source conventions. For the momentum block the input one-form is
the source after division by the **fixed** canonical density rho0;
using a variable density changes the operator and is outside this lemma.

The argument first applies to compact smooth radial tests and finite
smooth harmonic sums. Such tests are dense in the physical dotted
weighted spaces: on compact subannuli the weight is positive and smooth,
and angular spectral truncation converges in every fixed finite Sobolev
norm; then use the defining compact-support closure. Conversely, finite
J_m-norm fields are approximated by angular truncation, compact s cutoffs
and smoothing. After angular truncation zeta is a smooth bounded radial
coefficient; cutoff error terms tend to zero by the finite derivative
norms and tail integrability. On compact s intervals the remaining
smoothing convergence is ordinary finite-order Sobolev convergence.

The two estimates therefore extend to the respective completions and
identify them continuously in both directions. They impose no extra
endpoint trace on a potential. The physical correction's later
zero-extension theorem still uses its own fixed multiplication/operator
map; the present lemma neither changes that map nor replaces its domain.
All finite sums, truncations and closures can be restricted to the
reflection-invariant subspace because the radial weights, connections
and angular spectral projections commute with the reflections.

The graph expression `A=1-D^2+zeta^2` is compatible with these variables.
This proof supplies the full physical-norm/J_m bridge and the required
relative commutator bound; a separate graph regularity/inverse argument
must still prove the desired uniform estimate for the actual transformed
scalar and vector operators. No inverse constant is inferred from norm
equivalence alone.
