# All-mode physical correction maps for the selected weights

13 November 2026 UTC. Bounded linear-map derivation. This is a bounded linear-map
estimate for arbitrary potentials in the actual dotted weighted spaces.
It does not assign a nonzero amplitude, prove a nonlinear contraction,
complete the harmonic Cauchy map or control spacetime evolution.

## 1. Exact maps, norms and evaluated result

Use the selected annulus and the background conventions of the accepted
`HDBLAST_ALL_MODE_INVERSE_20261113/WEIGHTED_NORM_BRIDGE.md`:

```
gamma0=dy^2+rho^2 kappa_S3, dy=-rho dz,
ell=10^-6, x=(z-3ell/2)/ell, r=1-x^2,
kappa0=10^-8, q=kappa0^2 r^2, w=exp(-1/r),
a=rho_y/rho, b=rho^-1, p=phi0_y,
ds=dy/q, D=partial_s, M0=sqrt(q rho^3), M=w M0.
```

All scalar, one-form and tensor norms use gamma0, full sphere measure,
and the full covariant derivative tensor, including all derivative slots.
The weighted norm is

```
W_m(P)^2=sum_(k=0)^m ||w q^k nabla^k P||_L2(gamma0)^2.
```

The dotted domain is the compact-smooth closure. No angular-mode cutoff
or restriction to the original second source is imposed in the proof.
The result holds before restriction to the reflection-invariant subspace.

The exact corrections are

```
h=A_H N=w^2 q^4 P_g* N,
tau=A_M Z=w^2 q^2 P_m* Z,
pi_correction=rho0 tau,  rho0=dvol_gamma0.
```

Here tau is an ordinary symmetric contravariant tensor. For estimates we
lower both indices using gamma0, which commutes with every covariant
derivative and preserves its norm. The fixed canonical density rho0 is
parallel for the induced background density connection. Thus the natural
background norm of pi_correction/rho0 equals the ordinary norm of tau;
coordinate partial derivatives of a density are not substituted for that
convention.

Define the literal high-factor sums

```
E0(h)=sum_(j=0)^14 ||w^-1 q^(j-2) nabla^j h||_2,
E1(tau)=sum_(j=0)^13 ||w^-1 q^(j-1) nabla^j tau||_2.
```

The evaluated constants are

```
C_H=8*17^14*14! =117431129563348386287300198400 <10^30,
C_M=8*17^13*13! =493408107409026833139916800 <10^27,

E0(A_H N) <=C_H W_16(N),
E1(A_M Z) <=C_M W_14(Z).                                  (1)
```

These constants also bound the corresponding square-sum norms, since
each is at most its literal sum. Consequently the physical maps satisfy

```
||A_H N||_H14(gamma0) <=10^256 W_16(N),
||A_M Z||_H13(gamma0) <=10^253 W_14(Z).                    (2)
```

The proof below charges all multiplier, adjoint and background derivatives
through the required orders. The factor `sup w q^-12<10^226` is used
only after those complete weighted output-jet estimates have been proved.

## 2. The actual adjoints and their coefficient jets

The retained canonical-constraint source fixes

```
P_g* N=(1/2)[Hess N-gamma0 Delta N-N Ric(gamma0)
                                      +N dphi0 tensor dphi0],
P_m* Z=Lie_(Zsharp)gamma0
```

after lowering the second expression's two output indices. The scalar
principal map on a two-tensor is

```
A(T)=(1/2)[T-gamma0 tr(T)],   ||A||op=3/2.                (3)
```

It is parallel. Its trace eigenvalue is -3/2 and its tracefree
eigenvalue is 1/2. For momentum the map `B_sym(T)=T+transpose(T)` is
parallel and has operator norm two, in the full Frobenius convention.

The background ODE `a_y=-b^2-p^2/3` and the warped metric give

```
Ric_yy=3b^2+p^2-3a^2,
Ric_tan=(3b^2+p^2/3-3a^2) Id_tan.
```

Hence the scalar zero-order coefficient is

```
B=(1/2)[-Ric+dphi0 tensor dphi0]
 =B_R dy^2+B_T gamma_tan,
B_R=(3/2)(a^2-b^2),
B_T=(3/2)(a^2-b^2)-p^2/6.                                (4)
```

The radial and tangential projectors are round-parallel in the reference
splitting. Their physical angular derivatives are accounted for by the
connection on all output slots below.

Put `d=q_y`, `e_w=D log w`, `beta=qa`, and `c=D log M`. The exact
polynomials are

```
r_y=2xb/ell,
e_w=kappa0^2 r_y,
d=2kappa0^2 r r_y,
beta=kappa0^2 r^2 a,
c=e_w+(1/2)d+(3/2)beta.                                  (5)
```

The accepted norm bridge and its retained scalar regularity source give
the polynomial D-vector field in x,b,a,p,eta. Each component has
coefficient l1 norm at most 10^-8 and degree at most nine. Therefore a
polynomial of coefficient norm A0 and degree d0 has kth derivative
coefficient norm at most

```
A0 product_(i=0)^(k-1) [10^-8(d0+8i)].                    (6)
```

The actual variables have absolute value at most one on the whole
annulus. In particular all D jets through order fourteen of
`e_w,d,beta,c` have absolute value at most one. This follows already
from their coefficient norms `2*10^-10`, `8*10^-10`, `4*10^-16`,
`6.000006*10^-10`, respectively, and degrees at most five.

The coefficient needed from (4) is **q^2 B**, not a bare B estimate.
Since `q^2=kappa0^4(1-x^2)^4`, its scalar coefficient degree is at most
ten and

```
||q^2 B_R||coeff <=48*10^-32 <10^-30,
||q^2 B_T||coeff <=(152/3)*10^-32 <10^-30.                 (7)
```

Every factor in (6) for d0<=10 through k=14 is less than one. Thus
every D derivative through fourteen of each scalar coefficient in (7)
is also below 10^-30. This charges all radial background and multiplier
derivatives required by the map. The angular connection conversion in
the next section charges the remaining coefficient-tensor derivatives.

## 3. Complete scaled covariant jet recurrences

Use the radially parallel fiber isometry I on every original and
derivative slot, as in the accepted norm bridge. In this section S_k
denotes input jets, with original input rank rP=0 for N and rP=1 for Z:

```
S_k=M q^k I(nabla^k P),
sum_(k<=m) ||S_k||_L2(ds dOmega)^2 =W_m(P)^2.
```

Let `alpha=q/rho`. The exact inverse input recurrences are

```
D S_k=(S_(k+1))_0+(c+kd)S_k,
alpha nabla_S S_k=(S_(k+1))_S-beta K_(rP+k)S_k.            (8)
```

The subscript on a new jet selects its new radial or angular derivative
index and has operator norm at most one. The angular connection map
acts on **every** tensor slot. For rank R, one angular direction is
the sum of R elementary radial/angular skew rotations, each of norm one.
Collecting three directions gives

```
||K_R||op <=sqrt(3)R <=2R.                                (9)
```

All slot maps, radial/tangential projections and permutations used below
are round-parallel when their angular direction indices carry their
induced connection.

Treat the two output maps together by taking

```
sigma=2, p_order=2, rP=0, T=h,   or
sigma=1, p_order=1, rP=1, T=tau.
```

In both cases the output has original rank two and
`rP+p_order=2`. Set

```
Y_j=M0 w^-1 q^(j-sigma) I(nabla^j T).                     (10)
```

Because M0 is the bare radial volume factor,

```
||Y_j||_L2(ds dOmega)
 =||w^-1 q^(j-sigma) nabla^j T||_L2(gamma0).              (11)
```

Their initial expressions are exact:

```
Y_0=A S_2+(q^2 B)S_0                     for h,
Y_0=B_sym S_1                            for tau.         (12)
```

The first initial expression has total coefficient-map norm at most
`3/2+3*10^-30<2`, by (3), (7), and
`|gamma_tan|=sqrt(3)<2`. The second has norm at most two.
No quadratic or nonlinear adjoint term is present in these fixed maps.

Differentiating (10), including its radial prefactor and every existing
tensor slot, gives

```
(Y_(j+1))_0=[D-c+2e_w+(sigma-j)d]Y_j,
(Y_(j+1))_S=alpha nabla_S Y_j+beta K_(j+2)Y_j.             (13)
```

In particular the derivative of the flat multiplier is not discarded.
When the radial operator in (13) acts on a coefficient times S_k,
the c from (8) cancels, leaving exactly

```
coefficient_derivative*S_k + coefficient*(S_(k+1))_0
       +[2e_w+(k+sigma-j)d]*coefficient*S_k.              (14)
```

The angular row similarly has its direct S_(k+1) part and the two
connection actions, one on the output and one from (8) on the input.
No radial/angular derivative is reordered or replaced by a scalar
component derivative.

## 4. Explicit finite ledger through all required orders

Expand Y_j in the S_k, k<=j+p_order. A term has a fixed round-parallel
slot map multiplied by a monomial in D derivatives of e_w,d,beta and
the two small coefficients (7). Assign each primitive derivative its
bound from (5)–(7), multiply by the map operator norm and the absolute
numerical coefficient, and sum these costs.

At stage j every coefficient monomial has at most j+1 scalar factors:
the scalar lower term starts with one; every subsequent multiplication
adds at most one. Differentiating a coefficient therefore costs at most
j+1 times its assigned bound. Replacing a primitive by its next D
derivative preserves that bound. Angular differentiation of these scalar
coefficients vanishes.

The radial contribution (14) costs at most

```
(j+1) +1 +2+|k+sigma-j| <=2j+8,                          (15)
```

since `0<=k<=j+p_order` and `sigma=p_order<=2` imply
`|k+sigma-j|<=j+4`. For the angular connection terms the input rank
is `rP+k<=j+2`, and the output rank is j+2. Equations (8)–(9) thus
give their combined cost at most `4(j+2)`. Including the direct
new angular input derivative, the angular contribution costs at most
`4j+9`. The norm of the full radial/angular tensor is at most the sum
of its two pieces. If D_j bounds the total coefficient cost, then

```
D_0<=2,
D_(j+1)<=(6j+17)D_j<=17(j+1)D_j,
D_j<=2*17^j*j!.                                         (16)
```

This ledger includes all Leibniz branches through the recursion, the
adjoint maps, multiplier derivatives, derivative-index connections,
background jets, and input/output slot conversion. At no point is a
fixed harmonic polynomial used to bound an arbitrary field.

With m=14 for scalar and m=13 for momentum, every input jet in the
expansion is among those in W_(m+p_order). Hence

```
||Y_j||_2 <=2*17^j*j! W_(m+p_order)(P),    0<=j<=m.        (17)
```

The descending ratio between successive terms in the coefficient sum is
at most 1/17. Therefore, also for the literal sum convention,

```
sum_(j=0)^m 2*17^j*j!
 <=2*17^m*m! sum_(r=0)^m 17^-r
 <(17/8)*17^m*m! <8*17^m*m!.                             (18)
```

This proves (1) and the stated two exact integers. It avoids adding a
second derivative-count factor when a later argument uses E0/E1 as
literal sums instead of square sums.

## 5. Physical norms and zero extension

Only now strip the weighted factors. Since 0<q<=10^-16<1, for every
scalar j<=14 and momentum j<=13 the residual multiplier satisfies

```
w q^(sigma-j) <=w q^-12.
```

Writing t=1/r>=1 gives

```
sup_A w q^-12
 =10^192 sup_(t>=1)t^24 exp(-t)
 <=10^192*24^24 <10^226.                                (19)
```

For each j multiply (11) by this supremum, then sum derivative squares.
Its square sum is bounded by the literal sum already controlled in (1).
Using `C_H<10^30` and `C_M<10^27` proves (2). There is exactly one
use of the flat-factor bound (19); none of the differentiated
coefficients was treated as an undifferentiated flat multiplier.

For compact smooth potentials, h and tau are compact smooth tensors in
the same annulus. Approximate any dotted potential by such tests. The
bounded maps (2) give convergence in physical H14/H13. Extending each
compact correction by zero produces convergence in the ambient initial
slice norm with exactly the same integral norms. Thus

```
A_H:dot H16_(q,w) -> H14_0(A),
A_M:dot H14_(q,w) -> H13_0(A).                            (20)
```

Their zero extensions have no distributional boundary layer. No
unweighted trace of N or Z is assumed. The fixed open shell and vertex
collars remain unperturbed because the closed annulus avoids both.
The fixed-density conversion to pi commutes with this construction.

## 6. A complete, deliberately loose relative-metric C0 guard

The separate weighted point-embedding lemma may give a much smaller
guard. A finite guard is already available from retained generic tensor
estimates, without any angular restriction. After the zero extension in
(20), view h on the original full collar `[0,6ell] x S3`.
`HDBLAST_INITIAL_TRANSPORT_20261106/INITIAL_TRANSPORT.md` proves for
the generic physical-to-product fiber isometry

```
||I(h)||_H3(dz^2+kappa_S3) <=10^3000 ||h||_H3(gamma0).
```

`HDBLAST_ORDERED_EVOLUTION_20261106/QUADRATIC_FORCING.md`, section 4,
proves the generic boundary-inclusive product-cylinder H3 embedding
with constant 10^7. The original tensor rank here is two, within its
stated budget. Since I is a pointwise isometry,

```
||gamma0^-1 A_H N||_op,infinity
 <=||A_H N||_gamma0,infinity
 <=10^3007 ||A_H N||_H14(gamma0)
 <=10^3263 W_16(N).                                     (21)
```

Thus one may take `C_A0=10^3263` as an explicit conservative fallback.
This is a constant assembled from two cited all-tensor estimates and
(2), not a new embedding inferred from the special finite source.

If `S_0` is an independently evaluated relative-metric C0 bound for the
actual affine seed, then the entire product ball is positive whenever

```
|e| S_0+r C_A0 <=1/4.                                   (22)
```

It yields eigenvalues between 3/4 and 5/4 relative to gamma0. This is
only a sufficient design guard; no amplitude or nonlinear solution is
assigned here. A tighter separately proved C_A0 can replace (21).

For completeness, the background metric has physical H14 norm below
220 on this annulus: its higher covariant derivatives vanish, its norm
is two in dimension four, and
`Vol(A)<40*131^4/10^6<12100`. Consequently, if S_H and S_M denote
the actual affine seed's physical H14 metric and normalized H13 momentum
norms, a product ball has the finite estimates

```
||gamma||_H14(gamma0) <220+|e|S_H+10^256 r,
||pi/rho0||_H13(gamma0) <=|e|S_M+10^253 r.               (23)
```

The seed constants and the complete uniform second derivative of the
nonlinear constraint map are separate inputs. The correction maps do
not by themselves certify the nonlinear ball, a parameter interval,
harmonic gauge completion, spacetime PDE control, or a physical origin.
