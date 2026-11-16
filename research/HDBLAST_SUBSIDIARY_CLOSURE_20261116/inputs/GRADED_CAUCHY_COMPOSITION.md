# Uniform graded canonical and harmonic Cauchy compositions

13 November 2026 UTC. Graded composition proof. All statements retain the
actual selected background, canonical normalization, lapse and relative-
harmonic gauge. They hold on an explicit whole physical data ball, not
only at amplitude zero. This note supplies compositions and norm bridges;
the later boundary-chart map and evolved solution tube are separate.

## 1. Fixed norms and the whole data ball

Use the accepted annulus, weights and fixed background

```
gamma0=dy^2+rho^2 kappa_S3, dy=-rho dz,
129<=rho<=131, a=rho_y/rho,
ell=10^-6, kappa_w=10^-8,
r=(z-ell/2)(5ell/2-z)/ell^2,
q=kappa_w^2 r^2, w=exp(-1/r).
```

All spatial derivatives, tensor slots and volumes use gamma0. Retain the
literal-sum physical graded norms from PR140:

```
E_d(T)=sum_(j=0)^(14-d)
           ||w^-1 q^(j+d-2) nabla0^j T||_L2(gamma0),
d=0,1,2.
```

Their accepted interfaces include graded products with constant 10^17,
`E_(d+1)(nabla0 T)<=E_d(T)`, grade inclusion
`E_(d+1)(T)<=E_d(T)`, and pointwise control
`|T|<=10^12 w E0(T)` for grade-zero tensors. Tensor products in this
application have input ranks below 25, including all subsequent fixed
contractions; the accepted rank restrictions are respected.

Write the five ordinary fixed-density physical inputs as

```
z=(h,psi,tau,omega,l),
gamma=gamma0+h, phi=phi0+psi,
pi_act=rho0 tau, varpi_act=rho0 omega,
rho0=sqrt(det gamma0), alpha=rho(1+l),
delta(z)=E0(h)+E0(psi)+E1(tau)+E1(omega)+E0(l).
```

The density rho0 is fixed. It is not sqrt(det gamma). The lapse alpha
is the registered affine lapse when l=e f; it is not the radial metric
root used by the later boundary chart. Put Q=10^100 and use the open ball

```
Bphys={z: Q delta(z)<1/1024}.                            (1)
```

For an arbitrary direction z_i write D_i=delta(z_i). All first and
second derivatives below are Fréchet derivatives in this five-field
space. The bounds are uniform throughout (1), including arbitrary
angular modes and arbitrary state/direction derivative allocations.

The registered lapse profile satisfies an exact useful identity:

```
h_seed=f(-4dy^2+2gamma_tan),
f=(1/2)tr_gamma0 h_seed,
E0(f)<=E0(h_seed)<=10^200.                              (2)
```

The trace norm is two in dimension four and commutes with every
covariant derivative. Thus l=e f, l'=f, l''=0 fit the same graded
description without a new harmonic-polynomial estimate. The full
metric temporal component still contains l^2; it is not affine.

## 2. Evaluated background multipliers involving rho

The accepted nonlinear constraint ledger covered a,b=1/rho, the scalar
background and Ric0. Positive powers of rho must be charged separately.
Let D_s=q D_y, d_q=q_y and beta=qa. The retained polynomial background
certificate gives, through all derivative orders needed here,

```
|D_s^j a|, |D_s^j d_q|, |D_s^j beta|<=1,
D_s log rho=beta.
```

For an integer k with |k|<=3 set
`R_j=q^j rho^-k nabla0^j(rho^k)`.
The radial scalar covariant-jet recursion is

```
(R_(j+1))_0=D_s R_j+(k beta-j d_q)R_j,
(R_(j+1))_tan=beta K_j R_j.                            (3)
```

All derivative slots are included in the rank-j rotation K_j; its
aggregate norm is at most 3j. These radial tensors are built from dy
and the angular metric, so their pure round-sphere rotation action
vanishes. A coefficient monomial at stage j contains at most j profile
factors. Differentiating those factors costs at most j, multiplication
in the first row costs at most |k|+j, and the angular action costs 3j.
The total next-step coefficient factor is at most
`|k|+5j<=5(j+1)`. Therefore

```
sup |q^j nabla0^j(rho^k)| <=131^3*5^j*j!, j<=14.        (4)
```

For k<0, rho^k<=1, so the same ceiling is valid. This includes all
needed inverse powers of rho as well as rho,rho^2,rho^3.

The remaining nonparallel multiplier in the harmonic mixed velocity is

```
B_rho=rho d rho=rho^2 a dy.
```

Apply the same recursion to `q^j rho^-2 nabla0^j B_rho`, whose initial
tensor is a dy and whose rank is j+1. The coefficient derivative costs
at most j+1, the radial multiplication costs 2+j, and the angular
connection costs 3(j+1). Thus the next-step factor is at most
`5j+6<=6(j+1)`. It follows that

```
sup |q^j nabla0^j B_rho|<=131^2*6^j*j!, j<=14.          (5)
```

The exact endpoint arithmetic is

```
131^3*5^14*14! =1196195873059687500000000000 <10^28,
131^2*6^14*14! =117238012873166550151987200 <10^28.
```

The sequences are increasing, so these bounds cover every lower order.
The accepted weighted Leibniz multiplier estimate now gives a common
module norm on every required E_d:

```
E_d(B tensor T)<=10^33 E_d(T),
B in {rho^k: |k|<=3} or B=B_rho.                       (6)
```

Indeed the finite exterior derivative sum is at most 2^15-1<10^5.
A fixed trace or index contraction is separately charged below. The
new multiplier bound 10^33 is below Q=10^100; no enlarged algebra
smallness assumption has been imposed. Parallel gamma0 and the fixed
density rho0 retain their exact parallel action. No coordinate density
derivative or varying-density factor enters (6).

## 3. The positive fractional determinant and matrix inverse

Raise one index of h with gamma0 and denote the resulting matrix by H.
Let

```
B=I+H, A=B^-1,
J=det(gamma)/det(gamma0)=det(I+H),
S=J^-1/2.
```

The pointwise bound and (1) make B positive everywhere. Hence S is the
positive real square root of the inverse determinant, with its branch
fixed on the entire ball.

The accepted unital E0 matrix algebra has product constant at most Q.
Its Neumann and logarithm series give, with L=log det(I+H),

```
E0(A-I)<=2 E0(h),
E0(DA[h_i])<=4 E0(h_i),
E0(D^2A[h_1,h_2])<=16Q E0(h_1)E0(h_2),

E0(L)<=4 E0(h),
E0(DL[h_i])<=4 E0(h_i),
E0(D^2L[h_1,h_2])<=8Q E0(h_1)E0(h_2).                 (7)
```

These are full matrix statements, not a scalar-conformal restriction.
The constant identity is a unital multiplier, not an E0 field.

The required new fractional composition is `S=exp(-L/2)=1+s`.
The exponential series and `Q E0(L)<=4/1024` give

```
E0(s)<=4 E0(h),
E0(DS[h_i])<=4 E0(h_i),
E0(D^2S[h_1,h_2])<=16Q E0(h_1)E0(h_2).                 (8)
```

For the last two lines use the exact identities

```
DS=-(1/2)S DL,
D^2S=S[(1/4)DL[h_1]DL[h_2]-(1/2)D^2L[h_1,h_2]].
```

The multiplier 1+s has norm below two. Thus the second bound costs
`(1/2)*2*4=4`, and the third costs
`2[(1/4)*4*4+(1/2)*8]=16`, with one Q for the bilinear product.
The derivative series converge uniformly on closed subballs of (1).
The graded algebra already includes spatial derivatives falling on
arbitrary state factors throughout these series. No unevaluated
spatial Faà di Bruno constant is invoked.

For later bookkeeping put e_d(T)=Q E_d(T). In these rescaled ideal
norms, graded multiplication has norm at most one. The unital factors
B,A,1+l,S have module norms at most two. Their first derivatives have
e0 norm at most `4 sigma_i`, and their second derivatives at most
`16 sigma_1 sigma_2`, where sigma_i=Q D_i. The first and second
derivatives of the affine factors B,1+l themselves are smaller.

## 4. Exact canonical conversion and background-subtracted Cauchy map

The actual half-action canonical normalization is

```
pi_act=sqrt(det gamma)(K^ij-gamma^ij tr_gamma K)/2,
varpi_act=sqrt(det gamma) Pi.
```

Consequently the exact physical fields are

```
T=gamma tau gamma-(tr_gamma tau)gamma/3,
K=2S T, Pi=S omega.                                    (9)
```

The two gamma factors act by the matrix B in the fixed gamma0
identification. The trace is `tr_gamma tau=gamma_ij tau^ij`.
The background momenta vanish, but no canonical-to-physical nonlinear
term is discarded. In particular affine omega does not make Pi affine.

For the retained zero-shift initial section and relative-harmonic
completion, the exact coordinate Cauchy identities are

```
g_tt=-rho^2(1+l)^2, g_ti=0, g_ij=gamma_ij,
phi=phi0+psi,

g_ij,t=4rho(1+l)S T_ij,
g_tt,t=(4/3)rho^3(1+l)^3 S tr_gamma tau,
phi_t=rho(1+l)S omega,

g_ti,t=rho^2(1+l)^2 C_i
          +gamma_ij rho (nabla0^j rho)-alpha partial_i alpha,
C_i=(gamma^-1)^jk[nabla0_j h_ki-(1/2)nabla0_i h_jk].    (10)
```

The companion completion proof derives (10) from the actual reference
connection and verifies the cancellation of its spatial Christoffel
terms. Here (10) is used as that exact formula interface.

Before estimating, subtract the background from every nonzero position
or velocity. The mixed velocity has the exact graded expression

```
g_ti,t=rho^2(1+l)^2 C_i
       +(B_rho)^j h_ij
       -(2l+l^2)(B_rho)_i
       -rho^2(1+l)nabla0_i l.                          (11)
```

Thus the nondecaying background terms cancel. Positions are represented
by `g_tt+rho^2=-rho^2(2l+l^2)`, h, and psi. All velocities in
(10)--(11) have zero background value. Every output term now contains
an E0 or E1 physical perturbation. No E norm is assigned to the full
background metric, lapse or scalar.

## 5. Uniform state, direction, and second-direction ledger

Let the graded output size be the literal sum

```
Pout=E0(g_tt+rho^2)+E0(h)+E0(psi),
Vout=E1(g_tt,t)+2E1(g_ti,t)+E1(g_ij,t)+E1(phi_t).
```

The factor two in the mixed vector block covers its two symmetric
spacetime metric slots. Let W(z) denote this background-subtracted
unstraightened Cauchy map. One may append K and Pi to the output
without changing the constants below.

Expand the finite formulas (9)--(11) in products of the unital factors
B,A,1+l,S, and the ideal factors
`h,psi,tau,omega,l,nabla0 h,nabla0 l`.
There are at most 32 monomials, including the optional K,Pi output.
Every monomial has at most six nonbackground factors. For example the
largest is `(1+l)^3 S tr(B tau)` in g_tt,t. Every monomial contains
at least one ideal factor, whose state norm is at most
s0=Q delta<1/1024. Each term has exactly one of the background
multipliers in (6), or the identity.

The numerical coefficient, fixed contraction and symmetric mixed-block
factor together are bounded by 16 per monomial. The largest displayed
numerical coefficient is four, a trace costs at most two, and the mixed
block costs two. The two summands defining C_i are treated separately,
so their coefficient does not require an additional hidden trace or
permutation count. The free tensor in `tr(B tau) B` also needs an
identity insertion: write it as `B[(tr(B tau)) I]`. Its Frobenius
insertion norm is two, separately from B's single-index module norm.
Together with the trace norm two and the coefficient 4/3 in g_ij,t,
the numerical/contraction cost is 16/3<16; this term has no mixed-block
factor. The corresponding K term costs only 8/3. All other index
permutations are isometries, and each monomial remains within the
common ceiling 16.

For a product of n<=6 nonbackground factors, use their module bound two,
first-direction bound 4 sigma_i, and second-direction bound
16 sigma_1 sigma_2. The ordinary product rule gives the following
uniform ceilings in the rescaled output ideal norm:

```
state:       2^(n-1) s0 <=32 s0,
first:       n*4*2^(n-1) sigma_i <=768 sigma_i,
second:     [n*16*2^(n-1)+n(n-1)*16*2^(n-2)]
                         sigma_1 sigma_2
                                      <=10752 sigma_1 sigma_2. (12)
```

For the second line of variations, ordered choices of two different
factors account for both parameter allocations. A second derivative
falling on A or S is charged by (7)--(8). State factors remaining
elsewhere in a differentiated monomial are retained, with their norm
bounded by two; they have not been replaced by zero. Grade-one factors
and their directions remain in E1, while the unital compositions are
grade zero. At most one total grade-one factor occurs in any velocity
monomial, so the graded product rule applies directly.

Combining (6), (12), 32 monomials and the coefficient ceiling 16, and
then removing the rescaling Q, yields

```
Pout+Vout <=16384*10^33 delta <=10^38 delta,
||DW[z_i]||graded <=393216*10^33 D_i <=10^39 D_i,
||D^2W[z_1,z_2]||graded
 <=5505024*10^33*Q D_1D_2 <=10^140 D_1D_2.             (13)
```

Here the output graded norm is the displayed literal block sum. The
same proof bounds the optional physical K,Pi outputs within the stated
32-monomial reserve. Convenient non-strict common roundings, also valid
at zero, are

```
||W(z)||graded <=10^40 delta(z),
||DW(z)|| <=10^40,
||D^2W(z)|| <=10^140,       z in Bphys.                 (14)
```

Thus a single looser C2 composition ceiling 10^200 is available if a
later interface prefers one number. Its larger exponent is only a
rounding; it is not used to redefine the smallness condition (1).

## 6. Physical H14/H13 conversion without a new Sobolev constant

For any tensor of grade d=0 or 1 and j<=14-d, factor its unweighted
physical derivative exactly as

```
nabla0^j T=[w q^(2-j-d)]
                 [w^-1 q^(j+d-2)nabla0^j T].
```

The worst exponent is -12 in both grades. The accepted flat multiplier
calculation gives

```
sup w q^-12
 =10^192 sup_(t>=1) t^24 exp(-t)
 <=10^192*24^24<10^226.
```

The same ceiling controls smaller exponents because q<=1 and w<=1.
The physical derivative-square norm is bounded by the literal sum,
so there is no extra derivative-count factor:

```
||T||_H^(14-d)(gamma0)<=10^226 E_d(T), d=0,1.           (15)
```

Apply (15) to the Cauchy blocks. Their literal block sum dominates the
full raw tensor/scalar Hilbert norm in the fixed positive fiber metric
`dt^2+gamma0`, including the mixed metric slots. Therefore, with

```
R(W)=||position||_H14 +||velocity||_H13,
```

the actual unstraightened map obeys on the entire ball (1)

```
R(W(z))<=10^266 delta(z),
R(DW[z_i])<=10^266 D_i,
R(D^2W[z_1,z_2])<=10^366 D_1D_2.                       (16)
```

This conversion is a scalar multiplier bound. No unweighted Sobolev
embedding, unstated chart-cover constant, or finite angular polynomial
estimate has entered it.

For a C2 curve with graded envelopes

```
sup delta(z(e))<=10^-6201,
sup delta(z'(e))<=10^201,
sup delta(z''(e))<=10^6432,
```

the usual first and second chain rules, with every term retained, give

```
R0<=10^-5935,
R1<=10^467,
R2<=10^6698+10^768<=10^6699<=10^6700.                  (17)
```

The companion selected-curve proof supplies those three input envelopes
on the already admitted interval. Equation (17) is stated conditionally
on those precise interfaces and has ample strict inclusion in (1).
The second term in R2 is the composition Hessian applied to two first
directions; it has not been dropped because the scalar sectors are affine.

## 7. Domain and scope

The polynomial, Neumann and exponential constructions are C2 on the
whole explicit open physical ball (1). Products are interpreted in the
accepted completed graded spaces; background multipliers and parallel
identities remain in their separate multiplier algebra. Compact smooth
approximants establish the physical zero-extension statements through
(15), while the uniform algebra estimates pass derivatives to the
limits. The fixed canonical density rho0 stays fixed throughout.

All angular modes and all original/derivative tensor slots are included.
The registered affine lapse and zero shift are retained. The physical
formula for K includes its half-action factor two, and the varying
volume factor is the actual J^-1/2. This is a completed initial Cauchy
composition interface. It does not quantify the boundary-chart flow,
the quasilinear evolution, an observer tube or a physical observable
remainder; those require their own estimates after (16)--(17).
