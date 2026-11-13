# Complete all-mode momentum H14 inverse bound

13 November 2026 UTC. Reviewed continuation of
[MOMENTUM_ALL_MODE.md](MOMENTUM_ALL_MODE.md). The latter's accepted
all-mode low-order proof is retained unchanged. This companion supplies
the missing high-order estimate, conditional on the same actual smooth
background and selected weighted isomorphism premises.

The completed numerical conclusion, with exactly the physical norms,
fixed canonical density, weights, annulus, reflection group and dotted
domain specified there, is

```
||L_M^-1 f||_(dot H14_qw) <=10^2700 ||f||_(dot H12_qw)    (1)
```

for every Gamma-invariant one-form source f. This is an operator norm
bound on the actual selected inverse. It is not restricted to the
actual second source's angular support, and it is not an admitted
nonlinear amplitude.

## 1. Fixed all-mode transformations and shared quantified inputs

Keep

```
ell=10^-6, kappa_w=10^-8,
r=(z-ell/2)(5ell/2-z)/ell^2,
q=kappa_w^2 r^2, w=exp(-1/r),
gamma0=dy^2+rho^2 kappa_S3, dy=-rho dz,
rho_y=a rho, ds=dy/q, D=d/ds,
M=w sqrt(qrho^3), b=D log M, j=qa, c=q/rho,
h=D log c=q_y-qa.
```

The coordinate s is a complete real line. For every radial/exact
scalar-harmonic block use zeta=c sqrt(lambda), lambda>=8. For every
coexact block use zeta=c sqrt(mu), mu>=9. The monopole has zeta=0.
Both coexact helicities, all multiplicities and all invariant angular
degrees are included. For a coefficient vector X of the physical
one-form set V=M X. For the normalized equation L_M Z=f the
transformed source is g=M f_block. Its weighted L2 norm is exactly
the flat ds L2 norm. An unnormalized physical canonical source would
have a different transformation; it is not substituted here.

Define the weighted derivative-array norm

```
J_m(V)^2=sum_(alpha+beta<=m) ||zeta^beta D^alpha V||_2^2,
P=-D^2+zeta^2, A=1+P.
```

Use the componentwise version for the two-vector radial/exact block
and the Hilbert square sum over all angular blocks. At each fixed
finite eigenvalue, A is the positive self-adjoint operator on L2(R)
associated with its closed quadratic form; its potential is smooth
and bounded. In particular A>=1. The aggregate operator is their
orthogonal direct sum.

Three already-derived quantitative inputs are used, with their actual
proofs rather than an unknown elliptic regularity constant:

1. [MOMENTUM_ALL_MODE.md](MOMENTUM_ALL_MODE.md), equation (15), gives
   `||Z||_w<=1.44*10^38 ||f||_w` for the same selected inverse.
2. [WEIGHTED_NORM_BRIDGE.md](WEIGHTED_NORM_BRIDGE.md) proves, for
   every m<=16 on the entire scalar, radial/exact and coexact spaces,
   `J_m(MZ)<=10^29||Z||_Hm(q,w)` and
   `||Z||_Hm(q,w)<=10^29 J_m(MZ)`. Its actual common constant is
   `13*8^16*16!<10^29`, and it includes all derivative-index connections.
3. The explicit scalar graph-norm proof in
   [SCALAR_ALL_MODE.md](SCALAR_ALL_MODE.md), sections 4--5, gives

```
J_m(V)<=10^(10m)||A^(m/2)V||_2, 0<=m<=16,              (2)
|D(zeta)/zeta|<=10^-5,
|D^2(zeta)/zeta|<=10^-5,
|(zeta^2)''|<=10^-10 zeta^2.                            (3)
```

The relative statements in (3) are interpreted directly at zeta=0.
They hold uniformly in the angular eigenvalue. The graph lemma is
scalar only in the sense that A acts identically on each component;
its proof therefore applies to the vector blocks by square summation.
No fractional power of the matrix momentum operator is used.

For clarity, the coefficient proof behind (3) uses the actual
polynomial ODE on the global variables x=(z-3ell/2)/ell, 1/rho, a,
the scalar gradient and the scalar background. On their unit box the
D-vector field has coefficient sum at most 10^-8 and degree at most
nine. The exact quantities b,j,c,q_y have polynomial degree at most
five and small coefficient norms. Consequently

```
max_(0<=k<=14) ||D^k b||_infty,
max_(0<=k<=14) ||D^k j||_infty,
max_(0<=k<=14) ||D^k c||_infty,
max_(0<=k<=14) ||D^k q_y||_infty <=1.                    (4)
```

The same argument gives the stronger relative polynomial recurrence
`D^k zeta^beta=zeta^beta p_(k,beta)`,
`||p_(k,beta)||_coeff<=(4*10^-6)^k` for k,beta<=32.
It is essential that this is a relative bound, not an absolute bound
on zeta or on the angular eigenvalue.

## 2. The exact normal operator with all angular terms retained

For a radial/exact block the exact transformed first-order adjoint is

```
T=B D+C,
B=[[2,0],[0,sqrt(2)],[0,0],[0,0]],
C=[[-2b,0],
   [sqrt(2)zeta,-sqrt(2)(b+j)],
   [sqrt(12)j,-sqrt(12)zeta/3],
   [0,2c sqrt(2(lambda-3)/3)]].
```

The selected normal equation is L V=g, with
`L=T^*T=-B^T B D^2+(C^TB-B^TC)D+C^TC-B^TC'`.
Writing it out gives

```
L=[[-4D^2+2zeta^2+v1, 2zeta D-2zeta(b+3j)],
   [-2zeta D-2zeta'-2zeta(b+3j), -2D^2+4zeta^2+v2]],

v1=4b^2+12j^2+4b',
v2=2(b+j)^2-8c^2+2b'+2j'.                              (5)
```

In particular the term -8c^2 survives the trace-free Hessian
coefficient lambda-3 and must not be omitted. The derivative of the
upper/lower mixed coefficient also fixes the -2zeta' term in the
lower-left entry.

Introduce the row operator G=(D,-zeta), with its flat formal adjoint.
Then the exact decomposition is

```
L0=2P I+2G^*G,
L=L0+R,
R=diag(v1,v2)+r1 zeta [[0,1],[1,0]],
r1=-2(h+b+3j).                                        (6)
```

Indeed G^*G has upper-right entry D zeta and lower-left entry
-zeta D, so the difference from (5) is exactly the symmetric
coefficient -2zeta(h+b+3j). Reversing the sign of the angular
component conjugates these matrices by diag(1,-1) and is equivalent;
the signs above use the retained positive gradient convention.

For a coexact mode the exact operator is simpler:

```
L=2P+v2, zeta=c sqrt(mu).                               (7)
```

For the radial monopole,

```
L=4P+v1, zeta=0.                                       (8)
```

Equations (5)--(8) describe every angular block of the actual
momentum operator. They do not freeze its coefficients or replace
its angular principal terms by bounded lower-order coefficients.

## 3. Uniform second-order matrix graph estimate

All calculations first concern compact smooth coefficients. Denote
`u=||V||_2`, `p=||PV||_2`. Integration by parts gives the identity

```
||AV||^2=||V||^2+2||DV||^2+2||zeta V||^2
          +||D^2V||^2+||zeta^2 V||^2+2||zeta DV||^2
          -integral (zeta^2)'' |V|^2.
```

By (3), the final term is absorbed into one of the two displayed
`||zeta V||^2` terms. Thus this particular order has the sharper
bound

```
J2(V)<=||AV||<=p+u.                                    (9)
```

In the two-vector block, expand the square of L0 from (6):

```
||L0 V||^2=4p^2+4||G^*G V||^2+8 Re<PV,G^*G V>,
<PV,G^*G V>=<P G V,G V>+<[G,P]V,G V>.
```

The first term in the last line is nonnegative. The commutator is

```
[G,P]=(2zeta zeta', -2zeta' D-zeta'').                  (10)
```

Using (3),
`||[G,P]V||<=6*10^-5 J2(V)` and
`||GV||<=sqrt(2)J1(V)<=sqrt(2)J2(V)`.
Since `8*6*sqrt(2)*10^-5<10^-3`, equations (9)--(10) imply

```
||L0 V||^2 >=4p^2-10^-3 J2(V)^2
            >=(4-.002)p^2-.002u^2 >=p^2-u^2,
p<=||L0 V||+u.                                        (11)
```

No bounded-frequency hypothesis is present: zeta remains in the
positive principal operator throughout.

For the lower-order part, (4), |h|<=2, and (5)--(6) give

```
|v1|,|v2|<=20, |r1|<=12,
||RV||<=20u+12||zeta V||.
```

Positivity of P yields `||zeta V||^2<=<PV,V><=p u`.
Young's inequality in the precise form
`12 sqrt(pu)<=p/4+144u` therefore gives

```
||RV||<=p/4+164u,
p<=||LV||+p/4+165u,
||AV||<=p+u<=(4/3)||LV||+221u
          <=256(||LV||+||V||).                         (12)
```

For the coexact and monopole operators (7)--(8), the same bound
follows directly from |v_i|<=20; their constants are smaller.
Thus (12) is a common evaluated all-mode base estimate.

## 4. Explicit normal-order and commutator ledger

To quantify the higher step, put

```
Q_N(a)=sum_(k=0)^N sup|D^k a|/k!.
```

Use operator norm for matrix coefficients and sum over normal-form
monomials `a_(alpha,beta)(s) zeta^beta D^alpha`. These seminorms are
submultiplicative, and
`Q_N(Da)<=(N+1)Q_(N+1)(a)`.
By (4), for N<=14,

```
Q_N(b),Q_N(j),Q_N(c),Q_N(q_y)<3, Q_N(h)<6.
```

The exact operator (5), or (7)--(8), has weighted order at most two.
Its coefficient sum in Q_12 is below 10^4. Here is an explicit
allocation: `Q_12(b')<=13 Q_13(b)<39`, and likewise for j'.
Each of v1,v2 then has Q_12 norm at most 300. The two off-diagonal
zeta coefficients in (5) together cost at most 60. The principal
coefficients cost at most 10. The sum is below 1000, and 10^4 is
the common rounded ceiling.

For beta<=16, left multiplication of a normal-form term by D is

```
D[a zeta^beta D^alpha]
 =(a'+beta h a)zeta^beta D^alpha
     +a zeta^beta D^(alpha+1).
```

At a derivative-budget step N<=13, its coefficient sum costs at most
`15+16*6+1=112<130`, passing from Q_(N+1) to Q_N.
Multiplication by zeta has cost one. Consequently left multiplication
by A=1-D^2+zeta^2 has cost at most

```
1+130^2+1=16902<20000.                                 (13)
```

The necessary initial jet budgets are available: expanding A^6
with two derivatives left over uses h through order thirteen;
expanding A^6 L differentiates coefficients of L at most twelve
times, requiring b,j through order thirteen. All are within (4).

Six steps in (13) give

```
||A^6 g||<=20000^6 J12(g)<10^26 J12(g).                 (14)
```

The normal-form coefficient sum of A^6 L is below
`10^4*20000^6<10^30`. For L A^6, first expand A^6 with the Q_2
coefficient budget below 10^26. Acting with L costs at most
`10^4*130^2`, since its derivative order is at most two and its
remaining zeta multiplications do not alter coefficients. Thus the
coefficient sum is below 10^35. The difference is in particular
below 10^40.

The commutator has weighted order at most thirteen. To see the
essential cancellation, in highest weighted order replace D and
zeta by commuting variables. The degree-twelve symbol of A^6 is
scalar and therefore commutes with the degree-two matrix symbol
of L. Every actual commutation `D zeta-zeta D=h zeta` lowers the
weighted order by one, as does differentiation of a coefficient.
Hence no order-fourteen term survives. We have proved the evaluated
estimate

```
||[L,A^6]V||<=10^40 J13(V).                             (15)
```

This is a bound for the full matrix commutator. The mixed principal
terms zeta D are included, and no factor lambda has been bounded
by a constant.

## 5. Fourteen derivatives by one scalar spectral absorption

Apply (12) to A^6 V, commute L through A^6, and use (14)--(15):

```
||A^7 V||
 <=256[10^26 J12(g)+10^40 J13(V)+||A^6V||].
```

Use (2) at order thirteen and A>=1. Rounding the finite prefactors
upwards gives

```
||A^7 V|| <=10^29 J12(g)+10^174||A^(13/2)V||.           (16)
```

For t>=1, splitting at t^(1/2)=epsilon^-1 proves
`t^(13/2)<=epsilon t^7+epsilon^-13`. The spectral theorem and the
triangle inequality therefore give

```
||A^(13/2)V||<=epsilon||A^7V||+epsilon^-13||V||.
```

Set K=10^174 and epsilon=(2K)^-1. The K epsilon term in (16)
is exactly 1/2 and can be absorbed. It follows that

```
||A^7 V|| <=2*10^29 J12(g)+2^14 K^14||V||
            <=10^2442[J12(g)+||V||],
J14(V) <=10^2582[J12(g)+||V||].                        (17)
```

For the arithmetic, `K^14=10^2436` and `2^14<10^5`; the final
graph conversion costs 10^140. These estimates hold with the same
constants in every angular block. Orthogonal direct sums and
Minkowski's inequality preserve (17) for the full all-mode arrays.

There is no regularity circularity in this a priori argument. It can
first be made for compact smooth potentials and finite angular sums.
Normal-ordering controls all expressions by J14. The norm bridge
identifies the original dotted H14 approximations with this weighted
array completion, and the retained continuous map L_M:H14->H12
allows passage to the actual selected solution. Equivalently the
finite operator expansions show that A^7 V and L A^6 V are defined
in L2 for the relevant approximants and converge in the indicated
array norms. The qualitative selected H14 solution was an explicit
premise; an unknown numerical high-order bound was not inserted.

## 6. Exact geometric norm and source closure

Apply the two directions of the norm bridge to (17), remembering
that g=M f is the normalized source. With V=MZ this gives

```
||Z||_H14(q,w)
 <=10^29 J14(V)
 <=10^2611[J12(Mf)+||Z||_w]
 <=10^2611[10^29||f||_H12(q,w)+1.44*10^38||f||_w]
 <=10^2650||f||_H12(q,w)
 <=10^2700||f||_H12(q,w).                              (18)
```

The penultimate rounding uses
`10^29+1.44*10^38<10^39`. The last fifty decades are a stated
rounding reserve. All inputs to (18) have been quantified, including
the all-mode gap, endpoint-degenerate principal terms, coefficient
jets, derivative-index and source norm conversions, matrix
commutator, and spectral absorption.

By uniqueness, the operator bounded here is precisely the retained
selected inverse on the invariant dotted space. The canonical
momentum correction is still
`delta pi=rho0 w^2 q^2 P_m^*Z`, with the same fixed density rho0.
The correction map and its physical zero extension have not been
changed. The estimates extend by density and remain valid at both
annulus endpoints without assigning ordinary endpoint traces to Z.

This closes the momentum block of the requested all-mode initial
inverse certificate under the retained analytic premises. The scalar
block is treated separately in SCALAR_ALL_MODE.md. Neither block
alone, nor their direct sum, evaluates the nonlinear product-ball
constants or supplies a nonzero amplitude radius or evolved tube.
