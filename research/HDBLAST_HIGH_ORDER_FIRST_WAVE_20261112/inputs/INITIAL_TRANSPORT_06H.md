# Quantitative initial-data transport through the exact boundary chart

6 November 2026. This derivation is conditional on the same smooth background
and selected initial canonical curve as PR #133. It bounds amplitude derivatives
at zero, not a finite-amplitude solution or an evolved norm. Public source bytes
are pinned in SOURCE_PINS.json; REPRODUCE.md describes the finite exact replay.

## Fixed norms and input interface

Work on the physical one-sided collar `Sigma=[0,6 ell] x S3`, `ell=10^-6`.
Set `q=dz^2+kappa`, and `gamma0=dy^2+rho^2 kappa=rho^2 q`, with
`dy=-rho dz`, `129<=rho<=131`. The source raw tensor connection is
`D=nabla_gamma0`, extended by `D dt=0`, with positive fiber metric
`dt^2+gamma0`; the scalar has its usual connection. Its measure is the full
`dvol_gamma0=rho^4 dz dmu_kappa`. All Hm norms are the square sum of
covariant derivative orders 0 through m, with full sphere measure.

The target fixed bundle is

```
E = R_v + Sym^2(R dt + T*S3)_H + R_w + (R dt + T*S3)_D.
```

Use the direct-sum positive fiber norm, the full Frobenius norm on its
symmetric tensor part, product connection `bar D`, and measure `dz dmu_kappa`.
This is a precise fixed-bundle Sobolev norm compatible with the energy
construction; symmetric off-diagonal components are counted twice in the
Frobenius square. It is not the independent-entry coordinate Euclidean norm.

Write `X_n=(g_n,phi_n)` and `dot X_n=(partial_t g_n,partial_t phi_n)` for
amplitude derivatives at zero in the fixed background identification. The
input interfaces are

```
W2 = ||X2||H14 + ||dot X2||H13 < 10^4000,
W1 = ||X1||H14 + ||dot X1||H13 < 10^1500,
C1 = max_{j<=14} ||D^j X1||infty
     + max_{j<=13} ||D^j dot X1||infty < 10^1500.
```

The first inequality is the rounded raw-wave estimate in
[GAUGE_CAUCHY_COMPLETION_06G.md](inputs/GAUGE_CAUCHY_COMPLETION_06G.md).
The latter two follow from the pinned 06e source engine and the complete first
harmonic jet, not from W2. Write `Y=n1²-n2²`, `w_velocity=ut` and
`S=ut Y/rho`. In the raw orthonormal splitting the first fields are

```
g_tt,1=-2rho²uY, g_ti,1=0,
gamma1=uY(-4dy²+2gamma_angular), phi1=vY,
g_tt,t,1=-2rho²utY,
gamma_t,1=2utY(-2dy²+gamma_angular),
g_ty,t,1=-6rho²(u_y+4au)Y, g_tA,t,1=0, phi_t,1=vtY.
```

The mixed velocity is compulsory: the first mixed metric values vanish
while their time velocities need not vanish. The exact global Taylor replay
uses ordinary proper-y derivatives and gives these strict decade ceilings:

| Profile | Maximum derivative order | Ceiling |
| --- | ---: | ---: |
| u | 14 | 10^86 |
| v | 14 | 10^97 |
| ut | 13 | 10^79 |
| vt | 13 | 10^90 |
| rho²u | 14 | 10^90 |
| rho²ut | 13 | 10^83 |
| rho²(u_y+4au) | 13 | 10^90 |

All are below 10^100. The source engine stores sixteen normalized orders;
`u_y` uses one extra stored derivative before the thirteen velocity orders.
It bounds all centers of the compact packet, not a sample grid.

For completeness, the full covariant derivative conversion uses the retained
warped identity `D_0=partial_y`, `D_A=rho^-1 Dbar_A+a K_A`. The background
proper-y jets of a and rho^-1 through eighteen are below 10^80. At each of
fourteen stages, 1000 bounds the derivative/slot expansion coefficient sum
at ranks at most sixteen, and a term has at most fourteen background factors.
The angular Y and fixed radial/sphere tensors have reference angular jets
below 10^200, by `Dbar²n_i=-n_i kappa` and the finite product rule. This
charges the derivative slots as well as the original tensor slots. Thus the
pointwise derivative bound before displayed coefficients and finite sums is
`1000^14 * 10^(14*80+100+200)=10^1462`.

The known packet is flat at its support endpoints. Its enclosing annulus
has geometric volume below `40ell*131^4<10^6`, hence square-root volume below
10^3. A factor 1000 covers the finitely many displayed terms, their coefficients
and symmetric mixed slots; the derivative-order square root is below four.
Both the sum C1 of the two supremum budgets and W1 are bounded by the same
conservative expression

```
4*10^6*1000^14*10^(14*80+100+200) < 10^1500.
```

This proves the first-data interfaces in the specified raw norm. It supplies
smooth multiplier bounds for this explicit initial packet only.

The output claim is

```
T1 = ||V1||H14(E,q) + ||partial_t V1||H13(E,q) < 10^5200,
T2 = ||V2||H14(E,q) + ||partial_t V2||H13(E,q) < 10^10000.
```

The same rounded ceilings cover two identical reflected copies. Restriction
to the physical collar can only reduce the global input norms. The first seed is supported in `ell<=z<=2ell`; the second correction is
flatly supported in `ell/2<=z<=5ell/2`. Both lie inside the physical collar
where the retained parameter is exactly theta=z. Since F_-z(0)=0, the
transformed jets vanish off those supports. No derivative of the cutoff
outside that collar is used. On the second reflected copy use its positive
inward coordinate; no negative-physical-z background is asserted.

## Raw physical bundle to the fixed product bundle

Let `I` be the fiber isometry that sends the radially parallel coframe
`(dt,dy,rho theta^A)` to `(dt,dz,theta^A)`. This is a fiber identification,
not the coordinate pullback `dy=-rho dz`. It acts on every covariant slot,
including derivative slots. For a rank-r tensor T, the warped connection
formula gives

```
partial_z I(T) = -rho I(D_y T),
bar D_A I(T) = rho I(D_{e_A} T) - rho a K_A I(T),
e_A=rho^-1 E_A, a=rho_y/rho.
```

The sign on the first line is essential. Each `K_A` mixes the radial and
angular slots; its operator norm is bounded by `3r`. Combining three angular
directions and the radial direction gives the deliberately enlarged bound

```
|bar D I(D^k T)| <= 2 rho |D^(k+1)T|
                    + 10(r+k) rho |a| |D^k T|.
```

The differentiated identity is used on `D^k T`, so all derivative slots are
charged. `r+k<=16` in the present comparison. The radial and angular splitting
maps, signs, contractions, and permutations are product-parallel maps.

PR133's background z-jet certificate gives ordinary derivatives of `rho,a`
through order 19 below `10^100`. Hence through order 14,

```
max |partial_z^j rho| < 10^210,
max |partial_z^j(rho a)| <= 2^14 10^200 < 10^210.
```

In an expansion of j product derivatives, each term has at most j
background factors (each an ordinary derivative of rho or rho a), and one
transported `D^k T`. At the next stage there are at most j ways to hit such
a factor, one leading derivative term, and slot mixing bounded by
`10(r+k)<=160`; the full-direction leading factor is at most 2. Thus 10^4
uniformly bounds the sum of all branch and slot coefficients at each of the
14 stages. No angular coordinate charts or partition of unity are involved.
Induction therefore gives

```
|bar D^j I(T)| <= (10^4 10^210)^j sum_{k=0}^j |D^k T|, j<=14.
```

Because `rho>=1`, changing source L2 measure to target measure costs at most
one. Applying Cauchy-Schwarz to the finite derivative sums and summing the
15 output orders costs less than 10^3. Since `14*214+3=2999`, the strict
common estimate is

```
||I(T)||H_m(q) <= 10^3000 ||T||H_m(gamma0), m<=14,
max_{j<=m}||bar D^j I(T)||infty
  <=10^3000 max_{j<=m}||D^j T||infty.
```

For the second inequality the at-most-15 summands are included in the same
10^3 reserve. Scalar, tensor, and component sums have additional reserve in
the rounded exponent 3000. Consequently, setting `Z_n=I(X_n)`,

```
||Z2||H14+||dot Z2||H13 < 10^7000,
||Z1||H14+||dot Z1||H13 < 10^4500,
max_{j<=14}||bar D^j Z1||infty
 +max_{j<=13}||bar D^j dot Z1||infty < 10^4500.
```

## Exact chart and the two coefficient maps

For a physical metric perturbation in the I frame define the linear map
`L(z)` into root-unprocessed normalized variables by

```
x     = delta g_yy,
H_tt  = delta g_tt/rho_b^2,
H_tA  = rho delta g_(t,phys A)/rho_b^2,
H_AB  = rho^2 delta g_(phys A,phys B)/rho_b^2,
w     = delta phi,
D_t   = delta g_yt/rho_b,
D_A   = rho delta g_(y,phys A)/rho_b.
```

Here `rho_b=rho(0)` is a fixed constant. In particular the cross components
have exactly the proper-y convention of the retained generator. The map
does not divide g_ya by rho(z) or replace the component root by a Schur lapse.
Introduce

```
r(x,H,w,D)=(sqrt(1+x)-1,H,w,D),
R1=Dr(0), R2=D^2r(0),
(R1 u)_v=u_v/2, (R2(u,v))_v=-u_v v_v/4,
```

with R1 the identity on all other components and R2 zero there. Let
`F=F_-z` be the actual inverse boundary flow. Since `F_-z(0)=0`,

```
A(z)=D F_-z(0) R1 L(z),
B(z)[P,Q]=D^2 F_-z(0)[R1 L P,R1 L Q]
             +D F_-z(0) R2(L P,L Q).
```

The amplitude chain rule is exactly

```
V1=A Z1,
V2=A Z2+B(Z1,Z1).
```

Second derivative means the coefficient of e^2/2; there is no extra
one-half in the second formula. These are all the amplitude Faà di Bruno
terms at orders one and two. For 14 spatial derivatives the full two-factor
and three-factor Leibniz sums below charge all remaining partitions; one
does not need an unproved nonlinear H14 composition theorem at nonzero e.

The maps L, R1, R2, A, and B are natural under orthogonal rotations of the
three angular tangent directions. The autonomous generator's operations
are matrix inversion, multiplication, and contraction with J4, so its flow
and derivatives commute with these rotations. Consequently A and B are
parallel in angular directions with the induced product connections on
all input/output slots. Their base derivatives at t=0 are ordinary z
derivatives; there is no hidden angular-frame derivative in this assertion.

PR133's background jets and Leibniz give, with room for every component and
symmetric multiplicity,

```
max_{j<=14} ||partial_z^j L||op < 10^210.
```

The 06f flow certificate supplies mixed z/state derivatives through total
order 18 below 10^126 in star operator norm. For the Frobenius direct sum,
`star<=|.|E<=sqrt(7) star`; hence for each required derivative (j<=14 and
state order 1 or 2), the E multilinear operator norm is below `3*10^126`.
The expressions below use the larger factor 10^6 for this norm conversion
and all fixed component sums.
No erroneous use of the independent-coordinate infinity Jacobian as 1.01
occurs. State derivative order two plus spatial order fourteen uses only
total order sixteen, within the available eighteen.

Using ||R1||<=1 and ||R2||<=1/4 and ordinary Leibniz sums,

```
max_{j<=14} ||bar D^j A||op <= 10^6*2^14*10^(126+210) < 10^350,
max_{j<=14} ||bar D^j B||op <= 10^6*3^14*10^(126+2*210) < 10^560.
```

The common coefficient budget `M=10^600` is therefore valid, in precisely
the fixed bundle used above.

## Initial time velocities

The chart uses the fixed q_dS coframe `(dt,cosh(t)theta^A)`. Relative to
the fixed initial raw frame, its time dependence is through powers of
`cosh(t)^-1`. The radial component root and the autonomous F_-z are
time-independent. The background metric and scalar have zero first time
derivative at t=0. Thus the perturbation is zero for all background times,
and the induced coefficient maps satisfy

```
partial_t A|t=0=0, partial_t B|t=0=0.
```

This is an even-frame identity at the initial section; it is not an
assumption that the perturbed metric is even in time. Since I is fixed at
the initial section, dot Z_n=I(dot X_n). Differentiating the exact amplitude
formulas gives

```
partial_t V1=A dot Z1,
partial_t V2=A dot Z2+2 B(Z1,dot Z1).
```

The factor two is mandatory. In a different moving frame or a nonzero
initial time the general formulas also contain A_t Z_n and B_t(Z1,Z1),
which cannot simply be discarded. The proof here explicitly establishes
their vanishing for this data section.

## Sobolev products and final arithmetic

For a coefficient map A with product covariant derivatives through m
bounded by M, tensorial Leibniz and the above Hm convention give

```
||A P||Hm <= 2^m sqrt(m+1) M ||P||Hm <= 10^6 M ||P||Hm,
```

for m<=14. For a bilinear coefficient B and a fixed smooth P,

```
||B(P,Q)||Hm <= 3^m sqrt(m+1) M
               max_{j<=m}||bar D^j P||infty ||Q||Hm
             <=10^8 M C_m(P)||Q||Hm.
```

The exact bound uses the multinomial sum over coefficient/P/Q derivative
allocations. It applies for arbitrary Hm Q, and therefore does not promote
an H14 estimate to W14,infinity. The infinity norm is used only on the
explicit smooth first packet. These finite products make no statement
about evolved X14 commutators or normal recovery.

The converted input budgets now yield

```
T1 <10^(6+600+4500)=10^5106 <10^5200,
T2 <10^(6+600+7000)
       +3*10^(8+600+4500+4500)
   =10^7606+3*10^9608 <10^9610 <10^10000.
```

The coefficient 3 charges the value quadratic once and the time cross
term twice; using sums of value and velocity input norms only improves
this conservative estimate. The margins 94 and 390 decades easily absorb
the sqrt(2) gain for two identical reflected copies (even a factor two
would fit). One must retain the same rounded constants only for that
specified identical union, not for an arbitrary number of copies.

## What this closes and what remains

This closes the raw-to-straightened H14/H13 initial first/second jet
interface, conditional on the analytic premises in the pinned sources. The
finite replay checks the displayed arithmetic, not those existence premises. It
uses the actual proper-y radial root, frozen-shell F_-z, full geometric
raw measure, fixed product target measure, frame conversions, tensor
derivative slots, and initial velocity cross term.

It does not evaluate the complete ordered PDE lower maps, the quadratic
first-wave forcing after evolution, high-order energy products and normal
recovery, an actual-family patch/bootstrap radius, an admitted nonzero
amplitude, or the proper-time/area-normalized observable remainder.
The large constants are analytic ceilings in fixed model normalization;
they are not physical blast energy or evidence of a cosmological origin.
