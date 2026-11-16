# Actual common coordinate-probe clocks and moving-field composition

16 November 2026. New proof at immutable public base
`5be99501c7f58e4400212f374abb6334c22cbca5`. Exact source identities are
recorded in `SOURCE_PINS.json`.
All retained background, canonical-family, coefficient-chart, frozen-linear
construction, norm and physical-reconstruction premises remain in force.
This proves an actual coordinate-probe clock lemma and finite-regularity
composition lemmas. It does not identify these accelerated probes with
intrinsic shell geodesics or ambient free fall.

## 1. Inputs and statement

Use the actual original-time physical curve `X_e=(g_e,phi_e)` on

```
[-T,0] x Sigma, Sigma=[0,6ell] x S3, ell=10^-6, T=13ell/8,
I=[-a_iter,a_iter], a_iter=10^(-10^500209).
```

The retained evolved-parameter theorem gives C2 dependence into `C X10`
on a qualitative open parameter interval J containing I. Its reduced
extension reconstructs these same fields and has the same strict state
margin. Full unreduced equations and shell laws are needed here only on
the original I; no new physical parameter interval is asserted.
`C X10` is the compatible-temporal-jet Banach space, with fixed positive
product Sobolev norms in four spatial dimensions. No parameter continuity
in the individually available X14 norm is used.

Fix `z_* in [0,6ell]` and the prescribed coordinate congruence

```
chi(t,n)=(t,z_*,n), n in S3.
```

At z_*=0 this is a curve in the actual identified physical shell.
At an interior z_* it is a curve in one bulk copy. Its label n is the
fixed round label, not a freely transported geodesic label. Define

```
b_e(t,n)=sqrt(-g_e,tt(t,z_*,n)),
f_e(t,n)=tau_e(t,n)=integral_0^t b_e(r,n) dr.           (1)
```

There is an open interval J_c containing I, with compact closure inside
J, and a common positive b_*, such that b>=b_* on the original time
slab for every e in closure(J_c) and every label. For every symbolic
choice `0<tau_*<b_* T`, the actual common inverse exists on

```
C_tau=[-tau_*,0] x S3,
t=T_e(tau,n), f_e(T_e(tau,n),n)=tau.                  (2)
```

It has the precise regularity

```
T_e in C7(C_tau) for every e,
e -> T_e is C2 into C5(C_tau),
e -> partial_e T_e is continuous into C6(C_tau).      (3)
```

Its second parameter derivative is continuous into C5. All displayed
norms and derivatives are uniformly bounded on compact parameter
subintervals of J_c, hence on the entire original I. Endpoint parameter
derivatives on I have their ordinary two-sided meaning from J_c. Clock
derivatives at tau=0 are one-sided time derivatives from the past. No
future physical solution is used to construct (2).

No numerical duration tau_*, calibration, effective new amplitude ceiling,
or historical observer duration is evaluated by this theorem. The clock
is the physical proper clock of the declared coordinate probe. These
probes generally accelerate and differ from the intrinsic normal
congruence constructed by a separate observer ODE argument.

## 2. Bounded linear joint-C7 embedding

For a compatible table in C X10, every mixed derivative with time
index j and spatial index count k, `j+k<=7`, obeys

```
||D_x^k partial_t^j X||C0(Sigma)
 <=C ||partial_t^j X||H^(k+3)
 <=C ||partial_t^j X||H^(10-j).                       (4)
```

The integer H3 reserve is above the continuity threshold in dimension
four. Every jet is strongly time continuous in its own Sobolev space,
so (4) gives joint continuity in `(t,x)`. Distributional compatibility
identifies these continuous entries as the classical mixed derivatives
of one field, including the available one-sided derivatives at the time
and radial boundaries. Fixed smooth coordinate/bundle changes and a
finite chart cover are bounded on these finite derivative norms.
Thus the embedding

```
C X10 -> joint C7([-T,0] x Sigma)
```

is bounded and linear. Applying it to the actual C2 curve proves

```
e -> (g_e,phi_e) is C2 into joint C7.                 (5)
```

Restriction to the fixed shell preserves C7. Equivalently the shell
value trace of the jth time jet is `H^(19/2-j)(S3)`; a further k
tangential derivatives with j+k<=7 leave at least H^(5/2), above the
three-dimensional continuity threshold. The critical eighth-order
trace does not justify joint C8 by this argument.

On the admitted inverse patch, the finite-dimensional inverse metric
is a smooth Nemytskii map on the C7 Banach algebra. Differentiation then
gives the safe fixed-base statements

```
g^-1 is C2(e;C7), Gamma[g] is C2(e;C6), Riem[g] is C2(e;C5),
h=iota^*g is C2(e;C7),
Gamma[h] is C2(e;C6), Riem[h] is C2(e;C5),            (6)
```

where `iota(t,n)=(t,0,n)` is the retained actual physical shell
embedding. All inverse denominators remain bounded away from zero on
compact admitted sets. The intrinsic shell connection and curvature in
(6) are not ambient bulk curvature. No derivative across the Israel
interface is asserted.

## 3. Uniform timelike margin and common clock image

The exact retained reconstruction has time coframe entry dt and

```
||rho(z)^-2 G-J5||op<1/40.
```

Consequently on the retained I,

```
g_e(partial_t,partial_t)<-(39/40)rho(z_*)^2<0.         (7)
```

The lapse b in (1) is this probe's clock rate. It is distinct from
the radial root `N=sqrt(g_yy)` and from the temporal ADM lapse alpha;
the latter equals this clock rate only in the appropriate zero-shift
description. The chart flow `F_z` in field space is not a spacetime
flow moving these probe points. Since the background coordinate change
`dy=-rho(z) dz` is parameter independent, fixing z_* and fixing its
background y value are the same point prescription. Neither makes y
perturbed proper distance.

By (5), b is an actual C2 parameter curve into C7 on every compact
admitted parameter neighborhood where the argument of its square root
stays positive. Compactness of I and the strict gap (7) allow a smaller
open J_c containing I, with compact closure in J, for which this holds
uniformly with some b_*>0. This step uses continuity of the reconstructed
reduced extension, not a new physical theorem outside I.

For t<=0, the clock is negative into the past, and

```
f_t=b>=b_*,
f_e(0,n)=0,
f_e(-T,n)<=-b_* T<-tau_*.
```

The intermediate value theorem and strict monotonicity prove a unique
solution of (2) for every `(e,tau,n)` in the asserted product domain.
Let b^* be the finite maximum rate on the same compact parameter and
spacetime set. The mean-value bound gives the common past-edge gap
`T_e(-tau_*,n)+T>=(b_* T-tau_*)/b^*>0`; the anchor value is exactly
`T_e(0,n)=0`.
In particular no labelwise clipping or parameter-dependent loss of a
clock endpoint is hidden in the definition.

Square root is smooth on the positive C7 cone. The integration operator
`b -> integral_0^t b` is bounded linear into C7 on this fixed compact
cylinder, since a positive time derivative differentiates the integral
and zero time derivatives only integrate spatial derivatives. Thus the
clock curve f is C2 into C7. At fixed original coordinate t and label,

```
b_epsilon=-(g_tt)_epsilon/(2b),
b_epsilonepsilon=-(g_tt)_epsilonepsilon/(2b)
                         -((g_tt)_epsilon)^2/(4b^3),
f_epsilon=integral_0^t b_epsilon,
f_epsilonepsilon=integral_0^t b_epsilonepsilon.        (8)
```

The inverse powers of the actual positive clock gap are genuine costs
in uniform bounds. The formulas are not evaluated physical constants.

## 4. Finite-regularity composition and inversion

For a compact parameter/label domain K and a source domain D, suppose
the images of gamma lie in one common compact interior subset of D.
For each nonnegative integer k, the map

```
Comp:(A,gamma) -> A o gamma,
C^(k+2)(D) x C^k(K;D) -> C^k(K)
```

is C2 on its open image-admitted domain. Indeed its differentials are

```
D Comp[(a,Y)]=a o gamma+(dA o gamma)Y,
D2 Comp[(a,Y),(c,Z)]
 =(da o gamma)Z+(dc o gamma)Y+(d2A o gamma)(Y,Z).      (9)
```

Differentiate each expression at most k times on K, keeping all slots.
There are finitely many chain/Leibniz terms, with at most k+2 source
derivatives. Products are bounded in the Ck Banach algebra. The ordinary
finite-dimensional Taylor formulas, with uniform continuity of the
highest derivatives on the common compact image, prove the Frechet
remainders and continuity of the second differential. This proves (9)
as a Banach-space statement, not only pointwise Gateaux calculus. Fixed
finite charts and transitions give the same result for the tensor
component arrays and compact manifolds used here.

We handle the t=0 endpoint before applying (9). A fixed bounded linear
C7 extension across t=0 can be constructed explicitly: for 0<t<delta,
put

```
E b(t,n)=sum_(l=1)^8 c_l b(-l t,n),
sum_(l=1)^8 c_l(-l)^j=1, j=0,...,7,
```

and let E b=b on the original side, with delta<T/8. The Vandermonde
system is invertible, so this matches all seven boundary derivatives.
It is a bounded extension on this fixed finite strip. Compactness and
the positive anchor gap allow a sufficiently small common delta for
which the actual extended b remains positive. Integration from zero
defines the extended f and keeps f(0,n)=0. This extends coefficient
functions for calculus only, and does not posit a future Einstein field.

In the extended cylinder the range of the inverse at any fixed e is
compactly inside the time source interval, including at tau=0. Apply
(9), with k=5, to

```
H(f,T)=f o (T,id)-tau, C7 x C5 -> C5.
```

Its T differential is multiplication by `b o (T,id)`, with bounded
inverse multiplication by `1/(b o (T,id))`. The Banach implicit-function
theorem proves C2 local dependence of T into C5. The unique monotone
inverse makes the local definitions agree on all of J_c. The inverse
in (2) is therefore the actual global curve in (3).

At tau=0 the equation gives T=0. Since the extended clock is increasing
and f(0,n)=0, negative tau has negative inverse t. The past-edge margin
already proved prevents an exit at -T. Thus every value of the final
inverse uses the original physical past interval, despite the auxiliary
extension used to justify differentiation. A direct anchored one-sided
implicit theorem gives the same inverse.

Classical finite-dimensional inversion gives C7 for each individual T.
Its continuity in C7 as e varies follows without derivative loss in the
continuity statement: uniform inverse convergence follows first from
the b_* gap; successively differentiate (2) through seven orders, solve
each leading inverse derivative using the same positive factor b, and
use previously obtained convergence and uniform continuity of the
highest source derivatives. The explicit first parameter formula below
then gives continuous C6 regularity, while its second formula is C5.
The inversion map is not asserted C2 into the unreduced C7 target.

## 5. Complete inverse, label and moving-origin formulas

All following right-hand quantities are evaluated at `(e,T_e(tau,n),n)`
unless another argument is shown. Write `U=T_epsilon`,
`V=T_epsilonepsilon`. Differentiating (2) once and twice proves

```
U=-f_epsilon/b,
V=-(f_epsilonepsilon+2b_epsilon U+b_t U^2)/b
 =-f_epsilonepsilon/b+2b_epsilon f_epsilon/b^2
                                -b_t f_epsilon^2/b^3,
T_tau=1/b,
T_epsilon,tau=-(b_epsilon+b_t U)/b^2,
T_A=-f_A/b,
T_AB=-(f_AB+b_A T_B+b_B T_A+b_t T_A T_B)/b,
T_epsilon,A=-(f_epsilon,A+b_epsilon T_A
                                  +b_A U+b_t U T_A)/b.         (10)
```

A and B are label indices. These identities include every parameter/
label allocation; in particular the two cross terms in the second
parameter derivative are retained. The general formula, for variables
y other than the inverted t, is

```
T_i=-F_i/F_t,
T_ij=-(F_ij+F_it T_j+F_jt T_i+F_tt T_i T_j)/F_t,
F(y,t)=f(y,t)-tau.                                   (11)
```

For y=tau, `F_tau=-1` must be included. A derivative at fixed actual
tau cannot be exchanged for a derivative at fixed original coordinate t.

For a different prescribed moving anchor/clock origin, write

```
f_e(t,n)=c_e(n)+integral_(a_e(n))^t b_e(r,n) dr.
```

At fixed t the exact additional terms are

```
f_epsilon=c_epsilon+integral_a^t b_epsilon-b(e,a)a_epsilon,
f_epsilonepsilon=c_epsilonepsilon+integral_a^t b_epsilonepsilon
 -2b_epsilon(e,a)a_epsilon-b_t(e,a)a_epsilon^2
                                        -b(e,a)a_epsilonepsilon. (12)
```

A moving spatial anchor adds its ordinary moving-evaluation terms.
Equations (11)--(12) are calculus identities under a common admitted
domain, not a declaration of a new experiment or common image for
arbitrary anchors. The actual lemma (1)--(3) has a=0 and c=0.

On each compact parameter interval define the modulus of continuity
of `partial_e^2 T` in C5. It tends to zero by (3). The Banach fundamental
theorem therefore gives uniform Lipschitz bounds for T and its first
derivative, and the uniform Peano remainder

```
||T_(e+h)-T_e-h T_epsilon-(h^2/2)T_epsilonepsilon||C5
 <=(h^2/2) sup_(theta in [0,1])
       ||T_epsilonepsilon(e+theta h)-T_epsilonepsilon(e)||C5
 =h^2 o(1).                                         (13)
```

Only two parameter derivatives have been used. Their two-sided endpoint
meaning follows from J_c; restricting secants to I gives the same limits.

## 6. Actual moving composition and tensor pullback

Let `chi_e(t,n)` now be any supplied sufficiently regular path map on a
common admitted source domain, and set

```
gamma_e(tau,n)=chi_e(T_e(tau,n),n).
```

For the actual coordinate probes chi is fixed as in Section 1. They
have actual unit future velocity `gamma_tau=b^-1 partial_t`, since
`g(gamma_tau,gamma_tau)=-1`; `gamma_r=-gamma_tau` is past directed for
r=-tau. For a different moving chi, the same inverse from Section 1
defines only a parameterized moving evaluation. For tau to be that
path's actual proper time, T must instead invert chi's own clock, with
rate `sqrt(-g(chi_t,chi_t))`, and a common image must be proved for
that clock. The regularity in (3) then applies only if that clock also
has the required C2 C7 bound and positive uniform gap. Such a bound
is not inferred for an arbitrary chi from its position alone.
The following statements are conditional on the declared path/clock
regularity and common image; they do not construct those paths.
The exact derivatives at fixed tau and label are

```
Y=gamma_epsilon=chi_epsilon+chi_t U,
Z=gamma_epsilonepsilon=chi_epsilonepsilon
                         +2chi_epsilon,t U+chi_tt U^2+chi_t V. (14)
```

For a scalar or one fixed-coordinate component A,

```
(A_e o gamma_e)_epsilon=A_epsilon+dA[Y],
(A_e o gamma_e)_epsilonepsilon
 =A_epsilonepsilon+2dA_epsilon[Y]+d2A[Y,Y]+dA[Z].      (15)
```

Every right-hand quantity is evaluated at the moved point. The partial
parameter derivatives of A there remain fixed-base field derivatives.
At the anchor tau=0, apply the same fixed bounded time-reflection
extension used in Section 4 to each source scalar/metric component
array in C7, and its C5 version to curvature component arrays, before
using the open-image composition lemma. These extensions act in the
fixed product bundle and match all available original past derivatives.
They preserve the C2 parameter curves by bounded linearity. Monotone
retention of the actual coordinate inverse ensures that all final
values and endpoint derivative jets are the original past fields.
No future physical field equation or geometric identity is imposed
on these auxiliary extended arrays. At an additional radial boundary
one may likewise use a fixed finite-order extension if a moved map
requires an open radial source domain; the actual coordinate probes
do not move in the radial direction.
At the actual coordinate probes, gamma=(T_e,z_*,n), so

```
O_epsilon=A_epsilon+A_t U,
O_epsilonepsilon=A_epsilonepsilon+2A_epsilon,t U
                                            +A_tt U^2+A_t V. (16)
```

The second clock term `A_t V` is generally present. A special constant
background scalar may simplify a derivative at zero; it cannot remove
the corresponding general nonzero-parameter terms.

If an independently moving timelike chi is parameterized by t, its
clock rate is instead `b=sqrt(S)`, `S=-g(chi_t,chi_t)`. Put
`L=chi_t`, `Y=chi_epsilon`, `Z=chi_epsilonepsilon`, and

```
H1=g_epsilon+dg[Y],
H2=g_epsilonepsilon+2dg_epsilon[Y]+d2g[Y,Y]+dg[Z].
```

The complete metric/path/velocity contributions are

```
S_epsilon=-H1(L,L)-2g(L_epsilon,L),
S_epsilonepsilon=-H2(L,L)-4H1(L_epsilon,L)
                    -2g(L_epsilonepsilon,L)-2g(L_epsilon,L_epsilon),
b_epsilon=S_epsilon/(2b),
b_epsilonepsilon=S_epsilonepsilon/(2b)-S_epsilon^2/(4b^3). (17)
```

For a covariant tensor pullback, derivatives of gamma enter each input
slot as well. For a metric on labeled cuts, put `P=D_n gamma`,
`P1=D_n Y`, `P2=D_n Z`, and `H=g_e o gamma_e`; let H1,H2 be its
complete moving-evaluation derivatives from (15). Then

```
k=P^T H P,
k_epsilon=P^T H1 P+P1^T H P+P^T H P1,
k_epsilonepsilon=P^T H2 P+2P1^T H1 P+2P^T H1 P1
                   +P2^T H P+P^T H P2+2P1^T H P1.   (18)
```

Thus tensor pullback is not merely component evaluation. Velocity,
screen-frame and any contraction slots of an observable likewise vary.

For the actual coordinate inverse map, (3), (5)--(6), and the proved
composition lemma give the following safe uniform statements on the
common clock strip:

| Operation | C2 parameter target |
|---|---|
| Scalar phi or metric component evaluated at `(T_e,z_*,n)` | C5 |
| Bulk one-sided or intrinsic shell curvature component at that point | C3 |
| Covariant metric pullback, including `D_n gamma` | C4 |
| Covariant curvature pullback with the same input slots | C3 |

For a general supplied map gamma in `C2(e;C^k)`, source C^(k+2)
regularity suffices for C2 Ck moving evaluation; full covariant pullback
also needs its derivative slots and may cost one more target order.
In particular a C2 C4 observer position map and C2 C3 label derivative
map suffice for C2 C3 curvature and metric cut pullbacks from (5)--(6).
This is the ready-to-use composition interface for the separate actual
intrinsic observer construction. No C2 C5 curvature pullback is inferred
from only fixed-base C5 curvature.

If a separate positive cut-metric gap is established, determinant and
area-density maps are smooth on that positive cone in the indicated
Ck algebra. Coordinate-probe equal-clock cuts are not automatically
spacelike, so this clock lemma does not supply that extra hypothesis.

## 7. Anchoring, proper-time units and scope

The selected initial fields equal the background on an open shell
collar. At the fixed physical shell anchor `(t,n)=(0,n)`,

```
h0=rho_b^2(-dt^2+kappa),
u_future(0,n)=rho_b^-1 partial_t,
partial_e anchor=partial_e^2 anchor=0,
partial_e u_future(0,n)=partial_e^2 u_future(0,n)=0.   (19)
```

These are actual collar facts for this selected cut. At a general fixed
nonzero-state cut the future shell normal is instead
`u^a=-h^(a0)/sqrt(-h^00)`, and its state and evaluation derivatives
belong in the initial data. A moving cut varies its embedding too.

The clocks in (1) use physical signed proper time tau. With physical
past duration `r=-tau`, an intrinsic shell unit-speed geodesic uses
`x_r=-u_future`; its background coordinate path has `t=-r/rho_b`.
The historical dimensionless convention is `q=r/rho_b`, so its
background path is `t=-q` and `h(x_q,x_q)=-rho_b^2`. Unit norm for
that q-velocity refers to the normalized metric `h/rho_b^2`.
Past proper duration is not the field equation's coordinate elapsed
time s. The retained relation `v(s)=J_E V(-s)` is a same-member tensor
time reflection; reconstructing original t introduces no extra reversal
of the physical outward radial normal.

An intrinsic shell timelike normal is distinct from the spacelike
outward bulk normal. Shell geodesics are not generally ambient free
fall; neither is a coordinate probe. A clock inverse changes the
parameter along the declared path and cannot change its acceleration
or confinement prescription. Proper-r geodesics already use their
proper clock and require no independent clock inversion unless an
observable is reexpressed on coordinate-time cuts.

This establishes actual common coordinate-probe clocks, their complete
two-parameter-order calculus and the finite moving-composition interface.
It adds no long intrinsic observer cylinder, shell-crossing transmission
theory, historical proper duration, effective amplitude, detector/noise
model or measured physical calibration.
