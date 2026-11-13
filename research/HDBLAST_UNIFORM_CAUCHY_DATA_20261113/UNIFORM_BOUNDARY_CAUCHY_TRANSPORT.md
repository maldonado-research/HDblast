# Uniform nonzero-state boundary transport of the initial Cauchy pair

13 November 2026 UTC. Uniform initial transport proof. This is an initial-section map
on a specified completed Cauchy-data ball. It is not an estimate on a
spacetime solution tube, a prescribed-time evolution theorem, or a new
amplitude admission. The actual canonical-to-harmonic completion must
independently place the whole selected curve in the ball stated below.

## 1. Exact input, output and domain

Use the fixed physical collar `0<=z<=6ell`, `ell=10^-6`, with
`dy=-rho dz`, `129<=rho<=131` and spatial background
`gamma0=dy^2+rho^2 kappa_S3`. The raw tensor connection is
`nabla_gamma0`, extended by `D dt=0`; its positive fiber metric is
`dt^2+gamma0` and its measure is `dvol_gamma0`.

Let `X=(g-g0,phi-phi0)|_(t=0)` and `Xdot=partial_t(g-g0,phi-phi0)|_(t=0)`
be the **fully relative-harmonic-completed** initial fields and velocities.
They are measured in physical spatial H14 and H13, respectively. No second
physical time derivative is assumed here. Define

```
R_j=||partial_e^j X||_H14(gamma0)
                  +||partial_e^j Xdot||_H13(gamma0), j=0,1,2.
```

For j=0 this means the background-subtracted values themselves. Parameter
derivatives hold the original chart and parameter fixed. All tensor slots,
including derivative slots, use the declared full Hilbert/Frobenius norms.

The sufficient completed initial ball for this lemma is

```
R_0<=10^-3228.                                           (1)
```

The domain can be taken open around this closed ball. This is a sufficient
data-norm criterion, not an amplitude derived by Taylor expansion using
the desired output Hessian. Actual whole-interval membership in (1) must
be proved directly from the canonical family and lapse/completion formulas.

The output is the actual straightened Cauchy pair
`C(X,Xdot)=(V(0),partial_t V(0))` in the fixed product bundle

```
E=R_v + Sym^2(R dt+T*S3)_H + R_w + (R dt+T*S3)_D,
q_sp=dz^2+kappa_S3, dmu=dz dmu_kappa.
```

Its norm is `||V||_H14+||V_t||_H13`. On the entire ball (1) the map
is C2 and the following bounds hold:

```
||C(X,Xdot)|| <=10^3653 R_0,
||partial_e C(X,Xdot)|| <=10^3661 R_1,
||partial_e^2 C(X,Xdot)||
                 <=10^3661 R_2+10^6884 R_1^2.            (2)
```

Equivalently its first Fréchet derivative has norm at most `10^3661`
and its bilinear second derivative has norm at most `10^6884`, in the
sum norm on the input Cauchy pair, throughout that ball. These are new
uniform nonzero-state bounds, not the old A/B coefficients evaluated only
at the zero background. No B_C for the actual family is assigned until
uniform actual R_1 and R_2 bounds have been inserted.

## 2. The full nonlinear point chart at the initial section

Let I be the fixed fiber isometry from the raw radially parallel coframe
`(dt,dy,rho theta)` to `(dt,dz,theta)`, acting on every tensor slot.
It is a fiber identification, not the coordinate pullback `dy=-rho dz`.
The generic linear estimate in
`HDBLAST_INITIAL_TRANSPORT_20261106/INITIAL_TRANSPORT.md` is

```
||I(T)||_Hm(q_sp) <=10^3000 ||T||_Hm(gamma0), m<=14.       (3)
```

That argument is a fixed-connection/fiber conversion for arbitrary tensors;
it is not restricted to zero parameter jets.

Write L(z) for its linear normalization into unrooted variables:

```
x=delta g_yy,
H_tt=delta g_tt/rho_b^2,
H_tA=rho delta g_(t,phys A)/rho_b^2,
H_AB=rho^2 delta g_(phys A,phys B)/rho_b^2,
w=delta phi,
D_t=delta g_yt/rho_b,
D_A=rho delta g_(y,phys A)/rho_b.
```

Here rho_b is the fixed background shell radius. The same retained source
proves `max_(j<=14)||partial_z^j L||op<10^210`, including symmetric
multiplicity and component norm conversions. The linear multiplier formula
`2^m sqrt(m+1)<10^6` gives

```
||L I T||_Hm(q_sp) <=A ||T||_Hm(gamma0),
A=10^3216,                        m<=14.                 (4)
```

Set `P=L I X` and `Q=L I Xdot`. For the unrooted state P define

```
r(P)=(sqrt(1+x)-1,H,w,D),
T(z,P)=F_-z(r(P)),
V=T(z,P).                                                (5)
```

This is precisely the proper-y radial component root `sqrt(g_yy)`.
It is neither the temporal ADM lapse nor an off-shell Schur-complement
lapse. F is the actual background-subtracted boundary-frozen flow of
`HDBLAST_BOUNDARY_FLOW_PATCH_20261106/BOUNDARY_FLOW_PROOF.md`. On the
audited collar the retained parameter is theta(z)=z. The global cutoff
outside that identity region is not differentiated here.

The target inner chart condition is satisfied on (1). Indeed (4) gives
`||P||H14<=10^-12`. The retained boundary-inclusive cylinder H3 embedding
has constant `E3=10^7` and therefore

```
sup|P|_E <=10^-5,  sup||P||_*<=10^-5.                    (6)
```

The radial root obeys `|sqrt(1+x)-1|<=|x|` for real x>-1.
Consequently `||r(P)||_*<=10^-5`. The actual flow certificate gives

```
||F_±z(r(P))||_* <=exp(400|z|)10^-5
                   <1.01*10^-5<10^-3,  |z|<=6ell.       (7)
```

Thus the proper-y root, flow, inverse flow and initial reconstructed metric
all remain in their actual inner chart patch. The derivative certificate
applies on the larger real input region `||u||_*<=.003`, so it has strict
room around every state used here. Equations (6)–(7) establish only the
initial state patch. They do not assert that an evolved solution stays in it.

## 3. Uniform mixed point coefficients through the required orders

The flow certificate proves for both signs of z, all `j+k<=18`,

```
||partial_z^j D_u^k F_±z||_* <10^126.
```

Since `||.||_*<=|.|_E<=sqrt(7)||.||_*`, a unit E direction is a
unit-or-smaller star direction and the output conversion costs less than
three. Hence the full E multilinear operator bound is below `10^127`.
No independent-entry infinity Jacobian is substituted for that operator
norm.

On the slightly larger unrooted real state patch `||P||_*<=10^-3`,
the derivatives of the scalar root through order eighteen obey

```
||D r||op<=1,
||D^k r||op<=17!<10^15,          2<=k<=18.                (8)
```

For example its kth scalar derivative is the falling factorial of 1/2
times `(1+x)^(1/2-k)`. The factorial is at most `(k-1)!/2` in absolute
value; `(1-10^-3)^(-35/2)<2` proves (8). Other blocks are linear and
have no higher derivatives.

Faà di Bruno for k labelled state directions has one term for every
set partition, at most `k^k<=18^18<10^23`. In each term there are at
most eighteen root derivative factors. Therefore, uniformly on that
patch and throughout `|z|<=6ell`,

```
||partial_z^j D_P^k T(z,P)||op <10^(127+23+18*15)
                              =M_T:=10^420,
                                              j+k<=18.  (9)
```

For k=0 the direct flow coefficient bound suffices and is smaller.
The explicit z dependence in (9) is differentiated at fixed P. The
angular maps in (5) are natural under orthogonal rotations of the sphere
tangent frame; the fixed background split, scalings and generator are
product-parallel in angular directions. Thus full product covariant
derivatives of the composed field have only the explicit z and state
chain branches charged below, with the induced connection on every slot.

## 4. Explicit Nemytskii bounds at H14/H13

Take any base P with `||P||H14<=10^-12`. For m=13 or 14 and
`1<=k<=3`, consider the kth state derivative of T applied to k arbitrary
Hm sections. At each of m spatial derivative steps, one may hit an
explicit z coefficient, create a derivative of the base P, or hit any
existing input/base derivative factor. At step j there are at most
`k+j+2` such branches. The full expansion consequently has at most

```
product_(j=0)^(m-1)(k+j+2)
 =(m+k+1)!/(k+1)! <=18!                                  (10)
```

terms. All derivative indices and permutations are included. Each
coefficient in a term is a mixed derivative from (9), of total order
at most m+k<=17. No derivative of an evolving spacetime solution occurs.

The sum of derivative orders assigned to all field factors is at most m.
At most one factor can have derivative order greater than m-3; two such
orders already exceed m for m=13,14. Put that factor in L2 and all
others in L-infinity using the proven tensor H3 embedding. If there is
no such factor, choose any input factor for L2. If there are n factors
from the base P, the resulting estimate for a term is at most

```
M_T E3^(k+n-1) ||P||Hm^n product_(i=1)^k ||A_i||Hm
 <=M_T E3^(k-1) product_i ||A_i||Hm,                      (11)
```

because `E3||P||Hm<=10^-5<1`. Only field derivative tensors of rank
at most thirteen are put into the embedding; high-rank coefficient maps
are controlled by (9), not passed to a tensor Sobolev embedding.

Sum the m+1 output derivative squares, using sqrt(m+1)<4. Equations
(10)–(11) yield

```
||D^k T(P)[A_1,...,A_k]||Hm
 <=4*18!*10^420*10^(7(k-1)) product_i ||A_i||Hm.
```

Since `4*18!<10^17`, the following common explicit ceilings hold in
both target orders:

```
C1=10^437,   C2=10^444,   C3=10^451.                    (12)
```

The notation D^k T in this display denotes the kth **state** derivative
of the point map, evaluated on fields; spatial derivatives use Hm.
For H13 output, arguments in H14 are simply bounded in H13. In
particular (12) applies when precisely one argument is a time velocity
having only H13 regularity. No H14 bound for that velocity is assumed.

These bounds also give a C2 function-space map on an open neighborhood
of the displayed closed base ball. For example take `||P||H14<2*10^-12`;
the estimates still use `E3||P||<1` and remain inside (9). For
continuity of the third state derivative in H13, the same argument with
m=13,k=4 uses at most 18! terms and total coefficient order seventeen;
it gives the finite bound `10^458`. Integral mean-value comparison
then supplies operator-norm continuity. The H14 second derivative uses
the already bounded third derivative. This invokes no third derivative
of the actual parameter curve.

## 5. Initial time velocity at arbitrary nonzero state

There is an exact even-frame identity at t=0, not an assumption that
the perturbed solution is time symmetric. Extend the raw fiber frame
`(dt,dy,rho theta)` independently of t. In the fixed de Sitter coframe
`(dt,cosh(t)theta)`, the time-dependent linear normalization is L_t(z):
its H_tA and D_A factors acquire sech(t), and its H_AB factor acquires
sech(t)^2. Other factors are unchanged. The subtracted background
metric is even in t and the background scalar is independent of t.

Thus the explicit derivative of this normalization at t=0 vanishes
for **every raw state**:

```
partial_t L_t(z)|0=0.
```

The radial component root is time independent. Expressed in the same
de Sitter coframe, the actual autonomous generator uses J4 and fixed
shell coefficients, so the flow has no additional explicit time
dependence. Therefore the complete initial velocity is exactly

```
V=T(P),             V_t=D T(P)[Q],
P=L I X,            Q=L I Xdot.                           (13)
```

This identity is valid at nonzero parameter. The original temporal lapse
remains the prescribed `alpha_e=rho(1+e u_profile Y)`, with zero initial
shift and the actual harmonic-completed time velocities supplied by the
preceding Cauchy map. It is not changed into the radial square root in
(5). No second physical time derivative of g or phi is needed to calculate
the chart transformation of the initial velocity.

For any C2 parameter curve of these completed data, primes denoting
parameter derivatives, the exact identities are

```
V'=T1[P'],
V''=T1[P'']+T2[P',P'],
(V_t)'=T1[Q']+T2[P',Q],
(V_t)''=T1[Q'']+2T2[P',Q']+T2[P'',Q]
                                      +T3[P',P',Q].      (14)
```

Here Tj is the jth state derivative evaluated at the same nonzero P.
The last two terms in (14) vanish at the retained zero-state background
because Q=0; they must be retained at a general nonzero parameter.
Omitting them would reuse the old zero-state formula incorrectly.
The largest needed coefficient counts are 14+2 for V'' and 13+3 for
(V_t)'', both within the flow's total eighteen-jet budget. No full
spacetime H14 norm is inferred from the spatial Cauchy-pair norm.

## 6. Uniform linear and quadratic coefficients

Set `p_j=||P^(j)||H14`, `q_j=||Q^(j)||H13`. Equation (4) gives
`p_j+q_j<=A R_j`, and (1) gives `p_0<=10^-12`, `q_0<=1`.
Since T(z,0)=0, integrate its first state derivative along the segment
from zero to P. Equations (12)–(14) imply

```
||V||H14+||V_t||H13 <=C1(p_0+q_0),
||V'||H14+||(V_t)'||H13 <=(C1+C2 q_0)(p_1+q_1),
||V''||H14+||(V_t)''||H13
 <=(C1+C2 q_0)(p_2+q_2)
                  +(C2+C3 q_0)(p_1+q_1)^2.              (15)
```

The two mixed quadratic terms are charged by
`p_1^2+2p_1q_1<=(p_1+q_1)^2`. With q_0<=1,

```
C1+C2<10^445,      C2+C3<10^452,
A C1=10^3653,
A*10^445=10^3661,
A^2*10^452=10^6884.
```

This proves (2). The same calculation polarized in two Cauchy-pair
directions proves the bilinear second-derivative bound stated there.
All constants are explicit consequences of the retained flow jets,
the actual root, fixed linear normalizations, the cylinder embedding
and finite branch counts. None is an unexplained uniform composition
constant.

## 7. Support, copies, and the precise remaining join

The actual selected canonical curve and its completed lapse/shift fields
equal the background on fixed open shell and outer collars. The completion
is local in the canonical fields and their spatial derivatives. Hence the
completed raw perturbations and velocities have support in the retained
annulus, lying wholly inside theta=z. Since T(z,0)=0, (13) preserves
their support and the original physical anchor. No derivative of a
global cutoff outside the audited collar enters the estimates.

The same bounds apply to the specified reflected two-copy family using
its positive inward coordinate on each copy and summed Hilbert norms.
Linear estimates sum in squares, and the nonlinear product sum uses
`sum_i a_i^2 b_i^2 <=(sum_i a_i^2)(sum_i b_i^2)`; no new negative-z
background extension is assumed.

To obtain a numerical B_C for the actual whole parameter interval, the
canonical-to-physical/relative-harmonic completion must independently
prove (1), and provide uniform finite bounds B_raw,1 and B_raw,2 for
R_1 and R_2. Then the concrete join is

```
B_C=10^3661 B_raw,2+10^6884 B_raw,1^2.                   (16)
```

The admission of (1) should use direct whole-ball canonical size and the
actual lapse formula, not a Taylor estimate containing this not-yet-joined
B_C. Equation (16) is an initial Cauchy Hessian statement only. A uniform
mixed spacetime solution norm, prescribed-slab continuation and proper
observer tube require their separate estimates.
