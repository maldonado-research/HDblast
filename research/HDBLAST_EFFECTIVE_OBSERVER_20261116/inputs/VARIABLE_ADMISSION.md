# A numerical admission bound for compatible smooth Picard iterates

14 November 2026 UTC. Joined mathematical calculation, conditional on the
companion energy and normal-recovery bounds at the exact constants below.
This establishes a uniform smooth-iteration bound. It does not establish
the finite-Sobolev nonlinear limit, top-norm continuity, evolved parameter
regularity or an observer tube.

## 1. Inputs and source spaces

Retain the actual normalized equation, same-state elapsed reflection,
`Sigma=[0,6ell] x S3`, `ell=10^-6`, `T=13ell/8`, the orthogonal mixed rows,
and the fixed positive product-bundle norms. Use

```
X_m(w;s)=sum_(j=0)^m ||partial_s^j w(s)||H^(m-j),
I14=sum_(j=0)^14 ||J_j||H^(14-j),
R=10^-12.
```

The companion estimates, for any smooth predecessor v with
`sup_s X14(v)<=R` and arbitrary smooth successor u solving the frozen
equation with the common initial jets, are

```
Y14(u;s)<=C_Y I14+integral_0^s[K_u X14(u;r)+K_v X14(v;r)]dr,
X14(u;s)<=N Y14(u;s)+L_u X13(u;s)+L_v X13(v;s),
C_Y=10^30, K_u=10^500101, K_v=10^500100,
N=10^110, L_u=L_v=10^500208.                              (1)
```

Only the predecessor is assumed to lie in the ball. No successor bound
is used to establish (1). The coefficients lie in their actual state
patch by the already evaluated embedding. The smooth frozen linear
construction on this compact coefficient interval is the retained
construction in `QUASILINEAR_VARIABLE_MIXED_EVOLUTION.md`, Section3;
no qualitative energy constant from that construction enters (1).

For a common collar-supported initial Cauchy pair of size delta, the
accepted nonlinear initial recursion and auxiliary heat lift give

```
I14<=C_J delta,              C_J=10^6503000,
sup_s X14(U^0;s)<=C_seed delta,
C_seed=10^6503255.                                      (2)
```

Here C_seed acts on the Cauchy-pair size delta. The smaller C_lift=10^255
acts on I14, and cannot replace C_seed in (2). The initial-jet guard is
`delta<=10^-6503020`. The selected original canonical family obeys
`d(a)<=10^4200 a` on its previous interval; no parameter is rescaled.

For smooth data, U^0 is smooth. Its nonlinear jets agree through J14 and
it vanishes on fixed open endpoint collars for all auxiliary times.
Smooth compatible approximants first smooth the Cauchy pair, recompute
its nonlinear jets, and then lift. Independently truncating all jets
does not preserve the nonlinear recursion.

## 2. Common finite jets and lower-norm integration

Freeze the coefficients and lower source at U^n and solve for U^(n+1)
with the same original Cauchy pair. The triangular recursion for the
frozen equation through temporal row twelve uses only predecessor jets
through thirteen. If U^n has the exact nonlinear jets J0 through J14,
the frozen recursion therefore gives the same jets for U^(n+1).
The auxiliary starting field supplies the induction base.

The homogeneous rows and higher smooth corner compatibility also hold.
The starting field has all-order zero initial boundary collars. For each
subsequent frozen equation, its actual higher initial jets are computed
from the smooth predecessor and compactly supported Cauchy pair; they
vanish near the boundary by induction. Higher *interior* jets of the
starting field are not asserted to equal those of the nonlinear PDE.
Only J0 through J14 are used in (1)--(2).

For any such iterate w, including U^0, the fundamental theorem of calculus
in each H^(13-j), `0<=j<=13`, gives

```
X13(w;s)<=I14+integral_0^s X14(w;r)dr.                    (3)
```

The seed proof supplies the required strong derivatives in these spaces;
smooth successors satisfy them directly. This is a relation between
actual temporal derivatives, without interpolation between unrelated
time jets.

## 3. Evaluated Volterra majorant

Insert (3) into (1). For x_n(s)=X14(U^n;s),

```
x_(n+1)(s)<=D0 delta+integral_0^s[H_u x_(n+1)+H_v x_n]dr,
D0=(N C_Y+L_u+L_v)C_J,
H_u=N K_u+L_u,                 H_v=N K_v+L_v.
```

The displayed numerical constants imply

```
D0=(10^140+2*10^500208)*10^6503000 <10^7003209,
C_seed<10^7003209,
H_u=10^500211+10^500208<10^500212,
H_v=10^500210+10^500208<10^500212.
```

Use the safe rounded constants

```
C0=10^7003210,                 H=10^500212.              (4)
```

If x_n(s)<=C0 delta exp(2Hs), the same function is a supersolution of
the scalar Volterra inequality for x_(n+1), since its initial value is
at least D0 delta and its derivative is twice H times itself. Scalar
integral comparison, or Gronwall applied to the positive difference,
gives

```
x_(n+1)(s)<=C0 delta exp(2Hs).                            (5)
```

Equation (2) supplies the base. This induction may be used only when
the resulting majorant is inside R, so that every frozen predecessor
meets the hypothesis of (1). The next section verifies that admission;
the ball condition is not assumed as the conclusion of the induction.

## 4. An explicit smaller interval in the original parameter

Define the positive exact number

```
a_iter=10^(-10^500209).                                  (6)
```

This notation denotes a nested power, not zero and not an expanded
decimal or a measured cosmological amplitude. It is strictly smaller
than the preceding seed interval `10^-6507485`, because
`10^500209>=10^7>6507485`.

Take `|epsilon|<=a_iter` and allow a factor two for supported smooth
approximants, so

```
delta<=2*10^4200*a_iter.                                 (7)
```

This is strictly inside the initial-jet guard. For example
`2*10^(4200-10^7)<10^-6503020` is already a much weaker sufficient
comparison. The fixed zero endpoint collars and symmetry constraints of the initial
lift are retained. Its permitted slightly enlarged smooth support annulus
inside the cutoff plateau may be used for approximants; exact preservation
of a nonsmooth support edge is not assumed. The estimate covers one copy and the
retained identical-two-copy Hilbert union when delta and I14 already
use the global union norms.

The exponential is bounded without evaluating an enormous decimal:

```
2HT=(13/4)*10^500206<10^500207,
e<10  =>  exp(2HT)<10^(10^500207).
```

Combining (4), (6), (7), and this inequality,

```
C0 delta exp(2HT)
 <2*10^(7007410-10^500209+10^500207)
 =2*10^(7007410-99*10^500207)
 <2*10^-15<R/4.                                        (8)
```

The penultimate step needs only `99*10^500207>7007425`, a finite
exponent comparison. In particular it follows from
`500207>=6` and `99*10^6=99000000>7007425`. No integer with hundreds
of thousands of digits is expanded in the replay.

Equations (5)--(8) close the smooth-iterate ball induction on the
entire prescribed elapsed slab. Each predecessor has X14<R/4 and lies
strictly in the actual point-state patch. Linear frozen existence is
used successively under its verified coefficient and corner hypotheses;
this does not invoke an unevaluated nonlinear precontact lifetime.

## 5. Exact scope and remaining work

The conclusion is a uniform X14 bound for the compatible smooth Picard
iterates and supported smooth approximants under the retained analytic
premises. It supplies evaluated numbers for the previously conditional
majorant. Its very small sufficient interval is a conservative
mathematical threshold in the registered normalization.

An actual finite-Sobolev nonlinear solution is not concluded here.
The lower-energy consecutive-iterate difference estimate and its
convergence/uniqueness argument remain to be evaluated. Passing to a
limit with continuous full X14 still requires the separate order-fifteen
smoothing/tame interface or a proved substitute. The current starting
field reproduces fifteen initial jets, not all higher nonlinear interior
jets. Initial/seed C2 regularity does not prove evolved C2 dependence;
the retained conservative parameter scheme targets actual C2 in X10,
not C2 in X14.

The artificial seed support extends to11ell/4, beyond the true initial
support5ell/2. The outer causal marginell/4 belongs to an eventual
solution only after finite propagation is established. Homogeneous outer
Dirichlet flux handles the iterates independently of that support
argument. No observer tube, finite nonlinear observable remainder,
physical energy scale or empirical higher-dimensional event is proved.
