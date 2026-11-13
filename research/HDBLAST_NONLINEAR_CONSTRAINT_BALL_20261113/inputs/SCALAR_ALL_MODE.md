# Scalar all-mode weighted inverse: a complete conditional quantitative argument

13 November 2026 UTC. Private derivation with independent review. The
background, canonical constraint normalization, selected weighted domain
and qualitative projected inverse are retained analytic premises. This
note evaluates an a priori constant; it does not supply a nonlinear
parameter neighborhood or an evolution theorem.

## 1. Precise operator, weights, and conditional conclusion

Write `gamma0=dy^2+rho(y)^2 kappa_S3`, with the original fixed volume
`dmu0=rho^3 |dy| dmu_kappa`. Probability or full geometric sphere measure
can be used consistently: a common multiplicative volume factor does not
alter an operator-norm inequality. Keep

```
ell=10^-6, kappa0=10^-8,
x=(z-3ell/2)/ell, r=1-x^2,
q=kappa0^2 r^2, w=exp(-1/r), dy=-rho dz,
A={ell/2<z<5ell/2} x S3.
```

The symbol `q` is the scaled derivative length, not its square root.
For scalar functions the retained norm is exactly

```
||N||_(H^m_qw)^2 = sum_(j=0)^m integral_A
                     w^2 q^(2j) |nabla_gamma0^j N|^2 dmu0.       (1)
```

The dotted space is the completion of compactly supported smooth
functions in (1). The reflection group is the retained coordinate group
`Gamma=(Z2)^4`. Its scalar invariant harmonics have degrees 0,2,4,...;
in particular the degree-one boost-lapse kernel is excluded.

With the canonical scalar normalization, put

```
P_g h=(div div h-Delta tr h-Ric:h+(dphi0 tensor dphi0):h)/2,
P_g^* N=(Hess N-gamma0 Delta N-N Ric+N dphi0 tensor dphi0)/2,
L_H N=w^-2 P_g(w^2 q^4 P_g^* N).                            (2)
```

The source of `L_H N=F` is the *normalized* source `F`. If the physical
constraint source is `f`, then `F=w^-2 f`. It is not legitimate to
substitute the unnormalized physical source into the range norm below.

The conditional estimate is

```
||L_H^-1 F||_(dot H16_qw) < 10^5000 ||F||_(dot H12_qw)      (3)
```

for every nonzero Gamma-invariant source, with `<=` for the zero source.
Equivalently one can use `<=10^5000` for all sources. The selected inverse
and its existence are the original qualitative ones. The proof below is
an independent all-mode a priori estimate and does not insert the
qualitative inverse norm as a numerical constant.

The complete argument uses the two-sided norm bridge stated in section 7.
That bridge is proved and independently reviewed in the companion private
`WEIGHTED_NORM_BRIDGE.md`; its constant is below `10^29`. The deliberately
loose budget used here is `10^1000` in each direction through order 16.

## 2. An all-mode scalar coercivity bound

All following estimates first concern a compact smooth radial coefficient
of a unit scalar harmonic `E`, `-Delta_kappa E=lambda E`. Set

```
N=n(y)E, Q=P_g^*N,
T=||N||_w, E_tail=||r^4 Q||_w,
I=||N||_(w,{r>=1/5}).
```

The monopole estimate already proved in the retained scalar coercivity
note is

```
T <=1.3*10^24 E_tail,                       lambda=0.     (4)
```

Its proof keeps the actual scalar-gradient lower bound
`|p|>=.00012`; this input is not removed or replaced by a vacuum estimate.

For every `lambda>=8`, the angular trace-free block is
`n[Hess E+(lambda/3)E kappa]/(2rho^2)`. Its squared angular norm is
`lambda(lambda-3)|n|^2/(6rho^4)`. Consequently, writing `|Q|_Omega`
for the complete angular L2 tensor norm at fixed y,

```
|n| <=6650 |Q|_Omega,
lambda |n| <=53200 |Q|_Omega.                            (5)
```

Indeed `129<=rho<=131`, `lambda(lambda-3)>=40`, and
`lambda/(lambda-3)<=8/5`. These inequalities are uniform in lambda.
The second is the estimate missing from a direct reuse of the old
three-mode tail proof.

On either endpoint tail `r<=1/4`, use `t=1/r>=4`, `k=t_y` and the exact
weighted radial measure `W(t)dt=exp(-2t)rho^3|dy/dt|dt`. The retained
calculation gives

```
-W_t/W>=2, |alpha|<.55, |k^-1|<4.73*10^-6,
k^-2 <=6*10^-9 r^4,
alpha=2/t+1/[2t^2(1-1/t)]+a/k,
n_tt+alpha n_t=(-2B+C_lambda n)/k^2,
C_lambda=2lambda/(3rho^2)-R_t,
R_t=3/rho^2+p^2/3-3a^2, |R_t|<.04.                     (6)
```

Here `B` is the angular-trace coefficient of Q. In particular
`2|B|<=2|Q|_Omega`. Unlike the previous finite-mode argument, no
lambda-independent bound is asserted for `C_lambda/k^2`. Instead (5)
gives

```
|C_lambda n|
 <=[2*53200/(3*129^2)+.04*6650]|Q|_Omega
 <270 |Q|_Omega.
```

Thus the right side of (6) is bounded by
`2*10^-6 r^4 |Q|_Omega`.

Take the retained cutoff `chi(t)`, zero for `t<=4`, one for `t>=5`,
with `|chi'|<=2`, `|chi''|<=6`. For `f=chi n`, weighted integration
by parts gives `||f||_W<=||f_t||_W<=||f_tt||_W`. Therefore

```
||chi n||_W
 <=3[2*10^-6 E_tail
       +||2chi' n_t||_W+7.1 I].                         (7)
```

The radial adjoint component gives the exact identity

```
n_y=a n+(lambda-3)n/(3a rho^2)-2A_rad/(3a).
```

Using (5), `.111<=a<=.112`, and `rho>=129`, its absolute value is below
`800 |Q|_Omega`. On the cutoff transition this implies
`|n_t|<.004 |Q|_Omega<=2.5 r^4|Q|_Omega`. Hence (7) yields, on each tail,

```
||chi n||_W <=32 E_tail+22 I.
```

The residual cutoff is supported in `r>=1/5`. Adding the two tails and
that residual gives

```
T <=64 E_tail+45 I,
I <=625*6650 E_tail,
T <=2*10^8 E_tail,                             lambda>=8. (8)
```

The largest constant is therefore still the monopole constant (4).
Orthogonality of the scalar harmonic representations and their natural
tensor derivatives permits square summation of these estimates. Since
`q^2=kappa0^4 r^4`, every Gamma-invariant compact smooth potential obeys

```
||N||_w <=C_gap ||q^2 P_g^*N||_w,
C_gap=1.3*10^56.                                          (9)
```

Finite harmonic truncation and compact approximation extend (9) to the
retained dotted domain. The degree-one kernel is not silently included.

Pairing (2) with N in the weighted measure gives

```
||q^2P_g^*N||_w^2=<N,L_HN>_w.
```

Together with (9), this proves the all-mode zero-order inverse estimate

```
||N||_w <=C_gap^2 ||L_H N||_w <=10^113 ||L_H N||_w         (10)
```

for every N in the domain. This step does not use a particular compact correction
or finite angular support of the actual second source.

## 3. Exact radial unitary operator for arbitrary lambda

Set

```
ds/dy=1/q, M=w sqrt(q rho^3), v=M n, g=M F_lambda,
D=d/ds, b=1/rho, h=q b, zeta=h sqrt(lambda),
d=q_y, c=D log M, eta_h=D log h=d-qa.
```

The coordinate s covers the complete real line. The transform is unitary
from the weighted radial L2 measure to `L2(ds)`. In particular
`g=sqrt(q rho^3)w^-1 f_lambda` for an unnormalized physical source f.

The exact transformed adjoint rows, valid for arbitrary scalar eigenvalue,
are

```
Arow=zeta^2/2+alpha,
Brow=-D^2/2+zeta^2/3+beta,
sqrt(lambda) Crow=zeta(D-e),

alpha=alpha1 D+alpha0,
alpha1=-3aq/2,
alpha0=3aqc/2+3q^2(a^2-b^2)/2,

beta=beta1 D+beta0,
beta1=c+d/2-aq,
beta0=-(c^2+dc-Dc)/2+aqc-q^2 R_t/2,
e=c+aq.
```

Consequently, in flat ds formal-adjoint pairing,

```
L_lambda=Arow* Arow+3 Brow* Brow
          +(1/2)[zeta(D-e)]*[zeta(D-e)]
          +zeta^4/6-h^2 zeta^2/2.                       (11)
```

This is also correct at lambda zero, where all zeta terms vanish.
Let `m=zeta^2`, `P=-D^2+m` and `Aop=1+P`. Direct noncommutative
composition gives

```
L_lambda=(3/4)P^2+R,

R=m''/4+(m alpha+alpha* m)/2+alpha*alpha
  +3[(-D^2/2+m/3)beta+beta*(-D^2/2+m/3)+beta*beta]
  +[(m e)'+e^2m-h^2m]/2.                                (12)
```

Every term of R has differential order at most three when both D and
multiplication by zeta are counted as order one. No bounded-lambda
assumption is used. The `m''/4` and `-h^2m/2` terms are essential.

The companion replay independently composes the generic differential
operators with arbitrary smooth real coefficient functions and checks
all five coefficients in (12). That replay validates the identity, not
the analytic norm estimates which follow.

## 4. Uniform coefficient ledger

The retained scalar regularity note supplies a convenient global polynomial
coordinate description on `|x|,|b|,|a|,|p|,|eta_phi|<=1`:

```
D x=-10^-10(1-x^2)^2 b,
D b=-q a b,
D a=q(-b^2-p^2/3),
D p=q[U'(eta_phi)-4ap],
D eta_phi=q p.
```

Each vector-field component has polynomial coefficient l1 norm at most
`10^-8` and degree at most 9. For a polynomial of degree d0 and
coefficient norm C0,

```
||D^j polynomial||_coeff
 <=C0 product_(i=0)^(j-1) [10^-8(d0+8i)].                (13)
```

This estimate uses the exact retained polynomial potential and the actual
background ODE, not independently assumed high ordinary y derivatives.

Here `q=kappa0^2(1-x^2)^2`,
`d=4*10^-10 x(1-x^2)b`, and
`c=d/2+3qa/2+2*10^-10 x b`. Thus `c,d,qa,eta_h` have coefficient norm
below `10^-8` and degree at most 5. For integers `0<=j,beta_power<=32`,

```
D^j(zeta^beta_power)=zeta^beta_power p_(j,beta_power),
||p_(j,beta_power)||_coeff <=(4*10^-6)^j.                (14)
```

Indeed `p_0=1`, `p_(j+1)=D p_j+beta_power eta_h p_j`, the degree is at
most `8j`, and every induction multiplier is at most
`10^-8(8j+beta_power)<4*10^-6`. The zero-power/zero-mode cases are read
directly as derivatives of a constant or of zero. In particular

```
|m''| <=16*10^-12 m <10^-10 m.                          (15)
```

This is a *relative* bound. `m''` is not asserted to be absolutely bounded
independently of lambda.

For completeness, a deliberately generous ledger for (12) is as follows.
The coefficient polynomials of alpha have degree at most 10, those of
beta at most 13, and their coefficient norms are at most one. Normal
ordering the products in (12), using
`D zeta^b=zeta^b(D+b eta_h)`, requires at most two derivatives of those
rows. Every resulting coefficient has degree at most 50. The sum of its
coefficient l1 norms over all monomials `zeta^b D^a`, `a+b<=3`, is below
`10^6`. Here is an explicit generous product ledger. The coefficient sums
of alpha and beta are at most two, those of their formal adjoints at
most three, and that of `B0=-D^2/2+m/3` is below one. Normalized
coefficient derivatives through the two composition steps have no larger
norm, by (13)--(14); Leibniz multiplication therefore costs at most
`2^order`, with order at most two. The cross term involving alpha costs
at most `(2+3*2)/2=4`, and `alpha*alpha` costs at most `3*2*2=12`.
The three beta products cost at most `3*(8+6+12)=78`. The `m''/4`
term costs at most `1/4`; the final three multiplication terms cost at
most `5/2`. Thus the fewer than twelve displayed summands cost below
1000 in total. The ceiling `10^6` leaves a further factor 1000.

If `C zeta^b D^a` is in normal form, left multiplication by Aop gives
five groups of coefficients, with total weight at most
`1+1+2+1+1=6`, provided the first two derivatives of `C zeta^b`, after
factoring out `zeta^b`, have coefficient norm at most `||C||_coeff`.
Equation (13), (14), and `degree(C)<200`, `b<=16` ensure exactly that.
Each multiplication increases polynomial degree by at most 16. Hence

```
||Aop^6 g||_2 <=6^6 J12(g)<=10^5 J12(g),
||Aop^6 Rv||_2 <=6^6*10^6 J15(v)<=10^11 J15(v),          (16)
```

where `J_m` is defined next. The maximal coefficient degree in the second
line is at most `50+6*16=146<200`; the highest weighted order is 15.
These estimates remain uniform as lambda tends to infinity.

## 5. A scalar graph-norm lemma with explicit constants

For each fixed finite lambda define

```
J_m(v)^2=sum_(a+b<=m) ||zeta^b D^a v||_2^2.
```

The operator `Aop=1-D^2+zeta^2` is the positive self-adjoint operator on
`L2(R)` associated with its closed quadratic form. For fixed lambda its
smooth potential is bounded, so compact smooth functions lie in the
domains of all powers needed here. It obeys `Aop>=1`.

The following uniform estimate holds for `0<=m<=16`:

```
J_m(v) <=10^(10m) ||Aop^(m/2)v||_2.                     (17)
```

The case m=0 is equality and
`J1(v)^2=<v,Aop v>`. For the base second-order estimate, integration by
parts gives

```
||Aop v||^2=||v''||^2+||(1+m)v||^2
             +2 integral (1+m)|v'|^2 - integral m''|v|^2.
```

Here m means `zeta^2` in this displayed identity. By (15) the last term
is absorbed in `2||zeta v||^2`, already present in the preceding square.
It follows directly that `J2(v)<=2||Aop v||`.

Here is an explicit induction, avoiding an unquantified elliptic estimate.
For `a+b<=j<=14`, let `T_ab=zeta^b D^a`. Then

```
[Aop,T_ab]
 =-(zeta^b)'' D^a-2(zeta^b)'D^(a+1)
   -zeta^b sum_(k=1)^a binom(a,k) m^(k) D^(a-k).         (18)
```

By (14), this has weighted order at most j+1 and sum of coefficient
norms below `2^j+3<10^5`. Thus
`||[Aop,T_ab]v||<=10^5 J_(j+1)(v)`.

To recover every raw monomial of order at most j+2 from the collection
`J2(T_ab v)`, choose an outer monomial `zeta^d D^c` of order at most two,
and leave an inner `T_ab` of order at most j. Normal ordering produces
only lower terms of order at most j+1, with coefficient sum at most three.
There are at most six outer monomials per inner monomial, and at most
153 target monomials. Consequently

```
J_(j+2)(v)
 <=6 sum_(a+b<=j) J2(T_ab v)+459 J_(j+1)(v)
 <=10^9 [J_j(Aop v)+J_(j+1)(v)].                        (19)
```

For the last inequality there are at most 120 inner monomials; use the
base estimate, (18), and Cauchy--Schwarz for their Aop v terms.
If (17) holds at j and j+1, the right side of (19) is at most

```
10^9[10^(10j)+10^(10j+10)] ||Aop^((j+2)/2)v||
 <=10^(10j+20) ||Aop^((j+2)/2)v||,
```

where `Aop>=1` compares successive fractional powers. This proves (17).

## 6. Uniform fourth-order regularity through the required order

For a compact smooth radial v, write `g=L_lambda v`. Applying Aop^6 to
(12), and using (16), (17), gives

```
(3/4)||Aop^8 v||
 <=10^5 J12(g)+(3/2)||Aop^7v||+(3/4)||Aop^6v||
                     +10^11 J15(v)
 <=10^5 J12(g)+10^162 ||Aop^(15/2)v||.                  (20)
```

The spectral theorem and a scalar split at `s^(1/2)=delta^-1` give

```
||Aop^(15/2)v|| <=delta ||Aop^8v||+delta^-15||v||.
```

Take `delta=10^-164`. In (20) the coefficient absorbed on the left is
`.01`, leaving at least `.74`. Thus

```
||Aop^8v|| <=10^2624 [J12(g)+||v||],
J16(v) <=10^2800 [J12(g)+||v||].                        (21)
```

This estimate includes arbitrarily large angular eigenvalues; they were
part of Aop, rather than bounded coefficients in a radial ODE.

## 7. Bridge to the literal retained tensor-derivative norm

The norm conversion must hold for the source as well as the potential.
The required companion statement, with very generous common budget, is

```
10^-1000 J_m(M n) <= ||n E||_(H^m_qw)
                    <=10^1000 J_m(M n),  0<=m<=16,     (22)
```

uniformly over scalar harmonic degree and multiplicity. All quantities
are integrated in their respective unitary measures. The same statement
holds after square summing harmonic components.

For clarity, the mechanism of the bridge is recorded here. Put
`W_j=M q^j nabla_gamma0^j(nE)` in the product radial/angular orthonormal
frame. The exact radial and angular recursions are

```
(W_(j+1))_rad=(D-c-jd)W_j,
(W_(j+1))_ang=h nabla_kappa W_j+qa S_j W_j,              (23)
```

where the invariant tensor-slot insertion has collected norm at most
`4j`. The reverse formulas solve (23) for `D W_j` and
`h nabla_kappa W_j`. Normal ordering charges derivatives of h through
`D log h=eta_h`; derivatives of c,d,qa have the polynomial bounds in
section 4. Thus every stage is a finite triangular conversion including
all derivative slots, rather than a componentwise angular restriction.

For a scalar harmonic, natural SO(4) lifts give
`||nabla_kappa^b E||<=lambda^(b/2)||E||`. Conversely, adjacent metric
traces of the ordered covariant derivative, with one gradient for odd b,
give `lambda^(b/2)||E||<=3^(b/4)||nabla_kappa^b E||`.
No interchange of derivative order is required. Distinct scalar harmonic
degrees remain orthogonal after these natural tensor derivatives, since
the derivatives intertwine the rotation action and its self-adjoint
Casimir. This justifies the all-mode square sums in (22).

The companion bridge proof evaluates the stage factors in (23) and the
number of derivative slots, obtaining the common constant
`13*8^16*16!<10^29`. The `10^1000` used in (22) is vastly
larger than this factorial bound. It is not inferred from the existing
finite-mode conversion `10^2000`, which was a different calculation.

Apply (22) to v and to `g=M F_lambda`, sum (21) in square over angular
modes, and use the unitary identity for v. This gives

```
||N||_(H16_qw)
 <=10^4800 [||L_HN||_(H12_qw)+||N||_w].                 (24)
```

Initially this is proved for compact smooth potentials with finite
harmonic support. Harmonic truncations commute with the operator and
converge in the retained norms. Compact approximation then extends (24)
to the dotted domain: L_H is a bounded differential operator from
H16_qw to H12_qw, as follows from the same normal forms and the reverse
bridge. No unweighted trace or endpoint boundary condition is imposed
on the potential.

Finally use (10). The factor in (24) is bounded by
`10^4800(1+10^113)<10^4914<10^5000`. This proves the conditional
estimate (3) for the retained selected inverse. Existence and uniqueness
are still those of the stated weighted inverse premise; the numerical
norm follows from the separate inequalities above.

## 8. Premise ledger and source pins

The premises retained, rather than reproved here, are: the actual smooth
Einstein--scalar background and polynomial potential on the annulus;
the fixed canonical constraint normalization; the original dotted weighted
domain and Gamma-invariant qualitative isomorphism; and the already
proved monopole identity and numerical bound (4). The extension from
three scalar modes to all allowed scalar modes, the quantitative graph
regularity estimate, and the explicit source normalization are the new
parts of this derivation. The companion norm bridge and this scalar
argument have each received independent internal mathematical review.

The accepted companion file is
`/workspace/hdblast-all-mode-20261113/WEIGHTED_NORM_BRIDGE.md`, SHA256
`17bbd84869346f90e8f02e7b343c49fa155d9f22138e23ff366d56769551cd63`.

All paths below are under
`research/ in the pinned base tree`.

| Source | SHA256 |
|---|---|
| `HDBLAST_WEIGHTED_REGULARITY_20261106/README.md` | `595b8d46ad8db6e198dde730061eb5c612150b659cd266202353fca8cc2fc09e` |
| `HDBLAST_WEIGHTED_REGULARITY_20261106/SCALAR_COERCIVITY.md` | `c8bb6264c981c7e4c4b5a4e580552c38dc1093fecc8e8474502650888775d942` |
| `HDBLAST_WEIGHTED_REGULARITY_20261106/SCALAR_REGULARITY.md` | `e830a68ceff19e04f03d6ae4d53f7342526c54bb5f1e7a82d1cef8d542f02816` |
| `HDBLAST_WEIGHTED_REGULARITY_20261106/inputs/PRIOR_SCALAR_OUTPUT.md` | `67c3fe77b0fc6c2d0ce07e7466a46148e7983f93c4f14dae0344f00a2c6dc925` |
| `HDBLAST_EXACT_LOCALIZED_PRECONTACT_20261014/EXACT_LOCALIZED_INITIAL_DATA_AND_PRECONTACT.md` | `95ccaa550309b72226f8d03d3a950e4e10760d0d98810f0e818d19057a8cfd0b` |
| `HDBLAST_WEIGHTED_REGULARITY_20261106/inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json` | `fd93a0de3d321ff5182115138eed2725f3c327083cf5148e3efc1aca3659d5c0` |

No finite-amplitude implicit-function radius follows from (3). In
particular this note does not evaluate the weighted nonlinear derivative
defect, uniform initial-family gauge completion, or the nonzero-state
quasilinear evolution constants recorded in the private next-step plan.
