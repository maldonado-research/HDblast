# Full covariant momentum H13 and weighted potential H14

6 November 2026. This proof applies to the newly specified member of
[the coercivity calculation](MOMENTUM_COERCIVITY.md), with the
conditional selected solution estimated in
[the complete-line regularity proof](MOMENTUM_REGULARITY.md).
It treats all covariant derivatives, including angular derivatives and
the connection on earlier derivative slots.

## Norms and exact reconstruction

Let `gamma0=dy²+rho² kappa_S3`, and let `<.>` mean probability angular
average. Use

```
||T||_H13,prob²=sum_(j=0)^13 integral rho³|dy| <|nabla^jT|²>,
||Z||_H14(q,w),prob²
 =sum_(j=0)^14 integral rho³|dy| w² <|q^j nabla^jZ|²>.
```

All derivative and tensor slots use the fixed metric. Geometric
norms use sphere area measure and are `sqrt(2pi²)` times these norms.
The complete-line derivative-square norm is
`||V||_H14,s²=sum_(k=0)^14||D_s^kV||_L2(ds)²`.

For each unit harmonic with `lambda=8,24`, use the orthonormal
sphere-averaged tensor sections of `F=R e0+TS3`

```
S0=E e0* tensor e0*,
Sm=(e0* tensor dE+dE tensor e0*)/sqrt(2lambda),
St=E kappa_S3/sqrt(3),
Sf=[Hess E+(lambda/3)E kappa_S3]/sqrt(lambda d_lambda).
```

Here the split bundle identifies a tangential orthonormal component
with its unit-sphere component. For `lambda=0`, keep only `S0,St`.
With `U=sqrt(qrho³)w`, `V=U(u,v)^T`, and `A=B D_s+C` exactly as in
[the normal equation](MOMENTUM_REGULARITY.md), the lowered physical
tensor is

```
tau2=w²q²P_m*Z=m sum_i (AV)_i S_i,
m=kappa r w rho^(-3/2),
alpha=D_s log m=kappa²[(1+r)r_y-(3/2)ar²].              (1)
```

Lowering or raising with `gamma0` commutes with the covariant
derivative and preserves every norm here. The sign in `alpha` differs
from `b0=D_slogU=kappa²[(1+r)r_y+(3/2)ar²]`.

[The source and coefficient calculation](SOURCE_JETS.md) bounds all
`D_s` derivatives through fourteen of `q/rho`, `qa`, `q_y`, `b0`,
and `alpha` by1, and those of `C` by `10^-10` in operator norm.
For the estimates below the common entrywise ceiling18 suffices.
These coefficient bounds hold on the whole annulus, including tails.

## Explicit angular bounds

Restricted ambient coordinates `n_i` on the unit sphere satisfy

```
Hess n_i=-n_i kappa_S3,
|bar_nabla^j n_i|<=2^j.
```

The latter follows by differentiation, since `bar_nabla kappa_S3=0`,
`|n_i|,|dn_i|<=1`, and `|kappa_S3|=sqrt(3)`. The actual unit
harmonics are

```
E0=1,
E2=sqrt(48)[(n1²+n2²)/2-1/4],
E4=sqrt(720/19)[(n1²-n2²)²-(n1²+n2²)/2+1/12].
```

Their polynomial degrees are at most four and their absolute
coefficient sums are below32. Assigning each derivative slot to one
of at most four factors in the product rule gives

```
|bar_nabla^jE_lambda|<=32*8^j.                          (2)
```

The reference connection preserves `e0*` and the sphere metric.
Using the displayed definitions and `8<=lambda<=24`,
`d_lambda>=10/3`, it follows that every section above satisfies

```
||bar_nabla^hS_i||_L2(S3,prob)<=1000*8^h.               (3)
```

For example the trace-free section has bound
`32*8^h[64+lambda/sqrt(3)]/sqrt(lambda d_lambda)<1000*8^h`.
Thus the angular constants are obtained from the actual harmonic
polynomials rather than an unspecified norm equivalence.

## Full connection recursion

Choose an orthonormal frame parallel along `y`, with
`e_a=rho^-1 bar_e_a`. On the split bundle,

```
nabla_0=D_y,
nabla_a=rho^-1 bar_nabla_a+a K_a.
```

For a covariant rank-N tensor, each slot of `K_a` is a rotation of
operator norm at most1. Hence

```
||(K_aT)_a||<=3N||T||.                                (4)
```

The tensor map `T -> (K_aT)_a` is parallel for the reference
connection: it uses only `e0*`, the tangential identity, and the
sphere metric. The projections onto radial and tangential derivative
outputs are parallel too. This property accounts for subsequent
angular derivatives of these insertions without additional curvature
or connection coefficients.

Let `W_j=q^j nabla^j tau2`, of rank `j+2`. Exact differentiation gives

```
(W_(j+1))_0=D_sW_j-j q_yW_j,
(W_(j+1))_tan=(q/rho)bar_nabla W_j+(qa)K_(j+2)W_j.       (5)
```

The rank `j+2` includes every previous derivative slot. Equations
(4)--(5) therefore control the full covariant derivative tensor.

## Quantified coefficient bound

For a scalar coefficient set

```
J_M(f)=sum_(d=0)^M sup|D_s^df|/d!.
```

The product rule gives `J_M(fg)<=J_M(f)J_M(g)` and
`J_M(D_sf)<=(M+1)J_(M+1)(f)`. The common jet ceiling18 gives
`J_M(f)<54` for the coefficients used in (1), (5).

Expand `W_j/m` as coefficient functions times `D_s^kV` times sphere
tensors obtained from `S_i` by reference derivatives and parallel
insertions. Give a base section angular size1000, multiply that size
by8 for each reference derivative and by `3N` for each `K_N`.
Equations (3)--(4) and parallelism justify this bound even after
further derivatives. Let `F_j` sum the coefficient seminorms
`J_(13-j)` times these angular sizes. The at most eight entries of
`C` and the two nonzero entries of `B` give `F_0<10^6`.

At step `0<=j<=12`, the coefficient derivative, the combined
`alpha-jq_y` term, the derivative of `V`, the angular derivative,
and the connection insertion give respectively the factors

```
13-j, 54(j+1), 1, 54*8, 54*3(j+2).
```

Their sum is `824+215j<=3404<4000`. Thus

```
F_j<=10^6*4000^j<10^53,  0<=j<=13,
||W_j/m||_L2(ds x S3,prob)<10^53||V||_H14,s.             (6)
```

Every derivative of `V` that appears has order at most `j+1<=14`.
The quotient in (6) denotes this explicit expression; it assumes no
endpoint quotient trace.

## Physical H13 and zero extension

The volume conversion is exact: `rho³|dy|=q rho³|ds|` and
`sqrt(qrho³)m=qw`. Consequently

```
||nabla^jtau2||_L2,prob
 <=sup(wq^(1-j))F_j||V||_H14,s.                        (7)
```

For `j=0,1` the multiplier is at most1. For `2<=j<=13`, putting
`t=1/r>=1` gives

```
wq^(1-j)=kappa^(-2(j-1))t^(2(j-1))exp(-t)
 <=10^192 sup_(t>=1)t^24 exp(-t)
 <10^192*24^24<10^226.                                (8)
```

Equations (6)--(8), the fourteen squared derivative norms, and the
factor `sqrt(2pi²)<5` give per-sector bounds
`||tau2||_H13,prob<10^280||V||_H14,s` and
`||tau2||_H13,geometric<10^281||V||_H14,s`.
Summing the three sectors by triangle inequality costs at most
`sqrt(3)` relative to their square-summed potential norm. In
particular the rounded conversion

```
||tau2||_H13,geometric<10^300||V||_H14,s,total            (9)
```

holds. The transformed bound `||V||_H14,s,total<10^54` then proves
`||tau2||_H13,geometric<10^354` on each interior copy. We use the
public rounded ceiling

```
||tau2||_H13,geometric<10^355,                          (10)
```

which also bounds the square-summed norm on both copies: that norm
costs at most an additional factor `sqrt(2)`.

Take compact smooth `V_n` on the complete line converging to `V` in
H14. Each reconstructed `tau_n` is smooth and compactly supported
inside the correction annulus. Applied to differences, (9) shows that
their zero extensions converge in full geometric H13. The limit
equals `tau2` locally and vanishes outside the annulus. This proves
H13 zero extension without assuming an ordinary endpoint trace.

## Full weighted potential H14

Use reference covectors `R0=E e0*`, `Rg=dE/sqrt(lambda)` for the
nonconstant sectors, and only `R0` for the monopole. Their angular
jets obey the same `1000*8^h` bound. Write
`Z=U^-1 sum_i V_iR_i`, and
`q^j nabla^jZ=U^-1 P_jV`. The exact split recursion now has radial
part

```
D_sP_j-(b0+jq_y)P_j
```

and tangential part

```
(q/rho)bar_nabla P_j+(qa)K_(j+1)P_j.                   (11)
```

Here the radial derivative also differentiates the displayed
derivatives of `V`. The rank is `j+1`. Repeat the coefficient budget
with `J_(14-j)`, initial size below `10^4`, and `0<=j<=13`.
The step factors sum to

```
(14-j)+54(j+1)+1+54*8+54*3(j+1)
 =663+215j<=3458<4000.
```

Thus `F_j,potential<=10^4*4000^j<10^55` through order14. The weighted
volume factor cancels exactly, since `w sqrt(qrho³)U^-1=1`.
After summing derivatives and sectors and converting the sphere
measure, this proves

```
||Z||_H14(q,w),geometric<10^60||V||_H14,s,total<10^114.  (12)
```

The same compact approximation proves membership in the associated
weighted dotted completion. The square-summed norm on both copies
costs at most `sqrt(2)` and is below `10^115`.

All physical estimates here concern the ordinary tensor
`tau2=pi2/rho0`. With the induced norm and connection on the fixed
density bundle, `nabla(rho0tau2)=rho0nabla tau2`; hence the invariant
canonical-density norm agrees. An unrelated coordinate-component
norm is not specified. The arithmetic bounds used above are included
in [the replay instructions](REPRODUCE.md).

The smooth branch, canonical source, and selected inverse remain
analytic premises. These bounds apply to the new weight member and
do not identify an unspecified historical member. They are initial
analytic ceilings; they give no evolved fields, uniform nonzero
amplitude interval, or nonlinear remainder estimate.
