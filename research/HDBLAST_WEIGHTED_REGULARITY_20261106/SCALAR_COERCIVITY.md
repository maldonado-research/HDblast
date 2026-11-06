# Scalar weighted coercivity for the specified initial-data member

6 November 2026. Conditional on the retained exact-initial-data theorem,
actual background and canonical constraint conventions in
[EXACT_INITIAL_THEOREM.md](inputs/EXACT_INITIAL_THEOREM.md),
[REGISTERED_EXACT_FAMILY_BINDING.md](inputs/REGISTERED_EXACT_FAMILY_BINDING.md),
[CANONICAL_CONSTRAINTS.md](inputs/CANONICAL_CONSTRAINTS.md), and
[ACTUAL_POTENTIAL_JET_BOUNDS.json](inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json).
The complete scalar adjoint and its compact comparison are established in
[PRIOR_SCALAR_OUTPUT.md](inputs/PRIOR_SCALAR_OUTPUT.md).

This note proves the numerical low-order potential estimate used by
[SCALAR_REGULARITY.md](SCALAR_REGULARITY.md). It concerns the actual
selected scalar normal potential for one bulk copy. Sphere measure is
normalized to probability throughout the proof. Full geometric sphere
measure multiplies each displayed norm by `sqrt(2 pi²)=pi sqrt2<4.45`;
identical reflected bulk copies contribute a further `sqrt2` to the
norm on their disjoint union. These factors are kept separate.

The weights below specify a new member of the admissible construction.
They do not identify the defining function or weight ratio of an earlier
member whose numerical choices were unspecified.

## Result

Use `ell=10^-6`, `rho in[129,131]`, `a in[.111,.112]`,
`|p| in[.00012,.00014]`, the actual background ODEs, and

```
r=(z−ell/2)(5ell/2−z)/ell²,
kappa=ell/100=10^-8, q=kappa²r², w=exp(−1/r).
```

The norm in this note is `||u||_w²=integral w²|u|² rho³|dy|`, with
unit probability-sphere norm for the angular harmonic. For each scalar
harmonic potential `N=n(y)E_l`, `lambda_l=0,8,24`, the following estimates
hold on the dotted weighted domain:

```
||N||_w <= C_l kappa^-4 ||q² P_g*N||_w,
C_0=1.3e24, C_2=1.54e8, C_4=4.34e7.                    (1)
```

For the specified kappa, these constants are respectively
`1.3e56`, `1.54e40`, `4.34e39`. They are deliberately loose. They prove
an effective positive gap on the three actual scalar channels. No
unweighted trace of the potential is assumed.

Set `Q=P_g*N`, `E=||r⁴Q||_w`, `T=||N||_w`, and
`I=||N||_(w,{r>=1/5})`. The proof consists of

```
T <=.4 E+37 I,                                         (2)
I <=625*6650 E  (lambda8),
I <=625*1873 E  (lambda24),                             (3)
I² <=9e20 E T  (lambda0).                              (4)
```

Since `q²=kappa⁴r⁴`, these imply(1). The constants do not involve an
unknown inverse or a boundary-value solution.

## Boundary-tail Hardy estimate

First prove everything for compact smooth potentials. On either component
of `{r<=1/4}`, put `t=1/r>=4`; t increases toward the annulus boundary.
The signed derivative `k=t_y` satisfies

```
|k|=2 t² sqrt(1−1/t)/(rho ell),
k^-2=rho²ell²r⁴/[4(1−r)].                              (5)
```

The transformed radial measure, including the weight, is

```
W(t)dt=exp(−2t) rho³ |dy/dt|dt,
rho³|dy/dt|=rho⁴ell/[2t²sqrt(1−1/t)].
```

Its logarithmic derivative obeys

```
(log[rho³|dy/dt|])_t=4a/k−2/t−1/[2t²(1−1/t)] <=0.
```

Indeed `|a/k|<=a_max*rho_max*ell/[2sqrt(3/4)t²]`; its fourfold
coefficient is less than8, whereas `2t>=8`. Hence `−W_t/W>=2`.
For compact H2 functions f with f(4)=f_t(4)=0, integration by parts and
Cauchy–Schwarz give

```
||f||_W <=||f_t||_W <=||f_tt||_W.                       (6)
```

For example `2 integral W f²<=integral(−W_t)f²
=2 integral W f f_t`; all endpoint terms vanish for these compact test
functions. Applying the same argument to f_t proves the second step.

Let B be the angular-trace coefficient of Q, so

```
B=[−n_yy−2a n_y+C_lambda n]/2,
C_lambda=2lambda/(3rho²)−R_t,
R_t=3/rho²+p²/3−3a².
```

In the tail coordinate this is

```
n_tt+alpha n_t−beta n=−2B/k²,
alpha=2/t+1/[2t²(1−1/t)]+a/k,
beta=C_lambda/k².                                     (7)
```

For all three lambda, the actual boxes give
`|alpha|<.55`, `|beta|<.01`, `|k^-1|<4.73e-6`.
The beta ceiling is extremely coarse; its direct bound is below1e-12.

Choose chi(t)=0 for t<=4, chi(t)=1 for t>=5, and
`chi=3(t−4)²−2(t−4)³` on[4,5]. This is C1 and piecewise C2,
with `|chi'|<=2`, `|chi''|<=6`; it is adequate for(6), or can be
approximated in H2 by smooth cutoffs. Put f=chi n. From(7),

```
f_tt+alpha f_t−beta f
=−2chi B/k²+2chi' n_t+(chi''+alpha chi')n.
```

Using(6) and absorbing .55+.01 on the left gives

```
||chi n||_W
 <=2.5{||2chi B/k²||_W
       +||2chi'n_t+(chi''+alpha chi')n||_W}.            (8)
```

For the radial coefficient A=Q_yy/E_l,

```
n_y=b_lambda n−2A/(3a),
b_lambda=a+(lambda−3)/(3a rho²), |b_lambda|<.116.
```

Thus on the transition4<=t<=5,

```
|n_t|<=2.85e-5|A|+5.5e-7|n|.                          (9)
```

In(8), the first term is at most `3e-8 E` by(5), after multiplication
by2.5. The transition A term is at most `.178125 E`, since there
`r>=1/5`, `|A|<=|Q|<=625|r⁴Q|`. The transition n coefficient is
less than18. Consequently each tail satisfies

```
||chi n||_W <=.181 E+18 I.                             (10)
```

The two tails are disjoint. The remainder
`(1−chi_left−chi_right)n` is supported where r>=1/5 and is bounded
in norm by I. Summing(10) for the two tails gives(2). This is not a
Dirichlet estimate on the original potential; the only imposed endpoint
values belong to the auxiliary cutoff f at the interior point t=4.

## Nonconstant modes

The angular trace-free block of the adjoint is
`n T_l/(2rho²)`, in the background orthonormal angular frame, where
`||T_l||²=2lambda(lambda−3)/3` for unit E_l. Hence

```
|n|<=sqrt[6rho⁴/(lambda(lambda−3))] ||Q||_sphere.
```

The bracketed gain is below6650 for lambda8 and1873 for lambda24.
Since r>=1/5 in I, this proves(3), including the squared spherical
normalization. Equations(2),(3) give

```
T<(.4+37*625*6650)E<1.54e8 E,
T<(.4+37*625*1873)E<4.34e7 E.
```

## Monopole interior estimate

The exact adjoint identity derived in [PRIOR_SCALAR_OUTPUT.md](inputs/PRIOR_SCALAR_OUTPUT.md) is

```
n=c A+b(div Q)_y,
c=−2rho²/3, b=−2a rho²/p².                            (11)
```

Choose eta(r)=0 for r<=1/6, eta(r)=1 for r>=1/5, and the same cubic
step on the intervening interval of length1/30. Then
`|eta_r|<=45`, `|eta_y|<700000`. Its support is compact inside A.
Write `l=(log w)_y`; on this support `|l|<560000`.

Multiply(11) by eta²w² n and integrate in the fixed background volume.
The adjoint of div is minus sym grad. The angular trace contributes
`−3ab eta²w² nB`. After substituting
`n_y=b0 n−2A/(3a)`, `b0=a−1/(a rho²)`, the exact identity is

```
integral eta²w² n²
=integral w²{eta²(c−b_y−2bl−bb0)nA
             −2b eta eta_y nA−3ab eta²nB
             +[2b/(3a)]eta² A²}.                     (12)
```

The last term is nonpositive because a>0 and b<0. Drop it for an
upper bound. Use `|A|<=|Q|`, `sqrt3|B|<=|Q|`, together with

```
|c|<=11441, |b|<2.67e11, |b_y|<1.545e12,
|b0|<.113, sqrt3<7/4.
```

The coefficient multiplying `|n||Q|` is at most

```
K=11441+1.545e12+2*(2.67e11)*560000
  +(2.67e11)*.113+2*(2.67e11)*700000
  +(7/4)*.112*(2.67e11).
```

Since r>=1/6 on this support, `|Q|<=1296|r⁴Q|` there, and
the exact rational arithmetic gives `1296K<9e20`. Weighted
Cauchy–Schwarz proves(4). This estimate involves T on the right;
it is nevertheless enough to close(2).

Indeed `T<=.4E+37sqrt(9e20 E T)`. Young's inequality gives

```
37sqrt(9e20 E T)<=T/2+(37²*9e20/2)E,
T<(.8+37²*9e20)E<1.3e24 E.
```

The mechanism retains the small-p conditioning in b and b_y. It does
not assert a scalar-gradient-independent monopole gap.

## Dotted closure and excluded kernel

The inequalities above hold for compact smooth potentials in each
harmonic channel. The fixed dotted weighted H16 domain is the closure
of such potentials at k12. Its norm controls both `||N||_w` and
`||q²P_g*N||_w`, because the second derivative weight is q² and the
background lower-order multipliers are smooth and bounded on the
closed annulus. Passing to the limit proves(1) on that domain. No
claim of finite unweighted potential norm or an unweighted zero trace
is made at any stage.

An estimate including lambda3 would be false: the degree-one boost
lapse `N=rho E_1` belongs to the dotted exponential weighted domain
and has `P_g*N=0`. The present source/reflection sector excludes it.

Equation(1) is an explicit low-order coercivity constant. The following
energy comparison evaluates its right side for the selected potential;
[SCALAR_REGULARITY.md](SCALAR_REGULARITY.md) supplies the distinct
high-order argument.

## Numerical bound for the actual selected potential

Write `h2=w²q⁴P_g*N` and

```
E_nat=||q²P_g*N||_w=||h2||_(mu_H^-1), mu_H=w²q⁴.
```

The selected solution satisfies the weighted normal equation and is
orthogonal, in the correction energy, to every finite-energy homogeneous
correction in the same physical extension class. To justify this statement,
choose compact smooth potentials N_j converging to N in its dotted domain.
For a homogeneous correction k, the distributional identity P_g k=0 gives
`integral (P_g*N_j):k=0`. Weighted Cauchy–Schwarz passes this identity to N,
so `||h2||_(mu_H^-1)<=||h_part||_(mu_H^-1)` for each admissible same-source
particular h_part. No unweighted potential trace is used.

The compact scalar-frozen particular in
[PRIOR_SCALAR_OUTPUT.md](inputs/PRIOR_SCALAR_OUTPUT.md) has probability
L2 norm less than5.58e39 and support where r>=3/4. Hence

```
E_nat <= kappa^-4 (4/3)^4 exp(4/3) ||h_part||
      < kappa^-4 *2.72*4.42*5.58e39.
```

Here `e<2.72` and `(4/3)^4 exp(1/3)<4.42` are the rational-safe
exponential bounds used in that comparison. The three orthogonal scalar
channels have coercivity constants at most `1.3e24 kappa^-4`. Therefore

```
||N||_w <1.3e24 kappa^-4 E_nat
        <8.72104896e128 <9e128, kappa=10^-8.                (13)
```

The norm in(13) is an auxiliary potential norm on one bulk copy with
probability sphere measure. It is not a physical field measurement.
For full geometric sphere measure it is multiplied by `pi sqrt2`, and
for two identical reflected bulk copies by another `sqrt2`.

The coefficient inequalities and exponential bounds are replayed as
specified in [REPRODUCE.md](REPRODUCE.md). Finite arithmetic checks do
not replace the continuous Hardy, compact pairing, and closure arguments
above. No evolved estimate or uniform nonlinear amplitude interval is
asserted here.
