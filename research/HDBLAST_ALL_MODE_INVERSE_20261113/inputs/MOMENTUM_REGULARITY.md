# Momentum normal equation and quantified H14 regularity

6 November 2026. Use the newly specified weight member, fixed-density
canonical source, and probability angular measure of
[the coercivity proof](MOMENTUM_COERCIVITY.md). The retained smooth
background branch and selected weighted inverse are premises from
[the initial-data theorem](inputs/EXACT_INITIAL_THEOREM.md).
The historical member with unspecified numerical weight parameters
is not identified by this calculation. Norms through (8) refer to one
interior copy; the final extension estimates state the two-copy factor.

## Exact unitary transformation

Radial background and weight coefficients preserve each scalar
radial/gradient harmonic sector under `P_m w²q²P_m*`. The adjoint
form also preserves its orthogonal complement. Uniqueness of the
retained selected inverse therefore leaves no unsourced sectors for
the actual source in degrees zero, two, and four. This sector
restriction follows from the operator and its source.

For a unit harmonic with `lambda=8,24`, write

```
Z=u E dy+(rho/sqrt(lambda))v dE,
X=(u,v)^T, k=sqrt(lambda)/rho, d_lambda=2(lambda-3)/3.
```

The normalized orthogonal tensor components of the lowered adjoint are

```
S X=(2u_y,
     sqrt(2)(v_y-av+ku),
     sqrt(12)(au-kv/3),
     2sqrt(d_lambda)v/rho)^T.                          (1)
```

Their Euclidean square is the angular mean of `|P_m*Z|²`. For the
monopole retain only `u` and the first and third rows; the second row
is absent when `lambda=0`.

Set `D_s=qD_y`, so `ds=dy/q` and the two endpoints become the two
infinite ends of a complete line. Define

```
U=sqrt(q rho³)w, V=UX,
b0=D_s log U=kappa²[(1+r)r_y+(3/2)ar²].
```

The volume factor `rho^(3/2)` and the factor `sqrt(q)` are included
in this unitary map. Direct substitution into (1) gives

```
||V||_L2(ds)=||Z||_w,
||A V||_L2(ds)=||qP_m*Z||_w,
A=B D_s+C,
B=[[2,0],[0,sqrt(2)],[0,0],[0,0]],
C=[[-2b0,0],
   [sqrt(2)qk,-sqrt(2)(b0+qa)],
   [sqrt(12)qa,-sqrt(12)qk/3],
   [0,2q sqrt(d_lambda)/rho]].                         (2)
```

The weak selected equation `P_m(w²q²P_m*Z)=s_M` is therefore
exactly

```
A*A V=G,
G=sqrt(q rho³)w^-1 f,
f=(r_lambda,sqrt(lambda)g_lambda/rho)^T,
A*A=-M D_s²+K1 D_s+K0,
M=B^T B=diag(4,2),
K1=C^T B-B^T C,
K0=C^T C-B^T C_s.                                    (3)
```

All transposes in (3) are Euclidean component transposes. In the
monopole `M=4` and `f=r0`. Source coefficients in (3) are those in
unit harmonics: multiply the coefficients of the public unnormalized
`H0,H2,H4` by their respective norms `1,1/sqrt(48),sqrt(19/720)`.
There is no change of source or canonical density prescription.

For clarity, the source transformation follows directly by testing:
`rho³dy=q rho³ds`, so the source pairing with the test potential
`U^-1 Phi` is `integral sqrt(qrho³)w^-1 f.Phi ds`.
Expanding `(-B^T D_s+C^T)(B D_s+C)` gives exactly the signs of
`K1,K0` in (3).

## Global regularity from an L2 distribution solution

Use the derivative-square norm

```
||V||_Hm²=sum_(j=0)^m ||D_s^j V||_L2(ds)².              (4)
```

Consider an L2 distribution solution of
`-MV''+K1V'+K0V=G`, with constant symmetric `M>=2I`, `G in H12`,
and

```
max_(i=0,1; 0<=h<=12) ||D_s^h Ki||_op,Linf
 <=epsilon<=10^-6.                                    (5)
```

It belongs to H14, with the explicit estimate

```
||V||_H14<10(||V||_L2+||G||_H12).                      (6)
```

To prove regularity first, the Fourier multiplier
`(I-MD_s²)^-1` gains two Sobolev orders. Rewrite the equation as

```
V=(I-MD_s²)^-1(G+V-K1V'-K0V).
```

Starting from L2, multiplication by `W1,infinity` coefficients acts
on H^-1 by duality, so the right side gives H1. If `V in Hr`, the
forcing belongs to `H^min(12,r-1)` by the product rule. Iteration
through `r=13` gives H14 on the complete line. Derivatives through
twelve of the coefficients suffice at the last step. No endpoint
trace enters this distribution argument.

For the constant, let `T=||D_s^14V||`. Fourier's inequality
`|xi|^j<=1+|xi|^14` gives `||D_s^jV||<=||V||+T` for `j<=13`.
Differentiate the equation twelve times. Since
`sum_h binom(12,h)=2^12` and `||M^-1||<=1/2`,

```
T<=(1/2)||D_s^12G||+4096epsilon(||V||+T).
```

Absorb `4096epsilon<=0.004096`; this gives
`T<0.503||G||_H12+0.004114||V||`. Finally
`||V||_H14<=sqrt(15)(||V||+T)`, proving (6).
For the Bessel convention with weight `(1+xi²)^14`, a conversion
factor below64 suffices because `max_j binom(14,j)=3432<64²`.

## Application to the actual canonical source

[The source-jet calculation](SOURCE_JETS.md) proves for each of the
three channels that `sum_(j=0)^12 ||D_s^jG_lambda||<10^17`.
Consequently the aggregate norm satisfies `||G||_H12<3*10^17`.
That calculation also gives whole-annulus operator-norm bounds
`||D_s^jC||<10^-10` through order fourteen. Thus, for `h<=12`,

```
||D_s^h K1||<=4*10^-10,
||D_s^h K0||<=2*10^-10+4096*10^-20<10^-6.               (7)
```

Here `||B||=2`, and the product derivatives of `C^TC` are bounded by
`sum binom(h,j)*10^-20<=4096*10^-20`. These are uniform coefficient
bounds on the complete annulus, including both tails.

[Coercivity and compact comparison](MOMENTUM_COERCIVITY.md) give
`||V||=||Z||_w<6*10^52` for the same selected solution. Equations
(6)--(7), summed over the three sectors, therefore yield

```
||V||_H14,total<10(6*10^52+3*10^17)<10^54.             (8)
```

The selected potential begins in the dotted weighted domain, so its
compact weak tests give (3); coercivity supplies its global L2
transform, and the preceding bootstrap supplies H14. This does not
substitute a compact particular for the selected normal potential.
Finite arithmetic joins are included in [the replay](REPRODUCE.md).

[The extension proof](MOMENTUM_EXTENSION.md) converts (8) to the full
geometric covariant weighted H14 bound `||Z||<10^114` and a physical
H13 zero extension of the ordinary momentum tensor
`tau2=pi2/rho0`. Its derived per-copy ceiling is below `10^354`;
`10^355` is the common rounded physical ceiling for either one copy
or the square-summed two-copy initial section. These are conditional
initial bounds. They provide no evolved estimate or uniform nonzero
amplitude interval.
