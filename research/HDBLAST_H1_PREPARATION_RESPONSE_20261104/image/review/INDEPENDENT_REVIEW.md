# Independent check of the actual-image H1 obstruction

Independent analytic review, 4 November 2026. No background evaluation, oracle query or PDE solve
was used. The following is an analytic construction within the retained
compact mean-preparation image, rather than an example chosen in an arbitrary
function space.

Let `I=[1,5/4]`, `r0=9/8`, and take integers `n>=16`. Define

```
psi(u)=u^3 exp(1-1/(1-u^2))  if |u|<1, and 0 otherwise,
z_n(r)=n^(-5/2) psi(n(r-r0)),
h_n=3b z_n-z_n', chi_n=ell rho^3 z_n.
```

The bump is smooth and flat at both endpoints. Its first three center jets
are `psi(0)=psi'(0)=psi''(0)=0`, `psi'''(0)=6`. All supports lie strictly
inside `I`. The exact background relation `rho_r=-ell a rho^2` implies

```
chi_n,r=ell rho^3(z_n'-3b z_n)=-ell rho^3 h_n,
chi_n,y=rho^2 h_n.
```

Thus the weighted moment of `h_n` is exactly zero, since the primitive
`chi_n` is compactly supported. This identity uses no numerical cancellation.
The associated local `zeta_n` is compactly supported and flat by the retained
formula involving `h_n,y-a h_n` and `chi_n`. The existing admissibility
theorem therefore supplies the homogeneous metric addition while retaining
all canonical momentum and scalar canonical second sectors. Each fixed
finite `n` is a separate admissible direction; no family-uniform nonlinear
parameter interval is asserted.

Write

```
A=3/(2g), c=-2b-g, d=3(ell^2-b^2+bg-b'),
Q2=A'+Ac, Q1=A'c+Ac'+Ad, Q0=A'd+Ad'.
```

Direct substitution, with all signs retained, gives

```
D_n=A(z_n''+c z_n'+d z_n),
D_n'=A z_n'''+Q2 z_n''+Q1 z_n'+Q0 z_n.
```

Only `b in C2` and `g in C1`, with `g>0`, are needed to make these
coefficients uniformly bounded on `I`. In particular a C4 background metric
and scalar suffice for these finite estimates, since `a` and `p` involve
one background derivative. Literal smoothness of each admissible datum also
uses the retained exact background's smoothness. This should be distinguished
from the lower C4 topology in which parameter dependence was retained; C4
alone is not a statement of C-infinity background regularity.

At the center all lower-order terms vanish, so one has the exact identity

```
D_n'(r0)=6A(r0) sqrt(n)=9 sqrt(n)/g(r0).
```

There is no asymptotic error to dominate. The derivative supremum diverges
because `g(r0)>0` is fixed.

For `B_j=||psi^(j)||L2(-1,1)`, change of variables gives the exact scaling

```
||z_n^(j)||L2(I)=n^(j-3) B_j, 0<=j<=3.
```

Consequently, putting lower-case coefficient names inside suprema,

```
||D_n||2 <= ||A||inf [B2/n + ||c||inf B1/n^2 + ||d||inf B0/n^3],
||D_n'||2 <= ||A||inf B3 + ||Q2||inf B2/n
             + ||Q1||inf B1/n^2 + ||Q0||inf B0/n^3.
```

The first norm tends to zero and the second is uniformly bounded. Their sum
is an explicit finite upper bound `K` on the standard H1 norm for `n>=1`.
For every positive budget `E`, multiplying the entire preparation by
`t=E/(2K)>0` gives `||t D_n||H1<=E/2`, while
`|(t D_n)'(r0)|=9t sqrt(n)/g(r0)` still diverges. An energy budget of the
form `||D||H1^2<=E` is handled by replacing `E` in this scaling by `sqrt(E)`.
Thus no positive H1 budget forces the declared finite `R_A` threshold.

Fully explicit continuous bump constants can be used without evaluating an
exponential. With `v=(1-u^2)^-1>=1` and `E0=exp(1-v)`, differentiation gives

```
|psi| <= E0,
|psi'| <= (3+2v^2)E0,
|psi''| <= (6+14v^2+8v^3+4v^4)E0,
|psi'''| <= (6+54v^2+96v^3+96v^4+48v^5+8v^6)E0.
```

For positive integers `m`, `E0 v^m<=m^m`, by differentiation and
`exp(1-m)<=1`. Therefore valid rational supremum caps for orders zero
through three are `(1,11,1302,550638)`, and valid rational L2 caps are twice
these numbers. These bounds are coarse but continuous, including the flat
endpoints. The subsequence `n=m^2` makes every scaling power an integer:
`z_m=m^-5 psi(m^2(r-r0))`, and the exact derivative becomes `9m/g(r0)`.

For comparison, an added H2 hypothesis does control the desired profile
norm. For arbitrary `f in H1(I)`, the elementary estimate

```
||f||inf <= 2||f||2+(1/2)||f'||2 <= (sqrt(17)/2)||f||H1
```

gives `R_A(D)<=sqrt(17)/2 ||D||H2`. In the actual flat-at-contact image,
`D(1)=D'(1)=0`, so the stronger one-sided FTC bounds give

```
||D'||inf <= (1/2)||D''||2,
||D||inf <= (1/(8 sqrt(3)))||D''||2,
R_A(D) <= (1/2)||D''||2 <= (1/2)||D||H2.
```

The compactly supported counterexample also vanishes near the target, but
that second endpoint condition must not be imposed on the full retained
causal-segment preparation image. H2 control is an additional regularity
assumption; it has not been identified with a physical energy or a physical
preparation-selection rule.

The accompanying checkpoint now supplies a separate direct H1 energy route. This review concerns the actual-image obstruction to the earlier supremum-norm bridge; it does not independently review that separate positive result.
