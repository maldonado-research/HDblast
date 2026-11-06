# Affine canonical initial constraint source and momentum bound

6 November 2026. Conditional analytic calculation for the retained
`ell=10^-6` packet and background. Parameter coefficients are amplitude
derivatives, so a second derivative is twice its quadratic Taylor
coefficient. The scalar and scalar momentum density are prescribed affinely.

For an ordinary spatial tensor or one-form `T`, use the probability-sphere
spatial norm

```
N(T)² = integral rho³ dy <|T|²_gamma0>,
```

where `<.>` is the round probability mean on `S3`. The geometric spatial
norm is `sqrt(2pi²) N(T)`. For a canonical metric momentum density write
`pi2=rho0 tau2`, with `rho0` the fixed background volume density; all
momentum norms below are norms of the ordinary tensor `tau2`.

The model, canonical reconstruction and selected affine prescription are
those of [the registered binding](inputs/REGISTERED_EXACT_FAMILY_BINDING.md)
and [canonical constraints](inputs/CANONICAL_CONSTRAINTS.md). The retained
branch box is recorded in
[the actual potential and background bounds](inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json).
The analytic weighted inverse and zero-extension conclusions remain
premises from [the exact initial theorem](inputs/EXACT_INITIAL_THEOREM.md).
Original source paths, immutable Git references and content hashes are
recorded in [SOURCE_PINS.json](SOURCE_PINS.json).

## Background and scalar-frozen radial Hamiltonian

On the time-symmetric initial section,

```
gamma0=dy²+rho(y)² kappa_S3, a=rho_y/rho, p=phi0_y,
rho_y=a rho, p_y=V'-4ap,
a_y=-rho^-2-p²/3, V=p²/2-6a²+6rho^-2,
dy=-rho dz.
```

The retained collar box contains the complete correction annulus
`ell/2<z<5ell/2` and gives

```
129<=rho<=131, -14/10^5<=p<=-12/10^5,
111/1000<=a<=112/1000,
-9/10^5<=eta=phi0-1<=-8/10^5.
```

The potential is
`V=-2/27+(14/9)eta²+(50/27)eta³-eta^4/6-(4/9)eta^5-(2/27)eta^6`.
In this box `|V'|<3/10^4`, `3<V''<4` and `|V'''|<12`.

Let `Delta E=-lambda E`, `T=Hess E+(lambda/3)E kappa`,
`d_lambda=2(lambda-3)/3`, and write a scalar-sector metric perturbation as

```
h_yy=A E, h_ya=B D_aE,
h_ab=rho²(C E kappa_ab+D T_ab).
```

With the scalar variation frozen, the Hamiltonian metric derivative is

```
P_g h=E{(3a/2)A_y+[(lambda+6)/(2rho²)-V]A
 -(3/2)C_yy-6aC_y+(lambda-3)C/rho²
 -(lambda/rho²)(B_y+2aB)+lambda(lambda-3)D/(3rho²)}.
```

For `lambda=0`, omit `B,D`. Its formal adjoint is

```
P_g* N=(Hess N-gamma0 Delta N-N Ric+N dphi0 tensor dphi0)/2.
```

The scalar-gradient term remains in this metric adjoint. The four
metric-adjoint kernel lapses are the degree-one boosts `rho n_i`; the
momentum-adjoint kernel consists of rotations. The invariant scalar
channels used below have no such kernel component. Their boost and
rotation source pairings vanish by parity and angular integration.

## Exact quadratic source for the actual affine curve

All derivatives in the following formulas are `y` derivatives. To avoid
confusing the velocity profile with the localization weight, write
`w=ut` for the first velocity profile from the registered binding. The
complete first canonical datum is

```
h_yy=-4uY, h_ab=2rho²uY kappa_ab,
phi1=vY, pi1/rho0=tY partial_y tensor partial_y,
varpi1/rho0=chi Y,
t=-3w/(2rho), chi=vt/rho,
v=-3(u_y+2au)/p, vt=-3(w_y+2aw)/p,
Y=n1²-n2².
```

Define

```
F=-12u u_yy-12u_y²-144a u u_y
  -(192a²+96a_y+168/rho²)u²
  -(v_y+4up)²-V''v²-(6w²+vt²)/rho²,
J=(18u²-v²)/rho²,
R=8u_y t+2chi(v_y+4up),
L=-8ut+2chi v.
```

The exact constraint Hessian evaluated on the affine first tangent is

```
D²H=F Y²+J|DY|²,
D²(M/rho0)=R Y²dy+L Y dY.
```

The required second-order source is their negative. These formulas retain
the canonical density prescription. At time symmetry the Hamiltonian
derivative has no momentum terms, so changing affine physical versus
canonical second momentum data does not change its Hessian. It does
change the momentum source and full second canonical equation.

One direct derivation uses
`gamma_e=(1-4euY)dy²+rho²(1+2euY)kappa`. Its scalar-curvature coefficient is

```
[e²]R(gamma_e)=
 [-12u u_yy-12u_y²-144a u u_y
  -(192a²+96a_y+168/rho²)u²]Y²+(18u²/rho²)|DY|².
```

Since `H` contains `R/2`, this coefficient contributes directly to `D²H`.
The two momentum contributions to `D²H` are `-6w²Y²/rho²` and
`-vt²Y²/rho²`. The remaining scalar-gradient contribution is
`-(v_y+4up)²Y²-v²|DY|²/rho²`.

Use

```
H0=1, H2=(n1²+n2²)/2-1/4,
H4=Y²-(n1²+n2²)/2+1/12,
Y²=1/6+H2+H4, |DY|²=4/3+4H2-4H4.
```

The eigenvalues are `0,8,24`, and the squared probability norms are
`1,1/48,19/720`. In this unnormalized basis the second-order source is

| Mode | Hamiltonian `f_l` | Momentum radial `r_l` | Momentum gradient `g_l` |
| --- | --- | --- | --- |
| 0 | `-F/6-4J/3` | `-R/6` | absent |
| 2 | `-F-4J` | `-R` | `-L/2` |
| 4 | `-F+4J` | `-R` | `-L/2` |

The angular norms therefore satisfy the exact identities

```
<|D²H|²>=3F²/40+2FJ/5+38J²/15,
<|D²(M/rho0)|²>=3R²/40+L²/(5rho²).
```

## Evaluated source and source derivative bounds

The actual bump and amplitude profile are

```
Q=exp(4-ell²/[(z-ell)(2ell-z)]) for ell<z<2ell, zero otherwise,
alpha=(rho/rho_b)^(3/2) p/p_b, u=alpha Q/rho², w=exp(-z)u.
```

The background box gives `4/5<alpha<6/5`. Elementary global bump bounds are

```
0<=Q<=1, |Q_z|<=16/ell, |Q_zz|<=736/ell²,
|Q_zzz|<35000/ell³.
```

For example put `b=1/[x(1-x)]>=4`, `x=(z-ell)/ell`. Then
`Q_xx=(b^4-6b³+6b²)exp(4-b)` and `|b_x|<=b²`; the third derivative is bounded
by `(b^6+10b^5+24b^4+12b³)exp(4-b)`. Maximizing the four power terms and
using `e^-1<1/2` gives the bound `34201<35000`.

Writing `B=alpha/rho²`, `h=V'/p` and `K=h-9a/2` gives

```
B_y=KB,
u_y=B(Q_y+KQ),
u_yy=B[Q_yy+2KQ_y+(K²+K_y)Q],
Q_y=-Q_z/rho, Q_yy=Q_zz/rho²+aQ_z/rho,
K_y=V''-h²+4ah+(9/2)(rho^-2+p²/3).
```

Differentiating these identities and using the coupled background
equations gives the following outward bounds by exact rational arithmetic:

| Quantity | Upper bound on absolute value |
| --- | --- |
| `u, u_y, u_yy, u_yyy` | `10^-4, 10, 4*10^6, 1.3*10^12` |
| `v, v_y, v_yy` | `260000, 1.1*10^11, 3.3*10^16` |
| `w, vt, vt_y` | `10^-4, 270000, 1.1*10^11` |
| `F, J` | `1.22*10^22, 4.1*10^6` |
| `F_y, J_y` | `7.3*10^27, 3.5*10^12` |
| `R, L` | `4.7*10^14, 1.1*10^9` |

Thus `|(f0)_y|<1.217*10^27` and
`|(f2)_y|,|(f4)_y|<7.31*10^27`. The Hamiltonian source angular norm is
everywhere below `3.35*10^21`. At the exact bump center `z=3ell/2`,
`Q=1,Q_z=0,Q_zz=-32/ell²`; the same boxes prove `F<-3*10^18` and an angular
source norm greater than `8*10^17`. This is a nonzero initial forcing
statement, not a lower bound on an evolved response or a spatial L2 norm.

All source terms vanish outside the seed support `ell<=z<=2ell`. Hence
the probability-sphere support volume obeys
`integral_seed rho³dy<=131^4 ell=(17.161)²`. In particular,

```
N(s_M)<2.22*10^15, s_M=-D²(M/rho0).
```

## Compact momentum comparison and one specified new weight member

For the normalized momentum derivative `P_m tau=-2 div_gamma0 tau` and a
scalar source `r E dy+g dE`, an explicit compact particular is

```
tau_yy=tau_ya=0,
tau_ab=rho²(C E kappa_ab+D T_ab),
C=r/(6a), D=(C+g/2)/d_lambda   (lambda=8,24).
```

For the monopole use `C=r/(6a)` and omit `D`. Direct divergence proves
`P_m tau=r E dy+g dE`, with the same smooth compact support. For a
unit-norm harmonic, the actual background box proves

```
<|tau_part|²>=3C²+lambda d_lambda D²
           <=36²(r²+lambda g²/rho²).
```

The bound follows from an exact two by two positive-matrix check. It
applies to the sum of the three orthogonal channels.

To give an absolute output bound, specify the following **new** member of
the admitted weighted construction:

```
a_w=2, kappa=ell/100,
r(z)=(z-ell/2)(5ell/2-z)/ell²,
d=kappa r, beta=kappa, qweight=d², wloc=exp(-1/r).
```

This choice does not identify the previously unspecified historical
weights. It keeps the background, first tangent, affine scalar and
scalar-density data, and correction annulus fixed. For this member
`W_M=wloc² qweight²=kappa^4 r^4 exp(-2/r)`. The seed has `3/4<=r<=1`, so

```
sqrt(sup_annulus W_M/inf_seed W_M)=(4/3)² exp(1/3)<2.49.
```

Conditional on the retained weighted inverse theorem, the selected output
`tau2=W_M P_m* Z` minimizes `integral W_M^-1|tau|² dmu0` among corrections
with this source. Comparing with the compact particular gives

```
N(tau2)<36*2.49 N(s_M)<2*10^17.
```

The geometric source and output norms are respectively below `10^16`
and `9*10^17`. These are absolute initial L2 ceilings for the normalized
momentum tensor of this specified new member. They do not bound high
spatial derivatives, evolved second fields, a common nonlinear cylinder,
an amplitude radius, or a nonlinear remainder.

The accompanying exact verifier checks the finite algebra and rational
inequalities. It does not prove the retained analytic existence,
weighted-domain or background branch theorems.
