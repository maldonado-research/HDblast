# Weighted momentum-potential coercivity

6 November 2026. This calculation uses the fixed-density canonical
momentum prescription and background of
[the preceding source calculation](inputs/PRIOR_ACTUAL_SOURCE.md),
with the retained weighted-inverse premises of
[the initial-data theorem](inputs/EXACT_INITIAL_THEOREM.md).
The numerical member used here is newly specified:

```
ell=10^-6, kappa=10^-8,
ell/2<z<5ell/2,
r=(z-ell/2)(5ell/2-z)/ell²=1-((z-3ell/2)/ell)²,
q=kappa²r², w=exp(-1/r).
```

It is distinct from the historical member whose defining function and
weight parameters were not numerically specified. All estimates below
are conditional on the retained smooth background branch and selected
weighted inverse. They concern initial second amplitude derivatives.

## Fixed metric, density, and component norms

Write

```
gamma0=dy²+rho² kappa_S3, dy=-rho dz, rho_y=a rho,
129<=rho<=131, 111/1000<=a<=112/1000.
```

The background box is included in
[the retained branch bounds](inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json).
Use probability angular measure, with angular mean denoted by `<.>`,
and define

```
||T||_w²=integral rho³|dy| w² <|T|²_gamma0>,
N(T)²=integral rho³|dy| <|T|²_gamma0>.
```

These norms are for one interior copy; the same estimates apply on
either copy. Square-summing the two copies costs at most `sqrt(2)`.

The canonical momentum density is `pi2=rho0 tau2`, where `rho0` is
the fixed background volume density; `rho` is the warp factor. The
linear momentum operator and its adjoint are

```
P_m tau=-2 div_gamma0 tau,
P_m*Z=(Lie_(Zsharp)gamma0)^(ij),
```

where both tensor slots in the second formula are raised with the
fixed `gamma0`. Lowering those slots preserves the tensor norm.

Let `E` have unit probability norm and satisfy
`Delta_S3 E=-lambda E`. For `lambda=8,24`, put

```
T=Hess E+(lambda/3)E kappa_S3,
d_lambda=2(lambda-3)/3,
Z=u E dy+(rho/sqrt(lambda))v dE,  X=(u,v).
```

Then `<|Z|²>=u²+v²`, and direct covariant differentiation gives

```
<|P_m*Z|²>=4u_y²
 +2(v_y-av+sqrt(lambda)u/rho)²
 +12(au-sqrt(lambda)v/(3rho))²
 +4d_lambda v²/rho².                                  (1)
```

For `lambda=0`, use `Z=u E dy`, omit `v`, and retain only
`4u_y²+12a²u²`. In particular the mixed second square in (1) is absent.
These are exactly the scalar sectors of degrees zero, two, and four;
rotational vector harmonics are not included.

Set `L=<|P_m*Z|²>^(1/2)`. The last square, `d_lambda>=10/3`, gives
`|v|<36L`. The trace square gives

```
|u|<=(1000/111)[L/3+5|v|/(3*129)]<8L.
```

Thus `|X|<37L`; the monopole has the stronger `|u|<4L`.
Also, if `j=v_y-av+sqrt(lambda)u/rho`, then

```
|(u_y,j)|<=sqrt(3/4)L,
|av-sqrt(lambda)u/rho|<=sqrt(a²+lambda/rho²)|X|<(1/5)|X|.
```

Minkowski's inequality therefore proves, uniformly in the three sectors,

```
|X_y|<=L+(1/5)|X|.                                    (2)
```

## Hardy estimate on the two tails

On either tail put `t=1/r` and `sigma=sqrt(1-1/t)`. For `t>=2`,

```
rho³|dy|=h(t)dt, h=ell rho⁴/(2t²sigma),
|qD_y|=[2kappa²sigma/(ell rho)]|D_t|,
|D_t|<=(10^-4 kappa^-2)|qD_y|.                         (3)
```

The last inequality uses `rho<=131`, `sigma>=1/sqrt(2)`, and
`131/sqrt(2)<100`. The actual radial measure is retained throughout.
Its logarithmic derivative satisfies

```
h'/h=4a y_t-2/t-sigma'/sigma<0,
|4a y_t|<=2(112/1000)ell*131/(t²sigma)
         <0.000042/t²<2/t.
```

Consequently `m=h exp(-2t)` obeys `m'<=-2m`. For a compactly
supported smooth function on this half-line that vanishes at `t=2`,
integration by parts and Cauchy--Schwarz give

```
integral m f² <=-(1/2)integral m' f²
               =integral m f f_t,
||f||_L2(mdt)<=||f_t||_L2(mdt).                        (4)
```

The boundary term at the infinite end vanishes for these compact
tests. No endpoint trace is assumed for the completed potential.

Choose `chi=0` on `t<=2`, `chi=t-2` on `2<t<3`, and `chi=1` on
`t>=3`. Smooth approximation permits `|chi'|<=2` with the same
bounds. Let

```
energy=||qP_m*Z||_w,  interior=||X||_(w,t<=3).
```

Apply (4) to `chi X` on both disjoint tails, then use (2)--(3).
Since `q<=kappa²/4` on their supports,

```
||chi X||_w<=10^-4 kappa^-2 energy
             +5*10^-6||chi X||_w+2 interior,
||chi X||_w<1.00001(10^-4 kappa^-2 energy+2 interior).  (5)
```

On `t<=3`, `q>=kappa²/9`, so the algebraic estimate gives

```
interior<=333 kappa^-2 energy.                        (6)
```

The interior and the region where `chi=1` cover the annulus. Thus
(5)--(6) prove

```
||Z||_w<1000 kappa^-2 ||qP_m*Z||_w
        =10^19 ||qP_m*Z||_w.                          (7)
```

Indeed `333+1.00001(666+0.0001)=999.006760001<1000`.

## Dotted closure and the selected potential

Start with smooth compact scalar-gradient potentials in these three
sectors, and complete them in `||qP_m*Z||_w`. Equation (7) makes a
Cauchy sequence Cauchy in the actual weighted L2 norm. On every compact
subinterval, `q,w` are bounded below, and (1)--(2) identify its limit
with the local first-order graph closure. Equation (7) passes to this
completion. It proves absence of a kernel in this dotted scalar-sector
domain without imposing ordinary endpoint traces. Rotational sectors
would contain Killing fields and require a different kernel treatment.

The exact canonical source is radial/gradient. Before converting the
public `H0,H2,H4` harmonics to unit norm, its coefficients are

```
r0=-R/6, r2=r4=-R, g2=g4=-Lsource/2,
s_M=r_lambda H_lambda dy+g_lambda dH_lambda.
```

The formulas for `R,Lsource` and their density prescription are those
of [the preceding source calculation](inputs/PRIOR_ACTUAL_SOURCE.md).
The squared harmonic norms are `1,1/48,19/720`; conversion to unit
harmonics multiplies each radial and gradient coefficient by the
corresponding harmonic norm. Orthogonality sums (7) over all three
sectors.

The same source calculation gives a compact particular with
`N(tau_part)<=36N(s_M)` and `N(s_M)<2.22*10^15`. It is supported in
`ell<=z<=2ell`, where `r>=3/4`. Compact comparison in the actual
dotted domain for the selected output `tau2=w²q²P_m*Z` gives

```
energy<=36N(s_M)/[kappa²(3/4)²exp(-4/3)]
       <36(2.22*10^15)(16/9)4*10^16
       =5.6832*10^33<6*10^33.
```

Here `exp(4/3)<4`, for example by its positive Taylor series and a
geometric bound on the tail. Combining with (7) establishes

```
||Z||_w<6*10^52.                                     (8)
```

This is a bound for the same selected normal potential, using a
competitor only in the variational energy comparison. Its higher
derivatives are obtained separately in
[the regularity proof](MOMENTUM_REGULARITY.md). The finite arithmetic
supporting this proof is described in [the replay instructions](REPRODUCE.md).
