# Proper-clock background shoot and compact Gaussian tube

The result is conditional on the retained actual regular background branch and
its authenticated wider-collar first-exit proof. No background point, launch
radius, or proper clock is chosen from a floating-point integration.

The stage is fixed by `STAGE_FREEZE.json`. The length is `ell=1/1000000`;
`T=5ell/4`, `Z=ell/8`. Let `b=rho a`. On `0<=z<=ell/4` the source
`interval_core.background` evaluates the actual shell interval plus the
retained ODE integral image and supplies coherent derivative enclosures.
These are enclosures of the same actual branch. Their use does not select
independently varying coefficient functions. In particular,

```
rho_z=-rho^2 a,  a_z=1/rho+rho p^2/3,
b_z=1-b^2+rho^2 p^2/3.
```

The source bounds give `14<b<15`, `-216<b_z<0`. The checker actually uses
the narrower source-derived `m<=b<=M`, retaining all shell uncertainty.
The complete source closure, the actual first-jet source, and both old
capsule manifests are authenticated before importing repository Python.

The coordinate-time ODE is

```
z'=-tanh chi, chi'=-b(z), tau'=rho(z) sech chi,
z(0)=zi, chi(0)=tau(0)=0.
```

Thus `-Mt<=chi(t)<=-mt` and `z` increases. For `u>=0`,
`u-u^3/3<=tanh u<=u`; the first inequality follows from
`(tanh u-u+u^3/3)'=u^2-tanh(u)^2>=0`.
Consequently

```
m t^2/2-M^3 t^4/12 <= z(t)-zi <= M t^2/2.
```

In particular the two launch values `Z-M T^2/2` and
`Z-m T^2/2+M^3 T^4/12` bracket a root. The checker widens the outward
rounded bracket by one exact grid unit to record strict opposite signs.
The entire bracket is positive and lies inside the compact label interval.
All launches in `[0,Z]` remain in the proved background collar through T:
a first exit at `ell/4` contradicts `Z+M T^2/2<ell/4`.
Any positive launch above Z cannot hit Z, since `a_z>0`, `a(0)>0` keep
`b>0` on the retained inward regular branch, so `z` increases. The
subsequent Jacobi proof therefore gives uniqueness among positive launches
on that branch that are defined through T.

The first coordinate-time shooting derivatives are

```
j=partial_zi z|t, k=partial_zi chi|t,
j'=-sech(chi)^2 k, k'=-b_z j, j(0)=1, k(0)=0.
```

Put `K=216`, `H=8ell/5`. Up to any first time `|j|=2`,
`|k(t)|<=2Kt` and `|j(t)-1|<=Kt^2`. Since `KH^2<1`, that first time
cannot exist. Hence `j>=1-KH^2>0`. The sign of `b_z` then implies
`k>=0`, `j<=1`. This is a strict analytic bootstrap, not a finite sample
or an uncertified variational integration.

The elementary inequality `sech u>=1-u^2/2` yields

```
rho_min [t-M^2 t^3/6] <= tau(t) <= rho_max t.
```

For `0<=u<=sqrt(2)`, the derivative of
`cosh(u)(1-u^2/2)` is nonpositive because `tanh u<=u`; outside that
range its value is nonpositive. This proves the stated sech inequality
without a numerical transcendental decision. The target clock enclosure
comes from the integral on the exact implicit shooting root, with no
choice of a representative radius or clock inside the enclosure. The
rapidity, clock, and launch bounds are all uniform over the retained
actual-background uncertainty.

The proper-time tube is

```
0<=tau<=1/5000, ell/16<=zi<=3ell/16, whole round sphere.
```

At coordinate horizon H, the lower clock bound exceeds `1/5000`.
The target upper clock is strictly below `1/5000`. Monotonicity of the
clock therefore supplies the entire proper-time tube, and the radial
first-exit inequality keeps every trajectory below `ell/4`. Its chart
enclosure is `0<=x<=8/5`, `1/32<=r<=1/4`, with strict margins
`3-x>=7/5` and `8-x-r>=123/20`. Initial time is exactly the anchored
hypersurface `t=0`; future clock monotonicity supplies `t>=0`.

At fixed proper time, Gauss' lemma gives `g(U,J)=0`. Writing
`v=z_t=-tanh chi` and using the change of independent variable,

```
partial_zi z|tau=j cosh(chi)^2,
partial_zi t|tau=-j sinh(chi) cosh(chi),
partial_zi tau|t=rho j sinh(chi),
det d(t,z)/d(tau,zi)=j cosh(chi)/rho,
gamma_zi_zi=rho^2 j^2 cosh(chi)^2,
gamma_AB=(rho cosh(t))^2 kappa_AB.
```

Thus the proper-time radial connecting-field norm is
`A=rho j cosh chi`; there is no reciprocal cosh factor. The conserved
energy also gives `rho cosh chi=rho(zi)`, so `A=rho(zi)j`.
The report encloses A, the determinant, and the radial and angular Gram
margins. The Gaussian frame has zero cross terms by Gauss' lemma and
angular equivariance. Its outward radial normal is
`sinh(chi)e0+cosh(chi)ey`, since initial zi increases inward.

The map is injective on the whole stated tube: equal image times and
radii first force equal labels by strict `j>0` at that coordinate time,
then equal proper clocks by strict `tau_t>0`. Angles are unchanged.
Thus neither focal points nor trajectory crossings occur in this tube.

For the specified exact C2 metric family with piecewise spacetime C4
regularity, the background construction extends to a uniform open tube
by ordinary ODE dependence on parameters and labels. The positive clock,
determinant, and Gram bounds and the strict spatial and upper-time domain
margins persist on a common small parameter neighborhood. A finite compact
sphere atlas and C1 stability of the background embedding preserve
injectivity. This conclusion is qualitative; no effective epsilon radius
is claimed. The anchor keeps initial t=0 and positive future t_tau ensures
t>=0 throughout the family tube. This is a finite-dimensional observer
argument and does not prove a new nonlinear PDE theorem.

All exponentials used in boost-only API rows have rational alternating
Taylor degree 28 enclosures after argument reduction, with outward
squaring. The inherited `expbox(q)` means `exp(-q)`, so negative rapidity
is converted with `expbox(-chi)` for `exp(chi)`. No float is a proof input.

## Selected exact-family domain restriction

The selected source family is supplied only through coordinate time
`3ell/2`. The preceding parameter-continuity argument is applied to the
smaller proper-time tube `tau<=1/5500`, whose strict coordinate-time margin
is proved in `FAMILY_DOMAIN_JOIN.md`. The frozen `tau<=1/5000` background
tube remains an auxiliary overcover and grants no extension of the
selected exact-family slab. The frozen report explicitly requires family
definition near any compact tube before applying its stability statement.
