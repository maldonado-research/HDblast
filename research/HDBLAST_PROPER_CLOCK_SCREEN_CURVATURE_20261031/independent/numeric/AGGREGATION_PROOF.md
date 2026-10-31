# Independent numerical aggregation proof

The actual first-metric jet boxes remain inherited source inputs. This review
does not reprove the existing oracle or shrink its errors. It independently
checks all arithmetic following those boxes, using only stdlib Fraction.

The center trajectory is parameterized monotonically by coordinate time t.
The positive weight is `d tau/dt=rho sech(chi)<=rho`. Hence, on a scaled-time
cell `[xl,xh]`, `integral_cell max(|SP|,|SQ|)d tau` is at most
`ell*(xh-xl)*rho_max*S_cell`. Summing these nonnegative cell bounds is an
upper bound for every partial source integral up to the observation clock.
No assumption of independence, common extremizer, or source sign enters.

Let `W=sup max(|P|,|Q|)` through the selected proper time T. The triangular
variation system gives

```
W <= source_integral + M T W + K T^2 W/2.
```

With positive denominator this proves

```
W <= source_integral/(1-MT-KT^2/2),
|Y| <= T W,
|X| <= T W/(1-sup|as| T).
```

Here `X=rho xi^t` and `Y=xi^y`; they are exactly the two position quantities
in the full tilt formula. The root of the shooting problem is enclosed, so
all source/path and rapidity boxes must hold for every enclosed root. The
rapidity is nonpositive and decreasing when rho*a>0. Endpoint cosh/sinh
boxes therefore give the whole-path bounds `[1,cosh_chi_max]` and
`[sinh_chi_min,0]` used in every cell.

Every subsequent endpoint can be reproduced by the independent interval
class. Addition, multiplication, reciprocal, integer power, and coercion
round outward to the exact grid 2^-112. Source caps, proper-length weights,
M, K, a*s, denominators, and position bounds are then exact rational scalars.
Point first-Weyl invariants and both quadratic terms retain the ordinary
physical derivative factors. The full tilt is evaluated as

```
-4/R {[-Xcap,Xcap](c At+s K)
       +[-Ycap,Ycap](s Ay+c K)}.
```

It encloses all correlations by inclusion. The final coefficient is the sum
of the complete TT term, complete boosted quadratic term, and this complete
tilt interval. The TT term must already include its radial jet, inherited
boundary history, complete feedback error, and exact-family second-amplitude
normalization. The sign is certified only if this full interval excludes
zero.

The old auxiliary background proper-time slab 0<=tau<=1/5000 is not itself a
valid exact-family observation domain if its coordinate-time image extends
past the PR112 family horizon. Restrict the actual congruence to
`0<=tau<=1/5500` and verify, using common rho/rapidity bounds,
`(1/5500)*cosh_chimax/rho_min < (3/2)ell`, with the target proper-time upper
bound below 1/5500. This is a strict domain condition, separate from the
arithmetic sign. Spatial, Jacobian, and no-caustic margins remain required.
