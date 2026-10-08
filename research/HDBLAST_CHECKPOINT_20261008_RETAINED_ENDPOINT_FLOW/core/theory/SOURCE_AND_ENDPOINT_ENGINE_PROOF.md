# Independent source and endpoint-engine proof

8 October 2026 UTC. Static review of the candidate `implementation/source/later_source.py`, immutable helper `source_algebra_baseline.py`, and `implementation/engine/endpoint_engine.py`. No function constructing a registered source, retained decoder, or physical endpoint evaluator was called in this review. The engine is reviewed as mathematical code; authorization and byte custody belong to the separate execution guard.

## Analytic source disks: the existing constant 64 is valid

For a real center `c in [-9/2,-7/2]` and complex `eta` with `|eta-c|<=1/8`, put `z=eta+4`, `q=1-z^2`, `L=-1/eta`. Then

```
|z|<=r=5/8,    |q|>=39/64,    |1/q|<=v=64/39,
|eta|>=27/8,  |L|<=ell=8/27.
```

All denominators avoid zero on and inside the disks. For any complex `w` of modulus `s<1`,

```
Re(1/(1-w)) >= 1/(1+s).
```

Indeed, writing `x=Re w`, the difference has numerator `(1-s)(s+x)>=0` over the positive denominator `(1-2x+s^2)(1+s)`. Set `w=z^2`. Thus

```
Re(1/q)>=1/(1+|z|^2)>=64/89,
|B(z)|=|exp(1-1/q)|<=exp(25/89)<=1/(1-25/89)=89/64.
```

The last inequality follows term by term from `1/n!<=1` in the positive exponential series. It does not evaluate the physical bump at any real point.

Direct analytic differentiation gives

```
B'=-2z B/q^2,
B''=[4z^2/q^4-2/q^2-8z^2/q^3]B.
```

The following exact rational constants therefore bound `B,B',B''`:

```
b0=89/64,
b1=2r v^2 b0=7120/1521,
b2=(4r^2 v^4+2v^2+8r^2 v^3)b0=98210432/2313441.
```

For `h=B`, the source `g=4L^2 h-2L h'-h''` has bound

```
4ell^2 b0+2ell b1+b2=951819044/20820969 <64.
```

For the signed source `h=zB`, bounds are `|h|<=r b0`, `|h'|<=b0+r b1`, and `|h''|<=2b1+r b2`, so

```
4ell^2 r b0+2ell(b0+r b1)+(2b1+r b2)
   =3227905133/83283876 <64.
```

Both bounds also imply `|Lg|<32`, although the later endpoint engine needs only `g`. These estimates cover every closed radius-1/8 disk for all 64 source-cell centers, not just the physical real segment.

Cauchy's estimate gives `|a_n|<=64/(1/8)^n` for the exact Taylor coefficients of `g` about any center. On half-width `H=1/128`, the remainder after degree24 is bounded by

```
64 sum_(n>=25)(H/(1/8))^n
   =64*(1/16)^25/(1-1/16)
   =1/(15*2^90).                                    (S1)
```

The corresponding `Lg` bound is `1/(15*2^91)`. This independently proves the constants already present in the helper; no guessed extrapolation of the prehistory endpoint-cap proof is used. The later interval is strictly interior to compact support, so there is no omitted endpoint cap in this source representation.

## Source coefficient construction

At a rational center `c`, use local coordinate `x=eta-c` and `z0=c+4`. The helper constructs coefficients of the reciprocal of

```
d(x)=1-z0^2-2z0 x-x^2.
```

The recurrence `v0=1/d0`, `vn=-(sum_(j>=1)d_j v_(n-j))/d0` is obtained by equating coefficients in `d(x)v(x)=1`. The bump exponent is `f(x)=1-v(x)`. Its normalized exponential `beta(x)=exp(f(x)-f(0))` satisfies `beta'=f' beta`; equating coefficients gives `beta0=1`, `beta_n=(sum_(j=1)^n j f_j beta_(n-j))/n`. The helper implements exactly these rational recurrences, retaining coefficients through degree26 because `g` contains `h''`.

At every center, `f(0)=-z0^2/(1-z0^2)` lies in `[-1/3,0]`. The even-degree200 alternating exponential sum is an upper bound for `exp(f(0))`; adding the degree201 term gives the lower bound. Absolute term magnitudes decrease because `|f(0)|<=1`. All sums and coefficient multiplications use exact `Fraction` arithmetic.

For the signed source, the coefficient of `x^n` in `(z0+x)beta(x)` is `z0 beta_n+beta_(n-1)`, with `beta_-1=0`. The geometry coefficients `L_n=-(-1)^n/c^(n+1)` are the exact Taylor series of `-1/(c+x)`. The convolution and factors `(n+1)` and `(n+1)(n+2)` in the helper therefore produce the exact normalized degree24 coefficients of `4L^2h-2Lh'-h''`.

Each exact physical Taylor coefficient lies in an interval obtained by multiplying its normalized rational coefficient by the enclosed `exp(f(0))`, swapping bounds when its sign is negative. The chosen coefficient is the downward512-bit dyadic rounding of the interval midpoint. If that point is outside the interval by a rounding unit, the implemented radius `max(chosen-lo,hi-chosen)` still bounds both endpoints and the whole interval. Consequently

```
|g(c+x)-sum chosen_n x^n|
 <= analytic_tail+sum coefficient_radius_n H^n       (S2)
```

for every real `|x|<=H`. The source certificate's complete uniform error includes (S1), the alternating-series exponential uncertainty, and chosen-coefficient rounding exactly once. The new `later_source.py` uses only the immutable helper's pure formal-algebra primitives; it does not pretend that the old scope-specific source authorization grants a new run.

## Endpoint moment formula used by EndpointFlow

Fix `t=-4` or `t=-7/2` and `D=t-a`, respectively `1/2` and `1`. The 64 contiguous source cells cover `[a,b]`; the first32 cells cover `[a,-4]` exactly. The engine validates every supplied row's center/half-width before selecting the prefix whose right endpoints are at most `t`. There is no partial cell and no extrapolated source.

For the chosen piecewise polynomial `p(s)`, define real moments

```
m_n(t)=integral_a^t p(s)[2(t-s)]^n ds.
```

Uniform absolute convergence of the exponential series on the compact domain gives the exact forced response

```
W_forced(t)=-sum_(n>=0) i^n k^n m_n/n!,
U_forced(t)=-sum_(n>=0) i^n k^n m_(n+1)/[2(n+1)!].   (E1)
```

These follow directly by expanding `E_k(t-s)` and `Phi_k(t-s)`. Their zero-momentum values are `-m0` and `-m1/2`, respectively. There is no division by momentum and no homogeneous-mode normalization assumption.

For one cell, put `y=2(t-s)`, `shift=t-center`, so `s-center=shift-y/2` and `ds=-dy/2`. The exact rational transformation

```
p(s)=sum_j d_j y^j,
d_j=sum_(m>=j) coefficient_m binom(m,j) shift^(m-j)(-1/2)^j
```

is the code's `polynomial_in_lag`. With `y_left=2(t-left)` and `y_right=2(t-right)`, its contribution to `m_n` is

```
(1/2) sum_j d_j [y_left^(n+j+1)-y_right^(n+j+1)]/(n+j+1). (E2)
```

The reversed order of the `y` endpoints and the `1/2` Jacobian are both essential. The code implements (E2); all `d_j` are first obtained exactly and their subsequent finite arithmetic is enclosed by Arb balls. The coefficient arrays use the cycle `1,i,-1,-i` and the negative signs in (E1). The static source reviewed contains all25 coefficients in every cell, rather than dropping small or high-degree terms.

## Kernel tails, source transport, incoming state and export

The code's exact polynomial L1 upper is

```
A=sum_cells (2H) sum_j |coefficient_j| H^j >= integral_a^t |p(s)|ds.
```

For `0<=k<=256` and `0<=t-s<=D<=1`, let `xmax=2*256*D<=512`. Truncating (E1) at degree `N=2048` gives, with

```
T=xmax^(N+1)/(N+1)! / [1-xmax/(N+2)],
```

the bounds `error_W<=A*T`, `error_U<=A*D*T`. Every omitted exponential term after the first has successive absolute ratio at most `xmax/(N+2)<1`. The `Phi` series has the still smaller factorial `(n+1)!` and factor at most `D`, so the code's `A*D*T` is conservative. The tail depends on the chosen polynomial's L1 majorant, which the engine computes, not on an unproved smallness assumption for the source coefficients. All tail comparisons are exact rational comparisons.

The reviewed snapshot additionally rounds each nonnegative tail upward to the fixed2^-512 dyadic grid using the exact integer ceiling of `tail*2^512`. This increases the bound by less than2^-512 (or leaves an exactly aligned bound unchanged), avoids enormous factorial denominators in serialized budgets, and preserves containment. This enlarged tail is what enters the target rectangle and is counted once. Its value need not be as small as the unrounded analytic tail.

If cell `j` has certified error `q_j` from (S2), then

```
source_error_W=sum_j 2H q_j,
source_error_U=sum_j 2H(t-center_j) q_j               (E3)
```

bound the true-source response discrepancy. The first is the integral of `q`; the second is the exact integral of `(t-s)q`. The earlier theorem supplies optional sharper `sin` kernels, but the implemented (E3) is already rigorous. Cell-source uncertainty is added once to the arithmetic and series-tail rectangles. There is no later-time cap. Prior incoming BD-cap/source uncertainty is not part of this same-saved-initial endpoint target and is added only if a separate total BD comparison is performed.

The full target is

```
W_target=E_k(D)W_saved(a)+W_forced,
U_target=U_saved(a)+Phi_k(D)W_saved(a)+U_forced.
```

The exact rational incoming complex components are converted to containing Arb balls. `E` uses a validated complex exponential; `Phi=D exp(ikD)sinc(kD)` uses Arb's documented unnormalized sinc `sin(x)/x`, which equals one at zero. Therefore incoming-state phase arithmetic is enclosed, including at very small momentum. The engine neither projects the canonical invariant nor assumes it vanishes.

The arithmetic premise is the soundness of the pinned Arb/python-flint ball operations at1024-bit precision. Ball evaluation includes rounding in moment powers, sums, polynomial coefficients, polynomial evaluation and incoming phase multiplication. Severe cancellation may widen an interval; it does not invalidate containment. Fixed precision is not a mathematical proof of small width. The export gate independently rejects a result whose complete L1 radius exceeds `10^-18`. A gate rejection must remain a rejected/insufficient-precision result, with no automatic degree/precision change unless newly registered.

Each real/imaginary ball component has certified lower/upper endpoints. Subtracting/adding the exact analytic error enlarges it, then flooring the lower and ceiling the upper to the2^-96 dyadic grid preserves containment. The code checks finite balls, ordering through the ball API, and a complete L1 export radius. A modulus error bound is safely added to both Cartesian components; this can overestimate L1 radius but never underestimates it. The dyadic export step is included in that final gate.

Exact later saved-minus-target subtraction provides interval defects at the observation times. Direct finite stress discrepancy uses those state intervals and the fixed `L(t)`, `k`, original weights and represented Pi. Complete contacts cancel for the mathematical comparison of two states evaluated with identical exact contacts. Stored full-stress arrays are not automatically certified by that cancellation; their direct evaluation/contact/accumulation discrepancies are additional targets.

## Exact canonical drift sharpens endpoint stress without projecting the data

For real `g`, the exact target conserves `c=Re U-Im W/(2k)`. Therefore the exact represented difference

```
dc=c_saved(t)-c_saved(a)
```

obeys `dc=Re D_U-Im D_W/(2k)`; it is not uncertain because of the enclosed target. With `Q=w/(4Pi^2)`, substituting `Re D_U=dc+Im D_W/(2k)` into the direct operators gives

```
mu D_R=Q[k(2k^2+3L^2)dc+(3L^2/2)Im D_W-kL Re D_W],
mu D_P=Q[k(2k^2/3-L^2)dc-(2k^2/3+L^2/2)Im D_W-kL Re D_W].
```

These affine functions attain their exact ranges on the supplied `D_W` component rectangle by endpoint sign selection. That rectangle can itself overestimate a correlated true target, so the resulting range is an enclosure, not necessarily the true attained range. The generic `D_U,D_W` direct-operator rectangle is independently valid. Intersecting the two stress intervals preserves containment. Likewise, the supplied `D_U.real` interval must intersect `dc+D_W.imag/(2k)`; failure is a consistency rejection. No identity is enforced by changing a saved state or narrowing an interval without a proved constraint.

Signed canonical sums can be accumulated exactly before taking absolute values. The quantity `2epsilon*dc` is a first-order canonical normalization drift; it is not the entire finite-epsilon nonlinear Wronskian defect. The fusion preserves the original coordinate data and does not turn a conserved or nearly conserved canonical diagnostic into an accuracy certificate.

## Scope of this review

The reviewed formulas are mathematically sound under the declared source and arithmetic premises. No static flaw was found in the lag transformation, factorial indexing, forced-response signs, source-error first moment, incoming propagation or zero-momentum branch. This statement is tied to the SHA256 file snapshot in the companion review receipt; later edits require reinspection.

This proof does not certify a numerical execution, authorize actual inputs, validate the hardware/compiler/Arb implementation independently, or reconstruct unsaved historical solver stages. It establishes endpoint target containment conditional on authenticated exact inputs and sound arithmetic. Numerical endpoint PASS is not the historical continuous pressure/contact gate and is not a continuum quadrature or ultraviolet result.
