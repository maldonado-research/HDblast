# Actual combined TT jet at the proper-clock background event

This certificate supplies only the
metric-linear second TT term. It does not supply the complete proper-clock
observer coefficient or its sign.

The event remains `(x,r)=(T,R)=(5/4,1/8)`, with `L=T-R=9/8` and
`ell=10^-6`. The exact-family TT identification and the conditional actual
first-field, regularity and constraint premises are inherited from the
authenticated postcontact and exact-family capsules. No new source-oracle campaign is run.

## Sources and exact replay

`check_combined_tt.py` authenticates closed capsule manifests and every
listed member before reading scientific data and again before writing the
success report. It imports no project numerical module. The authenticated
PR111 producer, normal/optimized reports, source roles, initial TT data and
reflected-kernel derivation remain frozen.

The checker treats all old F/G ranges as **AUTHENTICATED_BOX_INPUTS**.
Authentication establishes retained bytes; it does not independently reprove
that these values enclose the physical source. Inherited actual-field errors
remain in every old source cell. The report explicitly records this scope.

Using integer/Fraction arithmetic, the checker independently verifies the 64
bulk-cap rectangles, eight full shell-history intervals, all 33 characteristic
intervals and their reflection split at `s=L`. Their union is covered by the
retained actual first-field domain. It reconstructs source-cap maxima,
operator-cap maxima, `sigma`, `q`, the full C1 norm and complete feedback tail.
It reproduces the existing actual Ux interval and positive TT factor on the
original `2^-112` outward lattice, including the outward conversion of the
nondyadic constant `ell`. Normal/optimized source reports must be canonically
equal.

## Endpoint signs and the combined expression

Write `Fplus(s)=F(s,T+R-s)`, `Fminus(s)=F(s,|s-L|)` and
`V=LG+K_N F`. Direct differentiation of the two integration endpoints gives

```
Vx = -G(L) + (1/2) integral_0^T [Fplus+Fminus] ds,
Vr = +G(L) + (1/2) integral_0^T [Fplus-Fminus] ds.
```

There is no derivative sign of the absolute-value fold in either endpoint
value. The fold appears inside the physical source argument. The plus/minus
sign in `Vr` comes from differentiating the upper/lower integration limits.
The inhomogeneous Neumann lift is `LG=-integral_0^(x-r)G`, so its time/radial
derivatives are `-G/+G`.

Set

```
A = (1/2) integral_0^T Fplus ds,
B = (1/2) integral_0^T Fminus ds.
```

For a single background rapidity chi,

```
cosh(chi)Vx-sinh(chi)Vr
  = -exp(chi)G(L)+exp(-chi)A+exp(chi)B.
```

The F boxes enclose the entire physical source on each cell. Multiplication
by the exact positive cell length encloses its integral without an uncharged
quadrature remainder. Exponential weights are positive and independent of
the integration variable. Thus weighted integrals can be aggregated as A/B
before multiplication. The calculation uses interval products, never
cancellation between independent source errors.

The pinned background report encloses chi uniformly over all coherent
shooting roots in its certified bracket. Its target is the same Pdagger.
Every selected background input pin is checked against the same repository,
and the actual-source capsule manifest agrees with the TT source capsule.
The apex rho boxes are identical. The background's slightly tighter
cosh(t) range nests in the old range and contains this checker's independent
rational Taylor enclosure. Consequently all source, q and coefficient-factor
premises refer to the same actual background data and event.

The checker independently encloses exp(chi) and exp(-chi). For
`y=|x|/2^k<=1/2`, it sums the positive exponential Taylor series through degree
64. The tail after degree 64 is bounded by
`(y^65/65!)/(1-y/66)`. Squaring and, for negative arguments, reciprocal
operations are rounded outward on the `2^-192` lattice. This bounds the full
chi interval by monotonicity; no numerical root endpoint is mistaken for the
actual root. Repeated interval occurrences deliberately overenclose shared
dependence.

## All feedback and the radial jet

The actual target is

```
U=LG+K_N(F+A_T Ux+B_T Ur+D_T U),
X(u)=max((2/T)sup|u|,sup|ux|,sup|ur|),
sigma=a_star+b_star+T*d_star/2,  q=T*sigma<1.
```

Both first derivatives of `K_N f` are bounded by `T sup|f|`. Thus the same
complete resolvent estimate that was previously charged to Ux also bounds
Ur:

```
X(U) <= (2*g_star+T*f_star)/(1-q),
d = q*(2*g_star+T*f_star)/(1-q),
|Ux-Vx| <= d,  |Ur-Vr| <= d.
```

These bounds sum every scalar feedback iterate. They do not replace U by
the finite W1 approximation. The radial interval is `Vr+[-d,d]`. The
combined error obeys

```
|cosh(chi)(Ux-Vx)-sinh(chi)(Ur-Vr)|
 <= (cosh(chi)+|sinh(chi)|)*d = exp(|chi|)*d.
```

The report uses the supremum of this gain over the entire chi interval.
Chi uncertainty is charged in both exponential weights and this feedback
gain.

The positive metric-linear coefficient factor is
`3/(2*ell*rho(P)^2*cosh(ell*T))`. Multiplying it by the combined actual jet
gives the second-derivative TT contribution. There is no additional factor
one-half. The complete boosted quadratic term and surviving angular-tilt
term must be added elsewhere before asserting a full-observer sign.

## Subprocess controls and replay

`run_controls.py` runs new normal and `-O` reports and requires canonical
equality. Seven genuine corrupted checker subprocesses exercise the Neumann
sign, reflected time sign, exponential weight, radial fold sign, feedback,
full q and second-derivative normalization. Exact nonzero-boost and polynomial
source fixtures catch these changes. Two isolated byte-corruption
subprocesses exercise the source-capsule and background-report pins. Each
must exit unsuccessfully at the expected witness without a success artifact.
The frozen source tree is never edited.

```
python -I -S -B check_combined_tt.py \
  --repo-root /path/to/retained-checkout \
  --background-report /path/to/BACKGROUND_NORMAL.json \
  --background-sha256 abd05a8a93fd2f543a112c344ca7f92cf95ce6d56e379cb80ba1e7d9dfb94882 \
  --output /fresh/path/TT_NORMAL.json
```

The checker is stdlib-only and all source paths are relative to `--repo-root`.
Only the final background report is pinned by its supplied digest.
