# Independent review of the physical second initial bounds

The quantitative bounds are accepted for the specified physical Cauchy second
jet, conditional on the authenticated actual background. I found no remaining
mathematical defect in the reductions, factorials, radial conversion, angular
factors, covariant connection terms or published rounded caps. This review
changes no producer file or frozen checkpoint.

The norm means the supremum of the absolute value of **each individual
component** in the adapted background orthonormal frame
`e_r=partial_y`, `e_i=rho^-1 E_i`. Each first spatial covariant derivative
index is also in that frame. It is not an ordinary coordinate-component
norm, a sum or Euclidean contraction of tensor components, an unrestricted
four-dimensional frame-rotation norm, a Sobolev or energy norm, or an evolved
second solution. Sphere-frame rotations preserving the distinguished radial
direction are allowed by the angular gradient estimates.

## Independent numerical conclusion

The standalone reviewer uses **unnormalized** derivative arrays `D_z^j`,
binomial Leibniz products and outward rational rounding to denominator
`2^96`. The producer uses normalized coefficients `D_z^j/j!` and a
`2^144` grid. No producer or other project program is imported.

I load the certified `rho_inv`, `p_inv` and `exp_minus_z` jets directly. The
producer instead recomputes the reciprocal jets and uses a slightly coarser
exponential value interval. Both procedures enclose the same actual functions;
the independent bounds are consequently a little sharper. All **17** published
rounded caps, including their directional subdivisions, pass strict rational
comparison. The principal bounds are:

| Field | Independent value bound, approximately | Published strict cap |
| --- | ---: | ---: |
| gamma2 | 1.0584259202820574e24 | 1.06e24 |
| K2 | 7.955430875643975e7 | 7.96e7 |
| phi2 | 5.120038469211462e33 | 5.13e33 |
| Pi2 | 1.99106307518752e11 | 2.00e11 |

| First covariant derivative | Independent component bound, approximately | Published strict cap |
| --- | ---: | ---: |
| nabla gamma2 | 4.242089672813432e29 | 4.25e29 |
| nabla K2 | 1.9594595981997273e13 | 1.96e13 |
| nabla phi2 | 3.050496305780613e39 | 3.06e39 |
| nabla Pi2 | 5.018291717701209e16 | 5.02e16 |

The exact reviewer fractions, input hashes and every directional comparison
are retained in `INDEPENDENT_INITIAL_NORMS_NORMAL.json`; floating values in
this table are explanatory only. Normal and optimized Python return identical
mathematical results.

## Actual input chain and derivative conventions

The reviewer pins these actual inputs independently of the producer's pin file:

- Static coefficient report: `3aa2c8d2e9ad652b72a34b4d8a6e6a6f14f72c32b8dcd1fbdc67556093265c22`.
- Exact bump report: `add927ba57117a967da7d6725821d95e894ef954405d503731b3055fc64377bd`.
- Explicit second-initial constraint derivation: `2cca9c9b87dfa5c2c21712ae81ea0c91d01b939b902e82d54dd62a6eae6bf0ed`.

The static report links to the retained actual potential report
`fd93a0de3d321ff5182115138eed2725f3c327083cf5148e3efc1aca3659d5c0`,
the analytic shell box and actual background theory note. I verify those
hash declarations, the registered polynomial potential, `delta`, `c`, the
whole-annulus domain, positive rho, negative p and the coefficient definitions.
In particular its `h` is `U'/p`; it is not `p_y/p`. The underlying smooth
actual background and its original existence proof remain retained hypotheses.
This review independently propagates the authenticated derivative enclosures;
it does not claim a new proof of that original background existence theorem.

For every input row j, the reviewer multiplies the normalized endpoints by
`j!` before performing any products. Its product and reciprocal recurrences are

\[
 (fg)^{(n)}=\sum_{j=0}^n\binom nj f^{(j)}g^{(n-j)},\qquad
 (1/f)^{(n)}=-f^{-1}\sum_{j=1}^n\binom nj f^{(j)}(1/f)^{(n-j)}.
\]

Differentiating an unnormalized array shifts it without another factor.
Exact polynomial and reciprocal derivative tests check these conventions.
The physical bump input is independently formed as `Q_j/ell^j`; only its
value interval is `[0,1]`, while higher derivatives are signed symmetric
intervals. The fourth derivative is retained for `phi2_y`.

The proper derivative is the operator `D_y=-b D_z`, `b=1/rho`. Applying
this operator twice automatically differentiates b. Since `b_z=a`, it gives

\[
 D_y^2f=b^2D_z^2f+abD_zf.
\]

It would be incorrect to use only `b^2D_z^2f`, or to change the first
derivative's sign. The final actual derivative budgets are `A:2`, all
`f_l:2`, all scalar profiles `zeta_l:1`, and `Pi2:2`. These are derivative
orders of the entire profiles, not orders of the bump alone. In particular
all requested first derivatives are explicitly represented.

All second Cauchy fields are second parameter derivatives in
`X=X0+epsilon X1+epsilon^2 X2/2`. No extra factorial is applied to those
fields or to the Hessian source when producing the component caps. Radial
Taylor factorials and the parameter expansion's factor two are different
conventions and have been kept separate.

## All-point reductions

Write `b=1/rho`, `j=v_y`, `E=e^-z`, and `d=3b/p`. The actual reconstruction
and background identities are

\[
 u_y=-pv/3-2au,\quad p_y=p(h-4a),\quad
 a_y=-b^2-p^2/3,\quad E_y=bE,\quad d_y=d(3a-h).
\]

Differentiating the first identity gives

\[
 u_{yy}=-pj/3+p(2a-h/3)v+(4a^2-2a_y)u.
\]

Substituting into the four original radial terms
`-24u_y^2-288auu_y-24uu_yy-(192a_y+384a^2+32p^2)u^2`
gives `8puj+8p(h+2a)uv-8p^2v^2/3+(-144a_y-32p^2)u^2`.
Combining the remaining original terms and substituting
`a_y=-b^2-p^2/3` yields exactly

\[
 A=-2j^2-8puj+8p(h+2a)uv+(16p^2-192b^2)u^2
   -(2U''+8p^2/3)v^2-b^2(12u_t^2+2v_t^2).
\]

This is an identity of the actual smooth functions, so its derivatives are
also identities. It is not a replacement of independently selected interval
values by differentiated interval relations.

For the momentum reduction, use
`u_t=Eu`, `u_yt=E(u_y+bu)`, `v_t=E(v-du)` and
`v_ty=E[j-d_yu-du_y+b(v-du)]`. Direct differentiation of
`F=b(18u u_t-v v_t)` gives

\[
 \frac{G-F_y}{p}=\frac{Eb}{p}\left[
 duj+(d(h-b)+10p)uv+(30a-6b)u^2+(2b-a)v^2\right].
\]

The `v v_y` terms cancel exactly. The same identities give
`F/3=(Eb/3)(18u^2-v^2+duv)`. The reviewer verifies all coefficients of
these three identities by a separate exact formal Laurent expansion,
including negative powers of p. Its zero residuals are all-point algebraic
checks, not generic-point tests.

## All angular and connection terms

The harmonics are `W0=1`, `W2=Z2=N-1/2`, and
`W4=H4=Y^2-N/2+1/12`. Since `0<=Y^2<=N^2<=N<=1`,

\[
 |W_0|\le1,\quad |W_2|\le1/2,\quad |W_4|\le7/12.
\]

The last bound uses the separate lower endpoint `-5/12` and upper endpoint
`7/12`. Further,

\[
 |\nabla N|^2=4N(1-N)\le1,\qquad
 |\nabla Y^2|^2=16Y^2(N-Y^2)\le4N^2\le4,
\]

so the gradient bounds are `0,1,5/2` for `W0,W2,W4`, and `2` for `Y^2`.
These are spherical gradient lengths and therefore bound every component
in every sphere orthonormal frame. They justify the precise weighted sums
used by both computations; no angular direction is omitted.

For `gamma2=H` times the angular projector, adapted-frame components are
`(0,H,H,H)`. The frame is parallel along y. In angular directions,
`nabla_ea er=a ea` and the radial part of `nabla_ea eb` is
`-a delta_ab er`. Thus the complete nonzero derivative component types are

\[
 (\nabla_{e_r}\gamma_2)_{ij}=H_y\delta_{ij},\quad
 (\nabla_{e_a}\gamma_2)_{bc}=\rho^{-1}(E_aH)\delta_{bc},\quad
 (\nabla_{e_a}\gamma_2)_{rb}=(\nabla_{e_a}\gamma_2)_{br}
                            =-aH\delta_{ab}.
\]

The component maximum is therefore bounded by the maximum of the radial,
angular tangential, and angular mixed estimates. The repeated mixed index
order is included with the same cap; a factor two would be appropriate for
some summed norms, not for this maximum of individual components.

For `K2=K(y)Y^2 gamma0`, metric compatibility gives
`nabla K2=d(KY^2) tensor gamma0`, with no projector connection term. Its
radial/angular bounds are `|K_y|` and `2|K|/rho_min`. The scalar `Pi2`
has the same two derivative estimates; scalar `phi2` uses the three harmonic
weights. Ordinary coordinate derivatives of angular metric components would
contain a further `2aH` term. That term cancels in radial covariant components
and must not be added to the reported adapted-frame bound.

## Scope and acceptance

The source and correction profiles are local differential expressions in the
exact compact flat bump and smooth background coefficients with nonzero p.
They and their first covariant derivatives extend by zero smoothly outside
the original annulus. Division by p is needed only near that annulus; no
uncertified remote background extension is introduced by the zero correction.

The reviewed caps concern the explicit physical Cauchy variables, with K a
covariant extrinsic curvature and Pi the unit-normal scalar velocity. A
conversion to coordinate velocities, spacetime metric derivatives, momentum
densities, a summed tensor norm, or a Sobolev norm needs its own factors and
cross terms. This review establishes no finite-amplitude epsilon interval,
gauge completion, second evolution, nonlinear remainder or invariant force.

The independent program accepts explicit `--background`, `--bump`,
`--initial-jet`, and `--producer` (or `--producer-report`) paths, and writes
only `--output`.
It uses only the Python standard library. Its normal and optimized outputs
are mathematically identical apart from the optimization flag. The analytic
reduction and connection proofs above, authenticated all-point derivative
inputs, and independent rational enclosures support acceptance of the
restricted quantitative first-spatial result.
