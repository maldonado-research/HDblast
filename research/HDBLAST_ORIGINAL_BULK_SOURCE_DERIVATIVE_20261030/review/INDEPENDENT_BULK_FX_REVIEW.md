# Independent internal review of the original bulk-source time derivative

30 October 2026. This is a mathematical and numerical review of the original
registered 23-amplitude forcing, not an outside specialist assessment, novelty
opinion or statement of physical evidence. The actual source is not smoothed
or replaced.

**Assessment: PASS**, conditional on the retained actual background, master,
gauge, harmonic-neighbourhood and initial/boundary premises. The complete
original forcing has the positive H-norm bound

```
sup ||F_x||_H <905405545954284 <9.055e14.
```

This is an upper bound only. It supplies no nonzero sign, physical force,
relative smallness, second-response C2 target-error theorem by itself, exact
nonlinear family or cosmological conclusion.

## Physical-coordinate derivatives

Set `L=1/(ell rho)`, with `x=t/ell`, `dy=-ell rho dr`. The background is
time independent except for the explicit hyperbolic factors. The radial
identity `rho_r=-ell rho^2 a` gives `L_r=a`. Thus, for a frame amplitude j,

```
j_y=-L j_r,
j_yy=L^2 j_rr+a L j_r,
j_yyt=(L^2 j_xrr+a L j_xr)/ell,
j_ytt=-L j_xxr/ell^2,
j_ttt=j_xxx/ell^3.
```

The `a L j_xr` term must be retained. These are ordinary derivatives of the
registered component amplitudes, not covariant derivatives of a source
tensor. The maps already encode the frame connection separately.

The independent identity checker applies `D_y=-L(r)D_r` and `D_t=D_x/ell`
directly to ten bivariate polynomial amplitudes with nonconstant L. It
checks all three resulting third-jet identities at the polynomial level,
not merely at one event. It also verifies the factorial conversion from
normalized Taylor coefficients for all mixed raw degrees through four.

The original maps involve field suffixes `'',y,t,yy,yt,tt`. Ordinary time
differentiation introduces exactly `yyt,ytt,ttt` for each of the seven metric
amplitudes and the scalar: 24 new entries. It introduces no pure `yyy`
entry. The frame amplitudes contain one raw master or gauge derivative, so
these new entries can involve raw degree-four jets. Existing background
coefficient jets through degree three suffice; no fourth background jet is
used in the third-field construction.

## Complete source chain rule and normalization

The review authenticates the collected source map with SHA256
`34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615`.
Its 2,233 original complete forcing terms yield 5,248 collected derivative
monomials across all 23 rows. The derivative substitutions are

```
D_t(field_suffix)=field_suffix_followed_by_t,
D_t c=c theta,
D_t theta=1-theta^2,
D_t(rho,a,ay,p,py,Vp,Vpp,Vppp,ell)=0,
F_x=ell D_t F.
```

Here the field suffixes are written in the canonical y-then-t order. The
other background variables are static because the background scalar and
warp factors depend only on y. Negative Laurent powers, including inverse
c and rho, are differentiated by their ordinary product rule.

The identity checker independently builds the derivative Laurent maps and
compares all 23 rows with direct dual-number evaluations at four generic
rational fixtures: 92 row identities. The fixture dual directions are
initialized separately from the symbolic differentiation code. Every
canonical ordinary-time derivative-row digest agrees with the producer.

It also checks the complete original and differentiated normalization
identities in all 20 metric and three scalar slots:

```
metric_forcing=-2 ell^2 rho^2 QEred,
scalar_forcing=ell^2 rho^2 Qphi.
```

These factors are already in the stored forcing maps. Multiplying by them
again would be wrong. Only the outer ell converts ordinary physical time
differentiation to the requested scaled derivative.

The actual QEred/Qphi specialization assumes `C1=0` as an identity on a
neighbourhood, including the derivatives needed here. Independent pointwise
choices from interval boxes would not prove that identity. It remains an
inherited actual-field premise of this review.

## Independent third-jet composition and numerical norm

The numeric checker imports no producer code. It authenticates the new
actual master/gauge reports, the retained master report through degree three,
the complete source map, and the independent identity checker before using
its symbolic operations.

The new composition is recomputed using **ordinary mixed derivatives**:

```
D_x^i D_r^j(fg)=sum_ab binom(i,a)binom(j,b) f_(a,b) g_(i-a,j-b).
```

An ordinary-derivative reciprocal recurrence follows from differentiating
`f(1/f)=1`. This differs from the producer's normalized Taylor convolution.
Raw derivative shifts use the authenticated degree-four bounds directly;
only coefficient jets through degree three enter. Reconstructing the
registered frame amplitudes by this independent calculus puts all 24 third
physical-derivative intervals inside the producer's outward-rounded boxes.
In particular, this validates the full raw-master/gauge-to-third-field
composition, rather than merely checking the producer's final H-norm sum.

The coefficient ledger for this computation uses the whole actual
background/time box. The lower ordinary field jets use a hull of all 81
retained source-cover cells. The independent check verifies the full clipped
cover geometry and input-symbol closure. These hulls may widen bounds and
do not require one independent interval selection to be an actual solution.
The background and lower-field interval certificates retain their existing
authenticated analytic premises; this review does not replace their proofs.

The independently differentiated complete maps are evaluated with exact
rational interval arithmetic on those lower boxes and the independently
composed third boxes. All 23 resulting intervals are contained in the
producer's dyadic enclosures. The positive Gram weights are independently
packed from the probability-normalized scalar/vector/tensor basis, including
both TT rows, factor two on mixed metric slots and factor three on J. The
reported full squared norm is reconstructed exactly; its strict integer
square-root ceiling is 905405545954284.

## Reproduction and controls

`check_bulk_derivative_identities.py` and
`check_bulk_derivative_numeric.py` use only the standard library. Their
normal and optimized reports are byte-identical. The identity audit includes
seven wrong-formula controls, and the numeric audit includes four, covering
the radial scale product term, physical derivative factors, Taylor factorials,
hyperbolic coefficient derivatives, scaled-time ell and norm packing.

Input roles and output metadata are independent of local file basenames,
so the same reports can replay in a checkout or CI. Exact report digests and
source hashes are recorded separately with the numerical artifacts.

The result closes the original full-source derivative input, subject to its
stated assumptions. A C2 error gate must still combine it with the correct
candidate residual derivative, higher initial errors and homogeneous shell
error conditions. Coordinate error control then requires the declared frame
and connection conversion before becoming a natural-equation statement.
