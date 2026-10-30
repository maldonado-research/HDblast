# The complete original bulk-source time derivative

For the fixed registered actual first field, the complete reduced 23-amplitude
bulk forcing satisfies the continuous bound

```
sup_D ||F_x||_H <905405545954284 <9.055e14,
D={0<=x<=3, r>=0, x+r<=8},   x=t/ell, ell=1e-6.
```

The norm is the probability-normalized angular S3 L2 norm of the positive
metric-plus-scalar bundle, including both transverse-traceless orientations
and the factor two on mixed metric components. This is an absolute upper
bound in the registered scaled coordinate and frame. It supplies neither a
sign nor a nonzero value, a physical force, an invariant observable, or a
relative-accuracy claim.

This is the derivative of the original actual source F. It is not the
derivative of the one-sided averaged source introduced in the separate
numerical-candidate checkpoint, and it is not a derivative of the response
quantity X(KF+LG). The actual source, first field, model and shell data are
unchanged.

## What is differentiated

The authenticated complete maps have SHA256
`34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615`.
Their first-metric frame-amplitude symbols are A,B,C,D,E,L,M, with ordinary
physical y/t derivatives through total order two, and scalar symbols s with
the same suffixes. The background time dependence is represented by
`c=cosh(t)` and `theta=tanh(t)`. All other stored background factors are
static in time. At fixed r, hence fixed physical y, the complete chain rule is

```
F_x=ell [sum_j F_(m_j) (m_j)_t
          +sum_j F_(s_j) (s_j)_t
          +F_c c theta +F_theta(1-theta^2)].
```

The implementation differentiates and collects every Laurent monomial before
interval evaluation. The temporal suffix map is
`''->t`, `y->yt`, `t->tt`, `yy->yyt`, `yt->ytt`, `tt->ttt`.
Every output row is retained. In particular, the new inputs are the ordinary
third field derivatives yyt,ytt,ttt; pure yyy is not needed for this time
derivative.

The stored physical metric forcing is already `-2 ell^2 rho^2 QEred` and
the scalar forcing is already `+ell^2 rho^2 Qphi`. Since rho and ell are
time independent, differentiating those maps gives the same factors times
the differentiated reduced/natural source maps. The executable checks the
resulting Laurent polynomial identities in all 20 metric rows and all three
scalar rows. There is exactly one additional outer ell converting the
ordinary physical t derivative to the scaled x derivative.

The retained actual-mode maps use the harmonic identity C1=0 on a
neighbourhood, including its derivatives. This is a stated conditional
premise of the actual first field, not a condition inferred from selecting
independent values inside a pointwise interval box.

## New actual field three-jets

The same registered amplitude reconstruction is used:

```
u=R vartheta W,
v=B vartheta (W_r/ell+k W),       vartheta=cosh(ell x)^(-3/2),
A=-4u-2 f_r/(ell rho),          Bmetric=f_x/(ell rho)+b_r/ell,
C=f/(rho c)-c d_r/ell,          D=-2(u+a f+b_x/ell),
E=c d_x/ell-b/c,                L=2(u+a f+theta b),
M=2d,                          s=v+p f.
```

Here f,b,d are the actual raw gauge amplitudes, W is the actual first master,
and the static coefficient B in v is distinct from Bmetric. The same
background ODE jet evaluator and temporal coefficient jets are retained.
The actual raw master and all three raw gauge amplitudes now have continuous
scaled mixed derivative caps through total order four. These are ordinary
derivative bounds, not normalized Taylor coefficients.

For a normalized raw Taylor jet, the coefficient of h_x^i h_r^j is the
ordinary derivative divided by i!j!. All product and inverse operations are
performed in that jet algebra, with outward rational intervals. The raw
master/gauge degree-four coefficients are needed because the reconstruction
already contains one raw derivative. Each static/temporal coefficient needs
only degree three: no derivative is applied to those coefficient jets inside
the displayed amplitude formula. Degree-four background slots can therefore
be set to formal zero without affecting any retained output through degree
three. This does not assert that actual fourth background derivatives vanish.

Let `Lrad=1/(ell rho)`. The static radial transformation obeys
`partial_y=-Lrad partial_r`, `partial_t=ell^-1 partial_x` and
`(Lrad)_r=a`. Hence for every reconstructed frame amplitude j,

```
j_yyt=(Lrad^2 j_xrr+a Lrad j_xr)/ell,
j_ytt=-Lrad j_xxr/ell^2,
j_ttt=j_xxx/ell^3.
```

The `a Lrad` product term is retained. These are ordinary derivatives of the
frame-amplitude functions, not covariant derivatives of a tensor. Tensor
connection terms have already been included in the authenticated source maps.
The bounds concern the physical side of the shell, with the retained
one-sided limits; no C4 regularity across a doubled-shell reflection is
asserted.

## Retaining the actual lower jets and covering the domain

The previously authenticated actual first-field adapter supplies all field
jets through order two. Its actual solution bounds include the retained C3
master/gauge approximation errors and the exact reflected comparison
expressions; the free comparisons alone are never substituted for the
actual first solution. The new actual fourth-order caps supply the higher
inputs separately. Their independent interval choices only widen a true
derivative enclosure and are not differentiated as arbitrary value functions.

The lower adapter is evaluated on all rectangles

```
[i/d,(i+1)/d] x [j/d,(j+1)/d],
0<=i<3d, 0<=j<8d, i+j<8d,
```

intersected with x+r<=8, with d=2. Their union is exactly D, including the
characteristic boundary, and there are `(39d^2+3d)/2=81` rectangles. The
global hull of those actual lower-jet boxes and the new whole-domain
third-jet boxes is substituted into the collected derivative maps. Taking
that hull loses correlations and sharpness, but preserves the all-point
upper bound. No small or diagonal cells are discarded.

The report records the actual raw caps, the normalized coefficient-jet
ledger, global lower/third input boxes, canonical differentiated-row hashes,
all 23 output intervals, positive Gram weights, and exact squared H-norm
upper bound. Its final integer upper bound is obtained by an exact rational
square-root ceiling, without floating-point decisions.

## Checks and scope

Before project imports the script authenticates the numerical source files,
their predecessor/input manifests, every closed dependency, the complete
source map and the master/gauge fourth-jet reports. All are read back after
the calculation. It checks all 23 differentiated forcing normalizations,
69 independent dual-number Laurent chain-rule fixtures, normalized Taylor
factorials, the radial product term, the gap-free cover and the full positive
23-component packing. Wrong outer-time scaling, a missing radial product
term and a third Taylor coefficient mistaken for an ordinary derivative are
rejected explicitly.

The independent review reconstructs the differentiated maps and higher-jet
composition using separate algebra. Normal and optimized execution produce
byte-identical reports, including when input reports are given different
filenames.

This closes the stated quantitative complete actual F_x input, conditional
on the authenticated smooth background/master/gauge, differentiated initial
and physical boundary equations, harmonic-neighbourhood and all-point
projection premises. A subsequent second-response C2 error estimate must
still supply its own differentiated residual, initial acceleration and
boundary ledger. No exact nonlinear family, finite-amplitude remainder,
physical detection, cosmological prediction or higher-dimensional Big Bang
origin follows from this derivative bound.
