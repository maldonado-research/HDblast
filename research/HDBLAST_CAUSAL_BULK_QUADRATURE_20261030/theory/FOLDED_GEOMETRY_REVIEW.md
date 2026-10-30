# Independent review of folded causal shell quadrature

The selected output is `(a,r)=(9/8,0)`. The correction obeys

```
C = K f + L(M PN C),       f = (La-L0) Bs + La(W0-Bs) + F.
```

Here `M` is the same fixed actual Robin matrix as in the prior reconstruction.
It is constant in shell time and has zero Dirichlet rows and columns. The
same fixed affine source `F` remains in the bulk forcing. A uniform genuine
component enclosure of `F` can be integrated as a symmetric error contribution;
the response quantity `X(KF+LG)` is not a pointwise source bound.

## Weighted source moments

The complete physical source is formed before reflecting its output. At the
shell, physical positive and reflected negative source regions give

```
h(t) = PN integral_{s,y>=0; s+y<=t} f(s,y) ds dy,
Z(t) = h(t) - integral_1^t M Z(v) dv.
```

The inherited zero-wedge premise removes `s+y<=1`; small source caps alone
would not justify removing this region. For every nonnegative integer `j`,
Fubini over the Volterra history gives

```
J^j h(a) = PN integral_{Omega_a} (a-s-y)^j/j! f(s,y) ds dy.
Z(a) = sum_{j>=0} (-M)^j J^j h(a).
```

The change of variables `alpha=s+y`, `theta=(y-s)/(s+y)` has positive Jacobian
`alpha/2`. Thus `1<=alpha<=a` and `-1<=theta<=1`. For a rectangular panel, its
exact positive weight is

```
(theta_hi-theta_lo)/(2 j!) *
 integral_{alpha_lo}^{alpha_hi} alpha (a-alpha)^j d alpha.
```

Multiplying this positive mass by an interval containing the source on the
whole panel encloses the panel integral. A midpoint source sample does not
supply that premise. The total mass is

```
A_j = d^(j+1)/j! * (a/(j+1) - d/(j+2)),   d=a-1=1/8.
```

In particular `A_0=17/128`. If `||f||_H<=R` and `||M||_H<=m`, truncation after
degree `N` has an H-norm tail bounded by

```
R m^(N+1) A_(N+1) / (1 - m d/(N+2)).
```

Indeed, `A_(j+1)/A_j <= d/(j+1)`; the first ratio after the first omitted term
has denominator `N+2`. The positive denominator must be checked. This argument
uses the complete residual cap `R`, including the affine source.

## Shell derivatives

Define the physical null-edge integral

```
p(a) = integral_0^a f(s,a-s) ds
     = (a/2) integral_-1^1 f(a(1-theta)/2,a(1+theta)/2) d theta.
```

Differentiating the physical reflected wave integral first gives
`(Kf)_x=PN p` and `(Kf)_r=PD p` at the shell. The inward-sign characteristic
lift contributes `-MZ` to the time derivative and `+MZ` to the radial derivative.
Consequently

```
PN C_x = PN p - M Z,     PN C_r = M Z,
PD C   = PD C_x = 0,     PD C_r = PD p.
```

The physical null edge has length parameter `a`, not `d`. A zero Dirichlet
trace does not set its radial derivative to zero.

## Exact independent fixtures

`check_folded_geometry.py` imports only standard-library modules and contains
27 named groups and seven wrong-formula controls. Its four-component fixture
has one Dirichlet output, three Robin outputs, a nonzero matrix `M` with
`M^3=0` and `M^2!=0`, and a source coupling that does not commute with parity.

For arbitrary physical-side polynomial sources, it compares five weighted
moments computed by expansion in `(s,y)` against exact rectangular integration
in `(alpha,theta)`. It also verifies continuous whole-panel interval ranges
against those exact integrals, over 64 panels per component and moment.

An independent polynomial Volterra reconstruction is compared with the folded
moments. The physical reflected kernel is integrated directly over its
physical positive and negative domains, then combined with the Robin lift.
Its shell derivatives are recovered by exact polynomial interpolation and
checked against the formulas above. No HDBLAST point evaluator or numerical
quadrature checker is imported. These manufactured sources need not be smooth
across `s+y=1`; the fixture validates integral geometry and algebra, without
proving the scientific source's regularity or support premises.

The wrong-formula controls reject omission of the Jacobian, moving the source
coupling through parity, the wrong lift sign, omission of the second Robin
moment, dropping the Dirichlet normal derivative, substituting `d` for `a` on
the edge, and retaining only the first omitted tail term.

## Concrete implementation review

The current `check_bulk_quadrature.py` implements the weights, matrix order,
tail, first derivatives, and intersection with prior analytic enclosures
consistently with this derivation. Its interval powers repeatedly use the
actual matrix box; repeated interval dependence broadens the range and does
not freeze the matrix to a midpoint. Projecting the final moment sum is
equivalent to projecting each source moment because `M PN=PN M=M`.

An enclosing Cartesian `(s,y)` rectangle for an `(alpha,theta)` panel can extend
beyond the causal triangle. This is harmless only when the coefficient and
profile evaluators enclose that larger rectangle in their certified domain.
The default fine theta mesh stays within the inherited bulk domain. A very
coarse mesh can extend outside it and must be rejected or subdivided by the
evaluator. The reviewed evaluator explicitly rejects a coefficient rectangle
outside the retained footprint.

The `interior_evaluator.py` source was also reviewed separately. Its prefix
recurrence is the normalized-primitive identity

```
U_j(b+z) = sum_{h=1}^j z^(j-h)/(j-h)! U_h(b)
           + integral_b^(b+z) (b+z-v)^(j-1)/(j-1)! u(v) dv.
```

Whole-panel nominal profile ranges bound the last integral by positive
weights. The uniform displacement smoothing error gives a primitive remainder
`alpha_i max(0,z-u_collar)^j/j!`; the velocity primitives are exact compact
kernel convolutions of P1 hinges. The displacement derivative uses its
separate registered smoothing bound. Negative primitives carry `(-1)^j P`;
the displacement derivative has the opposite reflection parity.

The complete frozen ridge functions use these normalized primitives. Their
matrix action occurs before output reflection. The interior reflected-kernel
branches agree with direct physical positive/reflected source integration;
inside the flat collar they reduce to a single branch across zero. The free
wave and kernel first derivatives retain separate correlated-argument
intervals outwardly, without treating the actual coefficients as constants.

For the boundary lift `Z-h`, homogeneous terms one through three are evaluated
without subtracting two large traces. Its remaining homogeneous H-norm tail is
`Hcap (m d)^4/4!/(1-m d/5)`. The actual affine source contributes at most
`d Gcap_i` componentwise from its first integral and
`m d^2 GH/2/(1-m d/3)` in H-norm from all later terms. The load is then
`M Z+G`, so lift time and radial derivatives have the required opposite signs.

Finally, the residual is assembled in cancellation-free form
`(La-L0)Bs+La(W0-Bs)+F`, with interval actual coefficient assembly on every
physical coefficient box. These observations found no formula, source-unit,
reflection-order or outward-enclosure gap. They do not replace authenticated
execution against the archived scientific input files.

For a row that is identically zero in the actual `M` box, every omitted Robin
history term is identically zero, since its power is at least one. Its history
tail component radius may therefore be set to zero. This does not remove the
zeroth source moment from that row's trace.

This review supports conditional reduced-system quadrature. It establishes
neither a full nonlinear spacetime nor a physical observable.
