# Independent direct metric checker

`direct_metric_checker.py` uses only the Python standard library. It does not
import inherited physics code or any other agent's derivation. It constructs
the metric and its first two coordinate derivatives, evaluates Christoffel
symbols and their coordinate derivatives, and forms Ricci and the scalar
Laplacian directly. All arithmetic is exact `fractions.Fraction` arithmetic.
The local background jets are free: no background Einstein or scalar equation
is substituted.

The perturbation parameter is a polynomial truncated after epsilon squared.
The quadratic output is the second Taylor derivative, namely twice the
epsilon-squared coefficient. The general scalar linear operator uses the
epsilon coefficient without this factor of two.

Run the complete checks, including deliberate corruptions:

```sh
python direct_metric_checker.py --include-controls --output normal.json
python -O direct_metric_checker.py --include-controls --output optimized.json
```

Each run compares 30 quadratic coefficients and 24 linear coefficients
exactly. The quadratic calculation uses three arbitrary background/field
two-jets and five rational angular points per jet. The linear calculation uses
the same three background jets, five independently assigned amplitude jets
`fyy`, `fyt`, `ftt`, `H`, and `w`, and three rational angular points per jet.
The off-diagonal base metric amplitude `fyt` is nonzero. Its candidate operator
is evaluated separately with two-dimensional covariant derivatives and
compared with the full five-dimensional Christoffel calculation.

There are also 15 exact checks of the explicit component mixed rows and the
forced relations in note equations (16), (17), and (18). After substituting
`N = d + G - 2K` and `T = G - K`, the temporal identity is checked as
`2(E_t + j_t) = LHS(17) - RHS(17)`. The radial identity is checked as
`E_y + j_y = RHS(18) - p w`, and setting `w = RHS(18)/p` is verified to make
that sourced mixed row vanish. This checks the mandatory completion and its
mixed constraints; it does not prove reduction of the remaining rows to a
complete master system or regularity where `p = 0`.

Angular scalar coefficients are extracted by an overdetermined exact fit to
`1`, `H = sin(theta)^2 - 1/2`, and `Q`. Mixed components are fitted to `DH` and
`DQ`. The quadratic angular trace-free tensor is fitted to the full primitive
pair `Y (Hess Y)^TF` and `(DY DY)^TF`; its scalar harmonic projection is then
calculated. The linear angular trace-free tensor is fitted directly to
`(Hess Q)^TF`. Any inconsistency or rank deficiency raises an exception.

The normalized sphere integrals give

```text
<Q,Q> = 19/720
<(Hess Q)^TF,(Hess Q)^TF> = 133/15
[Y (Hess Y)^TF]_(Hess Q)^TF = 1/6
[(DY DY)^TF]_(Hess Q)^TF = 1/3
[|DY|^2]_Q = -4
```

For the last three identities, the checker uses `Delta Y = -8 Y` and
`Delta Q = -24 Q`, the latter also checked at every angular point. The
gradient-square coefficient is `8 - 24/2 = -4`. Two integrations by parts give
the contraction of `(DY DY)^TF` with `(Hess Q)^TF` as `112 <Q,Q>`, while the
Hessian norm is `336 <Q,Q>`. The identity
`(Hess Q)^TF = 2Y(Hess Y)^TF + 2(DY DY)^TF - (Hess H)^TF/2`, together with
orthogonality of the `l=2` and `l=4` scalar-derived tensors, supplies the `1/6`
coefficient. The measure is normalized: `x = sin(theta)^2` is uniform on
`[0,1]`, and phi is uniform on the circle.

There are nine negative controls. They reverse the Ricci sign, drop Ricci
connection products, drop the scalar gradient square, change the potential
normalization, omit the second-Taylor factor, corrupt the angular projection,
or omit one of three terms of the linear operator. Every control is required
to fail. For a standalone deliberately corrupted run, for example:

```sh
python -O direct_metric_checker.py --control omit-d2-factor --output wrong-half.json
```

This command must return exit status 1 and a JSON report with `status: FAIL`.
All substantive checks use explicit exceptions, so Python optimization cannot
disable them.

`receipt.json` freezes the script SHA-256, the final four executed commands and
their observed exit codes, and SHA-256 hashes of their JSON reports. It also
records that the ordinary and optimized reports are identical apart from the
optimization indicator. The two explicit corrupted-command reports are
`wrong-half-optimized.json` and `wrong-projection.json`.

The reports establish exact equality at the stated free local jets and
angular points. They are reproducible independent checks of the displayed
formulas; finite local samples are not a replacement for an algebraic proof
for all possible jets.
