# What the complete differentiated bulk source still needs

This is a structural audit of the **complete** registered bulk-forcing maps.
It does not provide a numerical upper bound on the actual `F_x`, a sign for
that derivative, or a second-derivative error theorem. It identifies a
specific high-jet input which is not removed merely by substituting the
principal wave-trace relations.

## Inputs and derivative convention

The authenticated `GEOMETRIC_COEFFICIENT_MAPS.json` contains all metric
forcing sectors `yy,y0,00,ya,0a,ab` and the scalar forcing, packed into 23
components. The stored metric forcing already has normalization
`−2 ell^2 rho^2 QEred`, and the scalar forcing has
`ell^2 rho^2 Qphi`. Each is a second parameter derivative, not an `eps^2`
Taylor coefficient. The source map SHA256 is
`34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615`.
Its actual-mode QEred and Qphi formulas assume the actual first harmonic
condition `C1=0` as an identity on a neighbourhood, including its derivative;
the general off-shell alternatives are distinct stored maps. This inherited
premise is retained here, not inferred from a point residual test.

Here `A,B,C,D,E,L,M` denote the **first** metric's frame amplitudes, not the
23 second-response labels. In particular `A=q_yy/Y`, `B=q_y0/Y`, and
`D=q_00/Y`. The background-frame signature is `(+,−,+,+,+)` and
`Y=n1^2−n2^2`. All suffix derivatives in the map are ordinary physical y,t
derivatives of frame amplitudes. The target derivative is the ordinary
scaled-time amplitude derivative, `F_x=ell partial_t F`, at fixed r (hence
fixed y). It is not a covariant derivative of the source tensor.

For complete coefficient maps `F=F(m,s,b)`, the ordinary chain rule is

```
F_x=ell [sum_j F_(m_j) (m_j)_t
         +sum_j F_(s_j) (s_j)_t
         +F_c c theta+F_theta(1-theta^2)].
```

Other background variables in this representation are static in time.
The highest new field inputs are ordinary third metric/scalar derivatives
with at least one t derivative. Lower field jets, static radial factors and
the known hyperbolic derivatives supply the remaining terms. A value cap for
F cannot be differentiated to obtain this expression's magnitude.

The actual first-field evaluator `numeric/actual_first_jets.py` constructs
its frame two-jets from C3 master/gauge data. Its source includes one raw
gauge derivative in several metric amplitudes and one master derivative in
the scalar amplitude. Thus the new ordinary field three-jets can involve
fourth raw gauge/master derivatives. Qualitative smoothness ensures these
derivatives exist under the retained premises, but does not provide the
missing quantitative constants.

## Null highest-jet test

At an event, hold the background and all raw gauge jets of order at most
three fixed. Consider a normalized fourth-jet variation of one raw gauge
amplitude with physical y,t covector

```
zeta_y=1,   zeta_t=sign*rho,   sign=+1 or−1,
delta partial_y^i partial_t^j xi = zeta_y^i zeta_t^j tau,
i+j=4.
```

The covector is null for the y,t principal metric. Every pair trace of this
symmetric fourth tensor vanishes:

```
rho^(-2) delta xi_(y^i t^(j+2))
   −delta xi_(y^(i+2) t^j)=0,     i+j=2.
```

These are the highest-order homogeneous prolonged gauge-wave relations.
Since lower jets are fixed, lower-coefficient prolongation terms do not
change in this highest-jet test. This does **not** construct a second actual
solution with the same global initial/boundary prescription. It tests the
information supplied by local jet values and their wave-trace relations.

The full lowered gauge reduction removes the other highest Ricci terms.
An independent indexed normal-frame calculation gives the remaining metric
principal term

```
delta partial_t QEred_AB
   =q^(CD) delta(partial_t nabla_C nabla_D q_AB)
     at this highest-jet order.
```

The null contraction of the low first metric is

```
J_sign=A+D−2 sign B.
```

It is not the Lorentz trace of q, and harmonic gauge does not generally set
it to zero. For the radial contravariant gauge amplitude f, the full
`metric_forcing` maps give exactly

```
delta(F_yy)_x =−4 sign ell^3 rho^3 J_sign Y^2 tau,
delta(F_y0)_x =−2      ell^3 rho^3 J_sign Y^2 tau.
```

All other metric sectors and the scalar forcing have zero fourth-f
sensitivity. The scalar expansion of `Y^2` has coefficients `(1/6,1/2,1)`
in the registered degree-zero/two/four basis. The displayed result is
therefore a statement about the collected **total** forcing coefficients,
after all Ricci, gauge-gradient and matter contributions have been included;
it is not a selected potential or Ricci summand.

## Complete directional sensitivity audit

The checker differentiates the stored Laurent maps exactly, using its own
standard-library rational polynomial operations. It checks all 23 output
components for both signs and the following input directions. In this table
the entries are the nonzero harmonic coefficient vectors; every omitted
output is checked to vanish. Each row is multiplied by `J_sign*tau`.

| Highest input | Output | Coefficient vector |
| --- | --- | --- |
| radial gauge f, total order four | yy | `sign ell^3 rho^3 (−2/3,−2,−4)` |
| radial gauge f, total order four | y0 | `ell^3 rho^3 (−1/3,−1,−2)` |
| time gauge b, total order four | y0 | `sign ell^3 rho^4 (1/3,1,2)` |
| time gauge b, total order four | 00 | `ell^3 rho^4 (2/3,2,4)` |
| sphere-gradient gauge d, total order four | ya | `sign ell^3 rho^4 c (−1/2,−1)` |
| sphere-gradient gauge d, total order four | 0a | `ell^3 rho^4 c (−1/2,−1)` |
| first scalar s, total order three | scalar | `sign ell^3 rho^3 (−1/3,−1,−2)` |

The scalar row uses a third-order null variation with all lower scalar jets
fixed. The scalar first field contains a first master derivative, so this
row can require fourth master jets after composing the master-to-scalar
adapter. It is not a fourth-gauge sensitivity of the scalar source.

The source maps are affine in their highest second field jets, which the
checker verifies before differentiating. It checks 216 exact component and
wave-trace identities and four wrong-formula controls, normally and with
optimization. The independent five-index check verifies eight generic frame
fixtures, null quartic traces and three wrong-formula controls without
importing the coefficient maps.

## A nonzero coefficient for the retained actual first field

This dependence is not merely a generic algebraic possibility. The existing
actual first-field two-jet report, SHA256
`0aaa7b2fe30f546db92ef8b26a7224ffa7f6812d3780c331661395aa1991d3b8`,
encloses the evolved-packet cell

```
1/4<=x<=1025/4096,   5/4<=r<=5121/4096.
```

Its actual `A,B,D` intervals prove continuously on that cell

```
−4.859e−4 < J_plus < −2.233e−4.
```

Using the same report's positive rho interval and `ell=1e−6`, the complete
mean yy coefficient's sensitivity to the normalized fourth-f null jet obeys

```
3.265e−16 < −(2/3) ell^3 rho^3 J_plus < 7.104e−16.
```

This is a nonzero **derivative-adapter coefficient** on an actual enclosed
first-field cell. It is not the value or sign of actual `F_x`: the unknown
highest jets and every lower term still enter that derivative. Nor is it a
force or an observable. Other high-jet channels or special relations of the
fixed actual solution could cancel in the evaluated total derivative; that
would require additional proved information, not deleting the coefficient.
The checker authenticates this inherited report and performs new exact
interval arithmetic; it does not rerun its C3 field certificate.

## Consequence for the next certificate

The prolonged wave traces leave two fourth-characteristic directions for
each gauge amplitude. Tracing them away does not control both null
directions. A practical route is therefore to certify the required fourth
characteristic jets from differentiated master/gauge equations and their
actual initial/boundary data, then use the complete chain rule. An alternative
must prove a data-specific identity controlling the displayed nonzero
contractions or the complete combined high-jet term.

This audit proves neither that C4 bounds are logically unavoidable under
every possible reformulation nor that the actual source derivative is
nonzero. It rules out the simple universal claim that principal wave-trace
substitution alone cancels every new highest-jet term in these maps. The
complete actual `F_x` bound remains open. Consequently the quantitative C2
error gate, constraint-derivative estimate and full natural Einstein residual
are not closed by this structural result.

The maps, actual first field and all physical background/source assumptions
remain conditional. No nonlinear spacetime family or higher-dimensional
origin of the Big Bang is established.
