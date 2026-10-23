# Exact angular projection of geometric bulk sources

The retained parent basis is used without a change of normalization. Write
`Y=n0²−n1²`, `Z=n0²+n1²−1/2`, and
`F=Y²−(n0²+n1²)/2+1/12`. Scalar slots use `(1,Z,F)`, angular covectors use
`(DZ,DF)`, and angular covariant two-tensors use
`(kappa,Z*kappa,F*kappa,HessTF(Z),HessTF(F),Tplus,Tminus)`.
`Tplus` and `Tminus` are the two parent chiral Weyl features, including their
parent signs. Neither is removed or identified with the other.

The full tensor indices are `0=y`, `1=time`, and `2..5` are four redundant
ambient tangent slots on the unit sphere. They are fixed background-frame
covariant components. `tensor_fields` performs no time-index raising and no
coordinate/frame conversion. The 23 packing order is exactly the parent
`linear_operator.LABELS`: K,L,N,J,w each at degrees 0,2,4; v,u each at degrees
2,4; H at degrees 2,4; and Tplus,Tminus.

## Proof method

Every angular polynomial is reduced by the monic relation
`n3²=1−n0²−n1²−n2²`. The resulting polynomial dictionary has powers of n3 at
most one. Coefficients may belong to the rational numbers or the explicitly
implemented Laurent ring in independent ordinary metric, scalar, and
background jet symbols. `angular_generic.py` retains the parent's exact
polynomial and intrinsic tensor calculus and changes only coefficient
coercion in `clean`. The replay compares its complete primitive matrices to
the pinned parent engine.

For each target basis, `witnesses` selects linearly independent **monomial
coefficients** of tensor components. Exact rational elimination certifies ranks
3,2,7. These are coefficient witnesses, not evaluations at sphere points.
The inverse witness matrix supplies candidate projection coefficients.
`project_tensor` then reconstructs the tensor and checks the complete reduced
polynomial difference in every component. It accepts only zero remainders.
Consequently every accepted result is an identity at every sphere point and
for every value of the independent coefficient symbols in their domain.

The scalar odd polynomial n0 has zero in all three selected scalar witnesses
but is outside the scalar basis. Its deliberate rejection demonstrates why
witness rank alone is insufficient for closure. Independent symbolic tensor
coefficients test both TT orientations, their difference, and Hessian
normalizations. Controls delete both TT features, delete one orientation,
change a Hessian factor, replace full reconstruction with witness acceptance,
or exchange the final packing slots; each must fail for its stated reason.

## Complete source export and evaluation

`check_complete_source_projection.py` supplies all 42 metric two-jet entries
and all 6 scalar two-jet entries as independent Laurent generators to the
complete indexed geometric source engine. It projects the complete Ricci,
matter, Einstein, reduced Einstein, gauge-gradient, and forcing tensors,
including general off-shell corrections. Every output undergoes the full
polynomial reconstruction just described. It additionally compares QC,
Qphi, and Einstein matter amplitudes to the pinned predecessor adapters and
checks the complete Einstein join, gauge-gradient sign, general PC1 correction,
and metric/scalar forcing factors.

Metric amplitudes A,B,C,D,E,L,M and suffixes `y,t,yy,yt,tt` are **ordinary
physical derivatives**, not Taylor coefficients. They multiply Y,Y,DY,Y,DY,
Y*kappa,Hess(Y) in yy,y0,ya,00,0a,ab respectively. Scalar inputs are
`s,sy,st,syy,syt,stt` for eta=sY. Products of ordinary derivatives therefore
carry their ordinary binomial factors in the geometric jet engine. Source
values are second parameter derivatives along the affine first-field path;
they are twice the eps² Taylor coefficient. PDE forcing is
`−2 ell² rho² QEred` in metric rows and `+ell² rho² Qphi` in scalar rows.

The actual harmonic formulas for Qphi and QEred assume C1=0 as an identity
on a neighbourhood. This includes the derivative identity nabla(C1)=0; a
zero residual at one point is insufficient. The retained general alternatives
keep the off-shell C1 corrections. QC means the Hessian of the **full
metric-lowered** relative-gauge covector; QC_background records the distinct
background-lowered upper-gauge Hessian. Their difference is retained explicitly
in the geometric engine and vanishes on the actual first harmonic field.

`GEOMETRIC_COEFFICIENT_MAPS.json` stores each coefficient as an exact rational
linear combination of Laurent monomials in the named input symbols.
`source_coefficients.evaluate_sources` evaluates the collected maps with
Fraction, float, or interval arithmetic. `evaluate_bulk_forcing23` returns the
complete forcing in the registered packing order. Interval users must evaluate
the collected maps: substituting intervals before angular cancellation loses
dependencies and cannot reproduce an exact polynomial zero test.

The export closes the algebraic adapter task. Bounds on actual evolved first
jets, source enclosures, and the second-response PDE remain separate analytic
or numerical assertions; they are not inferred from algebraic closure.
