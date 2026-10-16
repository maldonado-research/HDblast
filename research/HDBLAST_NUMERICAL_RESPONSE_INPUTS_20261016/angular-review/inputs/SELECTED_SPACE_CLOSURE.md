# Closure of the selected quadratic angular feature space

This statement applies to the retained aligned state `Y=n1²-n2²`, actual first
master and gauge solutions, fixed background/potential/tension, and the selected
second initial prescription. It is not a restriction on arbitrary perturbations.

Let V be the real degree-two scalar SO(4) representation `(1,1)`, and let v
denote the matrix `diag(1,-1,0,0)`. The symmetric square is multiplicity-free:

    Sym² V = R00 + R11 + R22 + R20 + R02,
    dimensions = 1,9,25,5,5.

In the complete 287-dimensional metric–scalar angular section space the same
representations have multiplicities `5,8,8,1,1`. All are real type. Choose
compatible real equivariant identifications between their copies, and put
`z_alpha = projection_alpha(v tensor v)`. Every equivariant map from Sym² V
to the full section space sends `v tensor v` into

    E_v = direct_sum_alpha span(z_alpha) tensor R^(multiplicity_alpha).

All five vectors are nonzero. Explicit scalar representatives are
`1`, `Z=n1²+n2²-1/2`, and `F=Y²-(n1²+n2²)/2+1/12`; the two remaining
representatives are the nonzero chiral Weyl-derived TT features T+ and T−.
Their exact norms, Gram matrices and identities are checked separately.
Consequently `dim E_v=5+8+8+1+1=23`.

Each actual first metric, scalar and gauge field is an equivariant linear
image of the same v; differentiating in time or radius changes radial/time
coefficients, not the angular representation. Polarizing any quadratic
expression in these fields gives an equivariant linear map from Sym² V.
This includes covariant quadratic Einstein/scalar and harmonic residuals,
physical shell Hessians, embedding/normal pullbacks, and the quadratic terms
of the selected initial harmonic completion. Scalar potentials and the shell
tension are true scalar functions. Radial/time operations, contractions,
intrinsic derivatives, and spectral projections are natural under SO(4).
No extra angular state or orientation-dependent free datum is added.

The chosen physical second initial prescription uses precisely the l=0,2,4
projections of Y² and |DY|², their equivariant scalar gradients and the stated
radial multipliers. Thus it satisfies the same hypothesis. A different choice
of independent homogeneous second data need not satisfy it.

Every homogeneous bulk coefficient, relative-harmonic term, shell and gauge
boundary row, and initial linear completion is SO(4)-equivariant. On each
isotypic component, Schur's lemma makes it identity on the irreducible factor
times a matrix on the multiplicity space. It therefore preserves E_v. The
Dirichlet/Robin projections, reflected free propagator, causal Duhamel integral
and boundary lift commute with the same rotations and preserve E_v as well.

The accepted contraction equation `U=b+TU` consequently has `b in E_v` and
`T(E_v) subset E_v`. Every convergent iterate lies in this closed finite angular
subspace. Uniqueness identifies its limit with the accepted full-space formal
second response. Its existing scaled time and radial derivatives stay in E_v.
This proves selected-state closure without assuming that a numerical full
metric quadratic-source adapter has already been implemented.

The new exact angular matrix engine verifies the complete homogeneous bulk and
boundary implementation and the achieved explicit source adapters. The full
quadratic Ricci/reduced metric-source adapter remains a computational task;
the abstract naturality argument does not supply its coefficient values or
approximation errors. Neither nonzero TT feature alone proves nonzero actual
TT forcing. Both features must remain available unless that forcing cancels
by a separate exact calculation.

The sharp bound `|U(n)|_F <= sqrt(321/19) ||U||_H` may now be applied to this
selected formal response and its scaled time/radial derivatives. It also applies
to an error `U-V` only when the reconstructed candidate V lies in E_v. The full
287-section estimate retains factor17. Stabilizer symmetry alone gives larger
47- or 31-dimensional spaces and does not justify this reduction or its sharp
constant. No stronger parity reduction is used in this checkpoint.
