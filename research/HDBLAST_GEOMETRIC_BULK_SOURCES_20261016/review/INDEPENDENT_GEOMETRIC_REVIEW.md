# Independent internal geometric-source review

16 October 2026. **Accepted within the stated quadratic-source and intrinsic-
coefficient scope. No blocking mathematical defect found.** This is an internal
independent mathematical review, not external peer review. The exact accepted
files and completed receipts are bound in REVIEW_RECEIPT.json.

The main result is a complete, executable, all-point angular adapter for the
Ricci, natural Einstein, registered reduced Einstein, scalar and relative-
harmonic quadratic sources. It uses the actual degree-two sphere harmonic,
all seven metric amplitudes and their ordinary physical two-jets, the scalar
two-jet, the correct stress/potential normalization, and both chiral TT basis
features. It supplies coefficient values rather than only an abstract closure
argument. This review also accepts the separately scoped four-zero theorem
for the intrinsic material-shell scalar's quadratic mean.

## Mathematical derivation and conventions

INDEPENDENT_DERIVATION.md derives the source formulas independently from the
exact inverse-metric and connection-difference identities. The Ricci Hessian
is proved for arbitrary smooth metric variations without a harmonic condition
or a background field equation. Covariant second derivatives include the
connection acting on the derivative slot. Time-index raising uses the Lorentz
sign. Every source is a second parameter derivative, not an epsilon-squared
coefficient. The scalar kinetic factor and potential cross factor are both
two, with the retained Einstein residual E=Ric-dphi*dphi-(2/3)Vg.

The implementation correctly distinguishes the second variation of the full
metric-lowered gauge covector from background lowering of the upper-vector
Hessian. The former is universally -2q^cd P_bcd. General reduced Einstein
jets include +2P^c_ab(C1)_c. The registered harmonic source is obtained only
with the actual first-field gauge identity. Its polynomial interval extension
may be evaluated on independent jet boxes because those boxes enclose the
actual correlated harmonic field; their arbitrary points need not satisfy
the gauge equation. A pointwise C1=0 witness alone does not justify dropping
derivatives of C1 when comparing the two lowering conventions.

The forcing convention is also consistent with the inherited PDE. Natural
second equations have residual forcing -Q. The metric reduced operator has
principal part -Box/2 and the scalar operator has principal part +Box. Since
rho^2 Box has principal part -partial_t^2+partial_z^2, the x=t/ell,r=z/ell
equation with principal part partial_x^2-partial_r^2 has metric forcing
-2 ell^2 rho^2 QEred and scalar forcing +ell^2 rho^2 Qphi.

## Independent coordinate calculation

independent_coordinate_oracle.py imports no producer geometry or angular
implementation. It constructs genuine smooth local coordinate families with
the rational stereographic S3 metric, Y=n0^2-n1^2, a curved nonnormal warped
background, all seven populated metric amplitudes, and arbitrary populated
ordinary y,t two-jets. It expands the coordinate inverse, Christoffel symbols,
Ricci, matter, gauge and varying-connection reduction directly through
epsilon squared. A separate indexed path constructs covariant two-jets and
checks the independently derived formulas.

check_frame_against_coordinate.py then transforms every producer q,dq,ddq,
deta,ddeta,dQC slot from its redundant ambient frame into those coordinates.
It compares the full QRic,QE,QEred_general,actual-lowered QC,background-lowered
QC and Qphi_general against the direct coordinate expansion. Both manufactured
families are nonharmonic, so the general cross terms are exercised rather than
silently discarded. Final normal and optimized Python runs are byte-identical.

The finite families verify implementation; the universal formulas are proved
by the preceding derivation. Negative controls independently detect omission
of the derivative-slot connection, wrong time raising, omission of the general
2PC1 term, confusion of gauge lowering conventions, a missing Ricci product
factor, and a missing scalar-stress factor. Producer tests add the independent
conformal metric identity and check both exact integer and Fraction inputs.

## All-point angular projection

The projector uses rational monomial coefficients, not collocation values.
Invertible coefficient witness matrices prove basis ranks 3,2,7. Every
candidate projection is then reconstructed in every tensor component, and
every coefficient of the reduced polynomial remainder must vanish modulo
n0^2+n1^2+n2^2+n3^2-1. Thus acceptance proves equality everywhere on S3.

All 42 metric two-jet entries and all six scalar two-jet entries are independent
Laurent generators during the complete source export. This proves the maps
for arbitrary permitted jets rather than selected substitutions. The final
coefficient payload agrees exactly with the reviewer-regenerated all-point
export. The final integer-exactness change preserves that formal payload and
is covered by final producer normal/optimized export and coordinate checks.
The parent normalization and 23-slot packing are retained. Both TT orientations
are present. Deliberate off-span, deleted-TT, wrong-Hessian-factor and swapped-
packing controls are rejected. No claim that every possible angular
perturbation lies in the selected 23-dimensional sector is made.

## Actual source intervals

check_independent_source_intervals.py evaluates the final collected maps on
two full closed actual-first-jet cells with a separate exact Fraction endpoint
interval implementation. It checks 94 coefficient enclosures per cell and
verifies that the producer's outward-rounded intervals contain every result.
All 13,126 monomial entries of the full map are quadratic in the first-field
jet variables. The actual-jet oracle is an authenticated analytic premise;
the independent evaluator does not replace it with comparator values or
differentiate a value interval.

On x in [1/4,1/4+1/4096], r in [5/4,5/4+1/4096], the natural and reduced
radial Einstein-source means are strictly negative. Their displayed scale is
approximately [-6.120e17,-5.030e17]. On the reflected cell x in
[3/2,3/2+1/4096], r in [1/4,1/4+1/4096], the scalar-source mean is strictly
negative, approximately [-1.2742e11,-1.0381e11]. The reduced radial Einstein
mean on that reflected cell contains zero. Exact rational intervals are in
the receipts; these decimal summaries are not decision inputs.

The separate numeric-review directory owns the independent audit of actual
jet reconstruction, complete causal-domain coverage, global coefficient caps
and numerical controls. This geometric review's cell evaluation does not
itself establish a new global cap or compute the second response. Angular
coefficient bounds should not be mislabeled as raw frame-component bounds.

## Intrinsic shell scalar and moving-shell scope

The reviewed intrinsic scalar is I=gamma^mu nu d_mu psi d_nu psi with
psi=X^*phi and gamma=X^*g on the same material shell. Since d psi0=0,
I0=I1=0 and I''=2 gamma0^mu nu d_mu eta_material d_nu eta_material. Joint
field/embedding gauge covariance is essential. In the registered fixed-shell
representative f|shell=0, eta_material=vY; an additional independent physical
shell displacement is not assumed away.

Consequently the normalized sphere mean is

    mean(I'') = (-vt^2+8v^2/c^2)/(3rho_b^2).

For H=c vt/v, H_t=c N/v^2 with N=v vtt-vt^2+theta v vt. The independent
check_independent_shell_zeros.py reads the actual jet boxes but recomputes all
interval products, powers, inverses and sums directly with Fraction endpoints.
It proves v nonzero and N<0 throughout each of [5/4,11/8] and [13/8,7/4],
with H above sqrt(8) at each left endpoint and below -sqrt(8) at each right
endpoint. Threshold decisions use rational signs and squares, not a floating
square root. The factorization v^2(8-H^2)/(3rho_b^2 c^2) proves exactly two
simple mean zeros in each bracket, positive mean between them, and negative
mean outside them within that bracket. The brackets are disjoint, so four
distinct simple zeros are certified. Additional zeros elsewhere are not
excluded. Independent normal/optimized receipts agree byte for byte.

These are conclusions about an intrinsic perturbative coefficient and its
mean on the chosen background foliation. They are not a physical force, a
finite-amplitude sign theorem, a nonlinear family, or a second-order waveform.
The retained smooth background/master/gauge existence and differentiated
error certificates remain conditional analytic premises.

## Receipt policy

REVIEW_RECEIPT.json lists the final accepted inputs and authoritative replay
reports. Earlier unqualified NORMAL/OPTIMIZED reports are retained as
development receipts when a FINAL counterpart exists. Their old hashes are
not silently attributed to the final source. No frozen parent file or
producer implementation was edited by this reviewer.
