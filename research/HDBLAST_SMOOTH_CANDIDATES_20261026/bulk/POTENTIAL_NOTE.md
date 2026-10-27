# A localized derivative bound for one scalar bulk-source term

This note is a conditional arithmetic consequence of an earlier actual-first-pair
2-jet enclosure on the small evolved-packet cell
`1/4<=x<=1025/4096`, `5/4<=r<=5121/4096`. It bounds one summand of the
scalar quadratic bulk source. Other summands contribute to the same harmonic
coefficient and may cancel it; no sign or bound for the complete scalar source,
its error, or the Einstein residual follows.

The direct affine-path Hessian of the natural scalar equation includes

    Qphi_potential = -V'''(phi0) eta^2,

where the registered common-gauge first scalar is `eta=s(t,y) Y` and the
background `phi0(y)` is static. Write its coefficient of the angular function
`Y²` as `P(t,y)=-V'''(phi0(y))s(t,y)^2`. The stored first-field jet symbols
`s` and `st` are ordinary first-scalar amplitude and its physical coordinate
time derivative. For `x=t/ell` at fixed r (so fixed y),

    P_x = -2 ell V'''(phi0) s st.                             (1)

The inherited cell report bounds all three factors in (1) by strictly positive
intervals. Exact rational endpoint multiplication yields

    -342000000 < P_x < -3520000.                              (2)

The retained forcing convention defines the scalar equation's potential
contribution as `Fpot=ell² rho² P Y²`. Since background rho is static, (1)
also gives `(Fpot)_x=ell² rho² P_x Y²`; using the same cell's positive rho
interval yields the coefficient bound

    -5.78 < coefficient_Y²((Fpot)_x) < -0.059.                 (3)

These are coefficients in the registered coordinate/parameter convention, not
observed forces or radiation. The `Y²` function itself vanishes at some angular
points. The other geometric and scalar terms in Qphi have not been bounded in
derivative and their sum may cancel this term. These inequalities therefore do
not determine whether the full forcing derivative is nonzero.

The checker authenticates `ACTUAL_FIRST_JETS_NORMAL.json` by SHA-256
`0aaa7b2fe30f546db92ef8b26a7224ffa7f6812d3780c331661395aa1991d3b8`, verifies
its selected cell and exact endpoints, and computes both intervals using Python
fractions. The starting report is itself conditional on the retained background,
master and first common-gauge existence/error certificates; this checker does
not rerun them. It does not differentiate a C0 error envelope or independently
construct a field solution. Run

    python ../checks/check_potential_derivative.py --repo-root /path/to/HDblast
    python -O ../checks/check_potential_derivative.py --repo-root /path/to/HDblast

The direct complete reduced Einstein source contains metric second derivatives.
Because `q=h+L_xi g0`, bounding an additional derivative of its generic direct
formula typically calls for fourth derivatives of xi. Likewise the scalar
source includes a second derivative of eta, which in turn contains a first
master derivative, so its derivative can need a fourth master derivative.
The inherited general certificates reach C3; the isolated computation above
uses only first scalar derivatives and does not close this gap. PDE-based
identities or special cancellations might reduce derivative requirements and
must be proved before using them; no claim that C4 is necessary for the exact
system is made.
