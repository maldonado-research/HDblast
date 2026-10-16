# Independent actual-source numerical review

16 October 2026. Accepted within the conditional scope stated below. This is an
internal mathematical and computational review, not external peer review.
The review writes only this directory and temporary outputs; producer and
frozen predecessor files were not edited. Final identities are recorded in
`OUTPUT_PINS.json` and `REVIEW_RECEIPT.json`.

The accepted result is a continuous enclosure of the actual common-gauge
first metric/scalar ordinary two-jets and their geometric quadratic sources.
It includes every cell of the stated causal domain, the evolved radial
profile, and the physical shell source rows. It does not compute a second
solution or a continuous residual for such a solution.

## Actual derivative inputs and factorial conventions

The retained actual smooth background, master and mixed-boundary gauge-wave
constructions remain conditional premises. The review read their pinned
analytic derivations and the independently reviewed differentiated error
certificate. The master error includes the true initial velocity, the dynamic
shell lift and the reflected Duhamel terms through order three. The gauge
error includes the complete residual, coefficient commutators, initial
acceleration and stabilized third initial jet, and radial residual recovery.
No derivative is inferred from a value-only error bound.

The local predecessor bytes and five numerical input hashes are authenticated.
The original workspace copies were also compared byte for byte with all
listed frozen origins. The replay programs authenticate only the local copies
and never dereference those historical absolute origin strings. The inherited
independent error review binds both its producer and report, and the master,
coefficient and shell pins agree across the retained reports.

For a common smooth field F, the producer stores

    T[i,j] = (partial_x^i partial_r^j F)/(i! j!), i+j<=3.

Thus multiplication is ordinary finite convolution and inversion follows from
T(F) T(1/F)=1. Certified actual-minus-comparison derivatives are divided by
the same i!j! before being added. This encloses the derivatives of the common
actual fields; treating the entries independently only widens the image.
The algebra is a derivative jet at each point of a cell, not a Taylor
approximation in the cell width, so no spatial Taylor remainder is omitted.

Reconstruction differentiates W or a gauge amplitude at most once before
forming metric/scalar sectors. Ordinary output order two therefore consumes
only input order three. The artificially zero top coefficients returned by
the fixed-order `dx`/`dr` operations are outside the used output order; no
fourth master or gauge derivative is supplied or silently assumed. Actual
background ODE jets are available through order four.

`check_independent_actual_jets.py` uses ordinary binomial-weighted derivatives,
rather than normalized coefficients. It independently forms the coordinate
Lie-derivative formulas, divides by the moving background frame factors, and
checks all 48 metric-sector/scalar ordinary two-jet entries, including mixed
yt derivatives, against the producer. Ten further mixed rational quotient
derivatives check factorial normalization through total order three. The
same program independently recovers 150 actual-background ODE derivatives
and reconstruction coefficients through order four at three manufactured
base states. These finite implementation checks supplement the analytic
jet composition argument above.

## Comparison functions, packet tails and elementary intervals

With Qo(s)=sign(s)Q(|s|) and Pe(s)=s Q(|s|), the gauge comparison is

    P=(6ell/130)/4 [Pe(r+x)+Pe(r-x)-r(Qo(r+x)+Qo(r-x))].

For r+x>=0, its first combined term is x Qo(r+x). Writing s=r-x,
the second combined term is -x Qo(s) on s>=0 and (x-2r)Qo(s) on
s<=0. Differentiating these linear prefactors produces all their x and r
product terms. The apparent branch junction lies in an open zero collar,
so taking the hull of both branch images is valid on crossing cells. The
independent check substitutes signed polynomial packets into the original
uncollected Pe-rQo formula and checks all mixed derivatives on both branches.

On the packet, let S=1/[(r-1)(2-r)] and v=3-2r. The differentiated packet
polynomials through order three are those of exp(4-S). Their retained global
caps and the inequalities |Q^(n)|<=S^(2n)exp(4-S) justify the endpoint tails.
For the tail threshold S>=4096/63>2n, this expression decreases with S.
Outside the packet and at its flat endpoints all derivatives are zero.
Intersecting a local derivative image with the authenticated global cap is
valid because both enclose the same actual derivative. Independent ordinary
exponential/Bell derivatives and a separate direct rational exponential
series support 80 signed-packet and uncollected-comparator checks.

`interval_core.I` encloses every operation by rational outward rounding on
G=2^112. The reduced exponential argument lies in [0,1], so the alternating
degree-28/29 Taylor bounds enclose exp(-q); repeated positive squaring
preserves inclusion. The cosh/tanh values use these exponential intervals.
For the square root, the integer lower and upper brackets explicitly square
below and above the cosh endpoints. Derivatives of cosh^(-3/2) follow its
ODE v_t=-(3/2)tanh(t)v. They are derivatives of this smooth function, not of
selected interval endpoints. Sixty independently checked differentiated
temporal identities supplement these all-point enclosure arguments.

## Physical derivatives and shell identities

The coordinate relation is dy=-rho dz and z=ell r. For L=1/(ell rho),
rho_y=a rho implies

    F_y=-L F_r,  F_yy=L^2 F_rr+a L F_r,
    F_yt=-L F_xr/ell,  F_tt=F_xx/ell^2.

These are ordinary derivatives of background-frame amplitudes. They include
derivatives of every rho, cosh and reconstruction factor in those amplitudes.
The separate reviewed geometric engine supplies the covariant frame
connection and the prior derivative-index connection in its second jet.
Holding rho fixed would omit the a L F_r term and is rejected by the
independent manufactured test.

At the physical shell, f=0 and b_r=d_r=0. Available time derivatives of these
identities give f_x=f_xx=f_xxx=0 and b_xr=b_xxr=d_xr=d_xxr=0.
The producer imposes exactly these tangential identities at r=0. It retains
the unconstrained f_r, b and d amplitudes. Hence B=C and their first two
time derivatives vanish, while the tangential transformed metric generally
changes. The scalar trace is vY. No extra Dirichlet condition is imposed on
b or d, and no normal differentiation of a shell-only identity is used.

The physical shell helper agrees with the affine unit-normal Hessian before
elimination of normal derivatives. Its elimination uses the actual first
physical rows only on the shell. It retains all five tangential/scalar
families and uses sigma'=-2p, sigma''=4(1+eta_background), sigma'''=4.
One-copy outward, two-copy outward/inward and coordinate-forcing adapters
are exact, including the different metric and scalar factors. The report
audit verifies these adapters in every shell time cell and the shell witness.

## Exterior improvement and complete causal cover

Removing the dynamic lift is valid only on a whole box with r.lo>=x.hi.
For any point with r>=x, every backward characteristic point (s,q) satisfies

    q>=r-x+s>=s.

Thus its backward cone stays in the same exterior restriction and never
meets the physical shell at positive time. The lift 1_(x>=r)y(x-r) and its
jets through three vanish there, including x=r by the compatible flat corner.
Every other inherited global coefficient/kernel/positive resolvent majorant
remains an upper bound on the restriction. `check_independent_exterior_errors.py`
uses an independent adjugate inverse and rederives all exterior errors,
source defects and residuals, then restores the lift and recovers the inherited
baseline exactly: 92 exact rational comparisons. The actual-jet checker
also instruments the exterior selector and rejects crossing boxes as exterior.

For integer d>=1 the closed boxes

    [i/d,(i+1)/d] x [k/d,(k+1)/d],
    0<=i<3d, 0<=k<8d, i+k<8d,

are intersected with D={0<=x<=3,r>=0,x+r<=8}. Their union is exactly D.
Interior points are covered by their floor-index cell. A diagonal grid point
has an adjacent included cell immediately below its r coordinate; the
endpoints (0,8) and (3,5) are included by the last r or x cells respectively.
Clipping gives no points outside D. This establishes coverage of the closed
boundary as well as the interior.

On a diagonal cell, `clip_causal=True` asserts an actual-field enclosure only
on D intersected with the box. The smooth comparisons and whole-collar
background are evaluated on the larger rectangle; that only widens their
images at the actual points. The actual C3 error bounds are used only at
those points of D. No PDE solution or coefficient parity extension outside
D is assumed. The 81 default d=2 cells include all such diagonal pieces,
the shell and reflected joins.

## Integer source arithmetic, exact maps and the complete report

The source-only `fast_dyadic.D` ring has integer endpoints l/G,h/G. Addition
and negation are exact in ticks. Multiplication gives tick endpoints
floor(min four integer products/G) and ceil(max/G). On either same-sign
denominator the reciprocal endpoints are floor(G^2/h),ceil(G^2/l).
The identity ceil(n/d)=-floor(-n/d) is valid for either sign of d, and Python
integer division is mathematical floor division. Recursive powers have the
same operation and rounding order as the Fraction reference. Arbitrary-size
integers do not overflow. These formulas prove equivalence for all inputs,
independently of finite tests.

D imports external intervals through their exact Fraction lo/hi endpoints.
The source adapter converts every interval-valued metric/scalar/background
input to D before evaluating the maps; remaining constants are exact integers
or Fractions. Scalar arithmetic on either side follows the reflected methods.
No mixed operation with a reference-I left operand and D right operand is
required. `check_independent_dyadic.py` supplements this proof with 20,859
exact signed/grid-edge, coercion, power, inverse and denominator checks in
normal and optimized Python. The producer also compares all 74 witness source
endpoints with its Fraction reference.

Universal geometric formulas, the full 23-slot projection and the final exact
map payload have a separate independent review in `review/REVIEW_RECEIPT.json`.
Its independent exact endpoint evaluator checks 94 coefficients, including
QRic, on both exterior and reflected cells and is bound by
`review/SOURCE_INTERVALS_FINAL_NORMAL.json`. The source report and this review
bind the final map hash and current engine hashes. By the all-input ring
equivalence just proved, that independent interval-extension review applies
equally to the source-only integer implementation.

QEred and Qphi are simplified using the actual harmonic identity C1=0 on a
neighborhood. This identity is established by the retained actual gauge PDE.
It is applied analytically before rectangular jet enclosures are substituted;
arbitrary other points of an interval box need not satisfy C1=0 for its
polynomial image still to enclose the actual solution's source.

`check_independent_source_report.py` checks all 81 whole causal cells, all
64 evolved radial cells and all six whole shell time cells. It verifies every
registered source coefficient against the map schema, all interval orders,
every global/profile cap maximum and strict-sign count, the exact partitions,
the source/shell witness signs, and every shell orientation adapter. Three
complete selected causal-cell records are replayed numerically, including
the clipped diagonal cell. This is a complete receipt audit with selected
numeric replay, not an independent reimplementation of every cell's source
arithmetic; the interval-extension and analytic composition proofs provide
the all-point justification.

The forcing H conversion is independently recomputed from the full reported
coefficient caps. The review contracts the ambient tangent basis tensors by
an unreduced monomial integral, using the exact Dirichlet half-Gamma formula
for probability-normalized S3 moments, without the producer's inner-product
or polynomial-reduction routines. The Gram matrices agree exactly and are
positive diagonal. Mixed metric blocks have multiplicity two. Taking the
positive metric-plus-scalar Frobenius norm gives

    F_H <= 17387513156573 < 171610068644000.

This improves only the retained bulk C0 forcing bound. No shell Green-map,
operator, source quadrature, response residual or response majorant is
updated by this comparison. The coefficient caps are in the registered
unnormalized basis; a mean-sector sign is not a sign at every angular point.

## Executions, controls and limits

The actual-jet, exterior-error, integer-ring and complete-report review
certificates pass under both normal and optimized Python with identical
complete output bytes. The producer's retained actual-jet controls were also
replayed in both modes and match. The public independent control harness runs
12 nested rejecting processes: factorial, scalar pullback, radial frame factor,
radial conversion, changed pinned input and coherently changed error report
with stale independent review, each under both interpreter modes. Every child
must fail for its specific reason and leave no output certificate. The harness
is one top-level positive job; its 12 children are nested rejections, not 12
additional top-level jobs. A further independent reviewer execution of the
producer's five source mutations in both modes also rejected; those ten
duplicate source controls are excluded from this review's public roster.

No unresolved numerical defect remains in the final pinned artifacts.
Finite manufactured tests do not replace the displayed analytic arguments,
the retained smooth actual-solution premises, or the separate universal
geometric/projection proof. This review establishes neither a numerical
second response nor its residual, a nonlinear family, an invariant force,
or the shell zero theorem, which has its own independent review.
