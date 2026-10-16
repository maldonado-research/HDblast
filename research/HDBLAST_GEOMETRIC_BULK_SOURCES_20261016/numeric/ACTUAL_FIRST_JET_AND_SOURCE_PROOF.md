# Actual evolved common-gauge jets and quadratic source values

This unpublished checkpoint evaluates the actual first common-gauge solution
through its ordinary metric/scalar two-jets, then encloses its geometric bulk
and physical shell quadratic source coefficients on complete cells.  It uses
the retained smooth background, master, and mixed-boundary gauge solutions as
authenticated conditional premises.  It does not perform a new first-field
Picard solve or construct the numerical second response.

## Authenticated input and the meaning of a jet interval

`PREDECESSOR_PINS.json` authenticates the copied inputs, the inherited producer
and independently reviewed error certificate, and the explicit retained master,
gauge, actual first-field and physical-shell proof premises. `INPUT_PINS.json`
is the original five-input manifest, resolved only against this directory.
No earlier directory is imported or written. The copied `interval_core.py`
contains the outward rational interval, normalized radial Taylor, actual
background ODE and explicit bump routines used by the preceding checkpoint.
Only these relevant routines are used here.

Put ell=10^-6, x=t/ell and r=z/ell. For the smooth actual first solution on

    D={0<=x<=3, r>=0, x+r<=8},

the inherited certificate supplies each actual derivative of W-Wbar and
U-(P,0,0), through total order three. Wbar and P are single explicit smooth
comparison functions, defined in the pinned proof. Actual first metric
two-jets use the third raw gauge jet; actual scalar two-jets use the third
master jet. Every needed mixed derivative is present. No derivative is taken
of a value error or of an independent interval selection.

`actual_first_jets.py` forms bivariate normalized Taylor jets through total
order three for the actual coefficients, enclosing each certified derivative
separately. Products and inverses use the finite Taylor algebra. Evaluating
these expressions on intervals loses dependencies and can widen the result,
but it still encloses the derivatives of the actual common functions. The
comparison functions remain explicitly distinguished in the API and reports.

The actual background values use the authenticated shell box and whole-collar
ODE image. Its derivatives are calculated from the actual ODE through order
four. Alpha is enclosed by its integrated negative logarithmic derivative.
These are actual background jets, including the inherited shell uncertainty;
they are not derivatives of independently selected background midpoint values.

All nonlinear decisions use exact fractions. Every basic interval operation
rounds outward to the grid with denominator 2^112. Exponentials use a rational
alternating Taylor enclosure after halving and then repeated squaring. The
single rational square-root bracket needed for cosh^(-3/2) is certified by
integer square-root inequalities. Its higher derivatives come from

    v_t=-(3/2)tanh(t) v.

For the large collected source maps, `fast_dyadic.py` implements the same
interval operations using integer endpoint ticks. Addition is exact on that
grid; products use floor/ceiling after division by the grid; reciprocals use
floor(GRID^2/upper_tick), ceiling(GRID^2/lower_tick), including negative
denominators. Its rounding order is identical to the reference Fraction
implementation. The source run checks 284 rational/signed arithmetic cases
and complete equality of every exterior witness source endpoint against the
reference implementation. This changes evaluation cost, not the enclosures.

The explicit packet is evaluated through its actual derivative polynomials,
with analytic flat endpoint tails. Each derivative enclosure is intersected
with its independently retained global derivative supremum. This intersection
uses derivative certificates, not a differentiated displacement cap.

## Cancellation, reconstruction and ordinary physical derivatives

Let Qo be the smooth odd packet and c0=6ell/130. The inherited gauge comparison
is

    P=c0/4[Pe(r+x)+Pe(r-x)-r(Qo(r+x)+Qo(r-x))].

The correlated Pe-rQo products are simplified before interval evaluation.
Since r+x>=0, its first characteristic term is x Qo(r+x). The other term is
-x Qo(r-x) if r>=x and (x-2r)Qo(r-x) if r<=x. There is an identically zero
collar around the branch junction. All differentiated linear prefactors are
retained. Twenty independent polynomial mixed-derivative controls check both
characteristic branches, including the reflected one.

The actual reconstruction is exactly the retained first-fields construction:

    u=R cosh(t)^(-3/2) W,
    v=B cosh(t)^(-3/2)(W_r/ell+k W),
    R=alpha/rho^2, B=3alpha/(rho^3 p).

`frame_amplitudes` retains every registered gauge term, the raw vector's rho
factors, and the scalar pullback sign +pf. Its seven outputs A,B,C,D,E,L,M
multiply respectively Y,Y,DY,Y,DY,Y*kappa,Hess(Y) in the yy,y0,ya,00,0a,ab
sectors. The scalar amplitude multiplies Y. These are background orthonormal
frame components. The returned suffixes y,t,yy,yt,tt mean ordinary physical
derivatives of these frame amplitudes, before adding the frame connection.
The separate geometric adapter adds that connection, including the prior
derivative index in the covariant second jet.

With Lr=(ell rho)^(-1), the ordinary conversion is

    F_y=-Lr F_r,
    F_yy=Lr^2 F_rr+a Lr F_r,
    F_yt=-Lr F_xr/ell, F_tt=F_xx/ell^2.

The a Lr F_r term follows from rho_y=a rho and is checked against a second
application of the variable-coefficient differential operator. In particular
rho is never held constant during radial conversion.

At r=0, the actual physical rows f=0, b_r=d_r=0 and all their available
tangential derivatives are imposed as exact identities. They yield B=C=0
and their tangential derivatives, and the shell scalar is exactly vY.

## Sharper actual errors before shell contact

`exterior_error_caps.py` treats the backward-closed restriction

    Dext={0<=x<=3, x<=r, x+r<=8}.

The exact retained dynamic lift is 1_(x>=r)y(x-r). The pinned master join
proof supplies vanishing corner jets through order three, so this lift and
its derivatives through three vanish on Dext, including x=r. Every backward
point (s,q) in the unreflected retarded cone obeys q>=r-x+s>=s, so the cone
stays in Dext and meets no positive-time shell. All inherited global kernel,
coefficient, residual product and positive matrix-resolvent bounds remain
valid after this restriction. Only the lift row is removed.

The gauge error computation retains the complete residual, all coefficient
commutators, the nonzero stabilized acceleration, its radial derivative,
the initial third jet A(0)R(0)+R_x(0), and residual_r in third radial recovery.
With the lift restored, the new implementation must reproduce exactly all
four inherited master and thirty inherited gauge derivative error caps.
The exterior master errors are approximately

    3.45909e-6, 1.30606e-6, 2.09983e-5, 3.35972e-4.

These values are derivative approximation errors for the actual solution.
They are not claims that the comparison satisfies the full master equation.

## Exact geometric source evaluation and complete cells

The sibling projection certificate supplies algebraically collected sparse
Laurent maps in independent ordinary first-pair jet symbols. Its all-point
angular reconstruction occurs before any interval substitution. Numerical
evaluation of these maps preserves those exact cancellations; no attempted
equality of dependent interval angular expressions is used as evidence.
QEred and Qphi use the actual first harmonic identity C1=0, including its
derivative, established by the retained first-gauge PDE. The reported QC is
the quadratic gauge Hessian lowered with the actual metric convention in
that exact geometric adapter.

`evaluate_actual_sources.py` retains QE, QEred, Qphi, QC and all registered
bulk PDE forcing rows. Source normalization is the second parameter
derivative along the affine actual-first-pair path. Thus its epsilon^2
coefficient is one half the displayed Hessian. The scaled reduced metric
forcing is -2 ell^2 rho^2 QEred and the scalar forcing ell^2 rho^2 Qphi.

The full causal cover is finite. For integer d>=1 take the boxes

    [i/d,(i+1)/d] x [k/d,(k+1)/d],
    0<=i<3d, 0<=k<8d, i+k<8d,

each intersected with D. Their union equals D. Diagonal boundary cells are
retained explicitly. Comparison functions are harmlessly enclosed on each
bounding box; actual solution error bounds are asserted only on its
intersection with D. No solution outside D is assumed. At default d=2
there are 81 such cells, including the shell and reflected joins. The
additional 64-cell radial profile covers all r in [1,2] at the whole time
cell x in [1/4,1/4+1/4096].

The global positive bundle forcing bound uses the exact probability sphere
Gram matrix in each scalar/vector/tensor block. A mixed metric block is
counted twice. For coefficient caps c_i and Gram entries G_ij, the squared
block norm is bounded by sum |G_ij| c_i c_j, with its block multiplicity.
Summing the metric and scalar blocks and rounding its square root upward
gives a C0 bound in the same positive Frobenius bundle norm as the retained
causal inverse. This is only a bulk forcing bound. It does not update the
operator, shell Green map, retarded quadrature or response error.

A narrow whole evolved cell has strictly negative yy and 00 Einstein
Hessian mean-sector coefficients and a strictly positive y0 coefficient.
The complete exact rational intervals, all other components and the scalar
source enclosure are retained. A sign claimed for a mode coefficient is a
claim in the registered unnormalized basis; it is not a claim that every
angular point has that sign. The scalar-source interval can cross zero.

## Physical shell source rows

`shell_source_values.py` uses the pinned actual on-shell physical first rows
only at the shell. For H=A Y, m00=D Y, m0a=E DY,
mab=L Y*kappa+M HessY and eta=sY, set sigma'=-2p,
sigma''=4(1+eta_background), sigma'''=4. Its one-copy outward affine
Hessians are exactly

    QK=A0*gamma-L0*m,
    A0=a H^2/4-sigma' H eta/6-sigma'' eta^2/6,
    L0=a H+sigma' eta/3,
    Qp=p H^2/4+sigma'' H eta/2+2 eta^2.

The implementation retains all tangential tensor families and projects the
fixed Y^2, Y DY and Y HessY products by exact rational coefficients. A
separate rational manufactured-jet check recovers these formulas from the
pre-elimination unit-normal affine Hessian. It detects omission of the
perturbed-normal normalization. All one/two-copy orientation and coordinate
forcing adapters are reported explicitly; no shell-only identity is
differentiated normally. The default six-cell shell cover includes all
0<=x<=3. A narrow evolved cell has strictly positive one-copy scalar shell
Hessian mean coefficient.

## What has not been certified

These calculations evaluate actual first-pair derivatives and actual
quadratic source enclosures. They do not evaluate retarded source quadrature,
solve the coupled numerical second evolution, reconstruct a second waveform,
or certify a full continuous second-response residual with all input,
quadrature, reconstruction and boundary errors. No zero candidate is presented
as a waveform, and no invariant force or finite-amplitude conclusion follows.

Normal and optimized Python reports must match as complete bytes. Direct
checkpoint formula/input controls and separate inherited nested controls are
enumerated in `NUMERIC_CONTROLS.json`. All controls require failure, the
specific intended stderr reason, and absence of an output certificate.
