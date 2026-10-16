# Full common-gauge quadratic bulk forcing: bound proof

16 October 2026. Actual numerical instantiation and portable exact reports
accompany this proof. Earlier trees remain frozen.

The strict uniform component bounds are:

| Second parameter source | Absolute bound |
| --- | ---: |
| Natural Einstein QE, each of all 15 symmetric components | 9.444e20 |
| Registered reduced Einstein QEred, each component | 9.444e20 |
| Scalar Qphi | 3.855e15 |
| Quadratic relative-harmonic vector QC, each component | 5971 |
| Covariant derivative of QC, each ordered component | 2.684e9 |

The metric covariant 0/1/2 jets are bounded componentwise by
(0.01771,4497,8.792e8); the scalar covariant 0/1/2 jets by
(176500,2.173e10,4.355e15). These are exact finite-decimal strict caps,
verified against the computed rational upper bounds. Equal displayed
natural/reduced caps are rounding only; the exact reduced bound includes
the positive derivative-QC majorant.

The actual first pair is q=h_long+L_xi g0 and eta=chi_long+p f Y,
with the *same actual solved first gauge vector* xi in both transformations.
Use the affine pair path (g0+eps q,phi0+eps eta), with second parameter
derivatives (not coefficients of eps squared). The natural residual is
E=Ric-dphi tensor dphi-(2/3)U(phi)g; the registered reduced residual is
Ered=E-nabla_(A C_B), with symmetrization containing one half.

Domain: ell=10^-6, x=t/ell, r=z/ell, 0<=x<=3, r>=0, x+r<=8,
times the whole unit S3. The component norm is the maximum of absolute
individual components in the adapted frame ey=partial_y,
e0=rho^-1 partial_t, ea=(rho cosh t)^-1 Ea. It includes derivative slots
when applied to covariant jets. It is not a Lorentzian contraction or a
supremum over boosts. Angular frames may be any sphere-orthonormal frames.

## 1. Actual fields and derivative input

Use the six frame sectors from the pinned first-fields note. With c=cosh t,
H=tanh t, a=rho_y/rho and qtime=c^-3/2, their scalar amplitudes are

    yy: -4u-2 rho^-1 f_z, multiplying Y;
    y0: rho^-1 f_t+b_z, multiplying Y;
    ya: f/(rho c)-c d_z, multiplying D_aY;
    00: -2(u+af+b_t), multiplying Y;
    0a: c d_t-b/c, multiplying D_aY;
    ab: 2(u+af+Hb)Y delta_ab + 2d Hess_abY.

The scalar is eta=(v+pf)Y. Here

    u=C_rec qtime w,
    v=B qtime(w_z+k w),
    C_rec=alpha/rho^2,
    B=3alpha/(rho^3 p),
    k=rho(5a/2-U'/p), alpha(0)=1.

All background values/jets are on the whole 8ell collar, from the frozen
normalized static report and its first-exit proof. B, k and U derivatives
are recomputed with exact interval Taylor arithmetic from those global
jets; no initial-annulus jet report is used. Master C3 and gauge C3 mixed
caps bound the actual fields. Improved actual solver C0/C1 caps may be
intersected with the higher-norm caps. They are not differentiated.

On S3 the individual covariant derivative component caps of Y through
order four are 1,2,4,8,16. Thus the three angular tensor types Y, DY,
Hess Y have derivative cap triples (1,2,4),(2,4,8),(4,8,16).
Sphere covariant differentiation is performed on every angular slot,
including the angular derivative slots. The sphere connection is not
held zero under repeated differentiation.

## 2. Baseline connection and ordered derivatives

Define a metric baseline connection D that parallelizes the radial/time
split and retains the intrinsic sphere connection on the angular block.
It is torsionful. Set K=nabla-D. Its nonzero frame components are exactly

    K^0_(0y)=a, K^y_(00)=a,
    K^a_(by)=a delta_ab, K^a_(b0)=beta delta_ab,
    K^0_(ab)=beta delta_ab, K^y_(ab)=-a delta_ab,
    beta=H/rho.

For a separated frame component F(y,t) T(n), write s=(rho c)^-1.
Its ordered baseline second derivatives are

    yy: Fyy T;
    y0: rho^-1(Fyt-a Ft)T; 0y: rho^-1 Fty T;
    00: rho^-2 Ftt T;
    ya: s(Fy-a F) D_aT; ay: s Fy D_aT;
    0a: s/rho(Ft-H F)D_aT; a0: s/rho Ft D_aT;
    ab: s^2 F D_aD_bT.

These include derivatives of rho^-1 and s. They are not all obtained by
commuting D or by treating frame scales as constants. Proper-radial
amplitude derivatives are obtained from z derivatives by

    Fy=-rho^-1 Fz,
    Fyy=rho^-2 Fzz+(a/rho)Fz,
    Fyt=-rho^-1 Fzt.

The scalar and metric amplitude jets through order two, with these ordered
rules, give the baseline maxima H0,D1h,D2h and E0,D1eta,D2eta.

## 3. Full covariant jets, including the derivative slot

Let K0 bound max_(A,B) sum_C |K^C_AB|, and K1 bound
max_(A,B,C) sum_E |(D_A K)^E_BC|. Valid choices are

    K0 <= |a|+|beta|,
    K1 <= max(|ay|+|a beta|, sech(t)^2/rho^2).

Intrinsic angular DK is zero because K uses parallel splitting tensors.
Radial/time coefficient derivatives are ay, beta_y=-a beta and
e0(beta)=sech(t)^2/rho^2.

For the metric rank-two tensor,

    H1 = D1h+2K0 H0,
    H2 = D2h+5K0 D1h+(2K1+6K0^2)H0

bound max|nabla q| and max|nabla nabla q|. The second formula follows by
first applying D to nabla q=Dq-K on each original slot, then converting
the resulting rank-three covariant tensor by subtracting K on **all
three slots**. The extra derivative slot contributes K0 H1. Dropping it
would incorrectly replace 5 by 4 and 6 by 4.

For eta, E1=D1eta and E2=D2eta+K0 E1 bound its gradient and Hessian.
This connection correction is also retained.

## 4. Universal dimension-five contraction bound

Let n=5. Raising an index in the adapted orthonormal frame changes at most
a sign and introduces no additional sum. The first connection variation
P and its covariant derivative satisfy

    P0=(3/2)H1, P1=(3/2)H2.

For J=-2q*P, where one index is contracted,

    J1=2n(H1 P0+H0 P1).

The actual first relative-harmonic identity is C1=DC[q]=0. It holds for
the actual solution throughout D and its C1 derivative is zero as well.
The enclosure boxes need not satisfy this identity for arbitrary independent
selections: the simplification is made analytically for the actual pair
before the remaining terms are bounded.

The full natural Einstein Hessian has the common component bound

    QE <= 2n J1 + 4n^2 P0^2 + 2E1^2
          +(2/3)(|U''|E0^2+2|U'|E0 H0).

The first two factors count the two Ricci derivative contractions and the
two quadratic connection products with their Taylor factor two. The bound
covers every one of the 15 independent symmetric Einstein components,
without projection to scalar angular sectors.

Since C1=0, the quadratic gauge vector is QC=-2q^(BC)P^A_BC. Therefore

    QC <= 2n^2 H0 P0,
    nabla QC <= 2n^2(H1 P0+H0 P1),
    QEred <= QE+nabla QC.

The last coefficient is one: the half in symmetrization bounds the
average of two terms, each bounded by nabla QC. Lowering QC uses the
parallel background metric. No fourth gauge derivative is used.

The background scalar gradient has only its radial component p. Its
Hessian is diagonal in this frame, with entries py,-ap,ap,ap,ap,
where py=U'-4ap. Let Phi2=max(|py|,|ap|). Then the actual simplified
scalar Hessian source has the bound

    Qphi <= 2n^2(H0^2 Phi2+H0 E2+H0 P0 |p|)+|U'''|E0^2.

Only n^2 terms are needed for its first contraction because the background
Hessian is diagonal; the gradient contraction also uses only its radial
component. All scalar pullback terms, including +pf, have been retained.

## 5. Scope and verification contract

These are absolute continuous component bounds for the direct affine-path
second derivatives QE,Qphi,QC,nabla QC and the registered QEred. The
forcing in the second equations has the opposite sign; its absolute
bounds are unchanged. The coefficient of eps squared is one half the
displayed second derivative. No derivative of QE, Qphi or QEred is bounded
here; nabla QC is available directly from the metric two-jet.

The exact checker instantiates every quantity from authenticated reports,
checks the frame and component-count identities, and runs in normal and
optimized Python with meaningful corruption controls. Independent review
uses a separate source expansion and audit. This is no second-order
solution, full boundary/gauge/constraint forcing certificate, nonlinear
amplitude/remainder estimate, or observational interpretation.

## 6. Authenticated implementation and replay

`check_common_gauge_bulk.py` reads copied reports from the frozen public
first-gauge capsule. Their exact hashes and lengths are checked against
the externally pinned manifest
`adb6e3a7e00a8e477af80a49bc6b7a0f4037fafd0dd42ab8175c8b55e3e4465c`.
The unpublished first-field/source note and completed handoff are separately
pinned for convention provenance. Programs from earlier producers are not
imported or executed by this checker.

Static interval jets are normalized by factorial. Product/inverse Taylor
arithmetic reconstructs B,k,C_rec and the actual potential jets globally.
Nonnegative normalized mixed derivative majorants then compose the actual
metric/scalar amplitudes, with the full two-variable Leibniz rule. Gauge
jets of order zero and one use the smaller of the actual comparator-solver
enclosure and the higher-norm report; order two and three use the certified
mixed-derivative report. The comparator error itself is never differentiated.
All interval operations round outward, and positive derivative majorants
round upward, using exact Fractions on a grid of denominator 10^40.

Direct exact indexed tensors test the complete Einstein sums in all 15
independent components, scalar contractions, five QC components and all
25 derivative-QC components. A nontrivial harmonic first/second-jet family
has both C1 and its derivative equal to zero in these tests. The independent
reviewer additionally checks the source formulas against direct affine
inverse/connection expansions on generic curved nonnormal finite jets;
neither finite algebra test is presented as a PDE existence proof.

Ordinary and optimized reports are byte-identical. Seven semantic controls
run in each mode: radial scale derivatives, the metric derivative-slot
connection, scalar pullback, the actual fourth sphere derivative, both
quadratic-gauge contraction sums, the scalar Hessian connection, and the
second-parameter scalar-stress factor. Every control raises its named
exception. CLI outputs must be new and outside the source/input tree.

The report also records the universal source formula constants, explicit
background coefficient caps and the normalized S3 Hilbert-Frobenius
conversion: a 5x5 component tensor with component cap q has Frobenius
norm at most 5q. This is a positive component norm, not a Lorentzian
contraction. The Einstein source cap and its Frobenius conversion must
not be interchanged silently.
