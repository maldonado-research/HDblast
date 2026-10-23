# Certified necessary enclosure of the retained regular branch

23 October 2026. Private bounded input-precision pilot. This proof sharpens the
same regular branch selected by the authenticated 11 October theorem. It does
not independently construct a finite Dirichlet endpoint map. Its existence and
uniqueness premise is the retained theorem, restricted to that theorem's tube.
No floating shooting residual is used. The old analytic detuning remainder is
used only to select an initial region, and is not the new error certificate.

## Result and finite scope

At delta=1/1000 and c=96714321425692/161839258930731,

    -0.00008405268308665527 <= beta <= -0.00008405268308663626.

Both endpoints are exact terminating rationals. The radius is 9.505e-18, below
the provisional 1e-11 target. This alone does not prove the initial/free gate:
that requires the separate physical-profile difference computation.

The cap was one retained root region, graph degree four, one comparison domain,
four rational Newton seed proposals, and one validated necessary-enclosure
contraction. No radial IVP or source panels were run, and no cap was expanded.
All evidence-producing arithmetic is exact rational arithmetic. The final rho
projection uses integer-square-root bounds checked by exact squaring. Informative
floating displays are not inputs to any acceptance inequality.

The complete old capsule (43 payloads plus MANIFEST.json) is checked against an
external manifest pin, an independent immutable-public-source readback, and all
44 exact members of the public ZIP. The independently retained public ZIP
SHA256 is 69ab054a4a533756f561691c28906527dbf7d4bb90eaff0ce4198a5cbfd661a8.
Authentication is repeated after arithmetic. See BRANCH_DEPENDENCY_AUTHENTICATION.json.

## Exact state and entry into the comparison domain

Use x=eta, h=rho^(-2) (so h itself is the physical quantity called H²),
z=p/eta, and a=rho_y/rho. The model is

    W(x)=1/3+x²+x³/3, f(x)=1+c+c x,
    U(x)=(W'(x))²/2-2W(x)²/3,
    U=-2/27+(14/9)x²+(50/27)x³-x⁴/6-(4/9)x⁵-(2/27)x⁶.

Let M(x)=U'(x)/x, with its polynomial continuous extension at zero. Away from
the cone the exact equations and constraint give

    x'=x z, h'=-2a h,
    z'=M(x)-4a z-z²,
    a²=h+x² z²/12-U(x)/6.

Primes here mean differentiation with respect to proper radial distance y.

The retained theorem puts beta in its displayed old analytic interval, which
is strictly negative and has |beta|<9/100000. The regular scalar cannot change
sign; its divergence equation and strict convexity imply x<0 and p<0 and
|x(y)|<=|beta|. The center scalar is nonzero: otherwise uniqueness of the regular
interior would give x=p=0 and beta=0. Thus z is defined on the entire segment.
The regular cone expansion gives z->0 as y->0.

Exact polynomial bounds on |x|<=epsilon=9/100000 verify 0<M<4, U>=U(0), and
U<0. The positive a branch starts at the cone, and its constraint prevents it
from reaching zero. The vector field for z points into [0,2] at z=0 and z=2,
so 0<=z<=2 everywhere along this segment. The retained theorem also directly
supplies positive a; the constraint argument makes the use explicit.

Since h'=-2ah<0 and h diverges at the regular cone, the trajectory crosses
h0=1/100 exactly once before its endpoint: the exact shell law and old beta
interval give 0<h_endpoint<h0. On the entire remaining segment, the constraint
and the exact polynomial bounds prove

    a_min=1/9-1/1000000 <= a <= a_max=3/20.

No center amplitude, center-start radius, numerical cone jets, or numerical
radial trajectory is needed for this comparison.

## Polynomial graph and exact defect

GRAPH_POLYNOMIALS.json contains the exact total-degree-four rational
polynomials Z(x,h) and A(x,h). They start as

    Z=14/9-(42/5)h+(25/23)x+...,
    A=1/9+(9/2)h-(7/27)x²+... .

These are approximating polynomials, not an assumed invariant manifold. Define

    B=2h Z_h-4Z,
    R=M-Z²-x Z_x Z+A B,
    D=h+x² Z²/12-U/6-A².

The producer independently reconstructs the exact model and each full
polynomial identity. R and D vanish through total degree four; their remaining
coefficients are retained exactly, not discarded. On |x|<=epsilon, 0<=h<=h0,
coefficientwise rational bounds prove

    0<Z<2, A>0, B<0,
    4a_min+inf Z-sup |x Z_x| > 1.88625 >= kappa=37/20.

The report contains all exact fractions. In particular A_min>0.05143; therefore
a+A>=a_min+A_min>0.

Let e=z-Z(x,h), with Z evaluated along the actual trajectory. Its exact error
identity is

    e'=R+B D/(a+A)-K e,
    K=4a+z+Z+x Z_x-B*x²*(z+Z)/(12*(a+A)).

To verify it, subtract the derivative of Z along x'=xz, h'=-2ah and use

    a-A=[D+x²*(z+Z)*e/12]/(a+A).

Since B<0 and z+Z>0, the final term in K is nonnegative. Thus K>=kappa on the
entire domain. This controls broad entry data; it does not replace z by its
endpoint approximation inside a damping coefficient.

## Explicit error barrier with cone memory

Write each coefficient of R and D as r_ij and d_ij. The report computes
Q=sup|B|/(a_min+A_min) exactly. Replace every |x|^i by epsilon^i and reduce
h^j for j>5 to h0^(j-5) h^5. The resulting nonnegative rational coefficients
C_j, j=0,...,5, satisfy

    |R+B D/(a+A)| <= sum_j C_j h^j.

Define the explicit comparison barrier

    E(h)=2(h/h0)^6 + sum_(j=0)^5 C_j h^j/(kappa-2j a_max).

Every denominator is positive. The smallest is kappa-10a_max=7/20. The memory
term also satisfies kappa-12a_max=1/20>0. Along the actual flow,

    (h^j)'+kappa h^j=(kappa-2j a)h^j,

so E'+kappa E dominates the residual majorant. At the unique entry h=h0,
0<=z<=2 and 0<Z<2 give |e|<=2<=E(h0). The scalar differential comparison for
Dini derivatives of |e| now proves |e(y)|<=E(h(y)) throughout the remaining
segment. This is the error certificate; there is no assumption that a formal
series converges or matches regular cone data. The first barrier term explicitly
retains and damps every possible regular-cone entry value in the stated range.

The old beta interval and exact shell polynomial give a rational upper bound
h_endpoint<=h_upper. Since E has nonnegative coefficients,

    |z_endpoint-Z(beta,h_endpoint)| <= E(h_upper) < 3.752255e-13.

## One validated root-enclosure contraction

At the shell, all quantities share the same beta parameter:

    h_s(beta)=delta*(W*f/9-W'*c/12)+delta²*(f²/36-c²/48),
    p_s(beta)=-W'(beta)-delta*c/2,
    a_s(beta)=W(beta)/3+delta*f(beta)/6,
    rho_s(beta)=h_s(beta)^(-1/2).

The metric constraint plus both shell laws imply the first identity exactly.
It is used only on the actual shell. Define the polynomial

    G(beta)=beta*Z(beta,h_s(beta))+W'(beta)+delta*c/2.

The exact scalar junction implies |G(beta_actual)|<=epsilon E(h_upper).
The producer evaluates the full polynomial chain derivative

    G'=Z+beta*(Z_x+Z_h*h_s')+W''

on the entire old beta interval and proves G'>3.55469. This derivative is that
of the comparison polynomial, not an uncomputed finite-Dirichlet Jacobian.

Four exact rational Newton seed proposals, rounded to 24 decimal places, give
one rational center b. They are only seed proposals and carry no existence
claim. The mean-value theorem and the already established branch existence give

    |beta_actual-b| <= (|G(b)|+epsilon E(h_upper))/inf G'.

All terms in this bound are exact rationals. The resulting interval is checked
to lie inside the old selection interval and widened outward to 20 decimal
places, giving the displayed certified endpoints. This is a necessary enclosure
of the retained existent root; it does not prove a new independent root exists.
Its scientific use is the sharper authenticated parameter for the same branch.

The report retains the exact one-parameter shell definitions. Its four numeric
outer projections are descriptive and must not be substituted by independently
selected midpoints. The local profile/jet evaluator must carry the common beta
parameter and separately prove its actual-minus-candidate error. No initial
precision result by itself supplies B1, a total U, an observable sign, or a
nonlinear family conclusion.
