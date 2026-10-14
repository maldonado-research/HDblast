# Root mathematical audit: shell-contact continuation

13 October 2026 (America/Los_Angeles). Internal independent reasoning; not
external peer review. This note is an audit of the proof's scope, not a
replacement for the constructive existence and junction derivations.

## 1. Exact geometric domain

Put b=r_b>0, D=cos(tau)+cos(chi), and use
T=b sin(tau)/D, R=b sin(chi)/D. The interior is the hemisphere
0<=chi<=pi/2 and q=b/D. On |tau|<pi/2, D>0 and

    s=R^2-T^2=b^2(cos(tau)-cos(chi))/D.

For each fixed tau this is increasing in chi, with exact range

    -b^2 tan^2(tau/2) <= s <= b^2.

Thus a compact cylinder |tau|<=tau_* lies in the actual analytically
continued background whenever tau_*<pi/2 and
b^2 tan^2(tau_*/2)<s_-. This is a bound on s, not the incorrect requirement
that every ambient time on the cylinder be smaller than sqrt(s_-).
At its center T=b tan(tau/2), whereas at the shell T=b tan(tau).

For collar distance x=pi/2-chi and initial support in [a,d], finite
propagation gives x<=d+|tau|. The stronger bound d+2 tau_*<pi/2
keeps the entire causal support in chi>|tau|, hence s>0. This family
can interact with the shell while remaining in the original positive-s
background region. It is a different family from the preceding near-cone
nonzero-flux pair; no long-time continuation of that pair is inferred.

## 2. Conformal transformation and boundary inertia

In five spacetime dimensions the conformal scalar weight is 3/2. With
psi=q^(-3/2)u, conformal covariance and scalar curvature R(cylinder)=12
give the exact bulk equation

    u_tt-Delta_(S4)u+Q u=0,  Q=9/4+q^2 V_c(s).

At chi=pi/2, q_b=b sec(tau) and
n_eta psi=q_b^(-5/2)[u_chi-(3q_b/(2b))u]. The even boundary law
n_eta psi-c Box_(Sigma,eta)psi+d_e psi=0 becomes

    u_chi + (mu u_t)_t - mu Delta_(S3)u + beta u=0,
    mu=c/q_b>0,
    beta=q_b(d_e-3/(2b))
         -mu[(3/2)sec^2(tau)+(3/4)tan^2(tau)].

The odd Robin law becomes u_chi+alpha u=0, with
alpha=q_b(d_o-3/(2b)). The signs of alpha and beta are not prescribed.
All coefficients are smooth and time-even on the compact cylinder;
mu' is odd. Positivity of mu controls an auxiliary boundary variable,
not an independently added physical shell field or gravitational energy.

## 3. Energy estimate and a Robin pitfall

For the even equation, choose constants m,k>0 and let z be the trace.
Use the positive control energy

    E=1/2 int_M (u_t^2+|grad u|^2+m^2 u^2)
      +1/2 int_boundary [mu(z_t^2+|grad z|^2)+k z^2].

Exact integration by parts and the boundary equation yield

    E'=int_M (m^2-Q)u u_t
       +(mu'/2) int_boundary (|grad z|^2-z_t^2)
       +(k-beta) int_boundary z z_t.

A valid finite coefficient is

    C=sup|m^2-Q|/m + sup|mu'|/mu
      +sup|k-beta|/sqrt(k mu),

so |E'|<=C E. The compact-slab positive lower bound on mu is essential;
this estimate is not uniform as |tau| tends to pi/2.

For odd Robin data, adding an arbitrary positive boundary z^2 term and
bounding its derivative by bulk H1 x L2 energy is unjustified: the trace
of u_t is not controlled by that norm. Instead retain the actual Robin
term alpha z^2, choose the bulk mass sufficiently large to obtain
coercivity via the epsilon trace inequality, and control alpha' z^2.
A direct characteristic estimate for each finite angular mode is another
valid route. None of these control energies is asserted to be the physical
Einstein-scalar Hamiltonian.

## 4. Constructive scope avoids a regularity gap

A formal energy estimate plus Galerkin compactness alone does not prove
smooth solutions with a dynamic boundary. In particular projected initial
accelerations can violate uniform compatibility bounds. The core claim
should use the explicit one-dimensional characteristic construction for
arbitrary smooth compact collar profiles and fixed spherical harmonics,
or finite sums of these. This is already an infinite-dimensional family.
An all-angular-mode Sobolev completion requires additional uniform estimates
and is not inferred here.

For lambda=l(l+2), v=cos(x)^(3/2)u_l satisfies

    v_tt-v_xx+[(lambda+3/4)sec^2(x)+q^2 V_c]v=0.

Writing p=v_t+v_x, m=v_t-v_x and z=v(t,0), the even boundary is the
ordinary differential equation

    mu z''+(mu'+1)z'+(mu lambda+beta)z=p(t,0),
    m(t,0)=2z'-p(t,0).

For odd data z'=p-alpha z and m=p-2alpha z. These equations coupled
to characteristic integrals provide a causal local construction. Smooth
initial profiles vanishing near the boundary meet every corner compatibility
condition. The coefficients remain bounded on a strictly interior collar.

## 5. Arrival versus actual activation

Finite propagation establishes zero shell trace before the distance a;
it does not alone prove nonzero trace afterward. The additional argument
uses zero initial velocity, time-even coefficients and uniqueness to obtain
an even solution. If its trace vanished on (0,tau_0), the boundary equation
would force both v and v_x to vanish on x=0 for |tau|<tau_0. For a single
angular harmonic the equation may be treated as a 1+1 hyperbolic equation
with x as time. Local uniqueness then gives v=0 in |tau|+x<tau_0,
contradicting an initial profile whose support has infimum a<tau_0.
This proves a nonzero trace somewhere in every (a,a+epsilon) in the
controlled interval. A smooth compact profile is flat at its support edge;
neither a nonzero value at exactly a nor nonzero trace at every later instant
is claimed.

## 6. Local non-gauge test for l>=2

In a spherical annulus write the background as
g=g_ab(z)dz^a dz^b+w(z)^2 q_ij dtheta^i dtheta^j. The scalar-harmonic
part of a covector gauge generator is xi_a=f_a(z)H_l,
xi_i=J(z)D_i H_l. Its angular trace-free Lie derivative is
2J(D_iD_j H_l)^TF. On the unit S3,

    ||(D_iD_j H_l)^TF||^2=(2/3)lambda(lambda-3)||H_l||^2.

This is positive for l>=2. The reconstructed metric has only angular trace
and no mixed orbit-angular component, so a putative pure gauge would first
require J=0, then f_a=0 from
(Lie_xi g)_ai=[partial_a J-2(partial_a w/w)J+f_a]D_i H_l.
Orthogonal angular scalar/vector/tensor sectors cannot cancel this scalar
projection. Thus a nonzero finite l>=2 reconstructed family is not pure
smooth gauge on the annulus. This argument does not transplant the previous
small-ball Cauchy positivity result to a distant collar and does not classify
a completed gauge quotient or finite-cut charges.

## 7. Physical and publication scope

The full junction derivation must separately establish induced metric and
scalar matching, Israel and scalar jumps, and the normal-force identity.
An auxiliary scalar boundary theorem alone is insufficient. The inherited
bulk residual identity then supplies every bulk linearized field equation.
Actual covariant current balance must retain the odd moving-shell term;
the auxiliary even boundary energy must not be substituted for it.

This is conditional linearized dynamics of a fixed model and background.
It establishes neither nonlinear integrability, stability, reheating,
observational discrimination, a higher-dimensional origin of the Big Bang,
external mathematical originality, nor journal acceptance.

## 8. Time-odd activation and thin-annulus physical pairing

The activation argument also applies to zero initial displacement and nonzero
initial velocity: time reversal gives an odd solution; identically zero trace
on a forward interval gives zero Cauchy data for the sideways 1+1 equation
on the corresponding symmetric interval. Uniqueness then forces the initial
velocity to vanish up to that depth, contradicting its support infimum.
Support merely contained in [a,d] is insufficient to identify a as the exact
arrival infimum; a must be the infimum of the chosen nonzero profile.

Independently integrating the exact physical Cauchy bilinear by parts gives
K(q,q)=integral[m{(r q_r)^2+(2+9sL)q^2}+6r^7 q^2/Omega],
where m=9r^3/(4 Omega f^2), and r(log m)_r=3-2sL. For an annulus I
of width w, m*=min_I(m r^2)>0 and M*=max_I m(-2-9sL)_+, radial
Dirichlet Poincare gives K(q,q)>=(pi^2 m*/w^2-M*) integral q^2.
Shrinking an annulus within any fixed regular shell collar makes this
strictly positive, without assuming the small-ball sign bound from the
prior capsule. The diagonal positivity proves weak nondegeneracy by
testing a kernel pair (q,v) against (0,q) and (v,0). K need not be symmetric.

The two time-parity partners (q,0),(0,q) in the same reflection parity both
reach the shell and have doubled total physical pairing 2K(q,q)>0. Full
current balance with the odd moving-shell cut term conserves this total
pairing. This does not imply a nonzero local shell Wronskian, positive
physical energy, completeness, or nonlinear integrability. Junction note
Section 7 independently derives these identities, with six positive checks
and a wrong integration-by-parts control added to its existing roster.
