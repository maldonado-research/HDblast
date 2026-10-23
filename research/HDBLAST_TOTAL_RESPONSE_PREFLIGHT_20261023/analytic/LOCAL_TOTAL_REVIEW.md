# Conditional total-mean review and a zero-wedge estimate

23 October 2026. Independent internal mathematical review. This note preserves
all retained capsules. It supplies an analytic refinement and a small exact
arithmetic check, not an evaluated continuous response, a new PDE campaign,
external peer review, or a novelty claim.

## Result and conditions

The retained total problem is consistent with the selected physical initial
and harmonic one-jet. Its five invariant means close exactly, but every
internal metric/scalar, connection and Robin coupling must remain. The exact
total coefficient is U=(I-T)^(-1)(U_free+S), S=KF+LG. In particular
U_hom=(I-T)^(-1)U_free, U=U_hom+W and W=(I-T)^(-1)S. Neither summand
separately supplies the physical constraint completion of U.

On the specified backward cone and compatible zero-wedge subspace, the
retained coefficient norm bounds give the improved contraction

    q_w = 3409/8000000 = 0.000426125.

This follows from the elementary zero-wedge estimate below. It is a separate
valid bound, not a correction to the accepted q_D=3609/2000000. The inherited
cap ||U||_X<1.03*10^33 and all retained physical hypotheses remain premises.
The exact two-term tail is consequently <1.870300*10^26 in X. This does not
bound input, coefficient, interpolation or quadrature errors.

At the provisional 10^29 X accuracy, one can omit evaluation of S in the
first reconstructed term and charge its entire norm to the defect ledger:

    ||S||_X <= (3/2)*17387513156573 + 3851400
              = 52162547172519/2 < 2.609*10^13.

This is a rigorous way to enclose exact total U at that accuracy; it cannot
resolve the much smaller auxiliary W. It reduces the necessary pilot to an
initial/free-wave representation and one full coupled feedback term, if the
initial-data accuracy gate succeeds.

## 1. Geometry and the zero-wedge lemma

Write

    D+ = {0<=x<=3/2, max(0,x-5/4)<=r<=7/4-x},
    X(V)=max(sup |V|_H/3, sup |V_x|_H, sup |V_r|_H).

Let V be C1 on this closed domain (one-sided at its boundary), and let V=0
on x+r<1. Continuity gives V=0 at the edge x+r=1. Then

    |V(x,r)|_H <= (x+r-1)_+ X(V) <= (3/4)X(V).       (A)

For x<=1 and x+r>=1, integrate V_r from (x,1-x) to (x,r).
Both endpoints and the whole segment lie in D+. For x>1, r<=7/4-x<3/4;
integrate along (1,0)->(1,r)->(x,r). The first segment lies in D+ because
r<=3/4. The second lies in D+ because its largest lower radial bound and
smallest upper radial bound occur at the final time x. Its total length is
r+x-1. The vector-valued fundamental theorem of calculus proves (A).
The remaining zero-wedge case is immediate.

The boundary history is exactly 0<=x<=5/4. There,

    |V(x,0)|_H <= (x-1)_+ X(V) <= (1/4)X(V).        (B)

The two lemmas hold for differences of compatible candidates and exact
solutions. They do not hold for arbitrary C1 functions with nonzero data
on the wedge. No hypothesis that V has zero initial data on the entire
initial footprint is needed: its nonzero initial packet starts at r=1.

The norm estimates from differentiated endpoints of the reflected kernel
are ||K f||_X <= (3/2)sup|f|_H and ||L g||_X<=sup|g|_H. Reflection is
an H-isometry. With the retained bounds ||A||,||B||<=10^-4,
||D0||<=10^-6, ||M||<=1/2000, (A),(B) therefore give

    ||T V||_X <= [(3/2)(2*10^-4+(3/4)*10^-6)
                    +(1/4)(1/2000)]X(V) = q_w X(V).             (C)

T preserves the wedge and the homogeneous Dirichlet condition. Thus the
Neumann-series and residual arguments apply on this closed subspace.
The sharper constant applies to the full retained Hilbert space as well
as its five-mean invariant subspace. It uses neither angular cancellation
nor any claim about the numerical approximation of coefficients.

At the target (3/2,1/4), (A) gives mean-component error <=3E/4 when a
compatible reconstruction has X error E. The generic factor 3 remains
valid; this separate improvement uses the exact wedge and still cannot
resolve a 5.58*10^7 auxiliary margin at E=10^29.

## 2. Initial/free-wave gate and complete two-term ledger

Let the candidate initial displacement be C1 and the candidate initial
velocity be continuous. They must use the same parity reflection and be
identically zero on a neighborhood of r=0 and on the required initial
zero interval. Suppose uniform H errors on 0<=r<=7/4 satisfy

    |U0-hatU0|<=delta0,
    |partial_r U0-partial_r hatU0|<=deltar,
    |V0-hatV0|<=deltav.

Differentiating the d'Alembert endpoints gives exactly

    e_free <= max((delta0+(3/2)deltav)/3, deltar+deltav).        (D)

The velocity integral has length at most 2x and prefactor 1/2; its value
error is at most x*deltav. Both first derivative errors are at most
deltar+deltav. Only one radial derivative of displacement and no radial
derivative of velocity is needed. Reflection is an isometry; the flat
initial collar removes all transition defects. An interval cell table
alone is not such a differentiable candidate, and its width alone is
not a lower bound proving that every possible approximation must fail.
Fixed shell uncertainty should be propagated through shared background
expressions before any accuracy claim or refined-grid decision.

For B0,B1 in the same compatible class, let e0>=||B0-b||_X and
e1>=||B1-TB0||_X. The exact identity is

    U-(B0+B1) = T^2 U + (I+T)(b-B0) + (TB0-B1).

It yields

    ||U-(B0+B1)||_X <= q_w^2*(1.03*10^33)+(1+q_w)e0+e1.       (E)

When B0 approximates U_free alone, charge ||S||_X to e0 once. If the
initial/free-wave reconstruction error is efree_rec, take
e0=efree_rec+52162547172519/2. B1 must approximate T B0 using all five
mean couplings, the actual coefficient functions and the physical Robin
trace. Errors of coefficient approximation times the potentially huge
B0 and its first derivatives belong in e1. Changing the background
functions in T without charging these errors changes the equation.

Alternatively a full residual r=||V-b-TV||_X gives r/(1-q_w). The residual
already includes the omitted iterations. Adding a tail again is unnecessary
double accounting, though it is a weaker upper bound rather than a logical
contradiction. Normal/optimized arithmetic reports distinguish these routes.

## 3. Mean closure and the registered initial one-jet

The invariant tensor means are z=(K0,L0,N0,J0,w0), with Gram matrix
diag(1,2,1,3,1). K is yy, L is y0, N is 00, J is the coefficient of the
sphere metric, and w is scalar. Only L0 is Dirichlet. Rotational averaging
commutes with all linear bulk/boundary maps. First-order nonmean modes
contribute quadratically to mean sources before projection. The retained
operator primitive matrices have no off-block entries to or from these
five slots. A scalar-only mean is not closed under the retained system.

The reviewed initial adapter gives, in this order,

    U0 = (0, 0, -u^2/3, f0, -3(f0)_y/(2p)),
    (V0)_K = ell*(2rho Fmom/3-8u ut)/6,
    (V0)_L = ell*rho*((-26u uy-2a u^2)/6-3(f0)_y/2-3a f0),
    (V0)_N = ell*(-8rho Fmom/3+44u ut)/6,
    (V0)_J = ell*(2rho Fmom/3+4u ut)/6,
    (V0)_w = ell*(rho*(Gmom-(Fmom)_y)/p+2u vt)/6.

Here Y^2 has sphere probability mean 1/6, f0 is the invariant angular
metric coefficient, and ell=10^-6 multiplies every ordinary t velocity.
These formulas follow directly from source_adapters.initial_frame_data
and linear_operator.pack_initial, and agree with the retained harmonic
completion after division by frame factors at t=0. In particular the
mixed coefficient contains the second variation of gamma*dot(beta),
and the lapse choice N2=0 gives k00=-2u^2Y^2 rather than zero.

The physical initial construction uses p bounded away from zero. For
lambda=0,8,24, it sets
f_l=-rho^2 S_l/[2(lambda-3)], zeta_l=-3(f_l)_y/(2p).
The linear Hamiltonian radial terms cancel by p_y+4ap=V', leaving
-S_l. The three denominators are -6,10,42. Its momentum choice
K2=(Fmom/3)Y^2 gamma0 and Pi2=(Gmom-(Fmom)_y)Y^2/p gives
-Fmom d_sphere(Y^2)-Gmom Y^2 dy because space has dimension four.
These identities, not finite point-jet checks, supply the all-point
initial constraint proof. The harmonic completion then enforces C2=0
initially without altering those physical data.

## 4. Physical subsidiary bridge and limits of a mean evaluation

The retained completion is an existing conditional theorem, not a missing
construction. The following joins are necessary and are consistent in the
reviewed texts and adapter.

1. Use the full lowered quadratic gauge covector:
   C2=dC[k]+QC. Lowering varies with the metric. Its lower-order cross
   terms vanish only because the actual transformed first field has
   C1=0 identically. A small interval for C1 is not that identity.
2. The explicit harmonic initial velocities give C2=0. Tangential
   derivatives then vanish on the initial slice. For X=grad_n C2,
   the reduced equations give E_nn=X_n, E_nI=X_I/2, E_IJ=0,
   tr E=-X_n. Hence trace-reversed normal projections are X_n/2 and
   X_I/2. The physical Hamiltonian and momentum constraints imply
   X=0. The reduced PDE determines the second time accelerations;
   they are not independent initial data. Lower-order frame/normal
   variations multiply already vanishing lower subsidiary jets.
3. The normal primal row uses GN=-2 tr_gamma QR-2QC_y and gives
   (C2)_y=0. The physical shell rows and exact Codazzi/matter identity
   give E2_y mu=0. Therefore the coordinate tangential covector obeys
   partial_y C2_mu-2a C2_mu=0. Since a tangent covector coordinate
   component has one rho frame factor, this becomes
   partial_y C2_frame_mu-a C2_frame_mu=0. A covariant two-tensor has
   two factors; that separate conversion explains the zero tangential
   homogeneous diagonal in the parallel-frame primal Robin row.
4. Differentiating the exact trace-reversed Bianchi identity twice,
   using the exact lower-order equations and scalar equation, gives
   (Box0+Ric0)C2=0. Qualitative smoothness of the retained actual
   fields is needed here. It follows from their smooth data and finite
   characteristic/trace-ODE construction; it is not inferred from
   numerical C3 source boxes. Those boxes alone support a C1 response
   estimate, not a quantitative higher-derivative bound.
5. In the finite angular subsidiary system, outgoing characteristics
   feed the finite Robin trace ODE. Incoming characteristics are then
   fixed by the trace ODE and normal Dirichlet row. Local uniqueness
   follows by a short-interval contraction and continues over the
   bounded cone, regardless of Robin energy sign. Zero initial and
   boundary data imply C2=0. Natural Einstein/scalar, copy-even shell,
   induced-match and retained force/conservation identities follow.

The full theorem concerns exact total U. A five-mean numerical candidate
and C1 error bound evaluate only its five means; they do not evaluate all
angular constraints. C1 error controls first-derivative gauge differences
after the physical factors 1/ell and 1/(ell rho), frame connections and
angular operators are included. It does not control pointwise second
derivative Einstein residuals. Such a claim needs quantitative C2 control
or an explicitly specified weak norm. No invariant observable or finite
amplitude nonlinear conclusion follows from this bridge.

## 5. Reproduction and input scope

Run check_local_total_review.py with --research-root pointing to a tree
containing the four selected capsule paths and their external pin files.
The checker authenticates every declared member of C,G,P,N and rejects
undeclared files before importing fractions or any scientific code. It
authenticates the same complete closures after arithmetic. It imports
no retained producer. Normal and -O reports have identical mathematical
content; output paths are not embedded. Simple wrong-formula distinctions
are exact synthetic examples, not numerical mutation subprocesses or
proofs of the retained analytic PDE premises.

Primary exact files are the follow-on plan and the four closures named
there: C/inputs/FULL_SECOND_RESPONSE_CONSTRUCTION.md;
C/inputs/FORCED_INCREMENT_BOUND.md; N/response/REFINED_RESPONSE_BOUND.md;
G/projection/inputs/source_adapters.py and linear_operator.py;
G/response-signs/inputs/INITIAL_PROFILE_PROOF.md;
P/inputs/EXPLICIT_COMPACT_SECOND_INITIAL_JET.md;
P/inputs/EXPLICIT_INITIAL_HARMONIC_COMPLETION.md; and
P/boundary-review/BOUNDARY_COMPATIBILITY_AUDIT.md. Their full manifests
and fixed external manifest-pin records are checked by the local checker.
No previous large producer campaign is rerun.

For a small standalone replay, check_portable_local_total.py authenticates
this note and EXPLICIT_PREMISES.json, then checks only the rational
implications and synthetic formula distinctions. The small premise file
explicitly assumes the retained PDE/operator/source facts. This portable
route neither reauthenticates the four upstream capsules nor replays their
scientific certificates. The full local checker and its separate receipts
record that larger byte-authentication join.
