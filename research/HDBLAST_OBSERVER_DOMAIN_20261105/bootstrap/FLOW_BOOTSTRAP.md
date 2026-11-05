# Conditional observer ODE continuation from fixed coefficient data

5 November 2026, UTC. This note closes a finite-dimensional logical gap in
the earlier observer derivative budget. It gives a sufficient condition for
existence over a prescribed past duration before applying parameter-flow
estimates. Every coefficient, radius and initial-data bound below is an
input. No numerical amplitude, Einstein–scalar existence theorem, physical
shell embedding theorem, or global injectivity theorem is established here.

## 1. Fixed neighborhood theorem

Use one fixed finite-dimensional ambient realization and one fixed scaled
norm on its state space R^d. Let N be the compact set of labels; for the
application N=S^3, represented with fixed smooth frames or finite local
trivializations. Fix a>0, T>=0 and r>0. A known reference solution

    X0_s(s,n)=F(0,s,X0(s,n)),       0<=s<=T, n in N,

is assumed continuous in (s,n) and continuously differentiable in s. It is
the background solution, not an assumed family of perturbed trajectories.
Set

    K_r = { (e,s,X): |e|<=a, 0<=s<=T,
                       ||X-X0(s,n)||<=r for some n in N }.

This is compact. Assume K_r is contained in an open regular coefficient
domain U. The vector field F is continuous in s and jointly C^2 in (e,X),
with those derivatives continuous on U; endpoint clock regularity can instead
be understood by a local continuous extension in s. The closed parameter
interval must have an open parameter neighborhood for the usual C^2
endpoint interpretation. No F_s bound is used. A fixed label dependence in
F can be included as an extra continuous external parameter with uniform
bounds; it does not affect the proof.

Require on ALL of K_r, including every fiber ball and every e in [-a,a],

    ||F_X|| <= L,             ||F_e|| <= P1,              (1)

where L,P1>=0. Bounds only on the reference curve, only at e=0, or only on
already existing perturbed trajectories do not supply (1). The fiber ball
is convex, so the whole state segment from X0(s,n) to X is available. The
whole parameter segment from 0 to e at X0(s,n) is available as well.

Let the initial data X_in(e,n) be continuous, equal X0(0,n) at e=0, and obey

    ||X_in(e,n)-X0(0,n)|| <= a v_disp,     |e|<=a,         (2)

uniformly in n, where v_disp>=0. Define

    phi_L(s) = integral_0^s exp(L(s-z)) dz
             = (exp(Ls)-1)/L for L>0;   phi_0(s)=s,
    B(s) = exp(Ls) v_disp + phi_L(s) P1.

If

    a B(T) < r,                                          (3)

then for every |e|<=a and every n the initial-value problem has a unique
solution on [0,T], and

    ||X_e(s,n)-X0(s,n)||
       <= a exp(Ls) v_disp + |e| phi_L(s) P1
       <= a B(s) <= a B(T) < r.                          (4)

Thus the perturbed family is obtained by local existence and compact
continuation. Its existence over [0,T] is not a premise. If the stronger
initial estimate ||X_in(e,n)-X0(0,n)||<=|e| v_disp is supplied, the first
right side in (4) improves to |e| B(s). Both variants yield (3).

### Proof

The inequality (3) implies a v_disp<r, so all initial values are interior
to the prescribed coefficient neighborhood. The local Picard–Lindelöf
theorem applies because F_X is continuous and locally bounded. Fix e,n and
consider a local solution up to its first possible exit from the radius-r
fiber or its maximal existence time. While it remains inside the fiber,
the fundamental theorem of calculus on the state and parameter segments
gives

    ||F(e,s,X_e)-F(0,s,X0)||
      <= ||F(e,s,X_e)-F(e,s,X0)||
         + ||F(e,s,X0)-F(0,s,X0)||
      <= L ||X_e-X0|| + |e| P1.

Subtract the integral equations and apply the integral Gronwall lemma.
This proves (4) on every such local interval. The function B is
nondecreasing because all inputs are nonnegative. A first exit at a time
<=T would, by continuity, have distance r, contradicting a B(T)<r.

To exclude a finite maximal existence time b<=T, choose the fixed number

    R=(r+a B(T))/2 < r.

The trajectory up to b is confined to the compact set K_R, which has a
positive distance from the complement of U in any fixed product norm.
Continuity bounds F on K_R. Therefore X_e(s,n) is uniformly Lipschitz in s
there and has a limit as s increases to b. The limit is still in K_R.
Local existence at that limit and uniqueness continue the solution if b<T;
if b=T the limit supplies its endpoint. Both alternatives exclude failure
to exist on [0,T]. This argument is uniform in the label and parameter:
the same compact K_R and its coefficient-domain margin work for all of
them. At T=0 the conclusion is the initial-data assertion alone. This
proves the theorem.

This proof uses boundedness of F on a compact set as well as the displayed
derivative bounds. That boundedness follows from the stated continuity and
compact-domain premise; it must not be assumed on a noncompact union of
charts or on a family whose regular domain degenerates with e.

## 2. Parameter differentiation comes after continuation

Assume now that X_in is C^2 in e, with continuous uniform derivatives, and
that on the same full K_r the additional bounds hold:

    ||F_XX||<=H,  ||F_Xe||<=C,  ||F_ee||<=P2,
    ||partial_e X_in||<=V0,  ||partial_e^2 X_in||<=W0.    (5)

All bounds are nonnegative. Standard C^2 dependence of finite-dimensional
ODEs on parameters, applied to the continued solution, gives

    V=partial_e X,        W=partial_e^2 X,
    V_s=F_X V+F_e,
    W_s=F_X W+F_XX[V,V]+2F_Xe V+F_ee.                   (6)

For all e,n,s in the prescribed interval,

    V_* = exp(LT) V0 + phi_L(T) P1,
    ||V||<=V_*,
    W_* = exp(LT) W0 + phi_L(T)[H V_*^2+2C V_*+P2],
    ||W||<=W_*.                                         (7)

Here is the parameter-regularity argument under the stated weaker clock
hypothesis. Around any solution point choose a closed local state ball and
a short clock interval inside U, with bounded F_X and enough room for the
integral equation. The Volterra map on the sup-norm space of continuous
paths is a uniform contraction when the interval length times the local
F_X bound is less than one. This yields local continuous dependence on the
initial state and e. Subtract two parameter integral equations and divide
by their parameter separation. The resulting difference-quotient equation
has averaged F_X and F_e along the state and parameter line segments. The
segments lie in that local coefficient neighborhood; continuity of the
derivatives and the continuous dependence just proved give uniform
convergence of the averaged coefficients. Integral Gronwall then gives
uniform convergence of these quotients to the unique solution of the
first variational equation in (6).

Repeat this difference-quotient argument for the coupled continuous
integral equations for (X,V). The derivatives of its right side use only
F_XX,F_Xe,F_ee and the already continuous V; they are continuous and
bounded on the local compact set. Equivalently subtract the V equations,
split each difference in its state and explicit-parameter variables, and
use their averaged first derivatives. The two mixed contributions, one
from differentiating F_X V and one from differentiating F_e, are equal by
the symmetry F_Xe=F_eX; they sum to 2F_Xe V. Gronwall gives convergence to
W and continuity of the second parameter derivative. Thus (6) is valid
without any derivative of F with respect to the clock. The derivatives
are understood as parameter derivatives of a continuous path at fixed s.

Because the continued solution stays a uniform positive distance inside
U, finitely many such short intervals cover [0,T]; uniqueness concatenates
the local derivatives. Endpoint parameters have local extensions in the
open parameter domain, and the strict margin permits their finite-duration
continuation; compactness of parameter and label gives uniform local
neighborhoods as needed. Continuity gives the closed-interval bounds.
Finally ordinary integral Gronwall, first for V and then for W, gives (7).
No preexisting global perturbed tube was used to establish (6)–(7).

The displacement input v_disp in (2) is NOT a derivative bound. A small
bounded function of e can oscillate rapidly. To reuse the earlier note's
single name v0 for both inputs, explicitly assume
||partial_e X_in||<=v0 and X_in(0)=X0(0); integration in e then supplies
||X_in(e)-X0(0)||<=|e|v0, and one may put v_disp=V0=v0. W0 is still an
independent second-variation input. If only (2) is known, V0 in (5) must be
provided separately.

## 3. Geodesic state, constraints and units

The augmented intrinsic state in the earlier note is X=(x,u,P,Q), with
u=partial_tau x future directed, P_A=partial_A x, Q_A=partial_A u.
For past proper duration s=-tau its field is

    x_s=-u,  u_s=Gamma[u,u],  P_s=-Q,
    Q_s=(partial_x Gamma)[P,u,u]+2 Gamma[u,Q].           (8)

The initial data are the actual section sigma_e, its future unit normal
u_e, and the label derivatives of BOTH. Bounds V0,W0 charge every one of
these components, including the variation of the normal and Q. A fixed
coordinate section does not generally have a parameter-independent unit
normal. The commutations partial_e P=partial_A(partial_e x) and their
second-variation versions are valid when the initial compatibility and
label regularity are supplied, and follow by differentiating the same
geodesic solution in its label. Label-Jacobi data cannot be independently
chosen without checking that compatibility.

The full physical shell fields must already be h_e=iota_e^*g_e and
varphi_e=iota_e^*phi_e on one fixed reference shell. This note does not
construct iota_e or supply missing parameter derivatives of that physical
embedding. Sufficient coefficient regularity is C^2 parameter dependence
in a uniform shell C^4 topology, bounded inverse metric, and corresponding
initial data. Because (8) contains partial_x Gamma, its F_XX can involve
four spatial derivatives of h_e. Bounds must include mixed parameter and
space derivatives, with the uniform norm comparisons needed to obtain
(1),(5); a bare regularity label does not evaluate them.

The ambient extension is auxiliary and parameter independent. Start with
a compact neighborhood of the KNOWN background state in a fixed reference
shell, its relevant jet bundles and label trivializations; embed that
finite-dimensional compact neighborhood smoothly into Euclidean space,
and extend the lifted field to an open neighborhood by fixed local
extensions and a partition of unity (or a tubular retraction). Choose r
and a so every entire closed fiber ball in K_r lies in the regular
extension domain. Include the embedding, retraction, frames and extension
constants in (1),(5). This does not posit a single physical shell chart or
a global physical embedding of a perturbed observer map.

Off-constraint points in those balls are allowed and their coefficients
must be bounded too. On the physical constraint manifold the extension
must agree with the intrinsic field and be tangent to it. Initial position
and velocity tangency, h_e(u,u)=-1, and the differentiated label constraints
must hold initially. Intrinsic local existence and uniqueness, or the
constraint evolution identities, then preserve those constraints on the
continued solution. An arbitrary ambient extension that changes the field
on the constraint manifold would not justify physical interpretation.
Metric compatibility preserves the unit-speed constraint for (8).

Use positive reference scales for x,u,P,Q before forming any product norm.
If a normed state unit is denoted [X], parameter unit [e] and clock unit
[s], then [L]=[s]^-1, [P1]=[X]/([e][s]),
[V0]=[X]/[e], [W0]=[X]/[e]^2,
[H]=1/([X][s]), [C]=1/([e][s]), [P2]=[X]/([e]^2[s]).
The radius r and a v_disp have the same scaled state unit; LT is
dimensionless. Projection constants converting V_*,W_* to coordinate
time or label derivatives carry the appropriate reference scales. An
unscaled sum of positions and velocities cannot be substituted here.

Finite shell/label atlases are an alternative. Choose chart chains and
switch times from background compactness, fixed independently of e.
Each closed background segment must have a chart-domain margin large
enough to contain its prescribed coefficient fiber balls. On segment j,
apply the theorem with duration T_j, initial displacement d_j, and guard

    D_j=exp(L_j T_j)d_j+a phi_{L_j}(T_j)P1_j < r_j.

At a parameter-independent chart transition Psi, require whole segment
derivative bounds ||D Psi||<=alpha, ||D^2 Psi||<=beta on its regular
transition neighborhood. Then

    d_{j+1}<=alpha D_j,
    V_{j+1}<=alpha V_{j,end},
    W_{j+1}<=alpha W_{j,end}+beta V_{j,end}^2.

State comparison uses the segment between the perturbed and background
endpoint in that transition neighborhood, or a separate Lipschitz bound
valid there. The transition domain must therefore contain that segment.
Prove one segment, transfer its endpoint, and prove the next; this is a
finite induction, not an assumption that the entire perturbed chart chain
already exists. Take maxima over all segments and chains for final bounds.
Taking the largest chartwise L while ignoring transition Hessians is not
a valid substitute for this induction. The single duration-wide theorem
uses the coherent fixed ambient realization to avoid these transitions.

## 4. Root horizon and the physical past slab

For the actual reference clock q=-tau/rho_b, rho_b>0 is the fixed exact
background radius. Let ell=10^-6 in the registered coordinate units, and

    T_q=25 ell/16,   T_coord=13 ell/8,
    s=rho_b q,      T_proper=rho_b T_q.                 (9)

The background coordinate time is t_0(q,n)=-q. Applying the proper-clock
theorem means using T=T_proper, NOT T=T_q or T_coord. Alternatively solve
the rescaled ODE

    X_q=G(e,q,X)=rho_b F(e,rho_b q,X).

In the same norm its constants are
L_q=rho_b L_s, P1_q=rho_b P1_s, H_q=rho_b H_s,
C_q=rho_b C_s, P2_q=rho_b P2_s. In particular

    L_s T_proper=L_q T_q,
    phi_{L_s}(T_proper)P1_s=phi_{L_q}(T_q)P1_q.

The same identities apply to second-variation forcing. rho_b is fixed,
so there are no derivatives of rho_b with respect to e. A family-dependent
radius used to redefine the clock would change these derivative formulas.

Suppose the coordinate-time projection and (7) give the uniform estimate

    sup_{|e|<=a, 0<=q<=T_q,n} |partial_e t_e(q,n)|<=Z_t,
    a Z_t < ell/32.                                    (10)

Because the full parameter interval has now been constructed, integration
in e legitimately gives |t_e(q,n)+q|<=a Z_t. Independently require the
actual anchor t_e(0,n)=0 and a uniform positive monotonicity bound

    -partial_q t_e(q,n) >= m_t > 0.                    (11)

This can be obtained from a sign bound on the time component of G over
the prescribed coefficient neighborhood, or proved on the constructed
flow from separate input. It does not follow just from position closeness.
Then q>0 implies t_e(q,n)<0, and (9),(10) imply

    t_e(q,n) >= -T_q-a Z_t > -T_coord+ell/32.

Thus the actual flow stays in [-T_coord,0]. The UNPERTURBED horizon margin
is T_coord-T_q=ell/16; (10) consumes less than half and leaves a PERTURBED
lower-bound margin greater than ell/32. It would be incorrect to claim an
unchanged ell/16 perturbed margin from (10).

A convenient sufficient replacement for the separate position and sign
inputs (10),(11) uses the proper velocity retained in (8). Suppose its
time projection has uniform parameter-derivative bound

    |partial_e u_e^t|<=U1,   u_0^t=1/rho_b,
    c=rho_b a U1 < 1/50.                                 (12)

For example, U1 is the fixed projection constant times V_* from (7).
The coefficient bounds determine that symbolic majorant before any
perturbed full-duration flow is assumed. The theorem first constructs the
flow in the auxiliary coefficient domain, and only then applies (7) and
integrates in parameter. Since t_q=-rho_b u_e^t, (12) gives

    |t_q+1|<=c,  -t_q>=1-c>49/50,
    |t_e(q,n)+q|<=c q<=c T_q<ell/32,

using t_e(0,n)=0 and (25/16)/50=1/32. Hence (12) supplies both required
clock conclusions. It does not require a separately computed time-position
projection budget. The same full-parameter and extension hypotheses still
apply. The sequence is coefficient nonexit, continuation, parameter-flow
regularity, then the physical-slab guard; there is no use of a presumed
physical full-duration observer flow to obtain these steps.

For the older retained coarse root bracket the past endpoint has
q=6145 ell/4096<T_q; its distance from T_q is 255 ell/4096. The full
connecting interval [0,T_q] is covered. This note supplies no new moment
endpoint signs and does not transfer the old threshold 106 to a different
fine linear bracket.

If the physical PDE supplies coefficients only on the closed slab
[-T_coord,0], whole ambient balls at the anchor t=0 can extend into t>0.
One must supply a regular extension of the COEFFICIENTS near that boundary,
with joint parameter/state regularity and all the stated uniform bounds.
A suitable one-sided extension theorem or an explicit construction is an
additional analytic input; it must respect the intrinsic data on the
physical side. No PDE solution beyond t=0 is asserted by this device.
Run the local/continuation theorem for that auxiliary extension, require
(11) on the relevant extended neighborhood or initially extend its sign
from a strict physical margin, and use (10),(11) to prove that the
constructed trajectory actually stays in the physical slab. Hence the
extended and physical equations agree along it. This avoids using an
assumed physical perturbed tube to justify the very first local solution.

## 5. What remains conditional

The retained source already asserts uniform first B13 and second B12
variation bounds, actual variation limits, and C^2 dependence into X10
(hence a lower piecewise spacetime C^4 topology). Its constants and common
parameter radius remain unevaluated. These are more than expansion-point
tangent statements. Their registered past-slab application and the
norm/pullback/extension comparisons needed for this particular observer
field must still be supplied. This finite-dimensional theorem does not
repair or replace any PDE premise.

For this observer bootstrap one still needs: a registered exact physical
past solution family and full physical shell pullbacks on the same fixed
shell; a fixed regular coefficient neighborhood of the known background
with compact margin; the full-ball uniform constants and initial section,
normal and Jacobi data; compatible auxiliary extensions and norm constants;
and values satisfying (3),(10),(11). None is numerically evaluated here.
The full initial/PDE dependence may make these constants functions of a;
all guards must then be checked at the same candidate a.

Continuation and C^2 label transport do not by themselves give embedded,
noncrossing spacelike cuts, global injectivity of the normal exponential
map, a positive induced cut metric, or a physical foliation. Those require
separate geometric hypotheses or a separately proved observer-domain
argument. The quotient/area/moment budgets require those inputs in
addition to (7). Positive local cut metric alone must not be promoted to
global injectivity. No effective a is claimed, and no legacy artifact was
changed or published.

## 6. Exact finite checks

The companion checker uses Python's standard library and Fraction. It
checks the split fundamental-theorem identity on polynomial examples,
formal exponential/phi coefficients, the inhomogeneous Gronwall majorant,
the parameter-chain coefficients, clock-rescaling identities, dimensions,
strict-guard edge examples and exact horizon arithmetic. These finite
algebra checks do not prove Picard–Lindelöf, analytic Gronwall, compact
continuation, a coefficient extension theorem, or the physical/PDE inputs.
The analytic proof above uses those standard finite-dimensional facts.

    python3 check_flow_bootstrap.py --output /absolute/fresh/normal.json
    python3 -O check_flow_bootstrap.py --output /absolute/fresh/optimized.json

The output file must be fresh and absolute with an existing non-symlink
parent. The checker uses exclusive creation, emits the same canonical
JSON bytes to the file and stdout, and does not depend on assert statements.
Normal and optimized reports must be byte-identical. See CHECK_REPORT.md
and MANIFEST.json for the local receipts and source hashes.
