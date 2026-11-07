# Explicit retained sixteen-field energy majorant after Run06f

6 November 2026. Conditional mathematical checkpoint. This continues the actual
retained mixed-wave theorem; it does not replace that theorem by an assertion
that no evolution theorem exists. All conclusions below are conditional on
the physical coefficient patch and boundary structure holding at the points
and frozen iterates under consideration. No frozen predecessor was edited.

The quantitative result is (10): all multipliers of the explicitly defined
first-jet and lower-order budgets are evaluated. Neither those solution
budgets nor the complete second-wave forcing are bounded numerically by the
retained pointwise metric patch. The result is an auxiliary analytic energy
estimate, with no interpretation as physical gravitational energy.

## 1. Actual retained system and fixed geometry

The sources are the prior public metric certificate and these
Run06e shell-evolution inputs:

* `QUASILINEAR_VARIABLE_MIXED_EVOLUTION.md`, sections 1–4;
* `LINEAR_VARIABLE_WAVE_AUDIT.md`, basic energy and normal recovery;
* `ANGULAR_UNIFORM_LINEAR_CONTACT.md`, sections 2–5;
* `EXACT_FINITE_REGULARITY_SHELL_EVOLUTION.md`, sections 1, 2 and 4;
* `FULL_FINITE_ANGULAR_SECOND_ORDER_CONTACT.md`, equations (9)–(10).

The actual wider-collar background box is also pinned through Run06f
`inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json`.

The coefficient convention is the inward conformal radial coordinate
`dy=−rho(z) dz`, unscaled de Sitter time, and

    ell = 10^−6,  T = 13 ell/8 = 13/8000000,
    M = [0,Z] × S^3,  Z = 6 ell,
    dmu_ref = dz dmu_kappa,  vol(S^3,kappa) = 2 pi^2.

Use the actual sixteen independent reduced coordinates after the exact
pointwise boundary flow: four radial cross components form E_D; one radial
root, ten independent symmetric tangential metric coordinates and one scalar
form E_N. Thus E=E_D orthogonal-sum E_N has ranks 4 and 12. Choose the fixed
positive fiber metric and split-preserving metric product connection of the
retained proof. In a symmetric-tensor matrix representation, off-diagonal
entries require the corresponding sqrt(2) normalization if the Hilbert norm
is the Frobenius norm. One must not sum all sixteen metric matrix slots as
independent tensor unknowns. The total number of independent fields is 16.

Choose a round left-invariant orthonormal angular frame E_A, A=1,2,3,
normalized by [E_A,E_B]=2 epsilon_ABC E_C. Each E_A has zero round
divergence. Put

    s(t)=sech(t), h(t)=tanh(t),
    X_z=partial_z, X_A=s E_A,
    q=D_t u=u_t, W_i=D_(X_i) u.

Every X_i is divergence-free for the fixed reference spatial volume. The
bundle product connection has no mixed time/angular curvature, so

    D_t W_z=D_(X_z)q,
    D_t W_A=D_(X_A)q−h W_A.                         (1)

This time dependence must be included even though the volume is fixed.

After changing the equation sign as in the retained proof, use the explicit
ordering

    P u = a u_tt − 2 b^i D_(X_i)u_t
                       − c^{ij}D_(X_i)D_(X_j)u,
    L u = P u + L_t q + L_w(W) + L_0 u = f,          (2)

where a=−g^{tt}, b^i=g^{ti}, c^{ij}=g^{ij} are components in the dual of
`(dz,dt,cosh(t)theta^A)`. L_w is the full map E^4→E; it is not a scalar
coefficient estimate applied separately to sixteen components. The scalar
principal coefficients multiply the identity on all of E. The bounded
matrix coefficients L_t,L_w,L_0 include the actual reference-connection
terms, the state-flow chain-rule terms, and any terms needed to write the
actual Hessian in the ordering (2). They need not be symmetric.

For example, if a symmetrized time/angular derivative is used originally,
rewriting it in the order D_(X_A)D_t uses (1) and creates a first-order
term. It belongs in L_w. For the fixed round spatial product connection,
the frame-Hessian correction itself cancels: its angular part is
`sum_AB c^{AB} s epsilon_ABC D_(X_C)`, which is zero by symmetry of c.
The full physical/reference connection difference does not thereby vanish.
No numerical bound for those actual lower-order maps is assumed here.

## 2. Separate positivity conclusions and energy

The Run06f physical patch is

    ||Gamma−J_5||op <= 1/40,   129 <= rho <= 131,
    g=rho^2 Gamma,  J_5=diag(1,−1,1,1,1).

The resolvent identity gives ||Gamma^−1−J_5||op<=1/39. Applying it to
the temporal coordinate vector and, separately, every purely spatial
vector proves

    a >= 38/(39*131^2) > lambda=1/18000,
    c >= lambda I,
    ||g^−1||op <= 40/(39*129^2) < mu=1/16000.        (3)

In particular |a|,|b|_2,||c||op<=mu. Positivity of the full spatial
inverse block is not inferred merely from the timelikeness of dt.

Use exactly the retained energy

    e[u]=1/2 (a|q|^2 + c^{ij}<W_i,W_j> + |u|^2),
    E[u](t)=integral_M e[u] dmu_ref.

For all sixteen fibers together,

    1/2 [lambda(|q|^2+|W|^2)+|u|^2] <= e[u],
    Q_fixed/72000 <= e[u] <= Q_fixed/2,             (4)
    Q_fixed=|q|^2+|D_z u|^2+|D_Omega u|^2+|u|^2.

The second line uses s(t)^2>1/2 on |t|<=T. Consequently the integrated
fixed-product first-order Hilbert norm is less than 269 sqrt(E[u]).
No extra factor sixteen occurs: (4) uses the full positive Hilbert sum.

## 3. Exact identity, including moving frame and all boundary fluxes

First suppose P u=F. Metricity, spatial integration by parts, and (1) give

    E' = integral <F,q> + integral <u,q>
       + integral (a_t/2 − sum_i X_i b^i)|q|^2
       + 1/2 integral c_t^{ij}<W_i,W_j>
       − integral (sum_i X_i c^{ij})<W_j,q>
       − h integral sum_ij c^{ij}1_(j angular)<W_i,W_j>
       + integral_boundary [b^i n_i |q|^2
                           + c^{ij}n_i<W_j,q>].     (5)

Coefficient derivatives in (5) are derivatives of the displayed component
functions. If they are instead expressed using a geometric covariant
derivative, the resulting frame-connection conversion terms must also be
included; they cannot be omitted from the budget. There is no b_t term:
the time-space principal term is integrated in its spatial index. No
spatial frame commutator is used in the base cancellation. The last bulk
term can equivalently be written with the symmetric matrix
`−h(c Pi+Pi c)/2`, Pi=diag(0,1,1,1).

At the shell every frozen coefficient field must obey all four exact
coordinate cross-trace identities g_za=0, a=t,Omega. Then inversion gives

    b^z=0, c^{zA}=0 at z=0.

For the homogeneous transformed rows u_D=0 and D_z u_N=0, split
preservation gives q_D=0 and W_(z,N)=0. The shell flux therefore equals

    −c^{zz}(<W_(z,D),q_D>+<W_(z,N),q_N>)=0.        (6)

This uses exact 4/12 orthogonality and the scalar principal symbol. It
does not require a symmetrizer of the old coupled Robin boundary matrix.
Interior b^z,c^{zA} are allowed and must remain in (2)–(5).

At z=Z, (5) also has a flux. For this estimate impose the retained
auxiliary homogeneous Dirichlet wall on all components, so q|_Z=0, or
separately establish that the supported solution vanishes in a neighborhood
of that wall. Merely naming Z=6ell does not eliminate its flux. The
conditional support arithmetic L+2T=23ell/4<6ell with L=5ell/2 and c*=2
is consistent with the wall construction. It is not an evolved support or
metric-patch proof for the actual family.

## 4. Evaluated first-jet majorant

At each time define nonnegative coefficient budgets

    A_t = ||partial_t a||_(Linf),
    B_d = ||sum_i X_i b^i||_(Linf),
    C_t = ||partial_t c||_(Linf,op),
    C_d = ||(sum_i X_i c^{ij})_j||_(Linf,l2),
    l_t = ||L_t||_(Linf,op),
    l_w = ||L_w:E^4→E||_(Linf,op),
    l_0 = ||L_0||_(Linf,op).                        (7)

The supremum is over the actual spatial cylinder where (3) and (6) hold.
These budgets include time and spatial coefficient derivatives. They are
not supplied by (3). A sufficient component-derivative bound is obtained
by replacing each divergence norm with its sum of component derivative
norms; that weakens the estimate without adding any unknown constant.

From (4), at each point and after integration,

    |<u,q>| <= E/sqrt(lambda) < 135 E,
    ||q||_2 <= sqrt(2/lambda) sqrt(E) < 190 sqrt(E),
    integral |W||q| <= E/lambda = 18000 E,
    integral |u||q| <= E/sqrt(lambda) < 135 E.

Here each inequality involving E denotes the corresponding integrated
quantity. The moving term is bounded by
`2 mu |h| E/lambda=(9/4)|h| E`. Put

    K(t)=135 + 18000 A_t + 36000 B_d
              + 18000 C_t + 18000 C_d
              + 36000 l_t + 18000 l_w + 135 l_0
              + (9/4)|tanh(t)|.                   (8)

Substitute F=f−L_tq−L_w(W)−L_0u in (5). With (6) and the stated outer
condition,

    E'(t) <= K(t) E(t) + 190 sqrt(E(t)) ||f(t)||_2. (9)

No derivative of the force and no coefficient b_t is required for this
zeroth-order estimate. A regularization at E=0, followed by variation of
constants, yields the explicit propagator

    Y(t) <= G(t,0)Y(0) + 95 integral_0^t G(t,s)||f(s)||_2 ds,
    Y=sqrt(E), G(t,s)=exp(1/2 integral_s^t K(r)dr). (10)

In particular it is bounded by G(t,0)[Y(0)+95||f||_(L1_tL2_x)].
For t<0 apply the same calculation in the reverse time direction and use
the absolute-value budgets. On the specified positive interval, the
moving-frame contribution to K is at most 117/32000000, and its
contribution to log G(T,0) is at most

    (9/16)T^2 = 1521/1024000000000000.             (11)

The constant 135 is a conservative consequence of the retained added
mass term |u|^2/2. It is not an instability claim.

A useful version of (8) uses jets of the actual normalized physical
metric, rather than inverse-coefficient derivatives. In the specified
global invariant moving coframe define

    U_t=||partial_t Gamma||_(Linf,op),
    U_x=||partial_z Gamma||_(Linf,op)
                            +sum_A||X_A Gamma||_(Linf,op),
    R_rho=||rho_z/rho||_(Linf).

These are component derivatives in that frame, with no angular connection
conversion suppressed. For any derivative Y other than partial_z,
`Y(g^−1)=−rho^−2 Gamma^−1(Y Gamma)Gamma^−1` as a matrix of components.
The partial_z derivative additionally contains `−2(rho_z/rho)g^−1`.
Since

    rho^−2 ||Gamma^−1||op^2
             <= 1600/25310961 < 1/15000,

block restriction and the triangle inequality give

    A_t,C_t <= U_t/15000,
    B_d,C_d <= U_x/15000 + R_rho/8000.

For C_d, use the norm of each differentiated row followed by the sum over
i; there is no unnecessary factor four. Consequently a sufficient actual
metric form of the majorant is

    K <= 135 +(12/5)U_t +(18/5)U_x +(27/4)R_rho
                +36000 l_t+18000 l_w+135 l_0 +(9/4)|tanh t|. (11a)

No numeric evolved U_t,U_x or actual lower-order budgets are thereby
provided. The background warp budget R_rho is distinct from a derivative
of an evolved unknown. This formula records an explicit route for using
future actual-family first-jet estimates with (10).

The actual retained background box does evaluate the warp budget. Write
`a_bg=rho_y/rho`, to distinguish this background variable from the wave
coefficient a=−g^{tt}. Its wider-collar bound is a_bg<=14/125 and
rho<=131. Since dy=−rho dz,

    R_rho=||a_bg rho||_(Linf) <= 1834/125 < 15,
    (27/4)R_rho <= 24759/250 = 99.036.

The wider box explicitly covers the rectangular enlargement with
0<=t<=3ell, 0<=z<=8ell, hence contains this candidate slab. Keeping
R_rho<15 and |tanh(t)|<=T gives the simple conditional version

    K <= 237 +(12/5)U_t +(18/5)U_x
                  +36000 l_t+18000 l_w+135 l_0.     (11b)

The numerical constant 237 is therefore an evaluated background/frame
contribution. The two actual evolved metric first-jet budgets and the
three complete transformed lower-order budgets still require estimates.

## 5. First and second retained waves

Write the actual transformed equation as P(V)V=F(x,V,DV), F(x,0,0)=0.
Let V_1=partial_epsilon V_epsilon|_0 and
V_2=partial_epsilon^2 V_epsilon|_0; this is the retained convention
V_epsilon=epsilon V_1+epsilon^2 V_2/2+o(epsilon^2).
Linearizing at the background gives the full background linear operator

    L_bg = P(0) − D F(0)[u,Du],
    L_bg V_1=0,
    L_bg V_2=Q_bg[V_1],                            (12)
    Q_bg[V_1]=D^2 F(0)[(V_1,DV_1),(V_1,DV_1)]
                              −2 P_V(0)[V_1] V_1.

Any fixed lower-order terms placed on the left are included consistently
in L_bg and the derivatives in (12). P_V(0)[V_1]V_1 denotes coefficient
variation followed by the second derivatives of V_1, not a derivative of
V_1 in field space. The factor two is required by the epsilon^2/2
convention. The potentially demanding term is V_1 D^2 V_1, precisely the
one retained in the parameter-variation proof. The nonlinear flow is
pointwise and invertible; its second field jets contribute to the full
source and must not be dropped.

Because the exact nonlinear transformed rows are homogeneous for every
epsilon, both V_1 and V_2 have exactly homogeneous 4D/12N rows. Thus, with
their actual completed reduced initial jets,

    Y_1(t) <= G_bg(t,0)Y_1(0),
    Y_2(t) <= G_bg(t,0)Y_2(0)
                   +95 integral_0^t G_bg(t,s)||Q_bg[V_1](s)||_2 ds. (13)

These are finite-time propagators once the integrals in them are bounded.
At the actual background Gamma_0=J_5 throughout the collar, so U_t=U_x=0
in (11b). In particular, if bar(l_t),bar(l_w),bar(l_0) are the actual
complete background lower-order suprema on 0<=t<=T, then

    G_bg(T,0) <= exp[3081/16000000
                         +(117/4000)bar(l_t)
                         +(117/8000)bar(l_w)
                         +(351/3200000)bar(l_0)].           (13a)

This evaluates all time and coercivity multipliers for the first/second
background wave propagator; those complete lower-order suprema and the
coupled quadratic forcing still require their own bounds.
The initial canonical H14/H13 estimates from Run06e do not
by themselves give numeric Y_j(0): the complete reduced time velocities,
initial harmonic gauge completion, and flow reconstruction must be bounded
in the same norm. The exact canonical second derivative also must not be
identified with the distinct earlier selected finite-angular correction.

Alternatively the earlier linear Robin conjugation produces a second
wave with quadratic Neumann data. Its source is the complete differentiated
physical/gauge source, including k_yy=2 nu_2+2 nu_1^2. The nonlinear-flow
description and that linear-conjugation description require their own
matching chain-rule terms; the quadratic boundary source cannot simply be
erased from one description without changing its bulk source and initial
data. The following lifting handles the latter description honestly.

## 6. Nonzero Neumann data: explicit retained lifting cost

Suppose the frozen operator has the same coefficient shell identities,
u_D=0, and −D_z u_N=g(t,Omega). Retain the published sufficient smooth
lift

    p(z)=z eta(z), eta(0)=1,
    W_D=0, W_N=−p(z)g(t,Omega), R=u−W.             (14)

Here eta is the retained smooth collar cutoff, supported away from the
outer wall. No bounds on its derivatives are asserted. Define

    P_j=||partial_z^j p||_(L2(0,Z)), j=0,1,2.

Then R has homogeneous rows and L R=f−L W. With angular covariant
derivatives in the fixed product bundle, introduce the explicitly defined
boundary norms

    g_0=||g||_2, g_t=||g_t||_2, g_tt=||g_tt||_2,
    g_O=sum_A||D_(E_A)g||_2,
    g_tO=sum_A||D_(E_A)g_t||_2,
    g_OO=sum_AB||D_(E_A)D_(E_B)g||_2.

All are on S^3 at the same time; the use of sums is deliberate so no
unrecorded finite-dimensional norm constant enters. Since s<=1 and each
component |a|,|b^i|,|c^{ij}| is at most 1/16000,

    ||L W||_2 <= (1/16000)[P_0 g_tt +2P_1 g_t+2P_0 g_tO
                                 +P_2 g_0+2P_1 g_O+P_0 g_OO]
                 + l_t P_0 g_t
                 + l_w(P_1 g_0+P_0 g_O) + l_0 P_0 g_0.       (15)

This follows by expanding the explicitly ordered operator (2). There is
no assumption that coefficients depend only on z; their uniform bounds
are taken before integrating in z. The same bound includes all twelve
Neumann components in their full fiber norm.

Likewise, the upper energy bound implies

    sqrt(E[W]) <= P_0 g_0
                    +(1/178)(P_0 g_t+P_1 g_0+P_0 g_O).       (16)

The energy square root is a norm for each fixed coefficient field. Apply
(10) to R and add (16):

    Y_u(t) <= Y_W(t)+G(t,0)[Y_u(0)+Y_W(0)]
                  +95 integral_0^t G(t,s)(||f(s)||_2+||L W(s)||_2)ds. (17)

Equations (15)–(17) display the lifting derivative cost instead of
claiming an arbitrary L2 Neumann-input estimate. They do not evaluate
P_j, the actual quadratic boundary data, or their time/angular jets.
Initial corner compatibility of R is still required. It is automatic for
the retained shell-flat initial sources when the actual constructed lift
uses those data; it is not automatic for arbitrary boundary inputs.

## 7. High-order statement that is actually justified quantitatively here

Let T consist of partial_t and boundary-preserving tangential angular
derivations in the retained fixed global split frame. Equivalently use
the split-preserving SO(4) lifts, including their known zeroth-order fiber
terms. Every word T^I preserves homogeneous boundary traces and the outer
Dirichlet wall. Put

    Y_m^2=sum_(|I|<=m−1) E[T^I u],
    F_m^2=sum_(|I|<=m−1) ||T^I f||_2^2,
    R_m^2=sum_(|I|<=m−1) ||[L,T^I]u||_2^2.

The sums are finite and include every ordered word used to define the
norm. Commuting the full operator and applying Cauchy–Schwarz to (9)
gives the fully evaluated, angular-cutoff independent statement

    Y_m(t) <= G(t,0)Y_m(0)
                 +95 integral_0^t G(t,s)[F_m(s)+R_m(s)]ds.   (18)

All principal coefficient derivatives, frame/connection commutators,
and lower-order matrix commutators are inside the explicitly defined
R_m. In a component frame, writing S=diag(1,s,s,s), the fixed principal
coefficients to commute are a, beta=S b, C=S c S. Therefore their time
jets include derivatives of s and may not be treated as constant.

Equation (18) is an exact commuted-energy estimate, not a closed numeric
H14 estimate. The retained theorem closes that estimate qualitatively
using finite Sobolev products and normal recovery. Solve the actual
uncommuted equation for c^{zz}D_z^2u and use

    1/c^{zz} <= 18000.                              (19)

Applying normal/tangential derivatives then recovers higher normals
without imposing false homogeneous boundary rows on D_z u. The retained
source cost is

    sum_(j=0)^(m−1)||partial_t^j f||_(L1_t H^(m−1−j))
      +sum_(j=0)^(m−2)||partial_t^j f(0)||_(H^(m−2−j)).       (20)

For the actual m=14 theorem, the finite tame indices are those stated in
the retained proof: a safe low index is 10 for orders 10,...,14, and 11
at the separate smoothing order 15. The theorem already supplies the
conditional X14 evolution and C2-parameter dependence into X10. A generic
bound requiring W13,infinity coefficient jets would not quantitatively
apply to X14 data in four spatial dimensions. No such stronger hypothesis
is silently substituted here. Evaluating the finite product, commutator,
normal-recovery and complete initial-jet constants at the stated Sobolev
indices is the remaining high-order quantitative task. An unspecified
constant C_m or the source residual R_m cannot be represented as a
completed numeric high-order propagator.

## 8. Exact obstructions to closing from the pointwise patch alone

First, timelike dt is not enough for (3). The inverse metric with z slot
+1, t/x block `[[-1,2],[2,−1]]`, and the remaining two slots +1,+1 has
one negative eigenvalue and g^{tt}=−1, but g^{xx}=−1. The full spatial
quadratic form in this energy is not positive. The Run06f patch excludes
this example.

Second, even the full Run06f pointwise patch does not bound (7). Let

    Gamma_N(t)=J_5+delta sin(Nt)e_t tensor e_t,
    0<delta<1/40, g_N=rho(z)^2 Gamma_N,
    a_N=rho(z)^−2/(1−delta sin(Nt)).

The pointwise patch holds for every N; all shell cross traces vanish.
At t=0, (a_N)_t/a_N=delta N. For the admissible frozen scalar-principal
operator with these coefficients, zero lower terms and initial u=0,
q a nonzero Neumann-sector section compactly supported away from both
ends, (5) gives exactly

    E'(0)/E(0)=delta N.                              (21)

Thus a numerical exponent uniform over all metrics in the pointwise
patch cannot follow from that patch. This is a coefficient-hypothesis
counterexample, not a claim that these oscillatory metrics solve the
Einstein–scalar equations or belong to the actual evolved family.

Third, first-order energy alone cannot bound the complete second forcing.
In an interior four-dimensional coordinate ball take
u_epsilon(x)=epsilon^−1 psi(x/epsilon), with a fixed nonzero smooth
compactly supported psi. Its H1 norm stays bounded while
`||(D u_epsilon)^2||_L2` scales as epsilon^−2. The same scaling occurs
for a nonzero u_epsilon D^2u_epsilon term. The actual high-regularity
control is therefore needed to turn (13) into a numeric second-wave bound.
This does not contradict the conditional finite-regularity theorem.

The completed output is the evaluated energy multiplier (8), propagators
(10), (13), and the explicit lifting estimate (15)–(17). A numeric
actual-family Gronwall exponent, high-order propagator, nonzero admitted
amplitude interval, or common proper-time/area observable domain has not
been obtained by the retained patch alone.
