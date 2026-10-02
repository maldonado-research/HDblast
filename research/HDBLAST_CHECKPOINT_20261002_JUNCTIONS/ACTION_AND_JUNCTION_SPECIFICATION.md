# Action and junction prerequisite for a common-action 5D quantum source

Prepared 2026-10-02 UTC against repository commit
`e17a01b428bb8049e919c42376ab0359e41d613c`. This is a bounded derivation,
executable algebra check and prospective experiment design. No new physical
mode, shell, bulk or radiation evolution is performed. No coupled solution,
physical matching measurement or new cosmological mechanism is established.
The independent review is a separate internal check, not external peer review.

## 1. Controlling inherited action and actual code

The controlling action is the published 22 September matter extension,
`hdblast/checkpoints/2026-09-22_scalar-profile-branch-and-matter_zenodo-22922928/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/matter/MATTER_EXTENSION_AND_RESIDUAL_VACUUM.md`,
lines 9–18 and 53–91. The coupled implementation inspected here is
`research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source/evolve_b1.py`,
`shell_matter.py`, `static_w.py` and `b1_symbolic.py`. File hashes are in
`evidence/SOURCE_PINS.json`.

Use signature (−++++), two identical retained interiors, and one shell action:

    S = sum_(a=1,2) { [1/(2 kappa5²)] int_Ma sqrt(−g)
                       [R5 − (nabla phi)² − 2 U(phi)]
                     + [1/kappa5²] int_Sigma sqrt(−h) K_a }
        − int_Sigma sqrt(−h) sigma_R(phi)/kappa5²
        + Gamma_Q,ref[h,x(phi); state] + DeltaGamma_loc[h,x(phi)]
        + S_other,matched.

`Gamma_Q,ref` denotes the causal, state-dependent quantum expectation-value
prescription (an in-in/closed-time-path formulation if expressed as an effective
action). An in-out vacuum action alone is insufficient to define causal
backreaction. Its local subtraction convention is the positive-reference
stress order four / variance order two prescription in the frozen smooth-FRW
checkpoint; physical mode functions and their state are additional data.
`DeltaGamma_loc` contains only finite changes relative to that reference
prescription. `S_other,matched` must be specified if any other matter exists.

The inherited conformal chart has

    ds5² = exp(2B(T,Z))(−dT²+dZ²) + exp(2A(T,Z)) dX²,
    shell Z=0, retained interior Z<0,
    n = exp(−B_b) partial_Z,
    K_munu = h_mu^A h_nu^B nabla_A n_B,
    d tau = exp(B_b) dT,  a = exp(A_b),
    H = exp(−B_b) A_T,b,  v = exp(−B_b) phi_T,b.

The normal points from each interior toward the shell. With this normal positive
pure tension gives positive spatial extrinsic curvature. Reversing the normal
requires reversing both K and n(phi); the equations below retain the inherited
normal. The doubling is in the bulk/GHY action, not in the shell quantum field.

The exact classical action family used by the code is

    W(phi) = 1−phi+phi³/3,
    U(phi) = W_phi²/2 − 2 W²/3,
    sigma_R(phi) = 2 W + delta[1+c phi+d phi²/2],
    c = 2/1.0357712571566784 − 4/3.

The earlier registered model has d=0. B1 uses the explicitly changed A1 tuned
model: delta=0.1 and d=−3.106933495673783. These branches must not be silently
identified. The specification is valid for the displayed family; the next
registration must select its existing branch explicitly. The values above are
in inherited model-length units. Restoring the length unit L0 means U/L0² and
sigma/L0 in the dimensionful action.

## 2. Common quantum stress and scalar source

For a single canonical minimally coupled real shell field chi,

    S_chi = −1/2 int sqrt(−h)[(D chi)² + x(phi_b) chi²],
    x(phi_b) = m0² + G²(phi_b−phi_star)²,
    x_phi = 2 G²(phi_b−phi_star).

Here phi is the dimensionless inherited bulk scalar; the canonical bulk scalar
is Phi=phi/kappa5. Thus G is dimension one and equals the inherited bar-g;
the code's dimensionless G is G/H0. Specify multiplicity separately; auxiliary
PV fields are regulators and do not add physical species.

For the causal quantum source in its chosen finite prescription define

    T_munu,Q = −2 delta Gamma_Q/(sqrt(−h) delta h^munu),
    S_x,Q = −delta Gamma_Q/(sqrt(−h) delta x) = Q_chi/2,
    J_phi,Q = −delta Gamma_Q/(sqrt(−h) delta phi_b)
            = x_phi Q_chi/2 = G²(phi_b−phi_star) Q_chi,
    T^mu_nu,Q = diag(−rho_Q,p_Q,p_Q,p_Q).

`Q_chi` is the renormalized field variance current with mass dimension two.
The old B1 symbol `Q` is a decay-transfer power with dimension five; they are
different observables. In this proposed free-field prerequisite there is no
decay-transfer term. The physical rho and J_phi have mass dimension four,
sigma has dimension one, and kappa5² has dimension minus three.

The particle part of Q_chi is sum n_i/omega_i. In that limited particle regime
the source reduces to the B1 derived j=G²(phi−phi_star)sum n_i/omega_i.
Full mode-function Q_chi and p_Q also retain vacuum polarization and coherent
interference. A particle-occupation replacement does not reproduce the common
quantum source during a general nonadiabatic interval. The B1 field-space
production force j_prod is a numerical insertion choice, not an additional
force to add to exact quantum modes: production and its work are already
contained in the quantum Ward identity.

The matched action implies

    D_mu T^mu_nu,Q = −J_phi,Q D_nu phi_b,
    dot(rho_Q)+3H(rho_Q+p_Q) = J_phi,Q v = dot(x) Q_chi/2.

The fixed comoving cutoff K version of the reference subtraction satisfies the
same algebraic identity before discretization. A varying momentum boundary
would require its boundary-flux term. Keep the subtraction reference r>0
constant under metric/phi variation; it is distinct from the physical mass x
and from a physical EFT/gravity cutoff.

## 3. Finite local action: variations and counting once

For the retained homogeneous sector choose the declared local basis

    DeltaGamma_loc = int sqrt(−h)[−V(x)+F(x) R4+alpha R4²],
    R4 = 6[dot(H)+2H²],  alpha independent of x.

This is a finite change relative to the frozen reference action, not a second
copy of that checkpoint's subtractions or its already matched coefficients.
F may include a separately matched constant induced Einstein coefficient.
The displayed basis is a declared truncation, not the complete shell EFT:
additional independent scalar-derivative, higher-curvature or boundary
operators require their own matching and paired variations if retained.
With D_mu the shell covariant derivative, its tensor and scalar variations are

    T_munu,loc = −V h_munu
                −2[F G_munu+(h_munu box−D_mu D_nu)F]
                −2 alpha E_munu^(R²),
    E_munu^(R²) = 2 R R_munu − h_munu R²/2
                  +2(h_munu box−D_mu D_nu)R,
    J_phi,loc = x_phi[V_x−F_x R].

The resulting FRW sources, with dots denoting proper-time derivatives, are

    rho_loc = V−6F H²−6H dot(F)
              +alpha[36 dot(H)²−216H² dot(H)−72H ddot(H)],
    p_loc = −V+2F[2 dot(H)+3H²]+2 ddot(F)+4H dot(F)
            +alpha[108 dot(H)²+216H² dot(H)
                    +144H ddot(H)+24 H^(3)],
    Q_chi,loc = 2 V_x−2 F_x R,
    dot(F) = F_x dot(x),
    ddot(F) = F_xx dot(x)²+F_x ddot(x).

They satisfy dot(rho_loc)+3H(rho_loc+p_loc)=J_phi,loc v for
nonlinear V and F, including at H=0. Pressure is obtained by action variation,
without division by H. A constant alpha contributes no J_phi,loc. Making alpha
depend on x introduces an additional source −alpha_x R² and derivative terms
in its metric variation; that would be a separately declared extension.

Constant beta C² and gamma E4 belong to the covariant four-dimensional basis.
The first variation of C² vanishes on exactly conformally flat FRW, and E4 is
topological with appropriate fixed boundary data. This does not fix the beta
coefficient or justify an inhomogeneous extension without matching. Brane
initial/final boundaries and state surface terms require a specified
variational/state prescription; discarding them is not automatic for arbitrary
initial data.

The total RHS source in the explicit-source representation is

    rho = rho_Q,ref + rho_loc + rho_other,
    p = p_Q,ref + p_loc + p_other,
    J_phi = J_phi,Q,ref + J_phi,loc + J_phi,other.

An exactly equivalent representation absorbs V into
sigma_eff=sigma_R+kappa5² V and moves F and alpha to the metric left side:

    K_munu−K h_munu
       +kappa5²[F G_munu+(h_munu box−D_mu D_nu)F]
       +kappa5² alpha E_munu^(R²)
       = −sigma_eff h_munu/2 + kappa5² T_munu,Q+other/2,
    2 n(phi)+sigma_eff,phi−kappa5² F_phi R
       = −kappa5² J_phi,Q+other.

Remove rho_loc,p_loc from that right side and remove V_phi−F_phi R from its
scalar source when using these moved terms. Adding V to the tension while
retaining its rho and J_phi on the source side counts it twice. Moving the
metric F term does not remove its scalar F_phi R variation. Both a correct
action and a doubled local action can conserve stress; a Ward check alone
cannot diagnose consistent double counting. A pinned operator ledger must
state each term's coefficient and its sole location.

The smooth-FRW benchmark's matched one-loop convention is
V(r)=V_x(r)=V_xx(r)=F(r)=F_x(r)=alpha(r)=beta(r)=gamma(r)=0 in its
covariant derivative basis. Those are explicit finite EFT choices. They can
define a reference model for a future registration but are not measurements
of HDBLAST's tension, induced gravity or higher-curvature coefficients.

## 4. All three inherited sourced junctions

Set lambda=sigma_R/kappa5² and S_munu=−lambda h_munu+T_munu. Variation of
the two bulk/GHY pieces and the one shell action gives

    K_munu−K h_munu = kappa5² S_munu/2,
    scalar boundary coefficient
        −2 n(phi)/kappa5² −sigma_R,phi/kappa5² −J_phi = 0.

Trace reversing the metric equation and using
K^mu_nu=diag(k0,ks,ks,ks) gives the complete conditions

    ks = n(A) = [sigma_R+kappa5² rho]/6,
    k0 = n(B) = [sigma_R−kappa5²(2rho+3p)]/6,
    w = n(phi) = −[sigma_R,phi+kappa5² J_phi]/2.

These equations remain valid when the sources depend on intrinsic curvature;
they then form implicit higher-order boundary equations. They do not permit
using the same Neumann value for A and B once rho+p is nonzero. Their coordinate
form is A_Z=exp(B)ks, B_Z=exp(B)k0 and phi_Z=exp(B)w.

The B1 adapter returns quantities already multiplied by kappa5²:
`rho_code=kappa5² rho`, `pr_code=kappa5² p`, `J0_code=kappa5² J_phi`.
Its `Y phi_T` friction term must be absent for the proposed exact free-field
source. A separate measured/derived interaction may be added only with its
own action and equal opposite energy transfer.

For H0=1/rho_b, s=H0 tau and b=kappa5² H0³,

    rho=H0⁴ rho_hat, p=H0⁴ p_hat, Q_chi=H0² Q_chi_hat,
    J_phi=H0⁴ G_hat²(phi−phi_star) Q_chi_hat,
    kappa5² rho=b H0 rho_hat,
    kappa5² J_phi=b H0 J_phi_hat,
    H0/M5=b^(1/3),  M5=kappa5^(−2/3).

Changing a physical gravity scale, b, G_hat or a multiplicity is a model
choice; a formal small-source expansion does not supply an empirical prior.

## 5. Correct mode clock and source derivatives

The shell conformal time eta is distinct from the bulk chart time T:

    d eta = exp(B_b−A_b) dT,
    u_etaeta+[k²+a²x−a_etaeta/a]u=0,
    u u_eta*−u_eta u* = i.

In the inherited time coordinate this becomes

    u_TT+(A_T−B_T)u_T
       +[exp(2(B−A)) k²+exp(2B)x−A_TT−2A_T²+A_T B_T]u=0,
    u u_T*−u_T u* = i exp(B−A).

All displayed shell quantities are evaluated at Z=0. Using T as eta silently
changes both the mode equation and its normalization. Source integration uses
one fixed comoving momentum labeling and the same modes for rho,p,Q_chi.

A replacement quantum adapter must accept the geometry, scalar jets, state
and selected EFT treatment explicitly. The B1 `totals(m,Ab,phib,...)` particle
interface alone cannot supply H, curvature and their derivatives. A numerical
source iteration must update all three junctions together. Full differentiated
boundary data require rho_T,p_T,J_phi,T from that same treatment.

The finite local R² pressure and the reference fourth-order subtraction need
a fourth derivative of a. Differentiating that pressure boundary condition can
require a fifth derivative. F-dependent pressure and the current also mix
metric and scalar accelerations. Setting the additional finite alpha to zero
does not remove the high jets from the reference renormalized source.
Choosing an enlarged constrained initial-value system or justified order
reduction is therefore a physical/numerical prerequisite, not a solver detail.
In an order-reduced calculation use the chosen leading equations consistently
in rho,p,J_phi and the ledger, report the resulting remainder order, and compare
against the controlled first-order response. Do not silently discard pressure
terms or evaluate only the scalar source on old trajectories.

## 6. Paired energy, Codazzi and bulk/Weyl ledgers

The unchanged bulk scalar stress is

    T_AB,bulk = [partial_A phi partial_B phi
                −g_AB((partial phi)²/2+U)]/kappa5²,
    T_n tau,bulk = w v/kappa5².

With the inherited outward normals the exact local shell budget is

    dot(lambda+rho)+3H(rho+p)
       = [sigma_R,phi/kappa5²+J_phi]v
       = −2w v/kappa5² = −2 T_n tau,bulk.

Thus each interior's outgoing physical energy flux is −T_n tau,bulk.
The shell receives their sum. This is a signed local transfer law; it is not a
positive globally conserved scalar-plus-gravity energy on a spacetime without
a chosen time symmetry. The relevant comoving integral versions are

    Delta(a³ rho_Q)+int p_Q d(a³) = int a³ J_phi,Q v d tau,
    Delta[a³(lambda+rho)]+int(p−lambda) d(a³)
         = −2 int a³ T_n tau,bulk d tau.

Evaluate the same endpoints and time nodes in both ledgers. Retain bare mode
work, subtraction work and their paired renormalized difference separately;
algebraic subtraction cancellation in one diagnostic does not establish its
numerical accuracy in every source. A quantum work gain is not compensated
bulk loss until the simultaneous sourced scalar junction and bulk evolution
are evaluated. No force divided by v is needed at a turning point.

The shell Codazzi relation is

    dot(ks)+H(ks−k0) = −w v/3,
    dot(ks) = −w v/3−kappa5² H(rho+p)/2.

For the inherited flat-slice conformal PDE momentum constraint

    M=−3A_TZ−3A_T A_Z+3A_T B_Z+3A_Z B_T−phi_T phi_Z,
    M_b/exp(2B_b)=−kappa5²/2
        [dot(rho)+3H(rho+p)−J_phi v],

provided the same differentiated A boundary condition is used. This relation
is a local compatibility check, not a proof of interior constraint control.

Let Wcal=−E_tau tau/3 denote the inherited projected Weyl scalar. The bulk
projection gives

    H²=ks²+v²/12−w²/12+U/6+Wcal,
    dot(H)=ks(k0−ks)−v²/3−2Wcal,
    dot(Wcal)+4H Wcal
       =[4ks w v−4H v²−v dot(v)+w dot(w)−U_phi v]/6.

Use the complete J_phi in w and dot(w)=−[sigma_R,phiphi v+
kappa5² dot(J_phi)]/2. Wcal is a bulk geometry variable, not a population of
thermal particles; its sign is unrestricted by these identities. The Weyl
balance is a consistency identity requiring the bulk normal dynamics. It does
not close an autonomous four-dimensional radiation model.

## 7. Evidence and preserved limitations

`code/verify_action_junctions.py` varies the lapse-dependent minisuperspace
action before setting N=1, derives rho,p and S_x, and checks nonlinear V/F
and R² exchange. It independently trace reverses the shell stress, computes
the bulk flux and Gauss scalar projection, and checks Codazzi, the signed
budget, Weyl balance and mode-clock transformation. Its output records 21
positive identities and 14 detected negative controls. It uses explicit
exceptions under `python -O`. This checks displayed algebra, not quantum
state admissibility, numeric cutoff convergence or 5D existence/stability.

The read-only inherited-code audit found two narrow issues to retain:

- An earlier Chat13 physics-review sentence has a negative spatial K for the
  same outward normal while printing positive A_Z. The later doubled action,
  numerical review and executable solver fix the positive sign used here.
- B1 `hat_rates` explicitly excludes dot(j_prod); consequently the scalar
  `gFt` returned by `neumann_t` omits −exp(B)partial_T(kappa5² j_prod)/2
  during production. Its current `constraints()` uses only `gAt`, leaving
  `gFt` unused. This is a helper/documentation completeness limitation,
  not evidence that B1's current momentum constraint or evolution is wrong.

Earlier scientific outcomes remain unchanged. The original smooth-FRW curved
pressure cutoff matrix failed at K192; the separately registered K384 follow-up
passed all 13 groups under its fixed gates. B1's registered aggregate is
INCONCLUSIVE at lambda_c=1; below its stated gravity cutoff the derived
particle channel did not reproduce the phenomenological-friction outcome,
and strong-source runs lost reliability. None of these outcomes is superseded
by the present algebra checks. `PROSPECTIVE_WEAK_SOURCE_PROTOCOL.md` records
the missing inputs before a new physical experiment may be registered.
