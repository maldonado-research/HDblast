# A bounded static de Sitter bridge for HDBLAST quantum sources

Prepared 2 October 2026 UTC. **Feasibility and conditional analytic response only.** The proposed next study below has not been registered or executed. No de Sitter renormalized expectation value, physical parameter choice, new mode evolution, new radial integration, or coupled quantum shell solution is reported here.

The inherited corrected +1 branch offers a simpler first quantum-source closure than time-dependent bulk evolution: choose a massive, minimally coupled shell scalar in its Euclidean de Sitter vacuum, derive density and scalar current from the same matched action, and feed their leading values into the two static junctions. This would test a vacuum-polarization response, rather than particle production or reheating. The simplification follows from de Sitter invariance; it does not turn the prescribed smooth-FRW benchmark into a solution of the shell model.

## 1. Pinned inputs and inherited evidence

The checkout was clean at `e17a01b428bb8049e919c42376ab0359e41d613c`. The primary action and corrected profiles are in

`hdblast/checkpoints/2026-09-22_scalar-profile-branch-and-matter_zenodo-22922928/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/`.

The archived response Jacobian and linear-stability evidence are in

`research/HDBLAST_CHECKPOINT_20260927/stability_gauge_invariant/`.

The local-action conventions are read from `research/HDBLAST_CHECKPOINT_20261002_SMOOTH_FRW/COMMON_ACTION_SMOOTH_FRW.md`. `SOURCE_PINS.json` freezes SHA-256 hashes and sizes of the 12 public action, profile, solver, Jacobian producer/output, stability review, and matching inputs at the source commit. The verifier compares actual files against those expected hashes before calculation. `STATIC_BRIDGE_CHECKS.json` records the verified source pin separately from any observed execution checkout HEAD. The new algebra script imports neither solver and contains no integration or optimization routine.

At registered detuning delta=0.001 and c=0.5975949350280132, the archived corrected branch has phi_b approximately 0.9999159473169134, warp radius approximately 129.92476284968, and H_b^2 approximately 5.924014794329e-5 in registered units. The 27 September checkpoint supplies independently verified floating-point evidence of linear scalar and tensor stability and a nonsingular static junction Jacobian. Those results concern the classical empty branch. They do not establish quantum stability, nonlinear attraction, or the response after adding a shell quantum sector.

## 2. The actual action and static junctions

Use signature (-++++), two identical interiors, one shell action, and the outward normal from each interior toward the shell. Write K= kappa_5^2 in the formulas below; K is a coupling, not a momentum cutoff or extrinsic curvature.

    S = sum_two_copies { (1/(2K)) int sqrt(-g)[R-(grad phi)^2-2U]
                         +(1/K) int_shell sqrt(-h) K_extrinsic }
        +int_shell sqrt(-h)[-sigma(phi)/K+L_chi].
    W(phi)=1-phi+phi^3/3,
    U=W_phi^2/2-2W^2/3,
    sigma=2W+delta(1+c phi).
    L_chi=-(D chi)^2/2-x(phi)chi^2/2,
    x(phi)=m0^2+gbar^2(phi-phi_star)^2,
    x_phi=2gbar^2(phi-phi_star).

Here phi is the registered geometric bulk scalar and chi is a canonical shell scalar. The source definitions are T_mn=-2 delta Gamma/(sqrt(-h) delta h^mn), S_x=-delta Gamma/(sqrt(-h) delta x), Q=2S_x, and

    J_phi = x_phi Q/2 = gbar^2(phi_b-phi_star) Q.

The gravitational normalization, m0, gbar, phi_star and finite matching are not fixed by the empty branch or by the illustrative smooth-FRW control parameters.

For a static de Sitter-sliced bulk write ds^2=dy^2+R(y)^2 gamma_dS, with unit de Sitter gamma_dS, shell y=y_b, H_b=1/R_b, q=R'/R, and w=phi'. On the chosen copy n=+partial_y. The radial equations are

    R''=-R(w^2/4+U/6),  phi''=U_phi-4qw,
    R'^2=1+R^2(w^2/12-U/6).

The regular cone has R~y and phi'=0. A de Sitter invariant state at constant phi_b has

    T_mn=-rho h_mn,  p=-rho,  rho and Q constant on the shell.

Consequently both gravitational junctions coincide and the complete static residuals are

    E_1(u)=q_b-sigma_b/6-K rho(H_b,phi_b)/6=0,
    E_2(u)=w_b+sigma_phi,b/2+K J_phi(H_b,phi_b)/2=0.

Both sources must be retained. The radial constraint then gives the exact static identity

    H_b^2 = (sigma_b+K rho)^2/36
            -(sigma_phi,b+K J_phi)^2/48+U_b/6.

That algebraic identity alone does not select the regular radial solution. The regular cone and the two residuals still have to agree. At phi=1, sigma_phi=delta c; the identically constant scalar fails the vacuum scalar junction. The corrected scalar profile must be the expansion point.

## 3. Conditional linear response from the archived Jacobian

Use u=(ell,y_b), ell=log10(abs(eta_h)), eta_h=phi_h-1<0. Define E_vac=(q-sigma/6,w+sigma_phi/2). The archive's `spectrum.py` formed raw derivatives of these residuals by symmetric finite differences of size 1e-6 in each coordinate, using its `gi_core.Background` integrations at the archived shooting values. This is the unscaled residual Jacobian: do not confuse it with a residual divided by delta in the earlier branch solver.

At delta=0.001 its matrix is

    J_vac = [[-6.938893903907228e-12, -5.924014057079319e-5],
             [-6.879322173491877e-4, -4.6460663028142276e-4]].

Its floating-point condition number is 16.9364181840799 and determinant is approximately -4.075e-8. Nonsingularity supports conditional local solvability if the exact residual map is smooth and its exact Jacobian is nonsingular. Finite differences do not certify those hypotheses or bound nonlinear errors. In particular the tiny first-column metric entry corresponds to a residual difference of approximately 1.39e-17 before division by 2h. It lies at the numerical floor. The inverse-matrix component dy_b/d(K J_phi), approximately -8.51e-5, is entirely driven by that entry and is explicitly **unresolved and numerical-floor dominated**, rather than a precise physical susceptibility. The condition number does not validate each near-zero matrix component.

For independently supplied small source values at the unperturbed shell,

    Delta u = K J_vac^{-1} (rho/6, -J_phi/2)^T + O(source^2).

`STATIC_BRIDGE_CHECKS.json` supplies the inverse response matrix and verifies its algebra against the archived matrix. This is a reanalysis of archived values, not a new background solution. To obtain accurately resolved changes of phi_b and H_b one also needs the endpoint output derivatives:

    Delta phi_b = phi_ell Delta ell + w_b Delta y_b,
    Delta H_b^2 = (H_b^2)_ell Delta ell - 2H_b^2 q_b Delta y_b.

The archive does not provide phi_ell and (H_b^2)_ell with an independently resolved error estimate. They must be derived and prospectively computed through variational equations or independently resolved differences. The moving-endpoint columns themselves have exact identities on a vacuum solution:

    (E_1)_y = -H_b^2,
    (E_2)_y = B_b w_b,  B_b=phi''_b/phi'_b+sigma_phiphi,b/2.

The archived differences agree with these endpoint identities at the declared absolute checks in the script.

For a genuinely self-consistent source, the Jacobian of the coupled static map is

    M_1a=(J_vac)_1a-K partial_a rho/6,
    M_2a=(J_vac)_2a+K partial_a J_phi/2.

These source derivatives include changes in H_b and phi_b. A formal loop-amplitude expansion Gamma_chi=epsilon Gamma_1 gives Delta u=epsilon K J_vac^{-1}(rho_1/6,-J_1/2) at first order; source feedback in M enters at the next order. This is a useful first response study but does not close a finite-amplitude coupled BVP. Small K rho relative to sigma alone is insufficient: H_b^2 is near a detuning cancellation. The protocol must also bound the predicted relative H_b^2 change and the scalar-profile response.

## 4. A leading small-detuning check and its limited interpretation

For this paragraph only, count delta, s=K rho and t=K J_phi as formal quantities of the same small order epsilon. Hold rho and J_phi as supplied shell values when writing the expansion; this is not a prescription for varying a quantum action. Put

    d=delta c,  nu=d+t,  tau=delta(1+c)+s.

The regular linear scalar profile on AdS+ has logarithmic slope 14/9 at large radius. The leading scalar match is (14/9+2)eta_b+nu/2=0, hence

    eta_b = -9nu/64+O(epsilon^2).
    H_b^2 = tau/27+tau^2/36-d nu/192+nu^2/384+O(epsilon^3)
           = tau/27+tau^2/36+(t^2-d^2)/384+O(epsilon^3).

The source-free result is the inherited correction -delta^2 c^2/384. The absence of a d t term in this **frozen-shell-value** expression is an algebraic cancellation; it does not mean the scalar source is physically negligible or that general quantum backreaction leaves H unchanged.

An instructive action-consistent example has only a local potential V(phi)=V0+J0(phi-1), with K V0=s0 and K J0=t. The density on the shifted shell is V0+J0 eta_b, so it must vary with eta_b. Setting tau0=delta(1+c)+s0 gives instead

    H_b^2 = tau0/27+tau0^2/36-(d+t)^2/384+O(epsilon^3).

Here the potential can be absorbed as sigma_eff=sigma+K V, provided both V and its current are removed from the source side. This example illustrates why frozen density and an action-paired potential are different problems. It is not the general curved-space quantum source: curvature terms and state-dependent contributions cannot in general be absorbed as a phi-only tension.

## 5. Exact de Sitter limits of the common local action

For Gamma_loc=int sqrt(-h)[-V(x)+F(x)R+alpha R^2+beta C^2+gamma E4], with alpha,beta,gamma constant, exact de Sitter and constant phi_b give

    rho_loc=V-6F H_b^2,  p_loc=-rho_loc,
    Q_loc=2V_x-24H_b^2 F_x,
    J_loc=V_phi-12H_b^2 F_phi=x_phi Q_loc/2.

Thus J_loc is generally **not** partial_phi rho_loc: the latter is V_phi-6H_b^2 F_phi at fixed H_b. A constant F still changes the branch despite making no scalar current. Constant alpha R^2 has zero metric first variation on four-dimensional exact de Sitter; beta C^2 vanishes on this conformally flat background and constant Euler is topological under the declared boundary convention. Those facts do not eliminate the curved-space trace anomaly or determine the quantum stress. Phi-dependent curvature-squared coefficients would also carry a scalar source and require a different retained EFT action.

Every retained local piece must appear once, either inside the matched quantum sources or explicitly in the shell action/junctions. If F is moved to the gravitational side, its scalar variation must be moved consistently as well. Subtracting it twice changes the static branch even though the static Ward identity remains 0=0.

## 6. An explicit state and calculational route for the next study

The candidate state is the Euclidean (Bunch-Davies) de Sitter invariant Hadamard vacuum for a free massive minimal scalar with x(phi_b)>0 and H_b^2>0. It is defined by the regular Euclidean Green function on S4 continued to Lorentzian de Sitter. The strictly positive gap is essential: the minimally coupled massless field has a zero-mode obstruction to this prescription. A static-patch temperature H/(2pi) does not by itself define a globally matched state or a radiation population. The state here is different from the exact static incoming plane wave of the prescribed FRW control; that control's state-specific decoupling argument cannot simply be cited as a de Sitter calculation.

For S4 radius 1/H, the scalar eigenvalues and degeneracies are

    lambda_l=H^2 l(l+3)+x,
    d_l=(l+1)(l+2)(2l+3)/6,  l=0,1,...,
    Vol(S4)=8pi^2/(3H^4).

The formal Euclidean determinant Gamma_E=1/2 sum_l d_l ln(lambda_l/mu^2) and coincident Q=(1/Vol)sum_l d_l/lambda_l both require covariant renormalization. One prospective implementation uses the inherited PV weights (1,-3,3,-1), regulator masses x+j Lambda^2, fixed positive reference r+j Lambda^2, and the same covariant local action subtraction. Its finite-Lambda combined log and Green-function sums have large-eigenvalue tails of order lambda^-3 and lambda^-4 respectively after the three PV moment cancellations. With d_l~l^3 and lambda~l^2, both combined sums converge at fixed Lambda. This does not justify interchanging the harmonic cutoff and Lambda limits or establish their numerical error bounds.

The Euclidean image of Gamma_sub from the inherited Lorentzian convention is

    Gamma_sub,E=Vol[C_V-12H^2 C_F-144H^4 C_alpha-24H^4 C_gamma],

with beta's term zero on S4. The C coefficients must be those of one fixed finite matching prescription and the reference r must be held fixed in all variations. At fixed Lambda, remove the harmonic cutoff using a registered tail bound; then study Lambda removal with Euclidean-state heavy-field asymptotics or an independently controlled limit. Matching about x=r to V(r)=V_x(r)=V_xx(r)=F(r)=F_x(r)=alpha(r)=beta(r)=gamma(r)=0 is the **control's finite convention**, not an established physical HDBLAST gravitational matching condition.

For a common renormalized Euclidean action Gamma_E=Vol W_E(H,x), its constant sources satisfy

    Q=2 partial_x W_E,
    rho=W_E-(H/4)partial_H W_E,
    J_phi=x_phi Q/2,
    partial_x rho=Q/2-(H/8)partial_H Q.

The last relation is a useful nonvacuous static pairing check. It catches a dropped F_phi current even though energy exchange at constant phi cannot. The metric derivative is at fixed x, r, mu and regulator masses; the scalar derivative is at fixed metric and reference r.

The inherited positive-reference continuum trace identity, **if the same prescription has been established for this Euclidean-state calculation**, specializes to

    -4rho=-x Q+A_r,
    A_r=[(x-r)^2/2-2(x-r)H^2+(29/15)H^4]/(16pi^2).

It should be checked against independently obtained metric variation and Q. It must not be imposed to manufacture density or confused with a finite-harmonic-cutoff identity. No value of Q or rho is evaluated by the present algebra checks.

## 7. Concrete prospective study and gates

The next bounded study should be **one-loop Euclidean-state source and linear response on the inherited corrected delta=0.001 branch**. Its intended deliverables are resolved rho, p=-rho and Q in one stated finite convention; J_phi from the declared mass law; independently resolved endpoint sensitivities; and the conditional first-order Delta phi_b and Delta H_b^2 from the archived regular-cone branch. It should end before dynamic bulk evolution or claims about heating.

Before physical sums or response integrations, a prospective registration must freeze:

1. The physical gravitational normalization, mass law parameters and positive mass gap, EFT cutoff/domain, and their origin. None is selected in this note. If only a dimensionless mathematical benchmark is intended, label it separately and do not present it as the physical shell.
2. Euclidean state, retained local shell operators, physical matching conditions or an explicitly adopted benchmark convention, fixed reference r, and exact placement of every counterterm. State whether the proposed result is a formal one-loop response or a finite-amplitude static closure.
3. Harmonic truncation/refinement sequence, tail treatment, arithmetic precision, PV limit order and regulator sequence, independent derivative method, and absolute plus relative acceptance tolerances. The characteristic vacuum cancellation and source scales must determine tolerances before results are known.
4. Variational or independent-difference sensitivity computation at the archived branch, with cone-start and integration-error checks. Include the exact endpoint Jacobian identities and a reproducible scale for the predicted Delta phi and Delta H^2.
5. Gates for de Sitter symmetry, common-action cross derivatives, independently computed metric stress, matched trace identity, subtraction independence under explicitly paired finite shifts, and controlled truncation/regulator changes. Analytic negative controls should drop F_phi, flip the scalar-current sign, freeze the density of a linear potential, and double-count a local term.
6. A small-response gate based on the actual conditional changes and matched EFT hierarchy. Failure ends in a bounded diagnostic report or a separately registered finite-amplitude BVP; it does not trigger unregistered mass/coupling tuning or a switch to a desired cosmological history.

This is a concrete proposal ready for the missing inputs to be fixed, not an active experimental registration. Even a successful static study would establish a vacuum-polarization bridge at one endpoint, with perturbative error limits. It would not connect the unstable original shell to this endpoint, establish nonlinear dynamics, count produced particles, transfer energy into radiation, or predict observations.

## 8. Checks executed

Run `python verify_static_bridge.py --repo /path/to/HDblast --output /path/to/replay-output`; the repository can also be detected from current/script ancestors. The source pin is checked even when the execution checkout has a later HEAD. `--output` keeps the archived evidence untouched.

The current run passes **33 algebra/archive assertions and detects five deliberate wrong-formula controls**. It verifies the action coefficients at the + vacuum, both source junctions and constraint, exact endpoint Jacobian identities, local F/R2/Euler static limits, static common-action integrability, PV moments/tail cancellations, the formal small-source expansion and potential-paired correction, and inversion/endpoint consistency of the archived Jacobian. It contains no chi mode solve or new bulk/radial solve. Reproducing algebra and archived matrix inversion is the extent of the completed result.
