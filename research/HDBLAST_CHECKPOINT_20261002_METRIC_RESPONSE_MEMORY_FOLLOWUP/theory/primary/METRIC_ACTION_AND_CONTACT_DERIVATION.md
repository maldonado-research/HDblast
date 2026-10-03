# Homogeneous metric response from the common minimal-scalar action

Prepared 2 October 2026. This is a prospective analytic derivation, with exact symbolic verification only. No physical metric profile, time history, momentum quadrature, mode solution, response value, or fitted gate was evaluated while preparing these files. The numerical experiment requires its own prior public registration. This is an established action/linear-response construction; no novelty or discovery claim is made.

The inherited public science head is `b807a4d549a40bc78b66122cb5719987091d2ce7`. The reference convention is the same positive physical reference mass `r`, order-two field-square subtraction, order-four stress subtraction, and common finite action as in the matched-stress and smooth-FRW checkpoints. Hold `r` fixed during every metric variation.

## 1. The physical metric perturbation and canonical action

For signature (-+++), start with

    S_chi = -1/2 ∫sqrt(-g)[g^(mu nu)∂_mu chi∂_nu chi + x chi²],
    ds² = a(eta)²[-deta²+dX²],
    u=a chi,
    S_u = 1/2 ∫[u′²-(∇u)²-(a²x-a″/a)u²],

where the integration-by-parts boundary term is fixed by the initial-state and no-boundary-variation convention. The canonical mode equation is

    u_k″ + [k²+a²x-a″/a]u_k = 0.

At the selected calibration, `a0=L=-1/eta`, `L′=L²`, `H=1`, `x=r=2`, and minimal coupling `xi=0`; the incoming regular Euclidean/BD modes become exact free conformal-time plane waves,

    u0,k=e^(-ik eta)/sqrt(2k).

Prescribe a homogeneous physical conformal perturbation, with the same incoming state,

    a=a0(1+epsilon h),    x=r=2 fixed,
    delta L=h′,          delta(a″/a)=2Lh′+h″,
    g=delta(a²x-a″/a)=4L²h-2Lh′-h″.

The first-order canonical mode variation `w=delta u` therefore obeys

    w″+k²w=-g u0,    w=w′=0 before the compact perturbation.

The exact source derivatives are

    g′=8L³h+2L²h′-2Lh″-h‴,
    g″=24L⁴h+12L³h′-2Lh‴-h⁗.

A shared canonical forcing does not make physical stress or field square identical to the mass-source channel. Their explicit metric operators and subtraction denominators change.

The exact curvature variation is

    R=6a″/a³,
    delta R=(6/a0²)[h″+2Lh′-4L²h],
    g=-a0² delta R/6.

For a de Sitter reference with general fixed `x` and `H`, the linear canonical potential is

    delta V=a0²delta x+2a0²(x-2H²)h-a0²delta R/6.

Consequently the simple free-plane-wave memory functional is special to `x=2H²`. It is not the retarded response at a shifted shell endpoint with a different `x/H²`.

## 2. Direct physical operators, before renormalization

Put `d0=u0′-Lu0` and `delta d=w′-Lw-h′u0`. The physical per-mode bare observables are

    Q_b=|u|²/a²,
    rho_b=[|u′-Lu|²+(k²+a²x)|u|²]/(2a⁴),
    p_b=[|u′-Lu|²-(k²/3+a²x)|u|²]/(2a⁴).

Their metric variations, normalized with the reference geometry, are

    a0² delta Q_b=2Re(u0* w)-2h|u0|²,
    a0⁴ delta rho_b=Re[d0* delta d]+(k²+2L²)Re[u0* w]
                    +2L²h|u0|²-4h a0⁴rho0,b,
    a0⁴ delta p_b=Re[d0* delta d]-(k²/3+2L²)Re[u0* w]
                    -2L²h|u0|²-4h a0⁴p0,b.

Complex conjugation is denoted by `*`. The reference bare scaled stresses are `k/2+3L²/(4k)` and `k/6-L²/(4k)`. The measure stays fixed at `k²dk/(2pi²)` because the cutoff is comoving. Physical-volume factors appear in the operators through `a^-2` and `a^-4`. Converting to physical momentum and forgetting its moving upper endpoint would change these expressions.

For the fixed quadratic mass law with `b=1`, the scalar current is `j=x_phi Q/2`, with `x_phi=2` at the unchanged scalar reference. Thus `a0²delta j=a0²delta Q`. No `Q0 h/4` scalar mass-law contact is present: that prior contact arose from varying `phi`, which is held fixed here.

## 3. Full fixed-reference adiabatic inventory

Let `omega=sqrt(k²+a²r)`, `D=a²(x-r)=0`, `T=a″/a`, and retain the full inherited inventory:

    U2=-T/(2omega)-omega″/(4omega²)+3omega′²/(8omega³),
    U4=-U2²/(2omega)-U2″/(4omega²)+omega″U2/(4omega³)
        +3omega′U2′/(4omega³)-3omega′²U2/(4omega⁴),
    b=-k²/3-a²r,  c=1-b/omega²,
    J2=Lomega′/omega²+omega′²/(4omega³),
    J4=LU2′/omega²-2Lomega′U2/omega³
        +omega′U2′/(2omega³)-3omega′²U2/(4omega⁴).

Define scaled subtraction numerators

    F_Q=1/(2omega)-U2/(2omega²),
    N_rho=omega/2+[L²/omega+J2]/4
          +[U2²/omega-L²U2/omega²+J4]/4,
    N_p=k²/(6omega)+[cU2+L²/omega+J2]/4
        +[cU4+bU2²/omega³-L²U2/omega²+J4]/4.

Physical subtractions are `F_Q/a²`, `N_rho/a⁴`, and `N_p/a⁴`. The metric direction is

    delta omega=2L²h/omega,
    delta L=h′,
    delta T=2Lh′+h″,
    delta b=-4L²h,
    delta D=0.

Differentiate every term, including denominators and both `U2` and `U4`, at fixed eta, comoving k, physical x, and physical r. Derivatives commute with this variation. The physical scaled variations are `delta F_Q-2h F_Q`, `delta N_rho-4h N_rho`, and `delta N_p-4h N_p`. A truncated WKB frequency is only subtraction data; it never defines the physical state.

`metric_contact_algebra.py` independently implements this inventory with exact first-order dual algebra and proves the retained Ward identities, finite-moment formulas, baseline primitives, metric-minus-mass contact primitives, continuum limits, and detection of eight symbolic omissions. `FORMULAS.json` is machine-readable export. This is a different algebraic implementation from a numeric dual-WKB producer.

## 4. Exact reduction to the mass-memory route plus local contacts

Use the previous scalar-memory response with canonical forcing `g` to compute `q_mass`, `R_mass`, and `P_mass`. Here `q_mass=a0²delta Q_mass`, `R_mass=a0⁴delta rho_mass`, and `P_mass=a0⁴delta p_mass`; the physical mass-source comparison is `delta x=g/a0²`.

Let `Q0K` be the physical renormalized reference field square and `R0K=a0⁴rho0K`, `P0K=a0⁴p0K` the scaled reference stresses at finite comoving K. Then the direct operator reduction gives

    q_metric=q_mass-2h L²Q0K-Cq,
    R_metric=R_mass-4h R0K+Er,
    P_metric=P_mass-4h P0K+Ep.

These formulas derive both stresses from their physical operator and full subtraction; neither stress is defined by a Ward or trace reconstruction. The local terms are exact finite-band primitives:

    Cq=∫dmu[delta F_Q,metric+g/(4omega³)],
    Er=∫dmu[(Lh′+h″/4)/k-(delta N_rho,metric-delta N_rho,mass)],
    Ep=∫dmu[-h″/(4k)-(delta N_p,metric-delta N_p,mass)].

All high-order terms reduce to

    I_n=∫_0^K k²dk/omega^n
       =(2L²)^((3-n)/2) Σ_j (-1)^j binom((n-5)/2,j) v^(2j+3)/(2j+3),
    v=K/sqrt(K²+2L²),   n=5,7,9,11,13,15.

The common measure factor `1/(2pi²)` is absent from `I_n` and all expressions exported with suffix `times_2pi2`. For the lower powers, all logarithms cancel in the combined primitives. Specifically,

    ∫k²[1/k-1/omega-L²/omega³]dk=L²v²/(1+v).

The exact baseline lower-order primitives are

    Rlead=L⁴v²(2v+3)/[4(1+v)²],
    Plead=-L⁴v²(4v+3)/[12(1+v)²].

Subtract the high-power baseline coefficients in `FORMULAS.json` times `I_n`, then divide by `2pi²`. This gives the full baseline without a cancellation-prone numerical integral. The exact physical field-square baseline is

    Q0K=[v²/(2(1+v))-v³/48-v⁵/16]/(2pi²).

The machine-readable density and pressure coefficient tables extend through `omega^-13` and `omega^-15`; no term is discarded. Each closed contact and baseline has a verified derivative equal to its original combined integrand and a verified zero lower endpoint. Those two properties fix each primitive uniquely.

For the continuum, the limits are

    Q0=1/(12pi²),
    R0=11L⁴/(960pi²),    P0=-R0,
    2pi²Cq=-L²h/2+Lh′/12+h″/24,
    2pi²Er=L[-21L²h′+7Lh″+3h‴]/120,
    2pi²Ep=-7L³h′/360+7L²h″/180+Lh‴/120-h⁗/120.

Evaluating these continuum expressions at finite K would lose precisely the cutoff artifacts that a finite-band calibration must check. The exact rational v formulas can instead be enclosed by directed interval evaluation on the actual K-to-infinity interval, combined with the separately derived scalar-memory omitted-band bound for forcing `g`. These are omitted-UV-band bounds, not certificates for numerical finite integrals or total floating-point error.

## 5. Independent conservation, trace, and the absence of a physical Ward drift

For exact modes at fixed physical x, the bare per-mode identity is

    rho_b′+3L(rho_b+p_b)=0.

The complete retained positive-reference subtraction grades separately obey the same identity for this fixed-mass metric direction, before integrating at any fixed comoving K. This was checked directly; no continuum-covariance assumption was used to remove a cutoff term.

Therefore the linear physical renormalized metric Ward equation is

    delta rhoK′+3L(delta rhoK+delta pK)+3h′(rho0K+p0K)=0,
    R_metric′-L R_metric+3L P_metric+3h′(R0K+P0K)=0.

The last term is generally nonzero at finite K. It vanishes in the removed-cutoff reference limit because `p0=-rho0`. The distinct fixed-reference drift in the second-order canonical-energy work functional is not a drift in this physical fixed-mass stress Ward equation.

A future numerical conservation check must independently integrate the direct stresses and the full finite-K baseline term over the saved history. It must not define pressure by solving this equation. A finite-K trace test must likewise use the differentiated finite-K anomaly, not substitute the continuum heat-kernel anomaly. Both diagnostics are supplementary to independent direct stress comparison.

## 6. Gauge and shifted-root boundaries

The declared perturbation changes both homogeneous lapse and spatial scale in the same conformal factor. It is not an arbitrary lapse response. Under an infinitesimal homogeneous time change `eta -> eta+epsilon T`, the temporal conformal perturbation is `LT+T′`, while the spatial one is `LT`. Their equality requires `T′=0`. A constant translation produces `h=LT`; it is not a nonzero compact perturbation with unchanged past. Hence generic compact `h` here is a physical metric deformation. A complete compact diffeomorphism test requires a separate lapse degree of freedom and its response.

Varying the underlying covariant action before imposing the conformal gauge defines the physical density operator, but this single response direction does not supply all independent metric kernels or the gravitational lapse constraint. It also does not supply inhomogeneous/vector/tensor response.

The actual shell stationary root has its own `x`, curvature, scalar-current derivatives, boundary data, and quantum state. If `x/H²` shifts away from two, its canonical modes and memory kernel change. Neither this analytic reduction nor a passing prospective special-point metric experiment can be called a shifted-root Hessian, coupled shell evolution, self-consistent Einstein solution, quantum stability proof, or heating/particle-yield calculation.

## 7. Development evidence and a bounded next experiment

Two symbolic development issues are preserved. V1 stalled while proving a rational primitive after expanding it in a needlessly complicated radical k-coordinate; it was explicitly interrupted, with original source and KeyboardInterrupt evidence retained. V2 caught an omitted `-2Lh‴` in a hand-written expected `g″` expression. The generated inventory always differentiated `g` correctly; the failed expected identity and its source are retained. The passing version proves the same primitive derivatives in the exact rational v-coordinate and corrects the expected derivative. Neither issue involved a physical source or experiment.

An actionable next registration can prescribe two compact smooth conformal profiles, fixed incoming BD data, finite K values, six observations, and two independent numerical routes. Route one uses scalar forcing memory plus these independently derived exact local contacts. Route two solves complex forced modes and directly varies the full physical operators and unsimplified W0/W2/W4 inventory. Include pure algebra, pre-source zeros, real/complex finite guards, source pins, zero initial forced modes, Wronskian diagnostics, local-contact/operator mutations, independently integrated Ward history with finite-K baselines, refinement checks, and properly scoped omitted-band comparisons. Freeze every source, setting, gate, and implementation before physical evaluation. The result would calibrate a genuinely new metric response channel while retaining explicit separation from the missing shifted-root and full metric/lapse channels.
