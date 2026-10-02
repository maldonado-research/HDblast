# A common finite-action convention for the smooth FRW control

Draft analytic work for prospective review, prepared 2026-10-01 America/Los_Angeles (checkpoint folder date 20261002 UTC). This document contains no new field evolution, cutoff result, fitted cosmology, discovery claim, or claim about the archived shell. Its starting prescriptions are the flat PV and FRW checkpoints in the public scientific stack pinned at `7e508d946e9ecdc542ca3059263d93132321838b`; the later merge `0949317d5e0310535caa50576a83d4487815dd7a` contains that stack. Flat conventions below are those of `HDBLAST_CHECKPOINT_20261002_PV/THEORY_AND_MATCHING.md` and `DIRECT_CUTOFF_APPENDIX.md` at the pinned commit.

## 1. Conventions and the selected finite scheme

Use signature (-+++), a minimally coupled real scalar with prescribed mass squared x, spatially flat FRW, cosmic time t, conformal time eta, and

    ds² = -dt² + a(t)² dX² = a(eta)²[-deta²+dX²],
    H = adot/a,  R = 6(Hdot+2H²),  dt = a deta,
    T_munu = -2 δGamma/(sqrt(-g) δg^munu),
    S = -δGamma/(sqrt(-g) δx),  Q = 2S.

A common local action has

    Gamma_loc = ∫sqrt(-g)[-V(x)+F(x)R+alpha R²+beta C²+gamma E4].

The constants alpha,beta,gamma are independent of x. Its homogeneous variations are

    rho_loc = V - 6F H² - 6H Fdot + alpha(36Hdot²-216H²Hdot-72H Hddot),
    p_loc = -V + 2F(2Hdot+3H²)+2Fddot+4H Fdot
            + alpha(108Hdot²+216H²Hdot+144H Hddot+24H'''),
    S_loc = V_x - F_x R,   Q_loc = 2V_x - 2F_x R.

Here H''' is the third cosmic-time derivative of H. Constant beta C² has zero first metric variation on conformally flat FRW; constant E4 is topological in four dimensions after specifying boundary terms. These properties do not determine beta away from FRW. For arbitrary smooth F(x), direct variation gives

    rhodot_loc+3H(rho_loc+p_loc) = xdot S_loc.

The selected benchmark convention sets the one-loop local potential and curvature response, in the covariant heat-kernel/PV derivative basis, to

    V(r)=V_x(r)=V_xx(r)=0,
    F(r)=F_x(r)=0,
    alpha(r)=beta(r)=0,

at a fixed positive reference r. This means local derivative-expansion coefficients about x=r, not the full state-dependent effective action, vanish to these orders. We also choose gamma(r)=0 and discard total derivative action terms with the standard no-boundary-variation convention. These are declared finite EFT matching conditions, not measured HDBLAST couplings. The beta and gamma choices specify a covariant extension invisible to spatially flat FRW. Positive-reference fourth-order adiabatic stress and second-order Q need no additional local terms in this convention; sections 3–5 establish the FRW coefficient bridge and state its decoupling assumptions.

The reference r stays constant when varying the metric or x. Letting r track x, the scale factor, a cutoff, or a trajectory is a different prescription and invalidates these formulas.

## 2. Complete computational observables

The exact canonical modes obey

    u_k''+[k²+a²x-a''/a]u_k=0,
    u_k u_k'^* - u_k' u_k^* = i.

Let L=a'/a, D=a²h, h=x-r, w=sqrt(k²+a²r), and define

    u2 = (D-a''/a)/(2w)-w''/(4w²)+3w'²/(8w³),
    u4 = -u2²/(2w)-u2''/(4w²)+w''u2/(4w³)
         +3w'u2'/(4w³)-3w'²u2/(4w⁴),
    b = -k²/3-a²r,   c = 1-b/w²,
    J2 = Lw'/w²+w'²/(4w³),
    J4 = Lu2'/w²-2Lw'u2/w³+w'u2'/(2w³)-3w'²u2/(4w⁴).

The subtraction integrands are

    rho0 = w/(2a⁴),
    rho2 = [(D+L²)/w+J2]/(4a⁴),
    rho4 = [u2²/w-(D+L²)u2/w²+J4]/(4a⁴),
    p0 = k²/(6a⁴w),
    p2 = [c u2+(L²-D)/w+J2]/(4a⁴),
    p4 = [c u4+b u2²/w³-(L²-D)u2/w²+J4]/(4a⁴),
    Q0 = 1/(2a²w),   Q2 = -u2/(2a²w²).

Bare integrands, with v_k=u_k'-L u_k, are

    rho_b = [|v_k|²+(k²+a²x)|u_k|²]/(2a⁴),
    p_b = [|v_k|²-(k²/3+a²x)|u_k|²]/(2a⁴),
    Q_b = |u_k|²/a².

Integrate rho_b-rho0-rho2-rho4, p_b-p0-p2-p4 and Q_b-Q0-Q2 against k²dk/(2pi²). Set S=Q/2. Use exact physical modes at every k. A truncated WKB frequency is solely subtraction data and need not be positive at low k. All denominators are positive for a>0 and r>0, including x=0.

The algebra in `derive_local_trace.py` verifies the three graded Ward identities identically for general symbolic FRW jets. This strengthens, without changing, the previous finite selection of exact rational-jet checks:

    rho0' + 3L(rho0+p0) = 0,
    rho2' + 3L(rho2+p2) = x' Q0/2,
    rho4' + 3L(rho4+p4) = x' Q2/2.

They imply an exact exchange identity at a fixed comoving cutoff K, before a numerical approximation:

    rhoK' + 3L(rhoK+pK) = x' QK/2.

A time-dependent physical cutoff would add moving-boundary terms. No such cutoff is used in this identity. For a=1 the expressions reduce exactly to the inherited direct flat prescription, including its x'' pressure subtraction.

## 3. A constructive PV action bridge

Take PV weights c_j=(1,-3,3,-1), j=0,1,2,3, physical masses y_j=x+jLambda², and positive reference masses r_j=r+jLambda². They satisfy

    Σc_j = Σc_j r_j = Σc_j r_j² = 0.

Each regulator is assigned the same geometrical background and its exact static incoming vacuum. The regulator fields are auxiliary subtractions, not matter species. Lambda is a regulator mass, not a physical gravity cutoff.

Define

    V_PV(x) = Σc_j y_j² ln(y_j/mu²)/(64pi²),
    F_PV(x) = Σc_j y_j[ln(y_j/mu²)-1]/(192pi²),
    C_V(x) = V_PV(r)+V_PV'(r)h+V_PV''(r)h²/2,
    C_F(x) = F_PV(r)+F_PV'(r)h,
    L_PV = Σc_j ln(r_j/mu²),
    C_alpha = -L_PV/(2304pi²),
    C_beta = -L_PV/(3840pi²),
    C_gamma = +L_PV/(11520pi²).

The scale mu cancels from every PV combination. Subtract from the PV one-loop action the local action

    Gamma_sub = ∫sqrt(-g)[-C_V+C_F R+C_alpha R²+C_beta C²+C_gamma E4].

Thus the renormalized action is Gamma_PV-Gamma_sub. Stress and current must both be obtained with this same sign and subtraction. The C_V and C_F terms are exactly the inherited flat matching convention. The additional constants implement the declared curvature-squared matching.

These curvature coefficients follow from the local heat-kernel coefficient for a minimal scalar, with h counted as order two:

    a2(r) = (h-R/6)²/2 + (Riem²-Ric²)/180 + box R/30 - box x/6,
    (Riem²-Ric²)/180 = C²/120-E4/360.

The Lorentzian local effective-action coefficient multiplying the geometric a2 logarithm is -ln(y/mu²)/(32pi²). Hence the R², C², E4 coefficients above. The total derivative terms are needed for local trace identities even though their constant-coefficient integrals can be discarded under the stated boundary convention. The C² coefficient is the declared covariant heat-kernel extension, not something inferred from conformally flat dynamics.

### Direct FRW coefficient calculation

This bridge need not rest on Ward consistency or a flat limit. At a fixed time use physical momentum p=k/a and Omega=sqrt(p²+r_j). All per-mode subtraction densities carry an overall a^-3 that cancels the measure transformation. With Hdot=d, Hddot=e, xdot=z,

    u2/a = (h-d-2H²)/(2Omega)
           -r_j(d+3H²)/(4Omega³)+5r_j²H²/(8Omega⁵).

The conformal derivative is w'=a²r_jH/Omega. In particular, the coefficient is d+3H² in w''/a³, not d+4H². Define I_n=∫_0^P p² dp/[2pi²(p²+r_j)^(n/2)]. As P tends to infinity,

    I_1 = [P²+r_j/2-r_j ln(2P)+r_j ln(r_j)/2]/(4pi²)+o(1),
    I_3 = [ln(2P)-ln(r_j)/2-1]/(2pi²)+o(1),
    I_n = sqrt(pi)Gamma((n-3)/2) r_j^((3-n)/2)
          /[8pi² Gamma(n/2)], n>3.

The summed order-two geometric energy is

    H² Σc_j[I_1/4+r_j I_3/4+r_j² I_5/16]
      = -6F_PV(r)H².

This fixes the previously undetermined constant F(r) matching relative to PV; it is not inferred from flat space. At order four, all convergent, nonlogarithmic coefficients are independent of r_j after integration and cancel by Σc_j=0. The surviving logarithmic coefficients are

    ∫Σc_j rho4 = L_PV/pi² [h²/64-hH²/32-Hz/32
                                 +(-d²+6H²d+2He)/64],
    ∫Σc_j p4 = L_PV/pi² [-h²/64+hd/48+hH²/32+xddot/96+Hz/48
                               -3d²/64-3H²d/32-He/16-H'''/96],
    ∫Σc_j Q2 = L_PV(h-d-2H²)/(16pi²).

Together with orders zero and two these are precisely

    ∫Σc_j rho_sub,j = rho[Gamma_sub],
    ∫Σc_j p_sub,j = p[Gamma_sub],
    ∫Σc_j Q_sub,j = 2C_V'(x)-2C_F'(x)R.

`derive_local_trace.py` checks these logarithmic coefficients symbolically, including pressure directly; it does not obtain pressure by dividing a Ward equation by H. The power-divergent and boundary terms cancel using the three PV moment identities. The continuum cutoff limit is taken at fixed Lambda first. At finite K, treating summed bare vacuum terms as their continuum covariant action can leave cutoff artifacts, just as in the flat checkpoint.

The consequence at every finite Lambda, after momentum cutoff removal, is the identity

    (T,S)_PV - (T,S)[Gamma_sub]
      = Σc_j (T,S)_positive-reference,j.

This is the common-action bridge in FRW. It is stronger than simply declaring zero local additions to an unnamed subtraction.

## 4. The heavy-field limit: scope and assumptions

To identify the physical positive-reference prescription with the regulator limit, the j>=1 terms on the right must vanish. This requires a state and background assumption; it does not follow from algebraic matching alone.

For the proposed control, a is positive and C-infinity, all derivatives of a and h are bounded on every compact conformal-time interval, and both profiles are exactly constant in the past. Prepare every mass with the corresponding exact past plane wave. For sufficiently large r_j, the canonical frequency is uniformly gapped by a positive multiple of sqrt(k²+r_j). The standard higher-order WKB/Volterra asymptotic argument about this gapped frequency, with exactly matched static initial data, yields uniform symbol estimates on compact intervals under its smooth-symbol hypotheses. In particular the order-six remainder of stress, after order-four subtraction, and the order-four remainder of Q, after order-two subtraction, are bounded in magnitude by

    C_I (k²+r_j)^(-5/2),

up to bounded powers of a absorbed in C_I. One way to obtain these bounds is to construct the mode parametrix through at least eighth order: its higher-order frequency is positive at large r_j; its residual is estimated using bounded profile derivatives; the Volterra Green kernel adds a factor O((k²+r_j)^(-1/2)); and the parametrix expansion of the quadratic observables leaves the stated order-six/order-four remainders. Exact static preparation removes the algebraic initial-state mismatch. The estimate has a finite constant dependent on the compact time interval and sufficiently many profile derivatives. Consequently

    ∫_0^infinity k²dk (k²+r_j)^(-5/2) = 1/(3r_j),

so all heavy renormalized stress/current terms vanish at least as O(r_j^-1). The state-dependent reflection can also be made smaller than any power for the C-infinity, gapped heavy-field scattering problem by repeated adiabatic integration by parts.

This is a standard asymptotic/decoupling argument with its hypotheses stated, not an explicit certified constant or a numerical PV-limit test. It is not uniform over arbitrary families with shrinking transition width, unbounded derivatives, nonstatic finite-time vacuum choices, arbitrary nonsmooth FRW histories, or arbitrarily growing time intervals. A fully formal operator-norm proof of the cited symbol estimates is outside this checkpoint. The finite-Lambda coefficient bridge above does not depend on that unprovided bound constant. Under the stated smooth, static-prepared asymptotic hypotheses, the argument identifies the physical positive-reference stress/current with the declared PV/heat-kernel finite convention. This conditional identification is not a certified uniform-bound theorem for arbitrary states or backgrounds.

One must not use a finite ratio K/Lambda as a substitute for the ordered limits. Nor should a heavy-mass local approximation to the physical field be used through x=0.

## 5. A pressure-sensitive local trace identity

For the physical bare modes,

    -rho_b+3p_b = (Q_bddot+3H Q_bdot)/2-xQ_b.

The complete subtractions give the continuum identity

    -rho+3p = (Qddot+3H Qdot)/2-xQ+A_r,
    A_r = a2(r)/(16pi²).

This is a matched, scheme-dependent trace identity, not a universal cosmological prediction. In FRW,

    Riem²-Ric² = -12H²(Hdot+H²),
    box R = -6H'''-24Hdot²-42H Hddot-72H²Hdot,
    box x = -xddot-3H xdot.

Writing d=Hdot,e=Hddot,f=H''',z=xdot,zz=xddot,

    2pi² A_r = [58H⁴-14H²d-60H²h-42He+15Hz
                  -9d²-30dh-6f+15h²+5zz]/240.

For a static flat geometry, H and all its derivatives vanish, and this reduces to h²/(32pi²)+xddot/(96pi²), exactly the prior flat formula. An isolated point with H=0 need not have vanishing curvature derivatives. This identity can detect pressure defects that an energy ledger cannot detect on a static interval.

There is also an exact finite-comoving-K remainder, avoiding the mistake of using the continuum anomaly in a finite-band test. Let P=K/a, v=P/sqrt(P²+r), and

    J_n(P)=∫_0^P p²dp/(p²+r)^(n/2)
      = r^((3-n)/2) Σ_{ell=0}^{(n-5)/2}
           (-1)^ell binom((n-5)/2,ell) v^(2ell+3)/(2ell+3),

for n=5,7,9,11,13. Define

    c5 = r[-18H²d-6H²h-9He+5Hz-3d²-4dh-f+3h²+zz]/16,
    c7 = -r²[-40H⁴-24H²d+50H²h+5He+10Hz+10dh+f]/32,
    c9 = 7r³[55H⁴+48H²d+10H²h+4He+3d²]/64,
    c11 = -231H²r⁴(3H²+d)/64,
    c13 = 1155H⁴r⁵/256,
    A_K = Σ_n c_n J_n(K/a)/(2pi²).

Then exactly

    -rhoK+3pK = (QKddot+3H QKdot)/2-xQK+A_K.

In conformal time the derivative term is (QK''+2LQK')/(2a²). The finite-K expression is derived by differentiating each fixed-comoving-k subtraction before integrating; no moving-cutoff derivative is omitted. For a=1 its factor is v³[h²/(32pi²)+xddot/(96pi²)]. For a numerical pressure test, estimate Q derivatives from independently resolved output; replacing them with the same mode equations used to compute pressure checks mainly algebra and can conceal correlated errors.

## 6. Analytic smooth background and stable derivatives

Let s=eta/T, f(s)=exp(-1/s) for s>0 and zero otherwise, and

    B(s)=f(s)/[f(s)+f(1-s)],
    a(eta)=exp[A B(s)],   x(eta)=r[2B(s)-1]².

Thus a=1,x=r for eta<=0 and a=exp(A),x=r for eta>=T. All endpoint derivatives vanish; x=0 at eta=T/2. Incoming data at eta=0 are exact:

    u_k(0)=1/sqrt(2sqrt(k²+r)),
    u_k'(0)=-i sqrt(k²+r) u_k(0).

The physical canonical frequency squared can temporarily become negative at some low k. The exact linear ODE remains regular; this does not justify replacing the state with WKB data or declaring a divergent subtraction.

Inside 0<s<1 use B=logistic(g),

    g=-1/s+1/(1-s),
    g^(n)=n![(1-s)^(-n-1)+(-1)^(n+1)s^(-n-1)].

Set P1=B(1-B), P2=P1(1-2B), P3=P1(1-6B+6B²), P4=P1(1-14B+36B²-24B³). Then

    B_s = P1 g1,
    B_ss = P1 g2+P2 g1²,
    B_sss = P1 g3+3P2 g1 g2+P3 g1³,
    B_ssss = P1 g4+P2(4g1g3+3g2²)+6P3 g1²g2+P4 g1⁴.

For numerical stability calculate P1 as exp(-|g|)/[1+exp(-|g|)]², even when the rounded value B is exactly 0 or 1. Do not obtain small endpoint derivatives by subtracting two nearly equal sampled values. At the endpoints and outside the transition use the exact constant branches. Products of exponentially small factors and large rational derivatives can be evaluated with log scaling if necessary. Any deliberate endpoint clipping must have an explicit error bound.

Let b_n=B^(n)(s)/T^n, ell_n=A b_n, y=2B-1. Then

    a'/a=ell1,
    a''/a=ell2+ell1²,
    a'''/a=ell3+3ell1ell2+ell1³,
    a''''/a=ell4+4ell1ell3+3ell2²+6ell1²ell2+ell1⁴,
    x'=4r y b1,
    x''=8r b1²+4r y b2,
    x'''=24r b1b2+4r y b3,
    x''''=24r b2²+32r b1b3+4r y b4.

Useful cosmic quantities are

    H=ell1/a,
    Hdot=(ell2-ell1²)/a²,
    Hddot=(ell3-4ell1ell2+2ell1³)/a³,
    H'''=(ell4-7ell1ell3-4ell2²+18ell1²ell2-6ell1⁴)/a⁴,
    xdot=x'/a,  xddot=(x''-ell1 x')/a²,
    R=6(ell2+ell1²)/a²,
    Rdot=6(ell3-2ell1³)/a³,
    Rddot=6(ell4-3ell1ell3-6ell1²ell2+6ell1⁴)/a⁴.

For w derivatives let Z=k²+ra² and

    Z1=2ra²ell1,
    Z2=2ra²(ell2+2ell1²),
    Z3=2ra²(ell3+6ell1ell2+4ell1³),
    Z4=2ra²(ell4+8ell1ell3+6ell2²+24ell1²ell2+8ell1⁴).

Then

    w'=Z1/(2w),
    w''=Z2/(2w)-Z1²/(4w³),
    w'''=Z3/(2w)-3Z1Z2/(4w³)+3Z1³/(8w⁵),
    w''''=Z4/(2w)-(4Z1Z3+3Z2²)/(4w³)
            +9Z1²Z2/(4w⁵)-15Z1⁴/(16w⁷).

For Delta=a²h-ell2-ell1²,

    Delta'=a²(h'+2ell1h)-ell3-2ell1ell2,
    Delta''=a²[h''+4ell1h'+(2ell2+4ell1²)h]-ell4-2ell2²-2ell1ell3.

Writing w_n for the nth derivative,

    u2'=Delta'/(2w)-Delta w1/(2w²)-w3/(4w²)
          +5w1w2/(4w³)-9w1³/(8w⁴),
    u2''=Delta''/(2w)-Delta'w1/w²-Delta w2/(2w²)+Delta w1²/w³
          -w4/(4w²)+(7w1w3+5w2²)/(4w³)
          -57w1²w2/(8w⁴)+9w1⁴/(2w⁵).

These supply u4 and the pressure subtraction without numerical differentiation of the background.

## 7. Prospective numerical requirements and risks

The inherited illustrative parameters r=4,T=1,A=0.2 are suitable for a first control; use A=0 with identical x for the flat reduction, and a static a=1,x=r case as a vacuum null check. The smooth bump has a short effective transition region, so rT² alone does not quantify its largest derivatives. These are chosen control parameters, not fitted shell data.

The prospective protocol uses K=24,48,96, with K=192 only under its registered refinement condition; momentum-grid and ODE-tolerance refinements must be separate. The precise grid, time range, tolerances, deterministic output nodes and acceptance thresholds belong in a prospective registration before evolution, selected by the numerical owner. Include the transition, the exact zero crossing, and a nonzero static future interval. Retain coherent pressure/current; final occupation energy provides a useful independent future-energy check but is not a replacement for the full stress.

At eta>=T, the exact out decomposition determines alpha_k,beta_k. With omega_f=sqrt(k²+exp(2A)r), the renormalized energy mode is omega_f|beta_k|²/a_f⁴. Pressure and Q generally also contain alpha beta* coherence oscillations. The static out vacuum has zero source in the selected scheme. The exact incoming state and smooth compact transition give a Hadamard-compatible ultraviolet structure; any numerical all-k claim still requires cutoff evidence.

The main numerical danger is subtractive cancellation: bare rho scales as k and local terms cancel through k^-3, while the remaining per-mode stress generally scales as k^-5. A relative mode error epsilon can contaminate an integrated energy by order epsilon K⁴. Therefore a larger K with unchanged precision can worsen the answer; Wronskian conservation alone does not bound correlated quadratic-observable errors. Use stable mode variables or rationalized differences, meaningful absolute errors, compensated summation where useful, and precision/tolerance cross-checks targeted at the high-k band. In the transition, local vacuum-polarization tails may still be O(K^-2) even though final scattering occupation decreases faster than any power; do not infer the stress tail from beta alone.

Recommended registered acceptance logic:

* Verify analytic background derivatives and exact constant branches before evolution, including the zero-crossing neighborhood.
* Compare cutoff, quadrature and independent-solver changes against an absolute floor proportional to r² for rho,p and r for Q as well as a relative tolerance where the result is nonzero. Do not divide by a source that crosses zero.
* Test the finite-K exchange identity with a separately resolved temporal energy derivative and an integrated work ledger.
* Test pressure through the finite-K trace identity above with resolved Q derivatives; report both absolute residuals and a scale-normalized residual. One should not accept merely because energy exchange passes.
* Compare static-future rho against the exact out-basis occupation integral; retain all coherence in p and Q.
* Include deliberate pressure and current omission tests in analytic verification, where they can expose missing terms without creating spurious physical runs.
* State whether tolerances support a finite-band benchmark or a continuum extrapolation. A fitted K^-2 term is an asymptotic estimate, not a certified error enclosure.

A failed gate should produce a bounded diagnostic report, not an automatic declaration of physical finiteness or a shell modification. Predetermined refinement branches may be registered; unregistered refinements must be labelled post hoc.

## 8. What this closes and what it does not

The local FRW coefficient bridge, source pairing, finite-K exchange and pressure-sensitive trace identities are analytic results for this explicit prescription. They close the previously unstated F(r) and R² finite convention in a controlled smooth FRW calculation. They preserve the declared flat matching. The covariant beta choice is supplied by the heat-kernel extension; FRW itself cannot test it. Decoupling relies on the stated smooth, gapped, static-prepared heavy-field assumptions, and no rigorous numerical remainder constant is supplied here.

The finite choices can always be changed by adding a common local action. For a shift deltaV,deltaF,deltaalpha the stress/current shift is exactly the variation in section 1; a nonconstant deltaF must shift Q by -2 deltaF_x R. Holding Q fixed while changing such a stress violates the exchange identity. Constant deltaF shifts FRW stress without changing Q; constant deltaalpha also shifts stress without Q. Ward consistency alone cannot select these coefficients.

This construction does not specify physical HDBLAST gravitational matching, solve junction conditions, evolve backreaction, resolve a singular or nonsmooth archived background, demonstrate net energy transfer into a self-consistent matter sector, thermalize produced quanta, or produce a sustained hot expanding universe. Those remain separate physical tasks. No account handoff or external publication is authorized by this draft itself.

## Sources and provenance

Primary project sources read directly: the pinned flat `THEORY_AND_MATCHING.md`, `FRW_ACTION_BRIDGE.md`, `DIRECT_CUTOFF_APPENDIX.md`; the inherited FRW `FRW_COUNTERTERMS.md`, `THEORY_REGULARITY.md`, `NEXT_EXPERIMENT.md`, and `code/verify_frw_counterterms.py`. The algebra in the present `derive_local_trace.py` is new direct symbolic derivation of the displayed local formulas.

Relevant literature provenance is inherited from the project's hash-verified literature reviews, not newly fetched primary full text in this turn: Ferreiro–Navarro-Salas, arXiv:1812.05564 (off-shell reference mass and adiabatic finite shifts); Ferreiro–Monin–Torrenti, arXiv:2311.08986 (adiabatic scales and gravitational couplings, constant physical mass); Molina-París–Anderson–Ramsey, gr-qc/9908037 (second-order Q, fourth-order stress and exchange with a mean field); Markkanen–Tranberg, arXiv:1303.0180 (local counterterms and common action); Junker–Schrohe, math-ph/0109010 (finite adiabatic order versus Hadamard regularity). The familiar scalar heat-kernel coefficients are stated explicitly above so the coefficient calculation does not depend on an inaccessible citation. No new literature-based priority or discovery assertion is made.
