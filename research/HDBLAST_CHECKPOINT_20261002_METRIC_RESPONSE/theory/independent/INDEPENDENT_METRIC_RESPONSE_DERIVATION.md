# Independent homogeneous metric-response derivation

Prepared 2 October 2026. This is an independent derivation by a separate task
agent, within the same AI-assisted project. It is not external peer review.
No new physical response, source, spectrum, mode, stress, or quadrature has
been evaluated. The companion scripts perform exact symbolic algebra only.

## Convention and perturbation

Use the inherited real minimally coupled scalar, signature (-+++), physical
mass squared x and positive reference r both fixed at 2, H=1, and incoming BD
state. In the same conformal coordinate eta<0, vary

    a(eta) = a0(eta) [1 + epsilon h(eta)],
    a0 = L = -1/eta,  L' = L^2,  M = a0^2 r = 2L^2.

Here delta means the derivative with respect to epsilon at zero; h is the
unit-amplitude perturbation. No perturbed state is selected at a later time.
For a compact h with a source-free past, delta v and delta v' start at zero.
All momentum cutoffs below are fixed comoving cutoffs. The reference r is
held fixed under metric differentiation; it does not track the scale factor,
curvature, mass law, cutoff, or stationary root.

The canonical equation and physical operators are

    v'' + [k^2 + a^2 x - a''/a] v = 0,
    Dv = v' - (a'/a) v,
    rho_b = [|Dv|^2 + (k^2+a^2 x)|v|^2] / (2a^4),
    p_b = [|Dv|^2 - (k^2/3+a^2 x)|v|^2] / (2a^4),
    Q_b = |v|^2/a^2.

The exact reference v=exp(-ik eta)/sqrt(2k) has a flat canonical frequency.
Its physical stress remains the minimal stress, including Dv and its mass.

The metric source in the forced canonical equation is

    sigma = 4L^2 h - h'' - 2Lh',
    delta v'' + k^2 delta v = -sigma v.

Equivalently delta R=6(h''+2Lh'-4L^2h)/a0^2 and
sigma=-a0^2 delta R/6. This simple relation is specific to x=R0/6 at the
reference. It does not turn the model into a conformally coupled theory.

## Physical operator and prefactor contacts

At fixed eta and k, delta L=h', delta(a^2 x)=2Mh, and

    delta(Dv) = delta v' - L delta v - h' v,
    delta Q_b = [2 Re(v* delta v)-2h|v|^2]/a0^2,
    delta rho_b = {2 Re[(Dv)* delta(Dv)]
                  +(k^2+M) 2Re(v* delta v)+2Mh|v|^2}/(2a0^4)
                 -4h rho0,b,
    delta p_b = {2 Re[(Dv)* delta(Dv)]
                -(k^2/3+M) 2Re(v* delta v)-2Mh|v|^2}/(2a0^4)
               -4h p0,b.

The stars denote complex conjugation. These terms are all required even when
the canonical forcing is propagated with the previously calibrated scalar
kernel. Missing h' in delta(Dv), the explicit mass variation, or the metric
prefactors changes the observable.

For an analytic reduction, put A=delta|v|^2 and B=A'. The fixed-state forced
solution gives A''+4k^2 A=-sigma/k and delta|v'|^2=-k^2 A. These analytic
identities must not be imposed as numerical mode projections. Direct modes
can evaluate each complex bilinear independently. With the momentum measure
k^2 dk/(2pi^2), the canonical bare variations are

    a0^4 delta rho_b = [3L^2 A-LB+(Lh'+2L^2h)/k]/2
                       -4h(2k^2+3L^2)/(4k),
    a0^4 delta p_b = [-(4k^2/3+L^2)A-LB+(Lh'-2L^2h)/k]/2
                     -4h(2k^2/3-L^2)/(4k),
    a0^2 delta Q_b = A-h/k.

## Full inherited subtraction variation

Use the full smooth-FRW fixed-r second-order Q and fourth-order stress
subtraction. The source mass offset x-r is zero throughout this metric
variation. Set omega=sqrt(k^2+M), C=a''/a, and

    U = -C/(2omega)-omega''/(4omega^2)
        +3omega'^2/(8omega^3),
    V = -U^2/(2omega)-U''/(4omega^2)+omega''U/(4omega^3)
        +3omega'U'/(4omega^3)-3omega'^2U/(4omega^4),
    b = -k^2/3-M,  c=1-b/omega^2,
    J2 = Lomega'/omega^2+omega'^2/(4omega^3),
    J4 = LU'/omega^2-2Lomega'U/omega^3
         +omega'U'/(2omega^3)-3omega'^2U/(4omega^4).

The canonical subtraction densities R,P,S are

    R0=omega/2,
    R2=[L^2/omega+J2]/4,
    R4=[U^2/omega-L^2U/omega^2+J4]/4,
    P0=k^2/(6omega),
    P2=[cU+L^2/omega+J2]/4,
    P4=[cV+bU^2/omega^3-L^2U/omega^2+J4]/4,
    S0=1/(2omega),  S2=-U/(2omega^2).

R=R0+R2+R4, P=P0+P2+P4, S=S0+S2. Their physical prefactors are a^-4,
a^-4, and a^-2. Thus physical metric subtraction variations are
(a0^-4)(delta R-4hR), (a0^-4)(delta P-4hP), and
(a0^-2)(delta S-2hS).

For z=delta omega=Mh/omega, direct differentiation gives

    delta C=h''+2Lh',
    delta U= -delta C/(2omega)+Cz/(2omega^2)
             -z''/(4omega^2)+omega''z/(2omega^3)
             +3omega'z'/(4omega^3)-9omega'^2z/(8omega^4).

Differentiate U, U', U'', V, b, c, J2, and J4 as well as omega, L, and M.
`metric_wkb.py` performs this with dual Taylor jets through a''''. It never
uses the scalar-response contact inventory and never computes a physical
mode. The pure audit compares its complete baseline and directional R,P,S,
U,V and the first two S and W2 derivatives with separately expanded chain rules.

## Exact finite-band bridge to the scalar kernel

A scalar mass source d_eff=sigma/a0^2 produces the same forced canonical
modes on fixed geometry. Denote its separately renormalized response by the
subscript scalar, and use the actual renormalized finite-band baselines
Q0,K, rho0,K and p0,K in the following formulas:

    delta Q_metric,K = delta Q_scalar,K - 2h Q0,K + Cq,K/a0^2,
    delta rho_metric,K = delta rho_scalar,K - 4h rho0,K
        +(h''+4Lh')Q0,K/(2a0^2)+Cr,K/a0^4,
    delta p_metric,K = delta p_scalar,K - 4h p0,K
        -h''Q0,K/(2a0^2)+Cp,K/a0^4.

Each C is an individually convergent rational omega integral. The complete
exact Laurent coefficients are in `FINITE_K_CONTACT_INVENTORY.json`.
For each coefficient c_n multiplying omega^-n, n=5,7,...,15, integrate
c_n I_n(K,M)/(2pi^2), where

    vK = K/sqrt(K^2+M),
    I_n(K,M) = M^((3-n)/2)
      sum[ell=0..(n-5)/2] (-1)^ell binomial((n-5)/2,ell)
        vK^(2ell+3)/(2ell+3).

These primitives come from the variable v=k/sqrt(k^2+M). The pure proof
checks each transformed derivative, its zero endpoint and continuum
endpoint. Replacing them with continuum contacts changes a finite-band
comparator. Finite-K rho0,K+p0,K and Q0,K derivatives generally do not vanish.

After cutoff removal, the inherited continuum baselines are

    Q0=1/(12pi^2), rho0=11/(960pi^2), p0=-rho0,

and the exact continuum bridge contacts reduce to

    Cq,infinity=(12L^2h-2Lh'-h'')/(48pi^2),
    Cr,infinity=-L(61L^2h'+3Lh''-3h''')/(240pi^2),
    Cp,infinity=-(7L^3h'-44L^2h''-3Lh'''+3h'''')/(720pi^2).

This is an exact first-order bridge about this particular reference, not a
scalar-kernel claim about arbitrary metric perturbations or shifted roots.

## Ward and trace identities

With x constant, the bare and all three subtraction grades conserve exactly
at each fixed comoving k. They produce the finite-band metric Ward identity

    delta rho_K' + 3L(delta rho_K+delta p_K)
      + 3h'(rho0,K+p0,K) = 0.

There is no mass-work RHS and no additional W4 subtraction drift in this
fixed-x=r case. The pure proof explicitly checks the graded subtraction Ward
identities and their full metric directional variation before specializing
M=2L^2. A time-dependent physical cutoff would instead add a moving-boundary
term. In the continuum, rho0+p0=0, and the last baseline contact vanishes.
Neither the Ward identity nor a trace relation defines the direct stress.

The inherited trace is T=-rho+3p=-xQ-(box Q)/2+A. At finite K, metric
variation of the d'Alembertian acts on the time-dependent baseline Q0,K:

    delta T_K = -2delta Q_K
        +(delta Q_K''+2Ldelta Q_K')/(2a0^2)
        -h(Q0,K''+2LQ0,K')/a0^2 + h'Q0,K'/a0^2 + delta A_K.

The two baseline derivative contacts may be dropped only after a justified
continuum limit. A constructive independent finite-band anomaly is obtained
from the full subtraction itself. Define the canonical local remainder

    Ctrace=(S''-2LS'-2L^2S)/2-MS+R-3P.

Then

    delta Ctrace=(delta S''-2Ldelta S'-2h'S'
                 -2L^2delta S-2h''S)/2
                 -2MhS-Mdelta S+delta R-3delta P,
    delta A_K = integral_0^K dmu (delta Ctrace-4h Ctrace)/a0^4.

Its rational inventory is recorded in `EXACT_PROOF.json`. In the continuum,

    delta A=(-116L^4h+94L^3h'+47L^2h''-3h'''')/(240pi^2 a0^4).

The proof independently varies the inherited cosmic-curvature anomaly and
integrates the complete directional subtraction trace to obtain the same
answer. Finite-K trace comparisons should retain the rational primitives.

## Gauge, state, and scope

The residual conformal-coordinate time shift eta -> eta+epsilon T, with T
constant, has h=LT and sigma=0. The state and temporal boundary transform
with the geometry; its mode change is a harmless constant phase. At fixed
comoving K the correct observable change is T times the derivative of the
finite-band baseline, rather than an artificial zero. The pure audit checks
this directly for R,P,S. For the continuum de Sitter baselines, all three
physical responses cancel exactly to zero.

A general compact h produces actual delta R and is not this residual gauge
shift. The other homogeneous sigma=0 solution h proportional to eta^4 has
zero linear scalar-curvature variation but can encode the radiation-like
homogeneous metric sector; sigma=0 alone does not identify a pure gauge mode.
A nonzero smooth compact h cannot solve sigma=0 everywhere with a zero past.

For the inherited mass law, metric variation at fixed scalar field gives
only delta j=(br/2)delta Q. It supplies no scalar mass-law variation contact.
Additional shell or bulk couplings must be specified in their own action.

The analytic route is ready for a separately registered bounded calibration:
one route uses scalar memory plus the finite-band bridge, and the other uses
actual forced complex modes plus the full dual WKB operator variation. Fix
compact profiles and their analytic jets, time range, state, precision,
cutoff and refinement ladders, gates, wall/memory budgets, controls, failure
preservation, and an immutable public source freeze before either route is
run. Check both full stresses directly, Q derivative jets, unchanged state,
Wronskian diagnostics without projection, independent finite-band Ward and
trace ledgers, and the coordinate-shift analytic control.

The contact tail can be enclosed from exact omitted I_n tails. A continuum
bound for the scalar part needs derivative budgets for sigma, not simply the
previous B-source budgets. Stronger precision and independent quadrature
refinement may be needed because metric contacts multiply UV-divergent bare
terms before subtraction. Tail bounds alone do not certify mode, time,
momentum-quadrature, rounding, or ledger errors.

This establishes the local homogeneous metric-response convention about the
exact benchmark. It does not evaluate the response at an actual shifted
stationary root, specify gravitational EFT couplings, supply bulk/shell
boundary conditions, solve constraints and gauge sectors, or establish
coupled evolution or stability. It establishes no higher-dimensional blast,
Big Bang origin, heating, particle production, discovery, or physical fit.

## Proof and preservation

The pinned Python 3.12.14 / SymPy 1.14.0 pure audit uses explicit runtime
checks that remain active under Python optimization. It proves 55 identities
and rejects 10 exact symbolic mutation witnesses. The module SHA256 is
recorded in the receipts. The inherited input files are faithful copies of
already public science and are hashed in each proof receipt.

The first proof attempted to test a negative odd exponent with `power % 2 ==
-1`. Python's modulo gives +1 for these odd negative integers. The failed
source and an actual failed replay are preserved under
`development/first-modulo-guard/`; correcting that algebra guard changes no
physical source, operator, prescription, gate, or result.

Before public registration, ordinary proof execution was made read-only.
`explore_contacts.py` now exports its inventory only to an explicit external
`--output` target (or an explicit `INVENTORY_OUTPUT` exec namespace). Earlier
actual proof receipts, logs, exit statuses, and support sources are preserved
under `development/pre-read-only-export/`. This preparation change alters no
formula, profile, operator, subtraction convention, numerical gate, or mode.

All six final normal/optimized proof runs were repeated with an audit hook
that rejects writes anywhere in the source tree and disables bytecode writes.
Their outputs were written to a separate directory. All 63 existing source
files retained identical bytes during the run. An explicit external inventory
export also matches the previous exact inventory byte for byte. The actual
verification receipt is in `evidence/read_only_preparation/READ_ONLY_REPLAY.json`.
