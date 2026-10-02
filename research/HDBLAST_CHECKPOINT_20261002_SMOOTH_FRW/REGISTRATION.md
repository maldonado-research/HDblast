# HDBLAST smooth curved-source control: prospective registration

Prepared 1 October 2026, America/Los_Angeles; folder uses 2 October UTC.
Author: Ricardo Maldonado, with disclosed OpenAI Codex/Astra-assisted derivation,
implementation and independent internal review. This is not external peer review.

## Chronology and scientific question

This registration and source are committed before any new physical mode evolution.
Analytic/symbolic derivations were completed before registration; their logs are
included and must not be described as prospective numerical outcomes. Numerical
success is not presumed. The inherited curved-readiness audit is pinned at
7e508d946e9ecdc542ca3059263d93132321838b and its subsequent scientific merge is
0949317d5e0310535caa50576a83d4487815dd7a. Main remains separately preserved.

Question: in an explicitly stated finite action convention, does the fully smooth,
exact-static-prepared scalar control give resolved stress, pressure and paired
mass current through x=0, with cutoff evidence and independently checked physical
mode evolution, exchange and pressure trace? This is a prescribed-background
benchmark. Its background is not a solution of the coupled HDBLAST shell equations.

COMMON_ACTION_SMOOTH_FRW.md defines the action sign, reference r, finite potential,
xR, R2 terms and the declared Weyl/Euler extension. Its explicit FRW coefficient
bridge is algebraic; heavy-field decoupling is conditional on the stated smooth
state/asymptotic hypotheses and has no certified uniform remainder constant.
These qualifications remain part of the scientific interpretation even if every
numerical gate passes. PROTOCOL_SOURCES.json records exact pre-execution bytes.
Independent review found no blocker to executing this qualified control.

## Nonadaptive reporting rules

Preserve this file and original source snapshot after results. Every failed run,
mechanical correction and later refinement must identify its actual source hash,
protocol hash, settings and order. No threshold, action choice, input or hypothesis
may be changed retrospectively. An implementation repair must preserve the first
failure and describe whether any physical evolution had already occurred. A new
scientific setting beyond the registered contingency needs a separate prospective
addendum. Record absolute errors, finite precision and cutoff limits without a
claim of a certified continuum answer. Include null controls and deliberate
incorrect-formula detections; an import or zero-test run is not validation.

# Prospective numerical matrix: exact incoming vacuum on a prescribed smooth FRW background

This protocol must be committed with the implementation before research execution.
No result or successful convergence is presumed. The matrix stops after its one
registered UV extension even if a gate fails. Failed directories are retained.

## Physical problem and prescription

Use conformal time, r=4, T=1, amplitudes A=0 and A=0.2, and the specified C-infinity
compact step B: a=exp(A B), x=r(2B-1)^2. Integrate eta in [0,1.25]. The exact
incoming vacuum at eta=0 has u=1/sqrt(2 sqrt(k^2+r)), u'=-i sqrt(k^2+r)u.
Evolve the complex physical equation u''+(k^2+a^2 x-a''/a)u=0 with DOP853.
There is no WKB replacement for the physical mode, no Wronskian normalization
after evolution, and no clipping of low-k modes. Minimal coupling is understood.

Per mode, z=u'-H u, rho=(|z|^2+(k^2+a^2 x)|u|^2)/(2a^4),
p=(|z|^2-(k^2/3+a^2 x)|u|^2)/(2a^4), Q=|u|^2/a^2.
Subtract the complete fourth-order rho,p and second-order Q expansions from the
archived exact-rational verifier, using long-double Taylor jets to order five.
The subtraction reference frequency is sqrt(k^2+a^2 r), always positive.
The momentum measure is k^2 dk/(2 pi^2), with a fixed comoving cutoff. The
paired mass-squared current is Q/2. No coherent term is discarded.

The action matching document supplies the common finite-action convention.
According to its constructive Pauli-Villars/DeWitt-Schwinger matching, the base
subtraction needs no additive local correction: V,Vx,Vxx,F,Fx and the R^2
coefficient vanish at x=r in the stated local derivative basis. The C^2 and Euler
extensions are explicitly declared there; spatially flat FRW cannot measure
their independent coefficients. Numerical exchange alone does not prove
covariance or establish this matching.

## Fixed matrix and numerical settings

For each amplitude, use cutoff K=24,48,96 on common Gauss-Legendre panels of
width 2. The main grid uses 8 nodes per panel. For each k set the initial state
above and solve all physical modes together. The array RMS error controller is
supplemented by max_step=min(1/80, phase_step/sqrt(kmax^2+exp(2A)r)).

1. Primary: 8 nodes/panel, rtol=2e-11, atol=2e-13, phase_step=0.2.
2. Tighter solver: 8 nodes/panel, rtol=2e-13, atol=2e-15, phase_step=0.1.
3. Refined momentum quadrature: 12 nodes/panel and the tighter solver settings.

Use 801 uniform output times on [0,1.25]. The nested 401-node grid tests time
quadrature and finite-difference diagnostics. Time sampling does not control
the adaptive ODE steps, whose independent refinement is specified above.

If and only if any tight-solver K48-to-K96 observable cutoff gate fails, run the
registered K192 extension at both 8 and 12 nodes/panel with tight ODE settings,
and an additional 8-node solve with phase_step=0.05 at the same tight tolerances.
Compare K96-to-K192, quadrature at K192, and K96 evaluated within the original
and extended solves. Compare the two phase step sizes at K192 to test the new
high-k band independently. Do not increase K again or relax tolerances. A failure is a
result requiring a separately registered follow-up, not a successful source.

Run an exact-static negative control B identically zero with A=0, K96,
8 nodes/panel and tight ODE settings. Its exact subtracted rho,p,Q vanish.

## Diagnostics and gates

All maxima below include every output time and apply independently to rho,p,Q.

- Observable solver, quadrature and final cutoff changes must each satisfy
  max|f1-f2| <= 2e-5+0.005 max(max|f1|,max|f2|). Report the absolute change,
  scale, normalized change and threshold; the absolute allowance accounts for
  floating-point cancellation of large vacuum terms. K24-to-K48 is reported
  but is not the final cutoff acceptance comparison. Report the ratio of the
  K48-to-K96 and K24-to-K48 maximum changes; do not infer an asymptotic power
  from this alone or gate a near-zero sign-changing shell by its ratio.
- Raw normalized Wronskian max|(u u'^*-u' u^*)/i-1| <= 2e-9 for primary and
  <=2e-10 for all tighter runs. No rescaling enforces this identity.
- Compute the renormalized conformal-energy ledger with cumulative Simpson:
  E=a^4 rho, L=int a^4[H(rho-3p)+x'Q/2]deta.
  max|E-E(0)-L| <= 1e-5+0.005 max(max|E-E(0)|,max|L|).
  The fine/coarse L difference has the same absolute-plus-relative rule.
- Independently accumulate an extra ODE state per mode,
  Lbare'=H[-|u'-Hu|^2+(k^2+2a^2 x)|u|^2]+a^2 x'|u|^2/2.
  Compare the renormalized E change with Lbare minus the analytic subtraction
  energy change under the same ledger tolerance. Subtraction cancels
  algebraically in this diagnostic: it checks mode integration against
  integrated work, while the separate Simpson ledger checks the paired source.
- With fourth-order five-point finite differences (one-sided at the edges),
  max|rho'+3H(rho+p)-x'Q/2| <= 2e-4+0.01 max(max|rho'|,
  max|3H(rho+p)|,max|x'Q/2|). Report 401- and 801-node results.
- Report the counterterm exchange residual from analytic Taylor derivatives as
  a numerical subtraction implementation diagnostic, not an independent
  physical solver test. Report the deliberately source-omitted ledger residual
  as a diagnostic; it is not independently forced to exceed a threshold.
- Independently test pressure via the exact finite-comoving-cutoff trace
  remainder derived in the finite-action matching document. In conformal time,
  compare -rho+3p against (Q''+2H Q')/(2a^2)-xQ+A_K. Integrate the independently
  derived analytic moments A_K, rather than using the infinite-cutoff anomaly.
  Direct physical-equation Q derivatives minus subtraction-jet derivatives
  must agree within 2e-6+5e-6 times the larger maximum side. This tests the
  subtraction and pressure implementation; the bare trace identity is built
  into those direct derivatives. Separately use five-point sampled Q' and Q''
  at 401 and 801 times; require the fine sampled test within
  2e-4+0.01 times the larger maximum side, and report both resolutions.
- In the exactly static future eta>=1, compare full rho with the exact
  occupation energy int omega_final |beta|^2/a_final^4 using the ordinary
  observable tolerance. Report full p,Q, their occupation-only estimates, and
  retained coherent differences. beta=(omega_final u-i u')/sqrt(2omega_final),
  whose magnitude is phase independent. No relative beta gate is used for
  negligible occupations.
- The exact-static control has max|rho|,max|p|,max|Q| <=2e-5 and satisfies the
  tight Wronskian gate. Its residual quantifies the floating-point floor.
- All arrays, samples and strict JSON values must be finite. Save full physical
  complex modes, subtracted mode quantities, weights, integrated observables,
  background values and ledgers as numeric NPZ without object arrays; load
  artifacts with allow_pickle=False. Keep immutable run directories, failed
  outputs, protocol SHA256 and source SHA256.

## Independent solver check

The independent implementation uses four real components and Radau with an
analytic Jacobian and a separately implemented compact-step background.
Parameters: rtol=2e-12, atol=2e-14,
max_step=min(1/100,0.1/sqrt(k^2+exp(2A)r)). Check both A values at
k=0,0.5,2,8,24,48 and times 0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1,1.25.
Compare physical u,u' against DOP853: maximum absolute differences <=2e-8.
Compare raw rho,p,Q using 2e-7+2e-8 times the larger absolute quantity.
Report each raw Wronskian, endpoint alpha,beta, and the separate source/current
identity test. This solver check is necessary in addition to the integrated
matrix gates. The independent review may further constrain research claims.

## Scope

Passing every gate demonstrates only this regulated/subtracted prescribed
smooth-background control at the documented tolerances and finite-action
convention. Cutoff agreement at finite K is numerical evidence, not proof of
the infinite-cutoff limit. It does not change the archived shell, establish
backreaction or an energy-consumption model, produce thermal radiation, or
demonstrate a new expanding cosmological solution.
