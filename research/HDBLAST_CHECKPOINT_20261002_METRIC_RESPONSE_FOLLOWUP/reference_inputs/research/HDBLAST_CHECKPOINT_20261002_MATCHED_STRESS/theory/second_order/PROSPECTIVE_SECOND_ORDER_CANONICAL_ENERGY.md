# Prospective quadratic canonical-energy calibration

Prepared 2 October 2026, after the public matched-stress input freeze
`4a5dad8dda6d0a57a9cf88cc8b4a9b8feeaea45d`.
**Analytic follow-up proposal only. This is not a public registration, an
executed numerical energy experiment, or a computed energy yield.** Its
candidate grids and gates require implementation, independent review, and a
new public source/input freeze before any physical evaluation. The current
stress experiment neither registers nor executes this proposal.

The next bounded advance should measure the positive coefficient of
post-pulse canonical excitation energy in three independent ways. This
resolves a concrete interpretive gap: first-order minimally coupled stress
is coherent polarization and can have either sign, while canonical excitation
energy first appears at second order. Their different perturbative orders
and operators must remain explicit.

## 1. Fixed model and expansion convention

Keep the inherited real, minimally coupled scalar, incoming BD state, positive
reference `r=2H²`, fixed de Sitter metric `a=-1/(H eta)`, and fixed comoving
cutoff K. Set H=1 for the numerical calibration, so r=2 and
`h=a'/a=-1/eta`, `h'=h²`. There is no metric or initial-density-matrix variation.
With `v=a chi`, the exact reference canonical equation is `v''+k²v=0`.

Prescribe an external source exactly linear in the amplitude:

    s(eta)=a² delta x=epsilon f(eta),  epsilon=10^-4,
    u=eta+4,  B(u)=exp(1-1/(1-u²)) for |u|<1, otherwise 0,
    f=B(u) or u B(u).

The source-free initial time is -6, support is (-5,-3), and the two post-pulse
observations are -2.5 and -1.5. These retain the preceding controls; they are
candidate inputs for a new registration, not permission to reuse an old run
as this experiment. A nonlinear mass law would introduce an order-epsilon²
source as well; that different experiment is outside this proposal.

Use coefficients without factorials:

    v=v0+epsilon v1+epsilon² v2+...,
    alpha=1+epsilon alpha1+epsilon² alpha2+...,
    beta=epsilon beta1+epsilon² beta2+...,
    Ecan(epsilon)=epsilon² E2+O(epsilon³).

Thus E2 is the quadratic coefficient, not the second amplitude derivative.
Report E2 throughout and label any multiplication by epsilon² explicitly.
With dimensionless a, `s` and epsilon have mass dimension two, `f` is
dimensionless, F below has mass dimension minus one, and Ecan has mass
dimension four. H=1 fixes these units rather than changing the dimensions.

## 2. Positive spectral quantity and independent mode route

Define `F(omega)=integral f(t) exp(-i omega t) dt`. The compact source is
extended by zero only for this Fourier calculation. This does not extend the
physical de Sitter chart through eta=0. After the pulse,

    v=(alpha exp(-ik eta)+beta exp(+ik eta))/sqrt(2k),
    alpha1=-i F(0)/(2k),  beta1=i F(2k)/(2k).

The free canonical Hamiltonian gives, per comoving volume,

    E2,K = (1/(2pi²)) integral_0^K k³ |beta1(k)|² dk
         = (1/(8pi²)) integral_0^K k |F(2k)|² dk,
    E2   = (1/(32pi²)) integral_0^infinity omega |F(omega)|² domega.

E2 is finite for every smooth compact f and strictly positive for real,
nonzero f. Vanishing would force F to vanish away from zero and then
everywhere by continuity and Fourier uniqueness. This is the usual positive
quadratic spectral functional, not a new particle-production theorem.

The independent mode producer evolves `v1=v0 u1` by

    u1'=w1,  w1'=2ik w1-f,  u1(-6)=w1(-6)=0,
    beta1=exp(-2ik eta) w1/(2ik),
    alpha1=u1-w1/(2ik).

It must preserve the complex amplitudes and verify the post-pulse extraction
at both observation times. To avoid a spurious low-k numerical division,
archive the bounded spectral amplitude
`G=-exp(-2ik eta) w1=2k beta1/i=F(2k)` and evaluate the energy as
`integral k |G|²/(8pi²)`. These are algebraically equivalent, not an extra
mode approximation. An absolute-square check alone misses a reversed phase
or a global sign error; complex G must be compared independently with F.

The first-order Wronskian condition is `2Re(alpha1)=0`. The next condition is
`2Re(alpha2)+|alpha1|²-|beta1|²=0`. A first-order truncation is not an exactly
normalized finite-amplitude state. E2 nevertheless needs only beta1.

## 3. Matched work requires a background drift

Here `q_M=a² delta Q_ren/epsilon` denotes only the induced linear response;
it does not contain the background Q0. Inherit the common smooth-FRW local
action and fixed-r PSAR prescription, including its variance and stress
matching. Let `M=a sqrt(r)`. For a constant auxiliary scale mu,

    q_M=q_mu-ln(M/mu) f/(8pi²),
    W2_M=(1/2) integral f' q_M deta,
    D=(1/(32pi²)) integral h f² deta,
    E2=W2_M-D.

The auxiliary decomposition changes no finite action: mu cancels in q_M.
For the stated Fourier convention the retarded multiplier is
`[ln(|omega|/mu)-1]/(8pi²)+i sgn(omega)/(16pi)`. Its absorptive part yields
the positive E2. The contact contribution to work is
`-integral ln(M/mu) (f²)'/(32pi²)=D`, using the vanishing endpoint jets.
A constant real contact makes zero compact work, but this contact depends
explicitly on the background and cannot be discarded as a total derivative.

Use an exact finite-K comparison first:

    A_K=asinh(K/M)-K/sqrt(K²+M²),  V=K/sqrt(K²+M²),
    q_M,K=[f A_K-integral_-6^eta
             f(t) {1-cos[2K(eta-t)]}/(eta-t) dt]/(8pi²),
    A_K'=-h V³,
    D_K=(1/(32pi²)) integral h f² V³ deta,
    (1/2) integral f' q_M,K deta - D_K = E2,K.

Evaluate `1-cos(z)` as `2sin²(z/2)` near zero and combine subtraction and
bare terms before using them. The removable time-kernel limit at t=eta is
zero. This route integrates the momentum kernel analytically and uses no
Fourier or forced-mode arrays. Using D rather than D_K at finite K changes
the tested identity. A physical-momentum cutoff would be a different
regulator and is not admitted.

At quadratic order the physical source work is instead

    (1/(2epsilon²)) integral a⁴ delta x' delta Q_ren deta
       = W2_M - integral h f q_M deta.

The physical Ward identity additionally includes expansion work,
`(a⁴rho)'=h a⁴(rho-3p)+(a⁴Q/2)x'`. Neither W2_M nor its unsubtracted source
piece alone is physical energy gain. The second-order induced Ward equation
for this source is `rho2'+3h(rho2+p2)=Q1 (delta x1)'/2`; an order-two mass-law
source would add `Q0 (delta x2)'/2` and must be registered separately.

## 4. Why this does not yet compute second-order physical stress

After the pulse, local subtraction terms and anomalies cancel between the
evolved and reference states: their geometry, mass, reference r, and mass jets
are identical. Let `Delta q=Delta<v²>` retain all coherence. The exact
post-pulse operator relation is

    a⁴ Delta rho=Delta Ecan-(h/2) Delta q'+(3h²/2) Delta q.

At first order, Ecan has no term; the variance and its derivative produce
coherent minimal stress. At second order, with the same coefficient convention,

    a⁴ rho2=E2-(h/2)q2'+(3h²/2)q2,
    a⁴ p2=E2/3+q2''/6-(h/2)q2'-(h²/2)q2.

These obey the source-free Ward identity when h'=h² and E2'=0. They do not
fix the sign of rho2 or p2 without q2. In particular p is not assumed positive,
negative, constant in sign, or equal to rho/3. The full q2 requires
`|v1|²+2Re(v0* v2)`; merely squaring v1 is incomplete.

There is a concrete infrared obstruction to replacing q2 with occupations.
For `F0=integral f`, the positive pulse has F0>0 and

    beta1 ~ i F0/(2k),  alpha1 ~ -i F0/(2k),
    Qocc^(2) = epsilon²/(2pi² a²) integral k |beta1|² dk
             ~ epsilon² F0²/(8pi² a²) integral dk/k.

This isolated occupation variance diverges, while Ecan is infrared finite.
In the full variance, `alpha1 beta1* ~ -F0²/(4k²)` cancels the leading
`|beta1|²`. Moreover beta2 must be retained. Equivalently,
`u1 -> F1-eta F0`, and the Volterra kernel `sin(k Delta)/k -> Delta`
keeps u1,u2 bounded on every fixed finite observation interval. The full
q2 integrand is proportional to `k[|u1|²+2Re(u2)]` and is integrable there.
Occupation and coherence must be combined with the same infrared regulator
before its removal; adding separately divergent integrals is not defined.
The signed uB pulse has F0=0 at this order, which is not a guarantee for
higher orders. The independent analytic audit is archived alongside this
proposal. Dephasing, a late-time limit, or a new state cannot silently be
interchanged with this cancellation.

The proposed experiment therefore computes E2 and the matching work identity,
not rho2, p2, a thermal spectrum, a detector response, or a radiation-fluid
equation of state. A later second-order physical stress calculation should
retain v2 and the complete inherited nonlinear subtractions during the pulse.

## 5. Constructive spectral endpoint bounds

For any artificial infrared lower cutoff kappa,

    0 <= E2,[0,kappa] <= ||f||_1² kappa²/(16pi²).

The preferred implementation includes k=0 through the continuous weighted
integrand `k|F(2k)|²`; Gauss nodes avoid evaluating a singular beta itself.
If a positive lower cutoff is introduced, its value and the above allowance
must be frozen and reported. Omitting it from the ledger is a failure.

Repeated integration by parts gives, for integer n>=2,

    |F(2k)| <= ||f^(n)||_1/(2k)^n,
    0 <= E2-E2,K <= ||f^(n)||_1² /
             [8pi² 2^(2n)(2n-2) K^(2n-2)].

For n=3 this is `||f'''||_1²/(2048pi² K⁴)`. The already proved rational
global derivative certificates give `||B'''||_1<=97` and
`||(uB)'''||_1<=96`. These are conservative analytic envelopes, not fits to
an evaluated spectrum. Archive their source/certificate hashes when adopted.
No new source-norm sampling is needed for these bounds.

Compare all three producers at the same finite K. Then apply this positive
spectral omitted-band bound to the continuum claim. Integrating the older
variance-tail bound against |f'| would give a valid but less sharp work bound;
it should not replace the exact finite-K drift in the comparison.

For Fourier uncertainty `|F-Fhat|<=eF`, the energy-integrand error is bounded
by `k[2|Fhat|eF+eF²]/(8pi²)`, with an analogous mode bound. Quadrature and
finite time-step allowances need separate ledgers. Nested refinement or
ordinary library error estimates are diagnostics, not interval-certified
continuum integration; only the displayed omitted-band bound is rigorous
under its analytic assumptions.

## 6. Candidate finite experiment, requiring a new public freeze

Implement the following bounded design, inspect all code and error ledgers,
and freeze exact sources, library versions, quadratures, grids, budgets,
validators, mutation thresholds, and stopping rules before running it. The
numbers below are a concrete candidate design, not measured feasibility or
an already accepted accuracy result.

| Setting | Candidate coarse | Candidate fine |
|---|---|---|
| Fixed comoving K | 64, 128, 256 | same |
| Fourier time panels on [-5,-3] | width 1/64, Gauss 12 | width 1/128, Gauss 16 |
| Spectral momentum panels | width 1/2, Gauss 16 | width 1/4, Gauss 20 |
| Independent forced-mode time panels | width 1/128, forcing Gauss 8 | width 1/256, forcing Gauss 8 |
| Mode momentum panels | width 1/2, Gauss 16 | width 1/4, Gauss 16 |
| Work outer time panels on [-5,-3] | width 1/32, Gauss 12 | width 1/64, Gauss 16 |
| Work inner time panels, clipped to [-5,eta] | maximum width 1/64, Gauss 12 | maximum width 1/128, Gauss 16 |

The mode route uses an exact homogeneous exponential update and integrates
its forcing independently; it must not import Fourier values. The work
route uses the regular finite-K time kernel above and must not import either
producer's amplitudes or energies. No time grid is inferred from observed
oscillations after running. Stream arrays in bounded blocks. Candidate total
budget is 900 wall seconds and 512 MiB peak resident memory for the complete
three-route experiment; exceeding it stops and archives an incomplete run.
There is no automatic K extension or tolerance repair.

Proposed scalar gates, all on E2 or the bounded amplitude G:

1. Each coarse/fine energy or corrected-work difference is at most 2e-10.
   Record `b_num=max(2e-10,4*|fine-coarse|)` for each producer, retaining
   the raw difference. This is a declared numerical allowance, not a proof.
2. At every source and K, fine finite-K energies from each independent pair
   agree within `2e-9+b_num,A+b_num,B`. Archive every gap and allowance.
   No averaging may conceal a failed pair.
3. At complex-amplitude sentinel momenta `(1/8,1/2,1,4,16,64)`, compare
   independently produced G and F with absolute error <=1e-9 at both
   observations. Post-pulse change in G must be <=1e-10. Check
   `|2k alpha1+iF0|<=1e-9` and the scaled first-order Wronskian condition
   `|Re(2k alpha1)|<=1e-10`. A vanishing/small amplitude does not receive a
   relative-error exemption.
4. The K256 spectral lower estimate must be strictly positive after its
   numerical allowance. Report the continuum enclosure obtained by adding
   the one-sided analytic UV tail; do not require a universal relative error.
5. E2,K must be nondecreasing in K within the pair's numerical allowances;
   all cross-cutoff gaps and actual tail ceilings must be shown. All finite-K
   spectral/mode/work comparisons use D_K, and the continuum statement uses D.
6. Offline mutation checks must reject the missing factor 1/2 in source work,
   the wrong frequency k instead of 2k, wrong Fourier phase/sign using the
   complex-amplitude gate, omission or reversal of D_K, substitution of D
   at finite K when resolvable, moving r, a nonretarded history, and a changed
   initial state. If a deliberately changed quantity is below its declared
   numerical discrimination budget, record that control as inconclusive;
   do not tighten a gate retrospectively to force rejection. Exact symbolic
   contacts and phase identities still must reject their corresponding
   algebraic mutations.

Save unrounded complex amplitudes, integrals before and after drift
subtraction, both refinement levels, all endpoint bounds, elapsed time,
source/input hashes, environment versions, start/finish timestamps and exit
codes. Preserve failed runs and source snapshots. A registered failure leads
to an explicit bounded repair proposal and a new freeze, never a silent
change to source width, r, state, subtraction, K ladder, or acceptance rule.

This design tests the quadratic coefficient directly. A nonlinear +/-epsilon
mode run with amplitude extrapolation would be useful separately, but is not
part of these candidate three routes and must not be added after registration.

## 7. Interpretation and progression

A complete pass would calibrate the positive canonical quadratic functional
and the background-dependent matched source-work contact for this prescribed
pulse and state. It would not establish a finite-amplitude yield, full
second-order minimal stress, damping of a dynamical shell, a bath, entropy
production, thermalization, or a stable coupled solution.

The next hierarchy remains explicit: obtain the full q2/rho2/p2 with coherence
if physical second-order stress is needed; derive responses at the actual
shifted endpoint and the metric kernels; then specify state-compatible
initial data, constraints, and the bulk/shell equations before coupled
evolution. A variance or canonical-energy number cannot be inserted into
bulk/gamma phenomenology as a substitute for those kernels. Apply a declared
species or loop multiplier only once. No pressure sign or equation of state
is inferred from the first-order pulse examples.

## 8. Literature and inherited scope

No new literature search or primary-paper reading was performed for this
proposal. It inherits the bounded review in the smooth-FRW checkpoint, pinned
in SOURCE_PINS.json. In particular, the reviewed PSAR analysis
[Ferreiro, Monin and Torrenti, 2311.08986](https://arxiv.org/abs/2311.08986)
motivates keeping the declared reference scales and matching all stress
terms, while the reviewed discussion by
[Armendariz-Picon and Diez-Tejedor, 2305.16293](https://arxiv.org/abs/2305.16293)
distinguishes occupation and coherent contributions to physical stress.
The inherited changing-mass Ward discussion in
[Molina-París, Anderson and Ramsey, gr-qc/9908037](https://arxiv.org/abs/gr-qc/9908037)
supports explicit source exchange, with conventions translated in the prior
review. Their models and approximations are not silently equated to this
one. Bogoliubov perturbation theory, causal work identities, Fourier
positivity and smooth-source tail estimates are established methods; this
proposal makes no novelty, completeness, or priority claim.
