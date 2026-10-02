# Prospective causal-response calibration

This registration must be published with its exact code and full SHA256 input
manifest before any new pulse response or forced-mode evaluation. Symbolic
derivation, analytic derivative-envelope design, syntax checks, and static code
review may precede publication. The parent checkpoint is the stationary quantum
package at commit `a109360b840bc1387221dd7d04d002176443a1a9`.

The immutable machine-readable choices are in `EXPERIMENT.json`. Geometry is
fixed de Sitter with H=1, r=2, xi=0, a=-1/eta. Incoming plane-wave BD initial
data remain fixed. eta_i=-6 and eta_f=-1.5. For u=eta+4, B=exp(1-1/(1-u²))
inside |u|<1 and zero outside. The prescribed s=a² delta x is epsilon B or
epsilon uB, epsilon=1e-4. Observations are exactly
[-5.5,-4.5,-4,-3.5,-2.5,-1.5]; finite comoving cutoffs are [64,128,256].
These amplitudes and grids cannot change in response to outputs. The geometry,
subtraction reference, mass law, and initial state cannot be fitted.

The common numerical comparison variable is y=a² delta Q/epsilon. Therefore
all response tolerances are absolute errors in this stated normalization,
including the independent forced-mode route. They are not relative errors
near zeros. The finite-K regulators are removed mathematically; K is not a
physical effective-theory cutoff.

## Primary numerical routes

1. Logarithmic memory:
   y=-(8pi²)^-1 integral f'(t)[ln(M(eta)(eta-t))+gamma_E+1]dt.
   For an observation during support, QUADPACK's `alg-logb` rule integrates
   the logarithmic endpoint. The remaining constant multiplies f(eta)
   analytically. Post-pulse integration has no singular endpoint. f=s/epsilon.
2. A separately implemented finite-part memory representation integrates
   [f(eta-tau)-f(eta)]/tau, with the exact endpoint derivative. Close source
   values use a rational exponent difference and expm1. Its local constant is
   ln(M(eta)(eta-eta_i))+gamma_E+1. This is an internal identity check, not a
   third independent physics derivation.
3. At each finite K, exchange the two finite momentum/time integrals and
   evaluate directly
   y_K=(8pi²)^-1{f0[asinh(K/M)-K/sqrt(K²+M²)]
   -integral f(eta-tau)[1-cos(2Ktau)]/tau dtau}.
   Use 2sin²(Ktau), panels no wider than pi/(2K), and sum with math.fsum.
   This is an independent finite-regulator numerical expression. It does not
   use the continuum subtraction formula to construct its answer.

Every primary quadrature requests epsabs=1e-13 in the normalized integral,
epsrel=1e-12, and QUADPACK limit=250. Finite-K epsabs is divided across the
precomputed panels. IntegrationWarnings and nonfinite values are fatal.
Individual normalized quadrature estimates must be <=1e-11; the two memory
representations must agree to 1e-10. A separately declared absolute 1e-10
allowance covers floating-point disagreement. No such estimate is a rigorous
quadrature enclosure.

The genuinely independent real-time mode implementation and its own frozen
manifest are in `independent/`. It solves the linear forced mode equation
with zero perturbation initial data using exact free propagation and local
Gauss quadrature. It must meet the separately registered coarse/fine 2e-10
gate and linear Wronskian/epsilon 1e-10 gate. Its finite-K normalized result
must agree with the primary finite-K result within 2e-9 plus the primary
quadrature estimate. The 2e-9 tolerance is an acceptance gate, not a theorem
that bounds mode integration error. The implementation retains raw mode
arrays and both refinement runs.

## Constructive combined tail and failure gates

The source is C-infinity with all initial derivatives zero. Repeated
integration by parts and the local subtraction yield

    |y-y_K| <= [3|f(eta)| M²/2 + A_normalized/4]/(16 pi² K²),
    A_normalized=|f''(eta)|+integral_{-6}^eta |f'''(t)|dt.

`theory/tail_bounds.py` computes the actual partial total variation, using
the exact critical-point polynomials, rational bisection for 256 steps,
60-digit directed interval arithmetic, and an outward-rounded binary64
upper bound on the complete expression. Details and rigorous fallback
ceilings are in the theory derivation. A_norm<=127 for positive_B and <=121
for signed_uB are independent coarse ceiling checks; the tight partial
values, not these fallback ceilings, enter the recorded bound.

The finite-K/continuum discrepancy must be at most this actual analytic tail
plus the two separately reported quadrature estimates and 1e-10 roundoff
allowance. At K=256 the normalized tail must be <=5e-6. Three-route continuum
agreement additionally allows the independently declared 2e-9 finite-K
cross-route gate. The tail theorem certifies only omitted momentum. The
complete response is not advertised as a rigorously certified enclosure.

Any failed gate stops this bounded calibration and preserves failed outputs.
Do not silently raise precision, cutoffs, budgets, tolerances, change state,
retune contacts, shift observations, or replace an implementation. A later
attempt requires a public prospective amendment. Each producer has a fixed
900-second wall-time budget. Dependencies are pinned in requirements.txt.

## Exact references and wrong-formula controls

Before support both sources give exactly zero. For positive_B after support,
the response is negative with |y|>0.004 at both selected observations; this
follows from evenness, a fixed lower area bound, and eta+4<=2.5. The signed
uB source also has strictly negative post-pulse response since its paired
integrand is 2u²B/(D²-u²)>0 for D=eta+4>1; no 0.004 margin is imposed on it.
Sign reversal, an instantaneous-only response, or a window excluding the
pulse must fail. An advanced and a time-symmetric response to a future source
must fail the exact pre-pulse zero.

The past-infinite stationary source is an independent analytic identity,
not an additional finite compact-pulse run:
Q_x=-(2gamma_E+ln2)/(16pi²). It must agree with
[psi(2)+psi(1)-ln2-1]/(16pi²) to 1e-14. Omission of +1, omission of gamma_E,
and a moving subtraction reference are rejected with finite >0.001 margins.
An exact separated-time kernel reference rejects reversed Kubo sign, missing
Wick factor two, and replacing sin(2kDelta) by sin(kDelta).

Two deliberately different state diagnostics use n(k)=0.01B((k-1)/0.25),
or a real beta(k) of the same shape, at Delta=0.1 and0.4. The occupation
changes the normalized kernel by
-(2pi²)^-1 integral n(k)sin(2kDelta)dk<0, with magnitude>1e-6. The initial
beta changes a² deltaQ by
(2pi²)^-1 integral k beta(k)cos(2kDelta)dk>1e-5. Compact support is
0.75<k<1.25. The sine/cosine have definite signs on both supports. A separate
mpmath implementation in the validator checks the QUADPACK evaluations.
Silently returning a vacuum kernel or omitting the initial-state boundary
term is rejected. These states are not used in the primary response.

Omitting the input a(t)² gives the wrong forcing t²f(t), independently
integrated as a mutation. Omitting the output a(eta)^-2 is a separate mutation.
Theory's frozen symbolic checks additionally verify dimensions, coordinate
rescaling, quadratic mass-law contact, reference matching, and control
identities. All checks use explicit exceptions and must pass under Python -O.

## Provenance and replay

The producer requires `--registration-sha256`, `--public-freeze-commit`, and
`--output`. It checks every FULL_REGISTRATION.json input and refuses to write
inside the frozen package or overwrite an existing output. It records package
versions, command, exact producer/experiment hashes, start time, elapsed time,
and failures. Source imports perform no numerical source or mode evaluations.
An active point and completed rows are written incrementally to
`partial_results.json`; completed finite-K calculations at the active point
are retained. The wall budget is checked before and after every bounded
quadrature and observation, as well as at completion.

The standalone validator imports no producer functions. It requires the same
registration pin plus explicit primary/mode result files and a fresh output
path. It binds the independent manifest and file list to the full registration,
requires all four coarse/fine source runs and their passed status, checks their
global resource/Wronskian limits, verifies raw archive hashes within the fresh
output directory, and checks the summary against both recorded raw runs.
Run it once normally and once under `python -O` with different output
files. Preserve all logs and exit status. Reproduction must use the frozen
inputs even if a later diagnostic corrects a defect.

This test calibrates one retarded scalar susceptibility on a prescribed
geometry. It does not derive a stress-current response, evolve a shell,
establish quantum stability or relaxation, calculate particles, radiation,
thermalization, or heating, or constitute a claim of discovery.
