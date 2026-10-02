# Prospective fixed-geometry matched-stress calibration

Baseline: public causal checkpoint e784b18128012825f126c5bd96a4bf2033ab63d9.
This registration explicitly reuses its two compact pulses, six observation
times and K=64,128,256. The new stress calculation starts only after the
complete implementation, numerical choices and checks are publicly frozen.
No new mode, memory-response or stress experiment has run before the freeze.
Pure algebra, exact critical-root isolation and budget design may precede it;
deferred pulse intervals are evaluated in the registered run.

## Model and observables

H=1, r=2, xi=0; a=-1/eta; incoming BD state and geometry fixed.
With u=eta+4, B=exp(1-1/(1-u²)) for |u|<1 and zero otherwise,
s=a²delta_x=epsilon f, epsilon=1e-4, f=B or uB.
Initial eta=-6, final eta=-1.5. All source derivatives vanish outside
-5<eta<-3. These are formal perturbative inputs, not measured parameters.
The original quadratic mass-law translation is evaluated separately at
b=1, delta_x=br delta_phi.

Outputs: q=a²deltaQ/epsilon and q',q'', R=a⁴delta_rho/epsilon,
P=a⁴delta_p/epsilon, physical Q0,K, anomaly a⁴deltaA_K/epsilon and
current a²delta_j/epsilon. Stress is the minimally coupled physical
stress. Fixed-r second-order variance and complete fourth-order stress
subtraction inherit one common finite prescription; no coefficient is fitted.

## Independent calculations

The primary evaluates regular logarithmic memory integrals of the first
three source derivatives and exact finite-band sine-history expressions,
retaining scale-derivative contacts. Both stresses use direct subtraction
reductions; the pressure contact comes from its own integrated W4 mismatch.
Neither stress is defined through trace or Ward.

The independent producer evolves complex forced variations with exact
local free propagation and eight-node local forcing. Evolution, physical
bilinears, complete W2/W4 subtraction and sums use longdouble/clongdouble
with at least63mantissa bits. No Wronskian projection is permitted.
Both stresses are calculated directly from the actual complex modes and
explicit mass-operator contacts. True u,w equations supply variance
derivatives. Mode and counterterms are combined before integration.
Coarse/fine settings refine time and momentum together: empirical combined
refinement, not a certified error enclosure.

Direct density/pressure are also archived at every time step. A separate
composite-Simpson integral of
a²Q0,K(f'-2Lf)/2+LR-3LP is compared with the independently calculated
endpoint R. Density is never defined through this ledger. The grids align
all observations with even Simpson endpoints. Q0,K is time dependent at
fixed comoving K; finite-band trace uses its actual anomaly and baseline.

## Errors and coverage

EXPERIMENT.json supplies the numerical settings, gates and resource stops.
Core coverage is36finite-band points timesfiveqjet/stress quantities:
180comparisons. Baseline, anomaly, current, trace, Ward, refinement and
Wronskian records remain distinct. Validators use explicit checks normally
and under optimization. Mutations identify algebraic, raw-state or numerical
scope; a small finite-contact mutation is not presented as necessarily
failing a looser main gate. An input diagnostic is not a physical
alternate-state experiment.

Omitted-band bounds use separate budgets
|f^(j+2)|+integral|f^(j+3)| for j=0,1,2, exact critical-root total
variations and directed60decimal interval arithmetic. Local tail differences
and moving-scale contacts are retained before propagation to R and P.
These bounds cover omitted momentum only. QUADPACK estimates, grid spread
and arithmetic allowances remain separate empirical evidence. No certified
total error or continuum stress-derivative/Ward residual is claimed.
Actual continuum gaps and bound sizes are reported even when accepted.

Pre-pulse response and unchanged initial data must pass. Analytic theory
predicts negative relative linear density after a nonzero positive compact
pulse; this coherent polarization diagnostic is distinct from total
reference energy and quadratic excitation energy. Pressure has no universal
post-pulse sign.

## Stops, provenance and scope

Publish FULL_REGISTRATION.json and all pins, verify the remote commit,
then run with those public pins. Output paths must be new and external.
Record actual source, manifest, registration, runtime, resources and raw
archive hashes; retain partial progress, exceptions and failures.
No silent change of thresholds, grids, sources, finite contacts or
implementations after outcomes. A correction needs a separate explicit
published amendment and preservation of the original failure.

This is a fixed-geometry homogeneous linear scalar-to-stress/current
response at x=r=2H². The actual shifted-root propagator, metric responses,
bulk/shell boundaries, matching and quantum initial data remain necessary
for coupled evolution or stability. No first-order stress result establishes
net heating, particle yield, radiation transfer, thermalization, a hot
Big Bang, extra dimensions or observational support for HDBLAST.
