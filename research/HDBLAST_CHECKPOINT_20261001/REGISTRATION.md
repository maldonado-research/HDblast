# Validation protocol and exploratory disclosure — 1 October 2026

Ricardo Maldonado's HDBLAST research program. These are model-internal analytical diagnostics and numerical controls, not observational evidence or external peer review.

## Baseline and scope
Private baseline commit: 8f67197b730d4e1c43554b86f224c29cc72629eb.
Public baseline commit: 6b3aa166fa73b4d8f13071f788b74cf38b2946a0.
The September 30 registered outcomes remain unchanged. No new coupled 5D evolution, physical reheating mechanism, thermalisation result, interval certification or external mathematical novelty is claimed.

## Disclosure before this file was committed
The previous checkpoint and audits were reviewed. Candidate single-cohort energy-envelope and frozen-shell duration formulae were derived exploratorily, and several approximate endpoint numbers were already calculated: conditional r lower budgets about 51, 160, 1327 and 4183; the Y=2 negative-vacuum share is about 0.285 of linear radiation; its quadratic-density correction is small. Independent reviewers had begun checking these proposals. These are retrospective/exploratory results, not preregistered discovery tests.

The archived positive sech-squared mass-pulse production formula and its integer-coupling reflectionless cases were known before this registration. The new mode-integration validation below has not been executed by this research round at registration time.

## C1: paid particle-production budget
Specialise to one real scalar degree of freedom, or explicitly retain occupation/degeneracy factors. State any pre-existing radiation separately. Bound the combined surviving-particle plus decay-radiation energy for monotone expansion, bounded particle mass, nonnegative one-way decay, and at most the specified number of produced particles. Multiple crossing cohorts require separately summed bounds. Recollapse and unbounded mass are outside scope.
Check Gaussian number and momentum moments by independent quadrature; check finite-bin decay histories and Jensen inequality. Preserve population corrections and the distinction between a fixed recorded history and a new backreacted trajectory. Label all applications conditional and post hoc.
Identical-species scaling and any species-cutoff prescription must be stated separately. A brane-localised species cutoff is not established here. Do not claim a universal species no-go or a realised physical transfer rate.

## C2: radiation dominance and duration
Use an absolute-energy diagnostic that includes Weyl, residual vacuum, quadratic density and any unresolved scalar contribution. It is a new diagnostic, not a replacement for B3's registered verdict.
Frozen-source-free continuation requires fixed positive shell tension and negligible or explicitly controlled scalar evolution. Verify the exact quadratic-density turnaround root against the Friedmann identity and a separately implemented numerical evolution, including a deliberately large quadratic-density control. The already-present quadratic term in B4's differential equation is not a newly discovered omission; only its simplified closed form omits it.
Apply the radiation Ward identity during rolling using endpoint quantities, without assuming nonzero j*v on a strictly frozen scalar.
No entire-trajectory duration claim follows from one saved plateau. Full trajectory reclassification and coupled source/gauge evolution remain future tests.

## C3: prospective mode-solver controls
Equation in dimensionless time: u'' + [omega_inf^2 + lambda*(lambda+1)*sech(t)^2]u = 0.
Incoming vacuum is imposed in the asymptotic past. Extract outgoing Bogoliubov coefficients from u and u'; check the Wronskian and the known exact occupation sin(pi*lambda)^2/sinh(pi*omega_inf)^2.
Predetermined cases: lambda = 0, 0.5, 1, 1.25, 2 and omega_inf = 0.5, 1, 2 (15 cases). Include lambda=0 free-vacuum and integer reflectionless controls.
Use at least three time steps, halving each time; retain the full table. Require nonzero occupations with exact value >=1e-8 to agree within max(2e-6 absolute, 0.005 relative), and exact-zero cases to have occupation <=1e-8. Require Wronskian |alpha|^2-|beta|^2 error <=1e-6 at the finest level. Check finite-domain sensitivity independently by enlarging the time window.
Failed controls must be reported; an implementation change requires a dated amendment before its new outcomes are used. Passing this prescribed-history benchmark validates the control only, not a self-consistent 5D particle-production solution.

## Evidence and review
Save executable scripts, input provenance, complete numerical outputs and a separate skeptical AUDIT.md. Report execution environment accurately. New physical mechanisms require action-derived junctions and renormalised stress/current with consistent energy ledgers and cutoff checks.

## Publishing
Keep raw/private archive material private. Public deliverables are newly prepared scientific summaries and reproducible controls, with conditional labels and links to the existing public baseline. Verify Zenodo concept membership before any versioning; do not guess a successor label.
