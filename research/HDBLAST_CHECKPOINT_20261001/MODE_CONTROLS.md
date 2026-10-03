# Registered exact positive-pulse mode controls

## Scope

This tests a particle-mode integrator on a prescribed exactly soluble mass pulse. It does not derive a HDBLAST trajectory, account for the energy driving the pulse, thermalize particles, or include their gravitational backreaction. The solvable pulse was already in the archive; external mathematical novelty is not claimed.

REGISTRATION.md was committed at ceaeb64adf76b901ded8fa7d1bf1c913a29d22e4 before this round's first integration. The implementing agent reported its first execution at 2026-10-01T21:08:52.164Z. No implementation amendment was needed. The executed source and complete output are saved under code/ and outputs/.

## Equation and exact answer

In dimensionless time with pulse width one,

u'' + [omega_inf² + lambda(lambda+1) sech²(t)] u = 0.

Start in the asymptotic incoming positive-frequency vacuum. Extract outgoing Bogoliubov coefficients from u and u'. For this positive pulse the exact occupation is

n_exact = sin²(pi lambda)/sinh²(pi omega_inf).

Free lambda=0 and integer lambda=1,2 are reflectionless controls. The occupation is not monotonically increasing with pulse strength. This fact matters when screening a production prescription by coupling alone.

## Predetermined tests

lambda = 0, 0.5, 1, 1.25, 2; omega_inf = 0.5, 1, 2.
Time domains [-12,12] and [-16,16], with RK4 steps 0.02, 0.01, 0.005.

This gives 90 occupation checks and 30 finest-resolution normalization checks. The acceptance rule is |n-n_exact|<=max(2e-6,0.005 n_exact) for nonzero exact occupations, n<=1e-8 for exact-zero cases, and finest | |alpha|²-|beta|²-1 |<=1e-6. All full rows are retained.

The implementing agent's unchanged first source passed every registered check. Reported finest-resolution results:

| Quantity | Value |
|---|---:|
| Maximum occupation absolute error, either time domain | 2.4384440300e-10 |
| Maximum relative error for nonzero occupations, either domain | 3.6252797525e-8 |
| Maximum null-case occupation, either domain | 1.0801124143e-19 |
| Maximum Wronskian drift, either domain | 1.5013545962e-10 |
| Maximum fine occupation change, T=12 to T=16 | 2.2970839120e-10 |
| Maximum T=16 fine occupation absolute error | 1.4136011806e-11 |
| Maximum T=16 fine relative nonzero error | 4.8877721869e-10 |

The absolute tolerance is loose compared with the smallest omega_inf=2 occupations: pass/fail alone would allow approximately 14–29 percent relative errors there. The actual observed errors above are much smaller. Step-difference ratios are generally near fourth-order expectations for omega_inf<=1, but the tiny omega_inf=2 differences are affected by roundoff; no uniform fourth-order convergence claim is made for all rows.

## Exploratory signed work identity

For E=(|u'|²+Omega²|u|²)/2, the exact derivative is E'=(Omega²)'|u|²/2. The work integral is signed; it is not a positive source automatically available to radiation. Finite-endpoint pulse tails are included in E. The implementing source reports finest work residual at most 1.8558287306e-10.

This check is exploratory, as distinguished from the occupation/normalization registration. A future physical source must use renormalized stress and scalar current with consistent subtraction and an energy ledger. Passing an oscillator identity alone is insufficient.

## Reproduction and independent review

Run node code/mode_controls.mjs and compare the full JSON with outputs/mode_controls.json. The coefficient extraction and energy identity were separately checked algebraically. The root replay and any GitHub runner execution are recorded in EXECUTION_RECORD.md.
