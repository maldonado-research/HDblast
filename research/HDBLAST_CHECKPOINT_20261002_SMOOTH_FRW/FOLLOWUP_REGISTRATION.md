# Separate prospective follow-up: curved pressure at K384

This protocol follows the completed original matrix, whose overall status
remains FAIL. The original A=0 control passed. The original A=0.2 control failed
the pressure cutoff gate: K96-to-K192 maximum change 0.02262428337688216 exceeds
its unchanged allowance 0.011232728419574517. Its other final gates passed.
That outcome is permanently retained; this is a new numerical experiment.

Commit this protocol and its wrapper before executing any new K384 modes.
Do not alter the frozen repaired core, physical model or earlier results.
The registered maximum here is K384. There is no conditional K768 extension,
no assumed asymptotic power, and no extrapolated acceptance.

## Frozen inputs and implementation

Use only A=0.2,r=4,T=1 with the identical C-infinity a,x profiles, exact incoming
vacuum, complete stress/current subtractions, finite-action convention and
eta interval [0,1.25] from the original registration. The core source is
smooth_frw_control.py SHA256
32f6632bb84cc10dd6cef22c2311941acdbe223d846eb7157a2f1e9459fb160a.

The new wrapper imports this exact source and changes only its in-process
cutoff bookkeeping tuple to (24,48,96,192,384), enabling the existing run_case
implementation to integrate and record K384. The core file remains unchanged.
Import hash, follow-up protocol hash, wrapper hash, original summary hash and
all raw arrays are recorded. The baseline is a lossless selected-column export
of the complete preserved run registered-57a97a2-primary-002, shipped in
followup_baseline. Its arrays.npz files contain exactly t,K192_rho,K192_p,K192_Q
copied unchanged from the original full arrays. BASELINE_MANIFEST.json records
both full-original and minimal-consumed file hashes, sizes and exact selected
column equality. This compact baseline permits reproducibility without
requiring the original full per-mode arrays. All NPZ loads use allow_pickle=False.
The wrapper enforces these consumed artifact hashes before any new mode run:

| Original artifact | SHA256 |
|---|---|
| summary.json | f6cc20865e61d8134c63ecc8a9fa6bf1695d9ffafd0eea0890ca9aeb471e0228 |
| A0.2_tight_K192/arrays.npz | 546c09231f504da70f3a28cb598fe40559d56d26af4b0e7b01a998926fc3ed0d |
| A0.2_quadrature_K192/arrays.npz | 62cc23be671106b228f3bb2cb99870b4d15df41d1e00676ed7e35c9c5cce0994 |

The baseline and new time arrays must match exactly before comparing curves.

## Fixed new matrix

At K384 run the following curved controls, all with 801 uniform output times:

1. Tight: eight Gauss-Legendre nodes per width-two panel, rtol=2e-13,
   atol=2e-15, phase_step=0.1.
2. Quadrature refinement: twelve nodes/panel, same solver tolerances and
   phase_step=0.1.
3. Solver refinement: eight nodes/panel, same tolerances, phase_step=0.05.

The core's max_step formula remains min(1/80,
phase_step/sqrt(kmax^2+exp(2A)r)). No mode normalization is enforced.

Also execute an exact-static B identically zero,A=0,K384 negative control
with eight nodes/panel and phase_step=0.05 at the same tight tolerances.
Its exact subtracted rho,p,Q are zero. This measures the new larger-cutoff
floating-point floor rather than inferring it from K96.

## Unchanged gates and comparisons

Compare the tight solve's K192 and K384 integrals separately for rho,p,Q.
Compare eight versus twelve nodes at K384 and phase_step=0.1 versus 0.05 at
K384. Compare each new tight/quadrature solve's K192 subset against its
corresponding preserved original K192 solve. Every observable comparison uses
the same maximum-over-time rule as before:

max|f1-f2| <= 2e-5+0.005 max(max|f1|,max|f2|).

The selected curves are from twelve-node K384. All of the original paired-work
ledger, five-point local exchange, finite-K direct trace and sampled trace,
future spectral energy, and 401/801 time-quadrature gates apply unchanged:

- Every new raw normalized Wronskian maximum <=2e-10.
- Integrated energy-minus-work residual <=1e-5+0.005 maximum ledger/energy
  change scale. Use both Simpson paired renormalized work and the independent
  ODE-accumulated bare work minus the analytic subtraction energy change;
  retain the distinction that subtraction cancels in the latter diagnostic.
- Nested time-grid ledger comparison <=1e-5+0.005 maximum ledger scale.
- Local rho'+3H(rho+p)-x'Q/2 residual <=2e-4+0.01 maximum individual term.
- Finite-K trace identity with direct physical-equation Q derivatives and
  analytic subtraction derivatives <=2e-6+5e-6 maximum side.
- Finite-K trace identity with 801 sampled Q derivatives <=2e-4+0.01 maximum
  side; report the 401-node result too.
- Static-future full rho versus exact beta occupation energy uses the ordinary
  2e-5+0.005 maximum scale allowance; pressure/Q coherence stays included.
- New exact-static control maximum |rho|,|p|,|Q| <=2e-5.
- All stored arrays and strict JSON values finite; preserve each run, failed
  outputs and full raw mode arrays.

The follow-up passes only if every listed gate passes. No gate, denominator,
sampling set, cutoff, tolerance or phase-step setting may be revised after
viewing its results. A failure remains a failure. The independent Radau mode
review and algebraic finite-action/trace review from the original experiment
remain applicable to the unchanged physical solver and subtraction code;
new high-k numerical reliability is tested by the explicit K384 solver,
quadrature, Wronskian and static controls here.

## Interpretation

A pass provides finite-cutoff numerical convergence evidence for this smooth
prescribed curved control at the existing tolerances. It is not a proof of
the infinite-cutoff limit. It does not change the original matrix's FAIL
classification, the archived shell, or establish backreaction, source
consumption, decay, thermal radiation or a new cosmological solution.
