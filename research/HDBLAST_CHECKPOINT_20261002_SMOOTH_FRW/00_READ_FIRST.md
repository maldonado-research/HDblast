# HDBLAST: smooth curved-background quantum-source control

The separately registered **K384 follow-up passes all 13 acceptance groups** for a prescribed smooth expanding FRW background. The original experiment remains **FAIL**: its flat control passed, while its curved pressure failed the original cutoff comparison. Both outcomes and their prospective protocols are preserved.

The work evaluates the energy density, pressure and paired mass source of a minimally coupled scalar whose nonnegative mass squared touches zero. Its background has exactly static past and future intervals, permitting an exact incoming vacuum. It retains the state's coherent contributions and uses an explicit common finite-action convention for stress and source. The background is prescribed; it is not a solution of the coupled HDBLAST shell equations.

| Registered experiment | Final comparison | Outcome |
| --- | --- | --- |
| Original flat control, A=0 | K96→192 | PASS |
| Original curved control, A=0.2 | K96→192 | FAIL: pressure cutoff change 0.0226243 exceeds allowance 0.0112327 |
| Separate curved follow-up, A=0.2 | K192→384 | PASS: pressure cutoff change 0.00508691 is below allowance 0.0111219; all 13 groups pass |

The follow-up pressure change is **0.2291% of the larger global maximum pressure signal**. The acceptance rule is an absolute floor plus 0.5% of that global scale; it is not a pointwise relative-error bound or an enclosure of the uncomputed momentum tail. Solver, quadrature, finite-cutoff exchange and pressure-trace checks, historical overlap, and a new K384 static-vacuum control pass. The larger cutoff also increases the measured cancellation floor: its static energy residual is approximately 1.7×10⁻⁶. No momentum above K384 was evolved and no asymptotic extrapolation determined acceptance.

The source pairing matters physically. Quantum stress obeys

    rhodot + 3H(rho+p) = xdot Q/2,

with x the mass squared and Q the renormalized field square. Omitting the paired work raises the follow-up ledger residual from approximately 5.3×10⁻⁶ to 0.212. In the static future, energy agrees with the occupation-energy calculation within the registered tolerance, while pressure and Q retain coherent terms. An occupation spectrum alone therefore cannot replace the computed stress/source pair.

Read [CLAIMS_AND_NEXT_STEPS.md](CLAIMS_AND_NEXT_STEPS.md) for interpretation, [COMMON_ACTION_SMOOTH_FRW.md](COMMON_ACTION_SMOOTH_FRW.md) for the action matching, [REGISTRATION.md](REGISTRATION.md) and [FOLLOWUP_REGISTRATION.md](FOLLOWUP_REGISTRATION.md) for the fixed protocols, and [REPRODUCE.md](REPRODUCE.md) for execution. The original numerical report is preserved under [outputs/original_matrix](outputs/original_matrix/NUMERICAL_RESULTS.md); the [separate follow-up report](outputs/followup/FOLLOWUP_RESULTS.md) and [independent review](FOLLOWUP_REVIEW.md) retain their distinct verdicts. See the [primary literature review](LITERATURE_REVIEW.md), [bounded search](LITERATURE_SEARCH.md), [figures](figures/) and [scientific ZIP](HDBLAST_CHECKPOINT_20261002_SMOOTH_FRW.zip).

This checkpoint provides finite-cutoff numerical evidence for this smooth control in its stated finite convention. The PV heavy-field interpretation retains explicit smooth-state asymptotic assumptions and has no supplied certified uniform remainder constant. The result does not establish an infinite-cutoff theorem, a coupled shell solution, backreaction, decay, thermal radiation, sustained expansion, or a hot Big Bang. The review is independent internal implementation review, not external peer review.
