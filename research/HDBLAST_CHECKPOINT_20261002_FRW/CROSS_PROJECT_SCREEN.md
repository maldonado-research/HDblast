# Transferable methods from the author's other research

Five public repositories were screened read-only at the commits below. All READMEs and selected technical notes/code were inspected. No screened project supplies a completed curved-brane scalar renormalized stress tensor or common-action backreaction solver. Their physical theories and parameter values are not being identified with HDBLAST.

| Repository | Screened commit | Useful method and limit |
|---|---|---|
| Unified-Theory-of-Everything | cb4fef5f93a03c2cb110d91feb0503be08715721 | Independent evolution in changed variables with implicit Radau and analytic Jacobian; explicit preparation, endpoint-basis, solver and UV-tail error separation. Fermion SU(2)/Pauli bounds cannot be copied to scalar SU(1,1). |
| dmde-research | 4e7274cb3ce5256433406b4b192876bf41ec740b | Test that the actual evolution consumes its source; refine inside the active time/momentum windows rather than just increasing total grid size. Its distributions/thresholds are not HDBLAST inputs. |
| AntiMatter | 4be30eb0bb17ed7c11a44525c53929ce5d64985e | Include the driver reservoir and complete potential drop in the energy ledger; particle production alone does not establish an autonomous solution. No model parameters transferred. |
| black-holes-hrf | 04436742ebc8c20f7e3ff201f1d9c98e8761a25f | Nuisance projection and separation limits may help future observable identifiability; no present curved-source implementation found. |
| TMD | 5f8be36b22b4606cbf58fe3bb466dfdb9a6ab70d | Count-data confidence methods require assumptions absent in deterministic mode integration; no immediate transfer adopted. |

The independently written four-real-component Radau smooth-step control in this checkpoint concretely applies the first method. Its scalar equation, normalization and analytic occupation are separately derived; no fermionic transition formula is imported.

Key sources:
- [TOE independent reviewer](https://github.com/maldonado-research/Unified-Theory-of-Everything/blob/cb4fef5f93a03c2cb110d91feb0503be08715721/checkpoint_2026_09_06/fermion_portal/independent_review.py) and [dynamic EFT review](https://github.com/maldonado-research/Unified-Theory-of-Everything/blob/cb4fef5f93a03c2cb110d91feb0503be08715721/checkpoint_2026_09_06/dynamic_eft_review.md).
- [TOE off-shell vertex review](https://github.com/maldonado-research/Unified-Theory-of-Everything/blob/cb4fef5f93a03c2cb110d91feb0503be08715721/checkpoint_2026_09_06/routed_loop/SEED_VERTEX_REVIEW.md): preserve common-action derivative terms and correlated matching choices.
- [DMDE source-consumption gate](https://github.com/maldonado-research/dmde-research/blob/4e7274cb3ce5256433406b4b192876bf41ec740b/provider/docs/DMDE_v0920_SOURCE_OPERATOR_CLOSURE_GATE.md) and [meaningful refinement implementation](https://github.com/maldonado-research/dmde-research/blob/4e7274cb3ce5256433406b4b192876bf41ec740b/provider/code/dmde_v0920_meaningful_refinement.py).
- [AntiMatter energy frontier](https://github.com/maldonado-research/AntiMatter/blob/4be30eb0bb17ed7c11a44525c53929ce5d64985e/research/AM1231/v1.23.1/v1.23.1_efficiency_energy_frontier.md).
- [HRF theorem](https://github.com/maldonado-research/black-holes-hrf/blob/04436742ebc8c20f7e3ff201f1d9c98e8761a25f/research/v72/GW250114_V72_THEOREM_AND_PROOF.md).
- [TMD extension](https://github.com/maldonado-research/TMD/blob/5f8be36b22b4606cbf58fe3bb466dfdb9a6ab70d/TMD_RESEARCH_EXTENSION_0_5_0.md).

A connected Drive DBlast 3 vector-store folder and a September 8 master handoff were also inspected. That handoff covers historical PTA covariance certificates, brane mappings and bulk-blast reheating questions; it supplies no replacement for the current archived shell's local high-order derivatives or state completion. It was used as historical context, not as new validated numerical input. Raw handoffs and private correspondence are not reproduced here.
