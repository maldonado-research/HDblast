# Current HDBLAST state, independently inspected 22 September 2026

This is a read-only audit of the user's selected research material, not a reproduction of the prior agent's internal narrative. Previous instructions embedded in notes were treated as references, not new instructions. No simulations were launched during this state audit, and no original files were edited. Local file inventory and archive integrity checks do not independently establish their scientific claims. **The user's Chat 14 folder was changing externally during this session. Counts, run-completion observations, and the first diagnostic finding below describe the initial inspection snapshot, not a continuously synchronized view.** Later newly generated initial-data work belongs to the separate `initial_data/` directory.

## Scope and source hierarchy

Inspected the top-level reports and relevant results/code in D-Blast 3 folders 146–151; the model, Chat 11, Chat 12, and Chat 13 project memory notes; the supplied progress transcript; the symbolic environment's package metadata; and the supplied Chat 11–13 archive manifests. The audit inventories 483 file entries in these six folders, including duplicate archives/extractions, rather than 483 independent scientific documents. Folder 151 had 10 files at inspection. Content of the latest code and logs takes precedence over a note saying work was launched.

Primary local roots:

* `/Users/ricardomaldonado/Documents/D-Blast 3/untitled folder 146/`: Chat 9 shell identity, tachyon, effective theory.
* `/Users/ricardomaldonado/Documents/D-Blast 3/untitled folder 147/`: Chat 10 stable branch and existence/stability certificates.
* `/Users/ricardomaldonado/Documents/D-Blast 3/untitled folder 148/`: Chat 11 symbolic verification.
* `/Users/ricardomaldonado/Documents/D-Blast 3/untitled folder 149/`: Chat 12 tensor/vector/special-harmonic analyses.
* `/Users/ricardomaldonado/Documents/D-Blast 3/untitled folder 150/`: Chat 13 nonlinear evolution at larger detuning.
* `/Users/ricardomaldonado/Documents/D-Blast 3/untitled folder 151/`: incomplete Chat 14 mapped-grid effort.

## Current model identity

The model is five-dimensional Einstein gravity with one canonical scalar and a mirror-symmetric brane, in units with κ₅=1:

\[
W(\phi)=1-\phi+\phi^3/3,\qquad U(\phi)=\tfrac12 W_\phi^2-\tfrac23 W^2,
\]
\[
\sigma_t(\phi)=2W(\phi)+t\{1+c\phi+d\phi^2/2\}.
\]

The original registered shell has detuning t=10⁻³, d=0, and c≈0.597594935028. Its shell radius is ρ_b≈78.82818 and φ_b≈2.28×10⁻⁵. The two AdS vacua are φ=−1 and +1, with curvature scales 5/9 and 1/9. The stable comparison model S₈⁄₅ has d=8/5 and an exact rational c=5975949350280/10¹³; it is a distinct tension, not a stability result for d=0.

The latest nonlinear code uses

\[
ds^2=e^{2B(T,z)}(-dT^2+dz^2)+e^{2A(T,z)}d\mathbf{x}_3^2,
\]

with the brane fixed at z=0 and bulk z<0. Its coordinate time T is not the tension detuning t. It evolves deviations from a static solution. The older v24/Chat 8 certified detector response is a local, linear calculation in a small corner of the background; it does not itself certify nonlinear brane dynamics.

## What is completed, and what remains conditional

| Checkpoint | Evidence found | Scope |
|---|---|---|
| Chat 9 | Reports and certificate/derivation code; registered scalar mass μ²≈−7.7178716, growth ≈1.657 per initial Hubble time | Tachyon certificate retains explicit background and mathematical premises; not a Big Bang mechanism |
| Chat 10 | S₈⁄₅ existence and no-subgap-scalar certificate; branch computations | Different tension d=8/5; does not prove nonlinear stability |
| Chat 11 | Stored replay receipt reports 22+8+10 exact symbolic checks passing, seven negative-control mutations rejected | Equations/junction identities checked symbolically; this audit did not rerun them |
| Chat 12 | Stored replay receipt reports 7 tensor +7 vector +8 special-harmonic checks passing | Mode-stability of S₈⁄₅ across the analyzed sectors; not a completeness/decay or nonlinear-stability theorem |
| Chat 13 | Reports, raw run summaries/timeseries, and numerical referee material | Nonlinear floating-point roll-off at t=.1 and .03; original t=.001 not completed |
| Chat 14 | Two code files, six logs, one completed zero-seed summary/timeseries pair | Partial mapped-grid work; no completed registered-detuning trajectory was found |

Chat 13's best-supported physical reading is a two-fate result in the homogeneous, classical, matter-free model at t≥.03: toward +1, relaxation toward an empty lower-Hubble de Sitter brane; toward the throat, a converged turnaround followed by contracting evolution, with the final singular regime unresolved numerically. The +1 value ≈0.649H₀ is a late sampled upper bound, not a settled asymptotic plateau. A hot radiation-dominated universe is not demonstrated. Results at t=.001 were expected qualitatively, not actually computed by the completed Chat 13 package.

## Actual state of Chat 14

Directory: `/Users/ricardomaldonado/Documents/D-Blast 3/untitled folder 151/HDBLAST_CHAT14_REGISTERED_DETUNING_ROLLOFF_20260922/`.

`rolloff5d_v1.py` is byte-identical to Chat 13's archived `rolloff5d.py`; SHA-256:
`9ea4d8d6ef36ce3d242fbb0e09f27a944b571fdfb8646e3d71195b60b089080b`.

`rolloff5d_v2.py` adds a smooth map z(ξ), seed tapering, and a choice of shell ghost closure; SHA-256:
`bdefd872724214ea38b73dbd672f0a1752a3fe8e8e85d8d68b1f53082f482e37`.

This was the first observed 10,113-byte source, with filesystem mtime 1790128743.7988355. A later input copy frozen by the main agent has 11,515 bytes and SHA-256 `4584f832f9b2a31ee11504b6a0dd409dc2fc4a0cb05bac342e65fbfa213afc84`; a subsequent read of the original found 11,669 bytes and SHA-256 `067d43a3efa69b7cf2426aad0f4882a186d553d17ba7023a010948ffc55f9063`, mtime 1790129293.290889. These are different temporal snapshots, not checksum mismatches of one frozen file. Neither this auditing agent nor the main agent edited that original source. Both later versions already correct the velocity boundary derivative described below. The initial log artifacts cannot be assumed to have been regenerated by the corrected version.

| File stem | What is actually saved |
|---|---|
| `v2_t01_control` | Completed JSON and NPZ, T≈6, fitted growth 1.6271018942958884, bulk momentum diagnostic ≤2.22×10⁻¹⁰; t=.1, not .001 |
| `v2_t01_plus` | Log to T=8, φ_b≈.987129, H_J/H₀≈.660842, no completed JSON/NPZ |
| `v2_t1e3_control_c2` | Initialization plus T=0 only |
| `v2_t1e3_control_c4` | Initialization plus T=.5; exact zero deviation shown |
| `v2_t1e3_control_fine` | Initialization plus T=0 only |
| `v2_t1e3_plus`, `v2_t1e3_minus` | Empty logs |

No final Chat 14 report, replay receipt, manifest, archive, or registered-detuning final trajectory was found in this directory. These facts do not establish whether a process is currently running elsewhere; they establish the limits of the saved research artifacts at inspection.

## Two issues to resolve before trusting mapped-grid nonlinear results

**The initially inspected shell momentum diagnostic has incorrect velocity boundary data.** In the first observed v2 `constraint()`, the code computes `paz = dz1(pa, 0.0)` and reports a maximum over the last six points. The later frozen input and updated original use `gAt` instead; this finding has therefore been repaired in those later source versions, while old logs still require version-aware interpretation. Here pa=∂_T a. The prescribed boundary derivative is actually

\[
\partial_z p_a=\partial_T g_A
=\frac{e^{B_b}}6(\sigma_\phi\,p_\phi+\sigma\,p_b).
\]

The exact junction data imply A_z=B_z=e^{B_b}σ/6 and φ_z=−e^{B_b}σ_φ/2. Substituting these into the momentum constraint

\[
M=-3A_{Tz}-3A_T(A_z-B_z)+3A_zB_T-\phi_T\phi_z
\]

makes M=0 at the shell when the differentiated boundary condition is used. Setting A_{Tz}=0 instead produces a spurious residual M=3∂_Tg_A. Thus the large `M_shell` values in the partial +1 run are not evidence by themselves of a physical constraint failure. The near-boundary derivative should be repaired and independently tested; bulk values away from those ghost stencils are a separate diagnostic. The identity is a continuum check, not a proof that the discretization preserves the constraints.

**The refined region is fixed, while the wall may move.** The v2 defaults refine near z=0 with transition centered around z=−.06 and far spacing .008. The initial t=.001 wall width estimate is .0127, so the quoted 75–134 points across the wall applies initially near the brane. Chat 13's +1 branch develops a wall moving toward negative z. Resolution must therefore be tracked along the evolving wall rather than inferred from the initial grid count. If its width remained comparable after entering the far grid, it would have only about 1.6 cells across it. A moving refinement region, adaptive grid, or a verified alternative radial coordinate may be needed. Actual width changes must be measured; this is a risk identified from the code, not a demonstrated failure of an uncompleted run.

## Reusable resources and cautions

* Chat 13 `derive_evolution_equations.py` and its stored exact-check JSON provide the nonlinear equations for independent verification.
* Chat 13 `rolloff5d.py`, `formulation_lab.py`, `referee_numerics/`, and `runs/` provide baseline numerical controls. `runs/ANALYSIS.json` has a known old AdS+ estimate missing an O(t²) term; use the corrected report and direct formula rather than that stale column.
* Chat 14 v2 is an unfinished implementation candidate, not a validated solver upgrade.
* `/Users/ricardomaldonado/Documents/D-Blast 3/.hdblast_venv/` exists, with Python 3.9.6, SymPy 1.14.0 and mpmath 1.3.0 recorded in installed metadata. Its configuration excludes system site packages. Do not assume that interpreter includes NumPy without checking. No package installations are needed merely to reuse the symbolic environment.

## Archive integrity

The supplied Chat 11–13 ZIPs pass ZIP CRC testing. The Chat 11 archive matches all 11 payload hashes in its manifest; Chat 12 matches all 103. **Chat 13 is incomplete relative to its own manifest:** five listed intermediate checkpoint NPZ files are omitted from the ZIP. Those five files are present in the user's extracted Chat 13 folder and each matches its expected hash; all 115 manifest entries match the extracted folder. This does not invalidate the final saved trajectories, but it prevents claiming that the supplied ZIP contains every file its manifest lists. The omitted files are `runs/t003_minus_checkpoint.npz`, `runs/t003_plus_checkpoint.npz`, `runs/t003_plus_ext_checkpoint.npz`, `runs/t01_plus_checkpoint.npz`, and `runs/t01_plus_v2_checkpoint.npz`.

Their full SHA-256 hashes and manifest-by-manifest verification appear in `STATE_INVENTORY.json`. A matching checksum establishes file integrity, not mathematical correctness, independent peer review, or external novelty.

Source changes and open questions should be carried forward as a new dated checkpoint; originals should remain unchanged.
