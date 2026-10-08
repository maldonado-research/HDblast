# HDBLAST research checkpoints (September 2026)

**Current release — 8 October 2026 UTC:** 98,304 endpoint comparisons, exact normal/optimized identity, independently accepted fresh portable replay and conditional bulk/incident models. Fresh normal/optimized portable replay and numerical/bulk-model CI passed; final-document-head CI and live site deployment are recorded separately in the final delivery receipt. Latest verified Zenodo publication is **[23244754](https://doi.org/10.5281/zenodo.23244754) / concept 17088132 / 31 files / 453,646,943 bytes**. [Publication proof](https://github.com/maldonado-research/HDblast/blob/main/hdblast/publication/2026.10.08-retained-endpoint-flow/README.md) · [Replay and CI evidence](https://github.com/maldonado-research/HDblast/blob/main/hdblast/continuation/2026.10.08-retained-endpoint-flow/README.md). Historical entries below retain their dated scope.

Each folder is a research package extracted unchanged from its delivered archive (exceptions
below). Each contains a report (`README.md`), the code, the numerical outputs, and internal
referee or verification material produced with AI assistants. None has had external peer review. The original delivered archives are
preserved in the author's private research archive.

| Date | Checkpoint | Main result | Labels |
|---|---|---|---|
| 16–17 Sept | [Chat 9: shell dynamics](2026-09-16_chat09_shell-dynamics/) | The registered shell has a tachyonic scalar mode (m² ≈ −7.7179 H²; existence certified, uniqueness numerical). A closed-form 4D effective theory. Shell-modulus inflation is not viable (n_s ≤ 0.93). | Certified (conditional) / Numerical / Negative |
| 18 Sept | [Chat 10: stable shell](2026-09-18_chat10_stable-shell/) | With a modified (curved) tension, d = 8/5, a shell exists whose scalar sector has no unstable mode; the 4D formula's two-branch structure was confirmed at 12 shells in 5D | Certified (conditional) / Numerical |
| 18 Sept | [Chat 11: symbolic verification](2026-09-18_chat11_symbolic-verification/) | All 15 linearized Einstein components, the scalar equation and the junction algebra; 7 negative controls rejected | Exact |
| 19 Sept | [Chat 12: all sectors](2026-09-19_chat12_all-sectors/) | No unstable mode in the tensor, vector and special-harmonic sectors (mode stability); the only instability found is the Chat 9 scalar mode | Exact + referee |
| 21 Sept | [Chat 13: nonlinear roll-off](2026-09-21_chat13_nonlinear-rolloff/) | Full 5D evolution: two fates, heading toward an empty de Sitter brane or reversal and collapse; no radiation branch in the homogeneous, classical roll-off | Numerical / Negative |
| 22 Sept | [Chat 14: registered detuning](2026-09-22_chat14_registered-detuning-rolloff/) | Same two fates at the registered δ = 0.001; growth rate recovered to five digits | Numerical / Negative |
| 22 Sept | [Registered-shell controls](2026-09-22_registered-shell-controls/) | Replacement solver calibrated; exact constraint-transport identity; initial-data obstruction and repair candidates | Numerical / Exact |
| 22 Sept → Zenodo 23 Sept | [Scalar-profile branch and matter](2026-09-22_scalar-profile-branch-and-matter_zenodo-22922928/) ([record 22922928](https://zenodo.org/records/22922928), apparently this package) | The constant φ = 1 endpoint fails the scalar junction; corrected static branch; balanced-disturbance evolutions; proposed matter extension | Exact / Numerical / Negative / Conditional |
| 27 Sept | [Stability, exact series, particle production, new mechanisms](../../research/HDBLAST_CHECKPOINT_20260927/) ([public summary](../../research/HDBLAST_CHECKPOINT_20260927/PUBLIC_SUMMARY.md)) | The corrected end state is linearly stable in the tested sectors (two methods); exact series through δ⁸; no radiation era from particle production or six other routes in the registered model (leftover vacuum energy is the obstacle); constraint stall in the evolutions solved; a tension-term model change is inconclusive because of dark radiation | Numerical / Exact / Negative / Inconclusive |

Notes:
- One 33 MB data file (`eft_landscape.npz` in Chat 9) is not included here because of its size. It
  remains in the original archive in the author's private research archive.
- Absolute paths such as `/Users/...` inside logs, provenance files and some replay commands refer to
  the author's original workstation and are kept for provenance (a private workspace prefix in two
  provenance files is replaced by `<author-workspace>`). To replay, run the scripts
  from their own folder with `python3` (numpy, scipy, sympy, mpmath installed).
- Folder names use the checkpoint date; the Chat numbers are the author's working-session labels.
