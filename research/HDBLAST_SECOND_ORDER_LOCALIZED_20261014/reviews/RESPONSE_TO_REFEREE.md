# Response to the referee report (14 October 2026)

Every point in `REFEREE_REPORT.md` was accepted. The changes are in
`manuscript/SECOND_ORDER_LOCALIZED.md`, `README.md`, `NEXT_CALCULATIONS.md`,
`basis/README_PROVENANCE.md` and `reproduce.py`.

| Referee point | Action |
|---|---|
| Lemma 1(a): wrong reason in the `λ` case | Fixed: `(1/Ω)″ < 0` for `Ω = c s^{−1/2}` and `= 0` for `1/Ω` linear in `s`; both contradict `(1/Ω)″ > 0` |
| Lemma 1(b): `l = 0` case where `Ω + 2sΩ′ = 0` | Adopted the referee's pointwise argument; noted that the `l = 2, 3, 4` ranks are sample-point checks |
| Lemma 1(c): "not needed below" is false | Rewritten: full KIDs are used for necessity in Theorem 2 |
| Theorem 2: vacuous support statement; `v` must solve the linearized constraints; theorem numbering | Apply Delay on `𝒜′ ⋐ 𝒜`; hypothesis added; numbering tied to arXiv:1003.0535v4 (the numbers match the text of that version) |
| Corollary 2.1: missing hypothesis | Added: first-order shell displacement and tilt vanish at `T = 0`; claim restricted to `T = 0` conditions |
| Corollary 2.3: last sentence needs hypotheses | Added `K = ωG`, `ω ≠ 0`, `B_rad(G, G) ≠ 0` |
| Abstract: "settle the initial-data question" too broad | Limited to corrections vanishing near the shell, and to fixed parity for the global tests |
| **Theorem 3 / Corollary 3.1: `m` is not invariant under first-order gauge changes of `v`** | **Main correction.** Title, abstract, Theorem 3, Corollary 3.1, README and next calculations now state that `m` is fixed for `v` in the reconstruction gauge and that a first-order gauge change shifts it. The physical reading "the packet must deform the geometry" is withdrawn as a gauge-invariant claim. The referee's formula for `m[L_ξ²γ₀]` is quoted |
| Theorem 3: cleaner adjoint proof | Adopted: `P_ψ*(Ωr) = 0`, `P_g*(Ωr) = −Ω(2sδ + xxᵀ)/(2s^{3/2})` |
| Corollary 3.2: "short radial wavelength" imprecise | Replaced by thin support via Wirtinger, or rapid oscillation |
| Section 4: bump interval must contain `supp f₀` | Fixed |
| "Checked exactly" before receipts existed | Section 6 and `RESULTS.md` now cite only completed receipts; the original `verify_cc.py` terminated without output and is reported as such |
| Missing `RESULTS.md`, `MANIFEST.json`; wrong path; digits | Added; path fixed; README states 40-digit evaluation and 70-digit differentiation |
| `reproduce.py` ignores receipt contents | Now asserts the pass/fail fields inside each receipt |
| Basis note §0 item 4 overclaim | Pointer added in `basis/README_PROVENANCE.md`; the copied note is left unchanged as a record |
