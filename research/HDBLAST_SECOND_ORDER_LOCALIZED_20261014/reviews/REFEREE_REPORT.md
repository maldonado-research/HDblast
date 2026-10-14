# Referee report: SECOND_ORDER_LOCALIZED.md (14 Oct 2026)

Adversarial internal review (AI referee, not external peer review). Scope: the
manuscript, README, NEXT_CALCULATIONS, receipts in `new/`, `independent/` and
`basis/`, and the shell-Noether checkpoint [1]. I edited no project file except
this one.

## Verdicts

| Item | Verdict |
|---|---|
| 1. Lemma 1(a) | CORRECT WITH FIX (λ-case reason is wrong, conclusion true) |
| 1. Lemma 1(b), SO(4) reduction | CORRECT |
| 1. Lemma 1(c), "not needed below" | WRONG (it is needed for necessity) |
| 2. Theorem 2 / Delay | CORRECT WITH FIX. Delay's theorem and lemma numbers UNVERIFIED |
| 3. Cor. 2.1 | CORRECT WITH FIX (missing hypothesis) |
| 3. Cor. 2.2 | CORRECT |
| 3. Cor. 2.3 equivariance | CORRECT. Final "necessary and sufficient" sentence needs hypotheses |
| 4. Theorem 3 identity, Ψ-independence, radial-gauge invariance | CORRECT |
| 4. Name "gauge-invariant monopole moment" | WRONG as phrased (counterexample below) |
| 4. Cor. 3.1 formula and sign | CORRECT. Physical reading WRONG as phrased |
| 4. Remark (Ωr, radial boost) | CORRECT |
| 4. Cor. 3.2 (c₂, cross term) | CORRECT WITH FIX (precise "short wavelength") |
| 4. Sec. 4 l = 0 construction | CORRECT WITH FIX (bump interval) |
| 5. Overclaims and housekeeping | several fixes, listed below |

## 1. Lemma 1

**(a)** The Liouville reduction and the two separated equations are right (I
re-derived both: `div Y = 4λ + 8b·x`, `Y(ln Ω²) = (4Ω′/Ω)(a·x + λs + s b·x)`). The
error is the sentence "Either case makes (1/Ω)″ = 0". For `Ω = c s^{−1/2}` one
gets `(1/Ω)″ = −s^{−3/2}/(4c) < 0`, which is not zero. The contradiction still
holds, because the background equation gives `(1/Ω)″ = (Φ′²/3)/Ω > 0`.
*Replace with:* "A nonzero λ forces `Ω = c s^{−1/2}` with `c > 0`, so
`(1/Ω)″ = −s^{−3/2}/(4c) < 0`. A nonzero `(a, b)` makes `1/Ω` linear in `s`, so
`(1/Ω)″ = 0`. Both contradict `(1/Ω)″ = (Φ′²/3)(1/Ω) > 0`." The same fix
applies where (b) relies on "would force Ω ∝ s^{−1/2}".

**(b)** Correct. For l ≥ 2: the Christoffel terms
`w_i∂_jN + w_j∂_iN − δ_ij w·∂N` with radial `w` are tangentially pure trace. The
flat Hessian differs from the sphere Hessian by a pure-trace term, and
`∫|(∇²Y)^TF|² = (2/3)λ(λ−3)∫Y²` (Bochner, Ric_{S³} = 2g). A cleaner version
needs no isotypic step: on each sphere the tangential trace-free part of
`P_g*N` is `½(Hess_{S_r}N)^TF`, so `N|_{S_r} ∈ span{1, x_a}`. The l = 1 and
l = 0 reductions match my independent code (below) exactly. Simplification for
l = 0: where `Ω + 2sΩ′ = 0`, the eliminated relation reads
`nΩ′(Ω + sΩ′) = 0` with `Ω′ = −Ω/2s ≠ 0`, so `n = 0` pointwise there as well.
The unique-continuation step is then unnecessary. The tensor-only rank
statements for l = 2, 3, 4 are checked at two sample points; say so.

**(c)** The statement "it is not needed below" is false. The "only if" half of
Theorem 2 allows arbitrary `w`, including `ψ ≠ 0` and `ϖ′ ≠ 0`. It therefore
needs `P_ψ*(Ωx_a) = 0` and `Y_J φ₀ = 0`, i.e. that the ten are full KIDs.
*Replace with:* "The ten are full KIDs, being the T = 0 data of the background
boosts and rotations. This inclusion is used for necessity in Theorem 2. By
(a)–(b) there are no others."

## 2. Theorem 2

The decoupling at `π₀ = ϖ₀ = 0` is correct. The densitization is consistent:
`P_M(√γ₀τ, 0) = −2√γ₀ div τ`, the adjoint of div is `−½L_Yγ₀`, and the condition
`∫Y·S_M d⁴x` (S_M is a covector density) is the correct pairing. The condition
`∫Ωx_a S_H dμ_γ₀` matches the non-densitized `H` and the `dμ_γ₀`-adjoint
`P_g*`. Both principal symbols are injective: in n = 4,
`tr(ξ⊗ξ − |ξ|²γ)N = −3|ξ|²N`. Lemma 1 covers every smooth kernel element on the
open annulus, so any weighted-space kernel is included.

I could not extract the text of arXiv:1003.0535 in this review. The arXiv
abstract confirms the scope (underdetermined elliptic operators,
divergence-free fields and TT tensors), but not the numbers "Thm 1.3/5.1,
Lemma 8.1". Cite the version (v4 = Comm. PDE 37 (2012)) and check the
numbering.

Required fixes:
- "supported in a compact subset of 𝒜̄" is vacuous. Apply Delay on `𝒜′` with
  `K ⊂ 𝒜′ ⋐ 𝒜` (Lemma 1 holds on every annulus). Then `w` is compactly
  supported in `𝒜`.
- State that `v` solves the linearized constraints.
- Necessity: cite (c) as corrected above.

## 3. Corollaries 2.1–2.3

**2.1.** Add the hypothesis: "the first-order shell displacement and tilt and
their time derivatives vanish at T = 0". Restrict the claim to the T = 0
(initial) junction and compatibility conditions.

**2.2.** Correct. Quadratic sources are even under `v ↦ −v`, and the boost
lapse is `Ωx_a` in both copies. The two-copy gap is honestly stated. Change
the abstract's "We settle the initial-data question" to "We settle it for
fixed reflection parity, and in general within corrections vanishing near the
shell".

**2.3.** The argument is sound. Under SO(4) = (SU(2)×SU(2))/Z₂, Sym₀(4) = (1,1)
and `S²(1,1) = (0,0)⊕(1,1)⊕(2,2)⊕(0,2)⊕(2,0)` contains no so(4) component. Also
`Λ²(1,1) = (1,0)⊕(0,1)⊕(1,2)⊕(2,1)`, so the bilinear map is antisymmetric with
two constants `c±`, and reflections force `c₊ = c₋`. The receipts agree: the
non-commuting witness has density `∝ s²Q_rot/(Φ′²Ω)` only in J₁₂. The last
sentence needs hypotheses. *Replace with:* "…for fixed-parity quadrupoles with
`K = ωG`, `ω ≠ 0` and `B_rad(G, G) ≠ 0` (e.g. thin support, as in [1]),
`[A, B] = 0` is necessary and sufficient." Without these, every pair passes
(ω = 0), or necessity fails.

## 4. Theorem 3, Remark and Corollaries 3.1–3.2

**Identity: correct.** Independently (own jet code, no project imports):
- `P_ψ*(Ωr) = 0`.
- `P_g*(Ωr) = −Ω(2sδ + xxᵀ)/(2s^{3/2})`.
- `⟨P_g*(Ωr), αδ + a xxᵀ⟩_γ₀ = −3(3α + sa)/(2Ω³√s)`.

This reproduces `−(3/2)∫Ω√s(3α + sa)` and `m = −(2/3)∫s^{3/2}Ω⁵f₀`. It also
gives a cleaner proof with no boundary bookkeeping: for compactly supported `w`,
`∫Ωr P_H[w] dμ = ∫⟨P_g*(Ωr), h⟩dμ` exactly. Suggest using it. Add that `Ωr`
spans the l = 0 kernel of `P_ψ*` (from `nΩ + 2nΩ′s − 2n′Ωs = 0`). I confirmed
the radial-gauge invariance symbolically, and checked the Sec. 4 `a`- and
`Ψ`-identities by hand. The Remark is correct: `(R∂_T + T∂_R)s = 0` and
`L_K g = 2RTΩ²dΩ²` vanishes at T = 0.

**Gauge invariance: overclaimed.** `m` is invariant only under linearized gauge
transformations of `w` at fixed `v`. Under a first-order gauge change
`x̃(ε) = exp(ελL_ξ)x(ε)` one has `w̃ = w + 2λL_ξv + λ²L_ξ²x₀`. I computed, for
`ξ = ζ(s)x`:
```
m[L_ξ²γ₀] = ∫[8s^{5/2}Ω³ζ′² + 8s^{3/2}Ω(sΩΩ″ − sΩ′² + ΩΩ′)ζ²] ds
```
In flat space this is `8∫s^{5/2}ζ′² > 0`. On the independent test background
`1/Ω = 1 + s²/5 + s³/7`, with the bump on (0.5, 1.5), it is `+6.34×10⁻³`. For a
velocity packet `h_v = 0`, so `L_ξv` has no metric part and
`m[w̃] = m[w] + λ²m[L_ξ²γ₀]`. Choosing `λ² = −m[w]/m[L_ξ²γ₀]` gives `m[w̃] = 0`
for the same physical packet in another gauge.

*Replace the last sentences of Theorem 3 with:* "m is invariant under
compactly supported linearized gauge transformations of w and does not involve
ψ. It depends on the gauge of v: it is fixed once v is given in the
reconstruction gauge (22) of [1]."

*Replace the last two sentences of Cor. 3.1 (and the matching abstract and
README lines) with:* "Hence, for velocity packets in the gauge (22) of [1],
m[w] < 0 for every compactly supported correction. The spherically symmetric
part of the second-order metric correction cannot be a compactly supported
gauge mode, and cannot be replaced by a scalar correction. This is a statement
about the gauge-fixed packet. A first-order gauge change of v shifts m (e.g.
by λ²m[L_ξ²γ₀])."

Rename the quantity to "the areal-radius-weighted monopole moment" and drop
"gauge-invariant" from the title, or qualify it as above.

**Cor. 3.1 formula: correct.** By hand:
`S_H = 6p²s²/Ω⁶ + 9(Dp+5p)²/(4Φ′²Ω⁶)`, with `⟨(x₁²−x₂²)²⟩_{S_r} = s²/6`. This
reproduces the stated density, which `compare_f0.json` also confirms. The tr B²
scaling holds.

**Cor. 3.2: correct with fix.** `c₂` comes only from `|∇δφ|²`, giving
`6s^{13/2}G″²/(Φ′²Ω)` ✓. The cross term vanishes ✓. "At most one ω²" ✓.
*Replace* "For profiles with short radial wavelength" *with* "For G supported
in a sufficiently thin interval (Wirtinger: `∫G″² ≥ (π/w)²∫G′²`,
`∫G′² ≥ (π/w)²∫G²`) or sufficiently rapidly oscillating".

**Sec. 4, l = 0.** The bump interval must contain the support. *Replace with:*
"…with `b ∈ C_c^∞((s₁, s₂))`, `b ≥ 0`, `b ≢ 0`, where `[s₁, s₂] ⊂ (r₁², r₂²)`
contains `supp f₀`." Otherwise `Ψ ≠ 0` below `s₁`.

## 5. Overclaims and housekeeping

1. "checked exactly" (abstract, §6) for `P_H[h_T] = S_H − f₀`: `new/fast_verify.json`
   does not exist yet; the job is still running, as is `basis` `verify_cc.py`,
   whose log is empty. What exists now: exact `P_M`, step 1 (no radial
   derivatives), the exact `linsolve` construction, and the 40-digit check at
   6 points (`indep_check`). Either wait for the receipt or state exactly this.
   `independent/indep_l0.json` is also pending; its log shows the weighted
   identity to 3×10⁻²⁵.
2. `RESULTS.md` (cited in §6 and README) and `MANIFEST.json` (README) are
   missing. The path `basis/results_snapshot/construct_cc_M.json` should be
   `basis/results/construct_cc_M.json`. README says "70 digits"; the code works
   at 40 digits, with 70 only for differentiation. Make them consistent.
3. `reproduce.py` checks only return codes and file existence, not the
   booleans inside the receipts. It should assert them.
4. `basis/SECOND_ORDER_LOCALIZED_CONSTRAINTS.md` §0 item 4 ("must deform the
   monopole metric (a localized backreaction)") carries the same gauge
   overclaim. It is copied unchanged, so add a pointer to the corrected
   statement.
5. The physics language "the packet must deform the spherically symmetric
   geometry" and "the scalar field cannot screen it" goes beyond the
   mathematics until item 4 is fixed. Everything else in "Not established" is
   appropriately cautious.

## Computations run (scratch `/private/tmp/claude-501/referee_so/`)

- `jet.py`, `check_l01.py`: independent radial-jet engine. Confirmed the l = 0
  `cd`, `cx`, `λ` and compatibility residual `−Φ′²Ω²n/(3(Ω+2sΩ′)²)`, and the
  l = 1 mixed component `∝ (Ωn′ − Ω′n)/Ω`. Confirmed `P_g*(Ωx₁) = 0`,
  `P_ψ*(Ωx₁) = 0`, `P_ψ*(Ωr) = 0` and `P_g*(Ωr)` as above.
- `gauge2.py`: linear radial gauge gives reduced density 0. Second-order gauge
  gives `m[L_ξ²γ₀]` as above.
- `gauge_num.py`: `m[L_ξ²γ₀] = +6.344×10⁻³` on the test background.
- By hand: Lemma 1(a) separation, the l ≥ 2 norm, the Sec. 4 l = 0 identities,
  the Cor. 3.1 density, `c₂`, and the SO(4) decompositions.
