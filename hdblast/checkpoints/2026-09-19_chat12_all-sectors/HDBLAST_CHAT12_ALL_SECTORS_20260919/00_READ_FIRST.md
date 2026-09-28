# HDBLAST Chat 12 — the stable shell has no unstable mode in any sector

19 September 2026 · Ricardo Maldonado's higher-dimensional blast research · prepared with Claude (Anthropic)

**Chat 10 proved that the stable shell S₈⁄₅ has no unstable *scalar* mode. This continuation covers everything that
was left out: the graviton (tensor) sector, the vector sector, and the five special scalar modes that the earlier
analysis had to exclude. None of them contains an unstable mode.** The checks are exact computer-algebra identities;
a skeptical referee re-ran them, broke them on purpose 22 times to confirm they detect errors, and replaced two of my
arguments with stronger ones.

This is a statement about perturbation *modes* of one shell in the registered five-dimensional model. It is not a proof
of full linear or non-linear stability, and it says nothing about whether a higher-dimensional blast caused the Big Bang.

## The result

For the shell S₈⁄₅ (`σ_t = 2W + t(1 + cφ + dφ²/2)`, `t = 1/1000`, `c = 5975949350280/10¹³`, `d = 8/5`), whose existence
is interval-certified (Chat 10), and with masses measured in units of the de Sitter slicing curvature:

| Sector | Finding | How it is established |
|---|---|---|
| Scalar, all harmonics with a non-vanishing trace-free Hessian | no normalisable mode with `m² < 9/4` | interval certificate (Chat 10), equations verified symbolically (Chat 11) |
| The five special scalar harmonics (`∇_μ∇_νY = −γ_μνY`) | no physical mode | every bulk perturbation is exactly a coordinate change; displacing the shell violates the scalar junction condition by `ζ Y φ' B`, and `B = +5.3187×10⁻⁴ ≠ 0` is certified |
| Vector | no mode | bulk equations plus the Israel condition force the gauge-invariant combination `σ_V = B_v − ρ²F_v'` to vanish; Killing harmonics are pure gauge |
| Tensor (all polarisations, including the helicity-1 and helicity-0 parts of massive gravitons) | one massless graviton, continuum from `9/4`, nothing else | factorisation `−∂_z² + V₁ = Q⁺Q` (no tachyon); exact identity `V₂ − 9/4 = ρ²(9φ'²/16 − U/8)` with `U < 0` along the certified bulk; second proof by a Riccati barrier; numerical spot check |

So: **S₈⁄₅ is mode-stable in every sector of the four-dimensional-covariant decomposition.**

The same tensor, vector and special-harmonic arguments apply to the unstable registered shell of Chat 9
(`U < 0` up to `φ = 3/20`; `B = −2.68×10⁻⁴ ≠ 0`). Its only unstable mode is the scalar one found in Chat 9
(`m² = −7.7179`), so that instability is now fully localised.

## What was checked, and by whom

**Orchestrator, exact symbolic scripts** (SymPy, 22 checks, about one minute:
[tensor](verify_tensor_sector.py),
[vector](verify_vector_sector.py),
[special harmonics](verify_special_harmonics.py)):

* the tensor radial equation `h'' + 4Hh' + m²h/ρ² = 0` and shell condition `h' = 0`; its Schrödinger form; both
  factorisations; the partner-potential identity; exact Sturm root counts showing `U < 0` on `[−1, 3/20]` and
  `U ≤ −1/6` on `[−1, 0]` for the registered potential;
* the vector gauge rule and field equations `(ρ²B_v)' = 0`, `B_v(□+3)V = 0`;
* two explicit special harmonics; the two-parameter gauge family solving the master equation at `μ² = −4` (the
  cone-regular member is the rigid translation of the apex); the displaced-shell computation — the trace Israel
  condition holds identically and the scalar-junction mismatch is exactly `εζYφ'B`.

**Numerical spot check** ([tensor_float_spotcheck.py](tensor_float_spotcheck.py)):
on both shells the shell mismatch vanishes at `m² = 0`, is negative for every `0 < m² < 9/4`, positive for every
`m² < 0`, with no sign change and no node.

**Mathematical referee** ([report](referee_math/REFEREE_MATH.md)): not refuted.
It re-ran everything, wrote 22 negative controls (all rejected), and closed three gaps in my arguments with exact proofs:

1. my tensor check used a single polarisation — the referee verified the covariant equation for a general traceless
   tensor, so it holds for all polarisations;
2. my vector exclusion relied on regularity at the cone, which is only heuristic because the cone is a horizon — the
   referee showed the Israel condition itself forces `σ_V(y_b) = 0`, with no regularity assumption;
3. my special-sector argument fixed a gauge first — the referee proved without gauge fixing that the non-gauge remainder
   vanishes wherever `φ' ≠ 0` and `ρ' ≠ 0`, both certified.

It also supplied the boundary-term details of the tensor argument and a second, elementary proof.

**Claims referee** ([report](referee_claims/REFEREE_CLAIMS.md)): acceptable as a
mode-stability statement only. Two wording errors of mine were corrected: `inf V₂` equals `9/4` (reached at the cone),
so the argument is pointwise positivity, not a gap; and the special sector needs `B ≠ 0`, not `B > 0`.

## Corrections to earlier wording

* Chats 9–11 described the excluded modes as "`μ² = −4`". The precise criterion is whether the trace-free Hessian of the
  harmonic vanishes. Infinitely many harmonics have `□Y = −4Y` without being special; they belong to the ordinary scalar
  sector, and the Chat 10 certificate does cover them (it holds for every real `μ² < 9/4`).
* The scalar/tensor split overlaps at `m² = 0` and `m² = 2`; both analyses exclude the overlapping configurations
  consistently.

## Not established

* One configuration is not covered by any sector: the constant harmonic (`Y = const`, `μ² = 0`). It is the linearised
  form of *uniqueness* of the shell solution, which is still open. It is not an instability.
* Completeness, self-adjointness and decay of the continuum: "no unstable mode" is not yet "linearly stable".
* The hand lemmas inside the interval certificates; that perturbations odd under the mirror symmetry are excluded
  (assumed); non-linear stability; stability against tunnelling (Chat 10 estimates `e^{−7270}`).
* Novelty: none claimed. The de Sitter graviton gap, the absence of vector modes and the translation modes are standard
  (Garriga–Sasaki; Langlois–Maartens–Wands; Gen–Sasaki; Garriga–Vilenkin; DeWolfe–Freedman–Gubser–Karch). The identity
  for `V₂` was not found stated in a bounded search but is a two-line consequence of the background equations. Thick
  de Sitter branes *with* a massive bound graviton below the gap exist in the literature, so its absence here is a
  property of this model.

## Reproduction

In ``:

    "/Users/ricardomaldonado/Documents/D-Blast 3/.hdblast_venv/bin/python" REPLAY.py

(about one minute; needs the SymPy environment installed in Chat 11). The numerical spot check uses the system
`python3` with numpy and needs folder 146 in place (about six minutes). Nothing was published or sent; no original
research file was modified.

## Next

The mathematics of the shells is now in good order: an unstable registered shell with a certified growing mode, a stable
shell with no unstable mode in any sector, a closed-form four-dimensional description linking them, and symbolic
verification underneath. The decisive open question for the hypothesis itself is unchanged and is physical, not
mathematical: **what does the five-dimensional, non-linear evolution of the roll-off produce?** That calculation
(a 1+1-dimensional evolution code with the shell as a boundary) is the recommended next step. It needs a multi-agent run.
