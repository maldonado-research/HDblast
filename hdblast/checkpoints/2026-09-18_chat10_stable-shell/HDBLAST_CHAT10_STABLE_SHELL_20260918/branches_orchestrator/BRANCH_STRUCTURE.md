# Two branches of static shells: how the unstable and the stable shell are connected

18 September 2026 (Chat 10, orchestrator). Status: 4D effective-theory prediction (t → 0) confirmed by floating-point
solutions of the full 5D junction problem at t = 10⁻³. Not interval-certified. Family studied:
`σ_t = 2W + t(1 + c★φ + dφ²/2)`.

## Result

For every curvature `d` there are **two** static de Sitter shells near the wall, one stable and one unstable. They meet
and exchange stability at `d₀ = 1.1134966` (a transcritical bifurcation of the modulus potential `V_E = V/f²`):

* `d < d₀`: the shell at the wall centre (`φ_b ≈ 0`, the registered type) is the **unstable hilltop**; a **stable** shell
  sits on the throat side (`φ_b < 0`). As `d` decreases the stable shell slides toward the throat and leaves the wall
  (`φ_b → −1`) somewhere between `d = 0.9` and `d = 0.5`. At the registered `d = 0` only the hilltop remains.
* `d > d₀`: the shell at the wall centre is **stable**; the unstable shell sits at `φ_b > 0` and is the top of a small
  barrier separating the stable shell from the long slope toward `φ = +1`.

The effective theory predicted the second branch before it was looked for in 5D. All twelve predicted shells exist in
the full 5D problem (`branches_5d_check.py`, junction residuals ≤ 10⁻¹⁴):

| d | 4D prediction: φ_b, m²/H² | 5D solution: φ_b, m²/H² |
|---|---|---|
| 0.90 | −0.24761, +1.7124 | −0.24772, +1.7137 |
| 0.90 | 0, −1.4802 | +0.00012, −1.48102 |
| 1.00 | −0.11641, +0.8496 | −0.11664, +0.85279 |
| 1.00 | 0, −0.7869 | +0.00022, −0.78958 |
| 1.05 | −0.06156, +0.4593 | −0.06197, +0.46536 |
| 1.05 | 0, −0.4402 | +0.00040, −0.44573 |
| 1.30 | 0, +1.2930 | −0.00014, +1.2949 |
| 1.30 | +0.14175, −1.1492 | +0.14185, −1.15069 |
| 1.60 | 0, (+3.37: above the continuum) | −0.00005, no bound state |
| 1.60 | +0.29302, −2.5371 | +0.29299, −2.53654 |
| 2.00 | 0, (+6.15: above the continuum) | −0.00003, no bound state |
| 2.00 | +0.41727, −3.8363 | +0.41717, −3.83449 |

Differences are O(t), as expected. This is a test of the effective theory well away from the point where it was derived:
it locates new exact non-linear static solutions of the 5D equations up to `φ_b = −0.25` and `+0.42` and gives their
spectra to about 0.1–1%.

## Reading

1. **The stable shell is a false vacuum with a very low barrier.** For `d = 8/5` the barrier top lies only 0.72% above
   the stable shell's vacuum energy. Classically the shell is stable; quantum mechanically it can cross by the
   Hawking–Moss process with exponent `B = 24π²(1/V_min − 1/V_top) = 7.27/t` (4D estimate, `METASTABILITY_EFT.json`),
   i.e. `e^{−7270}` per Hubble four-volume at `t = 10⁻³`: irrelevant in practice. Coleman–De Luccia bubbles on the shell
   exist only when the barrier top has `m² < −4H²`, which needs `d ≳ 2.1`.
2. **Changing the tension slowly does not produce a blast.** If `d` drifts down through `d₀`, the shell follows the stable
   branch continuously toward the throat (a smooth, second-order-like change). The violent roll-off found in Chat 9
   (`m² = −7.7H²`) requires the shell to be *placed* on the hilltop, or the tension to change abruptly (a quench).
3. Beyond the barrier the slope toward `φ = +1` is steep (no slow roll) and ends in another de Sitter state, as found in
   Chat 9. Nothing here changes the negative conclusion about a hot Big Bang from this modulus alone.

## Not established

Interval certification of these shells (in progress for `d = 8/5`); the fate of the stable branch below `d ≈ 0.9`
(it leaves the tabulated range toward `φ_b → −1`); any 5D instanton; non-linear time evolution.
