# HDBLAST Chat 11 — the perturbation equations are now verified symbolically

18 September 2026 · Ricardo Maldonado's higher-dimensional blast research · prepared with Claude (Anthropic)

**Both recent certificates — the instability of the registered shell (Chat 9) and the stability of the shell S₈⁄₅
(Chat 10) — rested on one premise that had only been checked numerically: the linearised perturbation equations.
That premise is now verified exactly, with a computer-algebra system.** Every one of the fifteen linearised Einstein
equations and the scalar equation was computed symbolically and shown to reduce to the equations the certificates use.

This strengthens the mathematics under the two results. It adds no new physics and does not bear on whether a
higher-dimensional blast caused the Big Bang.

## What was installed

With your permission, SymPy 1.14.0 (and its dependency mpmath 1.3.0) from the Python Package Index, into an isolated
environment at `D-Blast 3/.hdblast_venv` (40 MB). Your system Python was not changed. To remove it, delete that folder.

## What was verified

All checks are exact symbolic identities, for an **arbitrary** bulk potential `U(φ)` (so they hold for the registered
`W = 1 − φ + φ³/3` in particular), on the background `ds² = dy² + ρ(y)²γ` with `γ` the unit de Sitter metric and a
general scalar harmonic `□Y = μ²Y`.

**1. The bulk equations** ([verify_linearised_equations.py](HDBLAST_CHAT11_SYMBOLIC_VERIFICATION_20260918/verify_linearised_equations.py), 22 checks)

* the background equations;
* the off-diagonal equation is exactly `(∂-terms of Y) × (ξ + 2ψ)`, so `ξ = −2ψ`;
* the mixed equations are exactly `(∂Y) × [−3(ψ' − Hξ) − φ'χ]`, the Codazzi relation;
* **sufficiency:** with those two constraints and the master equation

      ψ'' + 2(ρ'/ρ − φ''/φ')ψ' + [ −(4/3)φ'² − 4(ρ'/ρ)(φ''/φ') + (2 + μ²)/ρ² ]ψ = 0,

  all fifteen components `R_AB − ∂_Aφ∂_Bφ − (2/3)U g_AB` and the scalar field equation vanish identically at first order;
* **necessity:** with only the constraints imposed, the `yy` equation equals `2Y ×` (master-equation residual), so the
  master equation is forced, not merely allowed;
* the three forms of the shell condition used in the project are the same expression:
  `χ' + 2φ'ψ + (σ_t''/2)χ = Bχ + 3(μ²+4)ψ/(ρ²φ') = (3ψ/φ')[(μ²+4)/ρ² + B·R]`.

**2. The boundary-condition derivation** ([verify_junction_logic.py](HDBLAST_CHAT11_SYMBOLIC_VERIFICATION_20260918/verify_junction_logic.py), 8 checks)

* the gauge-transformation rules (Lie derivative of the background) for all components;
* the first-order extrinsic curvature in Gaussian-normal gauge, `δK^μ_ν = ψ'Y δ^μ_ν + E'∇^μ∇_νY`;
* the trace-free Hessian of a homogeneous harmonic vanishes only when `T'' = T'`, which forces `μ² = −4`: those are the
  only harmonics for which the shell may bend in longitudinal gauge;
* the trace Israel condition is identically the Codazzi relation (given `φ' = −σ_t'/2`), and the scalar junction becomes
  `χ' + 2φ'ψ + (σ_t''/2)χ = 0`.

**3. The equations the interval certificates actually integrate**
([verify_certificate_odes.py](HDBLAST_CHAT11_SYMBOLIC_VERIFICATION_20260918/verify_certificate_odes.py), 10 checks)

The background system `(φ, s, H, w)`, the identity `G' = U₂ − G² + 4HG`, the first-order pair `(ψ, X)`, the Riccati
equation `R' = R² + (2G − 6H)R + (μ²+4)w − (2/3)s²`, its cone form and indicial root `r₊ = −5/2 − ν`, the formula
`B = G − 4H + σ_t''/2`, and the mismatch `b = (μ²+4)w + BR` are all exact consequences of the verified master system.

**4. The checks have teeth** ([negative_controls.py](HDBLAST_CHAT11_SYMBOLIC_VERIFICATION_20260918/negative_controls.py))

Seven deliberately wrong versions — `(2+μ²) → (3+μ²)`, `4/3 → 5/3`, a wrong Codazzi coefficient, `ξ = −ψ`, `4H → 3H`,
`μ²+4 → μ²+3` in the shell condition, and a wrong de Sitter metric — are each rejected (between 2 and 10 failed checks).

## What this changes

| Premise under the Chat 9 and Chat 10 certificates | Before | Now |
|---|---|---|
| Linearised bulk equations, constraints, master equation | hand-derived in two gauges; numerical check to 10⁻⁹ | **exact symbolic identity, necessary and sufficient** |
| Shell boundary condition algebra, no-bending argument, `μ² = −4` exception | hand-derived | **exact symbolic identities** |
| Certificate ODEs (Riccati form, mismatch `b`) | re-derived by hand by a referee | **exact symbolic identities** |
| Archived M462 root enclosure (Chat 9 only) | inherited | inherited (unchanged) |
| Hand lemmas of the interval proofs (parity, cone barrier, monotone comparison, counting theorem) | not machine-checked | not machine-checked (unchanged) |
| Self-adjointness / completeness (from "no eigenvalue" to "mode-stable") | cited | cited (unchanged) |

The two safe statements therefore tighten to:

* *Chat 9:* the registered shell has a normalisable scalar perturbation with `−7.71788 < m²/H² < −7.71786`, conditional
  only on the archived M462 enclosure and standard analysis lemmas.
* *Chat 10:* S₈⁄₅ exists and has no scalar bound state below `9H²/4`, conditional only on the hand lemmas listed above.

Still not addressed: the two starting points of the junction analysis (that Gaussian-normal coordinates with the shell
at fixed position exist, and that the perturbed Israel conditions are the first-order expansion of
`K^μ_ν = σ_t(φ)δ^μ_ν/6`, `φ' = −σ_t'(φ)/2` — both standard), the vector sector, the `μ² = −4` harmonics, and non-linear
stability. A proof-assistant formalisation has not been attempted.

## Reproduction

In `HDBLAST_CHAT11_SYMBOLIC_VERIFICATION_20260918/`:

    "/Users/ricardomaldonado/Documents/D-Blast 3/.hdblast_venv/bin/python" REPLAY.py

runs the three verification scripts and the negative controls (about 30 seconds; result in `REPLAY_RESULT.json`).
Any Python ≥ 3.9 with SymPy ≥ 1.12 works. Nothing was published or sent; no original research file was modified.

## Next

The decisive open calculation is unchanged: the five-dimensional non-linear evolution of the Chat 9 roll-off. With a
computer-algebra system now available, two further items become practical: the vector sector and `μ² = −4` harmonics for
S₈⁄₅ (to upgrade "scalar mode stability" to linear stability), and a symbolic derivation of the four-derivative terms
of the four-dimensional effective action.
