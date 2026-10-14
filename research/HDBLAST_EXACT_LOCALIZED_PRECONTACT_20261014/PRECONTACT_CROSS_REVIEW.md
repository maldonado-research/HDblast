# Independent review of the exact short precontact corollary

14 October 2026. Internal AI-assisted review, not external peer review.

Reviewed the precontact section of `referee/FINITE_REGULARITY_EXACT_INITIAL_DATA.md`, including its fixed-`k≥12` initial-data theorem. The substantive short-time argument passes review with the conservative four-order parameter-loss statement below. The compact four-torus auxiliary-extension correction and four-order parameter margin are adopted in the final revised proof and were verified on readback.

Final reviewed theorem SHA-256: `95ccaa550309b72226f8d03d3a950e4e10760d0d98810f0e818d19057a8cfd0b`. Acceptance is for that exact finite-regularity statement, with its explicit hypotheses and limits.

## Regularity and parameter variations

The exact initial metric is `H^{k+2}` and its canonical momentum is `H^{k+1}`. Converting momentum to extrinsic curvature is algebraic in the uniformly positive metric and its determinant, so the same velocity order is retained. Initial gauge lapse/shift time derivatives use first spatial derivatives and the canonical fields; they therefore fit `H^{s−1}` when fields lie in `H^s`, `s=k+2≥14`.

For a reduced system `g^{αβ}(u)∂_α∂_βu=F(u,∂u,x)`, the first parameter variation `w` has a source containing `w∂²u`. With `u∈H^s`, this fits `H^{s−2}`, yielding the energy class `w∈H^{s−1}`, `∂_tw∈H^{s−2}`. The second variation source includes `w∂²w`, which fits `H^{s−3}`; thus the second variation has energy regularity `H^{s−2}`, with velocity `H^{s−3}`.

Difference-quotient convergence must be checked in lower orders. Comparing the first quotient to `w` introduces `(g_{ε+h}−g_ε)∂²w`, whose source has only `H^{s−3}` regularity. The comparison energy estimate therefore gives convergence at least in `H^{s−2}`. The second quotient comparison introduces second spatial derivatives of that first comparison and second derivatives of the second variation; its source is controlled in `H^{s−4}`, giving convergence in a further reduced order. To avoid claiming an unquoted optimal tame/Bona–Smith flow result, the accepted corollary conservatively states actual C² parameter dependence in fields `H^{s−4}` and velocities `H^{s−5}`. The available high-regularity products and homogeneous linear energy estimates cover these comparisons and continuity on the same short interval.

These are derivative-loss estimates along the smooth initial parameter curve. No same-topology smoothness of the full quasilinear flow map is asserted. With `s≥14`, four additional time/spatial differentiations after the four-order parameter margin leave `H^{s−8}`, hence at least `H^6`; this is above the four-dimensional continuous-embedding threshold `n/2=2`. The stated C² parameter dependence in a piecewise spacetime C⁴ norm remains adequately covered.

## Gauge and constraint propagation

The background-adapted fixed reference connection is essential: using an unrelated harmonic coordinate condition could destroy equality with the original-coordinate background collars. With the selected reference, the untouched initial gauge values and derivatives agree with the background. Exact Hamiltonian/momentum constraints and the local gauge construction give both zero gauge deviation and zero normal derivative where required.

The scalar equation makes the stress-energy divergence vanish; the contracted Bianchi identity then yields homogeneous gauge propagation. Recovering the full equations is restricted to the domain of dependence of constraint-satisfying initial data.

## Uniform actual cones and patching

The radii `R_±` enclose the **whole exact correction annulus**, rather than only the original first-order packet. Continuous dependence and regular Lorentzian geometry give a common coordinate causal-speed bound `c>1` for actual perturbed metrics; the conformal background's speed one is only a reference value.

The conditions

```
2cτ < B−R_+,
cτ < R_-/2
```

produce both an exterior overlap `R_++cτ < R < B−cτ` and a protected inner region. On the exterior overlap, finite propagation and uniqueness identify the evolved solution with the auxiliary and actual backgrounds. Patching to the actual background there introduces no new interface. The arbitrary auxiliary extension cannot affect the physical solution on this protected region.

The revised proof specifies the auxiliary construction: embed the controlled ball in a large flat coordinate chart of `T^4`, extend the reduced initial fields and reference connection smoothly, and use the reduced Cauchy problem on that complete compact manifold. This avoids treating a nondecaying background/potential as an `H^s(R^4)` datum. Constraints outside the ball are neither needed nor claimed.

## Shell and remaining scope

The corollary explicitly assumes the actual background is smooth on the needed initial spacetime neighborhood and near its timelike shell; this is not inferred from the localized constraint solver. Both evolved copies equal that same background on a common open shell collar. Their background shell embeddings consequently satisfy the exact material junction and normal-force equations on this sufficiently short interval.

This proves a genuine short precontact continuation in the stated finite regularity class, conditional on the initial theorem and the regular-background neighborhood assumption. It gives no through-contact evolution theorem, no long contact-time bound computed solely from background cones, no spatial C-infinity family, and no external novelty conclusion.

Primary-source locations supporting the Cauchy and gauge machinery are recorded in `CAUCHY_SOURCE_LEDGER.md`.

