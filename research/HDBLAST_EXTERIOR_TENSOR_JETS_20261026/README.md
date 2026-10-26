# Exterior tensor coefficient two-jets

26 October 2026. Conditional mathematics for the retained Einstein–scalar shell model.

At the exterior point `(x,r)=(1/4,5/4)`, this checkpoint encloses the actual total registered tensor coefficient and all of its first and second coordinate derivatives. The two retained tensor amplitudes remain exactly equal, so write `u=T+=T−`. Both the coefficient and its scaled time derivative are strictly negative. The three second-derivative intervals contain zero: their signs are **inconclusive** at the declared maximum budget.

| Scaled coefficient derivative | Outward enclosure |
| --- | --- |
| `u` | [-0.000207982516, -0.000010433243] |
| `u_x` | [-0.000280539069, -0.000046459086] |
| `u_r` | [-0.000127797148, 0.000106282835] |
| `u_xx` | [-0.000848299839, 0.003927176089] |
| `u_xr` | [-0.001583379872, 0.003192085893] |
| `u_rr` | [-0.000848295221, 0.003927170544] |

Here `x=t/ell`, `dy=−ell rho dr`, and `ell=10^−6`. These are ordinary derivatives of a registered-gauge metric amplitude, not covariant tensor derivatives, curvature components, detector strain, or a gravitational-wave observable. All displayed endpoints round outward from the exact rationals in `compute/TT_TWOJET_NORMAL.json` and `RESULTS_SUMMARY.json`. The six intervals are separate enclosures; they do not assert that every Cartesian combination of endpoints is jointly attainable.

The new value interval is a subset of the previous frozen tensor-value interval. The earlier full23 global approximation and all frozen manuscripts/capsules remain unchanged. This is a pointwise tensor two-jet certificate, not a full-domain C2 approximation or a complete natural Einstein/scalar residual certificate.

## Why the derivative calculation closes

The full tensor source has 78 terms. Its only second physical metric derivatives occur in `Myy`, `Myt`, and `Mtt`, while `M=2d_gauge` contains no derivative in its reconstruction. One coordinate radial derivative of the entire source therefore needs gauge derivatives through order three and master derivatives through order two. Existing authenticated actual C3 error bounds suffice. This is a property of this exact tensor row; no general reduction of Einstein–scalar derivative requirements is claimed.

For `u_xx−u_rr=A u_x+B u_r+D u+F`, spatial differentiation gives `v=u_r` with the same highest-order operator and extra forcing `F_r+A_r u_x+B_r u_r+D_r u`. Both initial tensor data vanish and the backward cone avoids the shell, so `v` also has zero initial data. A triangular resolvent bounds the complete feedback of `(u,v)`. Characteristic integrals enclose first derivatives and radial/mixed second derivatives; the equation recovers `u_xx`. The existing C1 error number is never differentiated.

The source, its radial derivative, actual input errors, coefficient derivative cross-couplings and all remaining scalar feedback are retained. The proof is in `theory/DERIVATION.md`; independent source and numerical reviews are in `review-inputs/` and `review-numeric/`. Ordinary `t,y` coefficient derivatives, including the nonconstant radial conversion term, are also enclosed in the primary report. They carry the same gauge and interpretation limits.

## Bounded campaign and verification

The declared first budget used 20 clipped cone cells and eight segments on each characteristic. No complete second derivative excluded zero, so the single declared refinement used 72 clipped cone cells and sixteen segments per characteristic. The diagnostic sign gate remained inconclusive and the campaign stopped. The primary finite two-jet enclosure succeeded. No further mesh or parameter search was performed.

Primary calculations run in normal and optimized Python with byte-identical scientific reports, including relocated repository/code/input paths. Independent exact algebra, physical-chain, polygon/characteristic, resolvent and wrong-formula controls accompany the package. Internal AI-assisted review is distinct from external specialist peer review and proof-assistant formalization.

## Replay and scope

Python 3.12 and the standard library suffice for the runtime replay. Keep the manifest SHA256 outside this capsule. The dedicated workflow supplies the external source, archive and entry-point pins.

```
python -I -S -B reproduce.py --manifest-sha256 EXTERNAL_PIN \
  --repo-root /path/to/HDblast --output-dir /tmp/fresh-tensor-jet-replay
```

The wrapper authenticates the new source closure and nine complete frozen upstream capsules before project execution and after replay. It regenerates the new bounded calculations and independent checks in both modes; it does not rerun old production campaigns. The SymPy fixture generator is retained as provenance; its pinned exported fixture is independently checked using standard-library exact arithmetic during replay.

Read `literature/LITERATURE_SCOPE.md` for the bounded primary-source comparison. No external novelty/significance assessment or observational link is established here. Retained smooth background/master/gauge, differentiated-error, formal coefficient-existence and subsidiary premises remain explicit. Nonlinear finite-amplitude realization, a well-defined physical observable, a cosmological prediction and the proposed Big Bang origin remain open. No new Zenodo version or journal submission accompanies this checkpoint. Private editorial history is excluded.
