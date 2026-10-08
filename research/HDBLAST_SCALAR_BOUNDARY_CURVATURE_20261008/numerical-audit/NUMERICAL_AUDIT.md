# Reproduction and exact conditional remainder audit

8 October 2026. Internal AI-assisted research audit, not external peer review.

The uploaded archive is intact: all 27 manifest payloads match their recorded
byte counts and SHA256 values. The separately uploaded assessment PDF is byte
identical to the embedded copy (58,491 bytes; SHA256
`c8f2c1f74b2a8886ec8db93977322d06c0a352de69852a7d4b13ffbbbf190581`).
The reference files were only read. Executed copies and new output are confined
to this audit directory.

## Reproduction

All nine portable exact algebra checks reproduce, with a byte-identical result
JSON under SymPy 1.14.0. All 48 final numerical integrations reproduce their
stated consistency checks. These represent **24 family/detuning configurations
at two numerical settings**, comprising 42 fitted final integrations and six
integrations of the analytically fixed zero-coupling branch. Intermediate
shooting trials are additional integrations, not included in the count.

The replay uses Python 3.12.14, NumPy 2.2.6 and SciPy 1.15.3. The attachment
records Python 3.12.14, NumPy 2.3.5 and SciPy 1.18.1. No dependency was installed
or changed. The source's historical `/private/tmp/aps-research-deps` directory
does not exist in this environment; reviewed copies of the three program files
retain their original bytes. Main programs were run with Python `-I -B`; the
matched program received only its reviewed audit directory as an explicit local
import path.

| Recalculated diagnostic | Maximum |
| --- | ---: |
| Original versus replayed H² | 9.44612 × 10⁻¹⁷ |
| Original versus replayed boundary scalar | 6.39680 × 10⁻¹⁸ |
| Standard versus refined H² | 6.38189 × 10⁻¹⁷ |
| Recorded junction residual | 1.11554 × 10⁻¹⁵ |
| Recorded sampled constraint residual | 3.74340 × 10⁻¹⁵ |
| H² versus boundary identity | 7.84150 × 10⁻¹⁷ |
| Observed absolute (H² − P₂) / δ³ | 0.0084864690 |
| Observed absolute (η − aδ) / δ² | 0.0639765950 |

Table values round the recorded maxima upward for presentation. The full values
and per-family results are in `REPLAY_AUDIT.json`. The independent postprocessor
recomputes the polynomial boundary identity and susceptibility using exact
fractions for the represented binary64 coefficients and saved scalar values. It
also recomputes the pass thresholds from the saved diagnostics. It cannot
recover an unsaved continuous trajectory from those diagnostics.

Across the seven nonzero-coupling families, coefficient-error halving ratios
range from 1.9995255 to 2.0003160. At δ = 0.0005 the relative coefficient
discrepancy ranges from 0.0264301% to 0.190904%. This reproduces the reported
third-order trend. For zero scalar coupling, relative coefficient error and a
halving ratio are deliberately omitted: the target coefficient is zero and
roundoff dominates. In the stated domain w > 4k, the leading scalar junction
inverse 1/[2(w − 2k)] is finite. No claim about the full shooting Jacobian's
condition number follows.

This is a cross-environment rerun of the **same solver**, with an independently
written result postprocessor. It is not an independent numerical integrator.
The root flags, sampled residuals and agreement between settings provide useful
diagnostics but do not establish existence, uniqueness, continuous error bounds
or stability.

## A computable remainder implication

Write the exact polynomial boundary identity as

    H² = δ A(η) + δ² B(η),
    A = W f/9 − W′ f′/12,  B = f²/36 − (f′)²/48.

Suppose an exact branch independently satisfies

    |η − aδ| ≤ Kη δ²,   0 < δ ≤ D,
    a = −f₁/[4(w − 2k)].

Set E = |a| + Kη D. If A(η) = Σ Aₙηⁿ and B(η) = Σ Bₙηⁿ,
the triangle inequality proves

    |H² − P₂(δ)| ≤ C_H δ³,
    C_H = |A₁|Kη
          + Σ(n≥2) |Aₙ| Eⁿ Dⁿ⁻²
          + Σ(n≥1) |Bₙ| Eⁿ Dⁿ⁻¹,

where P₂ is the second-order expression in the attachment. Indeed, subtracting
P₂ leaves δA₁(η − aδ), the higher-order A terms, and the nonconstant B terms;
|η| ≤ Eδ bounds each of them. This implication covers the entire stated
interval, rather than a list of sampled detunings.

`conditional_remainder.py` evaluates this bound using exact rational arithmetic
for the eight polynomial models. The accompanying example files choose D =
0.002 and **assume Kη = 1**. That value is not inferred or certified from the
observed maximum 0.06398: a finite collection of floating-point solutions cannot
prove a uniform scalar estimate. These files certify only the conditional
algebraic implication. A separately proved scalar bound can be substituted
without changing the implementation or its proof.

The implementation controls cover exact envelope points, detect omission of
the scalar feedback coefficient in all seven nonzero-coupling models, and
reject invalid domain inputs. The normal and optimized-Python results are byte
identical. Sampling is used to check the implementation; the displayed
triangle-inequality argument proves the continuum implication. No claim of
external mathematical novelty is attached to this elementary bound.

## Files and reproduction

- `INPUT_INTEGRITY.json`: all attachment payload hashes and sizes.
- `STANDALONE_PDF_IDENTITY.json`: independent attachment comparison.
- `replay/`: the three unchanged sources and their freshly computed results.
- `FAMILY_REPLAY.log`, `MATCHED_REPLAY.log`: source program completion logs.
- `audit_results.py`, `REPLAY_AUDIT.json`: independent arithmetic diagnostics.
- `conditional_remainder.py`, `CONDITIONAL_REMAINDER.json`,
  `CONDITIONAL_REMAINDER_OPTIMIZED.json`: exact conditional implication.

From this directory, with the stated scientific dependencies available:

```sh
python -I -B audit_results.py --reference ../reference/HDBLAST-APS-research-review-20261008/general-scalar-response --replay replay --out REPLAY_AUDIT.json
python -I -B conditional_remainder.py --replay-dir replay --out CONDITIONAL_REMAINDER.json
python -I -O -B conditional_remainder.py --replay-dir replay --out CONDITIONAL_REMAINDER_OPTIMIZED.json
```

The numerical replay has already been completed; the commands above recompute
its diagnostics and polynomial certificate without further ODE integrations.
Files here do not establish APS readiness, originality, a Big Bang mechanism,
observational evidence, or a certified solution of the full boundary problem.
