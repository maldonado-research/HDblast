# Prospective high-precision continuum followup

Proposed checkpoint:
`HDBLAST_CHECKPOINT_20261002_METRIC_RESPONSE_FOLLOWUP`.
This document supplies a scientific rationale; it is not itself a completed
registration or authority to perform a physical calculation before the new
public freeze. No new source, response, mode or root is evaluated here.

The original metric calibration stopped on an explicit continuum quadrature
roundoff warning before completing its first nonzero row. A bounded followup
can change the continuum arithmetic and quadrature implementation while
retaining the same physical question and every original numerical gate.
The original failed experiment and all its outputs remain separate and
immutable.

## Preserve the physical comparison

Retain H=1, x=r=2, xi=0, fixed physical scalar and reference mass, fixed
incoming BD state, the same conformal metric pulse `a=a0(1+epsilon h)`,
epsilon=1e-4, h=B or uB on `[-5,-3]`, initial eta=-6, final eta=-1.5,
the six observation times and K=64,128,256. Keep the exact canonical source

```
g=4L²h-2Lh'-h'', L=-1/eta,
```

the direct metric observable/subtraction contacts, reference baselines,
normalizations, mass-law current convention, independent mode algorithm,
finite-K primary algorithm, stopping rules and existing acceptance gates.
Changes should be restricted to the continuum logarithmic-memory evaluation
and its declared precision/refinement evidence. Record any unavoidable
administrative paths, manifest pins or provenance changes separately.

The proposed continuum calculation uses independently executed **50- and
70-decimal-digit** evaluations of the same mathematical memory formula.
Its physical source jets must also be evaluated at the declared precision;
converting already rounded binary64 derivative values to mpmath would not
restore the missing digits. Exact integer derivative polynomials and
analytic background jets permit this without changing the source profile.

## Retain the matched local constant

For derivative order n=1,2,3, the continuum history expression is

```
I_n(eta)=integral_-5^min(eta,-3) g^(n)(t) ln(eta-t) dt
       +[ln(sqrt(2)/(-eta))+gamma_E+1] g^(n-1)(eta).
```

Before the pulse it is zero. When eta lies inside the pulse, the logarithmic
endpoint is integrable. A high-precision quadrature may use a frozen
coordinate transformation or exact endpoint subtraction, for example

```
T=eta+5,
integral_-5^eta g^(n)(t)ln(eta-t)dt
 =g^(n)(eta)T(ln T-1)
  +integral_-5^eta [g^(n)(t)-g^(n)(eta)]ln(eta-t)dt.
```

The second integrand tends to zero as `(eta-t)ln(eta-t)`. This equality is
analytic and does not alter the matching. The `gamma_E+1` local constant
and scale-dependent contact must remain intact. After the pulse the upper
endpoint is separated from observation, so the same expression has no
logarithmic endpoint singularity.

The implementation must freeze its actual high-precision quadrature rule,
panel nodes or endpoint transformation, endpoint branch definitions,
source-polynomial evaluation, precision contexts, maximum work and timeout,
reported refinement estimate, final rounding and failure behavior before
execution. Pure symbolic equality and synthetic tests can be completed
before that freeze; physical source sampling and response quadrature cannot.

## Precision comparison is evidence, not a total error certificate

Keep both 50- and 70-digit histories and the normalized response differences.
Higher precision and agreement can address the observed numerical barrier,
but agreement between two evaluations of the same formula is not an
independent physical route or an interval-certified total error bound. Both
can share a source/contact mistake or quadrature bias. The separately coded
forced-mode calculation, normal/optimized validators, raw mode reconstruction,
Ward/trace checks, negative controls and refinement evidence remain required.

Apply the original primary quadrature/refinement and cross-route gates
without relaxing them. Freeze exactly how the new continuum precision
comparison contributes to the reported quadrature estimate; do not present
a precision difference as a rigorous error bound. Any new warning, budget
overrun, nonfinite value or failed registered gate must stop the followup
with all attempted data retained. Do not add unregistered precisions,
switching profiles, states, cutoffs or tolerances after observing results.

The existing rigorous UV envelopes remain applicable because the canonical
source and finite matching are unchanged. Their global bounds are loose,
especially for q'' and pressure; a completed finite-K calibration would not
automatically certify a precise continuum stress or its sign.

## Publication and classification

Before the first new physical call, publish the distinct followup checkpoint,
exact producer and inherited source pins, full experiment, unchanged gates,
independent manifest and stopping rules, then verify the remote frozen tree.
The followup must explicitly reference the original freeze, warning and
immutable failure record. After execution, report each original/followup
status separately; a followup PASS does not relabel the original FAIL.

A successful followup could validate this homogeneous special-point metric
calibration on prescribed geometry. Actual shifted-root numerical response,
the full metric/state response matrix, bulk gravitational and scalar
constraint matching, dynamical shell displacement, physical EFT validity,
coupled evolution, stability and heating would still require separate work.
