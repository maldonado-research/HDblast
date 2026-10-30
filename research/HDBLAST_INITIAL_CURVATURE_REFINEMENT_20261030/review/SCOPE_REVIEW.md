# Scientific scope and transfer review

30 October 2026. This is an internal review of the initial-curvature checkpoint,
not an external specialist assessment, originality finding or journal decision.

## Mathematical assessment

The initial-join argument separates the actual target-minus-candidate error
from a nominal curvature jump. The nominal initial displacement is C1, so its
distributional second derivative has no point mass. The even positive kernel
has half mass 1/2 and half first moment 35 delta/256. Integrating the sidewise
third-derivative bounds therefore gives the stated trace-average enclosure.
The two inner kernel footprints contain no other join or outer continuation.
The actual analytic packet is separately evaluated with signed derivatives;
this is essential at the first join, whose nominal curvature traces are zero.
Intersections of two valid enclosures preserve the same actual target.

The pointwise H lower bound is valid: every component lower bound holds for
the same actual vector, and the diagonal Gram weights are positive. Correlation
among interval coordinates does not weaken this sum-of-squares lower bound.
These are initial-point coordinate errors of the fixed reconstruction. They
are not physical forces or obstructions to the underlying nonlinear model.

The continuous refinement keeps the same actual and nominal functions. Ten
exact dyadic leaves replace the single original cutoff interval, while 47
authenticated old cells remain. For each telescoping channel, the maximum of
the cellwise H caps is a uniform cap; summing the three maxima is justified by
the triangle inequality. There is no requirement that the maxima occur at the
same radial point. This aggregation can be sharper than packing independent
global component maxima. The complete retained smoothing, inner-jump and P1
velocity charges remain necessary and are retained.

Replacing only the initial error cap in the PR102 differentiated theorem is
sound. In particular, the absolute initial curvature of the candidate remains
inside the differentiated forcing ledger; it is not replaced by the smaller
error cap. The source, coefficient, physical boundary inverse, C1 error and
intermediate correction estimates are unchanged. The nonzero initial
accelerations of the corrections and the original identity r1(0)=0 remain in
the calculation.

## Independent numerical checks in this review

The exact refinement report with SHA256
`915acf8abf3c3bb2c843eac509315ebdbd568d7c8b10ff2682c900d672e7d4d7`
gives an initial radial-curvature error upper cap approximately
`1.1845751478843721e28`, strictly below `1.185e28`. Its old/new upper-cap ratio
is approximately 524990.3296.

The transfer report with SHA256
`fa13efe588fa19ab6a09c517b67bc277cd3dcaeee9d428143bf479ef98128bb0`
gives caps approximately `1.1865241670e28` for xx and xr and
`1.1865284999e28` for rr. All are strictly below `1.187e28`; all old/new cap
ratios exceed 524000. These ratios compare two upper certificates for the
same true error. They do not demonstrate a decrease of the actual error.

The independent delta-identity checker was run in isolated normal and optimized
Python. Both returned the same 29-group PASS report, SHA256
`adc9e9342a2fff1f8b5e198e58b1c7d593a48eac58a00684bde6211e5595029c`.
It imports no producer or transfer checker code and verifies the exact identity

```
new second-error cap = old second-error cap + (new J - old J)/(1-q)
```

for all three derivatives, rebuilding the initial acceleration and max branch
independently. It also checks that all other initial, source, boundary, C1 and
correction derivative ledgers are retained exactly.

The initial contribution's share of the differentiated upper majorant is now
approximately 0.9995988191. The old statement `>0.999999999` applies to the
previous coarser bound and must not be repeated for this refined ledger. Neither
share is a decomposition or lower bound on the unknown true error.

## Publication limits

There is no mathematical blocker in the reviewed join, subdivision or monotone
transfer arguments. Acceptance remains conditional on the authenticated model,
initial-data, actual-source, regularity and physical-boundary premises of PR102.
The new numerical bounds are large absolute scaled-coordinate bounds. No
nonzero physical reference scale or demonstrated relative accuracy follows.
The work does not independently complete all physical initial or subsidiary
constraints, construct a finite-amplitude nonlinear family, or provide a
cosmological observable or evidence for a higher-dimensional Big Bang origin.

A higher-order cutoff would define a new candidate and would require new
initial and transfer estimates. The current reconstruction, original sources,
point tables and frozen public manuscripts are preserved by this checkpoint.


## Final staged publication review

The final capsule README and the new top-level README paragraph were reviewed
against all seven staged reference reports and the three derivation notes.
The published numerical thresholds, upper-bound reduction wording, unchanged
candidate/source statements, initial-point lower bounds and physical limits
agree with the reviewed reports. The text does not claim actual-error reduction,
physical relative accuracy, external certification, a nonlinear family or Big
Bang evidence. It does not repeat the old initial-majorant fraction for the
new ledger.

Every relative documentation link in the capsule README resolves to a staged
file. The bounded methods/literature note distinguishes full text from
abstract-only review, identifies the versioned primary sources, and states that
no inspected nonlinear Dirichlet-wall theorem is imported for the material-shell
boundary law. The manuscript and private editorial history are outside the
new public capsule. Earlier README checkpoint paragraphs remain historical
entries beneath the new update.

All seven staged reference report statuses are PASS. In particular, the final
join diagnostic report is
`e1d0d6e3d8650cbf402ccb4f5a87fd20c20db651e3bd323c68f0bdea0219dd35`,
its independent ODE audit is
`c6b8bacc8b2cf36d20bff6625d520bd38a631dfe48f7745d95c070160086d823`,
and the complete-domain initial audit is
`ed4204f9b8ce39b0cebbe75fd8fc6ffd88b4977d8810a382ebcf8061ba0c4364`.
The refined initial and C2 report pins above remain unchanged. The displayed
join H lower bounds were checked directly against the exact report values:
`0.6942779566... >0.6942` and `5.8967317546e20 >5.896e20`.

Final scientific/publication scope assessment: **PASS within the stated
conditional mathematical scope.** Full checkpoint replay, GitHub CI, merge and
remote byte readback remain separate publication steps; this review does not
assert their completion.
