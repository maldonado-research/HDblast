# Conditional realization of the postcontact tensor coefficient

This checkpoint connects the previously enclosed tensor response to the
second derivative of a **specified exact Einstein–scalar family**, conditional
on the retained analytic initial-data, shell-evolution and regularity
theorems. It also defines a covariant functional relative to explicitly
chosen reference data on that family's diffeomorphism orbit.

The numerical interval is inherited from the frozen postcontact calculation:

```
h = J''(0) in [-0.404311163330, -0.014601815851]
at the reference sphere (x,r)=(5/4,1/8).
```

The new result is the identification with an exact family and the complete
reference prescription. It is not a new numerical evaluation of `h`.

Read [the theorem and proof](DERIVATION.md),
[the anchored reference construction](PUBLIC_REFERENCE_LEMMA.md), and
[the independent review](independent/INDEPENDENT_TRANSFER_REVIEW.md).

## What the proof establishes

| Step | Status and scope |
|---|---|
| Analytic existence and regularity | The earlier weighted initial-data and finite-interval shell theorems remain explicit inputs, together with the actual background/master/gauge premises. |
| Exact family selection | Apply the weighted construction anew to the complete current packet, including its velocity and canonical density, on a specified annulus and controlled domain. |
| Tensor identification | The selected curve has the same zero initial TT displacement and velocity. The full angular projection closes independently of a finite angular ansatz, and uniqueness identifies its forced tensor response with the enclosed `U`. |
| Reference prescription | A physical-to-fixed-background wave map, prescribed initial anchors, observer, screen, inverse events, frame, angular templates and measure define one complete functional. Covariance is proved on the known family's diffeomorphism orbit. |
| Nonlinear sign | `J(0)=J'(0)=0`, `J''(0)=h<0`, and `J` is C2. Hence `J(epsilon)<0` for every sufficiently small nonzero parameter of either sign. |

More explicitly, there exists a positive but **uncomputed** `delta` such that

```
0 < |epsilon| < delta  implies
J(epsilon) <= -0.00365045396275 epsilon^2 < 0.
```

This is an existential small-parameter result. No specified nonzero amplitude,
effective parameter radius, or quantitative Taylor remainder is certified.
The quantity is an oriented angular projection of magnetic screen curvature,
not a pointwise scalar curvature or physical force.

The nonlinear completion is essential. The
[completion-dependence example](independent/COMPLETION_SCOPE_COUNTEREXAMPLE.md)
constructs another permitted exact initial curve with the same first
perturbation and different initial tensor data. The theorem concerns the
selected weighted-range curve, rather than every completion of that
quadrupole.

The reference metric, anchors, clock-rate profile, frame, orientation,
templates and measure are part of the functional. Changing them can change
the result. The time label is a wave-map label; it is not proper time, and
the observer need not be geodesic. The proof covers the selected family's
orbit, not an arbitrary nearby metric's map boundary problem. These
distinctions are retained throughout the reference lemma.

## Reproduce

From a complete repository checkout, with Python 3.12.14:

```sh
python -m pip install -r research/HDBLAST_EXACT_FAMILY_TT_TRANSFER_20261031/requirements.txt
python -I -S -B research/HDBLAST_EXACT_FAMILY_TT_TRANSFER_20261031/checks/run_checkpoint.py --repo-root . --output-dir /tmp/hdblast-exact-family-fresh
```

The output directory must not exist and must be outside the source checkout.
The runner authenticates this closed capsule, the inherited source closures
and selected source pins before scientific imports and again afterward.
It replays the new normal/optimized reports and deliberate corruption
controls, then replays the unchanged postcontact checkpoint once. New and
inherited checks are reported separately. No scientific network access is
required after checkout and dependency installation.

The new checks reproduce four report pairs and reject 66 real corrupted
subprocesses. The unchanged postcontact replay separately reproduces nine
report pairs and rejects its 68 corrupted subprocesses. Eight distinct
in-process wrong-formula fixtures are evaluated in the new canonical reports
in both modes; wrapper baselines may repeat them and are not counted as new
fixtures. These are verification counts, not numbers of new theorems.

Algebra, source authentication, reproducible enclosures and internal
AI-assisted review support the stated checks. They do not numerically prove
the inherited analytic theorems or replace external specialist review.
The [primary comparison](literature/PUBLIC_GENERIC_TENSOR_COMPARISON.md)
identifies established generic tensor-perturbation results; this checkpoint
claims a model-specific application and coefficient match, not a new
universal tensor equation or an exhaustive novelty assessment.

The cosmological origin hypothesis, detector interpretation, radiation
energy, stability and a proper-time/geodesic observation remain unestablished.
No new Zenodo version or journal submission accompanies this checkpoint.
Earlier capsules, candidates, manuscripts and websites are preserved.
