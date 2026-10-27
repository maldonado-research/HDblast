# Conditional clock anchor and shell-lift regularity

This checkpoint advances two prerequisites for interpreting and testing the
second-order response. It preserves the earlier screen-curvature, observer-frame
and full-angular-response capsules unchanged.

## Event identification and gauge

[The clock and event construction](anchor/ANCHOR.md) starts with additional
marked shell data, propagates a shell clock and angular labels along a
nonvanishing spacelike scalar-gradient flow, and uses five scalar labels to
identify a local event. It requires a regular, transverse, caustic-free tube,
a timelike clock gradient and a nonzero background scalar derivative; none of
these is established here for the perturbed shell-to-target region. An inherited
conditional background enclosure does give `|phi0,y| >= 3/25000` along the
registered target's background shell-to-target segment. The marked data must
be transformed together with the metric and scalar field.

[The gauge and moving-event argument](anchor/GAUGE.md) combines the established
second-order pullback rule with the model-specific first-Weyl structure from
[the observer-frame proof](../HDBLAST_OBSERVER_FRAME_STRUCTURE_20261026/).
Assuming the same background frame, event, orientation and smoothness, the
leading screen-curvature coefficient is unchanged by the corresponding smooth
perturbative identification change and first-order movement of the target.
An independently rebuilt smooth frame with that same background frame gives the
same leading components for this specific registered first metric. A sphere-wide
projection additionally requires consistently matched charts, frames, basis and
measure over the whole sphere. If all these geometric conditions and the
original response/error assumptions hold, the prior bound carries over; no
new numerical source or response calculation is made here.

This is not a construction of a preferred observer from the symmetric background,
an observed signal, arbitrary-background-frame independence, or a proof of a
finite-amplitude exact solution. H1=0 alone is not enough for independently
rebuilt frames: the controls include an explicit local metric counterexample.

## Shell regularity

[The shell-lift note](shell/NOTE.md) derives the C² characteristic join for
the frozen lift `Lg=-integral_0^((x-r)_+) g(s) ds`. For C¹ boundary source g,
the existing g(0)=0 gives a C¹ join; g'(0)=0 is also necessary and sufficient
for a C² join. Its second derivatives depend on g'. A source family with a
zero initial collar makes the lift error tend to zero in C¹ while a second
derivative remains nonzero on the retained domain. Thus a C¹ error certificate
alone does not bound a pointwise second-derivative residual. The differentiated
source includes both shell feedback and the affine boundary source.

This is a local result for the lift; it does not calculate the actual derivative
error, solve the coupled problem, or certify the pointwise Einstein residual.
The existing full-angular construction includes a **conditional formal**
second-order completion and subsidiary premises; this checkpoint does not
replace or independently verify them.

## Reproduce the algebra checks

Python 3 and its standard library are sufficient. From this directory run:

```sh
python checks/check_anchor_composition.py
python checks/check_label_jacobian.py
python checks/check_inherited_background.py
python checks/check_gauge_counterexample.py
python shell/checks.py
```

Each program uses explicit failure checks and may also be run with `python -O`.
The first exercises both gauge generators, inverse event jets, and an unmatched
moving-event control. The second checks the five-label Jacobian and its
nonzero-scalar-gradient condition. The inherited-background check authenticates the cited input bytes and checks
domain containment; it does not re-prove its interval bounds. The next check
constructs a generic local metric
with vanishing first screen curvature but nonzero change under an independently
rebuilt frame; it is a counterexample to a broader claim, not an HDBLAST
solution. The final program tests the lift identities, corner conditions,
finite-domain source family, and physical-coordinate chain rule with incorrect
formula controls. These are exact toy or algebra tests, not evaluations of the
actual field family. Internal AI-assisted reviews are not external specialist
or proof-assistant review.

## Related literature and limits

The second-order identification formula is established mathematics; see
[Bruni et al.](https://arxiv.org/abs/gr-qc/9609040), especially section 4.
[Mars](https://arxiv.org/abs/gr-qc/0507005) develops first- and second-order
moving-hypersurface perturbations, including ambient and intrinsic gauge
freedoms; its framework is relevant to a future complete shell calculation.
[Fröb and Lima](https://arxiv.org/abs/2303.16218) discuss relational coordinate
construction in cosmological spacetimes. Their free-falling synchronous
observers are not automatically the accelerated observer used here. These
focused comparisons do not establish novelty or solve this model.

Next, establish a perturbed scalar-gradient bound and a regular tube for the
specific target, derive numerical derivative bounds for the total shell source
and bulk terms, and evaluate the complete compatible second-order residual.
None of those results or a higher-dimensional origin of the Big Bang is
established by this checkpoint.
