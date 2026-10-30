# Initial-data point diagnostic

This diagnostic keeps the registered nominal initial functions and all frozen
inputs unchanged. It authenticates the initial producer and report, then uses
the producer's `Context`, which authenticates its seven complete scientific
upstream roots before imports and again after evaluation.

The old rational-packet/cutoff/exponential-approximation channel has a global
radial-derivative H upper cap approximately `1.5334374403e27`. Its dominant
continuous cell is report row 7,

```
263/256 <= r <= 33/32.
```

The scalar-mean `w:0` component dominates that cell, with component cap about
`1.3504368865e27`. The next-ranked old cell is `[143/128,287/256]`, whose
channel radial H upper cap is only about `1.3031971824e23`. This ranking uses
the existing exact rational continuous-cell report.

Nine rational point intervals across the dominant cell are much narrower.
After summing **all three** paired difference channels, the largest reported
point H upper bound is below `4.031e16`. At the specific point `r=1055/1024`,
the complete actual-minus-nominal initial radial scalar-mean coefficient obeys
the strict bounds

```
-3.5495e16 < (u_actual,r-u_nominal,r)_w:0 < -3.5494e16,
4.03e16 < ||u_actual,r-u_nominal,r||_H < 4.031e16.
```

The actual first-field/background inputs and their error bounds remain
conditional inherited premises. The lower bound above concerns the complete
three-channel sum. A lower bound for only the cutoff channel would not suffice:
the other channels can in principle cancel it.

The diagnostic also evaluates `r=573/512`, the midpoint of the next-ranked
old cell. There the summed radial H interval upper bound remains about
`5.96785e22` and its interval-derived lower bound is zero. This illustrates
that the paired rational-exponential approximation enclosure still leaves
unresolved cancellation even at that point.

Point samples do not replace a continuous cell certificate or establish a
global maximum. In particular, they do not prove a continuous `4.031e16`
global cap. Their role is to identify the cutoff-transition interval bound as
an efficient target for a continuous refinement. The observed nonzero
actual-minus-nominal point difference also prevents identifying the two
registered functions exactly. These unsmoothed nominal-data differences are
not a lower bound for the later smoothed reconstruction's actual response
error.

An exact algebraic explanation for potential interval overestimation is the
cutoff polynomial

```
S(t)=126t^5-420t^6+540t^7-315t^8+70t^9,
S(1-t)=1-S(t),
S'(t)=630t^4(1-t)^4.
```

The dominant cell has `3/4<=t<=1`, where a complementary-polynomial evaluation
can avoid cancellation in the alternating power expansion. A continuous
refinement, a proven algebraically identical evaluation, or both must be
checked separately before changing the certified global error bound.

`check_initial_points.py` emits every exact component interval and H cap, plus
lower **squared** norm bounds, in a deterministic JSON report. Run normal and
optimized Python and require byte-identical output. It uses explicit failure
checks and imports no new field solver. This is a numerical diagnostic of
registered coordinate coefficients, not a physical measurement, external
specialist review, nonlinear solution, or evidence for a Big Bang origin.
