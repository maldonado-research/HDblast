# Figure provenance

The left panel compares upper endpoints of two sufficient mathematical ranges.
It does not estimate the maximum range where solutions actually exist. The
right panel plots the rigorous analytic enclosure
`q - K_H*delta <= (H_metric^2-H^2)/delta^2 <= q + K_H*delta`, using the exact
rational constants in the primary certificate. Zero is a limiting coefficient;
the theorem applies to positive detuning only.

`controlled-response.csv` retains exact rational coordinates and separately
labeled display floats. `figure-data.json` authenticates both certificate inputs.
These are mathematical bounds, not observed data or sampled ODE trajectories.
The comparator is the separately defined scalar-independent detuning theory.

To regenerate, install `matplotlib==3.10.1` in addition to the capsule's exact
check dependencies, then run:

```sh
python -I -B figures/make_figure.py --output /tmp/hdblast-response-figure
```

All decisions in the mathematical certificate use exact arithmetic. Rendering
uses floating-point coordinates solely for display.
