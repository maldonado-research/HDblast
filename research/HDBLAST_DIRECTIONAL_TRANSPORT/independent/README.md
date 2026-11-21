# Independent check of the angular deficit

This calculation uses the pinned HDblast trajectory generator, source normalization, shell derivative, electric acceleration and birth-window cutoffs. It integrates the directional deficit in the reverse order from a calculation that averages over direction first. Shared archived kinematics make it a check of the new phase-space calculation on the specified numerical trajectory.

Divide all speeds by the required source per shell area. For a momentum direction with cosine x, the instantaneous missing fraction is max(0, 1 − |c − d x|). The electric acceleration d varies linearly between the two endpoints of an accepted wall step. At fixed x the radial velocity then ranges uniformly between c − d₀x and c − d₁x. Its time-averaged missing fraction is the area under a unit-height triangular function on this velocity interval divided by the interval width. The code computes that area exactly with clipped trapezoids, including the constant-velocity limit.

It then integrates this time average over x with eight-point Gauss–Legendre quadrature. The integration splits at x = 0 and at every angle where an endpoint velocity is −1, 0 or 1. On each resulting interval, the integrated tent has terms proportional to x, a constant and 1/x. Geometric panels bound the ratio of the nonzero angular endpoints by two. Three simple controls test a static field, a field reversing sign and a vanishing field. A fourth checks a source with all angles below capacity. [Quadrature and controls](reverse.hpp)

Run the following from this directory with Python 3 and a C++17 compiler. Supply an absolute or relative path to the archived HDblast coherent-decay study. The program verifies the input digests and recompiles the two published mass histories. It accumulates the missing fraction with the published comoving birth weights for both time steps and both action choices. The optional comparison also evaluates each reported peak state by the trapezoid rule split at the angular breakpoints, where the instantaneous deficit is piecewise linear. [Pinned HDblast source](https://github.com/maldonado-research/HDblast/tree/dc99ff5c5894abcff3e11e02e4b97e19da31de00/research/HDBLAST_COHERENT_DECAY_20261120)

```sh
python -I -B check.py --source /path/to/HDBLAST_COHERENT_DECAY_20261120 --compare ../replay/results.json
```

For the half-step runs and the central 99.998% of the births, the fractions of the total requested number missing at the published action choices are 0.008204183312 for mass 4 × 10¹⁰ and 0.00002899262089 for mass 4 × 10¹². At the lower action required by the reference drag these fractions are 0.000558575797 and 0.00000579363032, respectively.

The field, shell and requested source are prescribed at each step. Occupied initial states, coherence over multiple crossings, collisions, finite spectral width and backreaction of a modified source require further input. [Runner](check.py) · [Results](CHECK.json)
