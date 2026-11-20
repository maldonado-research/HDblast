# How does the birth-state count change at larger mass?

These are conditional tests of two **cycle-averaged** parent trajectories. The resolved drive, the subsequent averaged parent, the field, matter and geometry evolve together according to the envelope model. Both examples use $Q=10^7$ and $L\nu=18$. Larger $mL$ opens magnetic branches that were inaccessible at $mL=400$. The code counts candidate quantum states before asking how many particles the fluid puts into them. [Averaged histories](../envelope/ENVELOPE.json) · [Lower-mass test and conventions](KINETICS.md)

| Daughter mass $mL$ | Hot radius factor | Required action for hot drag throughout | Largest action allowed by collisionless count in the tightest sampled bin | Largest action in that bin with label clearance capped at $18/L$ |
|---|---:|---:|---:|---:|
| 4e+10 | 2.2528 | 2.08e+47 | 3.22e+49 | 4.78e+50 |
| 4e+08 | 2.2529 | 2.08e+41 | 3.22e+41 | 4.78e+42 |

The windows divide the 1st to 99th percentiles of produced daughter number using exactly the same native-step measurement and padding as the lower-mass study. For each interval let $[b_0,b_1]$ enclose the absolute magnetic shift and let $k_+$ exceed the largest sharp birth momentum in unit-sphere coordinates. In each interior $(\ell,j)$ block the necessary condition for a birth is
$$
(\ell-|b|)^2+4j|b|\le (k_++1/2)^2 .
$$
The code allows both positive-energy spin states in a block whenever either could be accessed. This allows intrablock mixing. It encloses all its possible $\ell$ values by expanding $[b_0,b_1]$ on both sides by $\sqrt{(k_++1/2)^2-4j b_0}$ and imposing $\ell\ge j+1$. Both unpaired extremal branches are counted separately. The bound grants every candidate state once in the collisionless test. A density-dependent collision frequency, broad formation bandwidth or additional nonperturbative sources would need their own analysis.

Bounds on the field and momentum are rounded outward to integers. The exact integer square root and triangular sums handle large quantum numbers directly. At most 4,096 groups of successive $j$ are needed, since the enclosing count per $j$ is nonincreasing. The first and last $j$ in each group enclose its sum. The lower sum checks numerical accuracy for this enclosing set; it is not a lower bound on the states that a scalar matrix element actually populates. For context the corresponding continuum state count for a momentum ball swept along a straight interval is $4k_+^3+3k_+^2(b_1-b_0)$ over both signs and all three flavors. This follows by multiplying the capsule volume in momentum space by the unit-sphere volume, two spins per species and $(2\pi)^{-3}$. It provides an independent normalization check at large branch count.

At $mL=4\times10^8$ the tightest collisionless window delivers 75% to 90% of the final number. Its enclosing count includes paired branch indices up to 694. The demand at the pointwise hot-drag floor is 0.6467 times the upper state count. The summed-count upper/continuum ratio there is 1.00039. Across all six windows the largest fractional width of the grouped-sum check is 1.1e-09. Its combined gauge/drag lower bound is 2.08e+41 and the mass-shift cap is 6.25e+39. Imposing the limited-clearance bound as well leaves no common action interval for these specified tests.

At $mL=4\times10^{10}$ the tightest collisionless window delivers 75% to 90% of the final number. Its enclosing count includes paired branch indices up to 69,433. The demand at the pointwise hot-drag floor is 0.006469 times the upper state count. The summed-count upper/continuum ratio there is 1.00015. Across all six windows the largest fractional width of the grouped-sum check is 2.2e-06. At the separate gauge, mass-shift and drag screen's example action $S_0=3.61e+48$, demand in that bin is 0.1121 times its collisionless upper count. Its combined gauge/drag lower bound is 2.08e+47 and the mass-shift cap is 6.25e+49. The three screening inequalities and all six measured collisionless upper counts admit the interval 2.08e+47 to 3.22e+49.

Passing an upper-count test only shows that this counting argument permits the specified yield. Angular and spin selection rules, the true resonance width, detailed balance, Pauli occupancy of individual modes, and a collision law in the evolving distribution can still restrict production. The shell widths, actions, charges and masses have no independent measured calibration here.

From the capsule root, run

```sh
python -I -B kinetics/many.py --envelope envelope/ENVELOPE.json --output /tmp/MANY.json
```

[Calculation](many.py) · [All bins and assumptions](MANY.json)
