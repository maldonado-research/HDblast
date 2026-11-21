# Which particles can still change their energy and direction?

Fresh particles have less time to collide. If too many unscattered particles end up in one energy-and-direction bin, even a perfect redistribution of every particle that did scatter cannot erase that excess. This calculation measures that constraint for the same 20 bins and the same early marker as the two-polarization study. The birth times and momenta come from its homogeneous calculation with both transverse currents. [Results](SURVIVAL.json)

| Prescribed collisionless history and pulse | Collisionless distance | Minimum distance with the registered collision-rate ceiling | Necessary multiple of that ceiling to get within 1 percentage point |
|---|---:|---:|---:|
| Separate interactions, pressure-fitted pulse | 47.16 pp | 9.74 pp | 10.60 |
| Shared interaction, pressure-fitted pulse | 55.70 pp | 15.42 pp | 15.47 |
| Separate interactions, smaller pulse | 61.91 pp | 12.81 pp | 6.11 |
| Shared interaction, smaller pulse | 61.74 pp | 18.50 pp | 12.47 |

The table uses 16,384 birth pairs. The distance is the fraction of particles that must change bins to match the moving thermal reference, expressed in percentage points. The first column is measured. The last two are optimistic constraints. A stronger collision mechanism can also change the wall, density and fields, in which case its history needs to be recomputed. The last column scales the mathematical upper envelope of the registered collision hazard on a fixed history. It does not specify a physical mediator with that strength. [Cross section and density normalization](../rate/README.md)

For a particle born at wall time $s$, let $H(s)$ bound its first collision optical depth through the marker $T$. In the classical binary collision model its probability of avoiding every collision is at least $e^{-H(s)}$. We use the uniform upper limit $m^2\sigma v_{\rm M}\le 1.118569\times 10^{-6}$ for this vector pole. With total daughter rest $D$, radius $a$, channel weight $b$ and action normalization $S_0$, the rate ceiling is this number times $b S_0D/[6(mL)^3a^3]$. The shared channel weight $1/3$ includes all allowed final flavors. The denominator six comes from the three equal flavor populations and their two charge signs. The target density seen by unscattered particles is taken from the homogeneous symmetric background. The rate integral uses the upper endpoint for $D$ and lower endpoint for $a$ on each recorded cell, retaining the time after the particle's own birth.

If particle $i$ would land in bin $j$ in the absence of collisions, let it contribute $e^{-H(s_i)}/N$ to $u_j$. For target bin probabilities $q_j$ the measured total variation must satisfy

$$
d_{\rm TV}\;\geq\;\sum_j\max(u_j-q_j,0).
$$

The reason is simple. The protected mass in a bin must remain there even if all other particles could be reassigned freely. Conversely, an unconstrained rest population can fill target deficits proportionately and attain precisely this bound. The script constructs that assignment as a check. We also solve the inequality with $H$ multiplied by a common factor to obtain the necessary increases in the table. Each input ensemble is independently checked against the predecessor's complete counts before weighting the particles by age. Three separately scrambled birth samples give early TV bounds differing by less than 0.05 percentage points.

The inequality applies to ensemble probabilities. The numerical table uses samples of the collisionless trajectories, fixed production and fields, a fixed measurement time, and a duration short enough for the round, homogeneous model to apply. Coupling the changed particle distribution back to the fields also changes trajectories that have never collided. The [live test](../live/README.md) evaluates that feedback for the separate-interaction model. The [collision cascade bound](../rate/README.md) uses which pairs actually approach the resonance to place a different constraint on fixed-field evolution.

~~~sh
python -I -B survival/cohorts.py --study /path/to/HDBLAST_POLARIZATION_20261121 \
    --particles /path/to/regenerated/ensembles --output /tmp/SURVIVAL.json
~~~
