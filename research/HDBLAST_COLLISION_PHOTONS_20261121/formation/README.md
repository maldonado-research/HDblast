# Does the photon have time to form between successive encounters?

The ordinary soft-photon calculation treats the motion on both sides of a collision as straight. Actual particles on this shell can collide again or bend in the electric and magnetic fields. We ask how much of the predicted photon energy lies in portions of the selected bands for which the simulation supplies well-separated, nearly inertial legs. The comparison also measures how much kinetic energy is available repeatedly to the registered encounters. [Collision records](../events/) · [Soft-photon calculation](../soft/)

## Time and energy in the emitting pair's frame

For each equal-mass elastic event define $s=(p_++p_-)^2/m^2$, $\gamma_*=\sqrt{s}/2$, $\beta_*^2=1-4/s$. The photon energy in the incoming pair's center-of-momentum frame is $\omega_*/m=\xi\gamma_*\beta_*^2$. The kinematic endpoint corresponds to $\xi=1$. If a particle has four-velocity $u$, the wave phase accumulates at proper rate $k\cdot u$. For a straight leg, its lab time for one radian of phase is $T=\gamma_{\rm wall}/(k\cdot u)$. The smallest possible $k\cdot u$ in the pair frame occurs for parallel momenta. Thus, using the shell time unit $L$ and $\mu=mL$,

$$
 T/L\ \leq\ \frac{\gamma_{\rm wall}(1+\beta_*)}{\mu\xi\beta_*^2}.
$$

This upper limit is valid for every photon direction on all four Born legs. It is conservative because different legs generally attain their individual maxima at different directions. For each incoming track we use the physical time since birth or its preceding encounter. For each outgoing track we use the time to its next encounter or to the end of the simulation. The last interval is censored at the endpoint. Dividing the right-hand side by each available time gives the minimum $\xi$ that passes a chosen spacing check. The archived Poisson times fall within the integrator steps. The momentum push and elastic update use the parent's operator splitting. [Event times and momenta](../events/README.md)

## Bending of the legs

As a local diagnostic we evaluate the magnitude of the electromagnetic proper acceleration $m\chi$ on each leg using its particle momentum and the interpolated wall field at the encounter. If that proper acceleration held fixed, the rapidity accumulated in one formation time would be at most $\chi(1+\beta_*)/(\xi\beta_*^2)$. A separate direction-independent estimate uses only the recorded fields throughout the corresponding straight-leg formation interval. The electromagnetic contribution to the rapidity path length is bounded by

$$
\int m\chi\,dt_{\rm proper}\ \leq\ |Q|\int (|\mathbf E|+|\mathbf B|)\,d(t_{\rm wall}/L).
$$

For each field cell we use the larger of its two endpoint values of $|\mathbf E|+|\mathbf B|$. Convexity bounds the linearly interpolated fields, including on partial cells. We invert the cumulative integral to find the longest available interval with a rapidity bound of 0.1. Combining this with the one-radian time inequality gives a second minimum $\xi$. The choices of ten radians of separation and 0.1 in rapidity express the adopted margins. The calculation archives a one-radian comparison and a local acceleration comparison as well. It gives a consistency test within a locally flat description of the recorded fields. Curvature-dependent photon propagation, dispersive effects, coherent emission across encounters and interference with radiation from the mean current require additional physics. [Normalization of the local force](https://github.com/maldonado-research/HDblast/blob/9467ffbc67e9540da6762fecd1c44e5e22c80522/research/HDBLAST_ENERGETIC_PARTICLES_20261121/radiation/README.md)

## Results

Each entry below is a fraction of the vacuum leading-soft **energy estimate** for the specified band at time 7.05. An event contributes according to its full four-leg interference coefficient, incoming shell energy and expansion weight. At leading order the energy per unit $\xi$ is constant within each event. The calculation integrates the portion of the window above the required lower cutoff. The last row extends into a region where the finite-energy calculation gives a separate test of the soft approximation.

| Photon energy as a fraction of the kinematic endpoint | Ten radians before another event or endpoint | Also local bending below 0.1 | Also integrated, direction-independent bending envelope below 0.1 |
|---|---:|---:|---:|
| 0.001%–0.01% | 83.32% | 0.000040% | 0% |
| 0.01%–0.1% | 98.12% | 1.35% | 0.033% |
| 0.1%–1% | 99.84% | 59.76% | 23.50% |
| 1%–10% | 99.98% | 97.11% | 78.04% |

For the 0.1%–1% band the integrated-envelope fraction is 22.57%–24.42% across the four late realizations. A failed margin leaves the emitted amount unresolved for that portion. Passing a margin gives geometric separation on these trajectories. These tests do not calculate the correction to the spectrum or imply that the failed portions are suppressed by the fractions in the table. At time 6.67 the same tests in the 0.1%–1% band give 99.87%, 18.35% and 0.66% for the finest run. [All eight event histories](FORMATION.json)

## Hard photons with recoil

More than 83% of the isolated-collision energy estimate at the late endpoint is carried by photons above a quarter of the kinematic maximum. We also test this hard band directly, with both recoil and the finite-energy cross section. For photon momentum $k$, the absolute virtuality of a radiating external line relative to its mass shell is $2p\cdot k$. Its one-radian time in a chosen shell frame is the inverse equal-momentum energy mismatch. For radiation before the scattering this is

$$
\frac{\sqrt{m^2+|\mathbf p-\mathbf k|^2}+E_p-\omega}{2p\cdot k}\leq\frac{E_p}{p\cdot k}.
$$

For radiation after scattering, replacing the minus signs with plus signs gives an upper limit $(E_p+\omega)/(p\cdot k)\leq E_{\rm pair}/(p\cdot k)$. Here $p$ denotes the real outgoing charge momentum. Each outgoing massive particle in the three-body final state has pair-frame energy at most $m\sqrt{s}/2$. Thus the same lower bound on $p\cdot k$ used above holds for all its allowed momenta after recoil. Use $\xi=0.25$ in these upper limits and the limits cover the entire harder band. On the outgoing side we use the shorter remaining time of the two original elastic tracks. These future intervals are a surrogate since a photon that takes substantial energy would change subsequent collisions.

Weighting each passing event by the hard-photon energy from the finite-energy calculation gives **99.99%** of the hard-band energy on events meeting the ten-radian time requirement, and **93.84%** on events also meeting the integrated rapidity margin on both sides. Across the four late realizations the latter fraction ranges from 93.23% to 93.84%. The result measures how much of the model prediction passes conservative time and force checks on the recorded history. It cannot bound the true radiation yield of the surrounding charged medium. [Finite-energy spectrum and weights](../hard/)

## Energy a stronger channel would have to carry

Sum the incident energy of both collision partners, with their numerical weight $w$ and the expansion factor $a$, and divide by the final conformal particle energy. On the finest late run this gives **4.6593**. The same energy can pass through several collisions, so the sum can exceed one. If an additional channel on this fixed event history were to take away conformal energy equal to a tenth of the final gas energy, it would have to remove an average of **2.15%** of the incoming pair energy per encounter, weighted by that energy. At 6.67 the corresponding target is **9.76%**. The tenth is an illustrative energy budget. Strong enough new reactions would change the event history, and accounting for emerging radiation would require propagation and possible reabsorption. [Data and implementation](measure.py)

From the parent study directory, after generating the tapes in EVENTS, run the following command. The same soft-photon coefficient is used by both analyses.

~~~sh
python -I -B formation/measure.py --events EVENTS --hard hard/HARD.json --output FORMATION.json
~~~
