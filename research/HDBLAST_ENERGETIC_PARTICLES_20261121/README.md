# What happens to the highest-energy particles?

In a published HDblast simulation, a pulse accelerates particles on an expanding shell. At the last simulated time, many particles still travel close to their original beam direction. Two measurements ask why the most energetic particles can retain that motion. They use the dedicated vector interaction, the recorded shell radius, and the particle source and incident fields from the primary study. [Primary study and reproduction instructions](../HDBLAST_AXIAL_CURRENT_20261121/README.md)

## Which particles does the interaction reach?

At shell time 7.05, the most energetic 10% of particles hold **36.37% of the gas energy**. They account for only **1.19%** of the expected collision participants at that instant. On the same snapshot, giving all pairs the same total cross section raises their share to **10.01%**. Each encounter counts both of its participants. The sums include every one of the 268,435,456 eligible positive-negative pairs in the largest sample.

| Shell time | Energy held by highest-energy tenth | Participants from this tenth, vector resonance | Participants from this tenth, constant cross section |
|---|---:|---:|---:|
| 6.67 | 31.12% | 0.934% | 8.704% |
| 7.05 | 36.37% | 1.192% | 10.010% |

These are instantaneous rate shares. The matched-count box experiment in the primary study follows what happens over the next eight collisions per particle on average. At the late start the initially highest-energy tenth received 2.49 collisions per particle under the resonance and 7.98 with a constant cross section. [Pair-rate calculation and complete results](pairs/README.md) · [Time evolution in the box](../HDBLAST_AXIAL_CURRENT_20261121/freeze/README.md)

## Could motion in the smooth field cool them by emitting radiation?

As a separate test, assume each charged particle emits independently according to the ordinary four-dimensional classical radiation law. Integrating its emission from birth to time 7.05 on the finest simulated trajectories gives energy equal to **1.48 parts in a quadrillion** of the final gas energy. Four sample and time-step choices give about 1.2–2.0 parts per quadrillion. Even an estimate that maximizes the smooth force over particle directions at every recorded midpoint gives less than four parts in ten trillion on these grids. The energy of each emitted contribution includes the appropriate shell expansion factor. The small predicted emission on this prescription gives a reason to look elsewhere for a substantial transfer of energy from these particles. [Derivation, analytic tests and full runs](radiation/README.md)

The emission calculation follows the original particle paths without damping them. Photons emitted during individual collisions, interference between particles and radiation propagating or being reabsorbed in five dimensions require their own calculation. The primary study already accounts for coherent waves radiated by the collective current. Adding the two energy budgets would require specifying which physical modes each counts. The conclusions apply to the prescribed benchmark and do not establish the origin of the hot Big Bang. [Classical radiation formula, section 2.1](https://arxiv.org/abs/1708.04511) · [HDblast coupling normalization](../HDBLAST_NEUTRAL_PARENT_DECAY_20261120/consistency/CONSISTENCY.md)
