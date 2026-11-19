# Can a second burst of particle production bring the expanding shell closer to chemical equilibrium?

HDblast asks whether energy released in a higher-dimensional space could have supplied the hot early universe. The previous study found that the simulated particles could undergo very few elastic collisions before the end of the recorded shell motion. Here we test an interaction that **creates more particles**. We let the shell's scalar field turn the original scalar-to-fermion interaction on a second time. We follow the new particles, the energy that supplies them, and the effect on expansion. [Previous collision study](https://doi.org/10.5281/zenodo.24375564) · [Registered experiment](https://github.com/maldonado-research/HDblast/issues/157#issuecomment-7017450151) · [Archived edition](https://doi.org/10.5281/zenodo.24376497)

**At the illustrative coupling $y=0.1$, any late pulse sequence with total instantaneous strength bounded by the original peak produces at most 22.3% of the particle count required for chemical equilibrium at the energy actually produced.** This bound includes the original particles, interference between production intervals, and arbitrarily high daughter frequencies in the specified leading-order massless Dirac model. It holds on the four prescribed shell histories, over the measured time remaining after the first pulse, with the original quantum state and its 24 modeled parent momenta. A count ratio of 100% is one necessary condition for equilibrium. The shape of the particle distribution and the rates of its reactions also matter. [Bound, assumptions and numerical controls](criteria/CRITERIA.md)

![On the quiet shell, the second pulse raises energy faster than particle count as successively higher frequencies are included. The original count is 12.23% of its equilibrium comparison, the two-pulse count at cutoff 1920 is 8.45%, and the quiet-shell ceiling for any allowed late pulse is 22.21%.](renewed.png)

## Why producing more particles can increase the shortfall

A massless gas in chemical equilibrium with zero chemical potential has particle number proportional to the three-quarter power of its energy at a given volume. More energetic particles raise the equilibrium count needed at the new energy. For example, putting the same energy into many soft particles makes more progress toward the count target than putting it into a few energetic ones. We compare a single massless Dirac species, with equal particle and antiparticle populations, on identical comoving volumes.

We kept the original smooth pulse of the local scalar-clock coupling. The second unit-height pulse occupies quiet-clock coordinates 2.055 to 2.145. Every wall uses the same function of the scalar field, including the walls perturbed by incoming disturbances. We end each calculation at shell proper time 2.16. The quantum calculation retains the relative phase of the two emission amplitudes. The table shows quiet-shell coefficients with the common factor $y^2$ removed from number and energy. The last column restores $y=0.1$.

| Source and maximum pair frequency | Pair count | Energy of both daughters, conformal units | Count as a share of chemical equilibrium |
|---|---:|---:|---:|
| Original pulse, 120 | 0.76982 | 27.81645 | 12.23% |
| Both pulses, 120 | 0.97628 | 44.85620 | 10.84% |
| Both pulses, 1920 | 1.20697 | 82.90746 | 8.45% |

At cutoff 120 the pair count increased by **26.8%** and the energy increased by **61.3%**. The cutoff also matters for short pulses. It includes only **30.9%** of the second pulse's energy found by extending to 1920, on the same quiet history. Changing the sign of the second pulse barely changes the total count, yet the two normalized single-daughter momentum distributions differ by 0.159 in total variation on the archived bins at cutoff 120. This is a phase-sensitive prediction for the model. [Pair amplitudes, spectra and frequency controls](emission/EMISSION.md)

## What rules out a better-timed burst?

We derived an upper limit using the inverse third moment of the pair mean energy, the available time and the parent field variance. Hölder's inequality bounds the achievable particle number at fixed energy. The triangle inequality allows the most favorable interference between the original and renewed production. A Fourier norm bound controls the entire high-frequency tail of the specified source. This combination covers an arbitrary number of late pulses with bounded combined strength. Sharp pulses that would inject infinite energy do not define a gas at finite energy.

The quiet upper limit is **22.21%**. The highest limit on the other three shells is **22.28%**. We also tested a different constraint that allows very large peaks provided the integrated squared coupling, weighted by the parent variance, stays within the budget of continuously running the unit-strength source over the rest of the interval. The largest ceiling in that case is **22.51%**. The squared-source budget measures driving strength. The energy transferred by the scalar clock is a separate, signed observable.

As a finite test of whether useful timing was overlooked, we varied two shorter pulses independently over positive and negative strengths. A rational polynomial calculation certifies the best achievable score for the stored spectral matrices. The best relative improvement over the original pulse was below 0.16 parts per million at cutoff 120 and below 0.003 parts per million at 1920. These are computational comparisons for this two-parameter family. The waveform-independent inequality gives the broader physical restriction. [Two-pulse optimization](family/FAMILY.md) · [Exploratory search protocol](https://github.com/maldonado-research/HDblast/issues/157#issuecomment-7017538828)

## Where the energy came from

In the native matter model with the archived static local correction and cutoff 120, the quiet shell gains **13.648** units of proper daughter energy per normalized comoving volume from the renewal. The change in the parent contribution is **−3.922**. Together the changes give **+9.726**, accounted for by **+9.971** of scalar-clock work, **−0.051** of mass-dependent work and **−0.194** of work from expansion. All figures are leading coefficients in $y^2$. Removing the local correction or changing the cutoff changes the separate work and parent terms. The code reports both conventions and checks the energy balance. [Matter, pressure, force and work](budget/BUDGET.md)

As a control, we froze the geometry and parent frequencies and applied the pulse to an initial vacuum. In that case the switch must pay the full energy of the emitted pair **and** the parent excitation. On the quiet shell it supplies **4.030** native units for the full second window. Halving the pulse duration at the same onset and peak raises that cost to **38.564**. The control has an analytic Fourier transform and a bound on its missing high-frequency work. [Stationary vacuum and excited-parent controls](stationary/STATIONARY.md)

## What gravity adds, and what remains open

Eighteen five-dimensional gravity evolutions test the finite-cutoff source directly. At the final time, renewal offsets **37.5%** of the original coupling's slowing of the expansion rate on the quiet fine grid and **34.3%** on the near-negative fine grid. The combined response remains negative on both at that endpoint. All compared trajectories agree before the second interaction starts, and an independent shell equation checks the endpoint differences. The gravity calculation uses cutoff 120. The source's higher-frequency energy and its cutoff-dependent work prevent extrapolating that gravitational response to infinite frequency. [Native gravity results](gravity/GRAVITY.md)

The mathematical restriction narrows what a reheating model can rely on: this short renewal of the existing interaction cannot supply chemical balance under the stated conditions. Changing the source population, spectrum, interaction or time available could change the result. A candidate in the D-Blast 3 archive involves an oscillating scalar after the initial release. Testing that candidate would require evolving its energy supply and motion with particle production, then calculating whether a hot, radiation-dominated phase and a distinguishing observable follow. [Primary literature and project archive](literature/README.md)

## Reproduce the study

From the study directory in the repository or in the standalone Zenodo ZIP, run the commands below. The first authenticates the study and 20 pinned historical inputs. The second regenerates the source spectra, work accounts, stationary controls, mathematical bounds, optimized pulse family and 18 native gravity trajectories. The reference runtime uses Python 3.12, NumPy 2.3.5 and SciPy 1.17.0. Matplotlib is needed only for the optional figure command.

~~~sh
python -I -B -S reproduce.py
python -I -B reproduce.py --verify --jobs 2
python -I -B figure.py --output /tmp/hdblast-renewed-pairs.png
~~~

[Historical input manifest](SOURCES.json) · [File checksums](SHA256SUMS)
