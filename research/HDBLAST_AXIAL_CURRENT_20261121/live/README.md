# What happens when the particles supply the axial current?

In the earlier benchmark the axial electric field followed a charged fluid, while individual particles felt that field and scattered. A fluid and a set of scattered particles can carry different currents. Here the same individual particles that gain energy supply the electric current in all three directions. Through the early part of the pulse both prescriptions give similar results. Later they predict different reversal times, different amounts of energy transferred to the gas and different pressures. This is a calculation on the dedicated-interaction wall history. [Primary runs](LIVE.json) · [Second reversal and later snapshots](CROSSINGS.json)

## Encounters before and after the field reverses

Each row follows the same source history from the onset of particle production to the stated wall proper time. The finest samples have 16,384 equal-weight charge pairs. The work percentage is an independent clock taken from the original fluid calculation. Contact counts refer to accepted elastic encounters per particle, including particles with zero encounters. The percentage ever hit is evaluated among particles born by the endpoint. Later snapshots include newly born particles. Samples for different endpoints draw independent representations of the full birth history, so an individual particle should not be followed across rows.

| Wall time | Reference work | Mean contacts | Ever hit | Largest eigenpressure / gas energy density | Distance from common thermal gas |
|---:|---:|---:|---:|---:|---:|
| 6.27855 | 1.0% | 0.074 | 7.0% | 0.996 | 49.5% |
| 6.30944 | 5.0% | 0.134 | 12.4% | 0.996 | 49.5% |
| 6.45 | 33.7% | 0.451 | 35.0% | 0.991 | 47.8% |
| 6.67 | 69.4% | 3.114 | 82.2% | 0.936 | 43.5% |
| 6.73 | 78.1% | 3.882 | 85.5% | 0.915 | 42.2% |
| 6.90 | 92.1% | 6.982 | 95.8% | 0.830 | 35.7% |
| 7.05 | 97.9% | 9.453 | 88.6% | 0.764 | 30.3% |

Pressure eigenvalues include both charges of one representative beam, measured in its zero-energy-flux frame. For an isotropic gas of highly relativistic particles each pressure approaches one third of its energy density. At time 7.05 the three computed ratios are 0.1164, 0.1194 and 0.7641. An observer looking at three equal rotated copies requires a different angular diagnostic. [Three-beam diagnostics](../triad/)

The thermal comparison has five energy bins and four angle bins relative to the wall's axial direction. The equilibrium temperature and bulk motion are constrained by the final mean energy and momentum of the measured population, using the pinned common relativistic gas solver. The statistic is half the sum of absolute differences in bin fractions. Thus at time 7.05 at least about 30% of the sampled population would need to move between these cells to reproduce their target probabilities. Comparing the charge-separated streams to the original drifting-fluid target instead gives 32.6%. Neither finite-bin comparison resolves variations within the cells. The spread in the common-gas statistic at 7.05 is 29.46–31.13% across the four tested sample and time-grid choices. [Histogram and equilibrium convention](https://github.com/maldonado-research/HDblast/tree/main/research/HDBLAST_SCATTERING_20261121/kinetic)

## Which particles receive the encounters?

At time 7.05, the fine run counts 9.45 accepted contacts per particle on average. Giving each particle weight in proportion to its final energy lowers that average to **7.43**. Weighting by its contribution to axial pressure lowers it to **6.28**. The 44.4% of particles with fewer than eight contacts carry 58.8% of the measured energy and 66.6% of the axial pressure. The most energetic quarter alone carries 72.5% of the axial pressure, yet it averages only 5.79 contacts. The axial weight here is the positive wall-frame contribution $p_z^2/\gamma$. These shares are calculated before combining rotated flavors. The final gamma includes rest energy.

Birth time matters as well. In the last birth-time quarter 39.6% never collide, whereas the two middle quarters average 12.73 and 11.44 contacts. These measured differences help explain why a population can exceed the earlier eight-contact box benchmark and still carry most of its pressure along one axis. Determining the separate causal effects of the drive and the energy dependence of the cross section requires controlled comparisons. The tagged replay reads counters already accumulated by the collision solver, reproduces its event and energy totals, and reports four sample or time-step settings at both 6.67 and 7.05. [Contact and birth cohorts](COHORTS.json) · [Tagged calculation](cohorts.py)

## Work and reversal times

The original axial electric field changes sign at times 6.47163 and 6.66085. The fine particle-current run through 7.05 puts these changes at 6.48246 and 6.72190. The second particle-current crossing spans about 6.718–6.723 across refinements. The original axial current reverses at 6.57481. At the 6.67 endpoint the measured particle-current reversal is near 6.581. In the final 0.01 time units near 7.05, the axial electric field's root-mean-square amplitude is of order one ten-thousandth of its value at first birth. Its late zero timings vary across particle realizations. The first two reversals are resolved more reliably.

At time 7.05, 613.388 conformal energy units were already stored in the axial field when particle creation began. About 0.00009 additional units entered through the wall's inward characteristic. The axial field transferred 501.140 to the simulated particles, emitted 111.990 and kept 0.259. The matter ledger adds transverse work, the birth energy and the work of expansion. Its final gas energy is 501.553. With the imposed fluid current, the fine run instead credits 681.580 units of axial work to the test particles while subtracting only 448.241 from the axial field for work on the fluid. The difference of 233.339 quantifies the energy error incurred by interpreting both calculations as one gas-field system.

## Equations and scope

Let $F=a^2E$, $G=a^2B$, $i=a^2 I_{\rm wall}$, and $d\eta=dt/a$, with $a(t)$ the wall radius. Each of the three homogeneous curl-2 components satisfies

$$
F'=-2G-\frac{2a}{K}(F-2i)-\frac{J}{K},\qquad G'=2F,
\qquad J=Q\sum_n w_n s_n v_n.
$$

A prime denotes a conformal-time derivative. Here $K=2.54$, $Q=10^7$, $s_n=\pm1$ and the sum of the weights of both signs is the accepted comoving rest population $D(t)$. The effective field convention is the aggregate symmetric three-flavor normalization. Each microscopic field is divided by $\sqrt{3}$ and its microscopic charge-to-mass ratio is multiplied by $\sqrt{3}$. Charge encounters sample one eligible opposite-sign flavor with density proportional to $D/(6a^3)$. The Maxwell work uses the full conformal $J$ defined above. In the energy ledger for each component,

$$
\Delta\left[\tfrac{K}{2}(F^2+G^2)\right]
=\int 2a[i^2-(F-i)^2]\,d\eta-\int FJ\,d\eta.
$$

We use the radial solver's measured incident wave, wall radius and production history. Each transverse component has amplitude $0.04/\sqrt{2}$ times its unit early or delayed wave. The wall front already contains stored electromagnetic energy at first birth. This model tests the subsequent redistribution of that energy with a specified supply of daughters. The kernel is the inherited opposite-charge elastic vector pole with its energy-dependent cross section, invariant Møller flux and spin-averaged angular distribution. Changing the daughter supply, modeling different scattering channels, or letting the matter stress alter wall gravity changes the physical problem. The companion radial study propagates the additional reflected wave back through the collar. [Field history](../fields/README.md) · [Radial return](../return/README.md)

The integrator splits electromagnetic leakage and rotation from the simultaneous particle-current update. It integrates leakage at fixed midpoint input, rotates the conformal fields analytically under their curl and matches the midpoint electric impulse to the particles' relativistic energy increment. It accounts for partial steps of new births and for both momentum signs. Collision tests check energy and momentum in individual encounters. For all 83 primary and extended simulations the maximum componentwise field imbalance is below $2.0\times10^{-14}$ and the matter imbalance below $8.4\times10^{-15}$ in units of the larger of one, initial total stored field energy, incident energy and final gas energy. These residuals test bookkeeping. Finite particle samples, input histories, split timesteps and the choice of collision model have separate uncertainties. Control modes restore the prescribed axial fluid current, restore all prescribed currents, or turn off collisions. [Numerical source](solver.cpp)

## Reproduction

Run the study driver for the full package, or supply the two pinned parent study directories and the wall histories produced by the field extractor. `dedicated_half.npz` below must cover the relevant endpoint. The primary input extends through 6.6805. For the other two commands generate the extended input with `fields/extract.py --cases dedicated_half --end 7.0505` and the normal extractor arguments. The code writes the sample momenta and field traces to transient snapshot directories beside the outputs.

~~~sh
python -I -B live/evolve.py --polarization "$polarization" --scattering "$scattering" \
    --history "$primary/dedicated_half.npz" --output "$results/LIVE.json" --jobs 2
python -I -B live/follow.py --polarization "$polarization" --scattering "$scattering" \
    --history "$extended/dedicated_half.npz" --output "$results/CROSSINGS.json" --jobs 2
python -I -B live/cohorts.py --polarization "$polarization" --scattering "$scattering" \
    --history "$extended/dedicated_half.npz" --output "$results/COHORTS.json" --jobs 2
~~~
