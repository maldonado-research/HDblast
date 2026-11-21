# Can a collision redirect a charged pair into a different field?

The three charged species in the shared-resonance model respond to fields pointing along different axes. A collision can convert an oppositely charged pair of one species into an oppositely charged pair of another. At wall time 6.67 this process lowered the largest pressure of a single species from **97.6% to 87.7%** of its rest-frame energy density. Its particle counts also moved closer to the counts expected for a single thermal gas. The remaining difference was **38.4 percentage points** in the 20-cell test. An isotropic, ultrarelativistic gas would put about one third of its energy density into pressure along each principal direction. [Recorded results](FLAVOR.json)

## Measured histories

Both calculations use the same shared-pole birth history, incident field at the shell, mass spectrum, total collision cross section, scattering angles and charge-force law. The control keeps the outgoing species equal to the incoming one. The physical model chooses each outgoing species with probability one third. In both runs all three components of the current respond to the particles and feed back into Maxwell's equations on the shell. The radius, particle births and radial incident fields follow their supplied histories. Times are in units of the model scale $L$.

| Wall time | Species changes enabled | At least one encounter, fraction of particles | Largest single-species pressure / energy density | Count distance to common thermal gas | Particle conformal energy |
|---|---|---:|---:|---:|---:|
| 6.30944 | No | 0.0405 | 0.9988 | 0.4986 | 71.994 |
| 6.30944 | Yes | 0.0399 | 0.9954 | 0.4959 | 72.029 |
| 6.45000 | No | 0.1440 | 0.9971 | 0.4926 | 493.883 |
| 6.45000 | Yes | 0.1437 | 0.9798 | 0.4696 | 494.287 |
| 6.67000 | No | 0.5417 | 0.9762 | 0.4659 | 500.599 |
| 6.67000 | Yes | 0.5421 | 0.8774 | 0.3840 | 502.002 |

The finest sample has 16,384 positive and 16,384 negative numerical particles in the representative species. The other two species are its cyclic copies. Each time uses its own sample of all births through that marker. The common-gas comparison fits a classical Maxwell–Jüttner gas to the measured mean energy and momentum of the two signs together. Each of 20 bins counts particles by rest-frame kinetic energy and direction relative to that species' beam. Half the sum of absolute differences gives the count distance. A separate comparison with the two prescribed, oppositely drifting reference streams changes from 0.5693 to 0.5404 at the last marker. These two comparisons use different targets. [Shared Maxwell evolution](../live/) · [Thermal targets](https://github.com/maldonado-research/HDblast/blob/6cb09d2f01c49aa4ca975d77775a97352652e399/research/HDBLAST_SCATTERING_20261121/kinetic/scatter.py)

Summing the three species at the last marker gives principal pressures divided by energy density of $(0.3247,0.3247,0.3506)$ for the control and $(0.3243,0.3243,0.3514)$ with species changes. Similar aggregate pressure can therefore accompany appreciably different single-species momentum distributions. The field supplies the extra 1.40 units of gas energy through its changed work on the particles. The collisions conserve energy individually. [Energy accounts and pressure eigenvalues](FLAVOR.json)

## How much mixing survives to wall time 7.05?

Extending the shared source from the same initial data covers 90.43% of the archived absolute electrical work at wall time 7.05. In the finest run, permitting species changes moved the common-gas distance from **0.3924 to 0.2718** and lowered the largest single-species pressure fraction from **0.8613 to 0.5640**. Only the latter measurement concerns pressure. The former compares 20 particle counts with a gas at the same measured mean energy and momentum. Along the three axes of a single species the pressure fractions with conversion were 0.2123, 0.2237 and 0.5640. The three-species aggregate gave 0.3254, 0.3254 and 0.3492. [Longer integration](EXTENDED.json)

| Observable at 7.05 | Species fixed | Species changes enabled |
|---|---:|---:|
| Mean encounters per particle | 4.881 | 4.300 |
| Particles hit at least once | 89.01% | 85.02% |
| Count distance to the inherited two-stream target | 0.4245 | 0.3974 |
| Particle conformal energy | 500.563 | 501.758 |
| First axial electric-field sign change | 6.47392 | 6.47383 |
| Second axial electric-field sign change | 6.68604 | 6.67900 |
| First axial current sign change | 6.56849 | 6.56531 |
| Second axial current sign change | 6.76216 | 6.74735 |

The collision count declined even as the representative distribution spread across more directions. A change of species places the outgoing particles in a different accelerating field for their subsequent motion. The four converting runs span 4.280–4.375 encounters per particle, 0.2665–0.2792 in common-gas count distance and 0.5640–0.5718 in their largest single-species pressure fraction. The unrestricted three-population calculation at the same wall time finds 4.338 encounters per particle and largest single-species pressure fractions between 0.568 and 0.571 across its three species. These results refer to the specified pulse and evolving fields on the supplied expanding shell. Measuring a hot equilibrium radiation era would also require particle-production and radiative processes beyond this elastic mass-degenerate gas. [Independent longer integration](../independent/LATE.json)

## Collision law and coordinates

The equal-width vector resonance has one third of the dedicated-pole cross section when summed over all final species, or one ninth into a particular final species. Incoming partners must have the same species and opposite charge. At total accepted comoving daughter rest $D$ the physical number density for that partner is $nL^3=S_0D/[6(mL)a^3]$. The three species receive equal shares of the normalized field energy. The field shared by this cyclic reduction has current $QD\langle s\mathbf{v}\rangle$ in local coordinates, where $s$ is the charge sign and the brackets average all live particles of that species. Each accepted pair uses the inherited vector differential cross section in the two-body center-of-momentum frame. The stochastic clock tests the live momenta and birth times. [Inclusive rate, angular law and units](https://github.com/maldonado-research/HDblast/tree/6cb09d2f01c49aa4ca975d77775a97352652e399/research/HDBLAST_SCATTERING_20261121/rate) · [Source normalization](../SOURCES.md)

Let $R(x,y,z)=(z,x,y)$ and let species $f=0,1,2$ use the physical basis $R^f$. We sample $j=0,1,2$ with equal probabilities. Each pair of outgoing local momenta is sent to species $g=(f+j)\bmod 3$ using $\mathbf{p}'_g=R^{-j}\mathbf{p}'_f$. It follows that $R^g\mathbf{p}'_g=R^f\mathbf{p}'_f$. Each daughter's physical momentum just after scattering is continuous through the change of species. In subsequent steps it feels the field appropriate to the new species. The equivalent events in the other two rotated populations use the same random draw. There is always one pair entering and one pair leaving each affected species across the three copies. Each physical event transfers one positive and one negative particle together and preserves net charge separately for each gauge field. This reduction enforces the exact cyclic symmetry of the prescribed source, equal masses and equal couplings. Its numerical sampling suppresses fluctuations between copies. [Species transformation](conversion.hpp) · [Runner](shared.py)

At the last marker the representative sample experienced 15,077 collisions. Of these, 4,976 retained the species, 4,999 went to the next axis and 5,102 went to the remaining axis. Each particle participated in 0.920 collisions and in 0.617 species-changing collisions on average. These means count repeated encounters. The control with the same total cross section experienced 0.960 collisions per particle. The effective rate changes as the momentum distributions evolve. [Collision counts](FLAVOR.json)

The unrestricted three-population solver samples the three fields separately and allows the particle counts per species to fluctuate. With 4,096 birth pairs per species at wall time 6.67, it finds a hit fraction of 0.539 and mean 0.914 encounters per particle. Its largest combined pressure divided by energy density is 0.3519 and its gas energy is 501.95. The corresponding cyclic calculation with 16,384 pairs per species gives 0.542, 0.920, 0.3514 and 502.00. Agreement on these measurements provides a check on the symmetry reduction. [Independent three-population calculation](../independent/)

## Controls and reproduction

A separate constant-rate check puts one opposite-charge pair at the exact center-of-momentum resonance. Its constant scattering rate predicts 0.800 events per trial and at least one hit in 55.07% of trials. Across 8,192 seeds the code finds 6,464 events and 4,477 trials with a hit. The three relative destinations receive 2,114, 2,198 and 2,152 of the events ($\chi^2=1.64$, two degrees of freedom). A second control explicitly computes global Cartesian momenta in all three bases for nonzero incident total momentum and tests each species' charge balance. Across all 24 expanding-shell runs the residuals of individual collisions after species assignment are below $3.4\times10^{-16}$ in momentum and $8.3\times10^{-16}$ in energy, normalized by incident pair energy. Field and particle energy accounts close to $1.9\times10^{-14}$ and $5.0\times10^{-15}$ of the stated energy scale. These checks test the integration and relabeling. [Executable controls](control.cpp) · [Full results](FLAVOR.json)

At the last marker the four sampled resolutions use 4,096, 8,192 or 16,384 birth pairs and 2,048 or 4,096 time steps. The 8,192-pair time-step comparison holds birth scramble and collision seed fixed. Each different population size uses its own scramble and collision seed. Across these runs the common-gas distance ranges from 0.4629 to 0.4659 with species fixed and from 0.3820 to 0.3875 with species changes. These empirical ranges describe sensitivity across these four settings. Spatially varying densities, symmetry-breaking fields, independently evolving shell geometry and the production of equilibrium radiation remain outside this reduction.

From the study directory, set `polar` and `scatter` to the predecessor study directories and create the input with `fields/extract.py`. The output directory also receives the particle and field arrays for angular-moment analysis. [Shared history generation](../fields/) · [Current solver](../live/)

```sh
python -I -B flavor/shared.py --polarization "$polar" --scattering "$scatter" \
    --history /tmp/axial-inputs/shared_half.npz --output /tmp/shared-study/FLAVOR.json --jobs 2
```

For the separate later endpoint, regenerate the shared history to wall time 7.0505 and specify `--endpoint`. The four converting tests use three independent collision and birth samples, including one fixed-sample time refinement. The species-fixed control matches the finest sample and time grid.

```sh
python -I -B fields/extract.py --study "$polar" --cases shared_half --end 7.0505 \
    --output-dir /tmp/later-input --output /tmp/later-input/FIELDS.json
python -I -B flavor/shared.py --polarization "$polar" --scattering "$scatter" \
    --history /tmp/later-input/shared_half.npz --endpoint 7.05 \
    --output /tmp/shared-later/EXTENDED.json --jobs 2
```
