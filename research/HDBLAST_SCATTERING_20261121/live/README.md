# What does the registered collision actually do while the particles form?

This calculation injects new charged pairs, accelerates them in two sideways electromagnetic modes and the axial field, and scatters eligible opposite charges as the wall expands. The two induced sideways fields respond to the evolving particle currents. The simulation starts with no daughters on the separate-interaction history from the preceding study. It uses that history's radius, axial field, birth law, external illumination and mean density. It follows the two charge signs for one representative beam direction. Other beam directions may compensate one another in a whole-shell pressure average. [All 66 runs and three absolute-rate controls](LIVE.json)

## Following the smaller pulse in time

The two added fields together use 0.16% of the reference pulse's emitted energy in addition to the reference axial source. Each marker denotes the accumulated fraction of the original history's whole-run absolute electrical work. The table averages three collision seeds, each with 16,384 birth pairs. We regenerate the population to each marker at its cumulative daughter mass. The 20-bin distance compares each charge with its moving thermal target. The pressure ratios in the last column use the combined stress tensor of this one beam's two signs, measured in its own energy rest frame.

| Reference work reached | Wall time | Fraction that has collided at least once | Encounters per particle | 20-bin distance | Three principal pressures / combined energy density |
|---|---:|---:|---:|---:|---|
| 0.1% | 6.26217 | 2.29% | 0.0231 | 59.48% | 0.0029, 0.0031, 0.9940 |
| 1% | 6.27855 | 7.53% | 0.0788 | 57.53% | 0.0025, 0.0026, 0.9950 |
| 5% | 6.30944 | 12.62% | 0.1360 | 55.87% | 0.0019, 0.0021, 0.9960 |

At the early 1% marker the collisionless distance was 61.91%. With collisions, the two sideways pressures of individual streams rose from 0.00000068 and 0.121 of their moving-gas target to **0.382 and 0.487**. Substantial lateral pressure can come from relatively few encounters. The combined gas remains dominated by motion along the beam. A single stationary ultrarelativistic thermal population would give three ratios near one third in the last column. Different species or beam directions, later epochs, and the emergence of a radiation-dominated cosmology require separate tests.

The high-cost control keeps the waves previously adjusted to fit the two sideways pressures at 1% work. At that marker, switching on collisions moved the pressure ratios from 0.982 and 1.011 to 1.324 and 1.428. The 20-bin distance fell from 47.16% to 44.68%. At 5% work its collisionless and colliding distances were 77.30% and 73.85%. That particular fit was specific to the earlier marker. The physical cost of this control is about 1,356 times the original emitted energy. [Predecessor pressure fit](https://github.com/maldonado-research/HDblast/pull/198)

## Collision law and energy accounting

The binary kernel $k=m^2\sigma v_{\mathrm M}$ is the registered pure vector resonance with pole mass $8300m$ and width parameter 0.30. Binary encounter rates use the invariant flux divided by the two wall-frame energies. At equal daughter masses the spin-averaged pure s-channel current contraction in the pair center-of-momentum frame gives

$$
\frac{d\sigma}{d\mu}\ \propto\ A+B\mu^2,\qquad A=3-2\beta_*^2,\quad B=\beta_*^4,
\qquad \beta_*^2=1-4m^2/s .
$$

The spatial current tensors are proportional to $\delta_{ij}-\beta_*^2 n_i n_j$ for the incoming and outgoing axes. Their contraction gives this angular law. A draw selects the new pair direction with this density and uniform azimuth. A boost back to the wall frame conserves each pair's total momentum and energy. The expected angular second moment is $(A/3+B/5)/(A+B/3)$, varying from one third at threshold to two fifths at high energy. An isotropic angular control changes only this differential law. The model specifies the resonant s-channel. A full same-flavor quantum-field calculation could also include crossed diagrams and interference. [Pair kernel and density](../rate/README.md) · [Collision operator](collide.hpp) · [Sources](../SOURCES.md)

For $N_+$ equal-weight positive particles the candidate event rate is $N_+ n_- k_{\max}/m^2$ in physical units with $c=1$. In the wall's dimensionless time the per-particle factor is $S_0 D k_{\max}/[6(mL)^3a^3]$. The six populations comprise three equally populated flavors and two signs. The simulation samples one flavor and applies its target density. A Poisson clock proposes opposite-charge partners uniformly, rejects pairs that have not yet formed, and accepts with probability $k/k_{\max}$. Repeated collisions are allowed. Within a numerical step the electromagnetic and expansion update precedes collisions, which act at the end-step radius. Halving the step checks this splitting. The calculation assumes mixing and factorized incident-pair statistics within this homogeneous cell. It does not resolve a correlated pair's separation or collisions in a spatially clumped gas.

The source, particle, induced-field and outgoing-radiation work are checked using the predecessor's energy ledger. Every collision is elastic. The changed distribution can draw a different amount of work from the imposed axial field. The axial field and production history here remain prescribed. The calculation includes transverse current feedback only. Tests at the 0.1%, 1% and 5% markers use the same source waveforms. The gauge kinetic coefficient after production is 2.54. The stored energies and field amplitudes in the original wave calibration recover this value at all three markers.

## Checks

- Switching collisions off recovers both published 4,096-pair early histograms. The absolute thermal integrals use the equivalent change of variable $x=u^2$ to control quadrature near zero kinetic energy at the other markers.
- Across the full simulations the relative error in total energy per collision is below $10^{-12}$. The matter and induced-field work residuals are below $3\times10^{-13}$ relative to their reference scales.
- On the lower-energy pulse the early three-seed mean changes from 57.57% to 57.53% in the count comparison when increasing from 8,192 to 16,384 pairs at 512 steps. Halving the step at 8,192 pairs gives 57.53%. Across the three finest seeds, 7.41–7.71% of the population collides. The three seeds measure numerical variation. They do not quantify uncertainty in the physical cross section.
- An independent two-particle test places opposite beams exactly on resonance at constant volume. Successive encounters rotate the pair without changing its collision rate. For a predicted mean of 0.8 collisions per run it measures 0.7975 in 8,192 trials. Beginning halfway through gives 0.4102 against 0.4. Two particles moving with identical velocities never scatter. [Control](control.cpp)
- The pressure tensor sums both charge signs before finding its principal pressures. The three pressures and rest-frame energy density are the spacelike and timelike eigenvalues of the wall-frame energy-momentum tensor with one index lowered. The code checks the massive-particle trace identity. The earlier two-pressure comparison treats each charge's drifting population separately.

~~~sh
python -I -B live/evolve.py --study /path/to/HDBLAST_POLARIZATION_20261121 \
    --fields /path/to/regenerated/fields/traces --output /tmp/LIVE.json --jobs 2
~~~
