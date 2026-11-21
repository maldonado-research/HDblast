# Three independently evolving charged flavors

Does the threefold cyclic reduction impose the appearance of isotropic pressure or conceal uneven changes in the populations? We tested the shared-resonance benchmark in one global orthonormal frame. The calculation carries separate lists of positive and negative particles for each of the three gauge flavors. Each flavor sources and responds to its own three-component Maxwell field. The three birth distributions use independent scrambled sequences. The collision scheduler can move individual pairs between flavors. This allows the flavor populations, their currents and their fields to differ during the run.

## Results

The table uses 4,096 birth pairs per flavor, or 24,576 particles. The listed collision count includes multiple hits on the same particle. The pressure entries are eigenvalues of the stress tensor in the appropriate energy rest frame, divided by the energy density there. TV20 uses the 20 bins and inherited sign- and flavor-dependent streaming-fluid target of the parent study.

| Wall time | Mean collisions per particle | Particles ever hit | TV20 | Largest pressure in each flavor | Three-flavor pressures |
|---|---:|---:|---:|---|---|
| 6.30944 | 0.046 | 4.50% | 0.594 | 0.996, 0.995, 0.995 | 0.323, 0.325, 0.352 |
| 6.45 | 0.156 | 14.36% | 0.543 | 0.980, 0.980, 0.980 | 0.324, 0.324, 0.352 |
| 6.67 | 0.914 | 53.87% | 0.544 | 0.884, 0.883, 0.884 | 0.324, 0.325, 0.352 |
| 7.05 | 4.338 | 85.05% | 0.400 | 0.568, 0.570, 0.571 | 0.325, 0.327, 0.348 |

At 6.67 the positive counts in the three flavors were 4,095, 4,058 and 4,135. Each flavor had equally many negative particles. The cyclic calculation with 16,384 birth pairs in the representative flavor gave 0.920 collisions per particle, a hit fraction of 54.21%, TV20 of 0.540 and a largest total pressure of 0.351. The unrestricted solver gave charged conformal energy of 501.95 against 502.00 in the cyclic calculation. In two additional unrestricted samples of 2,048 birth pairs per flavor the collision count ranged from 0.883 to 0.918. The variation gives a check on the numerical importance of tying the three orientations together at this benchmark. [Early results](INDEPENDENT.json) · [Cyclic calculation](../flavor/)

At 7.05, about 90.43% of the shared reference fluid's absolute electrical-work integral has elapsed. The separate flavor pressures are closer together. With the same inclusive cross section and conversions suppressed, the largest pressure is still 0.862–0.869 in the 2,048-pair-per-flavor control. Switching off collisions gives a largest pressure close to unity. All three independent conversion runs give a largest combined pressure between 0.347 and 0.350. [Late results](LATE.json)

## What survives after mixing the three flavors?

A gas can have equal pressure along three axes while its fast particles still prefer particular directions. We measured the angular distribution directly in the global frame, pooling all the particles and both charge signs at their actual statistical weights. Choose six separate spherical caps around the positive and negative coordinate axes at a fixed half-angle $\alpha$. At least two of the three coordinate axes make an angle of at least $45°$ with either pole of any trial symmetry axis. For a distribution uniform in azimuth around that trial axis, the probability on a single latitude inside each of those four caps is bounded by $\arcsin(\sqrt2\sin\alpha)/\pi$. This follows by maximizing the angular width of the intersection of a latitude with a cap over all possible latitudes. Summing over the four caps gives a ceiling valid for any distribution over latitude. The observed probability has a lower bound given by the least-massive pair of opposite-cap pairs. Their difference lower-bounds angular total variation against every axisymmetric distribution. Boosted isotropic thermal particle distributions are among these targets.

At 7.05, each union of two opposite pairs of 2° caps contained at least 20.58% of the particle energy. The ceiling for the union selected by a trial symmetry axis is 6.29%. This gives a gap of at least 14.29 percentage points. The other two runs give 14.68 and 14.84 percentage points. The same cap statistic weighted by particle count becomes inconclusive at this late marker. A broad angular moment can still distinguish the patterns when the narrow peaks spread out. We computed the degree-four spherical-harmonic moment of the full particle list and bounded its distance to the family of all axisymmetric moments. The test uses the actual global directions and occupancies. It gives a gap of at least **7.24 percentage points by number and 14.48 percentage points by energy** for the largest run. The corresponding ranges across three independent collision sequences are 7.17–7.47 and 14.48–14.90 percentage points. These bounds apply to the finite modeled samples. They measure the residual angular structure after most particles have experienced a collision. [Angular bound and derivation](../triad/) · [Measured caps, moments and energy ledgers](LATE.json)

## Normalization and checks

Let $N$ denote the final number of birth pairs per flavor at the requested marker, $D$ the marker's total accepted comoving rest, $S_0$ the action normalization, and $\mu=mL$ the daughter mass. Each simulated particle represents a share $D/(6N)$ of the total comoving matter. The normalized per-flavor Maxwell equation uses signed current $3Q\sum_{i\in f}(D/6N)s_i\mathbf v_i$. The initial fields and incoming characteristics are cyclic rotations of the canonical triad. They subsequently evolve independently. We average the three gauge-field energies for comparison with the total matter. The discrete field work enters the matter ledger with the opposite sign.

The inclusive shared-pole candidate rate summed over flavors is

$$
\lambda = {D\over N}{S_0\over18\mu^3 a^3}K_{\max}
                 \sum_f n_{f,+}n_{f,-}.
$$

Here $n_{f,s}$ is the instantaneous candidate count for a sign and flavor, and $K_{\max}$ is the dedicated kernel's rejection envelope. The $1/18$ includes the shared inclusive width factor of $1/3$ and the designated opposite-charge, one-flavor target normalization $1/6$. A candidate flavor is chosen in proportion to its pair count, then a particle of each sign is chosen uniformly within that flavor. A birth-time check discards candidates sampled before either particle was born. After each accepted collision the two four-momenta are scattered in the global frame with the verified parent differential kernel. Both particles enter one uniformly drawn outgoing flavor. The event rate immediately reflects the changed populations. Fields and currents, charge per flavor, global particle number, and momentum and energy per collision are tracked. During the collisionless cyclic control, identically sampled initial births on separate lists yield three matching rotated fields to better than $10^{-12}$.

The maximum scaled field and matter energy residual through the late marker is $1.3\times10^{-12}$. The largest relative pair-energy error is below $8\times10^{-16}$. The calculation uses three interacting flavor populations on the given homogeneous background with the inherited source, birth law, radial incident characteristic and cross section. It tests the cyclic statistical closure for this trajectory and supplies a finite-time diagnostic for equilibration. The origin of the parent source, reaction network in an inhomogeneous wall and simultaneous gravitational response remain open. [Ensemble solver](ensemble.cpp) · [Driver and checks](run.py)

~~~sh
python -I -B independent/run.py --polarization /path/to/HDBLAST_POLARIZATION_20261121 \
    --scattering /path/to/HDBLAST_SCATTERING_20261121 --history /tmp/axial-inputs/shared_half.npz \
    --output /tmp/INDEPENDENT.json --jobs 2
python -I -B fields/extract.py --study /path/to/HDBLAST_POLARIZATION_20261121 \
    --end 7.0505 --cases shared_half --output-dir /tmp/axial-late --output /tmp/SHARED_FIELDS.json
python -I -B independent/run.py --polarization /path/to/HDBLAST_POLARIZATION_20261121 \
    --scattering /path/to/HDBLAST_SCATTERING_20261121 --history /tmp/axial-late/shared_half.npz \
    --output /tmp/LATE.json --extended --jobs 2
~~~
