# How much can a pattern on the whole sphere change the particle energy?

We use a unit round three-sphere for comoving positions. Its physical radius is the recorded $a(t)$. Births, the radial background, the axial electric and magnetic fields, and the gauge coefficient follow the original history. Two calculations use the Maxwell perturbations driven from the outer emitting surface. One follows collisionless test particles through a specific smooth global field. The other bounds the energy for a whole space of possible angular patterns.

## Source patterns and particle motion

Take any real superposition of the two families of divergence-free vector harmonics with curl eigenvalue $\pm\kappa$, $2\le\kappa\le40$. Each family contains $\kappa^2-1$ independent modes. Normalize them by mean squared norm on the unit sphere. Every launched coefficient multiplies the same original **time** waveform. Write $\epsilon$ for the square root of the sum of the squared coefficients relative to the original launch amplitude. Its square measures the mean squared additional launch field relative to the reference. The change in total incident power also contains an interference term when the original polarization overlaps the added field.

The bulk calculation supplies signed coefficients $E_\kappa(t)$ and $B_\kappa(t)$ for a unit launch. The two helicities have the same electric transfer and opposite magnetic sign. Fix the unit Hopf vector along the reference electric field. It is a Killing vector, so a free particle conserves the component of its **comoving** momentum along it along every geodesic. Let $D(t)$ be the comoving number of daughters up to time $t$, in the common rest-mass normalization, $Q$ the absolute charge-to-mass ratio in the model, and $T$ the comparison time. The sign of the axial field $X$ stays fixed. Project the particle equation along this Hopf vector, sign each particle's momentum by its charge and the direction of the reference force, and average over births. The mean Lorentz factor obeys

$$
\langle\gamma\rangle_T\ge {Q\over a(T)D(T)}
 \int_{t_0}^T a(t)\,[D(t)|X(t)|-\epsilon L(t)]\,dt-p_*.
$$

Here $p_*$ bounds the magnitude at birth in units of $mc$. Every particle has $\gamma\ge |p_\parallel|/(mc)$, including pairs streaming in opposite directions. An upper bound $L$ on the contribution from the additional field gives a lower bound on the final energy.

**Arbitrary birth positions.** The sum of the dyadics $Y\otimes Y$ of any orthonormal basis in either helicity space is $(\kappa^2-1)g/3$ at every point. Rotational symmetry fixes the tensor and its integrated trace gives the coefficient. The projected additional Lorentz force tests the electric field along the axis and the magnetic field against a perpendicular vector of length at most one. Their cross term vanishes in the sum over a complete helicity space. Cauchy–Schwarz in the launch coefficients therefore gives

$$
L_{\rm point}=D\sqrt{\sum_\kappa {d_\kappa\over3}
      (E_\kappa^2+B_\kappa^2)},\qquad d_\kappa=2(\kappa^2-1).
$$

**Uniform birth positions with subsequent motion.** Every cohort starts uniformly on the sphere in the reference model. By time $t$ a particle born at $s$ can travel comoving distance at most $\int_s^t du/a(u)$. The sum of these maximal displacements over particles then present is $J(t)=\int_{t_0}^t D(u)\,du/a(u)$. Compare each field's magnitude at the later position with its magnitude at the birth position, evaluating the field at time $t$ at both points. The change is bounded by displacement times the maximum directional covariant derivative. For a fixed unit direction $w$ one has $\sum_Y|\nabla_wY|^2=d_\kappa(\kappa^2-2)/3$ over both families. The rough Laplacian eigenvalue is $\kappa^2-2$ and rotational symmetry shares the integrated gradient norm equally across the three directions. Define $G_E^2=\sum_\kappa d_\kappa(\kappa^2-2)E_\kappa^2/3$, with the same expression for $G_B$. Orthogonality and the uniform birth measure bound the initial-position mean field magnitudes by their spatial RMS. Thus another upper bound is

$$
L_{\rm travel}=D(\max_\kappa |E_\kappa|+\max_\kappa |B_\kappa|)
                   +J(G_E+G_B).
$$

Using the smaller of these bounds permits correlations between the fields, velocities and **later** positions. Integrals use endpoint rectangles chosen to weaken the exclusion for the recorded linear interpolation and the monotonically increasing $a,D$. The geometric identities are given in the conventions of the unit sphere. [Lindblom, Taylor and Zhang](https://arxiv.org/abs/1709.08020)

### The target frame

The earlier fluid has positive and negative species drifting in opposite directions with the same proper temperature. Its **wall-frame** mean energy includes that drift. The energy bound below tests an additional, stationary-gas interpretation. These energy thresholds do not exclude the original moving-fluid energies. The transverse momentum counts, in contrast, use the moving-fluid reference directly, as shown below.

| At 1% of whole-run absolute electrical work, finer time step | Separate interactions | Shared interaction |
|---|---:|---:|
| Lower bound on mean energy in units of stationary thermal mean, 1% launch RMS, arbitrary birth positions | 2.15 | 3.84 |
| Same energy bound, initially uniform positions and finite travel | 7.79 | 14.01 |
| Necessary launch RMS for energy to reach the stationary value with arbitrary positions | 1.19% | 1.26% |
| Same stationary-energy screen with uniform births and finite travel | 15.99% | 16.56% |

The corresponding uniform-birth floors on the original time grid are 15.97% and 16.54%. On a doubled radial grid the first is 15.95% using the same comparison knots. The largest allowed travel by this marker is 0.00261 or 0.00308 of the comoving radius. The energy columns compare with an isotropic gas **at rest on the wall** at the specified temperature. A larger perturbation can change the source, background, trajectories and births. Additional angular modes, different time waveforms, collisions and fields induced by currents require their own evolution. Reaching the quoted amplitude only removes this particular exclusion. [Calculation](all_patterns.py) · [All markers, grids and bounds](ALL_PATTERNS.json)

## A test that also applies to the moving thermal streams

An axial Lorentz boost leaves every transverse momentum unchanged. The equal-time number weight under a boost acquires a factor $\Gamma(1+\beta v_\parallel)$ in units with light speed one. For an isotropic rest gas the term odd in $v_\parallel$ integrates to zero over the cylinder $p_\perp\le p$, and also in the total number. Thus the normalized thermal transverse distribution $F(p)$ is identical in the rest and axial-stream frames at the same proper temperature. It remains a direct test of the original benchmark.

For each birth cohort the inherited axial fields can only rotate transverse comoving momentum. The triangle inequality bounds the population RMS of the **additional** transverse impulse by $\epsilon c_s$, where
$$
c_s={Q\over a(T)}\int_s^T a(t)\,[M(t)+\ell(s,t)G(t)]\,dt,\quad
M=\max|E_\kappa|+\max|B_\kappa|,\quad G=G_E+G_B .
$$
The spatial-RMS and geodesic-distance estimates apply to the whole cohort in its birth measure. A second bound uses the pointwise supremum of each field. Its integrand is at most $\sqrt{\sum d_\kappa E_\kappa^2/3}+\sqrt{\sum d_\kappa B_\kappa^2/3}$. Either gives a valid impulse bound. At cutoff $p>p_*$, any particle above $p$ must have picked up at least $p-p_*$ of transverse impulse. Markov's inequality for its squared magnitude therefore gives a lower bound on the predicted count below $p$:
$$
F_{\rm particles}(p)\ \ge\ {1\over D(T)}\int \left[1-\left({\epsilon c_s\over p-p_*}\right)^2\right]_+\,dD(s).
$$
For arbitrary birth positions the pointwise impulse bound permits counting every cohort whose entire support remains below the cutoff. The code gives each birth the whole time cell and uses upper endpoint rectangles for the impulse. It checks 420 trial cutoffs and records an explicit witness where the lower count exceeds the moving-gas prediction. Additional cutoffs can strengthen an exclusion.

With uniform births and **1% launch RMS**, the separate-interaction case gives at least 18.30% of particles below $666.5\,mc$ versus 10.21% in the reference gas. The shared case gives at least 4.79% below $207.4\,mc$ versus 2.51%. These are minimum discrepancies of **8.09 and 2.28 percentage points** on the time-refined histories. To approach the reference counts within one percentage point for every cutoff, the necessary launch RMS is approximately **2.70% or 1.49%**, allowing all 44,200 coefficients and both curl signs. Original time steps give 2.66% or 1.47%. Doubling the radial resolution gives 2.66% on the separate-interaction history. These tests use freely moving collisionless particles, identical time profiles at the outer source and fixed original birth, geometry and scalar histories. [Boost-invariant count witnesses](ALL_PATTERNS.json)

## Following particles through one smooth spherical wave

Write a point of the sphere as two complex coordinates $(z_1,z_2)$ with $|z_1|^2+|z_2|^2=1$. For $\kappa=2,5,7,10$ we choose the one-form

$$
Y_\kappa=\sqrt{\kappa-1}\,{\rm Re}\,[z_1^{\kappa-2}(z_1\,dz_2-z_2\,dz_1)].
$$

Its norm has spatial RMS one and maximum $\sqrt{\kappa-1}$ because $|z_1|^2$ is uniform on $[0,1]$ under the volume measure. It is everywhere perpendicular to the chosen Hopf axis. The oriented orthonormal particle frame has $e_1=(-\bar z_2,\bar z_1)$, $e_2=i e_1$, $e_3=-i(z_1,z_2)$. With this orientation $\mathrm{curl}\,Y_\kappa=-\kappa Y_\kappa$ and $\mathrm{curl}\,e_3=-2e_3$. The brackets of the frame vectors are twice the next vector in cyclic order. Writing $f=z_1^{\kappa-2}$ gives $e_2f=i e_1f$ and $e_3f=-i(\kappa-2)f$, which yield the curl identity and zero divergence. The code also checks these identities by differentiation along great circles. Its magnetic sign then satisfies Faraday's law with the archived convention $d(a^2 B_\kappa)/dt=+\kappa a E_\kappa$. The additional electric and magnetic fields both point along $Y_\kappa$ at each position.

The calculation integrates both signs of each emitted pair. Positions drift along great circles and the comoving momentum receives symmetric electric kicks and a magnetic rotation. The sampled birth time, angular acceptance and mean axial fields come from the original history. Scrambled Sobol points sample the initial positions, birth times and angles. Changes of scramble, time step, number of samples and native background resolution are reported for selected cases. In the zero-launch control the guide field preserves the transverse comoving momentum magnitude. The maximum sphere-normalization error is recorded for every run.

For example, with mode 5 and launch RMS fractions 4% and 2%, the finer spherical runs give 18.4- and 20.2-percentage-point discrepancies in the transverse cumulative counts. Mean energies are only 2.83% and 0.81% above the corresponding moving-fluid predictions. We resolve pressure at each position along the local harmonic and across it. Along the field pattern the two runs give 3.91 and 3.28 times the predicted lateral pressure. Across it they give only 0.000363 and 0.000354 of the predicted lateral pressure. A particle stream in local equilibrium has equal pressure in these two directions even after its longitudinal boost. Measuring directions relative to the local wave avoids hiding anisotropy when the wave orientation varies around the sphere. These are properties of the tested paired stream. A sum over several rotated streams can have different stress symmetry. The particle-current field is evolved separately in the local calculation. [Spherical trajectories](round.py) · [Frame-matched spherical results](SPHERE.json) · [Local particle-current calculation](../orbits/README.md)
