# Timed pulses through the finite five-dimensional collar

A second pulse can arrive while charged particles are forming. At equal emitted energy, a short pulse arriving at wall time 6.260 gives a peak electric coefficient of 0.00209 for angular wavenumber 7. The original pulse leaves a coefficient of 0.997 at first birth. A pulse arriving earlier can deposit a much stronger field before birth. Its angular concentration and time dependence then matter for how many particles it can accelerate sideways.

## Launch and normalization

The radial Maxwell and wall equations are those of HDblast PR #197. Each linear channel propagates on the archived spherical mean. The field histories use the wall proper time and orthonormal physical electric and magnetic components. The coefficients multiply normalized curl eigenforms on the unit three-sphere. The primary experiments keep five angular wavenumbers, 2, 5, 7, 40 and 258, and seven time profiles. The unforced probe stays exactly zero and the original wavenumber-2 probe reproduces the reference electric field up to the onset of daughter production. Each channel checks its wall energy against the sum of net boundary work and work done by the evolving gauge coefficient.

The incoming characteristic is prescribed at the moving outer tube using its proper time τ. The reference emitted energy in the solver's angular normalization is

$$
U_0=2\int r_{\rm tube}^3 i_0^2\,d\tau=0.97803744.
$$

Each profile is rescaled to this same energy. The measured dimensionless source matrix for superpositions is

$$
W_{ij}=\frac{2}{U_0}\int r_{\rm tube}^3 i_i i_j\,d\tau .
$$

For real coefficient vectors in orthogonal angular channels, the total additional emitted energy in units of the original is the sum of their quadratic forms in W. The square root of that total is denoted by ε. The tube area changes during expansion, so the emitted energy matters when comparing launches at different times. Wave energy reaching the wall and work transferred by the wall's changing gauge coefficient are separate ledger entries. Cross terms within a single angular channel must be included when combining overlapping profiles.

| Source pattern | Outer tube prescription | Measured behavior |
| --- | --- | --- |
| Early | Original compact pulse | Stored wavenumber-7 field survives until first birth |
| Prepare | Original pulse delayed by 4.989 proper-time units | Stores a larger wavenumber-40 field before release |
| Release | Compact pulse centered at 6.2100 with width 0.008 | Drives the wall before the first birth |
| Birth | Center 6.2446, width 0.006 | Incoming wavenumber-7 wall peak at 6.260002 |
| Late birth | Center 6.2546, width 0.006 | Incoming wavenumber-7 wall peak at 6.270016 |
| Ramp | Center 6.2446, width 0.026 | Tests a broader wave during births |
| Bipolar | Center 6.2486, width 0.022 | Reverses sign halfway through its compact support |

The dedicated half-step history begins sampling births at 6.251625. The late positive pulses have a sine-squared envelope. The bipolar pulse multiplies that envelope by a full sine period and has zero unweighted signed area. The seven-profile source matrix has eigenvalues from 0.211 to 2.135, so the independent drivers remain distinguishable at launch.

## Can a pulse supply two independent transverse impulses?

For a cohort born at wall time s and assessed at the 1% work marker T, define the fixed-position direct electric impulse
$$
J_i(s)=\frac{Q}{a(T)}\int_s^T a(t)E_i(t)\,dt .
$$

Use Q = 10⁷ in the archived dimensionless units. The population measure is normalized from the reference daughter production history D(s). The axial reference momentum gives a Lorentz factor γ₀(s). Weight the cohorts by their population divided by γ₀ to form the leading transverse-stress diagnostic. With a Cholesky factorization W = LLᵀ, the first two singular values of the weighted matrix of impulses multiplied on the right by L⁻ᵀ are denoted by s₁ and s₂. At unit total launch energy, the greatest attainable weaker eigenvalue of this two-axis stress matrix is

$$
\lambda_*=(s_1^{-2}+s_2^{-2})^{-1}.
$$

For a pair of orthogonal output vectors with equal squared length λ, the least required input energy is λ times the sum of their inverse squared gains. The two leading singular vectors minimize that sum. Conversely, if both output eigenvalues reach λ, projecting and scaling the outputs to equal lengths can only decrease the required input energy. The solver constructs the optimizing drives and verifies their total emitted energy and the complete two-by-two output matrix.

At wavenumber 7 the weaker-to-stronger singular amplitude ratio is 0.000344 for the dedicated history and 0.000303 for the shared history. Matching the two target transverse stresses per lab-frame particle in this calculation would formally require ε = 28.1 and 17.6 respectively. The calculation holds a particle's angular position fixed, omits the birth kick and magnetic deflection, and weights stress using the unperturbed axial Lorentz factor. Those amplitudes characterize the response matrix. The linearized estimates do not establish physical thresholds at large fields. The companion orbit calculation includes angular motion and magnetic forces.

## Where do the high-wavenumber fields act?

A second scan applied eight disjoint preparation delays at 12 wavenumbers from 7 to 515. The strongest electric coefficient at first birth in the sampled grid is 12.446, found at wavenumber 72 for a delay of 5.085. Doubling the radial mesh gives 12.514. A large coefficient measures a volume-normalized mode and can coexist with a very weak field over most of the wall.

For the ten sampled wavenumbers at or above 40, the following bound covers all mixtures of eight pulses in both real and imaginary quadratures, or 160 launch coefficients. It uses the highest-weight horizontal harmonics centered on the same Hopf circle. At each birth the radius coordinate u = |z₁|² is uniform from zero to one for a uniform population on the sphere. The mode magnitude is √(κ−1) |z₁|^(κ−2). Until T, travel on the unit sphere is bounded by ℓ = ∫dt/a from first birth, less than 0.00308 in these histories. The distance to the circle is arccos √u, so |z₁| at any later point is bounded by

$$
r(u)=\cos\bigl(\max[0,\arccos\sqrt{u}-\ell]\bigr).
$$

Write the time-dependent field coefficients at wavenumber κ as a vector Eκ across launch profiles. Maximizing the projection onto a unit direction by Cauchy–Schwarz gives the local bound
$$
 |E(t,u)|\le \epsilon\left[\sum_\kappa (\kappa-1)r(u)^{2\kappa-4}
 E_\kappa^\mathsf{T} W_\kappa^{-1} E_\kappa\right]^{1/2}.
$$
The real and imaginary quadratures are perpendicular with equal lengths at a point. A projection onto a single unit direction therefore counts their squared lengths only once. The same estimate applies to B. The axial background preserves the transverse comoving momentum magnitude in free propagation. Bounding the added electric and magnetic forces by their magnitudes, writing the unit-amplitude bounds as e and b, gives
$$
 |p_\perp(T)|/(mc) \le p_*+\frac{\epsilon Q}{a(T)}\int_s^T a(t)\,[e(t,u)+b(t,u)]\,dt .
$$
The birth kick p★ bounds the initial transverse momentum in the same units. This force bound allows arbitrary particle directions, finite angular displacement, both charge signs and nonlinear test-particle deflection. The fields and the birth history are prescribed.

The numerical implementation uses the larger endpoint dual norm on each cell of the piecewise-linear field history, the largest radius in each volume bin, the largest expansion factor on each time cell, and the earliest possible birth in that cell. Each choice increases the permitted transverse momentum. Summing the resulting guaranteed low-momentum fraction over birth and volume cells gives a lower bound on the discrepancy with the moving thermal target. The thermal transverse cumulative fraction at p in units of mc is
$$
F_{\rm th}(p)=1-\frac{(1+p^2)K_2(\sqrt{1+p^2}/\theta)}{K_2(1/\theta)} .
$$
An axial boost preserves transverse momenta. Its odd-in-axial-momentum contribution to the transformed number measure integrates to zero.

| Additional launch amplitude ε | Dedicated minimum CDF discrepancy | Shared minimum CDF discrepancy |
| --- | ---: | ---: |
| 0.02 | 69.11 percentage points | 66.84 percentage points |
| 0.1 | 63.32 percentage points | 61.23 percentage points |
| 0.5 | 58.00 percentage points | 56.07 percentage points |
| 1.0 | 55.85 percentage points | 53.99 percentage points |

For example, with ε = 1 in the dedicated history, at least 57.40% of the particles have transverse momentum below 242.64 mc. The moving thermal target puts 1.55% below that cutoff. The bound survives a doubled radial mesh with a gap of 55.84 percentage points. The archived product also includes the 128-volume-bin comparison to the main 256-bin calculation. These bounds apply to the specified numerical histories and harmonic family. Other members of the angular eigenspaces can concentrate on different circles. The same calculation also covers groups of modes rotated onto distinct circles. Assign each wavenumber and its real and imaginary components to one circle. If a mixture uses at most M circles, then the probability of any birth-radius coordinate exceeding u is at most M(1−u), for any relative placements. The q-quantile of the largest radius squared is bounded by 1−(1−q)/M. Substituting that upper quantile for u in the calculation preserves the force and CDF bounds. An S³ isometry preserves the curl eigenvalue and angular norm, and distinct wavenumbers remain orthogonal in emitted energy.

With two arbitrarily placed circles and ε = 1, the CDF gaps are at least 14.11 and 10.33 percentage points. With ε = 0.02 they are 41.18 and 36.55 points. With three circles and ε = 0.02 they are still 14.39 and 7.38 points. This union bound gives no positive guaranteed gap at ε = 0.02 when all ten wavenumbers can have distinct centers. A self-consistent angular plasma, collisions, a changing geometry and a different birth history can also change the conclusion.

## Reproduction

The top-level reproduction runner builds the required parent and compact histories. With those inputs already prepared, the individual commands are

```sh
python fields/launch.py --parent parent --history history --output-dir output/primary --jobs 2
python fields/launch.py --parent parent --history history --output-dir output/scan --scan --jobs 2
OPENBLAS_NUM_THREADS=1 python fields/impulse.py --fields output/primary --scan output/scan --history history --output-dir output/impulse
```

The launch commands write compressed field histories and their source and wall ledgers. The impulse command writes the response metrics, localized distribution bounds and optimal two-axis waveform histories. FIELDS.json, STORAGE_SCAN.json and IMPULSE.json hold the compact reference results.
