# Does retaining the magnetic orbit and the three field directions change the heating test?

Each charged particle now carries both transverse momentum components as well as the longitudinal component. The three species accelerate along three different axes. Scattering with a shared mediator can change the outgoing species, so a momentum transverse to one field can become longitudinal to another. The particle representation keeps the physical triad throughout the replay. The earlier isotropic-azimuth computation and this extension use the same prescribed geometry, field boundary signal and accepted birth history. Both allow the electromagnetic field on the wall to respond to the measured particle current. [Driven solver](run.py) · [Three-dimensional replay](oriented.py)

## Heating and the surviving particle population

At the reference radius-doubling time, random particle energy accounts for **97.52%** of the added wall energy in the dedicated-channel replay and **97.73%** in the shared-channel replay. At this high resolution, **4.57% of the shared-channel representatives have never scattered**. The two runs give 4.30% and 4.83%. All dedicated-channel representatives have scattered by that time. A large random energy fraction can coexist with a smaller population retaining the narrow range of transverse momenta at birth. [Results and individual-seed ranges](ORIENTED.json)

| Interaction and marker | Instantaneous random energy fraction | Mean over preceding $0.25L$ | Largest energy-distribution gap | Particles with zero encounters |
|---|---:|---:|---:|---:|
| Dedicated, reference first hot | 90.97% | 90.63% | 9.67 points | 9.90% |
| Dedicated, reference radius doubling | 97.52% | 97.51% | 0.78 points | 0.00% |
| Dedicated, turnaround | 97.83% | 97.83% | 1.03 points | 0.00% |
| Shared, reference first hot | 89.39% | 90.62% | 27.84 points | 34.88% |
| Shared, reference radius doubling | 97.73% | 97.70% | 4.45 points | 4.57% |
| Shared, turnaround | 97.85% | 97.84% | 4.40 points | 4.36% |

The table uses two seeds with 4096 final equal-weight representatives, the halved source integration step and the finer particle integration. The quoted heat fraction uses the rest frame of the energy flow (Landau frame). The energy-distribution gap compares the energy CDF in the particle-number rest frame to the Maxwell–Jüttner distribution with the same mean rest energy. The sampled pressure and temperatures use the same definitions as the gyrotropic analysis. At shared-channel doubling the measured rest-frame longitudinal share of pressure is 33.57%, close to the isotropic value of one third, despite the residual energy-distribution gap. The fraction within $p_\perp^2/m^2\le0.1025$ agrees with the fraction that has never scattered to the archived precision. A small transverse momentum alone does not specify the particle's longitudinal energy. [Spectral estimators](spectra.hpp)

The clocks come from the earlier fluid solutions. For the dedicated channel the reference first-hot, doubling and turnaround times are 8.44044, 25.74105 and 38.44394 in units of $L$. For the shared channel they are 8.75689, 29.66306 and 38.24116. The early half-work marker uses the integrated **absolute electrical work of the reference fluid up to its first-hot time**. The fraction of particles having no encounters at that marker is 32.80% in the dedicated run and 53.26% in the shared run. Changes to the current and electromagnetic field can change the actual time distribution of work in the replay. The old first-hot time and doubling are comparison times. The instantaneous shared replay falls below 90% at the former marker in both high-resolution seeds; its preceding finite-window mean exceeds 90% in both. A continuous kinetic hot-growth interval would require threshold tracking along a self-consistently evolved geometry.

Holding the source history and 2048 representatives fixed, halving the particle step changes the shared zero-encounter fraction at reference doubling from **0.20% to 4.47%**. With 4096 representatives and the half-step source it is **4.57%**. This numerical sensitivity is much larger than the change in the total heat fraction. The gyrotropic calculation at 4096 representatives gives **4.69%** in the narrow tube and **97.69%** instantaneous random energy at doubling. The shared reference first-hot **instantaneous** heat fraction is more variable between seeds, resolutions and orbital treatments than the radius-doubling fraction. Two seeds demonstrate the spread seen here. They give limited precision for estimates of its tails. [Gyrotropic comparison](DRIVEN.json)

## Discrete axes and charge symmetry

The calculation uses one representative distribution for a positive charge with its field along $z$. Cyclic coordinate permutations generate the distributions for the two other flavors. The three field magnitudes, masses and couplings are equal. A particle that changes to the $x$-flavor is represented by $(p_y,p_z,p_x)$ in its new coordinate system. For a change to the $y$-flavor the assignment is $(p_z,p_x,p_y)$. The shared resonance gives each of the three final flavors probability $1/3$ and uses the inclusive incident rate. Dedicated resonances keep the incident flavor.

The physical negative charge is represented by the proper rotation $R=\operatorname{diag}(1,-1,-1)$ of a positive representative. For the force $q(\mathbf E+\mathbf v\times\mathbf B)$, with both fields along $z$, $R\hat z=-\hat z$ and $R(\mathbf v\times\hat z)=-(R\mathbf v)\times\hat z$. The rotation therefore transports the positive trajectory into the opposite charge's trajectory even when the magnetic field is active. Incoming particles sample different representative indices. Each pair updates both representatives with the two outgoing momenta and contributes to the charge-conjugate ensemble with half the full proposal count. This uses molecular chaos and the model's equality of the charge-conjugate distributions.

Independent rotations by integer multiples of $90^\circ$ about the incident axis sample the exact discrete symmetry of the equal three-axis preparation. Neighboring flavor axes may be permuted or sign-reversed by a rotation. Charge conjugation of the affected flavor and the parity-even incoming/outgoing vector cross section preserve the collision statistics. The four directions give zero transverse current and a diagonal pressure tensor with equal transverse diagonal entries when averaged over copies. They preserve the real and imaginary parts of $\langle (\hat n_x+i\hat n_y)^4\rangle$. These two measured moments can detect angular structure which a continuous azimuthal average removes. The separate first-product calculation gives a nonzero real fourth moment of 0.14431 for three higher-energy, split input beams. [First-collision angular test](../box/AZIMUTH.json)

## Magnetic rotation and collisions

On each free-flight interval let $q_\perp=a p_\perp/m$ be the comoving transverse magnitude, $q_z$ the longitudinal canonical momentum, $I$ the comoving electric impulse, and $E=\sqrt{1+(q_\perp^2+(q_z+I)^2)/a^2}$. The positive particle's transverse azimuth evolves as $d\phi/dt=-Q B/E$ with the solver's charge-to-mass ratio $Q=10^7$. The field is $B=B_{\rm ref}+2Z/a^2$, including the response of the wall's magnetic potential. The magnetic force does no work. Newly created particles sample a uniform azimuth independently of their birth polar angle. Their transport to the operator midpoint preserves that uniform birth azimuth.

During each electrical half step the integrator fixes $a,B$ at the midpoint and linearly interpolates the parallel impulse between the step endpoints. With $z_{0,1}=(q_z+I_{0,1})/a$ and $b=\sqrt{1+(q_\perp/a)^2}$ it integrates the rotation using

$$
\overline{E^{-1}}=
\frac{\operatorname{asinh}(z_1/b)-\operatorname{asinh}(z_0/b)}{z_1-z_0}.
$$

The coincident-endpoint limit is $1/E_0$. The implementation switches between cancellation-resistant forms of this expression. During the frozen-coefficient collision step, pair proposals occur at their actual continuous times. Each candidate's transverse momentum is rotated from the midpoint to the candidate time before evaluating its flux and scattering angle. If a collision occurs, both products are rotated back to midpoint representatives using their outgoing energies and destination fields. For constant field and momenta between events, this accounts for arbitrarily many turns between collisions. Step refinement tests the split integration when fields and longitudinal momenta vary. A second calculation evaluates all pair azimuths at the midpoint and a third suppresses gyromotion. The latter retains the magnetic field in the electromagnetic equations and energy accounting. [Orbit implementation](oriented.hpp) · [Independent rotation, basis and quadrature checks](ORIENTED_CHECK.json)

The largest individual free-flight half-step phase in the magnetic test suite is 48.78 radians when rounded upward. Event rephasing and the analytic inverse-energy integral handle full rotations. They still assume slowly varying field coefficients and a midpoint split. Across all 20 productions and controls the largest relative per-collision energy defect stayed below $1.0\times10^{-15}$, the global particle-energy ledger residual below $1.7\times10^{-12}$, and the largest relative discrepancy between the cached transverse magnitude and the two rotated components below $3.4\times10^{-14}$. These tests include one-quarter turns, charge-conjugate gyromotion, cyclic axis labels, and inverse-energy quadrature. [Independent checks](oriented_check.py)

After producing the native traces and the gyrotropic archive, reproduce the three-dimensional results from the study root with

```sh
python -I -B driven/oriented_check.py --output /tmp/ORIENTED_CHECK.json
python -I -B driven/oriented.py --native /tmp/native --gyro /tmp/DRIVEN.json --output /tmp/ORIENTED.json --jobs 2
```
