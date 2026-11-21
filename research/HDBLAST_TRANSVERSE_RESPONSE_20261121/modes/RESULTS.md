# Can the changing plasma amplify a transverse electromagnetic disturbance before the early thermal deadline?

We followed electromagnetic waves and the particles born during the actual archived startup. At the 1% absolute-work marker, the field produced by the particle current is far too small to meet the earlier sideways-impulse requirement when driven by a weak initial wave or by ordinary Poisson fluctuations in the source. A prescribed incoming wave can itself exert a sideways force. The size of that pre-existing wave is a separate physical input. [Uniform bound on Poisson birth fluctuations](NOISE_BOUND.md) · [All numerical runs](MODES.json)

| From first birth to 1% of whole-run absolute work | Dedicated background | Shared background |
|---|---:|---:|
| Elapsed proper wall time, in units of $L$ | 0.02692 | 0.03186 |
| Required early RMS extra transverse impulse, in $mc$ | 541.65 | 674.10 |
| Smallest necessary initial RMS electromagnetic seed in the sampled modes, counting the entire force | 0.00401 | 0.00414 |
| Smallest necessary seed when crediting only fields produced by the perturbed current | 0.714 | 0.505 |
| Largest terminal total-field amplitude divided by initial seed amplitude on this scan | 0.9825 | 0.9797 |
| Largest impulse upper estimate from a unit coherent fractional sinusoid in the birth rate, in $mc$ | 66.29 | 116.45 |

The static-guide magnitude, axial electric acceleration, ongoing births and scale factor all evolve during each run. The field entries are in the preceding wall solver's units, where the early parallel electric field reaches 0.207 and $Q=10^7$. Initial amplitude means the spatial RMS of $\sqrt{|\delta\mathbf B|^2+|\delta\mathbf E|^2}$. The quoted initial seeds are **necessary** values obtained by dividing the previous impulse requirement by an upper bound on the linear impulse per seed. A seed above that limit still has to demonstrate actual turning and the desired distribution. In particular, extrapolating the current-produced result to seed fields larger than the original parallel field exceeds a small-perturbation approximation. The scan uses 27 comoving wavenumbers from 16 to 32,768 in each channel. The displayed extrema are over these modes. The Poisson comparison separately covers a band with a finite cutoff. Its three-dimensional extension includes oblique and longitudinal modes and the permitted transverse birth momenta. [Earlier impulse benchmark](https://github.com/maldonado-research/HDblast/blob/6a62a0798713dbe28d0d670376896893cd8d79a1/research/HDBLAST_KINETIC_TEST_20261121/ballistic/KICK.md)

## The transverse field and particle equations

We take a local plane tangent to the shell, with $z$ along one stream's archived accelerating field and guide field, and $x$ along a comoving wavevector $k$. The three equal sign-paired streams are normalized so the zero-wavenumber response of the symmetric projection uses the preceding $J=Q(D/a^3)\langle v_z\rangle$ and gauge coefficient $K=2.54$. Equivalently, giving each gauge direction one third of the total kinetic coefficient and matter gives the same ratio of current to field coefficient. We treat the short-wavelength tangential perturbation as a local four-dimensional Maxwell field propagating on the prescribed scale factor. Its coupling to a perturbed five-dimensional bulk has not been derived here.

Every new neutral pair puts its two charges at the same point. At birth their axial perturbation momenta vanish. In the baseline both have zero initial momentum. The archived electric impulse then gives the positive charge $p_z/m=[I(t)-I(t_b)]/a(t)$, with $\dot I=QaX$ and negative charge its mirror. Set $v=p_z/(m\gamma)$, $\gamma=\sqrt{1+(p_z/m)^2}$ and $\Omega=QB_0/\gamma$. Perturbations are Eulerian fields sampled on unperturbed cohort orbits at linear order. In the charge-symmetric filamentation sector both signs develop the same $x$ velocity and spatial density modulation. Their $y$ and $z$ momenta have opposite signs. The two perturbed charge densities cancel, while their transverse and axial currents can drive a wave. A separate electrostatic charge-imbalance sector is possible for different initial data.

Let $r=a_0/a$, $\widetilde f=(a/a_0)^2 f$, $P_i=a\delta p_i/(a_0mQ)$, $d=\delta n/(Qn)$, and $S=Q^2/(K a a_0^2)$. Use cosine amplitudes for $B_y,E_y,P_x,P_y$ and sine amplitudes for $E_z,B_z,P_z,d$. The other quadrature obeys the same real system with shifted phase. For the positive charge of each axial cohort,

$$
 \dot P_x=-rv\widetilde B_y+\Omega P_y,\qquad
 \dot P_y=r\widetilde E_y-\Omega P_x,\qquad
 \dot P_z=r\widetilde E_z,\qquad
 \dot d=(k/a)rP_x/\gamma.
$$

With $dD_b$ the increment in comoving rest from the whole sign and direction population, the transverse Maxwell equations in this effective normalization read

$$
 \begin{aligned}
 \dot{\widetilde B}_y&=(k/a)\widetilde E_z, &
 \dot{\widetilde E}_z&=-(k/a)\widetilde B_y-S\int dD_b[rP_z/\gamma^3+v d],\\
 \dot{\widetilde B}_z&=(k/a)\widetilde E_y, &
 \dot{\widetilde E}_y&=-(k/a)\widetilde B_z-S\int dD_b\,rP_y/\gamma .
 \end{aligned}
$$

Only already-born cohorts enter each integral. Compression changes the axial current through $vd$, even when the number of positive and negative charges at each point is equal. The guide field links the two wave polarizations through $P_x,P_y$. A full Maxwell–Vlasov calculation would additionally evolve the angular spread of the injected distribution, collisions, the sources and the bulk perturbation consistently. This cold probe imposes the archived background; the kinetic background current predicted by the ballistic cohorts can differ from the archived current. The labels dedicated and shared describe the two supplied histories. The perturbation calculation uses collisionless propagation in both.

## Distinguishing externally supplied force and field produced by the current

Each of four independent electromagnetic initial data defines one column of a field transfer matrix $T(t)$. Without the charge current, the solution $T_0(t)$ rotates the two vacuum polarizations through angle $\int k\,dt/a$ and redshifts their physical amplitudes by $(a_0/a)^2$. The source-induced fields are $(T-T_0)f(0)$ for exactly the same initial seed and boundary prescription. This difference includes the phase shift due to the plasma. It can be nonzero even if the energy of the field decreases.

For a tagged unscattered daughter the added transverse force in units of $m/L$ is $Q(-v\delta B_y,\delta E_y)$ at the spatial wave phase. The parallel guide rotates the transverse momentum. As in the earlier kick target, its rotation contributes no new transverse magnitude. Since $|v|\le1$, the largest singular value of the physical $(B_y,E_y)$ two-row block bounds the RMS force for any fixed initial polarization or phase in that mode. Integrate this upper envelope from each cohort's birth with weight $a(t)/a_*$, and then take its root mean square over the comoving number already present. Use the physical $T$ for the full columns, and $T-T_0$ for the source-induced columns. The procedure gives the seed every advantage by allowing the most effective instantaneous polarization at each time. In the numerical output, `terminal_field_transfer` is the full physical 4×4 matrix. The separately accumulated 4×4 momentum Gram matrix gives the exact optimal terminal *net* transverse response within the discretized linear model. The net response at large $k$ is smaller than the absolute-force integral because oscillations can undo a previous deflection. [Implementation](evolve.cpp)

A fifth input supplies an initial spatial source modulation $\delta n/n=\sin(kx)$ with zero fields. It produces a current through the differing density of pairs accelerated along $z$. The largest scanned response per unit modulation occurs at $k=1792$ in the dedicated example and $k=1280$ in the shared example. The entries in the table are linear transfer coefficients; order-unity density modulation is outside the strict small-amplitude limit.

For the shot-noise experiments, assign independent Fourier cosine/sine coefficients to the pairs born in disjoint time intervals at uniform random positions within one wavelength cube. A bin with $N$ expected neutral pairs has variance $2/N$ in the fractional cosine density coefficient. Propagating one covariance-factor column per independent bin determines the field variance exactly for that discretized source. This uses one pair for the co-located two charges. With one source column per each of 512 cohorts at $k=32768$, the RMS absolute-force envelopes are $2.83\times10^{-12}$ and $3.76\times10^{-12}$ in $mc$. Net RMS transverse momenta there are $2.12\times10^{-12}$ and $2.81\times10^{-12}$. Coarse bins can underestimate high-frequency noise. Increasing from 256 bins and cohorts to 512 changes the absolute-force estimates at this mode by less than 0.01%. The independent continuous-source ceiling and an explicit bound for the whole band in a common box do not depend on this bin convergence. [Continuous shot-noise proof](NOISE_BOUND.md)

## Numerical and physical scope

The solver advances each vacuum half step by exact phase rotation at midpoint scale factor, predicts the midpoint current response using the preceding current, advances each cohort through a time-centered gyro rotation with the applied fields, and inserts the updated midpoint current into Maxwell's equations. Close to birth it integrates $1/\gamma$, $1/\gamma^3$ and the first moment of $1/\gamma^3$ over each approximately affine interval in $p_z$. This resolves the rapid transition from nonrelativistic newborn to relativistic stream. Each source time bin contributes the actual rest increment of the pinned history, and its cohort is activated at the midpoint.

At $k=32,256,2048,8192$ we varied the cohort count from 256 to 1,024, the number of field steps from 2,048 to 8,192, and the time resolution of the original wall solver. Between the baseline (512 cohorts, 4,096 steps) and finest checks, the current-induced impulse ceiling changed by under 0.08% across these probes. Setting each birth's initial axial momentum to either extremal value $\pm0.32016m$ changed the tested current-induced impulse and coherent-birth response by less than 0.03%. Removing the guide made a small difference during this early interval. The full spherical birth distribution carries transverse momenta as well. A collisionless birth with $p_\perp/m\le0.32016$ can traverse under 0.0064 rad of unperturbed Fourier phase by the deadline even at $k=32768$, given the archived longitudinal acceleration. This bounds free spatial dephasing before scattering. It does not control all finite-pitch-angle contributions to the current.

Direct matrix exponentiation of five independently specified constant backgrounds checks all four polarizations against vacuum, plasma at rest, static growing beams and a strong-guide stable example. The largest fine-grid relative field-matrix error over the control suite is $1.4\times10^{-6}$. The finite-time propagator here tests energy exchange over an accelerating pulse and ongoing injection. A growing normal mode on a time-independent reference is a different measurement. [Independent matrix checks](MODE_CHECK.json) · [Checks implementation](check.py)

The transfer-matrix simulations use strictly transverse wavevectors. The analytical Poisson bound separately treats the full three-dimensional band. Curvature, transverse variation of the original driving field, a wavenumber-dependent bulk boundary condition, collisions and nonlinear saturation are outside the local tangent calculations. Independent Poisson birth positions and specified initial wave seeds are test inputs. A spatially correlated source or externally supplied transverse field requires an explicit model for its origin, spectrum and energy.

From the study root regenerate the products with

```sh
python -I -B modes/run.py --history history --output /tmp/MODES.json --jobs 2
python -I -B modes/check.py --history history --output /tmp/MODE_CHECK.json
python -I -B modes/noise_bound.py --history history --output /tmp/NOISE_BOUND.json
```
