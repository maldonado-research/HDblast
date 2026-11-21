# Can the smooth wall fields cool the scattered charges?

In the tested dedicated-interaction pulse, smooth-field radiation would carry away roughly **one to two parts in a quadrillion** of the final charged-particle energy by model time 7.05, if each daughter obeys the ordinary local four-dimensional classical radiation law. A deliberately generous force estimate puts the fraction below **3.65 parts in ten trillion** across the four late collision runs. These results use the simulated charges, masses, fields and velocities. They narrow one possible cooling mechanism for this benchmark. [Recorded results](RADIATION.json)

## Calculation

Write $\tau=t_{\rm proper}/L$ for the shell observer’s dimensionless proper time, $a$ for the wall scale factor, $\mu=mL$, $\mathbf p$ for spatial momentum in units of the daughter mass and $\gamma=\sqrt{1+|\mathbf p|^2}$. The physical electric and magnetic fields in the local orthonormal shell frame are $\mathbf E=\mathbf F/a^2$ and $\mathbf B=\mathbf G/a^2$, using the stored conformal Maxwell fields. The dimensionless force coefficient $Q=10^7$ is the coefficient in the single-flavor Lorentz equation. With $\mu=4\times10^{12}$, $K=2.54$ and $S_0=1.0415592372880226\times10^{54}$, the assumed canonically normalized Abelian interaction gives $\alpha=3Q^2\mu^2/(4\pi S_0K)=1.443822089\times10^{-16}$. This coupling controls photon emission by each charge. The phenomenological resonant two-body scattering kernel has its own parameters. [Project scattering normalization](https://github.com/maldonado-research/HDblast/tree/main/research/HDBLAST_SCATTERING_20261121)

The Lorentz-invariant force measure and local classical radiation rate are

$$
\chi^2=\left(\frac{Q}{\mu}\right)^2
   \left(\big|\gamma\mathbf E+\mathbf p\mathbin{\times}\mathbf B\big|^2
           -(\mathbf p\cdot\mathbf E)^2\right),\qquad
\frac{d(E_{\rm emitted}/m)}{d\tau}=\frac{2}{3}\alpha\mu\,\chi^2.
$$

The local power is measured per unit shell-observer time. The derivative needs no extra time-dilation factor. For a nonzero momentum the solver evaluates the force invariant as $(Q/\mu)^2[E_\parallel^2+|\gamma\mathbf E_\perp+\mathbf p\times\mathbf B|^2]$. This avoids cancellation when a relativistic particle moves parallel to the electric field. Checks cover a charge at rest, constant parallel electric acceleration up to $|\mathbf p|=10^8$, pure transverse magnetic acceleration and crossed fields. At each endpoint, the script independently evaluates the original four-vector expression in extended precision. [Classical normalization in Ridgers et al., section 2.1](https://arxiv.org/abs/1708.04511) · [Implementation](smooth_radiation.py)

The integration follows each existing kinetic particle. On each time step it uses the converged midpoint electric field, the magnetic field at the particle push, and the momentum at the middle of the electric kick after the first magnetic half-rotation. Only the elapsed portion of a birth step contributes. The kinetic solver continues to exchange energy with the wall field and to scatter particles with its original kernel. The instrumentation accumulates the diagnostic on those trajectories. The local emitted energy is converted to conformal energy via $\sum_i w_i\int a(\tau)(2/3)\alpha\mu\chi_i^2\,d\tau$ with $w_i=D/(2N)$ for the representative charge-symmetric population. The comparison energy is the actual endpoint $\sum_i w_i a_f \gamma_{i,f}$. This is the weight convention of the archived gas ledger, including the rotated representative populations. The conformal weighting corresponds to freely redshifting emission on the shell.

An orientation-independent envelope replaces the force invariant with $(Q/\mu)^2\gamma_{\max}^2(|\mathbf E|+|\mathbf B|)^2$ on each particle step. The upper $\gamma$ is the larger value on either side of the electric kick. The magnetic half-rotations preserve momentum magnitude, and energy during a fixed electric kick is convex. The envelope therefore bounds the sampled estimate for every step and allows any momentum direction during that particle push. Both estimates use the solved fields and simulated time grid. A continuum bound would require independent bounds on the fields between force evaluations.

## Results at wall time 7.05

| Birth pairs | Time steps | Mean collisions per particle | Radiated fraction of final gas energy | Orientation-independent envelope | Largest sampled $\chi$ |
|---:|---:|---:|---:|---:|---:|
| 4,096 | 4,096 | 9.50 | $1.24\times10^{-15}$ | $3.62\times10^{-13}$ | 0.0248 |
| 8,192 | 4,096 | 9.55 | $1.95\times10^{-15}$ | $3.64\times10^{-13}$ | 0.0408 |
| 8,192 | 8,192 | 9.66 | $1.40\times10^{-15}$ | $3.64\times10^{-13}$ | 0.0316 |
| 16,384 | 8,192 | 9.45 | $1.48\times10^{-15}$ | $3.63\times10^{-13}$ | 0.0408 |

In the largest late run the integrated emission is $7.42\times10^{-13}$ conformal energy units against a final charged-particle energy of 501.55. The 10% of particles ending with the most energy contributed 22.9% of the calculated radiation. Even the largest ratio for an individual particle, measured against its own final conformal energy, was $2.75\times10^{-12}$. With the same smooth driving and collisions switched off in the 8,192-pair control, the total fraction was $2.60\times10^{-22}$. A separate direct four-force calculation on the largest endpoint snapshots gives $\chi_{\max}=0.00781950$ at 6.67 and 0.00113965 at 7.05, with emission rates relative to the gas energy of $1.8273\times10^{-15}$ and $1.3556\times10^{-16}$ per unit wall time, respectively. Scattering redirects some daughter momenta into larger smooth transverse forces. The Poisson event realization and its time discretization affect the tiny estimate. The orientation-independent envelope varies by less than one percent over the four late collision runs. [Radiation, endpoints and per-particle checks](RADIATION.json)

The diagnostic sums incoherent single-particle emission under a local 3+1-dimensional electromagnetic law, using the leading classical term with no recoil fed back. The strongest sampled forces have $\chi\simeq0.04$, with a late orientation-independent maximum envelope around 0.065. Spectral, quantum and radiation-reaction corrections have not been evolved. Annihilation, bremsstrahlung during individual close encounters, collective and stimulated emission, absorption and bulk wave propagation require their own dynamics. The large-scale radiation carried by the solved mean current already enters the Maxwell energy balance. Adding this diagnostic to that ledger as an independent loss would require resolving coherence, interference and which emitted modes the wall can retain.

## Reproduce

The study uses Python with NumPy and SciPy and a C++17 compiler. From a checkout containing the axial-current study and its two predecessor studies, generate the extended native dedicated history and run the postprocessor. Set the postprocessor path to its location in your checkout. Each instrumented run checks its final charged energy and event count against the archived trajectories.

~~~sh
python -I -B research/HDBLAST_AXIAL_CURRENT_20261121/fields/extract.py \
    --study research/HDBLAST_POLARIZATION_20261121 --cases dedicated_half --end 7.0505 \
    --output-dir /tmp/hdblast-radiation-input --output /tmp/hdblast-radiation-input/FIELDS.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -I -B /path/to/smooth_radiation.py \
    --study research/HDBLAST_AXIAL_CURRENT_20261121 \
    --polarization research/HDBLAST_POLARIZATION_20261121 \
    --scattering research/HDBLAST_SCATTERING_20261121 \
    --history /tmp/hdblast-radiation-input/dedicated_half.npz --output /tmp/RADIATION.json --jobs 2
~~~
