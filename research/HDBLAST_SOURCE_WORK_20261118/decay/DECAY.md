# A field that loses energy into a massless reservoir

On the archived B1 wall trajectories, replacing the local crossing event by continuous quantum-field evolution changes the predicted daughter radiation. The quiet rolling solution produces radiation in its own right. In the specified momentum band and for the nominal decay strength, the near positive packet reduces the final daughter density by **0.178%** relative to quiet. The near negative packet raises it by **1.020%**. Moving the compensating part of the positive packet beyond the readout's causal past gives **+0.092%**. [Numerical results](CHECK.json) · [Prior wall histories](https://doi.org/10.5281/zenodo.24355584)

### Energy, force and pressure

Use dimensionless proper time $s=H_0\tau$, normalized $a_{*,0}=1$ at the historical quiet crossing and $d\nu=k^2dk/(2\pi^2)$ for one real scalar with $0\le k\le8$. The signed quantity $\mu=100(\phi-0.5)$ gives mode energy $E=(k^2/a^2+\mu^2)^{1/2}$. In physical field variables the canonical pair is $(\chi,\Pi=a^3\dot\chi)$. Define the Gaussian moments $u=\langle\chi^2\rangle$, $v=\langle\Pi^2\rangle$, $c=\langle\chi\Pi+\Pi\chi\rangle/2$. At each time the oscillator ground has $u_0=(2a^3E)^{-1},v_0=a^3E/2,c_0=0$. Initialize the field in that ground state at $s=0.8$ and the massless reservoir empty. This preparation differs from the adiabatic state used for the preceding release's rolling particle comparisons.

Subtracting $E/2$ from the finite-mode Hamiltonian gives the nonnegative comoving oscillator energy $U_k=E n_k=[(v-v_0)/a^3+a^3 E^2 (u-u_0)]/2$. For this specified prescription the comoving mass work rate is $S_k=a^3\mu\dot\mu(u-u_0)$, and the comoving pressure is $P_k=(v-v_0)/(2a^3)-a^3(\mu^2/2+k^2/(6a^2))(u-u_0)$. These are conjugate derivatives of the same subtracted energy before taking expectation values. Subtracting only the vacuum energy from the output would give an inconsistent balance. The band and its subtraction define an effective model. A covariantly renormalized stress over all momenta requires additional physics.

### An explicit decay model

Assume a zero-temperature, memoryless amplitude-damping channel in the instantaneous physical oscillator, with $\gamma=y^2\mu^2/(8\pi E)$ in units of $H_0$. This takes the scaling of the time-dilated decay width from the original B1 gas. The Gaussian channel drives $u,v,c$ toward $u_0,v_0,0$ at rate $\gamma$. In the covariance equations, the Hamiltonian terms are $2c/a^3$, $-2a^3E^2c$, $v/a^3-a^3E^2u$ respectively. The dissipative energy loss per comoving mode is exactly $\gamma U_k$. Feed it to a homogeneous massless component with physical pressure $\rho_r/3$. Writing $Z=a^4\rho_r$, the two balances are

$$
\dot U_k=S_k-3hP_k-\gamma U_k,\qquad \dot Z=a\int \gamma U_k\,d\nu,\quad h=\dot a/a.
$$

The summed energy is $\int U_k\,d\nu+Z/a$. Its change is scalar work minus the pressure work of both components. The loss rate cancels between them. This construction quantitatively tests one decay model on the source histories. A microscopic treatment of a rapidly changing mass can select a different dissipative basis. The reservoir evolves energy and pressure only. Its phase-space distribution and temperature are unspecified. [Related open-system literature](../literature/SOURCES.md)

### Numerical checks at $s=2.10$

| Fixed shell history | Scalar work | Field energy | Daughter energy | Daughter density |
|---|---:|---:|---:|---:|
| Quiet | 156.2444 | 83.5836 | 69.6923 | 36.92618 |
| Near positive | 155.9431 | 83.4162 | 69.5638 | 36.86050 |
| Near negative | 157.7523 | 84.3658 | 70.3909 | 37.30270 |
| Distant positive | 156.3322 | 83.6150 | 69.7472 | 36.96008 |

The first three values are per unit normalized comoving volume in $H_0^4$ units. The last is a physical density in $H_0^4$. On the quiet history the field and radiation expansion-work terms are +0.4961 and −3.4646 in the same comoving units. Integrated decay transfers 73.1569 to the daughters before their redshifting. The part of the scalar work obtained from occupations alone is 125.0441. Phase correlations account for the remaining 31.2003.

An adaptive eighth-order Runge–Kutta solver integrates the covariance and each work independently. The largest absolute residual among the three balances at the reported cuts is below $8\times10^{-6}$ across all runs. Increasing the Gauss momentum order from 24 to 48, changing the history resolution, and varying $y$ between 0.5 and 2 preserve the three response signs. $y=0$ transfers zero radiation. All comparisons retain the physical momentum dependence of $E$ and stay on the published screened band. A second integrator composes exact oscillator rotations and positive Gaussian-loss channels on midpoint-held intervals. It recovers the primary outputs to better than nine parts per million and reproduces the analytic decay of an initially occupied static mode. [Main solver](reservoir.py) · [Independent split evolution](split.py) · [Energy and convergence tests](check.py)
