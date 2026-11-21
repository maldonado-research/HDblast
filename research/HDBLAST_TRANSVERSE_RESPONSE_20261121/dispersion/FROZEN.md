# Is the early particle stream unstable if its properties are held fixed?

At the point when 1% of the archived run's absolute electrical work has occurred, a local electromagnetic perturbation can grow in the cold, collisionless model used here. The rate depends strongly on how the particles' axial speeds are distributed. For the dedicated interaction, replacing the newborn particles of different ages with a single beam at exactly the same density and charge current overestimates the instantaneous growth rate by a factor of **7.30** at comoving wave number 8,192. The shared interaction gives a factor of **7.51**. The background has very little time to realize those instantaneous rates. [Numerical results](FROZEN.json)

| At 1% electrical work, finer wall time step | Dedicated channel | Shared channel |
|---|---:|---:|
| Age-resolved growth rate at comoving wave number 8,192, in units of $1/L$ | 13.71 | 13.70 |
| Single speed giving the same density and current, same wave number | 100.09 | 102.93 |
| Comoving wave number at which one wavelength fits inside the light-travel distance since the first birth | 2,416 | 2,043 |
| Integrated ceiling on the instantaneous growth rate, allowing every wave number | 0.4964 | 0.6025 |

Comoving wave numbers here denote local plane waves with physical wave number $k=\kappa/a$. A pre-existing wave need not fit into the light-travel distance *since injection*. An extended pattern built causally after injection would have to respect its coherence scale. No spatial harmonics on the full closed wall were solved. The growth-rate integral in the table constrains the number of e-foldings one would infer by integrating instantaneous normal-mode rates. Evolving eigenvectors, electromagnetic wave propagation, birth of new particles and forced transient responses require a time-dependent calculation.

## Charge symmetry and the guide field

Consider one gauge axis with the inherited full-current wall normalization, $c=m=1$, field coefficient $K=2.54$ and charge-to-mass ratio $Q$. A cohort of total wall-frame rest density $\rho_j$ has opposite charges moving at opposite velocities $\pm v_j\hat z$. The pair is electrically neutral while its two contributions to the axial electric current add. Its homogeneous current $J=Q\sum_j\rho_j v_j$ belongs to the accelerating background maintained by the wall and bulk equations. The perturbation considered here has nonzero wave number along $x$, perpendicular to the guide field $B\hat z$. A beam of comoving positive and negative charges would have a different current response.

For a fixed snapshot, define $\alpha_j=Q^2\rho_j/(K\gamma_j)$ and $\Omega_j=QB/\gamma_j$. In the charge-symmetric filamentation sector, the two signs have the same transverse displacement $\xi_j$ along $x$. This displacement changes their number densities equally. Their electric charge perturbations and $x$ currents cancel, but the axial currents from the density perturbation add. The guide field rotates the transverse motions oppositely and couples them to the $y$ polarization. Use vector-potential polarizations $A_y$ and $A_z=iZ$ and Fourier convention $\exp(ikx)$. Set $\chi_j=\gamma_j\xi_j/Q$. Conserved transverse canonical momentum gives $\delta p_{yj}=-Q(A_y+\Omega_j\chi_j)$ for the positive sign when evaluating a growing mode. The perturbation equations reduce to

$$
\ddot\chi_j=-\Omega_j^2\chi_j-\Omega_j A_y-kv_j Z,\qquad
\ddot A_y=-(k^2+\sum_j\alpha_j)A_y-\sum_j\alpha_j\Omega_j\chi_j,
$$

$$
\ddot Z=-(k^2+\omega_\parallel^2) Z-k\sum_j\alpha_j v_j\chi_j,
\qquad \omega_\parallel^2=\sum_j\alpha_j/\gamma_j^2 .
$$

The scalar potential and the $x$-electric polarization decouple by charge symmetry in this cold axial sector. The gauge kinetic factor and matter weights reproduce the wall's longitudinal fixed-momentum susceptibility. An alternative transverse field theory with different kinetic coefficients would require its own normalization.

In the variables $(A_y,Z,\sqrt{\alpha_j}\chi_j)$ the squared-frequency matrix $V$ is real and symmetric. The cohort block is diagonal with entries $\Omega_j^2$. The two field diagonal entries are $k^2+\sum\alpha_j$ and $k^2+\omega_\parallel^2$, and the field-cohort entries are $\sqrt{\alpha_j}\Omega_j$ and $\sqrt{\alpha_j}k v_j$. For $B\ne0$ and $k>0$, eliminating the cohort block and the $y$ polarization proves that exactly one growing mode exists if

$$
k^2\big(KB^2-P_{zz}\big)+ KB^2\omega_\parallel^2-J^2/K < 0,\qquad P_{zz}=\sum_j\rho_j\gamma_jv_j^2.
$$

If the expression is positive, all modes in this frozen sector are oscillatory. Equality admits a zero mode. Both examples have $P_{zz}/(KB^2)>628$ and $(J/KB)^2/\omega_\parallel^2>1,241$ at the 1% marker on the finer mesh. The pressure exceeds the guide-field tension and the current term also favors growth. For the numerical rate the program computes the negative eigenvalue $-\lambda^2$ by eliminating the same block from $V+\lambda^2 I$. The zero-guide limit can be handled in that equation for $\lambda>0$.

The quadratic form of $V$ can be written as a sum of $k^2 A_y^2$, $\omega_\parallel^2 Z^2$, $\sum_j\alpha_j(A_y+\Omega_j\chi_j)^2$ and $(kZ+\sum_j\alpha_jv_j\chi_j)^2$, minus $(\sum_j\alpha_jv_j\chi_j)^2$. Cauchy–Schwarz bounds the final square by $(\sum_j\alpha_jv_j^2)\sum_j\alpha_j\chi_j^2$. Hence every frozen growth rate obeys $\lambda\le\sqrt{\sum_j\alpha_jv_j^2}$ for any $k$ or axial guide strength. This bound avoids optimizing over a hand-picked wave-number band.

## Distribution and numerical checks

Each cohort starts with zero axial and transverse momentum, then receives the axial electric impulse of the recorded history. For a birth at $t_b$, its physical axial momentum is $[I(t)-I(t_b)]/a(t)$. The weight is $dD/a^3$, where $D$ is accumulated comoving rest in the parent run. We use the recorded array of all accepted steps from the first birth. Within each cell we linearly interpolate $I$ and $D$. The moment integrals then have elementary primitives over a uniform axial momentum interval. The dispersion integral uses rapidity, with $dp/\gamma=dr$, to resolve recently born particles even when their momenta change by many rest-mass units in a step. Density and mean velocity alone cannot recover the transverse inertia of this mixture.

At the 1% marker the two wall time resolutions give integrated rate ceilings differing by at most $1.1\times10^{-5}$. Doubling the time quadrature within a fixed history changes these integrals by less than $3\times10^{-5}$. Assigning every birth a common axial momentum of $+0.32016$ or $-0.32016$ instead of zero changes either ceiling by less than 0.06%. This sensitivity check samples two alternatives. The calculation omits the actual transverse birth spread, subsequent elastic encounters, spatial gradients of the background, and backreaction of the perturbations. Direct diagonalization for 75 synthetic multistream matrices agrees with the reduced dispersion calculation to within $4\times10^{-13}$ in growth rate. [Reproduction script](frozen.py)

From the study directory run `OPENBLAS_NUM_THREADS=1 python dispersion/frozen.py --history history --output /tmp/FROZEN.json`.
