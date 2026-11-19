# Which Gaussian state minimizes sampled energy on an uneven bounce?

A short-wavelength completion for an entire quantum field requires control over all directions and wavelengths. This calculation instead solves the energy-minimization problem without a diagonal-mode approximation in one specified finite sector of the actual inhomogeneous HDblast spacetime. One real spherical harmonic of angular degree two and the reflection-even radial labels $n=2,4,\ldots,N$ describe the daughter field. The background is prescribed. The two trajectories have initial parent amplitude $A=+0.001$ and $-0.001$ on the potential whose minimum lies at $\phi_m=-0.5$. They first cross at the two poles or at the equatorial surface, respectively.

## Which quantity is minimized?

In the same original canonical field variables as the archived energy study, write $y=(q,p)$, with Hamiltonian $H=y^T G(t)y/2$ and $G=\mathrm{diag}(K,M^{-1})$. The matrices $M,K$ are the mass and stiffness obtained by integrating the evolving three-metric and $Z=f(\Theta_*)/f(\Theta)$ against the regular round-sphere basis. Off-diagonal entries mix radial labels. Write $U(t)$ for the symplectic map carrying canonical variables at the bounce to time $t$. In a fixed symplectic scaling of the initial variables the same formulae hold with $U(0)$ equal to the scaling matrix. Let $V(t)$ be the whole-slice proper volume. We minimize the mean of the spatially averaged original field energy,

$$
F_T=\int_{-T}^T w_T(t)\, \frac{U(t)^TG(t)U(t)}{V(t)}dt,\qquad
\langle\bar\rho\rangle_T=\tfrac12\operatorname{Tr}(F_T C).
$$

Here $C$ is the symmetrized covariance of the canonical variables at the bounce and $\hbar=1$. The normalized smooth weight is proportional to $\exp[-1/(1-t^2/T^2)]$ on $|t|<T$ and zero elsewhere, in the Einstein-frame geodesic proper time. These are the same window functions and density convention as in the homogeneous calculation. The global spatial average makes no preference for a polar or an equatorial observer.

The background extends by time reflection about its time-symmetric initial slice. Configuration and momentum transform under $R=\mathrm{diag}(I,-I)$, so $U(-t)=RU(t)R$ in a scaling which commutes with $R$. Since $V,G$ are even, the two halves of the integral cancel the mixed configuration–momentum block of $F_T$. This gives $F_T=\mathrm{diag}(D,B)$ with positive symmetric matrices. Both blocks retain coupling between different radial labels. Defining $L=(B^{1/2}DB^{1/2})^{1/2}$, the minimizing centered pure Gaussian has

$$
 C_{qq}=\tfrac12B^{1/2}L^{-1}B^{1/2},\qquad
 C_{pp}=\tfrac12B^{-1/2}LB^{-1/2},\qquad C_{qp}=0 .
$$

This follows by a canonical transformation bringing the sampled positive harmonic Hamiltonian to independent oscillators. Its ground state uniquely minimizes the sampled energy over the finite-sector states. For an independent check, the code uses the general Williamson expression $F_T^{-1/2}|iF_T^{1/2}JF_T^{1/2}|F_T^{-1/2}/2$ with canonical symplectic matrix $J$. It checks the purity equation $CJC=J/4$, the nonnegative cost of evaluating a state in the other window's objective, and the optimal improvement over the best product of independent oscillators in the *original* round radial labels. Constant rescalings used internally preserve commutators and are undone before reporting physical field and momentum correlations. The nonzero time interval specifies a preparation that can be varied between runs.

## Following differences through first passage

We vary $T=0.1,0.2,0.3$ and compare each of the other preparations with $T=0.2$ on the same geometry. For $\Delta C=C_T-C_{0.2}$, $\Delta E(t)=\mathrm{Tr}[U^T G U\Delta C]/2$ is a total energy difference in the one angular block, in inverse model-time units. Its division by the spatial volume gives the difference in that block's contribution to the mean energy density. The same common subtraction would cancel between the two states. An indefinite covariance difference permits either sign. Its energy cannot be interpreted as the density of produced radiation.

During the interval $[x_*-0.005,x_*+0.005]$ we directly integrate $U^T\dot G U$ with the partial derivatives of $G$ from the frozen work calculation. Varying $Z$ at fixed metric gives the direct coupling work. Varying the metric at fixed $Z$ gives the geometry work. The source sign is positive for work into the daughter. We evaluate the endpoint energy change and check the two work terms against it as a full matrix identity, before contraction with any chosen state difference. Setting $Z=1$ on the same geometry makes the coupling work vanish exactly. That control has its own minimizing states. The energy difference is between quantum preparations, and the background is the previously published single-field geometry.

At $N=160$ there are 80 radial oscillators in the chosen angular block. In the table, a positive difference means the indicated preparation has a larger value than the reference preparation. The two work entries have units of inverse model time. The last column is the difference in spatial mean density at first passage in inverse model time to the fourth power.

| First crossing | Preparation minus reference | Direct coupling work | Metric work | Mean density difference at crossing |
|---|---|---:|---:|---:|
| Poles | $T=0.10$ minus $0.20$ | $+1.7422\times10^{-5}$ | $-1.7188\times10^{-5}$ | $-2.5128\times10^{-5}$ |
| Poles | $T=0.30$ minus $0.20$ | $-2.0952\times10^{-5}$ | $+3.2803\times10^{-5}$ | $+2.5791\times10^{-5}$ |
| Equator | $T=0.10$ minus $0.20$ | $-3.9945\times10^{-5}$ | $-1.7693\times10^{-5}$ | $-4.3385\times10^{-5}$ |
| Equator | $T=0.30$ minus $0.20$ | $+5.1201\times10^{-5}$ | $+3.6466\times10^{-5}$ | $+4.8280\times10^{-5}$ |

The signs of both work terms in all four rows persist from $N=32$ through $160$ at the sampled cutoffs. At $N=16$ the narrow-window equatorial coupling work is positive and the narrow-window polar metric work is positive. From $N=128$ to $160$, the largest fractional shift of an individual work term is 1.9 percent, for the narrow-window equatorial coupling work. The polar narrow-window terms nearly cancel. Their residual at $N=160$ is $+2.34\times10^{-7}$ and was $+1.57\times10^{-6}$ at $N=96$. This cancellation makes its size especially sensitive to which harmonics are retained. These finite-cutoff comparisons give no error bound for the omitted radial or angular modes.

For $T=0.20$, the largest absolute correlation between distinct bounce field coefficients at $N=160$ is $4.43\times10^{-4}$ for the polar seed and $4.71\times10^{-4}$ for the equatorial seed. Allowing the modes to correlate lowers the sampled density by $3.04\times10^{-5}$ and $3.55\times10^{-5}$ relative to the optimal independent-oscillator preparation, respectively. An exactly round closed de Sitter reference using independent analytic Klein–Gordon mode functions at $N=16$ reproduces the full covariance to $2.4\times10^{-14}$ and the sampled density to $1.6\times10^{-15}$.

The calculation evaluates background spectral degrees $64$ and $96$, the latter with 256 spatial nodes. At fixed daughter cutoffs $48$ and $160$ we refine the spacetime, its passage bracket, the radial quadrature, the proper-time quadrature and the evolution step together. The largest of the monitored signed work, energy-change and density shifts at $N=160$ is below $2\times10^{-13}$. Over the seven reported radial truncations the absolute full-matrix endpoint work discrepancy stays below $1.5\times10^{-10}$. It is measured before contraction with the covariance difference. Each cutoff uses a freshly minimized state in its own finite space. Establishing an inhomogeneous all-sector state with the local short-distance behavior required by quantum field theory remains a separate problem.

[Full results and controls](PATCH.json) · [Calculation](patch.py)
