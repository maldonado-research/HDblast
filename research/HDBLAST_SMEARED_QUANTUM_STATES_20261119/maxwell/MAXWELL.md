# A transverse electromagnetic field through first passage

Free four-dimensional Maxwell theory has no direct coupling to a universal Weyl factor multiplying the spacetime metric. Its waves can still respond to an evolving, uneven geometry. We propagated one transverse sector on each of the two archived first-passage spacetimes, using the same symmetric, smooth energy-preparation windows as in the scalar calculations.

For the central window width $T=0.20$, the lowest measured occupation and the energy of one additional normalized launch excitation are

| First contact | Lowest occupation at $x_*-0.005$ | Lowest occupation at $x_*+0.005$ | Energy after / before for the added excitation | Bounds on the energy ratio of any classical Maxwell wave |
|---|---:|---:|---:|---:|
| Poles, $A_{\rm seed}=+0.001$ | $4.339\times10^{-8}$ | $6.455\times10^{-8}$ | 0.995751 | 0.995449 to 0.997550 |
| Equator, $A_{\rm seed}=-0.001$ | $3.483\times10^{-6}$ | $5.003\times10^{-6}$ | 0.996487 | 0.994459 to 1.012624 |

The occupations use the instantaneous Maxwell Hamiltonian on the named slices and depend on the prepared state and readout prescription. The classical bounds concern the whole electromagnetic field at all wavelengths. They follow from its local stress tensor. In particular, on the polar simulated geometry, the coefficients imply that the whole-sphere energy of every nonempty, source-free classical Maxwell test field declines across this interval. The equatorial bounds permit either sign. The simulation treats the background as prescribed. It gives no absolute vacuum stress or evidence of a thermal photon gas. [Numerical results](MAXWELL.json) · [Replay](maxwell.py)

## Axial field and the round limit

Take the spherical spacetime metric $ds^2=-N^2 dt^2+A^2 d\chi^2+R^2d\Omega^2$. An axial Maxwell potential with angular components $u(t,\chi)\epsilon_A{}^B D_BY_{\ell m}/\sqrt{\ell(\ell+1)}$ and vanishing radial and temporal components is transverse for $\ell\geq1$. With a unit-normalized real angular harmonic its reduced action is

$$
S=\frac12\int dt\,d\chi\left[\frac{A}{N}\dot u^2-\frac{N}{A}(u')^2-\frac{NA\ell(\ell+1)}{R^2}u^2\right].
$$

Each coefficient is unchanged when $N,A,R$ are multiplied by the same local positive function. Electric-magnetic duality maps axial solutions to polar solutions of the same vacuum dynamics. The table includes one axial polarization in a single real $\ell=2,m$ block with even radial parity. Neither angular multiplicities nor both polarizations have been included in its occupations.

For the archived unit-lapse metric, $A=a\exp(m+\sin^2\chi\,d)$ and $R=a\exp(m)\sin\chi$. Write $L=m+\sin^2\chi\,d$. The proper-time Hamiltonian in a fixed radial basis is $(p^T M^{-1}p+q^TKq)/(2a)$, with

$$
M_{ij}=\int e^L u_i u_j\,d\chi,\qquad
K_{ij}=\int\! \left[e^{-L}u'_iu'_j+\ell(\ell+1)e^{L-2m}u_i u_j/\sin^2\chi\right]d\chi .
$$

For a round spatial metric, the regular functions $u_n\propto\sin^{\ell+1}\!\chi\,C_{n-\ell}^{\ell+1}(\cos\chi)$ have $M=I$ and squared conformal frequencies $(n+1)^2$ for $n\geq\ell$. With conformal time $d\eta=dt/a$ the Hamiltonian is constant and every positive-frequency solution evolves by a phase. A round universe of any prescribed scale history has zero excitation in this preparation and readout. The quadrature and ODE on the round controls gave total occupation below $10^{-26}$. The same conclusion holds under a common conformal rescaling of the entire spacetime metric when comparing the conformally transported modes.

The uneven backgrounds have a measurable local obstruction to four-dimensional conformal flatness. Let $\delta$ be the difference of radial and tangential eigenvalues of the electric Weyl tensor for the archived normal observers. Spherical symmetry gives a vanishing magnetic part and $C_{abcd}C^{abcd}=16\delta^2/3$. In the single-parent Einstein–scalar background the spatial expression is

$$
\delta= {}^{(3)}R_{\hat r\hat r}-{}^{(3)}R_{\hat\theta\hat\theta}
+H_t(H_r-H_t)-\tfrac12(D_{\hat r}\Theta)^2.
$$

Using the relative proper-time derivative of the directional expansion rates gives independently
$\delta=-\partial_t(H_r-H_t)-(H_r-H_t)(H_r+H_t)+\tfrac12(D_{\hat r}\Theta)^2$.
At first passage the proper-volume mean of $C_{abcd}C^{abcd}/h^4$ is 0.0042882 and 0.22510, respectively, where $h$ is the reference hilltop expansion rate. The two calculations of the gap at the sampled times agree within $3\times10^{-13}h^2$. A nonzero invariant establishes that no coordinate change and single Weyl rescaling can reduce these geometries to a flat metric in the region where it is nonzero. The invariant alone gives no particle count.

## Gaussian preparation and measurement

We minimize the average of the whole-sphere physical energy divided by its evolving proper volume, using the weight $\exp[-1/(1-(t/T)^2)]$ for $|t|<T$ and zero outside. The weight is normalized to unit time integral. Time reflection of the archived solution makes the averaged launch quadratic form block diagonal, with field block $C$ and momentum block $B$. Both are positive. Its centered pure Gaussian minimizer has field matrix $\Omega=B^{-1/2}(B^{1/2}CB^{1/2})^{1/2}B^{-1/2}$. The launch mode columns are $Q=\Omega^{-1/2}/\sqrt{2}$ and $P=-i\Omega^{1/2}/\sqrt{2}$. We evolve these columns continuously to all measurements, retaining each distinct window preparation $T=0.10,0.20,0.30$. The modulus metric factor $f$ does not enter this Maxwell mode action.

At each measurement time, diagonalize $Kv_i=\omega_i^2 Mv_i$ with $v_i^TMv_j=\delta_{ij}$, using conformal frequencies $\omega_i$ and proper frequencies $\omega_i/a$. The absolute squares of $\sqrt{\omega_i/2}\,v_i^TMQ-i\,v_i^TP/\sqrt{2\omega_i}$, summed over launch columns, define the occupation in that instantaneous eigenmode. A separate positive-energy diagnostic adds one quantum in the first specified launch column to the same background state. The energy difference is $(P_1^\dagger M^{-1}P_1+Q_1^\dagger KQ_1)/a$. Diagonalizing an instantaneous Hamiltonian does not select a unique vacuum in the full inhomogeneous quantum theory. No mode sum over the full three-dimensional ultraviolet spectrum or semiclassical stress calculation is used here.

At $T=0.20$ the after-crossing excitation energies in this one sector for radial labels through 24 are $5.4552\times10^{-7}$ and $4.2487\times10^{-5}$ in inverse model-time units. The sum extends over more modes as the cutoff increases. Raising it from 24 to 36 changes the equatorial after-crossing energy to $4.2323\times10^{-5}$. The low-label occupations through 8 change by at most $4.7\times10^{-15}$ across windows and sampled times. At fixed cutoff 24, simultaneous time-step, geometry, angular and smearing-quadrature refinement changes those low occupations by less than $2\times10^{-17}$. Results at widths 0.10 and 0.30 and at radial cutoff 12 are also included in the JSON.

## Maxwell energy conservation on a closed slice

Write $\mathcal E=(|\mathbf E|^2+|\mathbf B|^2)/2$ for the Maxwell energy density measured by the normal observers. With $H_i$ the local principal expansion rates and $D=H_r+2H_t$, the diagonal Maxwell pressures are $p_i=\mathcal E-E_i^2-B_i^2$. Spatial flux integrates out on the closed slice. The total classical energy $E=\int\mathcal E\,dV$ obeys

$$
\dot E=-\int\sum_i H_i p_i\,dV
=\int\sum_i (-D+2H_i)\,\frac{E_i^2+B_i^2}{2}\,dV .
$$

Each directional energy contribution is nonnegative. Consequently $E(t_+)/E(t_-)$ lies between the exponentials of the time-integrated spatial minima and maxima of the coefficients $-D+2H_i$. For an isotropic, homogeneous FLRW geometry this gives exactly $E(t_+)/E(t_-)=a(t_-)/a(t_+)$. On the uneven geometries we included both poles and the equator in extremum searches and polished interior extrema. Refining spacetime, extrema and time integration changes the reported ratios by less than $2\times10^{-13}$. The inequalities hold for exact classical solutions. The tabulated coefficients are evaluations on the numerical background and carry no interval certificate for the continuum Einstein evolution. Positivity also covers a single positive added excitation relative to the same Gaussian state. It does not provide a bound for general subtractions between different quantum states.


## How much initially seeded light can remain important?

The same stress-tensor inequality covers the entire simulated interval from the bounce to first passage. For a positive pre-existing source-free Maxwell test field, divide the whole-sphere energy by the simultaneous total energy of the parent scalar, $r=E_\gamma/E_\phi$. Both integrals use the physical volume of their own slice. In the code's canonical normalization $E_\phi=\int [\dot\Theta^2/2+|D\Theta|^2/2+3u(\Theta)]\,dV$, including the full potential.

| First contact | Bounds on $E_\gamma(x_*)/E_\gamma(0)$ | Round reference $a(0)/a(x_*)$ | $E_\phi(x_*)/E_\phi(0)$ | Bounds on $r(x_*)/r(0)$ |
|---|---:|---:|---:|---:|
| Poles | 0.891408 to 0.902657 | 0.891540 | 1.367316 | 0.651940 to 0.660167 |
| Equator | 0.847386 to 0.950309 | 0.848452 | 1.455034 | 0.582382 to 0.653118 |

For instance, a classical photon seed carrying 0.1% of the parent energy at the bounce would carry at most about 0.0661% of the parent energy at polar first contact or 0.0654% at equatorial first contact. Even where the short equatorial bracket permits growth, the longer bounce-to-contact bound requires a net loss. These bounds compare endpoint energies. They apply on the two specified Einstein–scalar geometries in the test-field approximation. They neither specify a photon seed nor include production by a charge current, a gauge-kinetic coupling, a quantum vacuum stress, or changes to the geometry caused by the photons. Joint geometry, angular-extremum and time-integration refinement changes the tabulated bounce ratios by less than $2\times10^{-10}$.
