# The short waves during and after the change in coupling

On each archived round background, a rolling field changes the electromagnetic coefficient $I$ from one constant value to another. We kept the parent motion fixed to isolate the response of the photon modes. The modes from harmonic 49 through 768 carried **24.09%** of all the outgoing photon energy on the positive-minimum trajectory and **42.18%** on the negative-minimum trajectory, for $\Delta\log I=0.1$. Here the denominator is the earlier full numerical sum through harmonic 1536.

During the transition, much of the short-wave stress follows the field's current speed and acceleration. After the transition these local terms vanish, while a pair-excited photon state still carries energy. Matching low-frequency response therefore leaves a separate outgoing radiation contribution to measure. This calculation gives a diagnostic for a prescribed trajectory. The gravity and parent equations in the coupled calculation must determine their own trajectories.

## Exact finite-band observables

We match the earlier smooth gauge law $\log I=\delta S(u)$, with $u=(f-0.2)/0.6$ and $f$ the fraction of the directed displacement from launch to the first minimum. The switch derivative $S_u$ is proportional to $\exp[-1/(u(1-u))]$ on $0<u<1$ and vanishes elsewhere. Define $s=(\log I)'$, $g=\partial_\Theta\log I$, $d\eta=dx/a$, and $x=H_0t$. Here primes mean conformal-time derivatives. In vacuum at launch the occupation $n_L$ and anomalous covariance $m_L$ vanish. Their equations are

$$ n_L'=2s\,\mathrm{Re}\,m_L,\qquad m_L'+2iL m_L=s(1+2n_L),\qquad n_L(n_L+1)=|m_L|^2. $$

The numerical solver propagates phase-removed Bogoliubov amplitudes so that the small outgoing occupation is evaluated as $|\beta_L|^2$. With spectral weights $w_L=L(L^2-1)/\pi^2$ and window $R_L$, set

$$D=\sum_L R_L w_Ln_L,\qquad J=2\sum_LR_L w_L\mathrm{Re}\,m_L,\qquad D'=sJ.$$

The comoving photon excess is $D=a_{\rm physical}^4\rho_\gamma$. In dimensionless units its density is $r_\gamma=(H_0/M_{\rm Pl})^2 D/(3a^4)$ and its finite Maxwell pressure is $r_\gamma/3$. The term it would add to the dimensionless parent acceleration $d^2\Theta/dx^2$ is $-(H_0/M_{\rm Pl})^2gJ/a^4$. All source comparisons use this scalar acceleration. After switch-off $m_L$ can oscillate, while $g=0$ and $D$ is constant. The table's outgoing energy is consequently measured relative to a free out-vacuum with constant $I$.

## Which part can a derivative expansion reproduce?

Write $\lambda=2L$. Iterating $(\partial_\eta+i\lambda)m=s(1+2n)$ and imposing the purity identity at each order gives

$$
 n_2=\frac{s^2}{\lambda^2},\qquad \mathrm{Re}\,m_2=\frac{s'}{\lambda^2},\qquad
 n_4=\frac{3s^4+(s')^2-2s\,s''}{\lambda^4},\qquad
 \mathrm{Re}\,m_4=\frac{6s^2s'-s'''}{\lambda^4}.
$$

Each order separately obeys $n_j'=2s\,\mathrm{Re}\,m_j$, including its nonlinear terms. These local contributions integrate to zero net comoving work across this compact switch. At fixed trajectory the second-order energy sum grows quadratically with the maximum harmonic, and the fourth-order sum grows logarithmically. We evaluated these expressions on the actual rolling trajectory using derivatives of the archived potential spline on each interval. Its third derivative is piecewise defined. The calculation uses a finite spectral range and compares the exact evolution with this local approximation inside that range. It does not construct a covariant continuum stress tensor, select finite gravitational counterterms, or evolve a fourth-order gravity equation.

Residuals are the root-mean-square of the exact quantity minus the sum through the indicated order divided by the RMS of that same exact quantity. We include the powers of $a$ and the coefficient $g$ in the reported stress and source norms. The time measure is uniform in conformal time on the switching interval. A small percentage of a large transient can still exceed the outgoing energy. At the final plateau the local photon energy is zero at both orders, so the residual equals the entire outgoing excess there.

| First minimum | Harmonics | Density residual after order 2 | Density residual after order 4 | Scalar acceleration residual after order 2 | Scalar acceleration residual after order 4 | Outgoing $D$ |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| $+0.5$ | 49–192 | 13.78% | 5.76% | 11.98% | 3.41% | 91.7433 |
| $-0.5$ | 49–192 | 19.20% | 10.24% | 15.10% | 6.37% | 283.5004 |
| $+0.5$ | 49–768 | 1.530% | 0.3945% | 1.430% | 0.2334% | 91.9513 |
| $-0.5$ | 49–768 | 2.123% | 0.7294% | 1.875% | 0.4425% | 285.5551 |
| $+0.5$ | 385–768 | 0.4201% | 0.00368% | 0.4415% | 0.00470% | 0.000702 |
| $-0.5$ | 385–768 | 0.5712% | 0.00826% | 0.6236% | 0.01216% | 0.01642 |

The columns use $\delta=0.1$. The raw JSON also gives $\delta=0.5$, the shell windows 49–96, 49–384 and 129–384, sampled instantaneous quantities, integrated work, and the two spectral-window controls.

## Changing the shape of the upper frequency window

We compared each hard upper cutoff $K=96,192$ with a cosine taper among the harmonics above 48. The weight decreases from 1 at 48 as $[1+\cos(\pi(L-48)/(E-48))]/2$, reaching zero at and above $E$. We choose $E$ by solving for the same low-frequency response coefficient $W_2=\sum R_L w_L/(2L)^2$ as the sharp window. This matches both local second-order expressions at every time along the same roll. The weighting is a spectral readout or regulator on the *same* solution. It does not change the parent or postulate an additional gauge interaction.

| Sharp upper harmonic | Taper edge at matched $W_2$ | Sharp/tapered outgoing $D$, positive minimum | Sharp/tapered outgoing $D$, negative minimum |
| ---: | ---: | ---: | ---: |
| 96 | 140.731 | 84.077 / 78.586 | 242.503 / 225.396 |
| 192 | 318.549 | 91.743 / 89.592 | 283.500 / 273.408 |

With $K=96$, the taper records 6.53% and 7.05% less outgoing band energy at the same leading calibration. With $K=192$ the differences are 2.35% and 3.56%. Equality of a static local coefficient constrains a frequency-weighted sum. It leaves the rapidly driven occupations and their distribution over frequencies to be determined.

## Resolution and controls

An artificial control uniformly stretches the conformal switch duration by four while keeping the same change in $I$ as a function of $u$. The phase in the mode equations grows fourfold and the $j$th derivative of $\log I$ decreases by $4^{-j}$. This controlled drive on the conformal cylinder tests separation of time scales. It is not a second solution of the parent Einstein equations. The outgoing $D$ for the 49–768 band drops from 91.9513 to 0.000770, or from 285.5551 to 0.007685, on the respective profiles. The local fourth-order density residual in this band drops to 0.000935% and 0.00215% of the exact RMS.

The primary calculation uses 1025 switch samples and integration tolerance $2\times10^{-11}$. Doubling the samples and tightening the tolerance to $3\times10^{-12}$ changed the reported outgoing $D$ by less than $2.5\times10^{-9}$ in every case and each normalized RMS residual by less than $10^{-10}$. Both tests use the same archived geometry. The integrated source reproduces the directly measured final band energy, and the two integrated local work terms are individually checked against zero. The nominal clock reconstructed from the parent equations is compared directly with the historical phase spline. Numerical convergence for a finite set of harmonics supplies no bound on the untested infinite tail.

[Replay code](band.py) · [All numerical results](BAND.json)
