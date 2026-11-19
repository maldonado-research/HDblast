# Counting particles from an oscillating Jordan metric

A massless real scalar can absorb energy from an oscillating scale factor when it is minimally coupled to curvature. We counted the outgoing particles mode by mode and compared their total energy with the decay rate used in the cosmology calculation. For a pulse lasting 32 oscillations, the measured energy was **0.610% above the narrow-line rate prediction**. Accounting for the finite duration of the pulse reduced the difference to **0.0511%**. A massless scalar with conformal curvature coupling gave exactly zero production in the same geometry. [Study protocol](https://github.com/maldonado-research/HDblast/issues/157#issuecomment-7018476842) · [Mode results](MODES.json) · [Occupation spectra](LINES.csv)

## A pulse with an unambiguous vacuum

Take the Einstein-frame scale factor $a_E=1$, the canonical modulus displacement $\sigma=M_P(\Theta-\Theta_m)$, and $t=m_\sigma\eta$. We set $m_\sigma=1$ in numerical values. The imposed displacement and corresponding local approximation to the Jordan metric are

$$
 \sigma(t)=\Phi f(t),\qquad f(t)=w(t)\cos(t),\qquad
 w(t)=\sin^p(\pi t/T),\quad T=2\pi N,\qquad
 a_J(t)=\exp[\epsilon f(t)],\qquad \epsilon=\alpha_m\Phi/M_P .
$$

Here $\alpha_m=(d\ln A/d\Theta)_m$. The pulse is supported on $0<t<T$. The scale factor is one on both exterior intervals, where the positive-frequency vacuum is exact. The primary window uses $p=4$. We also tested $p=6$. For the actual archived coupling, we evolve the same sinusoidal field displacement with $a_J=A(\Theta)/A(\Theta_m)=\sqrt{f_m/f(\Theta)}$, using the archive's positive Planck factor in this last expression. In that expression alone $f(\Theta)$ denotes the archived background function. The switching window and stationary Einstein volume make a controlled test of the local particle-production rate. On the cosmological trajectories the modulus amplitude and $a_E$ evolve under their own equations. [Background dynamics](../background/) · [Width derivation](../widths/)

The canonical variable $u=a_J\chi$ for one massless real scalar with curvature coupling $\xi$ obeys

$$
 u_k''+[k^2-V(t)]u_k=0,\qquad V(t)=(1-6\xi)\,a_J''/a_J .
$$

Primes here are derivatives in units of $m_\sigma$. Evolving a normalized incoming $e^{-ikt}$ mode gives an outgoing $\alpha_k e^{-ikt}+\beta_k e^{ikt}$. The momentum occupation is $|\beta_k|^2$. The full-space radial integrals $n=\int k^2|\beta_k|^2dk/(2\pi^2)$ and $\rho=\int k^3|\beta_k|^2dk/(2\pi^2)$ count single particles in units of $m_\sigma^3$ and their energy in units of $m_\sigma^4$. The pair count is $n/2$. The direct oscillator equation and a first-order evolution for $\alpha,\beta$ give the same outgoing coefficients.

The local perturbation is $V_1=(1-6\xi)\epsilon f''$. First-order scattering gives

$$
 \beta_{k,1}=\frac{1}{2ik}\int_0^T V_1(t)e^{-2ikt}\,dt
            =2ik(1-6\xi)\epsilon\int_0^T f(t)e^{-2ikt}\,dt .
$$

The code evaluates the Fourier transform analytically as a finite sum of shifted sinc functions. An independent particle-count integral follows from Parseval's identity, using the negative-frequency conjugate for a real drive.

$$
 n_1=\frac{1}{16\pi}\int_0^T V_1(t)^2\,dt.
$$

The temporal integral is evaluated by orthogonality of the finite Fourier series. It agrees with the independent frequency integration through $k=4m_\sigma$ to within $1.1\times10^{-13}$ relative in the primary tangent cases.

## Comparison with the decay width

In flat space the linear scalar vertex evaluated on the mass shell has magnitude $g=(1-6\xi)\alpha_m m_\sigma^2/M_P$. Two identical real, massless daughters give $\Gamma=g^2/(32\pi m_\sigma)$. A harmonic classical condensate with peak amplitude $\Phi$ has cycle-averaged energy density $\rho_{\rm osc}=m_\sigma^2\Phi^2/2$. The long weak-pulse prediction is therefore

$$
 \rho_{\rm out}\simeq\Gamma\rho_{\rm osc}\int w^2\,d\eta .
$$

The normalization uses the exposure $\int w^2d\eta=2\pi N\binom{2p}{p}/(4^p m_\sigma)$. For finite duration, derivatives of the envelope and the finite width of the spectral line can change the outgoing energy. The explicit comparison below uses $\epsilon=0.001$, $p=4$ and $\xi=0$.

| Cycles | Exact outgoing energy in $m_\sigma^4$ | Exact / narrow-line estimate | Finite-window Born / narrow-line estimate | Mean daughter momentum in $m_\sigma$ | Energy-weighted RMS offset from $m_\sigma/2$ |
|---:|---:|---:|---:|---:|---:|
| 8 | $7.45417\times10^{-8}$ | 1.090438 | 1.090402 | 0.517369 | 0.050684 |
| 16 | $1.39798\times10^{-7}$ | 1.022523 | 1.022391 | 0.504433 | 0.024085 |
| 32 | $2.75105\times10^{-7}$ | 1.006098 | 1.005585 | 0.501114 | 0.011871 |

The line narrows about $k=m_\sigma/2$ as the duration increases. The largest occupancy at 32 cycles is about 0.00142. Doubling $|\epsilon|$ at 16 cycles increases the integrated energy by a factor of 4.002. At $\xi=1/12$ the energy is about one quarter of its minimally coupled value, as expected from the squared improvement factor. With $\epsilon=0$ or $\xi=1/6$, $V$ vanishes and the mode equations preserve the vacuum exactly. Long strong oscillations may amplify previously produced quanta; the small-amplitude, small-occupation tests here establish the perturbative regime.

On the two archival branches $\phi_m=+0.5$ and $-0.5$, the spline gives $\alpha_m=0.2999504$ and $0.1259528$ respectively. At $\epsilon=0.001$, the corresponding field amplitudes are $0.003334 M_P$ and $0.007939 M_P$. Repeating the 16-cycle test with the full tabulated matter coupling yields outgoing energies **0.000038% below** and **0.0932% above** the tangent metric result. These runs retain the nonlinear field dependence of $A$. Their modulus paths are prescribed and their Einstein-frame scale factors are stationary.

## A sign check from reciprocal metrics

For a minimally coupled massless scalar the occupation is exactly the same for $a_J$ and $1/a_J$ whenever both have constant exteriors and the same conformal clock. To see this, set $W=(\ln a_J)'$ and $D=\partial_t-W$. The two Schrödinger operators factor as

$$
 L_a=D^\dagger D=-\partial_t^2+W^2+W',\qquad
 L_{1/a}=D D^\dagger=-\partial_t^2+W^2-W'.
$$

Applying $D/(-ik)$ maps a mode of $L_a$ to a mode of $L_{1/a}$ with the same incoming normalization. Where $W=0$, the transformation preserves the coefficient of $e^{-ikt}$ and reverses the sign of the coefficient of $e^{ikt}$. Thus the two metrics have identical $|\beta_k|^2$ at every $k>0$. The exponential tangent pulses with opposite $\epsilon$ are precisely reciprocal. The computed $\alpha$ differences and $\beta$ sums both stayed below $1.1\times10^{-14}$ over the 16-cycle momentum grid. This factorization supplies an exact check of the numerical mode phases.

## Reproduction checks

The energy was integrated to $k/m_\sigma=1,2,4$. For the 8-, 16- and 32-cycle minimal-scalar tests, the largest change from 1 to 4 is $1.7\times10^{-8}$ of the integrated energy. The first-order solver preserves $|\alpha|^2-|\beta|^2=1$ to $1.3\times10^{-14}$ in the primary tangent cases. An independent time-domain solution at six selected momenta gives complex $\beta$ values within $3.1\times10^{-15}$ of the first-order solution for the reference pulse. The JSON reports the tighter ODE and quadrature run, both envelope controls, sign reversal, the actual-coupling comparisons, particle counts and all cutoff readouts.

From the study root, run the following command.

~~~sh
OPENBLAS_NUM_THREADS=1 python -IB modes/modes.py --archive archive --output /tmp/hdblast-modulus-modes
~~~
