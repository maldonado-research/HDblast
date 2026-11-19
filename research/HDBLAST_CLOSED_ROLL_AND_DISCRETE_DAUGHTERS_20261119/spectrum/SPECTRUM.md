# A discrete daughter spectrum on an externally held round sphere

A finite sphere can enhance or suppress particle production from a weak homogeneous oscillation. The result depends on where its discrete frequencies fall within the finite pulse bandwidth. At 16 cycles and radius `a Htop=1`, the positive-minimum benchmark produces 2.296 times the flat continuum energy; at radius 1.1 it produces 0.3039 times that energy. The negative-minimum benchmark gives 0.9353 and 0.9984, respectively. These are vacuum-subtracted energies in one real, massless, minimally coupled daughter field, for the same imposed pulse and single-particle occupation convention.

This calculation holds a **round Einstein-frame S3 static** and prescribes the homogeneous modulus displacement externally. It does not evolve the spatial metric, the background Einstein equations, or the modulus response to the daughter. The published homogeneous radii below are used as numerical benchmarks for this static experiment. A volume radius of an inhomogeneous geometry does not determine its Laplace spectrum and is not substituted for a round radius. The [dynamic spectator calculation](../dynamic/DYNAMIC.md) separately treats the evolving geometry.

The replay produces [SPECTRUM.json](SPECTRUM.json), containing 36 cases and controls, and [SPECTRUM.csv](SPECTRUM.csv), containing each discrete level through `omega/m_sigma=4`. The [independent audit](../checks/spectral.py) uses a cutoff-free continuum expression and a separate mode formulation.

## Canonical daughter and the trace interaction

The archival four-dimensional frame relation is `g_E=f g_J`. Define `A=f^(-1/2)`, its value `A_m` at a potential minimum, and `r=A/A_m`. A constant local rescaling `chi=A_m chi_J` gives the minimally coupled daughter action on the static Einstein sphere

$$
 S_\chi=\frac12\int dt_E\,a^3d\Omega_3\,r^2
 \left[\dot\chi^2-a^{-2}(\nabla\chi)^2\right].
$$

The canonical modulus is `sigma=Mpl(Theta-Theta_m)`, with reduced Planck mass. The same trace coupling used in the archived widths is

$$
 \alpha_m=-\frac12\left.\frac{d\ln f}{d\Theta}\right|_m,
 \qquad \ln r=\alpha_m\frac{\sigma}{M_{\rm Pl}}+O(\sigma^2/M_{\rm Pl}^2).
$$

Expanding the action gives `L_int=(alpha_m/Mpl) sigma T`, where the minimally coupled scalar trace in signature `(-+++)` is `T=dot(chi)^2-a^(-2)(grad chi)^2`. There is no extra factor of two. On a harmonic displacement, the canonical oscillator has a positive mass-squared insertion `alpha_m m_sigma^2 sigma/Mpl`; equivalently its equation contains `-ddot(r)/r`. The sign follows directly from the canonical transformation rather than from a rate convention.

Let real harmonics obey `integral_S3 Y_nI Y_mJ dOmega3=delta_nm delta_IJ`. Their frequencies and multiplicities are

$$
 -\Delta_{S^3}Y_{nI}=n(n+2)Y_{nI},\qquad
 d_n=(n+1)^2,\qquad \omega_n=\frac{\sqrt{n(n+2)}}a.
$$

For `q_nI=a^(3/2) r chi_nI`, integration by parts gives

$$
 \ddot q_n+\left[\omega_n^2-\frac{\ddot r}{r}\right]q_n=0.
$$

The static minimal spectrum is `sqrt(n(n+2))/a`; replacing it by `(n+1)/a` would introduce the conformal curvature shift and change the problem.

## The pulse and its Fourier transform

The solver uses `tau=m_sigma t_E`, `mu=m_sigma/Htop`, and `T=2 pi N`. Its prescribed local coupling is

$$
 \sigma=\Phi f(\tau),\qquad
 f(\tau)=\sin^4(\pi\tau/T)\cos\tau\quad(0<\tau<T),\qquad
 r=e^{\epsilon f},\quad\epsilon=\alpha_m\Phi/M_{\rm Pl}=10^{-3}.
$$

Both exteriors have `r=1`, giving unambiguous static in and out vacua. The exponential keeps the local logarithmic coupling exact in this benchmark. It is the tangent approximation to the archival `A(Theta)`, rather than an assertion that the full archival function is exponential. At this amplitude the canonical excursions are `Phi/Mpl=0.00333415` and `0.00793621` on the positive and negative minima. Their different masses and the chosen sphere radii determine the dimensionless level spacing.

The effective exposure is

$$
 \tau_{\rm eff}=\int_0^T\sin^8(\pi\tau/T)d\tau
 =\frac{35T}{128}.
$$

With `k=omega/m_sigma`, the dimensionless driving potential is `V=epsilon f''+epsilon^2(f')^2`. The incoming mode evolves to `alpha exp(-ik tau)+beta exp(ik tau)`. Its occupation is `|beta|^2`. The first-order coefficient is

$$
 \beta_1(k)=\frac1{2ik}\int_0^T\epsilon f''e^{-2ik\tau}d\tau
 =2i\epsilon k\widetilde f(2k).
$$

The code evaluates the exact compact Fourier transform as a finite sum of shifted sinc functions:

$$
 \widetilde f(\nu)=\frac12\sum_{j=0}^4c_j\sum_{s=\pm1}
 D\left(\nu-s-\frac{j-2}{N}\right),\quad
 c_j=\frac{(-1)^{j-2}}{16}\binom4j,
$$
$$
 D(z)=T e^{-izT/2}\operatorname{sinc}(zT/2),\quad
 \operatorname{sinc}x=\frac{\sin x}{x}.
$$

The spectral line centers on the pair condition `2 omega=m_sigma` and narrows with duration. At that condition, `N>=3` integer cycles give `f_tilde(1)=3T/16` and

$$
 |\beta_1(1/2)|^2=(3\pi\epsilon N/8)^2.
$$

The `sin^4` window has three continuous derivatives at the support endpoints, making the outgoing particle energy converge in the ultraviolet. Abruptly truncating `sin^2(pi tau/T) cos tau` instead leaves jumps in `f''`; it generally gives `|beta|^2` proportional to `k^-4` at high frequency and logarithmically divergent particle energy. A `sin^2` window with an integer-cycle sine carrier is a different, smoother-endpoint protocol.

## Particle count, energy, and continuum reference

The sphere volume is `V_S3=2 pi^2 a^3`. The produced particle and energy densities are

$$
 n_{S^3}=\frac1{2\pi^2a^3}\sum_{n\ge1}(n+1)^2|\beta_n|^2,
 \qquad
 \rho_{S^3}=\frac1{2\pi^2a^3}\sum_{n\ge1}(n+1)^2\omega_n|\beta_n|^2.
$$

The occupation counts individual daughters. The number of pairs is half the particle count. The comparison uses the same oscillator occupation in the flat continuum:

$$
 n_\infty=\frac1{2\pi^2}\int_0^\infty k^2|\beta_k|^2dk,
 \qquad
 \rho_\infty=\frac1{2\pi^2}\int_0^\infty k^3|\beta_k|^2dk.
$$

In these integrals `k` denotes a dimensional frequency; the numerical output divides it by `m_sigma`. JSON densities are provided in both `m_sigma` and `Htop` units. “Yield” here is particle number per physical volume, with no entropy or temperature assigned. Vacuum and Casimir energies are excluded.

The zero mode has `omega_0=0`, so a massless minimally coupled scalar has no preferred static oscillator vacuum for that mode. It is excluded from every particle sum. In the original variable its canonical momentum is constant, and its endpoint transfer is `(chi,P) -> (chi+P integral r^(-2)dt/a^3,P)`. A zero-mode quantum state must be supplied separately. A constant initial field with zero momentum returns unchanged.

## Calibrated sphere benchmarks

The archival masses and trace couplings are read from `gravity/previous/prior/archive/BLAST_REHEATING_RESULTS.json` and checked against `modulus/widths/WIDTHS.json`. The round homogeneous radii come from `modulus/background/BACKGROUND.json`, with the width off. They do not describe the static pulse experiment's time evolution.

| Minimum | `mu=m_sigma/Htop` | `alpha_m` | Published first minimum crossing radius | Published first maximum radius |
|---|---:|---:|---:|---:|
| `phi_m=+0.5` | 7.683603799 | 0.2999263092 | 1.110835249 | 1.208508434 |
| `phi_m=-0.5` | 16.54253876 | 0.1260047724 | 1.136351847 | 1.301152792 |

All radii below are in `Htop^-1`. These tables give the **full oscillator** result divided by the full flat continuum result for the identical pulse. The Born ratios differ only slightly; the largest full/Born energy difference among the 36 cases is 0.0195%.

| Minimum | Radius | Role | Energy ratio, 4 cycles | 8 cycles | 16 cycles | Particle ratio, 16 cycles |
|---|---:|---|---:|---:|---:|---:|
| `+0.5` | 1 | Unit | 1.026873 | 1.129393 | 2.295960 | 2.297666 |
| `+0.5` | 1.1 | Representative | 1.022266 | 0.941055 | 0.303897 | 0.333287 |
| `+0.5` | 1.110835 | Published first crossing | 1.021839 | 0.956564 | 0.214379 | 0.234269 |
| `+0.5` | 1.208508 | Published first maximum | 1.018486 | 1.062311 | 1.260857 | 1.205631 |
| `-0.5` | 1 | Unit | 1.005865 | 1.006858 | 0.935307 | 0.940917 |
| `-0.5` | 1.1 | Representative | 1.004850 | 1.005671 | 0.998411 | 1.004003 |
| `-0.5` | 1.136352 | Published first crossing | 1.004546 | 1.005315 | 0.985128 | 0.983755 |
| `-0.5` | 1.301153 | Published first maximum | 1.003469 | 1.004057 | 1.008441 | 1.008121 |

The positive branch at radius 1 is close to a level, with fractional pair detuning `2 omega/m_sigma-1=+0.00811636`. At radius 1.1 the closest level has detuning `-0.0835306`. Neither published radius was selected to be exactly resonant: the positive branch's first-crossing and maximum radii give detunings `-0.09247` and `+0.05517`; the negative branch gives `-0.04839` and `+0.01786`.

Two additional radii are constructed independently from the spectrum. For a chosen integer level `n`, exact resonance is `a_res=2 sqrt(n(n+2))/mu`. The midpoint of neighboring frequencies lies at `a_mid=[sqrt(n(n+2))+sqrt((n+1)(n+3))]/mu`.

| Minimum | Level `n` | Exact resonant radius | 16-cycle energy ratio | Frequency-midpoint radius | 16-cycle energy ratio |
|---|---:|---:|---:|---:|---:|
| `+0.5` | 3 | 1.00811636 | 2.216600 | 1.14164695 | 0.201640 |
| `-0.5` | 7 | 0.95961739 | 1.088365 | 1.02049184 | 0.958273 |

The same finite pulse can therefore yield enhancement or suppression. The more widely spaced positive-branch spectrum resolves these 16-cycle lines more strongly. A weak source does not by itself justify replacing a finite spectrum by a continuum.

As the radius increases, the Born energy ratios approach unity. At `a=10,20,40` the positive branch gives `1.00033351,1.00008339,1.00002085`; the negative branch gives `1.00007196,1.00001799,1.00000450`. A smooth curved-space state count itself differs slightly from the flat count: writing `j=n+1` gives `j^2 dj/r^3=k^2 sqrt(1+1/(rk)^2) dk`, where `r=mu a`. This explains the positive correction proportional to inverse radius squared after the individual lines become dense.

## Width and independent controls

The same interaction has the flat-space identical-real-scalar width

$$
 \Gamma=\frac{\alpha_m^2m_\sigma^3}{32\pi M_{\rm Pl}^2}.
$$

For `rho_osc=m_sigma^2 Phi^2/2`, its narrow-line pulse prediction is

$$
 \frac{\rho_\Gamma}{m_\sigma^4}
 =\frac{\epsilon^2\tau_{\rm eff}}{64\pi}.
$$

The flat finite-window Born energy divided by this prediction is `1.37500014,1.09040179,1.02239118` at 4, 8, and 16 cycles. This duration effect exists before introducing a discrete sphere spectrum; the reported sphere/continuum ratios keep both pulse shapes identical. The model includes neither damping nor Bose enhancement from a populated incoming state.

The replay checks the following independent statements:

- Zero drive produces exactly zero `beta`. Direct second-order oscillator evolution agrees with the Bogoliubov solver to `4.4e-15` in the complex coefficient. The maximum `|alpha|^2-|beta|^2-1` residual is `1.1e-14`.
- Reciprocal pulses `r` and `1/r` have identical occupations. Factorization of `-d^2/dtau^2+r''/r` maps the two scattering problems and gives equal outgoing `alpha` and opposite `beta`; the largest numerical complex-coefficient discrepancy is `4.9e-15`.
- Doubling `epsilon` increases the representative 8-cycle energies by `4.000414` and `4.000398`. Halving it gives `0.2499935` and `0.2499938`. The maximum occupation over all primary cases stays below `3.6e-4`.
- Direct time quadrature of the pulse agrees with its analytic Fourier transform to `5.1e-17` in `beta_1`. The exact resonant Born occupation matches `(3 pi epsilon N/8)^2`.
- The continuum Born particle count obeys the independent Parseval identity `n_1/m_sigma^3=integral(epsilon f'')^2 d tau/(16 pi)`. Its largest relative discrepancy is `4.5e-11`.
- Increasing the continuum quadrature order from 8 to 12 produces no resolved change at 16 cycles. Increasing the exact frequency cutoff from 4 to 8 changes the 4-cycle continuum energy by `3.7e-10` relative; the two representative sphere energies change by at most `3.7e-10`.
- Analytic Born energy tail bounds above `omega/m_sigma=4` are at most `1.12e-16 m_sigma^4` for the representative sphere cases and `1.03e-16 m_sigma^4` for the continuum. The bounds are for the linear calculation; the full calculation has the separately measured cutoff test.

From the capsule root:

```sh
OPENBLAS_NUM_THREADS=1 python -IB spectrum/spectrum.py --gravity gravity --output /tmp/hdblast-spectrum
```

All numerical inputs are pinned local capsule files. The pulse is prescribed for a finite duration with exact static exteriors; using this result as a local rate in an evolving condensate requires an additional dynamical approximation.
