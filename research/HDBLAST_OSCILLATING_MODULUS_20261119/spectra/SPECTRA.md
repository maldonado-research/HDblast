# How many quanta carry the radiation?

We follow one real, massless, minimally coupled scalar produced in the two-body decays of a spatially homogeneous, oscillation-averaged parent of mass $m$. The parent is nonrelativistic, its rest width is constant, and neither daughters nor parents undergo additional interactions. The flat expanding universe initially contains only the parent. This control specifies $H_i=m/50$ and $a_i=1$. A second calculation at detuning $10^{-3}$ begins at the late minimum crossing of each resolved **undamped** flat history. Switching on decay at that crossing isolates late emission. The curvature study separately carries the small founder radiation density of the harmonic-gated histories through the fluid equations. No initial momentum distribution is assigned to that fluid.

## Momentum as a record of decay time

In each rest-frame decay two daughters emerge with physical momentum $m/2$. At scale factor $a$ they enter the comoving momentum bin $k=am/2$. The expansion is monotone in this test. Each bin therefore receives daughters at a single time. With $N$ the number per comoving volume and $Q=a^4\rho_r$, the source identities are

$$
\dot N=2\Gamma a^3\rho_m/m,\qquad \dot Q=\Gamma a^4\rho_m,\qquad
\frac{dN}{d\ln k}=\frac{2\Gamma a^3\rho_m}{mH}.
$$

The third identity divided by the phase-space measure $k^3/(2\pi^2)$ gives the occupation of each newly crossed bin,

$$
f(k)=\frac{32\pi^2\Gamma\rho_m}{m^4H}
=3\pi\alpha^2\,\Omega_m\frac{H}{m},\qquad
\Omega_m=\frac{\rho_m}{3M_{\rm Pl}^2H^2}.
$$

Here we used the width for this particular daughter, $\Gamma=\alpha^2m^3/(32\pi M_{\rm Pl}^2)$, with the reduced Planck mass. Subsequent collisionless motion preserves $f$. Since $H$ falls on these expanding flat histories and $\Omega_m\le1$, $f$ stays below 0.01696 and 0.002994 for the two $m/50$ launches. These low occupations support the dilute approximation once oscillations are established. Coherent amplification in the earlier anharmonic roll requires a different calculation. The injection formula with $1/H$ cannot be applied through a turnaround.

For a stable integration over many decades let $u=\Gamma(t-t_i)$, $\gamma=\Gamma/H_i$, $b=\gamma^{2/3}a$ and $P=\int_0^u b(v)e^{-v}\,dv$. The expanding Friedmann solution satisfies $du/d\ln a=b^2/\sqrt{be^{-u}+P}$. The comoving number below the bin born at $u$ is proportional to $1-e^{-u}$, and comoving radiation energy below it is proportional to $P$. This gives the full momentum distributions at any observation time. We independently integrated $p^3 f/(2\pi^2)$ and $p^4 f/(2\pi^2)$ over $d\ln p$ to recover the analytic number density and the Friedmann radiation density. An independent RK45 run checks the DOP853 milestones.

## Fixed-energy Bose comparison

Define $T_*=(30\rho_r/\pi^2)^{1/4}$ solely to specify the equilibrium distribution having this energy. Its number density for one real scalar at zero chemical potential is $n_* = \zeta(3)T_*^3/\pi^2$. At equality $\rho_r=\rho_m$. At the second event $\rho_r=9\rho_m$. The entries give $100 n/n_*$.

| Detuning | Equality, $+0.5$ | Equality, $-0.5$ | 90% radiation, $+0.5$ | 90% radiation, $-0.5$ |
|---|---:|---:|---:|---:|
| $10^{-3}$ | 8.744% | 5.390% | 7.180% | 4.426% |
| $10^{-10}$ | 0.1555% | 0.09587% | 0.1277% | 0.07871% |
| $10^{-20}$ | 0.000492% | 0.000303% | 0.000404% | 0.000249% |

The two delayed-start checks at $10^{-3}$ give 8.724% and 5.368% at equality. At $H_i=m/50$ the mean momentum at equality is about $0.268m$ and the median by energy is $0.348m$ on both branches. Radiation can carry most of the energy while the gas contains relatively few energetic quanta. At fixed event the number fraction is proportional to $\sqrt{M_{\rm Pl}\Gamma}/m$ times a dimensionless history factor. In the small-width limit this factor approaches a constant. Here $m\propto\sqrt{\text{detuning}}$ and $\Gamma\propto m^3$, giving the observed fourth-root dependence on detuning.

There is also an all-time ceiling in this production model. Among isotropic distributions with fixed energy and $0\le f\le f_{\max}$, filling the lowest momenta to the cap maximizes number. Replacing any occupied higher-momentum element by lower-momentum quanta can carry the same energy with more particles. For the maximizing step distribution of radius $p_F$, $n=f_{\max}p_F^3/(6\pi^2)$ and $\rho=f_{\max}p_F^4/(8\pi^2)$. Eliminating $p_F$ gives

$$
\frac n{n_*}\le \frac{(8\pi^4/30)^{3/4}}{6\zeta(3)}\, f_{\max}^{1/4}.
$$

Even at arbitrarily late times this gives ceilings of 57.57% and 37.32% for the $H_i=m/50$ controls. The cap follows from this source and the specified empty initial state. Bose statistics alone impose no upper occupation. Additional species, a populated initial bath, broad resonant production or reactions that change particle number require a different calculation.

[Epochs, quantiles and numerical checks](SPECTRA.json) · [Momentum distributions](SPECTRA.csv) · [Calculation](spectra.py)
