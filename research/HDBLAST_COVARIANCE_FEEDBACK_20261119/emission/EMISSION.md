# Weak-probe fermion pairs following a change in the parent covariance

We evaluate a switched scalar–Dirac Yukawa interaction on the returned five-dimensional shell trajectories. The scalar covariance and the geometry already influence each other during each native run. The daughter readout is the coefficient at order $y^2$ as $y$ tends to zero. Feeding the emitted pairs back would change that perturbative order and requires a stress and energy balance for the interaction and its switch.

## Preparation and units

At the common pre-pulse fork the positive-frequency parent solution is fixed by $u_i=(2a_i^3E_i)^{-1/2}$ and $a_i^3\dot u_i=-i a_i^3 E_i u_i$. The instantaneous frequency is $E_i^2=(q/a_i)^2+100^2(\phi_i-1/2)^2$ in units of $H_0$. Write $I=a_i^3E_i$, $X=\langle\chi^2\rangle$, $P=\langle\Pi^2\rangle$ and $C=\langle\{\chi,\Pi\}/2\rangle$. The incoming occupation and the two-mode correlation of opposite momentum are
$$
n=\tfrac12(IX+P/I-1),\qquad c=\tfrac12(IX-P/I)+iC .
$$
We use the actual native Gauss–Legendre momentum nodes and weights $q^2 dq/(2\pi^2)$, covering $0<q<8$. The uniform canonical momentum shift changes the occupation spectrum. Each state uses its own recorded covariance. The massless Dirac daughters start in the conformal vacuum. The two proper-time windows are $1.35\le H_0\tau\le2.05$ and $1.45\le H_0\tau\le1.95$. Each uses a $\sin^4$ switch compactly supported on its interval, which starts after the fork at $H_0\tau\simeq1.30372$.

For the conformally rescaled parent mode $v=a u$, define $F_\pm(\omega)=\int d\eta\,f v e^{\pm i\omega\eta}$. For positive total pair energy $\omega\ge q$ the inclusive differential pair count per squared coupling is
$$
\frac{q^2dq}{2\pi^2}\,\frac{(\omega^2-q^2)d\omega}{8\pi^2}
\left[ |F_-|^2+n(|F_+|^2+|F_-|^2)+2\operatorname{Re}(c F_+F_-)\right].
$$
The first term is emission on the same wall from the instantaneous parent vacuum at the fork. The last two terms give the excess above this reference. A vacuum evolved from the earlier field preparation can have nonzero occupation and phase at the fork. We check this distinction quantitatively by continuing the original vacuum mode through the uninterrupted retained-state history and recovering its count. Energy weights the integrand by $\omega/a_f$ at $H_0\tau_f=2.10$. Number and energy are per unit comoving volume, in $y^2 H_0^3$ and $y^2 H_0^4$ respectively with the native scale-factor convention. Dividing by $a_f^3$ gives physical densities.

## What the comparison measures

For each preparation use its covariance to read both the returned altered-state geometry and the returned retained-state geometry. Repeat each experiment for a quiet shell and three incident pulses. The difference between those two pulse-minus-quiet responses isolates the effect of the responding wall on the weak-probe readout. The JSON field named `geometry_effect` includes changes in both expansion and the scalar-dependent parent mass. For a given preparation, units in the reported contrast tables are percentage points of the quiet, before-subtraction count or energy on the retained-state geometry. This shared positive denominator also applies to the matched-vacuum excess. The geometry of the reference experiment evolves under the native equations throughout.

## Returned-wall results

| Fine quiet shell, wide switch and full resolved spectrum | Total pair count on own wall | Change from the retained preparation on its own wall |
|---|---:|---:|
| Retained correlation | 0.76981630 | — |
| Erased correlation | 0.77525214 | +0.70612% |
| Reflected correlation | 0.78068680 | +1.41209% |
| Uniform canonical impulse | 0.77513041 | +0.69031% |

All four positive total counts are below their respective switched fork-vacuum references on the quiet fine wall. The reflected state gives $0.78070653$ on the retained wall. The switch from that wall to its own wall changes its count by −0.002528%. For each preparation, nearby positive, nearby negative and distant-compensated positive fronts change the total count in the directions decrease, increase and increase. The vacuum-subtracted excess changes in the opposite directions. Both patterns hold across both mesh resolutions, both switching windows and both total-energy selections. Energy and physical total densities were also evaluated.

The fine, wide reflected-state total-count contrasts on their own walls are −0.169868%, +0.989009% and +0.093876%, normalized as above. On retained-state walls the same preparation gives −0.169909%, +0.989114% and +0.093856%. Across both grids, all preparations, windows and energy selections, the largest change in this comoving *number* contrast due to wall response is 0.000112 percentage point. For the corresponding vacuum-subtracted number contrasts the wall-dependent changes can be smaller than the numerical comparison uncertainty. [Numerical results](EMISSION.json)

### A combination that vanishes on a single wall

At order $y^2$ the readout is linear in the covariance. Since the erased covariance is exactly halfway between the retained and reflected ones, using one background for all of them gives $N_{\rm retain}-2N_{\rm erase}+N_{\rm reflect}=0$ for each pulse. On their separately evolved quiet walls the fine wide-window result divided by the retained count is $-1.52927\times10^{-6}$, or −1.52927 parts per million. The coarse value is −1.55706 parts per million. The independent methods differ by less than 0.018 parts per million across all the reported midpoint observables. The two-grid agreement is an empirical check. The effect includes the changing mass and scale factor of the responding wall. Other nonlinear dynamics can also violate this affine relation.

### Separate contributions from the mass and scale-factor histories

For a diagnostic on the reflected quiet state, combine either scale-factor history with either mass history on the same proper-time clock. For each input, average the effect of changing it first and changing it second. This symmetric attribution exactly adds up to the difference between the two computed physical walls. The two combinations using a scale factor and mass from different native runs serve as readout controls. On the fine wide, full-spectrum total pair count, the responding wall gives −25.2765 parts per million relative to the retained wall. The scale factor contributes +0.9432 parts per million and the mass history contributes −26.2197 parts per million. The corresponding coarse terms are +0.9582 and −26.8409 parts per million. This comparison identifies the dominant change for this particular readout and state. [Four combinations and independent integration](BACKGROUND_SPLIT.json)

### Preparations with matched initial source derivatives

This follow-up uses the two additional quiet native evolutions with opposite nonzero correlations that also match the initial first pressure and force derivatives of the erased state. Using its quiet total as the common denominator, the difference between the positive and negative correlation profiles is −0.145662% in the fine wide-window full count on the responding walls. On the common erased-state history it is −0.145838%. The wall response thus shifts the difference by +0.000176 percentage point. The coarse shift is +0.000177 percentage point. Matching the initial sources and their first derivatives still allows the later emitted counts to differ. The response difference is mostly due to the different parent state entering the readout. [Supplemental quiet emission and both quadratures](SMOOTH.json)

## Reproduction

The calculation authenticates the published source-matching helper and its two moving-parent dependencies before importing them. The primary propagation uses proper time, DOP853, a complex Fourier transform shifted to the light-cone threshold, and a separate 96-node frequency integral for $\omega\le60$. The independent propagation uses a symplectic transfer matrix with 16,384 conformal steps and direct Gaussian frequency quadrature through $\omega=768$. The full Fourier calculation has a separate, much higher Nyquist frequency. All emitted-pair counts are restricted to total comoving momentum at most 8. The input wall and fork data reside in the adjacent native directory. From this directory run

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -I -B check.py /path/to/HDblast ../native --full
~~~
