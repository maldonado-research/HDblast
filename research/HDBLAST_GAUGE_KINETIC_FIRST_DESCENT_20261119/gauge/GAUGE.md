# Photons from a specified gauge coupling on the two rolling spheres

Free Maxwell waves on a homogeneous expanding universe obey the same conformal-time oscillator equation as waves on a fixed cylinder. A change in the electromagnetic kinetic coefficient supplies a way to excite those oscillators. This study specifies one such change and computes the excess over the free outgoing vacuum. For a log change of 0.1, the photon energy density at the first potential minimum is **0.773%** or **1.228%** of the original parent-field density, depending on the rolling branch. The size and duration of the coupling change are assumptions of the calculation. The earlier parent-field fit determines the trajectory once the new interaction profile and its parameters have been chosen.

## Pulse and initial state

On the archived round closed backgrounds, let $q=\Theta-\Theta_{\rm top}$, $x=H_0t$, $a=H_0a_{\rm phys}$, and $d\eta=dx/a$. We integrate the unseeded parent equations from the turning point at the bounce. The code checks the parent position, velocity, scale factor and expansion rate at the previously archived first minimum. Define the directed field fraction $s=(q-q_0)/(q_{\min}-q_0)$, which grows from 0 to 1 on either first descent. The tested action and pulse are

$$
S_\gamma=-\frac14\int d^4x\,\sqrt{-g}\, I^2(\Theta)F_{\mu\nu}F^{\mu\nu},
\qquad \ln I=\delta S(u),\qquad u=(s-0.2)/0.6.
$$

The dimensionless function $S$ rises from 0 to 1. Its derivative on $0<u<1$ is $b(u)=\exp[-1/(u(1-u))]/\int_0^1\exp[-1/(v(1-v))]\,dv$ and vanishes elsewhere. Thus $I=1$ at launch and $I=e^\delta$ after switch-off. This is a smooth, positive, explicitly added coupling. A universal conformal metric factor drops out of the four-dimensional Maxwell action and cannot supply this term. We launch the conformal photon vacuum on the pre-switch plateau of the expanding solution. There is an unambiguous free in-basis on that plateau and a free out-basis from switch-off through first passage. The ancestor geometry can also be continued to earlier contracting times. A state originating earlier on that continuation would need its own initial conditions and history of $I$.

## Mode propagation and physical normalization

For an orthonormal transverse harmonic on the unit three-sphere let $L=2,3,\ldots$ be its conformal frequency. Each $L$ has $2(L^2-1)$ real modes including both polarizations. One oscillator has action $\tfrac12\int I^2(Q_\eta^2-L^2Q^2)d\eta$ and momentum $P=I^2Q_\eta$. Define $b_L=(I\sqrt L\,Q+iP/(I\sqrt L))/\sqrt 2$ and remove the free phase by writing $b_L(\eta)=e^{-iL(\eta-\eta_0)}(\alpha_L b_{L,0}+\beta_L b_{L,0}^{\dagger})$. In the switch coordinate the evolution is

$$
\frac{d\alpha_L}{du}=\delta b(u)e^{2iL(\eta-\eta_0)}\beta_L^*,\qquad
\frac{d\beta_L}{du}=\delta b(u)e^{2iL(\eta-\eta_0)}\alpha_L^*,
\qquad (\alpha_L,\beta_L)_{\rm in}=(1,0).
$$

Here $\eta_0$ is conformal time at switch-on. Changing $\delta$ to $-\delta$ leaves $\alpha_L$ fixed and flips $\beta_L$ for these source-free linear equations. The spectrum is therefore identical for reciprocal $I$ profiles. An independent evolution directly in $Q,P$ checks selected shells. The first-order (Born) amplitude follows by setting $\alpha_L=1$ in the second equation. After the switch the occupation $N_L=|\beta_L|^2$ is time independent.

The photon density in excess of the free out-vacuum uses the physical volume $2\pi^2(a/H_0)^3$ and photon energy $H_0L/a$. Consequently

$$
\rho_\gamma=\frac{H_0^4}{a^4}D,\qquad
D=\frac1{\pi^2}\sum_{L=2}^\infty L(L^2-1)N_L,\qquad
\frac{\rho_\gamma}{\rho_\phi}=\frac{(H_0/M_{\rm Pl})^2D}{3a^4[U+q_x^2/6]} .
$$

The reduced Planck mass is used and $H_0/M_{\rm Pl}=0.0088134414$ is the archived calibration. $U$ is the fitted dimensionless parent potential. The free out-vacuum subtraction compares states of the same final Hamiltonian and geometry. It removes the common zero-point term. A local absolute stress tensor during the pulse would require additional renormalization. After switch-off the isotropic excess has radiation pressure $p_\gamma=\rho_\gamma/3$ and redshifts as $a^{-4}$. Occupancy of free photon modes says nothing by itself about thermal equilibrium or about the origin of $I$ in a higher-dimensional theory.

## Scan at the first minimum

| $\lvert\delta\rvert$ | $D$, $\phi_m=+0.5$ | $\rho_\gamma/\rho_\phi$, $+0.5$ | $D$, $\phi_m=-0.5$ | $\rho_\gamma/\rho_\phi$, $-0.5$ |
|---:|---:|---:|---:|---:|
| 0 | 0 | 0 | 0 | 0 |
| 0.01 | 3.802 | 0.00007699 | 6.745 | 0.00012236 |
| 0.10 | 381.697 | 0.007730 | 676.998 | 0.012281 |
| 0.50 | 10,510.5 | 0.21285 | 18,502.7 | 0.33564 |
| 1.00 | 56,548.9 | 1.14521 | 97,554.8 | 1.76966 |

Ratios at switch-off for the primary $\delta=0.1$ case are 0.007632 and 0.012113 in branch order. The ratio can rise slightly after switch-off even though the photon density falls because the parent density also changes. The Born integrated energy underestimates the exact primary by 0.4034% and 0.3711%. At $|\delta|=1$ those errors grow to 32.8% and 30.9%. The zero-coupling controls give exactly zero in every shell. The direct negative-amplitude test recovers $\beta(-0.1)=-\beta(+0.1)$.

The unit-calibrated conformal switch intervals are $[0.3361853,0.4654661]$ and $[0.3748579,0.4938455]$. In the primary scan the most energetic single shells are $L=35$ and 40. Half the cumulative energy lies at or below $L=38$ and 45. The 90% thresholds are $L=61$ and 84, and the 99% thresholds are $L=115$ and 152. At first passage those 99% shells have physical momenta $0.91$ and $1.18\,M_{\rm Pl}$. Any physical ultraviolet cutoff chosen for the effective interaction can be applied directly to the supplied cumulative spectrum. The fitted potential and the chosen gauge interaction have no independently established domain extending to all those momenta.

With the original background total energy budget held fixed, assigning half of it to radiation at first minimum would require $D=24{,}689$ or $27{,}563$. Requiring equality with a remnant whose density is the entire original parent density instead gives $D=49{,}379$ or $55{,}126$. For the first definition to be realized dynamically, the parent, geometry and photons would have to exchange energy consistently. The primary ratios are small compared with either benchmark. The large-$|\delta|$ entries show where backreaction becomes essential. In particular, the $|\delta|=1$ prescribed backgrounds still contain their full original parent energy even as the computed photon energy exceeds it. Those numbers cannot represent a self-consistent radiation-dominated solution. Small final photon fractions are a necessary test of the spectator approximation. They do not bound a transient renormalized stress or scalar force during the switch.

## Duration sensitivity

A supplementary test keeps the same endpoints of $I$, the same $|\delta|=0.1$, the same midpoint in field space and the same family of smooth profiles. Extending the support to 15%–85% of the field displacement lowers the crossing fractions to **0.426%** and **0.699%**. Shortening the support to 25%–75% raises them to **1.573%** and **2.426%**. The narrower pulse also moves the median energy shells from 38 and 45 to 45 and 52. Thus changing just the assumed switch width can approximately double or halve the primary yield. An independent derivation or constraint on the gauge profile would substantially change how predictive the calculation is.

## Numerical checks

The primary product sums complete shells $2\le L\le1536$ and saves one occupation, Born occupation and cumulative energy per shell. Index zero in those arrays is $L=2$. The background phase is interpolated from 2,049 samples and checked on their interleaved midpoints. The two-dimensional oscillator evolution preserves $|\alpha|^2-|\beta|^2=1$ to below $2\times10^{-12}$ across the registered scan. Direct evolution of the canonical field agrees in primary occupations within $7\times10^{-14}$ in the checked shells. A separate run with 4,097 background points, $L_{\max}=2048$, relative solver tolerance $10^{-12}$ and smaller blocks changes the primary $D$ by less than $2\times10^{-10}$ on each branch. These are finite-resolution comparisons on the calibrated interpolant. An infinite-tail error bound for the exact fitted-roll spectrum is a separate question.

Run the calculation from the capsule root with the command below.

```bash
OPENBLAS_NUM_THREADS=1 python -IB gauge/gauge.py --output /tmp/GAUGE.json
```
