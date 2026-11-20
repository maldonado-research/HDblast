# Continuous parent decay and electromagnetic heating on a gravitating wall

A neutral population can release charged particles through a decay threshold controlled by the wall scalar. On a selected trajectory the daughter particles acquire 91.49% of the added wall energy as internal motion while the wall expands. The decisive energy transfer is electric work. The example assumes a pre-existing wall and parent population, a deliberately chosen scalar-dependent parent mass, a large charge-to-mass ratio and a collisional fluid closure. The cold parent abundance and the scalar's work on its mass enter the gravitational simulation.

## Charge and energy supply

There are two copies of the bulk and a closed wall initially at rest with $a_0=2$. The principal finite source pulse has dose 0.123, and the post-exit electromagnetic coefficient is 13,200.04. The late wall law is inherited from the preceding study. The previous proper-time conversion switch has been replaced by the evolution below. The neutral parents are at rest on the wall. Each decay produces one positive and one negative daughter of equal fixed mass $m$. The microscopic benchmark assigns equal branching fractions to three daughter flavors, each coupled to its own gauge field. At birth, each has energy $M/2$ and momentum magnitude $m\sqrt{R^2-1}$ in a uniformly distributed direction, where $R=M/(2m)$. Equal populations across the three independent Abelian fields retain isotropic total stress. A physical identification of the parent, daughter flavors and interactions remains open. [Preceding charged-particle calculation](https://doi.org/10.5281/zenodo.24405438) · [Underlying finite-source wave and wall](https://doi.org/10.5281/zenodo.24387565)

The parent mass ratio is $R(\psi)=R_i+(R_f-R_i)S[(\psi-\psi_*)/w+1/2]$. The infinitely differentiable step is zero for $x\leq0$, one for $x\geq1$ and $S(x)=[1+\exp(1/x-1/(1-x))]^{-1}$ between them. In the primary case $R_i=0.95$, $R_f=1.05$, $w=0.03$, $\psi_*=0.77$. The threshold $R=1$ is crossed at $\psi_*$ for these symmetric endpoints. This choice targets the part of the neutral trajectory with a strong electric field. The activation time is determined by the dynamical scalar and the wall, and shifts when the incoming pulse changes.

A Dirac fermion pair coupled to a scalar parent through a scalar Yukawa interaction gives the vacuum rest-frame width $g^2M(1-4m^2/M^2)^{3/2}/(8\pi)$ for an open channel. For three equal-mass channels set $g^2=\sum_i y_i^2$, where $y_i$ is the Yukawa coupling of one flavor. In proper wall time measured in units of the bulk curvature length $L$, the total width gives
$$
 \dot f=\widehat\Gamma(1-f),\quad
 \widehat\Gamma=\kappa R[1-R^{-2}]_+^{3/2},\quad
 \kappa=\frac{g^2mL}{4\pi},\qquad f(0)=0.
$$
Here $f$ is the cumulative fraction of parents which have decayed. This local on-shell width neglects inverse reactions, expansion and curvature corrections to individual transitions, final-state blocking, and radiative corrections. The one-species width and the equal-branching total width share this form when the effective coupling is matched. A choice of $\kappa$ fixes a product involving $mL$. For example, $\kappa=100$ and $g^2/(4\pi)<1$ require $mL>100$. A prediction for a physical interaction strength needs a scale and charge normalization; $Q=|q|/(\sqrt3m)$ refers to the model's electric-field units. The bulk normalization and electromagnetic coefficient affect the charge measured by real matter.

Write $D_0=10^{-4}$ for the *initial parent rest energy* in tension units and $B_*=(D_0/R_i)(2/a)^3$ for the daughter rest mass density that complete decay would supply. With $f$ decayed, the parent energy density is $B_*R(1-f)$ and the daughters contribute rest mass $B_*f$. The parent pressure is zero. Keeping the wall initially stationary gives $\mu_0=4-16D_0(2+D_0)=3.99679984$ on each side. The scalar acceleration contains $-B_*(1-f)R_{,\psi}$. This is the force that balances the work when surviving parents become heavier. The pressure, mass and temporal junction, reflected bulk waves and geometry evolve with the total matter stress.

## Daughter transport and exact balances

The charged matter consists of six tilted perfect-fluid streams, balanced in sign and orientation. Their internal pressure follows the relativistic Taub–Mathews closure $h(\Theta)=5\Theta/2+\sqrt{1+9\Theta^2/4}$, with $h$ the specific enthalpy including rest mass and $\Theta=p_0/(n_0m)$ in one stream's rest frame. The proper electric field of one member of the triad is $X/\sqrt3$. The summed current doing work on the wall is $J=B_*Qfv$ with $v$ the signed positive-stream velocity. The gas receives energy $JX$ per proper wall volume and time. Pulse, geometry and three gauge fields respond to the resulting current. A coefficient $\nu$ relaxes the relative momenta of the positive and negative components. The corresponding energy is conserved in the summed gas and the field except for tracked fluxes and expansion work. The model assumes rapid within-stream relaxation when it incorporates new particles into one temperature and drift per stream. The decay rate, this relaxation and the finite inter-stream drag are distinct assumptions. [Taub–Mathews equation of state](https://arxiv.org/html/0704.1679#S2.E13)

Define $E=\rho_g/B_*$, $Z=fh\gamma v$ (the magnitude with sign of the directed momentum in equivalent one-stream units) and $P=p_0^{\rm sum}/B_*$, with proper pressures summed across streams. The isotropic *wall-frame* gas pressure is $p_g=B_*(P+Zv/3)$. Distribution of newborn particles over directions gives zero mean injected momentum. The conservative equations are
$$
 \dot Z=Q f X-(H+\nu)Z,\qquad
 \dot E=Q f v X-H(3P+Zv)+R\dot f .
$$
For stable evolution through zero abundance and the ultrarelativistic limit, the code advances $U=E-C$, $C=\sqrt{f^2+Z^2}$, with $v_c=Z/C$. In the zero-density limit these have continuous fluxes. The stable velocity difference for positive abundance is $v-v_c=-v_c(U+P)/(E+P)$. It gives
$$
 \dot U=Q fX(v-v_c)-H[Z(v-v_c)+3P]+\nu Zv_c+(R-f/C)\dot f.
$$
Use $f/C=1,v_c=0$ for a new, initially empty gas. The primitive inversion is one monotone cubic. Write $z=f/E$, $A=1-(Z/E)^2$, and $p=P/E$. The Taub equality for the proper energy and pressure gives
$$
 3p^3+(3A+2+z^2)p^2+(A+2z^2)p-A(A-z^2)=0.
$$
All nonconstant coefficients are positive for a populated physical state. The cubic has a unique nonnegative root for $A\ge z^2$. In floating point the difference $A-z^2$ is computed as $U(2C+U)/E^2$. Then $v=Z/(E+P)$, $\gamma=(1-v^2)^{-1/2}$ and $\Theta=P\gamma/f$. Random energy is measured conservatively as $B_* f(h-1-\Theta)$, the proper per-particle random energy above rest mass times lab number density. The larger $\rho_g-B_* f\gamma$ additionally includes boosting the internal energy and pressure. Neither statistic recovers a spectrum or measures a common photon temperature.

Three separately integrated ledgers compare to the instantaneous material and field state. Here $q_b=I^2-O^2$ is the net bulk flux at one wall face and $B_m$ denotes the magnetic field on the wall.
$$
\begin{aligned}
 (a^3\rho_g)^{\cdot}+3H a^3p_g &=a^3 JX+a^3 B_*R\dot f,\\
 (a^3\rho_{\rm parent})^{\cdot}&=a^3B_*(1-f)R_{,\psi}\dot\psi-a^3B_*R\dot f,\\
 (a^4 e_{\rm EM})^{\cdot}&=a^4[2q_b+\dot K(B_m^2-X^2)/2-JX].
\end{aligned}
$$
Since $a^3B_*=8D_0/R_i$ is constant, the total directly injected birth kinetic energy is at most $8D_0(R_f-1)/R_i=4.21053\times10^{-5}$ when $R_f\ge1$, independent of the positive reaction rate. Complete rest-plus-kinetic injection is bounded by $8D_0 R_f/R_i=0.00088421$. During expansion with vanishing electric work, positive pressure further reduces the comoving excess energy. The mass-change work has already entered the parent balance and must not be added again to direct injection.

## Results and controls

At 64 early radial elements, the wall bounces at radius **1.881040**. At time **10.341103** it reaches the decay threshold at radius **9.51647** with 99.887% of the surface electromagnetic energy in the electric field. Half the parents decay by **10.73140** and 90% by **11.27212**. Intrinsic daughter random energy crosses 90% of the total added wall energy at **11.26695**, radius **9.68672**. It reaches **91.4888%** during expansion and the lab thermal excess reaches 91.4890%. At the intrinsic peak $\Theta=1760.77$, $f=0.996752$ and $v=-0.001474$. The pressure-to-density ratio of the daughters is 0.33333332. The wall expands to radius **9.76159** at time **12.42354** and begins contracting. The single expanding interval above 90% gives **0.773%** growth, or 0.0077001 e-folds. It is far shorter than the radiation era of standard cosmology.

At the intrinsic peak, comoving electrical work on the daughters totals **4.513107**, direct decay energy including rest totals **0.000877375**, and direct birth kinetic energy totals **0.0000380044**. The growing parent mass drew **0.0000802464** of comoving energy from the scalar. The daughters have $a^3\rho_g=4.433822$ and their accumulated pressure work is **0.0801627**. Accordingly the gas energy and expansion work agree with the electric and decay terms. The measured electrical work is over 5,100 times all the directly injected energy and over 118,000 times the kinetic part. This comparison uses the same proper-volume normalization for all four amounts. The separate electromagnetic ledger, weighted by $a^4$, records work of 43.27274 and residual coherent energy of 0.08563 at the peak.

The table shows the peak intrinsic random fraction of all added energy in expanding post-exit states. Each control shown uses 24 early radial elements and six late ones, with the same parameters as the primary unless otherwise specified.

| Change | Peak intrinsic fraction |
|---|---:|
| Reference ($\psi_*=0.77, w=0.03,\kappa=100,Q=10^8,\nu=6$) | 91.49% |
| Set $\kappa=0$, or keep mass below threshold | 0% |
| Set daughter charge and drag to zero | 0.000762% |
| Set explicit inter-stream drag to zero | 4.67% |
| Reduce $Q$ by a factor of ten | 59.45% |
| Reduce $\kappa$ to 0.1 or 10 | 2.86% or 88.15% |
| Raise $\kappa$ to 1,000 | 91.66% |
| Set $\psi_*$ to 0.42, 0.64, 0.76 or 0.79 | 1.45%, 53.23%, 91.61% or 88.79% |
| Raise $\psi_*$ to 0.85 | 0% on the reported trajectory |
| Increase $w$ from 0.03 to 0.12 | 90.66% |
| Change pulse dose from 0.12300 to 0.12295 or 0.12305 | 70.03% or 69.51% |

The threshold is an explicit function of an evolving scalar. The chosen mass law, rate and scalar trajectory still control when charges become available. In particular, changing the pulse by around four parts in ten thousand can move the threshold well away from the favorable electric phase in this parameter region. The scans illustrate dependence on supplied physics and do not estimate a probability of a hot universe.

## What the single-fluid calculation misses

A gas born continuously from two-body decays already has random momentum at birth. Mixing populations born at different times also spreads their momenta. In the run with $\nu=0$, the one-fluid closure still approximates each tilted charge class by one internally equilibrated component. We computed an independent collisionless current on that run's **prescribed** $A=a^2B_m/2$ and scale factor to test that approximation. For birth at time $b$ with $p_b/m=\sqrt{R_b^2-1}$ and $\mu$ uniform on $[-1,1]$, the impulse on a positive charge and its transverse redshift give
$$
 p_\parallel/m=\frac{a_b\sqrt{R_b^2-1}\,\mu+Q(A-A_b)}{a},\qquad
 (p_\perp/m)^2=\frac{a_b^2}{a^2}(R_b^2-1)(1-\mu^2).
$$
The magnetic force rotates transverse motion about the invariant field axis. It preserves its magnitude after redshift. Curvature on the round three-sphere introduces no local spatial acceleration term proportional to the symmetric product $p_i p_j$ in the invariant orthonormal frame, whose spatial connection is antisymmetric in those indices. The test integrates over $df$ and all launch angles, retaining the velocities of particles of different ages. At time 12 the fluid's current divided by $Q B_*f$ is **−0.99931**, while the resolved collisionless particles on its field yield **−0.20553**. Halving the temporal resolution changes this normalized kinetic current by at most 0.00363 over eight samples; halving the angular resolution changes it by less than $2\times10^{-15}$. The discrepancy means the collisionless charges would source a different electric field. This comparison fixes the field to isolate the current response; it cannot establish energetics of an alternative fully coupled trajectory. Anisotropic momentum dispersion or a suppressed mean current can coexist with an unthermalized set of particle beams. [Collisionless calculation](cohorts.py) · [Momentum and pressure diagnostics](COHORTS.json)

The closure models an ultrarelativistic, stable charged population with a prescribed relaxation rate. At these energies a realistic visible-sector model would need cross sections, screening, radiation, other species, the quantum state and an energy-dependent collision kernel. The parent decay law and its effect on gravity are now calculable within the model. An observational temperature, a photon blackbody spectrum, sufficient duration and the origin of the parent population remain open.

## Reproduce

The full radial Maxwell and Einstein constraints are evaluated at every wall time step. The primary evolution refines the initial pulse on 64 elements and uses six thereafter. Compared with 48 early elements, the peak intrinsic fraction changes by about $1.4\times10^{-6}$ and the 90% duration by about $7.2\times10^{-6}$ e-folds. Increasing late resolution from six to ten elements and reducing the time step on the 24-element early grid changes the intrinsic fraction by about $10^{-7}$. The primary maximum Friedmann error is below $5.4\times10^{-11}$; parent/gas comoving energy errors are below $1.9\times10^{-17}$ and $8.2\times10^{-11}$. The discrepancy from the exact number law $f=1-\exp(-\int\widehat\Gamma\,d\tau)$ is below $9.4\times10^{-14}$. The minimum bulk mass is 0.09654 and minimum lapse 0.99778. The monitored tube mass discrepancy is below $1.7\times10^{-6}$. A separate build with floating-point contraction disabled and different optimization reproduces the reported peak, threshold, duration and gas energy figures.

Run from this directory with Python 3, a C++17 compiler and NumPy for the angular quadrature.

~~~sh
python -I -B decay.py --output /tmp/decay-results.json --jobs 2
python -I -B cohorts.py --output /tmp/collisionless-packets.json
~~~

[Parameters and controls](cases.json) · [Coupled solver](decay_forward.cpp) · [All numerical records](DECAY.json) · [Source normalization and runner](decay.py)
