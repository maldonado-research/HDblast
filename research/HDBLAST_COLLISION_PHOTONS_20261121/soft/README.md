# Soft photons from the recorded charged-particle collisions

At model time 7.05 the largest dedicated-interaction run had recorded 154,876 elastic collisions. For each event we calculated the leading low-energy photon emission expected if the daughters carry ordinary four-dimensional Abelian charge. In the decade from 0.1% to 1% of each pair’s available center-of-momentum photon-energy scale, the emission amounts to **$6.90\times10^{-17}$ of the final charged-gas energy**. Four particle samples and time grids give $6.90$–$6.97\times10^{-17}$. The calculation includes interference between radiation from both charges before and after each encounter. [Event-level calculation](soft_photons.py) · [Results and controls](SOFT.json)

## Numerical result

For dimensionless squared pair energy $s$ in daughter-mass units, define $\Omega=(s-4)/(2\sqrt{s})$. This is the kinematic photon-energy endpoint for an on-shell equal-mass outgoing pair. Both limits of each tabulated band refer to photon energy $\omega$ in that individual pair’s center-of-momentum frame. The reported fractions refer to conformal energy in the shell frame. They count photons potentially emitted at the events, prior to any absorption, propagation or transfer into other modes.

| Wall time | Accepted events in largest sample | $10^{-5}<\omega/\Omega<10^{-4}$ | $10^{-4}<\omega/\Omega<10^{-3}$ | $10^{-3}<\omega/\Omega<10^{-2}$ | Expected photons per charged particle per decade |
|---|---:|---:|---:|---:|---:|
| 6.67 | 51,015 | $1.52\times10^{-19}$ | $1.52\times10^{-18}$ | $1.52\times10^{-17}$ | $1.15\times10^{-14}$ |
| 7.05 | 154,876 | $6.90\times10^{-19}$ | $6.90\times10^{-18}$ | $6.90\times10^{-17}$ | $3.52\times10^{-14}$ |

The median pair at 7.05 has $\sqrt{s}=9{,}295$ daughter masses and a boost of 1.92 from its center-of-momentum frame to the wall. For the upper listed decade, the conformal emitted energy is $3.4602\times10^{-14}$ against a final gas energy of 501.55. This is 4.66% of the separate single-particle smooth-field radiation estimate on the same trajectory. The two estimates condition on how the charge accelerates and how contributions from different times combine. A prediction of their joint signal would require an emission model that resolves interference between successive accelerations. [Recorded event histories](../events/) · [Smooth-field calculation](https://github.com/maldonado-research/HDblast/tree/9467ffbc67e9540da6762fecd1c44e5e22c80522/research/HDBLAST_ENERGETIC_PARTICLES_20261121/radiation)

## Four charged legs and a finite mass

Use $m=c=\hbar=1$ for the local kinematics, metric signature $(+,-,-,-)$, $\hat k=(1,\mathbf n)$, and $\alpha=3Q^2\mu^2/(4\pi S_0K)=1.4438220891\times10^{-16}$ for the coupling to the putative photon. In the following expressions $u_i$ are four-velocities. The signed incoming-positive, incoming-negative, outgoing-positive, outgoing-negative factors are $c=(-1,+1,+1,-1)$. The leading external-charge current is

$$
j^\nu(\mathbf n)=\sum_{i=1}^4 c_i {u_i^\nu\over u_i\cdot\hat k},
\qquad
{dN\over d\ln\omega\,d\Omega_{\mathbf n}}=
{\alpha\over4\pi^2}(-j^2).
$$

The sum of the four $c_i$ vanishes. The polarization sum can therefore use the transverse spatial current. For $g_{ij}=u_i\cdot u_j$ let $L(g)=g\,\mathrm{acosh}(g)/\sqrt{g^2-1}$ with $L(1)=1$. The angular average of $g_{ij}/[(u_i\cdot\hat k)(u_j\cdot\hat k)]$ is $L(g_{ij})$. Evaluating the integral in one particle’s rest frame gives $\mathrm{atanh}(\beta_{ij})/\beta_{ij}$, an equivalent form. Consequently $B=dN/d\ln\omega=-(\alpha/\pi)\sum_{ij}c_i c_j L(g_{ij})$. These are applications of the leading soft-photon theorem conditional on the elastic process. [Low, 1958](https://doi.org/10.1103/PhysRev.110.974) · [Weinberg, 1965](https://doi.org/10.1103/PhysRev.140.B516)

Writing $k_*^2=s/4-1$ and $z=\cos\theta_*$ for scattering of the positive charge gives

$$
B={4\alpha\over\pi}\left[
 L(1+2k_*^2)-1+L(1+k_*^2(1-z))-L(1+k_*^2(1+z))\right].
$$

For forward scattering $B=0$. Backscattering gives $8\alpha[L(1+2k_*^2)-1]/\pi$. At small $k_*$, $B\pi/\alpha=(16/3)k_*^2(1-z)+O(k_*^4)$. The reflection identity $[B(z)+B(-z)]/2=4\alpha[L(1+2k_*^2)-1]/\pi$ explains why suppressing incoming/outgoing interference can barely alter the sum for an even elastic angular law, even though it would give an incorrect answer for a particular event. The largest late ensemble's energy integral differs by just 0.00144% from that reflection average. The implementation retains the interference event by event.

## Counting energy in the wall frame

Each equal-mass pair is back-to-back before and after the underlying elastic event in its center-of-momentum frame. In this frame the spatial current for either charge dipole is proportional to $2\beta\hat{\mathbf v}/[1-\beta^2(\hat{\mathbf v}\cdot\mathbf n)^2]$. The full radiation intensity satisfies $I(\mathbf n)=I(-\mathbf n)$. Full-solid-angle emission in a photon energy interval independent of direction therefore carries zero *mean* three-momentum there. A boost into the wall frame multiplies its mean energy by $\Gamma_{\rm cm}=(E_++E_-)/\sqrt{s}$. Directional selections or photon cuts in the wall frame require a different integral.

An event occurring at wall radius $a_e$, with macroparticle weight $w=D/(2N)$, contributes

$$
\Delta U_e=w a_e(E_++E_-){1-4/s\over2}B\, (\xi_{\rm hi}-\xi_{\rm lo}),\qquad
\Delta N_{\gamma,e}=B\ln(\xi_{\rm hi}/\xi_{\rm lo}),
\qquad \xi=\omega/\Omega.
$$

The energy sum uses each accepted event once. The expected number of photons per charged particle is $\sum_e\Delta N_{\gamma,e}/(2N)$. At leading order every decade contributes the same expected photon count and carries ten times the energy of the preceding decade. The calculation uses positive infrared cutoffs so the stated counts are finite.

## Validation and interpretation

An independent two-dimensional angular integration directly sums both photon polarizations. Four overlapping angular grids resolve the narrow cones around the four charged legs. These integrals reproduce the closed expression at five controls from just above threshold through $\sqrt{s}=8300$, including a scattering cosine of 0.99. Their maximum relative discrepancy is $5.2\times10^{-7}$. The calculated odd angular moments are below $7\times10^{-17}$ of the zeroth moment. Computing all six independent Minkowski products directly from the incoming and outgoing wall-frame momenta separately reproduces the center-of-momentum coefficient on every event of every run, with a maximum relative difference below $3.0\times10^{-12}$.

At a band’s upper edge, an initial-state photon would change the pair's squared invariant mass by $2\sqrt{s}\omega=\xi(s-4)$. For $\xi\le 0.01$ this shift is at most 0.0334 times the magnitude of the elastic model’s complex resonance denominator on the sampled events. This checks the resonance scale. A quantitative error estimate for finite-energy photons also requires the radiative matrix element, including its interference terms and recoil. The output supplies explicitly labeled wider-band extrapolations for comparison with such a calculation. The physical soft-band results assume photon emission from independent localized events, the archived elastic distributions, fixed electromagnetic coupling, and no depletion or backreaction. Overlapping formation times, medium response and a nonlocal wall gauge field can change their interpretation. [Formation intervals](../formation/) · [Finite-energy calculation](../hard/)

~~~sh
OPENBLAS_NUM_THREADS=1 python -I -B soft/soft_photons.py \
    --artifacts /path/to/regenerated/dedicated/event/tapes --output /tmp/SOFT.json
~~~
