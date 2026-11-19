# How much of the emitted spectrum must survive?

The test follows the original massless particles that have experienced zero collisions since the shell's emission ends. A surviving particle keeps its comoving momentum $k$. Every other particle makes a nonnegative contribution to the final momentum distribution. This gives information about the final spectrum even when the momenta of all scattered particles are unknown. [Spectral-test protocol](https://github.com/maldonado-research/HDblast/issues/157#issuecomment-7016431913) · [Evolving-target test](https://github.com/maldonado-research/HDblast/issues/157#issuecomment-7016466688)

## A bound for any momentum interval

Let $P_0$ be the normalized one-sign number distribution, $S(k)$ a lower bound on survival, and $Q$ the hypothetical equilibrium distribution with the same energy and particle count. At every native endpoint the signs have equal energy by charge symmetry. A source-free, elastic gas conserves the total comoving energy of both signs. Thus the comparison uses the initial comoving moments, which are constant in time. For disjoint bins $B_i$ the total variation distance obeys

$$
 d_{\rm TV}(P_f,Q) \;\ge\;
 \sum_i\left[\int_{B_i} S(k)\,dP_0(k)-Q(B_i)\right]_+ .
$$

Indeed, the union of bins with positive excess is one of the events over which total variation takes a supremum. Its final probability exceeds its equilibrium probability by at least the displayed sum. Coupling initial and final momenta of the same original particle also gives $d_{\rm TV}(P_f,Q)\ge d_{\rm TV}(P_0,Q)-F$, where $F$ is the fraction that collided at least once. The latter result can use the stronger dynamic bound in the [collision calculation](../envelope/ENVELOPE.md).

For the classical contact benchmark $\sigma=\kappa\hat s$, the relative flux is $1-\cos\theta$. The isotropic loss rate involves only the opposite sign's first energy moment. Even when elastic collisions redistribute the targets, that moment is conserved by charge symmetry and free redshift. The original tag survival is therefore $S=\exp(-Ak)$, with $A=(8/3)y^2(\kappa H_0^4)M_1\int ds/a^5$. Pauli blocking of final states can only raise this survival at the same cross section. For classical equilibrium the normalized momentum density is $k^2e^{-k/T}/(2T^3)$, with $T=M_1/(3M_0)$ in comoving units. A finite amount of scattering can be necessary even to permit a desired spectral agreement.

On the quiet wall this test excludes distance ≤0.10 whenever $A\langle k\rangle<0.27427$, and distance ≤0.01 whenever $A\langle k\rangle<0.52718$ (values rounded; see numerical lower endpoints). Passing those necessary exposures allows many outcomes. No sufficiency or actual approach rate follows.

## Pauli bound with evolving collision partners

Here the benchmark allows any isotropic elastic cross section bounded by $16\pi/\hat s$, for one distinguishable channel, on homogeneous, statistically independent incident particles. The Pauli cap $0\le f\le1$ for each spin state and the conserved number $n=y^2 M_0$ yield

$$
 k_F=(6\pi^2 n/g)^{1/3},\qquad
 I_{\rm target}(t)=\frac{g}{2\pi^2}\int_0^\infty k f(k,t)\,dk
 \;\le\; \frac{g k_F^2}{4\pi^2}=I_{\rm max}.
$$

For fixed number, moving occupation from higher to lower momenta increases the average of $1/k$. The completely filled sphere up to $k_F$ attains the maximum. This applies even if the target distribution at a later time is otherwise unknown. With $B=8\pi I_{\rm max}\int ds/a$, the first-collision hazard of an original tag integrates to at most $B/k$. Jensen's inequality gives $F\le 1-\exp[-B M_{-1}/M_0]$ by particle count. Weighting by initial particle energy gives $1-\exp[-B M_0/M_1]$. Units of momentum are $H_0$ and $g=2$ for each sign. On the quiet wall these conservative bounds are 0.26704% and 0.23363% at $y=0.1$. The smaller fraction from tracking the scattered population separately is in the [dynamic-target calculation](../envelope/ENVELOPE.md).

For the fermionic equilibrium shape at the strictly conserved count and energy, we solve for the fugacity and comoving temperature of the isotropic massless Fermi–Dirac distribution. The quiet fugacity is approximately 0.000174. The quiet initial distribution differs from this equilibrium by at least 0.27891 in total variation on quarter-unit bins. Using the conservative Pauli survival probability directly leaves a distance of at least 0.27709 at shell time 2.16. A larger distance is compatible with the bound. These figures refer to probabilities over single-particle momentum. A temperature associated with an equilibrium reference describes no measurement of the shell.

## Numerical controls

The reference bins have widths 1, 1/2 and 1/4 over $0\le k\le64$, with an overflow bin retaining the equilibrium tail. Source counts use a piecewise angular integral. Within every bin, subintervals of width 1/64 bracket the survival weights at both ends. A second source and time resolution checks the integrated masses. The code independently evaluates the filled-sphere moments and solves for Fermi–Dirac equilibrium using the convergent fugacity expansion. It reports both ends of the numerical bracket for the exposure threshold of each finite-bin test. All bounds concern this finite source, recorded time interval and declared transport assumptions. [Calculation](survival.py) · [All four histories](SURVIVAL.json)

## How much energy would a larger target population require?

We can relax the assumption that the emitted antiparticles are the only targets and ask for a **necessary energy budget** for a proposed replacement. Keep the original freely redshifting tags until their first encounter, the measured geometry, and the same one isotropic elastic channel. The total target pool has two fermionic spin states, $g=2$, and a homogeneous, isotropic one-body distribution satisfying Pauli occupation. Its number and momentum distribution may vary in time. Let $L=\int I(\eta)\,d\eta$ be its total inverse-momentum column. A tag with initial comoving momentum $k$ experiences its first collision with probability no larger than $1-\exp(-8\pi L/k)$. Thus a necessary column for ten percent to collide solves

$$
 0.10=\int [1-e^{-8\pi L_{10}/k}]\,dP_0(k).
$$

For any target inverse moment $I$, the smallest permitted comoving energy density is $E_{\min}=2\pi^2 I^2/g$. To see this, the energy cost per unit inverse moment at momentum $k$ is $k^2$. Filling the lowest available states first minimizes the cost. A filled sphere of radius $r$ has $I=gr^2/(4\pi^2)$ and $E=gr^4/(8\pi^2)$, attaining the bound. Averaging over the conformal interval and applying the mean-square inequality gives

$$
 \overline E_{\rm target} \ge \frac{2\pi^2}{g}
       \left(\frac{L_{10}}{\Delta\eta}\right)^2 .
$$

The quiet wall needs $L_{10}=0.0670490$ in native units. The available interval is 0.0889485. The required mean comoving energy density in opposite-sign targets is therefore at least 5.608 native units. The emitted antiparticles supply 0.13908 at $y=0.1$, giving a ratio of **40.32**. The three pulse ratios are 40.54, 39.17 and 40.24. An independent momentum histogram brackets the required column using tag survival at both edges of each cell. Its conservative endpoints still require **more than 39 times** the emitted one-sign energy on every wall.

This is the total target energy averaged over the observed conformal interval. The count threshold concerns original particles with an encounter; ten percent need not suffice to make a thermal spectrum. The comparison fixes the shell motion so it can specify what a proposed source would have to supply. Introducing that much material and allowing it to gravitate would require a new shell evolution. Changing the target degeneracy, interaction channels, angular law or time interval changes the budget. [Exploratory energy-budget test](https://github.com/maldonado-research/HDblast/issues/157#issuecomment-7016666861)
