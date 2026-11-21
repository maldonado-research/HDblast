# A homogeneous two-stream Boltzmann experiment

Consider one positive and one negative population in a spatially uniform, locally flat box without an electric field, expansion, injection, escape, or a heat bath. The target number density $n$ counts **one sign in the box frame**. There are $N$ simulation particles per sign, each with the same statistical weight. Only unlike charges collide. Each daughter has mass $m$, and throughout the equations below its four-momentum is in units of $m$. The calculation uses a classical dilute gas and the flavor-specific elastic vector channel of the preceding study. A collision changes the momenta of its two incident particles in three dimensions. Each particle retains its charge and rest mass. [Earlier microscopic model](https://github.com/maldonado-research/HDblast/pull/194)

## Relativistic collision law

The pole mass in daughter units is $M=8300$ and its fractional width is $w=0.30$. For $s=(p_++p_-)^2$, put $\beta^2=1-4/s$ and $g=w(s+2)$. The total and transport cross sections coincide for this parity-even elastic channel. The probability of a collision uses the invariant flux divided by the incoming box-frame energies.

$$
\widehat\sigma=m^2\sigma=\frac{12\pi g^2}{s[(s-M^2)^2+\beta^2g^2]},\qquad
K=v_{\rm M}\widehat\sigma=\frac{s\beta}{2E_+E_-}\widehat\sigma .
$$

In the pair center-of-momentum frame, let $c$ be the cosine between the incoming and outgoing positive particle. The spin-averaged vector currents of equal-mass spin-$1/2$ particles have spatial tensors proportional to $\delta_{ij}-\beta^2\hat k_i\hat k_j$. Contracting the incoming and outgoing tensors gives $A+Bc^2$, with $A=3-2\beta^2$ and $B=\beta^4$. Thus the normalized density in $c$ is $(A+Bc^2)/[2(A+B/3)]$ and the azimuth is uniform. Its expected second moment is $(A/3+B/5)/(A+B/3)$. It is isotropic at threshold and proportional to $1+c^2$ in the ultrarelativistic limit. The solver samples this law by rejection, constructs the two outgoing momenta back to back at $|\boldsymbol q|=\sqrt{s/4-1}$, then boosts them to the box frame. Recomputing the mass-shell energies provides an independent energy check at every event. [Particle solver](particles.cpp)

## Initial distributions and clocks

The **drift-scale** fixture sets the dimensionless incident energy $E_0=35.71061$ to the maximum bulk Lorentz factor in the preceding half-step, flavor-specific trajectory. The earlier gas has its own temperature. Replacing it by cold beams at that speed tests the effect of exceptionally narrow birth momenta. The **hot-energy** fixture sets $E_0=6616.03485$ to the mean particle energy of a stationary Maxwell–Jüttner gas at the earlier hot-radius doubling temperature $T/m=2205.34488$. A cold beam with this energy moves far faster than the actual bulk drift there. The two fixtures match different properties of the previous trajectory. [Archived history](../prior/boost/EVOLUTION.json)

Each collimated monoenergetic preparation has positive particles moving along $+z$ and negative particles along $-z$ at $E_0$. In the split beams half of each sign starts with $0.15(E_0-1)$ of kinetic energy and half with $1.85(E_0-1)$. These preparations match mean energy, number and zero total momentum. Their per-sign directed momenta can differ. A third control samples each sign independently from a stationary Maxwell–Jüttner gas, boosts away the small random total momentum, and rescales the momenta to match the exact starting energy. All preparations for the same energy use the same reference clock

$$
\tau=\frac{n t}{m^2}K_0,\qquad K_0=K(E_0,+p_0; E_0,-p_0),\qquad p_0=\sqrt{E_0^2-1}.
$$

One unit is the starting mean encounter time in the monoenergetic beam. After the distribution evolves, the actual mean rate is different. Measured collisions per particle count accepted pair events divided by $N$. For any collimated positive/negative beam energy mixtures, summing the initial $K_{ij}$ with population weights gives their exact event rate. The loss rate of initial mean directed momentum weights each pair by $(p_i+p_j)/(2\langle p\rangle)$ as well, since the outgoing mean pair difference vanishes. These benchmarks appear in the JSON next to independent random-pair rate estimates.

For a connection to the wall clock, let $a$ be the radius in units of $L$, $D$ be the comoving daughter rest energy in solver units, $\mu=mL$, and $S_0$ the action. The previous model has $nL^3=S_0 D/(6a^3\mu)$ in the wall frame, counting one opposite sign and flavor. Its reference frequency is $L\,d\tau/dt=S_0DK_0/(6a^3\mu^3)$. Between doubling and the first turnaround $D$ is nondecreasing and $a$ grows. Using $D_{\rm double}/a_{\rm turn}^3$ and $D_{\rm end}/a_{\rm double}^3$ gives lower and upper bounds on the integrated reference column. This comparison holds the reference cross section fixed. Actual trajectories in the box retain their full energy-dependent cross section.

## Event sampling and numerical controls

Each ordered positive–negative pair in the homogeneous finite ensemble has rate $K_{ij}/(NK_0)$ per unit $\tau$. The continuous-time algorithm proposes pairs at exponentially distributed times using an exact rejection bound. After each accepted collision the momenta and the candidate rate update immediately. It advances from collision to collision, so the reporting times do not discretize the physics. The RNG for independent sampled rates is separate from the one for dynamics.

For the hot-energy cases a uniform proposal uses $K\le C_u$, where $C_u=\max(108\pi w^2/M^2,48\pi/(M^2-8))$. This follows by using $v_{\rm M}\le2$ and the off-resonance denominator $(s-M^2)^2\ge M^4/4$ below $s=M^2/2$, and the channel bound $\widehat\sigma\le12\pi/(s-4)$ above. For the lower-energy cases a more efficient proposal samples each sign with probability proportional to its energy. Its positive weight is updated after every event. Using $s\le4E_iE_j$ gives $K\le C_w E_iE_j$, with

$$
C_w=\frac{96\pi}{M^4}\max\left(\frac{2.25w^2}{0.99^2},
 \frac{w^2(1+q)}{q}\left(1+\frac{200}{M^2}\right)^2\right),
\quad q=w^2(1-400/M^2).
$$

For $s<0.01M^2$ the off-resonance denominator is at least $0.99^2M^4$. Elsewhere it is at least $(s-M^2)^2+qs^2\ge M^4q/(1+q)$. The code checks the selected bound on every candidate pair. In both samplers the proposal rates and the acceptance fractions multiply to the desired rate $K_{ij}/(NK_0)$.

The standard comparison uses six seeded ensembles of 4096 particles per sign. Four independent reruns at twice the particle count test finite-size sensitivity. A separate discrete-time scheme shuffles the opposite sign each step, proposes one distinct partner per particle, and scatters with probability $1-\exp(-K_{ij}\Delta\tau/K_0)$. This scheme limits each particle to one encounter per step. The three batch sizes 0.1, 0.05 and 0.025 test the resulting bias against the continuous-time history. Reported standard errors describe variation over randomized ensembles. They do not incorporate model uncertainty or the systematic finite-size error.

## What the observables measure

The scaled charge-odd momentum is $(\langle p_z\rangle_+-\langle p_z\rangle_-)/(2\langle p\rangle_{\rm initial})$ for the collimated fixtures. For the isotropic initialization the scale is $p_0$. The pressure anisotropy is $[P_{zz}-(P_{xx}+P_{yy})/2]/{\rm tr}(P)$ with $P_{ij}\propto\sum p_i p_j/E$. Both collimated fixtures start at 1. The fourth Legendre moment of the particle directions provides another angular check. The collision probe also estimates the instantaneous rate of removing charge-odd momentum in units of its initial momentum per unit $\tau$. It can become negative within sampling noise when the total current is small.

The predicted isotropic energy density is proportional to $E\sqrt{E^2-1}\exp(-E/\Theta)$, with $\Theta$ fixed by $E_0=3\Theta+K_1(1/\Theta)/K_2(1/\Theta)$. The supremum distance between its CDF and the empirical particle-energy CDF is computed on the actual sorted particle energies, for the pooled ensemble and for each sign. The largest per-sign distance is archived. The CDF and second moment are checked with independent quadrature. Near equilibrium the energy CDF must fluctuate in a finite ensemble. Comparing it with the same-size thermal control and with doubled particle count helps interpret that floor. Momentum relaxation, angular isotropy and a close energy distribution answer distinct questions. The probes do not test all higher correlations in the distribution.

Even a good approach to classical kinetic equilibrium in this box cannot establish radiation production, chemical equilibrium, release from a gravitational source, or how large a cosmological region a bulk pulse could heat. Spatial gradients, gravitational expansion, electromagnetic acceleration, other particle-production channels and quantum occupation factors remain absent.

## A check on fixed directions when three flavors share the pole

A shared resonance can scatter into a flavor accelerated along a different field axis. For this check populate the three Cartesian axes with the same positive/negative beams and condition on the first collision that produces a chosen flavor whose own axis is $z$. Define $H_4(\hat{n})=n_x^4-6n_x^2n_y^2+n_y^4$. The three equally likely input axes give $H_4(\hat{x})=H_4(\hat{y})=1$ and $H_4(\hat{z})=0$. Since the products of each incoming beam pair are axially distributed about that input axis, the mean contribution equals $H_4(\text{axis})\langle P_4(\cos\theta_{\rm lab})\rangle$. The average includes the center-of-momentum angle, its boost to the lab and the initial pair's $K$.

The mean $H_4$ vanishes for equal monoenergetic head-on beams. It is **0.02069** and **0.14431** for the lower- and higher-energy split beams. Their low-high pairs have a moving center-of-momentum frame. Randomizing every outgoing azimuth continuously about the destination axis would set both values to zero. Keeping a fixed three-axis basis preserves the harmonic. These are initial product distributions in the specified beam fixtures. Their effect on a wall trajectory requires a separate driven evolution. Reproduce with `python -I -B box/azimuth.py --input box/BOX.json --output /tmp/AZIMUTH.json`. [Angular quadrature](azimuth.py) · [Values and order refinement](AZIMUTH.json)
