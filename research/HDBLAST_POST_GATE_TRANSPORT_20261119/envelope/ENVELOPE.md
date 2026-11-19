# A collision bound for one isotropic elastic channel

On the four modeled walls the available post-emission time is short. For an illustrative production coupling of 0.1, **fewer than 71 particles per million** can undergo a first collision with an independently sampled antiparticle by shell time 2.16 when the target spectrum redshifts freely. Allowing already scattered particles to become new targets still gives fewer than 74 per million in the closed elastic gas below. The bound includes the entire emitted spectrum with total pair frequencies up to 120. It uses the leading production density and allows the elastic cross section to saturate a single partial wave's unitarity ceiling at every collision energy. [Numerical results](ENVELOPE.json) · [Calculation](envelope.py)

## Counting and weighting the particles

Let $p$ be a daughter's comoving momentum in units of $H_0$ and $dP$ the ordered pair number per comoving volume with $y^2 H_0^3$ factored out. For an unpolarized scalar producing a massless Dirac particle and antiparticle, the spin-summed intensity is isotropic in the pair rest frame. If the pair has total comoving energy $\omega$ and momentum $q$, the momentum of one daughter is uniform on $[(\omega-q)/2,(\omega+q)/2]$. Each pair supplies one target antiparticle. Its distribution equals the particle's after angular averaging. For either charge species the relevant moments are
$$
N=\int dP,\qquad M_1=\int \frac{\omega}{2}\,dP,\qquad
M_{-1}=\int\frac{1}{q}\log\frac{\omega+q}{\omega-q}\,dP .
$$
Here $M_1$ is half the total conformal energy coefficient. The logarithm in $M_{-1}$ comes from averaging $1/p$ over each interval. At a fixed native $q>0$, $dP/d\omega$ vanishes at least linearly in $\omega-q$, as follows from the vacuum spectral factor $\omega^2-q^2$ and the finite-time source integral. The product with the endpoint logarithm is integrable. These statements use all 24 archived parent momenta. [Source spectrum with its original units](https://github.com/maldonado-research/HDblast/blob/49756c922760880a9f5da0e4da6f2d8cce36738a/research/HDBLAST_QUANTUM_MEMORY_20261119/quantum/SPECTRUM.md)

## Cross section and expansion

Assume a homogeneous gas of independently sampled particles and antiparticles on the collisionless wall background. Each tagged particle scatters elastically on antiparticles drawn from this population. Take an isotropic cross section in the collision center-of-momentum frame, so that the total and transport cross sections agree. We impose
$$
\sigma_{\rm tot}=\sigma_{\rm tr}\le\frac{16\pi}{\hat s},\qquad
\hat s=2 H_0^2 p k(1-\cos\theta)/a^2.
$$
The coefficient is the distinguishable two-body elastic $J=0$ ceiling. For a single allowed helicity channel $k_{\rm cm}=\sqrt{\hat s}/2$ and $\sigma_0=(4\pi/k_{\rm cm}^2)\,|(S_0-1)/(2i)|^2$, with $|S_0|\le 1$. In a set of coupled spin channels, each unit-normalized incoming state's summed elastic probability obeys the same single-$J$ ceiling, and spin averaging can reduce it. The choice permits all initial spins to reach the ceiling as a conservative bound. Higher partial waves, strong forward scattering and additional targets require their own cross sections. The production interaction switches off at the gate. The interaction operating afterward is an explicitly chosen extension.

The massless relative flux factor is $1-\cos\theta$. Its product with the ceiling has a finite limit at forward incidence. Multiplying by the target physical density $y^2H_0^3 dP/a^3$ gives
$$
\frac{\Gamma(p)}{H_0}\le\frac{8\pi y^2 M_{-1}}{p\,a},\qquad
{\cal T}(p)\le\frac{8\pi y^2 M_{-1}}{p}\,\Delta\eta,\qquad
\Delta\eta=\int_{s_{\rm off}}^{s_f}\frac{ds}{a(s)},\quad s=H_0\tau.
$$
Using the full daughter distribution for tags and targets then yields
$$
\langle{\cal T}\rangle_N\le 8\pi y^2\Delta\eta\,\frac{M_{-1}^2}{N},
\qquad
\langle{\cal T}\rangle_E\le 8\pi y^2\Delta\eta\,\frac{N M_{-1}}{M_1}.
$$
Their ordering follows from $M_{-1}M_1\ge N^2$. Replacing $\Delta\eta$ by $1/[a(H/H_0)]$ gives the respective instantaneous bounds on $\langle\Gamma/H\rangle$. The cumulative first-collision probability in this test-particle description is $1-e^{-{\cal T}(p)}$. Its concavity makes $1-\exp[-\langle{\cal T}\rangle]$ an upper bound on the weighted fraction that has scattered. The energy weighting measures the share of the initial comoving energy carried by particles with a first collision. It assigns no temperature or energy exchange per collision. We set final-state blocking to zero in the collision estimate, which allows the largest rate for these fixed incoming distributions.

| Wall | Gate closes at shell time | Upper bound on mean depth per $y^2$ at 2.10 | Same bound at 2.16 | At $y=0.1$, particle fraction at 2.16, per million | Energy fraction at 2.16, per million |
|---|---:|---:|---:|---:|---:|
| Quiet | 2.050000 | 0.003149 | 0.006890 | 68.90 | 60.27 |
| Near positive | 2.050089 | 0.003131 | 0.006858 | 68.57 | 59.97 |
| Near negative | 2.049257 | 0.003246 | 0.007046 | 70.46 | 61.70 |
| Distant | 2.049887 | 0.003155 | 0.006897 | 68.96 | 60.33 |

For the quiet wall the full one-species moments $(N,M_1,M_{-1})$ are $(0.76981625,13.90822508,0.04871048)$. If both tag and target must originate from the upper pair band 60–120, the ceiling on the number-weighted depth becomes $0.0000378654\,y^2$. If the **tag alone** comes from that upper band and all emitted antiparticles can scatter, the corresponding ceiling is $0.00336892\,y^2$. Adding the lower-frequency targets substantially changes the estimate. Because all terms in the ceiling use the same factorized kernel, the low–low, mixed and high–high contributions to that ceiling can be calculated separately. Actual cross sections can change those shares.

The calculation compares two independent frequency coordinates at the lower-band threshold. One uses $\omega=q+(60-q)u^2$ and another ordinary Gauss integration in $\omega$. The largest relative difference across the three moments and four walls is below $4\times10^{-10}$. We also varied the source-time quadrature and checked the energy against both archived cutoff energies. These checks address numerical precision on the archived momentum grid. They supply no estimate for omitted cosmological species, production above the cutoff, scattering during emission or a later expansion history. A same-pair final-state interaction is a correlation of the produced pair and cannot be estimated by multiplying independent one-body distributions. Full thermalization would require a specified collision dynamics and its competition with expansion over the relevant epoch.

## Allowing previously scattered particles to become new targets

The fixed-target result prompted an exploratory question. Could a small group that has already collided redistribute to low momenta, become a more efficient target, and substantially increase the collision count? For this test the homogeneous, closed population contains the same two massless charge species with two spin states each. Elastic scattering preserves the number of each charge. We allow arbitrary redistribution consistent with Pauli occupation and retain the same isotropic single-wave cross-section ceiling and independent-particle kinetic description.

Use comoving momenta throughout. Write $n=y^2 N$, $I_0=y^2 M_{-1}$, and $E_0=y^2 M_1$ for either species with the production coupling now included. Let $F(\eta)$ be the fraction of the **original particles** of a given charge that have experienced their first collision by conformal time $\eta$. Charge symmetry gives equal $F$ for the particle and antiparticle. Equivalently, one may bound the larger of the two fractions for equal initial marginals and densities.

An ensemble of fermions with number $m$ and $g$ spin states obeys a bound on its inverse momentum moment $I=\int d n(p)/p$. The largest $I$ comes from filling the lowest momenta first with occupation one in each spin state, since moving a particle from high to available low momentum increases $1/p$. Thus
$$
I\le B(m)=\frac{3m}{2p_F(m)},\qquad
p_F(m)=\left(\frac{6\pi^2 m}{g}\right)^{1/3}.
$$
Before its first collision a particle retains its initial comoving momentum. Its surviving number density at each momentum is bounded by the original density there, so the surviving inverse moment $I_{\rm surv}\le I_0$. Particles that have already scattered form a positive component of the total occupancy, subject to the same Pauli bound. Their number is $nF$ irrespective of the number of later collisions. Together these facts give $I_{\rm target}\le I_0+B(n)F^{2/3}$ at each time.

The loss rate of *uncollided* tags per comoving number $n$ is bounded by $8\pi I_{\rm surv} I_{\rm target}/n$ per unit conformal time. Set $x=8\pi I_0^2(\eta-\eta_{\rm off})/n$ and $C=B(n)/I_0$. Then
$$
\frac{dF}{dx}\le 1+C F^{2/3},\qquad F(0)=0.
$$
Because the reciprocal of the right-hand side is positive, integration bounds $F(x)$ by the solution $U$ of
$$
x=\int_0^{U}\frac{dv}{1+C v^{2/3}}
 =\frac{3}{C}\left[U^{1/3}-\frac{\arctan(\sqrt C\,U^{1/3})}{\sqrt C}\right].
$$
This comparison carries the initial spectrum forward only for tags before their first event. The scattered targets can continue interacting throughout the interval. It also bounds the fraction of the *original* energy carried by particles that ever scatter. In fact $I_{\rm target}\le I_0(1+C U^{2/3})$ and $dU/dx=1+C U^{2/3}$ imply $8\pi\int I_{\rm target}\,d\eta\le n U/I_0$. The energy-weighted initial mean of $1/p$ is $n/E_0$, giving an upper first-event energy fraction $1-\exp[-U n^2/(I_0E_0)]$.

| Wall | Evolving-target particle fraction at 2.16, per million | Evolving-target initial energy fraction at 2.16, per million |
|---|---:|---:|
| Quiet | 71.65 | 62.67 |
| Near positive | 71.31 | 62.36 |
| Near negative | 73.30 | 64.18 |
| Distant | 71.72 | 62.73 |

All rows use $y=0.1$. On the quiet wall $C=38.8079$ and $x=0.0000689027$. Using the full Pauli bound on every target without distinguishing survivors gives a weaker particle-fraction ceiling of 0.0026704. The tighter ceiling is 0.000071650. A numerical initial-value solution and the inverse integral agree to better than $2\times10^{-13}$ in absolute fraction on the four histories at both endpoints.

### How long would ten percent require?

The same comparison bounds the exposure needed to reach a **chosen collision fraction** $f$. A necessary conformal column since switch-off is
$$
\Delta\eta_{\min}(f)=\frac{n}{8\pi I_0^2}\int_0^f\frac{dv}{1+C v^{2/3}}.
$$
For ten percent at $y=0.1$, the quiet shell gives $\Delta\eta_{\min}=26.491$ in the native dimensionless conformal clock, compared with $0.0889485$ available through $s=2.16$. It requires at least 26.402 additional conformal units, or about **298 times the observed total column** from switch-off. Across the four histories the necessary total ranges from 26.177 to 26.571, or 292–299 times the respective observed durations. The inequality deliberately allows independent rate maxima to be saturated at each step. It supplies a necessary condition; the actual time can be longer, or the chosen fraction may never be reached. Converting this column to a proper time would require the subsequent scale factor. A long-lived microscopic interaction would also have to continue satisfying the same assumptions.

## Replay

From the study directory, use the pinned inputs from the preceding release.

~~~sh
OPENBLAS_NUM_THREADS=1 python -IB envelope/envelope.py --repo ../.. --output /tmp/ENVELOPE.json
~~~
