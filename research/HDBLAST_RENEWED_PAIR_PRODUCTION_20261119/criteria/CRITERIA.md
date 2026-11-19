# Can any bounded late pulse close the chemical particle deficit?

For the illustrative production coupling $y=0.1$, the quiet shell's original pulse supplies **12.23%** of the number of massless Dirac particles in a gas at chemical equilibrium with the same energy, using the archive's cutoff at pair frequency 120. Over the remaining measured interval, any sequence of renewed pulses with a combined peak of at most one has an upper limit of **20.91%** at this cutoff. Extending the same leading pair spectrum to arbitrarily high daughter frequencies gives a **22.21% upper limit** on the quiet shell. The three disturbed walls give all-frequency ceilings of 22.20–22.28%. These ceilings apply to the retained 24 parent momenta, their prescribed Gaussian preparation and the recorded time available for the source. [Registered waveform test](https://github.com/maldonado-research/HDblast/issues/157#issuecomment-7017450151) · [Numerical results](CRITERIA.json) · [Replay](criteria.py)

## Counting the energy that additional particles bring

Let $N=\int dP$ be the emitted pair number with the factor $y^2$ removed. Each pair contributes one fermion and one antifermion. Write $E=\int e\,dP$ for the mean conformal energy carried by either sign, with $e=\omega/2$ for a pair of total comoving energy $\omega$. Angular averaging in the scalar decay gives that mean. All moments are per normalized comoving volume in native $H_0$ units. The total fermion plus antifermion energy in this normalization is $2E$.

A charge-symmetric gas of a single massless Dirac flavor in pair-changing equilibrium has zero chemical potential. For each charge, $g=2$ and
$$
n_{\rm eq}(\rho)=A\rho^{3/4},\qquad
A=\frac{3g\zeta(3)}{4\pi^2}\left(\frac{240}{7g\pi^2}\right)^{3/4}=0.2764102451 .
$$
For the source with comoving number $y^2 N$ and energy $y^2 E$, the count ratio is $R=|y|^{1/2}N/(A E^{3/4})$. The common redshift cancels. After any change from $(N,E)$ to $(N',E')$, this ratio improves precisely when $N'/N>(E'/E)^{3/4}$ for positive final moments. For infinitesimal positive additions the energy per extra particle must be below $4E/(3N)$ to improve the ratio. On the original quiet shell at cutoff 120 this threshold is 24.09 in comoving energy units. Reaching $R=1$ is a necessary number test for this specific equilibrium distribution. Momentum-dependent occupations and equilibration rates also enter a thermal interpretation.

## A waveform-independent inequality

Use one auxiliary spectral moment
$$
L=\int e^{-3}dP=\int (2/\omega)^3dP.
$$
Here $e$ is the mean daughter energy at fixed pair momentum and frequency. An individual daughter's inverse cubic momentum moment would give a different quantity. Hölder's inequality with exponents $4/3$ and $4$, applied to $e^{3/4}$ and $e^{-3/4}$ against the same positive measure $dP$, gives
$$
N^4\le E^3 L,\qquad R\le\frac{|y|^{1/2}}{A}L^{1/4}.
$$
Thus chemical number balance requires $L\ge A^4/y^2=0.58373609$ at $y=0.1$. This test already accounts for any energy deposited by a new pulse. The archived quiet pulse gives $L=0.00032033766$ for $q<\omega<120$.

The abstract amplitude at each native $q,\omega$ is $\int ds\,F(s)e^{i\omega\eta(s)}W_q^{1/2}u_q(s)$, where $u_q=(h_q,k_q)^T$ is the free parent transfer row and $W_q$ the ordered canonical covariance at the fork. Its squared norm, times the native $q$ weight $w_q$ and $K_q(\omega)d\omega=(\omega^2-q^2)d\omega/(8\pi^2)$, gives $dP$. The squared norm of the *instantaneous* real parent row is $X_q(s)=\langle\chi_q(s)^2\rangle$. The Hilbert-space norm of the integrated amplitude with an additional weight $(2/\omega)^3$ is $\sqrt L$. These formulas propagate the original phase even when the switch changes sign.

Let $A_0$ denote the original pulse and $G$ a later pulse or sum of pulses. The supports are disjoint. For a finite splitting frequency $\Lambda>\max q$, define
$$
J_q(\Lambda)=\int_q^\Lambda (2/\omega)^3 K_q(\omega)\,d\omega
=\frac{1}{\pi^2}\left[\log(\Lambda/q)-\tfrac12(1-q^2/\Lambda^2)\right],\qquad
B(\Lambda)=\int_{\rm late} ds\,\sqrt{\sum_q w_q X_q(s)J_q(\Lambda)} .
$$
If $|G(s)|\le b$ throughout that interval, the triangle inequality for the weighted amplitude gives
$$
L_{<\Lambda}(A_0+G)\le\left[\sqrt{L_{<\Lambda}(A_0)}+b\,B(\Lambda)\right]^2 .
$$
The early norm is directly integrated. The late bound allows constructive interference, any waveform, and any number of subpulses under the single peak constraint. Its maximizer need not be dynamically attainable. On the quiet wall $B(120)=0.01550565$ for the entire post-gate window. The upper bound on $R$ becomes 0.209095 for $b=1$.

## Controlling all daughter frequencies

The original gravity evolution used a finite comoving cutoff. To examine the emission sector outside it, assume that the same free massless Dirac spectral kernel and prescribed parent row apply to all $\omega\ge q$. Restrict to sources that give a finite positive outgoing number and energy for the chemical comparison. Smooth compact pulses meet this requirement in the specified histories. A discontinuous ultraviolet drive may inject unbounded energy.

For $\omega>\Lambda$, $(2/\omega)^3K_q(\omega)\le1/(\pi^2\Lambda)$. Since $ds=a\,d\eta$, the conformal-time transform acts on the vector $a(A_0+G)W_q^{1/2}u_q$. Plancherel over the entire frequency axis bounds its positive-frequency squared norm. With the Fourier convention above the factor is $2\pi$. The vector can be complex because its covariance is quantum ordered. Pointwise disjoint time supports remove the cross term from the **full-line** time norm:
$$
L_{>\Lambda}(A_0+G)\le \frac{2}{\pi\Lambda}(T_A+b^2 T_G),\qquad
T_A=\int_{A_0}ds\,a A_0^2\sum_q w_q X_q,\qquad
T_G=\int_{\rm late}ds\,a\sum_q w_q X_q .
$$
The actual positive-band production retains cross-pulse interference. Summing the low- and high-frequency bounds therefore supplies an all-frequency ceiling without dropping the coherent contribution. At $b=1$ on the quiet wall the two time integrals are 0.0529791 and 0.0112959 in the full post interval. At split 240, the resulting upper bound on $L$ is 0.00142023.

| Fixed wall | $R$ ceiling through 120, entire post interval | All-frequency $R$ ceiling, entire post interval | All-frequency $R$ ceiling, late scalar-clock window |
|---|---:|---:|---:|
| Quiet | 0.209095 | 0.222093 | 0.211879 |
| Near positive | 0.208992 | 0.221988 | 0.211636 |
| Near negative | 0.209730 | 0.222747 | 0.212890 |
| Distant | 0.209162 | 0.222162 | 0.211882 |

The full post interval starts at each wall's original switch-off and ends at shell proper time $s=2.16$. The shorter window follows the quiet scalar-clock markers 2.055–2.145 on each trajectory. Each entry allows **every** admissible late waveform in its chosen support. We checked splits 24, 48, 120, 240 and 480. Each gives its own valid bound on the same continuum spectrum. The tightest count ceiling among the five uses 240 in all rows. At fixed split 120 the largest continuum ceiling is 0.22413. The values are numerical evaluations of the analytic inequalities on the retained grid.

Solving the quadratic bound in $b$ gives a necessary increase in peak strength if all other conditions are held fixed. For the entire post interval it requires at least a 42.8-fold increase on the most permissive of the four histories. Such a change would take $|yb|$ above four. The leading-order calculation supplies no prediction for that regime. Lengthening the interval, increasing the parent population, imposing a new energy supply, or changing the reaction can alter the bound through the displayed integrals. The result concerns an isolated gas, the specified parent spectrum and finite measured times. A completion of the parent momentum grid and regulator-independent gravitational evolution require separate calculations.

## Replacing the peak constraint with an integrated source budget

The preceding bound constrains the aggregate peak at every time. As a subsequent check we also allowed the waveform to concentrate its strength in short intervals. Define the actual weighted squared source integral and the reference budget
$$
T[G]=\int_J ds\,a G^2\sum_q w_qX_q,\qquad T_{\rm unit}=T[1]
$$
on the same available interval $J$. The reference describes a unit-magnitude coupling active throughout the interval. For arbitrary real $G$ define
$$
D_\Lambda=\int_J ds\,\frac{\sum_q w_q X_q J_q(\Lambda)}{a\sum_q w_q X_q}.
$$
Cauchy–Schwarz applied to the late amplitude gives $\int_J |G|\sqrt{\sum_q w_q X_q J_q}\,ds\le\sqrt{D_\Lambda T[G]}$. The low-frequency contribution and rigorous ultraviolet tail now yield
$$
L_{\rm all}\le\left[\sqrt{L_{<\Lambda}(A_0)}+\sqrt{D_\Lambda T[G]}\right]^2+
\frac{2}{\pi\Lambda}\big(T_A+T[G]\big).
$$
At $T[G]\le T_{\rm unit}$ the largest all-frequency count ceiling on the four histories is 0.225019, even with no separate peak bound in the leading-source functional. For the quiet shell $T_{\rm unit}=0.01129585$ and $R\le0.224330$. Inverting this inequality in $\sqrt{T[G]}$, and taking the strongest necessary condition from the five frequency splits, gives **$T[G]\ge19.86069$** on the quiet history for chemical number balance. This is at least 1,758 times the available reference budget. The corresponding factors on the positive, negative and distant histories are 1,761, 1,735 and 1,755 after rounding. Each result is a necessary condition. A waveform at that strength need not yield a chemical population.

These numbers follow from the squared norm of the time-dependent source that appears in the Fourier amplitude. With dimensions restored, $T$ has energy units (here $H_0$), while $D_\Lambda$ has inverse energy units. The signed work done by the clock depends on the dynamical source response and is accounted for separately. The large formal resource shortfall does not establish what a nonperturbative source would produce.

## Replay and checks

Run from the study directory using the archived inputs:
~~~sh
OPENBLAS_NUM_THREADS=1 python -IB criteria/criteria.py --repo ../.. --output /tmp/CRITERIA.json
~~~
The replay checks the early spectrum with two independent frequency coordinates, refines the source-time integration, and verifies the closed-form frequency measure, the canonical transfer determinant, the earlier cutoff energy and the measured finite-band tails against Plancherel. It evaluates each proposed temporal support on the recorded shell history. [Quantum emission measure](https://github.com/maldonado-research/HDblast/blob/49756c922760880a9f5da0e4da6f2d8cce36738a/research/HDBLAST_QUANTUM_MEMORY_20261119/quantum/SPECTRUM.md)
