# Pair multiplication and interference after a second pulse

On the quiet shell, opening the scalar interaction a second time raises the number of pairs by **26.8%** and their energy by **61.3%** in the tested band $q<\omega<120$. As the added particles carry more energy on average, the final population is farther below the particle count of a gas in chemical equilibrium at the same energy. These are leading-order source predictions on the four prescribed expanding walls. [Registered source comparison](https://github.com/maldonado-research/HDblast/issues/157#issuecomment-7017450151) · [Numerical results](EMISSION.json) · [Downloadable daughter spectra](SPECTRA.csv)

## Two windows on one scalar clock

For every wall we map its scalar $\phi$ to the proper time $u(\phi)$ at which the quiet wall reaches the same value. Each gate is $\sin^4[\pi(u-u_l)/(u_r-u_l)]$ within its window and zero outside. Gate A uses the earlier interval $[1.35,2.05]$, and B uses $[2.055,2.145]$. Two controls repeat the unit-peak shape separately over the first and second halves of B. The half windows preserve **peak strength**. Their shorter durations change both the duration-integrated strength and the Fourier bandwidth. The common readout is proper time $s=2.16$, after all four recorded histories have crossed both gates.

For a parent mode at comoving momentum $q$, write its two free solutions as $h,k$, with initial variables $r=(\chi_0,\Pi_0)$. Their ordered quantum covariance is $W=\langle r r^T\rangle$. In the Fourier amplitudes of each gate $i$, both oscillations have the same conformal-clock origin.

$$
 v_i(\omega,q)=\int ds\,F_i[\phi(s)]\,(h(s),k(s))^T e^{i\omega\eta(s)},
 \qquad d\eta=ds/a.
$$

Define matrices for one-sign comoving number, both-sign conformal energy and an inverse moment by

$$
 G^{(m)}_{ij}=\sum_q w_q\int_q^\Lambda d\omega\,
 \frac{\omega^2-q^2}{8\pi^2}\,m(\omega)\,v_i^\dagger Wv_j ,
 \qquad m\in\{1,\omega,(\omega/2)^{-3}\}.
$$

For signed real gate coefficients $c$, the pair count and total conformal energy are $N=c^T G^{(1)}c$ and $E=c^T G^{(\omega)}c$. Each contributes one particle and one antiparticle. The initial daughter field is the conformal vacuum. The complete matrices in the machine-readable results span A, B and the two half gates. A random sign for B averages away the off-diagonal contribution, producing the incoherent count $G_{AA}+G_{BB}$. The coherent sign difference measures $4\operatorname{Re}G_{AB}$. Each signed case uses the same geometry, prepared parent and real coupling.

For both charge signs, $n=2|y|^2 H_0^3N/a_f^3$ and $\rho=|y|^2 H_0^4 E/a_f^4$. For a charge-symmetric massless Dirac gas in chemical equilibrium at the same energy, $\rho_{\rm eq}=7\pi^2 T^4/60$ and $n_{\rm eq}=3\zeta(3)T^3/\pi^2$. The ratio $n/n_{\rm eq}(\rho)$ therefore scales as $\sqrt{|y|}$, independent of subsequent collisionless redshift. This ratio compares particle counts. Assigning an actual temperature would require a separate equilibration test. [Earlier one-pulse particle distribution and kinetic test](https://doi.org/10.5281/zenodo.24375564)

## Number grows more slowly than the chemical target

The table gives pair count in units of $y^2H_0^3$ per normalized comoving volume and conformal energy in units of $y^2H_0^4$. The final two columns show fractions at the illustrative $y=0.1$ for $q<\omega<120$.

| Wall | A pairs | A+B pairs | A energy | A+B energy | A fraction of equilibrium count | A+B fraction |
|---|---:|---:|---:|---:|---:|---:|
| Quiet | 0.769816 | 0.976278 | 27.81645 | 44.85620 | 12.229% | 10.837% |
| Near positive | 0.766784 | 0.973443 | 27.71383 | 44.76872 | 12.214% | 10.822% |
| Near negative | 0.782563 | 0.989464 | 28.26106 | 45.34480 | 12.284% | 10.895% |
| Far | 0.769900 | 0.976689 | 27.82320 | 44.89047 | 12.228% | 10.836% |

On the quiet wall, B alone gives 0.206260 pairs and 17.03615 energy units. An individual daughter has mean proper energy $33.08H_0$ from B alone, compared with $14.47H_0$ from A or $18.40H_0$ from A+B. The high average energy of the new particles increases the reference equilibrium number faster than it increases the actual count.

## Interference changes the daughter spectrum

The integrated quiet count for A+B differs from the incoherent result by only **0.0207%**. We also counted particles in 64 unit-wide comoving momentum bins. A scalar emitting massless Dirac daughters is isotropic in the pair rest frame. For each $(q,\omega)$, an individual daughter has uniform comoving momentum magnitude on $[(\omega-q)/2,(\omega+q)/2]$. We integrated the accepted angles exactly across each bin to obtain one-sign number and energy and the per-spin mean occupation. The CSV lists proper momentum edges at the readout in addition to comoving edges.

After inspecting those spectra we calculated the total variation distance, $\frac12\sum_b|P_b-Q_b|$, with each binned particle distribution normalized to unit total. This treats the count as given and tests the *shape* of the distribution. Opposite signs for B yield a distance of **0.15874** on the quiet wall, with a range **0.15858–0.15878** across the four walls. The distance for positive B versus the incoherent prediction is **0.07935** on the quiet wall. The alternating signs of the interference pattern largely cancel when all momenta are simply added. For an illustration selected after inspecting the bins, $25<k<26$ holds 2.54% of the quiet incoherent count, and the signed count contrast in that bin is 42.6%. These are predictions for a hypothetical particle measurement given a controllable change to the source. They provide no identification of cosmic dimensionality.

The earlier occupation-norm bound applies to each real two-gate combination by summing the positive time-domain norm integrals over both disjoint gates. In the A+B example the pointwise per-spin bound is $0.018915 y^2$ on quiet and at most $0.019047 y^2$ across the four walls. These bounds concern the calculated initial emitted sector. The density matrix contains pair correlations beyond a one-particle spectrum.

## Energy above the gravity cutoff

The short pulse has considerable power at frequencies above the cutoff 120 used by the coupled gravity calculation. The following comparison extends the upper limit of the **source integral** with the same fixed wall, parent momenta and unit gates.

| Quiet cutoff $\Lambda$ | B pairs | B conformal energy | First half-gate energy | A+B conformal energy |
|---|---:|---:|---:|---:|
| 120 | 0.206260 | 17.03615 | 14.31669 | 44.85620 |
| 240 | 0.430656 | 53.47201 | 106.37220 | 81.28958 |
| 480 | 0.436975 | 55.08959 | 204.69816 | 82.90709 |
| 960 | 0.436975 | 55.08996 | 205.18897 | 82.90745 |
| 1920 | 0.436975 | 55.08996 | 205.18981 | 82.90746 |

Cutoff 120 includes **47.20% of B's pair count and 30.92% of B's energy** relative to 1920 on the quiet wall. At 1920 the A+B chemical count fraction is **8.452%** at $y=0.1$. For the first half-gate, the shorter unit-peak pulse at 120 includes only **6.98%** of its computed energy at 1920. Frequencies beyond 1920 are outside this numerical comparison. The large-band results retain the same prescribed wall and free parent modes. A gravity calculation including these additional frequencies would need its own stresses, vacuum terms and cutoff matching.

## Independent checks and replay

We propagated the two real canonical parent modes through both gates from the same fork. The maximum canonical determinant error was below $3\times10^{-11}$. The calculation integrates positive square-root amplitudes of the ordered covariance, splitting the proper-time quadrature at native trajectory knots. The pair-frequency integration uses ten-unit panels above 10 and a logarithmic map over $q<\omega<10$ to resolve the inverse cube near threshold. Raising the quadrature from 12 to 16 nodes per frequency panel and from four to six time nodes per native interval changes the tested counts, energies and chemical fractions by at most $5.5\times10^{-10}$ relative across all profiles and cutoffs. The direct A-only number and energy below 60 and 120 agree with the previous archive to $1.7\times10^{-12}$ relative. The binned number and energy use independent spectral interpolation and piecewise angular antiderivatives. Raising the source interpolation degree from 80 to 112, along with the time refinement, changes those bin integrals in normalized $L^1$ by at most $1.2\times10^{-10}$ across the 52 wall–source combinations.

The command below runs from the study root.

~~~sh
OPENBLAS_NUM_THREADS=1 python -IB emission/emission.py --repo ../.. --output /tmp/hdblast-repump-emission
~~~

[Calculation](emission.py) · [Gram matrices, strengths and convergence tests](EMISSION.json) · [Binned spectra](SPECTRA.csv)
