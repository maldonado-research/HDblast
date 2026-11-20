# Quantum pairs on the closed wall

During expansion the spatial gauge potential of the published wall passes through zero. Charged particles in a given angular shell then have exactly equal instantaneous energy, whatever their alignment with the electric field. Starting in vacuum, the model produces very different populations in those directions. This tests whether a specified collisionless field can have a thermal distribution.

## Field and state

Take one massless complex scalar confined to the wall, with conformal curvature coupling $R|\Phi|^2/6$ and fixed charge $q$ under one of the three Abelian fields. The original equal triplet was normalized so its sum has electromagnetic energy $K(X^2+B^2)/2$. One invariant direction therefore has potential $A/\sqrt3$, where $A=a^2 B/2$. Along the dense forward history $d\eta=dt/a$ and $A'=a^2 X$. The calculation uses temporal gauge and the curvature units of the input. The chosen charges describe a test species.

On a unit round three-sphere, ordinary scalar harmonics have degree $n\geq0$ and eigenvalue $n(n+2)$. They are the matrix entries of the spin-$n/2$ representations of $SU(2)$. The eigenvalues of a unit invariant momentum component are $k=-n,-n+2,\ldots,n$. Left and right group actions commute, leaving $n+1$ states at each $k$. In the covariant spatial Laplacian one invariant momentum is shifted by $s=qA/\sqrt3$. After rescaling by $a$, the conformal curvature term supplies the remaining $+1$ in the oscillator frequency. With $j=n+1$ the exact test equation is

$$u_{jk}''+\omega_{jk}^2 u_{jk}=0,\qquad \omega_{jk}^2=j^2-2ks+s^2.$$

Here $\omega_{jk}^2=(s-k)^2+j^2-k^2>0$. The initial empty field has $u(0)=1/\sqrt{2j}$ and $u'(0)=-i\sqrt{j/2}$. The neutral conformally coupled massless scalar stays in its vacuum on the expanding background. The instantaneous occupation per degenerate state is $N_{jk}=|\sqrt\omega\,u-iu'/\sqrt\omega|^2/2$. The antiparticle is paired with the conjugate label. Each two-mode squeezed state has probability $1/(1+N_{jk})$ of retaining its vacuum. For a finite angular band the probability of at least one pair is $1-\exp[-\sum_{jk} j\log(1+N_{jk})]$. The conformal excitation energy of both charges is $\sum_{jk}2j\omega_{jk}N_{jk}$. A physical energy density requires division by $2\pi^2a^4$ and a choice of physical units before comparison with the wall tension.

## Equal energies, unequal populations

At wall time **10.48052**, radius **9.53159**, the electromagnetic share of the wall matter is **91.65%** and $A=0$, while $A'=-0.0810418$. Each state in the $n=2$ shell has conformal energy $3$. Its angular populations for $q=1$ are

| $k$ | $-2$ | $0$ | $+2$ |
|---:|---:|---:|---:|
| Occupation per state | $1.15874\times10^{-5}$ | $2.17376\times10^{-10}$ | $1.05460\times10^{-5}$ |

The largest occupation exceeds the smallest by **53,306 times**. A Bose–Einstein population at a common temperature and chemical potential in the wall rest frame assigns equal weights to these degenerate states. The $k=0$ channel first responds to a term quadratic in $qA$, and its weak-charge occupation begins at order $q^4$. The other two channels respond already at order $q^2$ in occupation. For $q=0.1$ the population ratio exceeds five million. Interactions among produced quanta are absent.

This angular test resolves a single charged species. Combining three identical species directed along orthogonal axes and their antiparticles can conceal the angular difference. For $n=2$, the sign-averaged distribution along one direction is an even function of a spin-one generator, of the form $c_0+c_2J_i^2$. The sum over the three axes is proportional to the identity. A second diagnostic survives that averaging. The mean occupation **per state** is $4.71038\times10^{-8}$ at energy 1 and $2.82394\times10^{-5}$ at energy 2. The higher level is **599.515 times more populated per state** in this instantaneous basis. Rotating the species preserves each shell trace. Every positive-temperature Bose–Einstein population falls with energy for fixed chemical potential, and averaging such distributions preserves that ordering.

The field continues to change at this instant. A finite-time particle count depends on the basis. To specify out-vacua for a second check, keep the measured history until the first expansion maximum at $t=12.42163$, $a=9.74401$. After that, prescribe a readout over conformal duration $L$, with derivative $A'(\eta)=A'_{\rm turn}[1+\cos(\pi(\eta-\eta_{\rm turn})/L)]/2$, then hold the potential constant. These smooth continuations are detector conventions.

For $q=1$, $0\leq n\leq32$ and one complex species on the entire sphere, the chance of at least one pair is **0.1621%** for $L=0.125$ and **0.1634%** for $L=0.25$. The expected numbers of pairs are 0.001622 and 0.001635. These probabilities refer only to the stated band. Shorter waves contribute appreciably to the energy.

## Short wavelengths and timing

Weak-charge perturbation theory permits a larger mode sum. At first order, up to a phase,

$$\beta_{jk}=-\frac{qk}{2\sqrt3 j^2}\int A'(\eta)e^{-2ij\eta}\,d\eta+O(q^2).$$

Write the integral as $F_j$. Since $\sum_k k^2=j(j^2-1)/3$, the total expected **particle plus antiparticle** number in shell $j$ is $q^2(1-j^{-2})|F_j|^2/18$ at leading order. Multiplying by $j$ gives the leading conformal excitation energy.

| Highest $j$ | Short smooth readout, energy divided by $q^2$ |
|---:|---:|
| 32 | 0.04864 |
| 512 | 3.29304 |
| 2,048 | 27.23586 |
| 8,192 | 44.23312 |
| 16,384 | 44.25366 |
| 32,768 | 44.25479 |

Half the resolved energy sits at $j\leq1{,}706$. Doubling the input trace spacing changes the final coefficient by $5.0\times10^{-7}$ relative. The long smooth readout gives 44.25499. For the interpolated short readout, the entire tail beyond $j=32{,}768$ contributes at most **0.183** to this leading coefficient. For $D=A'$ continued by zero, $|F_j|\leq\mathrm{TV}(D')/(4j^2)$ by two integrations by parts. Each cubic segment and the analytic continuation give an explicit total variation of $336{,}153$. Thus the energy-coefficient tail is bounded by $\mathrm{TV}(D')^2/(576J^2)$. A sharper numerical Parseval estimate from the independent integral of $|D'|^2$ and the resolved spectrum bounds the remaining coefficient by about $0.000217$. These results describe the weak-charge limit on the supplied trajectory. A fundamental ultraviolet cutoff remains unspecified.

Most of this energy follows the short incoming pulse. Split the electric history $D$ into $wD$ and $(1-w)D$, where $w$ is one through wall time 1, zero after wall time 2, and a decreasing quintic smoothstep in conformal time between those instants. This window lies before the scalar gives its energy back to the classical waves. Through $j=32{,}768$ the short-readout energy coefficients of the two pieces are **44.25423** and **0.00023157**. Their interference contributes **0.00032714**. The three terms reproduce the full coefficient. The partition is a diagnostic of a coherent amplitude.

Freezing the potential abruptly at turnaround changes the high-frequency result. Here $D_f\ne0$ and the transform contains the endpoint term $F_j\sim-D_f e^{-2ij\eta_f}/(2ij)$. The instantaneous-basis energy grows logarithmically with upper mode number, with coefficient $q^2D_f^2/72$. Subtracting the last octave of the smooth sum from the last octave of the raw sum recovers this predicted logarithmic increment to within 1%. This ultraviolet tail prevents interpreting the abrupt readout as a total heat density.

## Reproduce

From the study directory run

~~~sh
python -I -B quantum/pairs.py --trace background/trace.csv.gz --output /tmp/PAIRS.json
~~~

The exact solver composes symplectic oscillator propagators and compares two time meshes, the reverse history and charge conjugation. Direct adaptive integration independently checks the homogeneous oscillator and the degenerate shell used for the population test. Halving the time step changes populations by at most $2.7\times10^{-11}$ across the exact crossing calculation. The empty-field maximum occupation is below $10^{-28}$. The normalized Wronskian is accurate to $5\times10^{-14}$ for the $q=1$ primary runs. Born and exact finite-band occupation sums agree at $q=0.01$ to $3.3\times10^{-7}$ relative. This calculation omits particle backreaction, collisions, annihilation and secondary emission. A temperature, plasma conductivity and matter energy fraction would require further physical input. [Program](pairs.py) · [Machine-readable data](PAIRS.json) · [Forward background](../background)
