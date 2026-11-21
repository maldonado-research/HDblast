# Where would the additional pressure relaxation come from?

The shell gas still favors the original beam axis at times 6.67 and 7.05. We held its measured momenta fixed and calculated the *initial* change of directional pressure due to elastic collisions. We compared the archived isolated resonance with the specified pointlike massive-vector exchange model, including same-sign scattering. The calculation transforms both outgoing particles back to the shell frame and measures their actual change in pressure.

At time 7.05 the extended Born prescription gives a pressure difference decreasing about 59 times faster than under the isolated resonance. **99.93% of the additional decrease comes from pairs with $\sqrt{s}/M>1.42749$.** At those energies the candidate model already violates an elastic unitarity constraint in one unlike-charge helicity channel. At time 6.67, 99.80% of the additional decrease comes from this energy region. These results locate the source of the apparent enhancement and prevent a quantitative inference about its physical relaxation time. [Per-batch results](STRESS.json) · [Elastic amplitudes and channel conventions](../crossed-exploration/)

## Pressure in the shell frame

Let $E=\sqrt{1+|p|^2}$ and $q(p)=(p_z^2-|p|^2/3)/E$, with momenta in units of the daughter mass. For both charge signs of one representative species, $Q=\langle q\rangle/\langle E\rangle=(P_{zz}-\operatorname{tr}P/3)/\rho$. The $z$ axis is the imposed primary field direction. The calculation uses the orthonormal shell frame. The eigenpressures of a fitted fluid in its rest frame or of a combination of three species are different observables. All rates below use shell proper time in units of $L$.

| Shell time | Initial $Q$ | Isolated-pole $dQ/d(t/L)$ | Full Born $dQ/d(t/L)$ | Full opposite-charge portion | Full same-charge portion |
|---|---:|---:|---:|---:|---:|
| 6.67 | 0.60175 | −0.4357 ± 0.0042 | −23.111 ± 0.028 | −10.696 | −12.414 |
| 7.05 | 0.43039 | −0.3595 ± 0.0057 | −21.139 ± 0.050 | −8.961 | −12.178 |

| Shell time | Pair-energy range | Isolated-pole derivative | Full Born derivative | Additional derivative, paired estimate |
|---|---|---:|---:|---:|
| 6.67 | Below cut | −0.09207 | −0.13639 | −0.04431 ± 0.00163 |
| 6.67 | Above cut | −0.34365 | −22.97436 | −22.63071 ± 0.02939 |
| 7.05 | Below cut | −0.01836 | −0.03278 | −0.01442 ± 0.00151 |
| 7.05 | Above cut | −0.34112 | −21.10590 | −20.76478 ± 0.04582 |

The cut is $s/m^2=140379347.50895345$. We apply this same kinematic division to all sign combinations. It comes from the first known crossing of the equal-helicity diagonal unlike-charge $J=0$ Born partial wave through $|\operatorname{Re}a_0|=1/2$ in this prescription. It is a sensitivity test. Restricting the integrals to lower energies does not supply a unitary theory in the omitted region or prove that all lower-energy channels are accurate. The smaller added derivatives below the cut are resolved on these finite samples.

The finest snapshots have 16,384 particles of each sign. A second ensemble with 8,192 particles per sign gave full derivatives of −23.045 and −21.338, respectively, with 99.79% and 99.91% of the corresponding added decreases above the cut. The uncertainty labels in the tables are standard errors across eight independent, scrambled quasi-Monte Carlo integrations of each fixed snapshot. The difference between snapshots measures another source of numerical variation.

## Collision generator

Let $\Delta q=q(p')+q(k')-q(p)-q(k)$ and $W=v_{\rm M}\,d(m^2\sigma)/d\Omega$. Define $I_{ij}$ as the mean over one particle of each specified sign of $\int W\Delta q\,d\Omega$. Equal sign densities give

$$
\left.\frac{dQ}{d(t/L)}\right|_{\rm coll}
=\frac{S_0 D}{6\mu^3a^3}\,
 \frac{I_{+-}+I_{++}/2+I_{--}/2}{2\langle E\rangle}.
$$

The halves on the same-sign terms ensure that unordered incoming pairs are counted once. The same-sign differential cross section already includes the separate $1/2$ for identical **outgoing** particles when integrating the full solid angle. Self-pairs in the discretized product have zero invariant flux. The isolated-pole comparator includes only opposite-sign pairs. Elastic collisions conserve the pair energy, so the denominator has zero collision derivative. In this frozen local measurement the expansion, Lorentz force, source births and change of the gas distribution at later times are excluded.

For total spatial momentum $P=p+k$, pair energy $E_T$, invariant mass $W_c=\sqrt s$ and outgoing CM momentum $k_*'$, the exact wall-frame momenta are
$$
p'=k_*'+\left[\tfrac12+\frac{P\cdot k_*'}{W_c(E_T+W_c)}\right]P,\qquad k'=P-p'.
$$
Each trial scattering uses these momenta and its own two outgoing energy denominators. This accounts for particles from a nearly co-moving pair even when their directions in the pair rest frame change substantially.

## Reproduction and checks

The calculation uses 16,384 sampled pairs for each of three sign combinations in each of eight randomized Sobol batches. Ten Gauss–Legendre nodes resolve the polar angle in a coordinate uniform in $\log(1+|t|/M^2)$ on each folded hemisphere. Eight evenly spaced azimuths share a randomized phase. The observable is invariant under exchange of the two outgoing momenta, which permits combining the two hemispheres with their respective differential rates. For rotationally randomized observer axes the expectation is zero. In the two fine completed-Born runs the results were −0.0008 ± 0.0230 and +0.0605 ± 0.0546. Explicit outgoing energy checks agreed to better than $1.4\times10^{-15}$ relative.

A separate control sends two beams into a head-on collision in their center-of-mass frame. The pressure loss for each angle is then $-Q_{\rm in}[1-P_2(\cos\theta)]$. The code compares its fully boosted generator with independent numerical integrals of this analytic expression and the total cross section at $s=M^2,20M^2,100M^2$. On a fixed 4,096-pair sample per sign combination, increasing the physical angular integration from 10 polar nodes and eight azimuths to 18 nodes and 24 azimuths changed the complete derivative by less than $8\times10^{-6}$ relative at both fine endpoints. The larger change for the isolated pole was $3.5\times10^{-4}$ relative. The event rates also agree with exhaustive pair counting within two batch standard errors on all four snapshots. [Generator and built-in controls](generator.py)

Run from this study directory, using the momenta from the axial-current replay and its extended dedicated history.

~~~sh
OPENBLAS_NUM_THREADS=1 python -I -B crossed-stress/generator.py \
  --kernel crossed-exploration/check.py --history OUT/extended/dedicated_half.npz \
  --snapshots OUT/live/snapshots/t6.67000-p14-n4096-r91283-particle-pole.npz \
  OUT/live/cross-snapshots/t7.05000-p14-n8192-r91283-particle-pole.npz \
  OUT/live/snapshots/t6.67000-p13-n4096-r70163-particle-pole.npz \
  OUT/live/cross-snapshots/t7.05000-p13-n8192-r70163-particle-pole.npz \
  --power 14 --batches 8 --polar 10 --azimuths 8 --output crossed-stress/STRESS.json
~~~
