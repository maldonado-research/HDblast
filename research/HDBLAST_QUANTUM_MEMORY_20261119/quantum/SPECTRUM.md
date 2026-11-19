# From driven work to quantum pair energy

On the fine quiet wall, the earlier driven oscillator calculation gave **2.57397203** units of conformal energy in the comoving band 60–120. The same parent fluctuations also enter the quantum calculation. Quantum emission samples an ordered product of the parent field, and that ordering contributes **−1.42145800** in the same units. The result is **1.15251403** units of conformal pair energy. At the final readout the scale factor is 1.23570620, giving **0.93267641** units of proper pair energy per comoving volume. An independent integral over this band agrees with the difference of the two archived cutoff energies to $7.0\times10^{-11}$ on the quiet wall. [Prior driven energy](https://github.com/maldonado-research/HDblast/blob/67b2e6ab85bb14d09f248415fdf31b7b66326c31/research/HDBLAST_DISPERSION_MATCHING_20261119/driven_check/ENERGY.md) · [Archived cutoff calculations](https://github.com/maldonado-research/HDblast/blob/67b2e6ab85bb14d09f248415fdf31b7b66326c31/research/HDBLAST_SCALAR_CLOCK_GRAVITY_20261119/source/SOURCES.json) · [Recomputed results](SPECTRUM.json)

## Which parent correlation does emission measure?

Consider one real parent mode with canonical fork variables $r=(\chi_0,\Pi_0)^T$, $[\chi_0,\Pi_0]=i$, zero mean, and symmetrized covariance $V=\left(\begin{smallmatrix}X&C\\ C&P\end{smallmatrix}\right)$. Evolve the two real solutions $h,k$ for unit initial coordinate or unit initial momentum, so $\chi(s)=h(s)\chi_0+k(s)\Pi_0$. The amplitude for creating a bath excitation of conformal frequency $\omega$ contains the parent operator

$$
Z=v_\chi\chi_0+v_\Pi\Pi_0,\qquad
(v_\chi,v_\Pi)=\int ds\,F[\phi(s)]\,(h,k)e^{i\omega\eta(s)},\qquad d\eta=ds/a .
$$

The conformal source is $aF\chi$, giving the proper-time integration factor $F\chi\,ds$. Starting the bath in its vacuum makes its excitation probability proportional to $\langle Z^\dagger Z\rangle$. The covariance and canonical commutator give three distinct real quantities,

$$
\begin{aligned}
S&=\tfrac12\langle Z^\dagger Z+ZZ^\dagger\rangle=v^\dagger Vv,\\
K&=\tfrac12\langle Z^\dagger Z-ZZ^\dagger\rangle=-\operatorname{Im}(v_\chi^*v_\Pi),\\
P_{\rm emit}&=S+K,\qquad P_{\rm reverse}=S-K .
\end{aligned}
$$

The classical ensemble of real drives measures $S$. On a fixed trajectory, $K$ depends only on the canonical commutator and free evolution, so it is identical for every admissible parent preparation. Both orderings have nonnegative spectral weight. The sign of their difference can depend on the waveform.

The spectral factor for conformal **energy** is $\omega(\omega^2-q^2)d\omega/(8\pi^2)$, multiplied by the archived weight of the parent momentum $q$. In the oscillator calculation the same factor is $g_j^2/2$. Proper energy per comoving volume at readout is conformal energy divided by $a_f$. The results are weak-coupling coefficients with $y^2$ factored out. [Vacuum spectral measure and Hamiltonian](https://github.com/maldonado-research/HDblast/blob/67b2e6ab85bb14d09f248415fdf31b7b66326c31/research/HDBLAST_CAUSAL_YUKAWA_ENERGY_20261119/theory/THEORY.md)

## Four trajectories and two independent integrations

The external scalar gate is calibrated to proper times 1.35–2.05 on the quiet trajectory. Each pulse history crosses the same scalar endpoints at its own proper times. The calculation reconstructs this mapping directly from the archived scalar histories. It integrates the canonical parent and conformal clock, then evaluates the Fourier integrals separately across the measured time intervals. All 24 native momenta are included.

| Fine retained wall | Symmetrized conformal energy | Commutator term | Ordered conformal energy | Proper band energy at $s=2.10$ |
|---|---:|---:|---:|---:|
| Quiet | 2.57397203 | −1.42145800 | 1.15251403 | 0.93267641 |
| Plus | 2.56313516 | −1.41446043 | 1.14867473 | 0.92959827 |
| Minus | 2.63300992 | −1.44722724 | 1.18578268 | 0.95964447 |
| Far | 2.57784191 | −1.42028855 | 1.15755335 | 0.93680362 |

For each arm, the direct band integral runs from 60 to 120 using 96 frequency nodes. Two additional integrals run separately from each $q$ to 60 and to 120 using 64 nodes per range, reproducing the original frequency discretization. The reconstructed absolute energies agree with the earlier archive within $3.4\times10^{-9}$. Across all four walls, the largest absolute difference between the direct band energy and the difference of the archived values is $2.2\times10^{-10}$. The agreement checks both the ordering sign and the final scale-factor conversion. [Numerical implementation](spectrum.py)

## The two signed-frequency branches

To check the ordering independently, let $I=a_0^3\sqrt{(q/a_0)^2+m_0^2}$ and $f=(h-iIk)/\sqrt{2I}$. The real parent field has the expansion $\chi=f\,a+f^*a^\dagger$. Define $A=\int ds\,F f e^{+i\omega\eta}$ and $B=\int ds\,F f e^{-i\omega\eta}$. With $n=\langle a^\dagger a\rangle$ and $c=\langle aa\rangle$ from the fork covariance, one obtains

$$
\langle Z^\dagger Z\rangle=n|A|^2+(n+1)|B|^2+2\operatorname{Re}(cAB),\qquad
K=\tfrac12(|B|^2-|A|^2).
$$

The independent canonical and branch calculations agree at every frequency node. On the quiet retained wall, the parent-annihilation branch contributes $0.00045676$, the parent-creation branch contributes $1.18327885$, and their interference contributes $-0.03122159$, in conformal units. Holding the wall fixed and setting the initial coordinate–momentum covariance $C$ to zero, as in the archived erase preparation, changes the total to $1.16017615$. Preparing the instantaneous fork vacuum gives $1.18315134$. Each of these three results contains the same $-1.42145800$ commutator contribution. The retained covariance determinant exceeds the pure one-mode value $1/4$ by $1.1$–$1.5\times10^{-12}$ across the archive. This matters for uses that distinguish exact purity from its numerical approximation.

A stationary control uses a parent and one bath oscillator of frequency 4, a constant scale factor of one, and the same sine-to-the-fourth pulse shape with duration $2\pi$. The pulse harmonics give exactly $B=0$ at this selected resonance and finite $A$. A stationary vacuum then gives zero emission, while the symmetrized calculation gives $0.34697828$ before the cancelling commutator is included. At parent occupation two, the emitted spectral weight becomes $1.38791312$. The cosmological-history calculation includes a time-dependent mass, gate and scale factor. A nonzero fork-vacuum branch there cannot by itself identify which of those time dependences supplied the energy.

## Replay

From the study directory, run

~~~sh
OPENBLAS_NUM_THREADS=1 python -IB quantum/spectrum.py --repo ../.. --output /tmp/SPECTRUM.json
~~~

The inputs are the native fork, the archived fine wall trajectories, the published scalar-clock source readouts, and the earlier quiet driven-work result. The quantum check uses only the specified high-frequency band and previously modeled backgrounds. An observation linking that band to the origin of the universe would require a separate physical prediction and data.
