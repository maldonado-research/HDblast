# Vacuum scalar production during the narrow source pulse

The HDblast heavy-source family transfers energy from a rolling scalar to a homogeneous neutral oscillator. The curvature of the same potential mixes their two quantum fluctuation fields during the pulse. On the published trajectory, modes with physical momenta between $m/8$ and $64m$ acquire **$0.00989942\,m^4$ of excitation energy per physical volume** in the leading, locally flat calculation. Here $m$ is the daughter mass used to set the time and momentum units. At the published $mL=4\times10^{12}$ action, the band carries **0.896% of the outgoing homogeneous parent energy**. Increasing the action lowers that fraction, but it also increases the number of fermion births requested by the fluid. At the highest action allowed by the published angle-integrated sharp-shell test, the scalar band still costs **0.737%**. [Classical source at the pinned revision](https://github.com/maldonado-research/HDblast/blob/dc99ff5c5894abcff3e11e02e4b97e19da31de00/research/HDBLAST_COHERENT_DECAY_20261120/envelope/ENVELOPE.md) · [Sharp-shell assumptions](https://github.com/maldonado-research/HDblast/blob/dc99ff5c5894abcff3e11e02e4b97e19da31de00/research/HDBLAST_COHERENT_DECAY_20261120/kinetics/SWEEP.md)

## Background, state and energy

Take both wall scalars to have the canonical two-derivative kinetic terms with the common action normalization $S_0=\sigma L^4$. Set $\hbar=c_{\rm light}=1$, $\mu=mL$ and use the source coordinates $s=\mu(\tau-\tau_s)$, $u=\mu(\psi-\psi_s)$ and $c=\mu\chi$. The leading fast potential is $U=m_0^2(c-h)^2/2$, with $m_0=1.9$ and $h(u)=0.008\,b(u)$ the published smooth bump of width $0.24$. The incoming data are $u=c=c'=0$ and $u'=0.12287824489895026$, with primes here denoting $d/ds$. Integrating both classical fields gives a duration $1.9671456754$ and outgoing parent energy $e_c=(c'^2+m_0^2c^2)/2=0.0002716758554102$ in wall-tension units. The driving scalar reacts to the homogeneous parent throughout this trajectory. [Pinned inner solver](https://github.com/maldonado-research/HDblast/blob/dc99ff5c5894abcff3e11e02e4b97e19da31de00/research/HDBLAST_COHERENT_DECAY_20261120/envelope/inner.py)

Write $q=p/m$ for physical momentum in the local frame. Each column of $F_q$ starts in one incoming scalar vacuum mode and satisfies

$$
F_q''+[q^2\mathbf{1}+U'']F_q=0,\qquad
U''=m_0^2\begin{pmatrix}
(h_{,u})^2-(c-h)h_{,uu} & -h_{,u}\\
-h_{,u}&1
\end{pmatrix}.
$$

Before and after the bump the vacuum frequencies in fast time are $\omega_j=(q,\sqrt{q^2+m_0^2})$. Use $F_{ji}(0)=\delta_{ji}/\sqrt{2\omega_j}$, $F'_{ji}(0)=-i\omega_j F_{ji}(0)$. The squared out-amplitudes are $|\beta_{ji}|^2=|\omega_jF_{ji}-iF'_{ji}|^2/(2\omega_j)$. This projection has zero excitation for a freely evolved vacuum. The coherent displacement stays in the homogeneous background. The two real fields each contribute one oscillator per momentum label in the continuum integral. In units of $m^4$ the outgoing vacuum-subtracted excitation density in a band is

$$
C_{[a,b]}=\frac1{2\pi^2}\int_a^b dq\,q^2\sum_{j,i}\omega_j|\beta_{ji}(q)|^2.
$$

| Upper momentum $p/m$ with lower limit $1/8$ | Cumulative energy $C$ |
|---:|---:|
| 1 | 0.00142498 |
| 2 | 0.00510905 |
| 4 | 0.00723442 |
| 8 | 0.00933511 |
| 16 | 0.00987451 |
| 32 | 0.00989920 |
| 64 | 0.00989942 |

For example, the $[1/8,2]$ band alone gives 0.462% of the parent at the higher published mass and action. The energy formula is positive across momenta. Extending a specified quadratic theory to additional momentum bands can only increase this out-energy. Each quoted band requires that theory to apply up to its stated cutoff. [Species-resolved energies, spectra and numerical checks](PULSE_MODES.json)

## Comparing two costs at the same action

The classical parent density is $\sigma e_c=S_0 e_c/L^4$. The ratio of the finite-band density to the classical parent density therefore equals $r=(C/e_c)\mu^4/S_0$, with $C/e_c=36.43834$. This follows independently from the locally rescaled action. Scaling wall time and all three local spatial coordinates by $\mu$ multiplies the inner action by $S_0/\mu^4$. Each canonically normalized fast field is consequently $\sqrt{S_0}\,(u,c)/\mu^2$ in units of $m$. Its homogeneous energy in units of $m^4$ is $(S_0/\mu^4)e_c$.

| $\mu$ | Published trial $S_0$ | $r$ at trial | Published half-step aggregate sweep ceiling | $r$ at that ceiling |
|---:|---:|---:|---:|---:|
| $4\times10^{10}$ | $3.609\times10^{48}$ | 0.002585% | $1.266\times10^{46}$ | 0.73687% |
| $4\times10^{12}$ | $1.042\times10^{54}$ | 0.89560% | $1.266\times10^{54}$ | 0.73688% |

For each fixed history, the aggregate peak demand-to-capacity ratio $D$ for sharp fermion births scales linearly with $S_0$. Its product with $r$ is independent of $S_0$ and equals the last column as a fraction. Thus requiring $D\le1$ sets a floor on this scalar energy fraction. At the smaller mass the published elastic drag floor already exceeds the sweep ceiling. At the larger mass, continuing a family with $S_0$ proportional to $\mu^3$ would make $r$ grow linearly with $\mu$ if the same fast trajectory persists. The numerical comparison uses the two published masses only. [Actions and half-step ceilings](https://github.com/maldonado-research/HDblast/blob/dc99ff5c5894abcff3e11e02e4b97e19da31de00/research/HDBLAST_COHERENT_DECAY_20261120/kinetics/SWEEP.json)

## Numerical and physical reach

The replay evolves the classical background with adaptive DOP853. It evolves the mode amplitudes in the interaction picture by fourth-order Runge–Kutta and integrates the spectrum with positive composite Simpson weights on each dyadic momentum interval. At 1,024 time steps and 128 momentum panels per interval, halving the time resolution changes $C$ by $6.4\times10^{-15}$ and halving the momentum resolution changes it by $1.8\times10^{-9}$. The maximum sampled errors in $\alpha^\dagger\alpha-\beta^\dagger\beta=\mathbf1$ and $\alpha^T\beta-\beta^T\alpha=0$ are $1.5\times10^{-13}$ and $3.2\times10^{-15}$. A separate DOP853 solve in the original field and velocity variables agrees with every complex $\beta_{ji}$ at all ten reported band edges to $1.6\times10^{-12}$ in absolute value. These are numerical convergence checks. [Replay code](pulse_modes.py) · [Receipt](PULSE_MODES.json)

This is a free linear fluctuation calculation on a prescribed classical trajectory. The short duration and $pL\gg1$ for the reported heavy masses motivate a locally flat wall and a continuum momentum density. The result depends on the incoming vacuum, the canonical potential and the chosen bandwidth. An evaluation of physical backreaction would need the renormalized stress and field forces during the pulse, the validity range and interactions of the two-scalar effective theory, and any coupled matter, gauge, metric or bulk modes. The finite out-energy establishes a cost in the stated model and leaves the total correction to the later hot-wall solution open.

From this directory run `python -I -B pulse_modes.py --output /tmp/PULSE_MODES.json` with NumPy and SciPy installed.
