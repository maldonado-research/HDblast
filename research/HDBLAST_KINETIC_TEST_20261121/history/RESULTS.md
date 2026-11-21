# When could the charged particles first have scattered?

Half of the wall's total absolute electrical work has occurred while most of the newly created particles have yet to scatter, even using the collision rate from the assumed thermal distribution. In the dedicated interaction the mean is only 0.250 encounters per particle at that point. In the shared interaction the mean is 0.428. These budgets require that at least **75.0%** and **57.2%**, respectively, have experienced zero encounters on average. By the start of the 90%-hot phase, the mean counts have risen to **63.7** and **33.5**. The count changes substantially between the main period of field work and the later moment when the hot condition is first met. [Numerical results](HISTORY.json)

| Quantity, half time step | Dedicated pole | Shared pole, all outgoing flavors |
|---|---:|---:|
| Mean encounters when half the absolute work is complete | 0.24979 | 0.42808 |
| Minimum zero-encounter share then, conditional on the mean | 75.02% | 57.19% |
| Zero-encounter share then under the population-averaged Poisson model | 78.24% | 65.79% |
| Absolute work completed when the mean first reaches one | 71.22% | 71.05% |
| Fraction of final particles already born by first hot time | 92.59% | 94.91% |
| Mean event count by first hot time | 63.73 | 33.50 |
| Zero-encounter share then under the same Poisson model | 0.1447% | 0.2467% |
| Share with fewer than five encounters then under that model | 0.7552% | 1.3652% |

“Half the work” means half of $\int a^3 |\mathbf J\cdot\mathbf E|\,d\tau$ over the complete numerical run. This includes energy flowing back to the field. At that time 29.88% and 35.59% of the final particles have formed in the dedicated and shared runs. The denominator for an encounter statistic is always **all daughters already born at the stated time**. The source continues afterwards. The calculation counts both charges. Every binary event gives one encounter to each member. It uses a single dedicated pole for each flavor or the inclusive $1/3$ event weight when three equal flavors share one pole. The target density per incoming particle is that of the opposite sign of its own flavor. [Parent wall and flavor calculation](https://github.com/maldonado-research/HDblast/tree/7fe36012da64e98add3915acd0356c733e717670/research/HDBLAST_DIRECTIONAL_BIRTH_20261121)

## Thermal encounter rate

Take two equal Maxwell–Jüttner distributions with temperature $\Theta=T/m$, proper density $n$, and signed wall-frame rapidities $+r$ and $-r$. Use $m=1$ inside the integrals, $\gamma=\cosh r$, and $f_\pm=n\exp[-\gamma(E\mp v p_z)/\Theta]/[4\pi\Theta K_2(1/\Theta)]$. The number of events per unit wall volume and wall time is

$$
\mathcal R=\int\frac{d^3p_1}{E_1}\frac{d^3p_2}{E_2}\,f_+(p_1)f_-(p_2)
\sqrt{(p_1\cdot p_2)^2-1}\,\widehat\sigma(s).
$$

Here $\widehat\sigma=m^2\sigma_{\rm tot}$ uses the isolated angular-even vector channel with $M/m=8300$ and $g^2/(12\pi)=0.30$. The streams are distinguishable. Dividing the event density by the positive stream's wall density $n\gamma$ gives its wall-time hazard $\Gamma$. A separate division by proper target density gives the dimensionless function $K=m^2\Gamma/n$. The opposite sign gives the same answer. In particular there is no enthalpy factor or transferred-momentum weighting in the event integral.

Parameterize the pair invariant mass by $\sqrt{s}=2\cosh\eta$. At fixed total longitudinal momentum $P_z$, integration over transverse total momentum gives a Boltzmann factor $\exp[-A E]$ with $A=\gamma/\Theta$, $E=\sqrt{s+P_z^2}$. The angular integral over the relative incoming momentum gives $\sinh(bE)/(bE)$, with $b=\sinh(r)\tanh(\eta)/\Theta$. Changing variables to $P_z=\sqrt{s}\sinh y$ leaves

$$
J=\int_0^\infty E e^{-AE}\frac{\sinh(bE)}{bE}\,dy
=\frac{K_0(x_-)-K_0(x_+)}{2b},
\qquad x_\pm=\frac{2\cosh(\eta\pm r)}{\Theta}.
$$

The result for one incoming particle is

$$
K(\Theta,r)=\frac{4}{\Theta\gamma^2 K_2(1/\Theta)^2}
\int_0^\infty\cosh\eta\,\sinh^3\eta\;\widehat\sigma(\eta) J\,d\eta .
$$

At zero relative drift $J=\sqrt{s}\,K_1(\sqrt{s}/\Theta)$, giving the usual invariant thermal average for equal masses. With two cold opposite beams $K=2\sinh|r|\,\widehat\sigma(|r|)$. The one-stream density in wall units is $nL^3=S_0D/(6a^3\mu\gamma)$, where $D$ is the comoving daughter rest energy summed over all six sign-and-axis populations, $\mu=mL$, and $S_0$ is the action scale. Thus the wall-time rate is $\Gamma L=(nL^3/\mu^2)K$ times the channel weight. [Rate implementation](event.py)

## Cohorts and the counting bound

Let $B(t)=\int^t\Gamma(t')\,dt'$ and let $dD$ be accepted comoving rest added by a birth cohort. Daughters have equal rest mass, so $dD/D(t)$ is exactly their share of the number present. Under an exchangeable local thermal and Poisson encounter prescription, the mean depth of that cohort is $b=B(t)-B(t_{\rm birth})$, and the share that have never scattered is

$$
F_0(t)=\frac{1}{D(t)}\int^t e^{-[B(t)-B(t_{\rm birth})]}\,dD(t_{\rm birth}).
$$

This equation makes a stronger assumption than the mean count. Given the modeled event density, the mean count for the entire population is $\bar N=D(t)^{-1}\int^t D\Gamma\,dt'$. For **any** distribution over nonnegative encounter counts with that same mean, $\Pr(N\ge1)\le\bar N$ and hence $\Pr(N=0)\ge\max(0,1-\bar N)$. This bound can be informative during the main work phase even if cohorts differ in energy or successive encounters are correlated. Its premise is the integrated collision count evaluated with the inherited thermal background. A genuinely nonthermal background could change that count. A separate diagnostic follows individual newly born particles in the electric field up to their first scatter. [Ballistic calculation](../ballistic/BALLISTIC.json)

## Resolution and interpretation

The four wall replays preserve the hot-onset time, hot fraction and expansion extrema of the pinned dedicated and shared models. Main and half time steps give 63.64 / 63.73 and 33.46 / 33.50 events per particle at hot onset, respectively. The zero-collision lower bounds at half absolute work agree to better than 0.001 percentage point between the time steps. The numerical reconstruction of total absolute electrical work agrees with the native Runge–Kutta ledger to better than one part in a million. The integrated participation count also agrees with the independently weighted birth-cohort count. Near and after radius doubling the cumulative rest has almost stopped changing at floating-point precision. Extremely small residual noncollision fractions there depend on late, very weak source increments and should be treated as resolution limited. [Native history exporter](export.py) · [Cohort integration](optical.py)

The analytic thermal event rate also agrees with an independent direct average over Lorentz-boosted particle pairs. Across five cold, resonant and hot probes, each using eight scrambled Sobol sequences of 131,072 pairs, the largest difference was 0.028%. The check also recovers the cold and hot constant-cross-section limits. Halving both lookup spacings and increasing invariant-mass quadrature order from 16 to 24 changed the reported hot-onset Poisson shares and mean encounter counts by at most 0.0017% relative. At the half-work checkpoint the change was under 0.00051%. [Independent and refinement checks](CHECK.json)

The diagnostic gives a time and a population against which a kinetic formation mechanism can be tested. An abundance of late elastic collisions alone does not establish that earlier field-driven particles had a thermal momentum distribution or that the daughter gas reached local equilibrium. The Poisson model assigns the thermal population mean hazard to every particle present at a given time. A ballistic daughter with a particular energy can sample the narrow range of invariant masses near the pole very differently. [Kernel checks](CHECK.json)
