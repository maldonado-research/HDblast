# A transverse-momentum test of thermal startup

A particle born in the original wall calculation receives little momentum perpendicular to its own electric-field axis. Cross-axis scattering can move it out of this narrow region. The finite interaction rate limits how fast that can happen. On the two tested wall histories there are many more particles still required inside this region than an exactly thermal stream can accommodate.

For example, at **1% of the accumulated absolute electric work**, even the largest collision rate allowed by the adopted vector interaction leaves at least **33.7%** of the dedicated-pole population in the narrow region. A Maxwell–Jüttner stream at the temperature of the published fluid solution puts **2.8 parts in 100 million** there. The corresponding shared-pole bounds are **56.4%** versus **6.3 parts in 100 million**. At these early markers only 2.76% and 3.33%, respectively, of the eventual daughter populations have been born. Fractions below are always relative to particles already present at the marker. [Numerical results](TUBE.json)

| Fraction of whole-run absolute electric work | Dedicated minimum for any target momenta | Dedicated Maxwell–Jüttner fraction | Shared minimum for any target momenta | Shared Maxwell–Jüttner fraction |
|---|---:|---:|---:|---:|
| 0.1% | 88.83% | $2.23\times10^{-7}$ | 94.33% | $4.26\times10^{-7}$ |
| 1% | 33.72% | $2.80\times10^{-8}$ | 56.45% | $6.31\times10^{-8}$ |
| 5% | 5.39% | $8.21\times10^{-9}$ | 12.01% | $1.76\times10^{-8}$ |
| 10% | 2.48% | $4.40\times10^{-9}$ | 5.06% | $8.90\times10^{-9}$ |
| 50% | 0.206% | $2.71\times10^{-10}$ | 0.431% | $9.10\times10^{-10}$ |
| Reference first hot crossing | 0.00634% | $1.13\times10^{-9}$ | 0.0140% | $1.47\times10^{-9}$ |

The minimum columns use the same birth and target density histories, the same parallel electric and magnetic fields and the sole specified binary pole. Targets may have arbitrary momenta. Those columns do not use the Maxwell–Jüttner event rate.

## Why the transverse population is compulsory

In units of daughter mass $m$, injection has $p_b=\sqrt{R_b^2-1}$, with $R_b=M_{\rm parent}/(2m)\le1.05$. Until scattering, the parallel electric force changes only longitudinal momentum, the parallel magnetic force rotates the transverse momentum, and expansion reduces its magnitude. Thus every unscattered particle satisfies

$$
p_\perp/m\le q_0=\sqrt{1.05^2-1}=0.3201562\ldots
$$

The transverse momenta of *all* particles can then differ from an exactly thermal marginal even if the first few fluid moments are similar. The probability difference for this single set bounds the total variation distance from below.

For an isotropic Maxwell–Jüttner distribution with $\Theta=T/m$, use cylindrical momentum coordinates and integrate over $p_z$. Writing $b=\sqrt{1+(p_\perp/m)^2}$ gives $\int_{-\infty}^{\infty}dp_z/m\,e^{-E/T}=2b K_1(b/\Theta)$. A longitudinal boost preserves the count fraction inside the transverse region. In the rest-frame measure its number-density Jacobian is proportional to $1+v p_z/E$, whose odd part integrates to zero at every $p_\perp$. With $b_0=\sqrt{1+q_0^2}$ the exact fraction is

$$
P_{\rm MJ}(p_\perp/m\le q_0)=
\frac{\int_1^{b_0} b^2 K_1(b/\Theta)\,db}{\Theta K_2(1/\Theta)}
=1-\frac{b_0^2 K_2(b_0/\Theta)}{K_2(1/\Theta)}
\sim \frac{q_0^2}{4\Theta^2}\quad (\Theta\gg1).
$$

The code integrates the positive numerator to retain accuracy where the fraction is below $10^{-9}$. An independent integral over particle energy and solid angle checks the result at every marker.

## Two ways to bound the surviving birth population

**From the largest allowed cross section.** The exact pair bound is $m^2\sigma v_{\rm M}<1.118569\times10^{-6}$ for this vector pole. Use the full possible opposite-sign target density in the wall frame and multiply by the appropriate inclusive channel fraction $b_{\rm ch}$. On each interval $[t_i,t_{i+1}]$ along an expanding and continuously producing trajectory, the cumulative accepted rest $D$ increases and the scale factor $a$ increases. The maximum exposure of any tag in that interval therefore satisfies

$$
B_i=\Delta t_i\,\frac{b_{\rm ch}S_0 D_{i+1}}{6(mL)^3 a_i^3}K_{\max},\qquad
P_{\rm unhit}(t_j)\ge
\frac{1}{D_j}\sum_{i<j}(D_{i+1}-D_i)
\exp\left[-\sum_{k=i}^{j-1}B_k\right].
$$

This is an ensemble probability bound for the classical binary-scattering model. It assigns every birth the largest possible exposure over its entire birth cell. It holds for any target velocity distribution on the specified density history. The tabulated endpoint histories and linear interpolation to the reported marker times are numerical inputs. Halving the wall integrator time step moves the independent 1%-work floor from 33.724% to 33.761% for the dedicated interaction and from 56.453% to 56.478% for the shared interaction. [Pair bound and characteristics](BALLISTIC.md)

**From the total number of encounters in the thermal calculation.** Write $C$ for total collision incidences divided by the number of particles born so far. A collision contributes two incidences, one for each participant. At most a fraction $C$ can have experienced any collision, since repeated collisions of the same particle also count. Consequently

$$
P(p_\perp/m\le q_0)\ge P_{\rm unhit}\ge\max(0,1-C).
$$

At half absolute electrical work the thermal encounter calculation yields $C=0.24979$ for the dedicated model and $0.42807$ for the shared model. These counts require at least **75.02%** and **57.19%**, respectively, in the transverse region. The thermal distributions used to compute the counts place only $2.71\times10^{-10}$ and $9.10\times10^{-10}$ there. The halved wall time step gives 75.02% and 57.19% as well. The inequality needs no assumption about a Poisson distribution of individual collision counts. It is a conditional incompatibility between the exact thermal shape, its own pair rate, and the given injection and streaming. A true nonthermal gas can have a different collision count. [Thermal event budgets](../history/HISTORY.json) · [Rate checks](../history/CHECK.json)

At the first reference hot crossing the thermal count exceeds one, so counting alone gives no positive floor. Continual cold injection and the finite cross-section ceiling still require a transverse fraction at least 56,000 times the dedicated equilibrium fraction and 95,000 times the shared equilibrium fraction. The guaranteed number fractions at this late marker are small.

These inequalities apply to the momentum distribution of this open wall model. They do not fix the error in its heat, pressure, electric current or radius. Fields with transverse electric components, scattering through additional channels, a broader birth spectrum or a different expansion history alter the test.

Reproduce with `python -I -B ballistic/tube.py --traces /path/to/native --history /path/to/HISTORY.json --output /tmp/TUBE.json`. The calculation uses the exact stage times from the thermal event record. Three additional round startup markers use the native absolute electrical-work integral.
