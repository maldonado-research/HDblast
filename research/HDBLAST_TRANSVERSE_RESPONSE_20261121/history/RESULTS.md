# How much sideways force do newly born particles need?

A particle created just before the measurement has little time to turn. Even a weak continuing source of cold particles can therefore prevent *exact* thermal equilibrium, while leaving a gas that is close to equilibrium at many measured energies. We quantified that distinction on the two early HDblast reference histories. The thermal comparator is the transverse momentum distribution expected at the recorded fluid temperature.

At 1% of whole-run absolute electrical work, the largest reference axial electric field since birth is 0.2071 in wall units. Let the **equivalent field** of an extra force be its peak transverse force per unit charge in these same units. The following are necessary lower bounds, allowing the new force its most favorable direction. [Focused results](FOCUS.json) · [Original impulse test](https://github.com/maldonado-research/HDblast/blob/6a62a0798713dbe28d0d670376896893cd8d79a1/research/HDBLAST_KINETIC_TEST_20261121/ballistic/KICK.md)

| Required agreement at every transverse momentum threshold | Dedicated interaction | Shared interaction |
|---|---:|---:|
| Within 1 percentage point | 0.0391 | 0.0219 |
| Within 0.1 percentage point | 0.126 | 0.0691 |
| Within 0.01 percentage point | 0.395 | 0.215 |
| Exact thermal marginal | 11.98 | 4.38 |

For exact agreement the force requirement is about **58 and 21 times** the original peak axial electric force. The exclusion at the strongest cutoff involves only about **1.13 in ten million** and **2.52 in ten million** of the population, respectively. Requiring agreement within one percentage point gives much smaller targets, about 19% and 11% of that axial peak.

Suppose the extra force per charge can reach 10% of the reference axial peak. At some transverse momentum threshold, the particle count must then differ from the thermal count by at least **3.12** and **1.11 percentage points**, respectively. Even with a force equal to the full original axial peak there remains a minimum difference of **0.0369** and **0.0108 percentage points**. Under that latter cap, an alternative way to attempt exact agreement would be to broaden the birth distribution. The allowed transverse birth momentum would have to extend to at least $18.57\,mc$ and $6.78\,mc$, respectively. For a cold parent at rest decaying to two equal-mass daughters, this requires accessible parent masses $M/(2m)\ge18.60$ and $6.86$, compared with $1.05$ in the benchmark. These are kinematic trade-offs on the specified histories. Changing the masses, birth rates, plasma density or electric pulse would require a new wall calculation. [Trade-offs and witnesses](FOCUS.json)

## Counting the youngest particles

Use wall time in units of $L$, with $c=1$. Write $D(t)$ for cumulative accepted comoving rest in all daughters and $a(t)$ for wall radius. Inherited birth kinematics require $p_{\perp,b}/m\le b=\sqrt{1.05^2-1}$. Let $G$ be an upper bound on the magnitude of an extra transverse force divided by $m/L$. During expansion a particle born at $s$ that retains its original flavor and avoids the selected elastic interaction must have

$$
p_\perp(t)/m\le b+G A(s,t),\qquad A(s,t)=\int_s^t a(v)\,dv/a(t).
$$

The original axial magnetic field can rotate the transverse momentum. It cannot enlarge its magnitude. For electromagnetic forces, $G/Q$ limits the equivalent field with $Q=10^7$ in the dimensionless wall conventions. A bound on $|\delta E_\perp|+|\delta B_\perp|$ is one sufficient bound on the transverse deflecting field per unit charge. A spatially localized source requires a bound on what tagged particles encounter along their trajectories.

Let $s$ be any cutoff birth time and $T$ the test time. Define $S(s)$ as a lower bound on the *number fraction among all daughters present at $T$* born after $s$ and still unscattered. On the recorded monotonically expanding, continuously producing histories, use the number density of possible opposite-sign targets and the global upper bound $k_*=1.118569\times10^{-6}$ (rounded upward here) on $m^2\sigma v_{\rm M}$. The inherited channel fraction is $b_{\rm ch}=1$ or $1/3$. For any incident and target momenta each cohort’s classical binary collision hazard is bounded by

$$
\Gamma_{\max}(v)=\frac{b_{\rm ch} S_0 D(v)}{6(mL)^3 a(v)^3}\,k_*.
$$

Each birth increment contributes at least $dD(u)\exp[-\int_u^T \Gamma_{\max}(v)dv]/D(T)$ to $S(s)$ when $u\ge s$. This bound uses the given density history along the trajectories and the molecular-chaos binary collision law. It permits an arbitrary distribution of target momenta. Spatial concentrations that raise the target density, new collisions or a new source can change the hazard.

Let $F_\Theta$ be the cumulative fraction below $p_\perp/m$ for a Maxwell–Jüttner distribution boosted along the field. A longitudinal boost leaves this transverse number marginal invariant. Its exact tail at $q=p_\perp/m$ is $(1+q^2)K_2(\sqrt{1+q^2}/\Theta)/K_2(1/\Theta)$. If the largest CDF discrepancy is to be at most $\epsilon$, then for *every* cutoff $s$ one needs $F_\Theta(b+GA(s,T))\ge S(s)-\epsilon$. This implies

$$
G\ge \max_s \frac{[F_\Theta^{-1}((S(s)-\epsilon)_+)-b]_+}{A(s,T)}.
$$

The code also evaluates the forward discrepancy $\max_s [S(s)-F_\Theta(b+GA(s,T))]_+$ and the required birth cap $\max_s [F_\Theta^{-1}(S(s))-GA(s,T)]_+$. Each discrepancy bounds the total variation distance from below. We compare a positive density integral to the Bessel CDF at each archived broad-window witness. [Calculation](budget.py)

## Time resolution and a separate short-time check

The four compact histories retain every original step from the final zero-population knot to the first step at or beyond 5% of whole-run absolute electrical work. They include the original 14 columns, the reconstructed guide field $B=2\int_0^t aX\,dt/a^2$ with $C(0)=0$, and the axial impulse $I=Q\int aX\,dt$ initialized to zero at the birth window’s left knot. Fractional work is measured from $\int |QD X u/\sqrt{1+u^2}|dt$ divided by the predecessor’s whole-run ledger. [Extraction](extract.py) · [Clocks and parameters](markers.json)

For each fully completed native cell the sum counts its accepted births and puts them at the earliest cell time for the hazard calculation. The density bound uses the largest recorded rest and smallest radius in each cell. We discard all births in the terminal bracketing cell, allow survival and impulse to be integrated to its right endpoint, and take the smaller endpoint temperature. The thermal transverse distribution stochastically broadens with temperature: the rest-frame energetic likelihood ratio increases with energy, and the uniform angular law is independent of that energy. The smaller temperature therefore gives a conservative required momentum among the bracketing temperatures. These are conservative enclosures under the smooth monotone history/interpolation used by the wall model. They do not constitute a verified continuum integration of the five-dimensional equations. [Broad-window results](BUDGET.json)

The original time steps resolve the one-percentage-point bound reasonably well, giving equivalent fields 0.0369/0.0209 on the native grid and 0.0379/0.0214 with the time step halved. Exact matching probes exceedingly recent births. We therefore reran the coupled source equations and shortened their time step within $0.0015L$ of the 1% marker. The adaptive cap is the larger of a minimum step and the remaining time divided by 64, with one step ending at the marker. We vary the minimum from $10^{-8}L$ to $2\times10^{-11}L$, repeat on the native and half-time background trajectories, and shift each finer background's marker by $\pm2\times10^{-5}L$. These focused runs count only births within the refined window. At the finest resolution the main native/half exact bounds are 12.002/11.981 and 4.384/4.380, respectively. The shifted markers retain bounds above 11.93 and 4.37. The maximizing recent age is about $2.7\times10^{-9}L$ in the dedicated case and $7.3\times10^{-9}L$ in the shared case. The thermal cutoffs are both near $0.64\,mc$. [Focused solver](focus.py) · [Resolution and timing tests](FOCUS.json)

A separate local estimate explains these scales. Write $\lambda=\dot D(T)/D(T)$. Over a very short age $h$ during smooth injection, $S(T-h)\simeq\lambda h$, $A\simeq h$, and in a hot thermal gas $F_\Theta(q)\simeq q^2/(4\Theta^2)$ for $q\ll\Theta$. Optimizing $(2\Theta\sqrt{\lambda h}-b)/h$ gives $h\simeq b^2/(\Theta^2\lambda)$ and $G\simeq\Theta^2\lambda/b$. The measured values $\lambda L\simeq41.98,34.59$ and $\Theta\simeq956.75,637.09$ give equivalent fields 12.003 and 4.385. This estimate uses smooth instantaneous injection. Microscopic intermittency, a formation-time model, collective effects on the birth process or changes to the incoming density are physically different assumptions, especially at such short ages.

From the study root, `python -I -B history/budget.py --history history --output /tmp/BUDGET.json` uses the four compact tables. The full study driver regenerates the wall histories and runs the focused source calculation with the supplied predecessor inputs.
