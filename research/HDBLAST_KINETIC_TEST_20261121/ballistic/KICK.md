# How much transverse momentum would an extra force have to supply?

At the point when the reference wall has exchanged 1% of its whole-run absolute electrical work, some daughters have had too little time to scatter through the selected elastic interaction. Suppose their transverse momentum distribution must already match the thermal distribution assumed in the earlier fluid calculation. In the dedicated-channel example, an additional force must then supply an RMS accumulated transverse impulse of **at least $541\,mc$**, averaged over all daughters already born. The shared-channel example requires **at least $674\,mc$**. Here $m$ is the daughter mass. These are optimistic requirements on a proposed mechanism acting on the specified histories. [Earlier transverse-momentum bounds](TUBE.md) · [All numerical inputs and results](KICK.json)

## From a collision probability to an impulse

Use units $c=1$ and let $b=0.32016$ bound $p_\perp/m$ for any particle that has experienced only expansion and the reference axial fields since birth. The parallel electric field adds longitudinal momentum. The parallel magnetic field rotates the transverse momentum, and expansion reduces its magnitude. Now allow another force with physical transverse component $\mathbf F_{\rm extra,\perp}$. For a daughter born at $t_b$ that has avoided the specified elastic interaction through $t$, define

$$
J=\frac{1}{a(t)m}\int_{t_b}^t a(t')|\mathbf F_{\rm extra,\perp}(t')|\,dt',
\qquad J\ge(p_\perp/m-b)_+ .
$$

The absolute value inside the integral grants the new force perfect directional efficiency. The original pole kernel has a global upper bound for any incident momenta. The fixed birth, density and expansion histories and a Markov binary collision law therefore imply an unhit probability of at least $f_0$ for an already-born daughter, even if target momenta vary at fixed density. This statement concerns the population distribution or expected fractions. An individual finite Monte Carlo sample can fluctuate. [Derivation of the any-target bound](BALLISTIC.md)

Let $F$ be the cumulative distribution of $p_\perp/m$ in the target stream at the marker. To minimize the required sum of squared impulses, give the unhit population exactly probability $f_0$, assign it the *lowest* eventual transverse momenta, and grant everyone else free redistribution. Since $(p_\perp/m-b)_+^2$ increases with $p_\perp$, every allowed allocation obeys the bound below.

$$
\sqrt{\mathbb E[J^2]}\ \ge\left[\int_0^{f_0}(F^{-1}(u)-b)_+^2\,du\right]^{1/2}.
$$

The expectation uses **all daughters present** in its denominator. The bound credits only those avoiding the original interaction. At least $f_0/2$ of all daughters must receive $J\ge(F^{-1}(f_0/2)-b)_+$, and the essential supremum of $J$ is at least $(F^{-1}(f_0)-b)_+$. This is an impulse bound. Inferring supplied energy requires the force direction and particle velocity along each path.

## The thermal transverse distribution

Write $p=p_\perp/m$, $w=\sqrt{1+p^2}$ and $\Theta=T/m$. Integrating the isotropic Maxwell–Jüttner momentum density over the longitudinal component gives

$$
F'(p)=\frac{p w K_1(w/\Theta)}{\Theta K_2(1/\Theta)},
\qquad 1-F(p)=\frac{w^2 K_2(w/\Theta)}{K_2(1/\Theta)} .
$$

The second identity follows from $d[x^2 K_2(x)]/dx=-x^2 K_1(x)$. A boost along the field leaves $p_\perp$ fixed. The equal-wall-time measure has a Jacobian proportional to $1+v p_z/E$. Its odd term integrates to zero at each $p_\perp$ for this isotropic rest distribution. Consequently the normalized transverse marginal above also applies in the wall frame.

| At 1% of the whole-run absolute work | Dedicated pole | Shared pole |
|---|---:|---:|
| Minimum probability of avoiding the pole interaction | 0.33724 | 0.56453 |
| Minimum RMS $J$ over all daughters already born | 540.96 | 673.71 |
| Fraction that must acquire at least the impulse in the next row | 0.16862 | 0.28226 |
| Corresponding minimum individual $J$ | 895.54 | 831.34 |
| Minimum essential supremum of $J$ | 1416.01 | 1445.63 |
| Minimum peak $L|F_{\rm extra,\perp}|/m$ | 52,604 | 45,384 |
| Equivalent lower bound on peak transverse electric field in wall units | 0.005260 | 0.004538 |

The peak-force estimate grants even the last-born particles the full elapsed time since the left endpoint of the first numerical birth interval. Because $a(t')/a(t)\le1$ during expansion, $J\le \sup(L|F_{\rm extra,\perp}|/m)\,\Delta t/L$. For a purely electric force the dimensionless acceleration is $Q|X_\perp|$ with $Q=10^7$. The original *parallel* field reached about 0.207 in the sampled interval. Its orientation is essential. A transverse magnetic force could instead redirect existing particle momentum while doing no work. Maintaining the same density history, collision kernel and target temperature as a new force acts would have to be tested in a coupled evolution.

At halved integration steps the RMS lower bounds become 541.65 and 674.10. Direct density integration verifies both target quantiles. Independently evaluating the truncated squared moment as $2\int_b^{F^{-1}(f_0)}(p-b)[f_0-F(p)]\,dp$ reproduces the result. The code uses scaled Bessel functions and verifies normalization. [Reproduction code](kick.py)

From the capsule root run this after generating the common native histories and the transverse-momentum bounds.

```sh
python -I -B ballistic/kick.py --tube ballistic/TUBE.json --native /tmp/hdblast-native --output /tmp/KICK.json
```
