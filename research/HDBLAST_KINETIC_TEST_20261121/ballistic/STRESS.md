# Do enough particles scatter to carry the fluid current and pressure?

At half of the total absolute electric work in the dedicated-pole history, **at least 45.33% of all particles born so far must have scattered** to account for the adopted fluid current and its pressure along a charged stream. The thermal collision integral gives **0.2498 collision incidences per particle**, enough to touch at most 24.98%. The joint current and pressure requirement therefore exceeds the available incidence count by a factor of at least **1.81**. The calculation gives the model every advantage by allowing collisions to select whichever particles need redirection most, and by charging only one incidence per affected particle. [Results](STRESS.json) · [Thermal collision counts](../history/HISTORY.json)

The shared vector model gives a different outcome. At its half-work stage at least 31.46% must have scattered by this test, while the encounter budget reaches 0.4281 incidences per particle. This leaves the shared model compatible with the **necessary moment bounds** at that stage. Its full momentum distribution still fails the transverse-distribution test. [Transverse number-distribution test](TUBE.md)

| Dedicated half-work marker | Native wall step | Half wall step |
|---|---:|---:|
| Thermal collision incidences per particle | 0.24979 | 0.24979 |
| Minimum collision participants to match both moments | 45.33% | 45.345% |
| Lowest attainable mean axial velocity with this collision count | $-0.47413c$ | $-0.47412c$ |
| Reference fluid velocity, including adjacent native endpoints | $[-0.54447,-0.54396]c$ | $[-0.54441,-0.54415]c$ |
| Minimum longitudinal stress from unhit particles, in units $n_{\rm lab}m$ | 24,891 | 24,894 |
| Upper reference fluid longitudinal stress in the same units | 21,920 | 21,916 |

The sign is relative to the positively charged stream's field axis. The negative stream is its reflection. The electrical current depends on the number-weighted axial velocity through $J=QD\langle v_z\rangle/a^3$. The stresses shown here refer to a single stream's longitudinal component. Summing the three equal perpendicular axes makes the total spatial stress isotropic.

## Bounding every possible assignment of collisions

Use the positive-stream characteristics $p_z/m=[I(t)-I(t_b)+a_b p_b\xi]/a(t)$ and $p_\perp/m\le q_0=0.320157$ throughout the expanding interval. All quantities below use $m=1$ and $c=1$. The cohort born in step $i$ has number weight $w_i=\Delta D_i/D(t)$. Bound its birth impulse $I(t_b)$ using both native endpoints. For $I'=QaX$ the largest absolute endpoint force times the cell duration is $J_i$. With linear interpolation of the native force it encloses the within-cell impulse change. Intersect $[I_i-J_i,I_i+J_i]$ and $[I_{i+1}-J_i,I_{i+1}+J_i]$. Together with the two extreme birth directions this gives a longitudinal interval $[z_-,z_+]$ for the cohort. Using $a_b/a\le1$ to bound its birth impulse enlarges this interval.

The function $v_z=z/\sqrt{1+z^2+p_\perp^2}$ increases with $z$. Minimize it at $z_-$ and maximize at $z_+$, allowing $0\le p_\perp\le q_0$. Call the resulting bounds $v_i^-,v_i^+$. Writing $z_0=\min_{z\in[z_-,z_+]}|z|$, the cohort's unhit longitudinal stress per particle is at least
$$
p_i^- = \frac{z_0^2}{\sqrt{1+z_0^2+q_0^2}}.
$$

Let $C$ be the mean total collision incidences per currently present particle. The fraction still unhit is at least $F=\max(0,1-C)$. For any nonnegative cohort cost $g_i$, define
$$
L_F(g)=\min_{0\le x_i\le w_i,\ \sum_i x_i\ge F}\sum_i x_ig_i .
$$
Sorting the cohorts by cost and keeping the cheapest $F$ of their number weight solves this linear program, including a fractional final cohort. The reconstructed trajectory and supplied incidence budget then obey
$$
-1+L_F(1+v^-)\ \le\ \langle v_z\rangle\ \le\ 1-L_F(1-v^+),\qquad
\frac{P_{zz}}{n_{\rm lab}m}\ge L_F(p^-).
$$
Every particle outside the unhit set may travel at either limiting speed in the current bound. Its stress is nonnegative. The test also bounds the maximum two-direction transverse pressure if the candidate total energy is retained.

For the **joint** test, keep the same unhit cohort allocation for both moments. Fix $\ell$ to the upper reference mean lab energy in units of $m$. For the dedicated stream at the marker, the mean velocity is negative. All particles contribute a nonnegative amount to
$$
G=\frac{P_{zz}}{n_{\rm lab}m}+\ell(1+\langle v_z\rangle).
$$
Its lower bound is $L_F(p^-+\ell(1+v^-))$. For a positive reference mean velocity we use $1-\langle v_z\rangle$ and the upper velocity bound. Compare with the **largest** $G$ allowed by the fluid's bounding native endpoints. Inverting the sorted cumulative costs gives the largest permitted unhit fraction. We discarded the whole positive contribution of every hit particle. A full collision solution may therefore require more participants than this lower bound.

With the fluid's Taub–Mathews specific enthalpy $h=5\Theta/2+\sqrt{1+9\Theta^2/4}$, $u=\gamma v$ and $\Theta=T/m$, the reference mean energy is $h\gamma-\Theta/\gamma$ and the reference axial stress is $(hu^2+\Theta)/\gamma$. The two-direction transverse stress is $2\Theta/\gamma$. The code forms separate endpoint bounds on $|u|$, $\Theta$, and signed velocity to allow for interpolation within the marker cell.

## Resolution and scope

The largest impulse any single input cell contributes at the dedicated native marker is 50.27 in physical $p_z/m$ units, falling to 25.14 on the half step. An additional sensitivity test widens **every** birth cell by $\pm500$ in these units and increases the entire collision budget by 1%. It still needs at least 44.86% of particles to scatter, retains a velocity gap above $0.058c$ and gives a longitudinal stress floor more than 10.8% above the fluid value. These perturbations are numerical sensitivity tests. The native and halved histories are numerical solutions of the parent model.

At 10% of the absolute electric work the respective current and pressure tests do not exclude either benchmark. At 90% and at the reference hot crossing the total incidence budgets exceed one per particle; a pure counting argument then permits every particle to have scattered. The mismatch found here is a specific failure of the dedicated-pole fluid closure early on its published history. An actual kinetic solution can change its collision rate, current, or electric field. Determining the resulting heat and radius requires evolving those quantities together.

Run `python -I -B ballistic/stress.py --traces /path/to/native --history /path/to/HISTORY.json --output /tmp/STRESS.json`.
