# How many particles can the inherited interaction reach during the pulse?

At the first 1% of total absolute electrical work, the two pressure-fitted source histories still differ from a drifting thermal gas in 46.78% and 55.70% of the counts redistributed among 20 bins. A collision can only change the path of a particle that has experienced one. We bounded the fraction of daughters that could have encountered even one partner by that measurement, allowing previously struck particles to take arbitrary momenta and seed additional encounters.

| Prescribed early/late source and interaction | Potential encounters per daughter if both populations continued freely | Maximum fraction ever struck from the comparison equation | Upper mean encounters per daughter | Remaining 20-bin discrepancy, lower estimate |
|---|---:|---:|---:|---:|
| Pressure fit, dedicated pole | 0.07534 | 23.06% | 0.386 | 23.72 percentage points |
| Pressure fit, shared pole | 0.02908 | 4.89% | 0.0687 | 50.81 points |
| 4% RMS, dedicated pole | 0.07708 | 23.64% | 0.396 | 38.25 points |
| 4% RMS, shared pole | 0.02932 | 4.94% | 0.0694 | 56.78 points |

Numbers are quadrature estimates of one-sided analytic inequalities on the prescribed model. The first column counts potential encounters on freely moving pairs. Actual struck particles may leave the original stream or collide again. Both effects are included conservatively in the next two columns. Each two-body event gives each participant one encounter. The fractions and means refer to daughters already born by the marker and average both charges. [Numerical results](RATE.json)

## The inherited event kernel

Time is measured in wall proper units $L$, with $c=\hbar=1$. A daughter's mass is $m$, $\mu=mL=4\times10^{12}$, the dimensionless shell action is $S_0=1.0415592372880226\times10^{54}$, and the total accepted comoving daughter rest is $D(t)$. One of the three flavor axes, with a specified sign, has wall-frame physical number density $n_{\rm lab}L^3=S_0D/(6\mu a^3)$ on the homogeneous symmetric history. The allowed incoming target has the same flavor and the opposite sign. This density counts the occupants of a common wall-frame volume.

For $p_i$ in units of $mc$ and $E_i=\sqrt{1+|p_i|^2}$, put $s=2(1+E_1E_2-p_1\cdot p_2)$ and $\beta=\sqrt{1-4/s}$. The isolated vector pole has $M/m=8300$ and $g_s=0.30(s+2)$. The dimensionless event kernel including the Møller flux is

$$K=m^2\sigma_{\rm tot} v_{\rm M}=
\frac{6\pi\beta g_s^2}{E_1 E_2[(s-8300^2)^2+\beta^2 g_s^2]}.$$

The dedicated rate uses the factor $b=1$. A single pole with three equal partial widths gives $1/9$ of the dedicated kernel into a particular final flavor and $b=1/3$ when summing all three possible final flavors. The shared rate already counts flavor-changing events. Combining that inclusive rate with a target density summed over all incoming flavors would double count channels. The reported kernel normalization follows the vector pole used in the two histories. The scalar-pole transport calculation in PR #193 uses a different interaction and supplies no additional rate in this component. [Event-rate derivation](https://github.com/maldonado-research/HDblast/blob/6a62a0798713dbe28d0d670376896893cd8d79a1/research/HDBLAST_KINETIC_TEST_20261121/history/RESULTS.md) · [Scalar transport case](https://github.com/maldonado-research/HDblast/pull/193)

The bound $K\le k_*=1.118569\times10^{-6}$ applies for any pair momenta. The inequality uses $E_1E_2\ge s/4$. Below $s=8300^2/2$ the pole denominator is bounded below by $(8300^2/2)^2$. Above that point a rational upper envelope is maximized analytically. The script recomputes the tighter unrounded constant. This cap belongs to the specified pole and need not bound a different interaction or a plasma collective field. [Two-body maximum and cross-section derivation](https://github.com/maldonado-research/HDblast/blob/6a62a0798713dbe28d0d670376896893cd8d79a1/research/HDBLAST_KINETIC_TEST_20261121/ballistic/BALLISTIC.md)

## A bound that allows subsequent collisions to spread

Let $T$ be the measurement time, $z(t)=D(t)/D(T)$ the born fraction, and $A(t)=b S_0D(T)/(6\mu^3 a(t)^3)$. For each sign let $F^0_\pm(t,p)$ denote the collisionless momentum counting measure normalized to the sign's final population. Its mass at time $t$ is $z$. Define an integrated first source and a maximum rate by

$$R_0(t)=A(t)\int K F^0_+(dp)F^0_-(dq),\qquad \lambda_*(t)=A(t)k_*z(t).$$

Use a spatially homogeneous classical binary Boltzmann gas with the prescribed births and fields. Tag a particle at its first encounter. Before being hit each daughter follows the collisionless characteristic of its birth cohort. Therefore the measure of untagged particles of either sign is dominated by $F^0_\pm$ at every time. Let $h_+,h_-$ be the tagged masses normalized per sign at $T$, and $h=(h_++h_-)/2$. An untagged pair that encounters each other adds two tags. An encounter of a tagged and untagged particle adds one tag. Encounters between two tagged particles add none. The rate of the first class is bounded by $R_0$, using positivity of the collision kernel and domination of both measures. The rate contribution to $h$ from the mixed class is bounded by $A k_* z(h_++h_-)/2=\lambda_*h$. This gives

$$
 h'(t)\le R_0(t)+\lambda_*(t) h(t),\quad
 h(T)\le Y(T)=\int_0^T R_0(s)\exp\!\left[\int_s^T\lambda_*(v)dv\right]ds.
$$

This first-hit argument also works after averaging over rotational copies of three equally supplied flavors under flavor conversion. Each untouched flavor has mass at most the collisionless mass per flavor. Summing mixed-tag encounter contributions then requires only the total tagged mass. A statement about one particular flavor after conversion requires the assumed flavor symmetry. The homogeneous kinetic treatment assumes molecular chaos and the given densities. It does not bound concentration of targets in clumps.

With the symmetric total densities held fixed, the rate of events of *all* classes per normalized sign is bounded by $R_0+2\lambda_*h$: charge every event involving at least one tag against the density of tagged positive or tagged negative participants and permit every partner with momentum-independent rate $k_*$. A two-tag event can be counted twice in this upper bound. Integrating with $Y'=R_0+\lambda_*Y$ gives a mean encounter cap $2Y(T)-\int_0^T R_0dt$. A simpler, independent upper bound on the number of encounters for a particle born at $t_b$ is $\int_{t_b}^T\lambda_*dt$. The age-weighted mean of this simpler cap is 1.258 and 0.611 on the two fine-time histories. The corresponding upper fractions that ever encounter any partner are 66.24% and 43.53%, from the inherited local collision hazard. [Numerical comparison equations](exposure.py)

Any assignment of new momenta to a fraction $h$ of particles can change a histogram's total variation distance by at most $h$. We use the previous 20-cell distances on the same prescribed fields with both transverse currents switched off. The table's last column subtracts $Y(T)$ from those distances. As a sensitivity test, multiplying the entire event kernel by $c$ multiplies both $R_0$ and $\lambda_*$ by $c$. Reaching a discrepancy below five percentage points would require at least $c\simeq 1.275$ and $3.023$ for the pressure fits, or $1.412$ and $3.125$ for the 4% RMS pair. These are necessary multipliers on frozen histories. The prescribed pole has $c=1$. A new mediator would generally change the shape of the event kernel as well as the histories.

## Evaluation and controls

The homogeneous two-polarization fields and continuous angular birth acceptance come from the versioned predecessor. In the invariant basis on the round wall the free trajectory's comoving momentum has constant components. We redshift with the evolving radius and advance the electric force and magnetic rotations with the inherited charge-to-mass normalization. A six-dimensional scrambled Sobol net covers both particles' birth ages, polar angles and azimuths jointly. Its first and last three coordinates are independent in the continuum quadrature. Only pairs for which both daughters have been born enter $R_0$. Integrating over the marker's full birth population automatically includes the factor $z^2$. The program checks the sampled live-pair mass against $z^2$, the pointwise kernel against its analytic maximum and the forward comparison equation against its integrated exponential solution. Exact controls include resting, parallel and head-on resonant pairs.

Each finer history uses 8,192, 16,384 and 32,768 independent pair quadrature nodes with different scrambles, and 192, 256 and 384 time intervals, respectively. A separate run on each original time resolution uses 16,384 nodes. For the fitted cases the reported first-hit caps differ by less than 0.004 percentage point across these settings. The tests measure numerical stability without certifying continuum quadrature error. They apply before the first 1% of absolute electrical work on spatially homogeneous prescribed histories. Medium corrections to the microscopic kernel, changes to the mean fields and scattering from other species require a new rate calculation.

Use the unpacked #198 study as the parent and its regenerated primary field traces.

~~~sh
python -I -B rate/exposure.py --parent /path/to/HDBLAST_POLARIZATION_20261121 \
  --fields /path/to/polarization-results/fields/traces --output /tmp/RATE.json --jobs 2
~~~

[Pair integrator](pair_clock.cpp) · [Floating-point inputs and results](RATE.json)

## Direction-independent enclosure using birth times

A second calculation encloses $R_0$ using only the fields, the maximum birth momentum, and when the particles form. It allows **every** birth momentum direction. Across 768 time and birth intervals the first-hit fraction in the pressure fits is bounded above by **31.42% and 6.873%**. Combining these upper bounds with the measured collisionless histogram distances leaves at least **15.36 and 48.82 percentage points**. The corresponding 4% source estimates are **31.82 and 54.94 points**. This bound on the encountered fraction avoids angular sampling. Its inputs are the piecewise linearly interpolated field and birth histories. The archived particle-count distances have finite sampling error. [Finite-cell enclosures at multiple resolutions](ENCLOSURE.json)

For the prescribed homogeneous modes, let $s_i=\pm1$ be a charge sign. The comoving momentum obeys $d(a\mathbf p)/dt=s_i Q a(\mathbf E+\mathbf v\times\mathbf B)$ in the inherited invariant frame. The axial electric field stays negative and exceeds the horizontal magnetic field in magnitude over the recorded birth interval. Hence $d[-s_i a p_z]/dt\ge Q a(-E_z-|B_\perp|)>0$. Also $d|a\mathbf p|/dt\le Q a|\mathbf E|$, since magnetic rotation does no work. Bounding the integrands on every original field interval gives positive cumulative lower and upper integrals. The extreme birth and observation times within each new integration cell then give a lower bound on forward $-s_i p_z$ and an upper bound on total momentum. The injected momentum magnitude never exceeds $\sqrt{1.05^2-1}\,mc$.

If both particles have guaranteed forward magnitudes $Z_i>0$ with opposite physical directions along the axis, their invariant $s$ satisfies

$$s\ge2\left[1+\sqrt{1+Z_1^2}\sqrt{1+Z_2^2}+Z_1Z_2\right].$$

The transverse inequality follows from $E_1E_2-|p_{\perp1}||p_{\perp2}|\ge\sqrt{1+p_{z1}^2}\sqrt{1+p_{z2}^2}$. A newborn cohort with an unresolved sign instead uses $s\ge4$. Upper bounds on the two momentum magnitudes give $s\le2(1+E_1^{\max}E_2^{\max}+p_1^{\max}p_2^{\max})$. The pole factor divided by $6\pi/(E_1E_2)$ is $\beta/(y^2+\beta^2)$, with $y=(s-M^2)/[0.30(s+2)]$. Monotonicity of $y$ gives its minimum absolute value over the allowed interval. Maximizing the remaining rational function over the allowed $\beta$ interval gives a pairwise upper bound. Below $s=M^2/2$ the pole factor is less than one. Above that point the bound $1/\beta\le(1-8/M^2)^{-1/2}$ also applies. The code takes the minimum of these bounds and the global kernel cap.

For each time cell we use the upper bound for every pair of birth cells, charge their entire rest increments including any newborns in the observation cell, use the smallest cell radius and the largest cell population, then integrate the positive comparison equation. At 512 cells the separate-interaction fitted upper fraction is about 35.54%. At 768 it is 31.42%. Finer numerical angular quadrature estimates the bound's comparison integral more tightly in the table above. The finite-cell calculation deals conservatively with birth and observation intervals and uses ordinary floating point with small outward margins. It does not certify the continuum gravitational evolution.

~~~sh
python -I -B rate/enclosure.py --parent /path/to/HDBLAST_POLARIZATION_20261121 \
  --fields /path/to/polarization-results/fields/traces --output /tmp/ENCLOSURE.json --jobs 2
~~~
