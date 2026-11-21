# Can the weaker local sideways pressure reach its thermal value?

Oppositely moving thermal streams have two equal lateral pressures. Particles at different places on the shell can move along different sideways directions and still be nearly one-dimensional within each small neighborhood. This calculation bounds the smaller horizontal pressure **at the particles' final positions** and averages it over particles in a specified region.

Take both spatial quadratures for each of the five angular modes 2, 5, 7, 40 and 258 and all seven tested source waveforms. This gives 70 real freely chosen source coefficients. With a 4% budget in launch amplitude, in the half of the wall satisfying $|z_1|^2\le 1/2$, the average smaller lateral pressure per particle is at most **0.62% or 2.64% of the moving-gas prediction** on the finer time histories. The original time steps give 0.62% and 2.65% after rounding upward. Doubling the radial mesh for the first interaction history gives 0.67%. This regional calculation requires uniform initial positions within each birth cohort. It allows subsequent clustering and correlations of fields with particle velocities. [All regional and whole-wall bounds](MIXTURES.json)

| RMS launch amplitude relative to original emission | Maximum weak pressure as a fraction of thermal in half the wall, separate interactions | Shared interaction |
|---|---:|---:|
| 2% | 0.00198 | 0.00859 |
| 4% | 0.00617 | 0.02637 |
| 8% | 0.02162 | 0.09170 |
| 16% | 0.07975 | 0.33611 |

The entries are conservative upward roundings of the finer-time results. At 4% amplitude the added field's RMS-square budget is 0.16% of the original emitted characteristic-square integral with the specified tube-flux normalization. The source coefficients may have either sign. They are charged for their full waveform overlap. The strongest new short-scale field can be large near $|z_1|=1$. Its highest-weight harmonic falls rapidly away from that part of the sphere. Each cohort can travel less than 0.0031 comoving radians by the marker. The 70-coefficient bound over the entire sphere, using the largest possible field at every position, permits substantial pressure from this short-scale mode and is inconclusive for the whole-wall average. These regional bounds apply to the ten selected spatial functions with a shared angular coordinate. Sources in other angular patterns or distributed around different angular centers require new tests.

For a second subset, the original pulse feeds one spatial mode and an independent timing is allowed in a second mode. The two-mode calculation covers every split of a specified amplitude budget. When both modes lie among 2, 5 and 7 and the total amplitude is 4%, the weaker lateral pressure averaged over all final particle positions is at most **6.4% or 26.6% of its thermal target** across all seven choices of the second timing. This whole-sphere result needs no spatial uniformity assumption at birth. [Two-mode results](WEAK.json) · [Two-mode calculation](weak_pressure.py)

## Why a global pressure can look thermal

In the chosen left-invariant frame write $S(x)$ for the local two-by-two matrix of horizontal pressure per particle and $q(x)=S_{11}-S_{22}+2iS_{12}$. With averages over the final particle-number measure,
$$
\lambda_{\min}(\langle S\rangle)-\langle\lambda_{\min}(S)\rangle
={\langle|q|\rangle-|\langle q\rangle|\over2}\ge0 .
$$
For a simple example, let the local stress be entirely along the horizontal vector proportional to $\bar z_1^{\kappa-2}$ and let its magnitude depend only on $|z_1|$. At every point the weaker eigenvalue vanishes. For $\kappa>2$, uniform integration over the phase of $z_1$ gives two equal diagonal entries in the averaged matrix. This example uses the same directions as the waves and explains how the orientation of the averaging axes affects the interpretation.

## Combining arbitrary source coefficients

Normalize the applied waveforms by their energy at the outer emitting surface. Their measured overlap matrix $G$ has unit diagonal and is positive definite. In one angular quadrature, write the coefficients as $b$ and use $G=LL^T$, $c=L^Tb$. The emitted perturbation budget summed over both quadratures and all angular degrees is then $\sum c^2\le\epsilon^2$. Use $L^{-1}E$ and $L^{-1}B$ for the corresponding temporal transfer vectors. This conversion retains destructive and constructive interference between waveforms.

For each angular degree $\kappa$, write $Y_{\kappa j}$ for the two perpendicular quadratures, $j=1,2$, of the same complex highest-weight angular function. In the tested frame the horizontal complex components are proportional to $\sqrt{\kappa-1}\,\bar z_1^{\kappa-2}$ and its right-angle rotation. During any paths ending in the chosen region, let $M_\kappa$ bound their common magnitude. For the whole sphere $M_\kappa=\sqrt{\kappa-1}$. For a region with $|z_1|\le R$, take $R_+=\min(1,R+\chi)$, $\chi=\int_{t_0}^T dt/a(t)$, and $M_\kappa=\sqrt{\kappa-1}\,R_+^{\kappa-2}$. The real frame components of each pair of quadratures can change by at most
$$
\Delta_\kappa(t)=\min\left(2M_\kappa,\sqrt{\kappa-1}(\kappa-2)R_+^{\kappa-3}\int_t^T {du\over a(u)}\right).
$$
This follows by factoring the difference of two complex powers and bounding the chord displacement by the path length. Both quadratures share the same complex difference, up to rotation.

Use a combined degree and temporal-profile index $m$. At cohort birth time $s$ define $I_m(s)=M_\kappa\int_s^T a(t)(L^{-1}E)_m(t)dt$. Choose a single unit coefficient vector $r$ for all cohorts within this test. At each observed point $x$ choose a horizontal unit $n(x)$ perpendicular to $\sum_{m,j} c_{mj}r_mY_{\kappa j}(x)/M_\kappa$. At zeros any horizontal unit is admissible. The coefficients $v_m=\sum_j c_{mj}\,n\cdot Y_{\kappa j}(x)/M_\kappa$ obey $v\cdot r=0$ and $\|v\|_2\le\epsilon$, by orthogonality of the two quadratures at the same angular degree. The static-position electric impulse along this direction is therefore bounded by $\epsilon\|I(s)-r(r\cdot I(s))\|_2$ times $Q/a(T)$. A leading eigenvector of the measured age-integral covariance supplies one choice of $r$. The inequality holds for every choice. The two-mode calculation also tests trial directions directly.

Changing position contributes a further upper impulse of
$$
{Q\epsilon\over a(T)}\int_s^T a(t)\left[\sum_m\Delta_\kappa(t)^2(L^{-1}E)_m(t)^2\right]^{1/2}dt .
$$
The projection uses left-invariant frame components along the path. Comoving momentum components are constant along a free geodesic in this frame. The tested direction is evaluated once at the final position and can vary across observation points. Its spatial derivatives never enter the bound.

## Magnetic forces and relativistic pressure

Write $H_E,H_B$ for upper bounds on the magnitudes of the actual horizontal electric and magnetic field. For either field take $\epsilon$ times the Euclidean norm of its $M_\kappa$-weighted whitened transfer vector. Let $X$ and $B_3$ be the given axial electric and magnetic fields and $p_*$ the maximum birth momentum in units of $mc$. The comoving momentum evolves by the Lorentz force in the left-invariant frame. The horizontal magnetic force is at most $H_B+|B_3|v_h$, where $v_h$ bounds the sideways speed in units of $c$. The magnetic guide preserves the horizontal momentum magnitude. The axial electric force has fixed sign for each charge. The horizontal fields can cancel that force by at most $H_B$. Thus at the marker a tagged particle has
$$
\gamma(T)\ge\sqrt{1+\left[{Q\over a(T)}\int_s^T a(t)(|X|-H_B)dt-p_*\right]_+^2}.
$$
Using a speed bound of one always works. When $|X|>H_B$ throughout, define $r_v=\max_t(H_E+H_B)/(|X|-H_B)$. The horizontal momentum at any intermediate time is at most $p_*+r_v g$, and the magnitude of the axial momentum is at least $[g-p_*]_+$, with $g$ the accumulated lower axial impulse at that time. Maximizing the resulting speed bound over $g\ge0$ gives $v_h\le h/\sqrt{1+h^2}$ with $h^2=r_v^2+p_*^2(1+r_v)^2$.

Add $p_*$, the bounds on the projected electric impulses and the integrated magnetic force. Square this sum and divide by the lower bound on $\gamma$. This bounds $p_n^2/\gamma$ for each cohort. Average with its birth weight $dD(s)/D(T)$. Conservative endpoint rectangles overestimate all positive forces, overestimate travel, and underestimate the signed axial impulse. The signed electric transfer uses exact products of piecewise-linear fields and radius. The full-space response test adds an upper force for the entire uncertain birth cell. It places every birth anywhere in its native cell. [Implementation](mixtures.py)

For a uniform cohort on $S^3$, $|z_1|^2$ is uniform on $[0,1]$. All particles born with $|z_1|\le R_-=\max(0,R-\chi)$ stay inside $|z_1|\le R$ until observation. Every particle ending inside was born with $|z_1|\le R_+$. Consequently the numerator of the conditional pressure is bounded using birth fraction $R_+^2$, while its denominator contains at least $R_-^2$ of the particle count. The region-free test uses factor one. The moving thermal reference has horizontal pressure per particle $\theta/\Gamma$ in either direction, where $\theta$ is its dimensionless proper temperature and $\Gamma$ its axial Lorentz factor. This comparison applies to the specified paired stream.

A regional average need not describe every point inside the region. The JSON also gives lower bounds on the fraction of **all particles** located where the smaller local pressure falls below a specified fraction $f$ of its benchmark. For a regional conditional pressure cap $C_R$ in units of thermal, that fraction is at least $R_-^2[1-C_R/f]_+$ by nonnegativity. The code chooses the largest of four such lower bounds. At a 4% launch, more than 70% and 60% of particles in the two cases must be at sites whose smaller local horizontal pressure per particle is below 90% of the moving thermal value, across the tested resolutions.

These inequalities use the fixed round background, uniform axial guide and original birth history, linear radial transfer of prescribed transverse signals, and collisionless trajectories. They allow arbitrary signs, phases and strengths of the specified 70 coefficients within the budget. Changes in birth distribution, additional current-generated fields, collisions or evolution of the source and geometry can change the bounds. A low pressure can exclude local equilibrium in these tests. Larger pressures would still leave the momentum distribution to be checked.
