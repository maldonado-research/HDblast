# How much transverse field must persist during startup?

A strong field for one instant may give a small impulse. This calculation finds a lower bound on transverse field strength integrated over the whole interval in which the daughters are born. It allows an arbitrary time profile. The same field envelope must satisfy all birth ages at once.

At the 1% marker of the whole run’s absolute electrical work, the finer birth history gives the following necessary exposures. “Peak” here means the largest strength anywhere the particles can reach at a specified time. The RMS then averages its square over time from the first birth cell to the marker. Fields use the units of the recorded wall gauge field. The comparison allows an electric force and a velocity-dependent magnetic force to cooperate with maximum efficiency. [Numerical product](EXPOSURE.json)

| Required transverse distribution at the marker | Dedicated RMS peak field | Shared RMS peak field | Dedicated integrated peak energy / recorded axial field energy | Shared integrated peak energy / recorded axial field energy |
|---|---:|---:|---:|---:|
| Exact thermal transverse distribution | 0.007645 | 0.004839 | 0.1428% | 0.05783% |
| Within one percentage point at every transverse-momentum cut | 0.006006 | 0.004026 | 0.08816% | 0.04003% |

The native wall steps give 0.1303% and 0.05400% for the exact target. The finer bounds increase because the most recently born cohorts have shorter times to receive a force. At one-percentage-point tolerance, the native bounds are 0.08688% and 0.03967%. These lower bounds assume the recorded birth rate, expansion and thermal target, and the specified limit on the number of elastic collisions. A source modified by the added fields would need new histories. The temperature and unhit counts use conservative marker-cell enclosures in the birth-age calculation. [Birth histories and peak force test](../history/BUDGET.json)

## Field exposure as a constraint on all particle ages

Use units with the speed of light and the daughter mass equal to one. Measure proper time in $L$. For a particle’s field axis let

$$
G(t)=\sup_{\mathbf x}\sqrt{|\delta\mathbf E_\perp(t,\mathbf x)|^2+
                                      |\delta\mathbf B_\perp(t,\mathbf x)|^2}.
$$

The supremum can also include the three symmetry-related species. The homogeneous parallel electric field leaves transverse momentum unchanged. A parallel magnetic field only rotates it. The force that can increase its magnitude is at most $Q(|\delta\mathbf E_\perp|+|v_\parallel|\,|\delta\mathbf B_\perp|)\le\sqrt{2}QG$. A daughter that avoids the binary interaction after birth at $s$ starts with $p_\perp\le b=\sqrt{1.05^2-1}$. Its maximum additional transverse momentum at $t_f$ is bounded by $\sqrt{2}Q\int_s^{t_f}(a/a_f)G\,dt$.

For a suffix of recent births starting at $s$, let $h_s$ be the lower bound on its unhit share among all particles present at $t_f$. If all those particles remain below a given transverse-momentum cut, the actual cumulative distribution reaches at least $h_s$ there. Write $F$ for the target thermal transverse cumulative distribution and $\epsilon$ for the maximum permitted difference in the two cumulative distributions. Each suffix therefore imposes

$$
\int_s^{t_f} (a/a_f)G(t)\,dt \ \ge\ 
g_s=\frac{[F^{-1}((h_s-\epsilon)_+)-b]_+}{\sqrt{2}Q}.
$$

The numerical enclosure grants births the full left endpoint of their cell and extends the available response time through the marker cell’s right endpoint. Its $w_i=\min(1,a_{i+1}/a_{\rm left\,marker})$ overestimates redshift weighting. The physical field may be extended by zero from the actual marker to that endpoint. In increasing lookback order define $x=\int_s^{t_{\rm end}}w^2 dt$ and $Y(x)=\int_s^{t_{\rm end}} w G\,dt$. Then $\int G^2 dt=\int (Y')^2 dx$, with $Y(0)=0$ and $Y(x_i)\ge g_i$. The minimum is the piecewise-linear least concave majorant of the required points $(x_i,g_i)$. Where constraints are inactive, constant slopes minimize squared slope. At an active point the slope can only decrease. The computed minimizing exposure equals the sum of squared segment slopes times their widths. It is an optimistic lower bound on the physical field energy needed to accelerate real particles.

An independent certificate assigns each active point weight $\lambda_i=m_i-m_{i+1}\ge0$, where $m_i$ is the incoming segment slope and the terminal outgoing slope is zero. For any admissible field, Cauchy–Schwarz gives

$$
\int (Y')^2dx \ \ge\ \frac{(\sum_i\lambda_i g_i)^2}
{\int [\sum_i\lambda_i\,\mathbf1_{[0,x_i]}]^2 dx}.
$$

The code verifies equality of this certificate with the constructive minimum and checks a separate small problem by numerical quadratic optimization. In the finer histories at the 1% marker, combining all recent-birth constraints gives 2.08 times the integrated bound obtained from the strongest single suffix, for either interaction. [Code and controls](exposure.py)

## Independent check using the RMS impulse over all births

At fixed birth history, particles born before $s$ constitute fraction $D(s)/D_f$ of all particles present at $t_f$. Squaring the allowed integrated impulse and averaging over birth times yields a positive kernel proportional to $D(\min(s,t))/D_f$. In discrete time an entrywise upper bound is

$$
M_{ij}=\sqrt{\Delta t_i\Delta t_j}\,w_iw_j d_{\min(i,j)},\qquad
d_i=\min(1,D_{i+1}/D_{\rm left\,marker}).
$$

The squared RMS impulse is bounded above by $2Q^2\lambda_{\max}(M)\int G^2dt$. Matching a thermal transverse target, with an allowed cumulative difference $\epsilon$, requires squared RMS impulse at least $\int_0^{(h_0-\epsilon)_+} [F^{-1}(u)-b]_+^2du$. This gives another field lower bound. At exact matching it gives time-RMS peak strengths 0.002383 and 0.002482 in the finer histories, below the all-suffix requirements. The script bounds the matrix norm on both sides using a positive trial vector, compares to direct diagonalization on an irregular example, and recovers the continuum norm $4/\pi^2$ for uniform birth over a unit interval. [Study data](EXPOSURE.json)

## Spatial and energetic interpretation

On the recorded interval the wall gauge factor is constant. The ratio in the table divides $\int G^2dt$ by the recorded $\int(X_0^2+B_0^2)dt$; the common Maxwell energy normalization cancels for the selected sector. This is a ratio of time-integrated field energies along the benchmark. A magnetic force can redirect momentum while doing zero work. Both the force target and its duration can matter for any proposed source.

When the perturbation contains one real spatial Fourier harmonic, each transverse field component is a linear combination of the same sine and cosine. Over one period the supremum of the sum of squares is at most twice its spatial mean. For that case the required *spatially averaged* and time-integrated field energy is at least half the tabulated peak-energy bound. General concentrated fields can have much smaller spatial means at a given peak. Tests of such fields require the spatial support, phases and particle paths to translate this result into a volume energy budget.

The numerical command from the study root is

~~~sh
python -I -B energy/exposure.py --history history --output /tmp/EXPOSURE.json
~~~
