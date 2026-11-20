# Checking the parent energy and its source

The coupled run excites an initially empty homogeneous parent field. An independent integration follows that field on the saved wall trajectory. A separate two-field example asks whether a finite source that loses energy can excite a much heavier parent while its oscillation action is approximately conserved after release. The two tests use different geometries and masses. Their energy budgets apply to their respective histories.

## Independent integration on the coupled history

The only trajectory inputs for advancing the parent are $a,H,\psi,\dot\psi$. Cubic interpolation uses the measured proper-time velocities. For the minimally coupled homogeneous parent use $P=a^3\dot\chi$ and the turning minimum $v(\psi)$. The canonical equations are

$$
 \dot\chi=P/a^3,\qquad \dot P=-a^3 M^2(\chi-v).
$$

The replay starts with both canonical coordinates zero. It freezes $a^3,M,v$ at each substep midpoint and applies the exact harmonic-oscillator rotation about $v$. This differs from the Runge–Kutta integrator that advances the full wall. The geometric and driving inputs for the replay include the effect of the parent on the original coupled solution. Midpoint quadrature independently computes the work from $-a^3M^2(\chi-v)\dot v$ and $a^3M\dot M(\chi-v)^2$, and the pressure contribution $3Ha^3p_\chi$. At each time the comoving energy equals the sum of both work terms less the pressure contribution. [Replay](trace.py) · [Results](TRACE.json)

| Epoch | Independent $J/J_{\rm target}$ | Coupled $J/J_{\rm target}$ |
| --- | ---: | ---: |
| Minimum finishes moving | 1.5141391 | 1.5141391 |
| Parent mass starts rising | 1.7544088 | 1.7544086 |
| First decay threshold | 1.6735161 | 1.6735162 |

Here $J=a^3\rho_\chi/M$ is an instantaneous classical oscillation action and $J_{\rm target}=2^3(10^{-4})/1.9$ matches the previous cold-parent normalization per unit coordinate volume. Event ratios agree within $1.1\times10^{-7}$ in relative terms. Subsampling the trajectory by two changes the independently calculated event and final ratios by at most $2.6\times10^{-6}$. The final drive work agrees to $2.3\times10^{-11}$ and the mass work to $6.0\times10^{-11}$ in units of comoving energy. The largest independently reconstructed energy residual is below $10^{-10}$.

The meaning of this action depends on the epoch. Between the end of the drive and the start of the mass rise, archived samples span **1.082 to 2.167** times the reference. After $v$ has returned to zero, the equation of motion gives

$$
 \frac{\dot J}{J}=-\left(3H+\frac{\dot M}{M}\right)\frac{p_\chi}{\rho_\chi}.
$$

For an expanding wall and nondecreasing mass, $|p_\chi|\le\rho_\chi$ bounds $J_f/J_i$ between $(a_i/a_f)^3 M_i/M_f$ and its reciprocal. Applying this to the brief final mass rise gives **1.6501 to 1.8653** times the reference at opening for any initial oscillator phase at the measured pre-rise action, on the same prescribed geometry and mass history. This bound needs no slow-variation assumption. At opening $|\dot M|/M^2=5.19$ and $H/M=0.172$. The rapidly changing frequency and coherent oscillations limit interpreting instantaneous action as a count of independently decaying particles. [Native trajectory](../coupled/PARENT.json)

## A finite source driving a heavy parent

Take two minimally coupled homogeneous fields on a prescribed expanding geometry. The local Lagrangian per unit coordinate volume is

$$
 L=a^3\left[\tfrac12\dot s^2+\tfrac12\dot\chi^2-\tfrac12 M^2(\chi-v(s))^2\right].
$$

The source coordinate $s$ rolls with only its supplied initial kinetic energy. Both $\chi$ and $\dot\chi$ start at zero. Give $v(s)=g[S(2s/w)-S(2s/w-1)]$ support on $0<s<w$, with the same compact smooth step as the coupled study. Writing $\delta=\chi-v(s)$, the scalar equations give $\ddot s+3H\dot s=M^2\delta v_{,s}$ and $\ddot\chi+3H\dot\chi=-M^2\delta$. The interaction is absent once $s$ has crossed the finite region. Source and parent receive opposite proper-time work per unit comoving volume, with $-a^3 M^2\delta v_{,s}\dot s$ delivered to the parent. The sum obeys $d(a^3\rho_{\rm tot})/d\tau=-3Ha^3p_{\rm tot}$. This construction supplies a finite precursor energy and specifies its interaction. It does not derive that source from a blast. [Calculation](driver.py) · [Full results](DRIVER.json)

For a numerical example take the archived wall expansion beginning at proper time 4.52 with $a=4.93445$. This is a short external geometry segment from the earlier release. Set $M=237.5$, corresponding to $mL=125$ and initial mass ratio $0.95$. Supply source kinetic energy $a_i^3\dot s_i^2/2=0.0016$, twice the former cold-parent energy target 0.0008. In the primary example $M w/\dot s_i=4$ and $M g/\dot s_i=0.5$. These are exploratory inputs. With $g$ fixed in these dimensionless units, the comparison varies the field width and retains the finite source reservoir. [Geometry data with archive digest](geometry.json) · [Earlier wall release](https://doi.org/10.5281/zenodo.24406309)

| Source choice | Parent action at exit / reference |
| --- | ---: |
| Primary | 1.44507 |
| Half the displacement amplitude | 0.33454 |
| Half the field width | 0.55148 |
| Twice the field width | 0.66728 |
| Zero displacement | 0 |

In the primary case, crossing takes 0.02115 proper time units. Net work into the parent is **72.12%** of the source's original comoving energy. At exit **72.25%** is parent energy, **26.15%** remains in source kinetic energy, and **1.60%** has gone into expansion work. The signed work and the final parent share differ because the parent pressure contributes during the motion. The total energy ledger closes within $8\times10^{-12}$ of the original source energy. Halving the oscillator time step shifts the reported action ratio by less than $2\times10^{-11}$. With $H/M=0.00168$ at exit, the subsequent free oscillation stays within **1.4394 to 1.4467** times the reference over one period. Prescribing the free source trajectory while omitting the force exerted on it produces a 64.84% mismatch in the combined budget. 

The heavy example uses the homogeneous mean and a constant post-pulse parent mass, and measures action per unit comoving volume. The earlier conformally coupled vacuum calculation counted fluctuations in a specified quantum state. The companion conformal calculation evaluates the coherent contribution of a source on that earlier history. Coupling the heavy example to the inflowing bulk, converting the coherent heavy field into individual particles, and heating daughters require additional dynamics. [Conformal calculation](../conformal/FORCED.md)

From the study root, run the following commands with Python 3.

~~~sh
python -I -B replay/trace.py --output /tmp/TRACE.json
python -I -B replay/driver.py --output /tmp/DRIVER.json
~~~
