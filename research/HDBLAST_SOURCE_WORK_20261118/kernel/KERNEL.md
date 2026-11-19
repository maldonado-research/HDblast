# How the shell history sets the late particle count

A single response calculated on the **quiet rolling shell** recovers the sign of the directly propagated particle excess in all four pulse histories. Its largest error on either bulk grid is about 0.33% of the small change in the count. The calculation takes a history of the mass and geometry as input. Predicting that history from the initial five-dimensional disturbance still requires the bulk dynamics. [Numerical results](KERNEL.json) · [Input study](https://doi.org/10.5281/zenodo.24355584)

## Observable and response functions

Use the common wall proper time $s=H_0\tau$, comoving momenta $0\le\kappa\le8$, and the quiet crossing scale factor as the unit scale factor. For the minimally coupled real field use $X=a^{3/2}\chi$ with canonical oscillator data $(X,X')$; primes in this section denote $d/ds$. In dimensionless units,

$$
X''+V_k X=0,\qquad
V_k=\kappa^2 a^{-2}+m^2-\tfrac94 h^2-\tfrac32 h',\qquad
m^2=10^4(\phi-0.5)^2,\qquad h=H/H_0.
$$

Start in the common second-order adiabatic state at $s_i=0.8$. Measure at $s_f=2.10$ with the prior study's adiabatic convention. Let $X_q$ be the quiet positive-frequency mode evolved forward; let $w_q$ solve the **same** quiet mode equation backward from the positive-frequency final analyzer $f_q$. For two functions $W[u,v]=uv'-u'v$. Set $\beta_q=W[w_q,X_q]$, constant in time, and $d\nu=\kappa^2d\kappa/(2\pi^2)$. The observable is the comoving band number $N=\int|\beta|^2d\nu$ in units of $a_{*,0}^3H_0^3$.

For a pulse history $p$, form the *difference of squared frequencies* $\Delta V_k=V_{p,k}-V_{q,k}$ at equal proper time. Differentiate the oscillator along the affine comparison $V_{\lambda,k}=V_{q,k}+\lambda\Delta V_k$ at $\lambda=0$. This defines a direction in the space of oscillator histories. Intermediate values of $\lambda$ need not solve the five-dimensional gravity equations. For initial data $z_i=(X_i,X_i')$ and final analyzer data $f=(f,f')$, the derivative is

$$
\delta\beta=W[w_q(s_i),\delta z_i] + W[\delta f,X_q(s_f)] -
  \int_{s_i}^{s_f}w_q(s)X_q(s)\Delta V_k(s)\,ds ,
\qquad \delta N=2\,\Re\int\overline{\beta_q}\,\delta\beta\,d\nu .
$$

The common initial data give zero initial term. Backward propagation computes the weight of an earlier perturbation in this specified final measurement. The integral has no dependence on pulse inputs after the readout. The moving endpoint analyzer in the prior study uses frequency samples extending $0.002$ on either side of $s_f$ for its derivative stencil. Its change is shown separately. A fixed quiet analyzer chosen for all arms removes that arm-dependent endpoint term. The proper-time and canonical-variable choices define the split between mass, geometry and analyzer in this comparison.

Only two scalar functions of time are needed for the bulk integral, one with an extra $\kappa^2$ in the momentum weight. For $j=0,2$ define

$$
K_j(s)=-2\,\Re\int \kappa^j\overline{\beta_q}\,w_q(s)X_q(s)\,d\nu .
$$

Their action is $\int [K_0(\Delta m^2-\tfrac94\Delta h^2-\tfrac32\Delta h')+K_2\Delta(a^{-2})]ds$. Here $\Delta h^2=h_p^2-h_q^2$. The supplied 2,049-sample table lets an independent calculation test new frequency histories against the same quiet reference. Linear interpolation of that table changes the fine-grid projection in the four examples by less than $1.8\times10^{-6}$ in $N$. [Quiet response table](RESPONSE.csv)

## Pulses on the simulated rolling shell

All entries are percentage changes relative to $N_q=3.51407719$ in the stated band. The mass, geometry and analyzer columns are signed contributions to the linear estimate, in percentage points.

| Fine-grid history | Direct change | Linear estimate | Mass history | Geometry history | Final analyzer |
|---|---:|---:|---:|---:|---:|
| Near compensator, positive front | −0.203329% | −0.203450% | −0.199378 | −0.002952 | −0.001121 |
| Near compensator, negative front | +0.991794% | +0.988551% | +0.982386 | −0.003286 | +0.009452 |
| Distant compensator, positive front | +0.051120% | +0.051129% | +0.054319 | −0.004700 | +0.001511 |
| Distant compensator, half-positive front | +0.026682% | +0.026685% | +0.028201 | −0.002261 | +0.000744 |

Keeping the output analyzer fixed to the quiet reference gives directly evolved fine-grid changes of −0.202209%, +0.982344%, +0.049658% and +0.025950%, in the same row order. This choice preserves all four signs without needing an arm-specific endpoint frequency derivative.

The near and distant full positive fronts give equal accumulated contributions through the quiet first crossing. Both reach $+0.00313129$ in $N$ at $s=1.48612$. By the proper time corresponding to quiet $T=1.25$, the accumulated contribution has fallen to approximately $+0.00041881$ for both. They agree there to $3.4\times10^{-12}$ in $N$. The nearby compensator first has causal access at $T=1.255$. The late mass-history contribution to $N$ is $-0.00700628$ for its positive front, compared with $+0.00190880$ when the compensator is distant. These partial integrals apportion the final measurement in time. A particle measurement performed at an earlier cut would have a different analyzer and response. The constructions include the metric needed by the bulk constraints. [Paired histories and causal geometry](https://github.com/maldonado-research/HDblast/tree/b5be2cb3c044457e20a9251658fd7e36f5e65864/research/HDBLAST_PULSE_PARTICLE_RESPONSE_20261118)

A finite-amplitude identity checks the response directly. Replace $X_q$ by the fully propagated $X_p$ in the source integral, use $W[f_p-f_q,X_p(s_f)]$ at the endpoint and the actual initial-data difference. The resulting complex number is exactly $\beta_p-\beta_q$ since $dW[w_q,X_p]/ds=-w_qX_p\Delta V_k$. Thus $\Delta N=2\Re\int\overline{\beta_q}(\beta_p-\beta_q)d\nu+\int|\beta_p-\beta_q|^2d\nu$. The positive quadratic term alone cannot predict the sign. Comparing the tangent with centered finite differences at $\lambda=\pm0.01$ changes $dN/d\lambda$ by less than $7.1\times10^{-8}$ on the fine grid. The solver also compares 32 and 48 Gauss–Legendre nodes, 4,096 and 8,192 time intervals, and both bulk meshes. Exact spectra in this note use Gauss–Legendre momentum integration and a tighter mode integration tolerance, which accounts for small differences from the preceding release. [Solver and checks](kernel.py)

## A timing test with identical crossing data

Keep the quiet geometry fixed and prescribe a small positive addition to the squared mass. Take $\Delta m^2/H_0^2=\sin^4[\pi(s-c+L/2)/L]$ inside a window of width $L=0.08$ centered at $c$, zero elsewhere. Its area in dimensionless units is $3L/8=0.030$. Every trial starts after the quiet crossing, leaves the full crossing neighborhood unchanged, and returns to the same outgoing mass trajectory. Applying the quiet kernel to 16 candidate centers between 1.60 and 1.90 selected the following two extrema. Direct integration of each model then tested the forecast.

| Window center | Linear change in $N$ | Direct change in $N$ | Direct change at half height |
|---|---:|---:|---:|
| 1.62 | +0.00802340 | +0.00803363 | +0.00401426 |
| 1.70 | −0.00550243 | −0.00549799 | −0.00275010 |

Moving a positive mass disturbance in time can reverse the change in the final count even though its area, shape, crossing data and final analyzer are the same. Relative to quiet, the two full-height counts change by +0.229% and −0.156%. A prescription that uses only derivatives of the mass trajectory at the crossing must predict identical results for these two inputs. The experiment is an externally prescribed drive of a free field. It establishes a counterexample in that model class and leaves the bulk realization and source energy open. The project had previously studied retarded mass-kick response on a stationary de Sitter comparison. Here the weights are calculated on the rolling B1 history, including its nonadiabatic crossing. [Earlier project response study](https://github.com/maldonado-research/HDblast/tree/b5be2cb3c044457e20a9251658fd7e36f5e65864/research/HDBLAST_CHECKPOINT_20261002_CAUSAL_RESPONSE)

Pass the preceding particle-study directory to the full NumPy/SciPy solver. The separate audit needs standard Python only and projects all 16 test windows on the supplied quiet table.

~~~sh
python -I -B kernel.py /path/to/HDBLAST_PULSE_PARTICLE_RESPONSE_20261118 --write-kernel regenerated.csv > regenerated.json
python -I -B audit_kernel.py
~~~
