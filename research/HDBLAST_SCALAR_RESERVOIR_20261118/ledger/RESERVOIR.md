# What supplied the energy in the B1 radiation calculation?

In this B1 simulation the shell scalar supplies energy to a massive particle population, $\chi$. The particles subsequently decay into radiation. On the fine grid, about **60.20% of the integrated scalar work** comes from the increasing energy of the existing particles as their mass changes. The other **39.80%** pays for particles added during the specified production ramp. Expansion takes about 19.91% of this net scalar input as pressure work. The radiation grows while the particle population loses energy through decay. This calculation follows a particular sourced model with detuning $0.1$, $\phi_*=0.5$, $G=100$, $y=1$ and cutoff $\lambda_c=1$. The original run recorded no qualifying radiation plateau. [B1 model and numerical limits](https://github.com/maldonado-research/HDblast/blob/a011dee2b490faa56ae81d3165b4b5fa5bbfe636/research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source/README.md)

## The two energy accounts

Use forward proper time $s=H_0\tau$, $h=H/H_0$, $a=\exp[A_b(s)-A_b(0)]$, $r_b=1/H_0$ and $b=\kappa_5^2H_0^3=6.49\times10^{-6}$. The comoving particle and radiation energies per unit reference volume in the common normalization are $E_\chi=a^3 b\widehat\rho_\chi$ and $E_r=a^3r_bR$, where $R=\kappa_5^2\rho_r$ is the stored solver variable. Hatted matter densities and rates use the solver's $H_0$ units. With dots denoting $d/ds$, the five rates entering the comparison are

$$
\begin{aligned}
W_m &= a^3 b\,\widehat j\,\dot\phi_b,&
W_f &= a^3 b\,\widehat j_{\rm prod}\,\dot\phi_b,\\
D &= a^3 b\,\frac{y^2\widehat m^2\widehat n}{8\pi},&
C_\chi &= 3h a^3 b\,\widehat p_\chi,&
C_r &=h E_r .
\end{aligned}
$$

The B1 particle energy and pressure moments follow the massive dispersion relation $\omega_i^2=p_i^2+\widehat m^2$. The shell scalar changes $\widehat m^2=G^2(\phi_b-\phi_*)^2$ in this zero-bare-mass case. Differentiating the energy sum with evolving occupations gives $\dot E_\chi=W_m+W_f-D-C_\chi$. The decay rate includes the factor $\widehat m/\omega_i$ from time dilation. Multiplying by the energy per decaying particle cancels $\omega_i$ and gives the stated $D$. The radiation equation is $\dot E_r=D-C_r$. Their sum is
$$
\frac{d}{ds}(E_\chi+E_r)=W_m+W_f-C_\chi-C_r .
$$
Here $3h(E_r/3)=C_r$ accounts for the pressure of radiation. The earlier $a^4 R$ diagnostic removes this dilution from the radiation-only equation. Using $a^3$ for both components allows their internal transfer to cancel exactly in the combined equation. [Microscopic law and recorded quantities](https://github.com/maldonado-research/HDblast/blob/a011dee2b490faa56ae81d3165b4b5fa5bbfe636/research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source/shell_matter.py) · [Previous radiation-only investigation](https://github.com/maldonado-research/HDblast/issues/157)

## Reconstructing the inputs

For mass work, use $\widehat\rho_\chi-3\widehat p_\chi=\widehat m^2\sum_i n_i/\omega_i$. Thus $\widehat j=(\widehat\rho_\chi-3\widehat p_\chi)/(\phi_b-\phi_*)$ away from the crossing; the original finite node sum has $\widehat j=0$ at exact crossing. For filling work, regenerate the 24 generalized Gauss–Laguerre momentum nodes and target populations from the crossing record, using the production factor and ramp width for the selected run. For each saved progress $u$, the production derivative is $30u^2(1-u)^2$ inside $0<u<1$ and zero outside. Multiplication by the particle energy, dilution relative to crossing, the direction of the crossing, the registered stall switch and the actual scalar speed gives the work of filling the cohort. The separately saved population density and field value give the decay input. These reconstructions agree with the solver's reported force and power columns to floating-point precision. They use no finite difference of the radiation density and no saved energy derivative. The particle moments, occupation progress, radiation state and geometry all come from the same simulation.

At cohort creation, the finite step crosses $\phi_*$ before its endpoint. The implementation initializes the cohort at the small nonzero ramp fraction corresponding to that endpoint. It reports an initial dimensionless density jump $r_b\Delta X$ of $2.8763\times10^{-12}$ on the coarse grid and $5.6828\times10^{-13}$ on the fine grid. Multiplying by the interpolated crossing value of $a^3$ gives approximately $1.90\times10^{-10}$ and $3.75\times10^{-11}$, respectively. The creation-step endpoint has a slightly different scale factor. The continuous rate columns omit this discrete positive addition; the reported residuals below are the direct, uncorrected rate-minus-gain values. The jump is far smaller than the output-cadence residual in either run. [Cohort insertion and recording](https://github.com/maldonado-research/HDblast/blob/a011dee2b490faa56ae81d3165b4b5fa5bbfe636/research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source/evolve_b1.py)

## Results on the same interval

Both grids cover $1\leq s\leq4$, inside their reported reliability bounds. This is the interval used for the earlier radiation-only reanalysis. Linearly interpolate each weighted quantity to the endpoints and integrate the piecewise-linear rates by trapezoids in $s$. The 155 points on each grid include the two endpoints.

| Quantity on $[1,4]$ | Coarse grid | Fine grid |
| --- | ---: | ---: |
| Scalar work through mass change | 0.06951064 | 0.06953286 |
| Scalar work through cohort filling | 0.04595534 | 0.04597003 |
| Internal decay transfer | 0.11381177 | 0.11384805 |
| Particle pressure work | 0.00049755 | 0.00049771 |
| Radiation pressure work | 0.02248979 | 0.02249600 |
| Gain in particle energy | 0.00115386 | 0.00115447 |
| Gain in radiation energy | 0.09132200 | 0.09135208 |
| Total energy gain | 0.09247585 | 0.09250655 |
| Net input minus total gain | $2.7880\times10^{-6}$ | $2.6292\times10^{-6}$ |

The fine-grid residual is 28.4 parts per million of the total gain. Two neighboring subintervals, $[1.375,1.75]$ and $[1.75,2.125]$, give residuals of about $-6.19\times10^{-5}$ and $+6.51\times10^{-5}$. Their cancellation improves the aggregate. Keeping every fourth output on the fine grid with the same endpoints gives a total residual of $-6.65\times10^{-5}$. The measured total gain changes by 0.0332% when the radial mesh spacing is halved. The saved time cadence is comparable on both grids; these figures do not provide a continuum error bound.

While decay continues, the total comoving particle-plus-radiation energy decreases in each of the five late subintervals after $s=2.125$. During those bins the external scalar work is less than the pressure work. Omitting mass work or filling work produces fine-grid residuals of $-0.06953$ or $-0.04597$. Omitting both pressure terms gives $+0.02300$. Counting decay again as an external source gives $+0.11385$. These are balance checks on the same saved trajectory, with the dynamics held fixed. [All eight intervals and controls](RESERVOIR.json)

The matter-side work term measures the energy delivered to these two components in the B1 equations. The available outputs do not independently measure the loss of an integrated bulk scalar or gravitational energy reservoir: that would require the corresponding bulk fields, fluxes and boundary history. The original runs reached radiation fractions of about 1.16% under the B1 classifier; production alone did not establish a radiation-dominated era. The optical shell study uses a different local family and amplitude normalization. Converting these units to a cosmic thermal history requires additional matter physics and a physical identification of the production epoch.

## Reproduce

The standard-library checker takes the two publicly distributable, derived series and prints the component balances and controls. Run it from this directory. An optional extractor regenerates the series from the privately archived NPZ files and their two public summaries using NumPy and SciPy. It verifies seven source documents by SHA-256. The public solver, summaries and B1 report match their pinned revision; the two NPZ inputs match their private archive revision. [Independent extraction from the saved particle moments](extract_reservoir.py) · [Input digests and crossing data](EXTRACTION.json) · [Coarse series](reservoir_coarse.csv) · [Fine series](reservoir_fine.csv)

~~~sh
python -I -B reservoir.py > /tmp/hdblast-reservoir.json
cmp RESERVOIR.json /tmp/hdblast-reservoir.json
~~~
