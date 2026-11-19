# Quantum field energy supplied by the rolling B1 mass

On each of four rolling shell histories, the time-dependent scalar mass does positive net work on a freely evolving quantum field between $H_0\tau=0.80$ and $2.10$. The sign of the **additional** work relative to the undisturbed rolling shell depends on the whole disturbance. For the positive front it is negative with nearby compensation and positive when the compensating signal cannot arrive by readout. This calculation identifies the source term in the quantum field's energy equation. The scalar and spacetime histories are prescribed by the earlier B1 simulation, which used its own particle and decay rule. A self-consistent quantum loss from that scalar, decay into daughters and a temperature would require additional dynamics. [Input histories and particle study](https://doi.org/10.5281/zenodo.24355584)

## Band, clock and state

The homogeneous shell metric has flat spatial slices. Use $s=H_0\tau$, $h=d\log a/ds$, $a=1$ at the quiet mass crossing and $m=100|\phi_b-0.5|$ in units of $H_0$. We use the previous experiment's fixed comoving band $0\le k\le8$, with measure $d\nu=k^2dk/(2\pi^2)$ for one minimally coupled real scalar. Every arm starts from the same second-order adiabatic state on the common prehistory at $s=0.80$. We interpolate $h$ and $\phi_b$ in proper time using the preceding calculation's cubic splines, and integrate $h$ to obtain $a$. The near and far results apply to the respective constraint-matched initial perturbations. The previous release checks the band's physical-energy ceiling on these histories. [Input digests](INPUTS.json)

Let $q_k=\chi_k$ and $p_k=a^3dq_k/ds$, with $q_kp_k^*-q_k^*p_k=i$. Write $u=k^2/a^2$ and $\omega^2=u+m^2$. The equations $dq_k/ds=p_k/a^3$ and $dp_k/ds=-a^3\omega^2 q_k$ propagate each state, including phase through the crossing. Introduce $C_k=a^3|q_k|^2$, $K_k=|p_k|^2/a^3$. A spatially averaged comoving energy and pressure per mode are

$$
E_k=\tfrac12(K_k+\omega^2 C_k),\qquad
P_k=\tfrac12[K_k-(m^2+u/3)C_k].
$$

In the units used here $E=\int E_kd\nu=a^3\rho/H_0^4$ and $P=\int P_kd\nu=a^3p/H_0^4$. Equality of scale-factor normalization makes the compared comoving volumes consistent across histories. The explicit mass derivative and the time dependence of the geometry separately give

$$
E_f-E_i=W_m+W_a,\qquad
W_m=\frac12\int_{s_i}^{s_f}\frac{d m^2}{ds}\left(\int C_k\,d\nu\right)ds,\qquad
W_a=-3\int_{s_i}^{s_f}hP\,ds.
$$

The identity also follows by differentiating $E_k$. If $R_k=\mathrm{Re}(q_k^*p_k)$, the evolution yields $C_k'=3hC_k+2R_k$ and $K_k'=-3hK_k-2\omega^2R_k$. The two $R_k$ terms cancel in $E_k'$. These formulas continue through $m=0$.

## Absolute and paired energies

The quiet field gains $200.68208$ in comoving energy. Its mass supplies $207.96875$ and its net expansion work is $-7.28667$. The initial comoving energy is $203.06853$ on all four prepared states. The difference in final comoving energy therefore equals the difference in total accumulated work. The following numbers include the vacuum contribution within the stated band.

| Pulse minus quiet, fine mesh | Mass work $\Delta W_m$ | Expansion work $\Delta W_a$ | Final energy $\Delta E_f$ | Final energy above instantaneous vacuum, difference |
|---|---:|---:|---:|---:|
| Near positive | −0.35753 | −0.00161 | −0.35914 | −0.36306 |
| Near negative | +1.79640 | +0.02107 | +1.81747 | +1.81061 |
| Far positive | +0.10627 | +0.00471 | +0.11099 | +0.10397 |

The mass term accounts for 99.55%, 98.84% and 95.75% of the respective signed changes in raw band energy. These fractions describe the probe's balance on these trajectories. The calculation does not measure how much energy came across the shell from the bulk. The background already rolls and does work in the quiet experiment. [Full fine-grid integrals](REFINED.json)

## How much work an occupation-only description misses

At each proper time, diagonalize the physical field Hamiltonian in $(q_k,p_k)$. Set $n_{H,k}=E_k/\omega-1/2$. After dropping the anomalous pair correlation, a state with that instantaneous occupation has $C_{\rm diag}=E_k/\omega^2$ and $P_{\rm diag}=E_k u/(3\omega^2)$. Its vacuum parts are $C_0=1/(2\omega)$ and $P_0=u/(6\omega)$. The covariance difference $C_k-C_{\rm diag}$ and pressure difference $P_k-P_{\rm diag}$ isolate the omitted coherence in this basis. Every quadrature node has $k>0$. The isolated $k=0$ point at $m=0$ has zero measure in the band integral. A second, finer momentum quadrature checks this region.

| Fine mesh, difference in mass work | Instantaneous vacuum | Diagonal excitations | Coherence | Sum |
|---|---:|---:|---:|---:|
| Near positive | +0.00401 | −0.23162 | −0.12991 | −0.35753 |
| Near negative | +0.00391 | +1.17799 | +0.61449 | +1.79640 |
| Far positive | +0.00623 | +0.07463 | +0.02541 | +0.10627 |

The coherence terms supply about 36.3%, 34.2% and 23.9% of the respective mass-work differences. An occupation-only force misses this contribution even with the exact instantaneous number as input. This is an algebraic decomposition of the same evolving state, with a specified basis. The occupation used in the preceding particle study projected instead on a finite-time second-order adiabatic mode. The present state evolution reproduces its late quiet and pulse occupations to the accuracy of the changed momentum integration.

Integrating the instantaneous vacuum parts alone gives $\int (\omega_f-\omega_i)d\nu/2$, including expansion work. Subtracting that quantity produces the last column in the first table. This finite-band comparison checks sensitivity to a particular reference energy. It supplies no ultraviolet renormalization or determination of the local covariant counterterms of a four-dimensional theory.

## A paired decomposition of the source

A positive field-energy difference may arise from a change of mass history, a change of variance, or both. At common proper time write $C=\int C_kd\nu$ and let a bar denote the arithmetic mean of pulse and quiet. The symmetric product identity gives

$$
\Delta W_m = \frac12\int \overline C\;\Delta (m^2)'\,ds
             +\frac12\int \overline{(m^2)'}\;\Delta C\,ds .
$$

| Fine mesh | Changed mass derivative, raw | Changed variance, raw | Changed mass derivative, above vacuum | Changed variance, above vacuum |
|---|---:|---:|---:|---:|
| Near positive | −0.18371 | −0.17382 | −0.19760 | −0.16393 |
| Near negative | +1.26980 | +0.52660 | +0.99955 | +0.79293 |
| Far positive | +0.16333 | −0.05706 | +0.05903 | +0.04101 |

The two descriptions agree on the sign of the net mass-work difference in each case. In the far case even the sign assigned to the change in variance depends on the chosen reference. The symmetric rule shares product interactions between its two factors and does not determine a unique causal partition between the incoming pulse and the background roll. The report also gives the paired expansion split and partial work accumulated at 33 common clock cuts.

## Validation and replay

The central comparison uses 96 Gauss–Legendre momentum nodes and 8192 elapsed proper-time intervals for Simpson integration of the work; DOP853 evolves the canonical modes with relative tolerance $6\times10^{-11}$ and maximum internal step 0.0025. The largest absolute raw work–energy residual is $1.3\times10^{-8}$ on the refined quadrature. The mode norm differs from unity by less than $1.7\times10^{-11}$. Reducing the momentum nodes and time intervals to 64 and 4096 gives paired energy differences within $5.3\times10^{-8}$. On the coarser bulk grid they are −0.36396, +1.82559 and +0.10723. The spread of 0.00376 for the far positive signal across the two bulk grids is larger than the quadrature error. These two bulk grids test sensitivity to resolution and do not establish a continuum error bound. [Default quadrature](LEDGER.json) · [Coarse bulk grid](COARSE.json)

An exact control uses $a(s)=\sqrt{1+2Hs}$, with $m=m_0/a$ and $H>0$. The conformal mode $a\chi$ has constant frequency $\Omega=\sqrt{k^2+m_0^2}$. Its normalized plane wave has $E_k=\Omega/(2a)+H^2/(4\Omega a^3)$. The integrated mass work is $m_0^2\Delta(1/a)/(2\Omega)$ and the expansion work is $k^2\Delta(1/a)/(2\Omega)+H^2\Delta(a^{-3})/(4\Omega)$. Applying the production energy, pressure and work functions to these waves recovers both terms within $1.1\times10^{-12}$ in the two-mode test.

Use Python 3.12, NumPy 2.3.5 and SciPy 1.17.0. Supply the extracted HDBLAST_PULSE_PARTICLE_RESPONSE_20261118 directory from the preceding release. The script validates its input hashes and raises an error on a failed control, norm or work balance.

~~~sh
OPENBLAS_NUM_THREADS=1 python -I -B ledger.py PREVIOUS > ledger-replay.json
OPENBLAS_NUM_THREADS=1 python -I -B ledger.py PREVIOUS --order 96 --steps 8192 > refined-replay.json
OPENBLAS_NUM_THREADS=1 python -I -B ledger.py PREVIOUS --grid coarse > coarse-replay.json
~~~
