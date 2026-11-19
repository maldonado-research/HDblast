# Fermion pairs from a population of moving parent scalars

A tagged parent with comoving momentum $q=8$ on the fine quiet shell generates $0.222940\,y^2$ excess pairs in the wider switching window, compared with $0.494140\,y^2$ from the parent at rest. For an isotropic population with one parent per momentum state up to $q=8$, the mean is $0.297910\,y^2$ pairs per parent. Near-positive, near-negative and distant-positive bulk disturbances change that mean by $-0.08919\%$, $+0.58212\%$ and $+0.06909\%$ relative to quiet. The three signs persist for a second parent distribution, a shorter window, both bulk meshes and a common daughter energy screen. [Numerical report](MOVING.json) · [Published rest calculation](https://github.com/maldonado-research/HDblast/tree/0c79f123c00918f16647df9c868b4ca889336dc5/research/HDBLAST_TAGGED_PARENT_20261119)

## Massive scalar and massless Dirac daughters

We evaluate the excess over a vacuum experiment that uses the same wall, switching and interactions. The wall is homogeneous and spatially flat. The real scalar parent is minimally coupled with $m/H_0=100|\phi-0.5|$. The fermion starts in the empty conformal vacuum. The Yukawa coupling is $y f(s)$, where $s=H_0\tau$ is the proper clock. We prepare the parent modes in the same instantaneous physical basis at the published initial time. With $\Pi=a^3\dot\chi$, the free evolution in units of $H_0$ is

$$
\dot\chi_q=\Pi_q/a^3,\qquad \dot\Pi_q=-a^3(q^2/a^2+m^2)\chi_q,\qquad d\eta/ds=1/a .
$$

At preparation $\chi_q=(2a_i^3 E_{q,i})^{-1/2}$, $\Pi_q=-i a_i^3 E_{q,i}\chi_q$. For an additional parent in the prepared mode write $v_q=a\chi_q$, $g=fv_q$, and $F(\omega)=\int g\,e^{i\omega\eta}\,d\eta$. The massless two-body invariant phase space is constant above the light cone. The spin-summed scalar Yukawa vertex contributes $2(\omega^2-q^2)$. The pair carries total comoving momentum $+\mathbf q$ in the annihilation channel and $-\mathbf q$ in the excess stimulated-creation channel. At leading order per extra parent,

$$
\frac{1}{y^2}\frac{dN_{\rm ann,stim}}{d\omega}
=\frac{\omega^2-q^2}{8\pi^2}\Theta(\omega-q)\,|F(\pm\omega)|^2 .
$$

Annihilation contributes $+1$ relative to vacuum. The matrix-element square for creating a parent in the occupied mode is $2$, compared with $1$ for the same switched vacuum. The two contributions add to the inclusive pair count and subtract for depletion of the initial mode. Multiplication by $\omega/a_f$ yields the pair energy at the common readout $s_f=2.10$. The massless daughters are observed after the coupling switches off. For both windows we use $f(s)=\sin^4[\pi(s-s_l)/(s_u-s_l)]$ within the support and zero outside. We evaluate $1.15\leq s\leq2.05$ and $1.25\leq s\leq1.95$.

## Finite momentum brings a threshold correction

The rest-particle time-domain count involves the derivative of the switched mode. At nonzero $q$, extending the polynomial $\omega^2-q^2$ over the whole real line includes an unphysical negative contribution from $|\omega|<q$. The exact inclusive identity is

$$
\frac{N}{y^2}
=\frac{1}{4\pi}\int (|g'|^2-q^2|g|^2)\,d\eta
+\frac{1}{8\pi^2}\int_{-q}^{q}(q^2-\omega^2)|F(\omega)|^2\,d\omega .
$$

The second term is positive. Omitting it at $q=8$ loses $4.35\%$ of the quiet wide-window count. For uniform occupation the corresponding shortfall is $2.15\%$. The code evaluates the local term with an independent proper-time integral using the parent canonical momentum and the switching derivative. The refined local-plus-threshold and spectral counts agree within $7.0\times10^{-13}$ in absolute units per squared coupling. This particular comparison largely checks the full-line Parseval normalization since both readouts use the same below-threshold subtraction; the independent moving-parent propagator and alternate frequency integration test the threshold itself.

## What is being averaged

The parent density matrix is number-diagonal in the specified initial basis, isotropic and free of anomalous pair coherence. The chosen occupations are $n_u(q)=1$ or $n_t(q)=[1-(q/8)^2]^2$ on $0\leq q\leq8$, with zero occupation outside. All results are linear in $n$ at this order, and the comoving parent density is $n_c=\int q^2 n(q)\,dq/(2\pi^2)$. Its values are $8.64607434$ and $1.97624556$ in units of $H_0^3$. For a single real scalar the number of excess Dirac pairs per unit final physical volume is $n_c\langle N\rangle/a_f^3$. Each pair contributes two particles. The final daughter energy density is $n_c\langle E_{\rm pair}\rangle/a_f^3$. Isotropy and massless free propagation imply that the **outgoing pair population** has pressure one third of this density and zero mean momentum flux. These formulas describe the normal-ordered free particle content after switching.

| Fine quiet population and window | Pairs per parent / $y^2$ | Daughter energy per parent / $y^2H_0$ | Daughter energy density / $y^2H_0^4$ |
|---|---:|---:|---:|
| Uniform, wide | 0.297910 | 6.600639 | 30.247249 |
| Tapered, wide | 0.365664 | 7.971689 | 8.349725 |
| Uniform, narrow | 0.244195 | 5.189158 | 23.779177 |
| Tapered, narrow | 0.300975 | 6.262174 | 6.559141 |

For the uniform wide-window quiet case, the count splits into $0.235501$ from annihilation and $0.062409$ from excess stimulated creation in units of $y^2$. Both pieces carry positive daughter energy. The signed depletion coefficient of the initial parent number is $0.173091$. In the same order the local energy density contrasts are $-0.07694\%$, $+0.71355\%$ and $+0.12542\%$. They differ from the contrasts in mean energy per parent because $a_f$ varies across the shell histories.

With the common cutoff $\omega\leq60$, the fine quiet wide uniform count and energy per parent are $0.296310$ and $6.516946$. The largest possible individual comoving daughter momentum equals $(\omega+q)/2$ by two-body kinematics. For all the sampled windows and histories with $q\leq8$, the largest such screened physical momentum is $44.459H_0$, below the historical model scale $M_5=53.611H_0$. The full-spectrum calculation has no such ceiling. The screen bounds external daughter momenta on this interval and makes no claim about other short-distance corrections.

Both parent distributions are prescribed extra matter on stored geometries. These leading coefficients require perturbative $y$. A full dynamical energy budget would also include the initial parent stress, interaction and switching work, scalar work, depletion, and the response of the geometry. Energy supplied by switching and by the varying mass can contribute to the observed pairs.

## Numerical checks

The code imports the published wall loader, switching function, initial normalizations and rest reference under their SHA-256 pins. It evolves each group of proper-time modes with DOP853. It constructs a monotonic conformal clock for the Fourier waveform using cubic Hermite interpolation. It integrates the timelike spectrum across the true threshold and the fixed energy screen using cubic interpolation below the screen and Simpson integration on the larger Fourier grid. The requested tagged momenta $0,2,4,6,8$ are separate evaluations from the Gauss–Legendre quadrature nodes used for the distributions.

The baseline uses 16 nodes in $q$, 1024 samples of the switched conformal waveform and 32-fold zero padding. The $q$-only check doubles the quadrature order. The final refinement doubles the waveform resolution and padding as well. Relative to the published rest study, the largest refined discrepancy across total number, energy and net depletion on both windows and both screen choices is $1.90\times10^{-9}$. The largest shift in the inclusive per-parent number or energy contrast across the simulations on the two bulk meshes is $0.00271$ percentage point. This difference between mesh resolutions is an empirical convergence diagnostic. [Resolution and contrast checks](CHECK.json)

From this directory, use a checkout of HDblast containing the earlier tagged rest study, the B1 wall archive and quiet rolling summaries. NumPy and SciPy are required.

~~~sh
OPENBLAS_NUM_THREADS=1 python -I -B moving.py /path/to/repository > moving-rerun.json
OPENBLAS_NUM_THREADS=1 python -I -B moving.py /path/to/repository --nodes 32 > q-refined.json
OPENBLAS_NUM_THREADS=1 python -I -B moving.py /path/to/repository --nodes 32 --samples 2048 --padding 64 > refined.json
python -I -B checks.py moving-rerun.json q-refined.json refined.json > checked.json
~~~
