# How large a transverse impulse can randomly located births generate?

For the first 1% of absolute electrical work, model each neutral pair as born at one location. Its daughters have opposite charges and opposite axial velocities. Birth locations are independent and uniform. The source rate, expansion, longitudinal fields and current normalization follow the corresponding archived wall history. Start with zero perturbing electromagnetic field. Under these assumptions the statistical fluctuations in where births occur give a very small source for transverse fields. An inequality covers the continuous injection history and all permitted three-dimensional Fourier directions in a specified local box. The tighter transverse-sector result and its finite-pitch generalization are given separately. [Numerical bounds](NOISE_BOUND.json) · [Linear field and particle equations](RESULTS.md)

| 1% whole-run work marker | Dedicated history | Shared history |
|---|---:|---:|
| Earlier necessary RMS sideways impulse to match the prescribed thermal spread | 541.65 $mc$ | 674.10 $mc$ |
| Conservative one-mode RMS impulse ceiling at any $0<k\le32768$, using each mode's own wavelength cube | $6.07\times10^{-10}\,mc$ | $1.25\times10^{-9}\,mc$ |
| Combined ceiling, 6,588,364 transverse wavevector pairs and cold axial birth cohorts, 1% margins | $4.38\times10^{-8}\,mc$ | $8.99\times10^{-8}\,mc$ |
| Combined ceiling, 17,990,702,160 wavevector pairs in all directions and finite birth momenta, 1% margins | $2.031\times10^{-4}\,mc$ | $4.518\times10^{-4}\,mc$ |

The table uses the half-step histories. The common patch has comoving side $2\pi/16$. The cutoff is comoving, so the shortest included physical wavelength is $2\pi a/32768$ in wall length units. The combined ceiling gives every mode full permission to add its RMS impulse in the same direction. For the transverse plane of cold streams the remaining deficit is over seven billionfold in each channel. Including the full three-dimensional band and finite birth momenta still leaves a deficit over a millionfold. This is a bound on a population and source-ensemble RMS. A single random realization can fluctuate. The previous impulse target tests recovery of the earlier model's thermal distribution. [Earlier target and its assumptions](https://github.com/maldonado-research/HDblast/blob/6a62a0798713dbe28d0d670376896893cd8d79a1/research/HDBLAST_KINETIC_TEST_20261121/ballistic/KICK.md)

## Bounding the linear response in time

Use the four-component comoving field $\widetilde f=(a/a_0)^2(B_y,E_z,B_z,E_y)$ of one spatial sine/cosine quadrature. Let $F(t)$ bound its source-ensemble root mean square Euclidean norm. Write $r=a_0/a$, $S=Q^2/(K a a_0^2)$ and $D$ for total comoving daughter rest in the effective wall field equation. A cold cohort born at $b$ has axial velocity $v_b$ and energy $\gamma_b m$ per particle; this subscript denotes a birth time. Its scaled perturbation momentum is $\mathbf P_b=(a/a_0)\delta\mathbf p/(mQ)$ in $c=1$ units.

Free Maxwell propagation preserves the Euclidean norm of the comoving field. The parallel guide field rotates $P_x,P_y$ and preserves their norm. A cohort is born with zero momentum perturbation. Its driving term is $r(-v_b\widetilde B_y,\widetilde E_y,\widetilde E_z)$, of norm at most $r\|\widetilde f\|$. Therefore $P(t)=\int_0^t rF\,dt'$ bounds the RMS norm of every cohort's momentum perturbation. Time here starts at the first birth interval's left endpoint.

Define $A(t)=\int^t dD_b/\gamma_b(t)$. Since the perturbation density divided by $Q$ has derivative $(k/a)rP_{x,b}/\gamma_b$, swapping the time integral and the integral over all births gives an upper bound on the density-induced current from the motions. It is $U(t)=\int^t (k/a)r A P\,dt'$ in the dimensionless units below. The direct velocity contributions to the current have coefficients $1/\gamma_b$ and $1/\gamma_b^3$, both bounded by $1/\gamma_b$. The triangle inequality for root mean squares and the two norm-preserving propagators then give

$$
 F(t)\le\int_0^t S\big[rA P+U+\eta\big],dt',
\qquad P(t)=\int_0^t rF\,dt',\qquad
 U(t)=\int_0^t(k/a)r A P\,dt'.
$$

Here $\eta$ bounds the current from the directly injected spatial modulation. The comparison includes the field responding to the particles and the induced density responding to that field. A single encounter between a tagged particle and a discrete partner is outside this mean-field calculation.

## Statistical normalization of one neutral pair

In a comoving volume $V$, the mean number of neutral pairs of a particular flavor born in an interval with comoving rest increment $dD_b$ is $dD_b/w_1$, where $w_1=6\mu/(S_0V)$, $\mu=mL$ is the daughter mass and $S_0$ is the action scale. A point at a uniform location contributes $2w_1\cos(\mathbf k\cdot\mathbf x)$ to the real Fourier coefficient of the pair density. The two charges in that pair share its location. Each nonzero mode has zero mean and variance $2w_1dD_b$ in either quadrature. Distinct birth intervals have independent increments. Thus the RMS of the source current $\int v_b\,dD_{b,\rm fluct}/Q$ at time $t$ is at most $\sqrt{2w_1D(t)}/Q$. This uses $|v_b|\le1$ and requires no independence of the current at different times. Conditional on an existing tagged Poisson birth, the other births still form the original independent Poisson process. We allow the tag's own contribution of up to $2w_1/Q$ in each quadrature. The deterministic background density entering the response coefficients is held fixed in this linear calculation. Set

$$\eta=(\sqrt{2w_1D_*}+2w_1)/Q.$$

The force on any unperturbed particle uses the two real spatial quadratures. At any phase their combined norm is bounded by the square root of the sum of the squared quadrature norms. Thus the redshift-weighted transverse impulse has RMS at most $\sqrt{2}Q P(t_*)$, including the two quadratures and granting particles the entire interval. Centered fluctuations in a fixed total number of independent uniform pairs have the same nonzero-mode variance.

## Bounding a continuously injected population

The tabulated $X(t)$ is negative and $a(t)$ increases throughout the interval. If $b_0=\sqrt{1.05^2-1}$ bounds the birth momentum in mass units, then $|p_z|/m\ge[\alpha(t-t_b)-b_0]_+$ with $\alpha=Q\min|X|a_0/\max a$. This follows by integrating the axial force in comoving momentum. If $T$ is the full time window and $\dot D_{\max}$ is the largest birth slope of the linearly interpolated history, then at every time

$$
 A(t)\le A_* =\dot D_{\max}\left[\min(T,b_0/\alpha)+
 \operatorname{asinh}([\alpha T-b_0]_+)/\alpha\right].
$$

Use $S_*=Q^2/(K a_0^3)$, $r\le1$, $a\ge a_0$, $k\le k_{\max}$. With $x=t/T$, $\beta=S_*A_*T^2$ and $\zeta=\beta k_{\max}T/a_0$, divide $F$, $P/T$ and $S_*TU$ by $S_*T\eta$. They are bounded componentwise by the solution at $x=1$ of the following positive linear system. Its initial state is $(0,0,0,1)$.

$$
 \frac{d}{dx}
 \begin{pmatrix}f\\p\\u\\1\end{pmatrix}=
 \begin{pmatrix}0&\beta&1&1\\1&0&0&0\\0&\zeta&0&0\\0&0&0&0\end{pmatrix}
 \begin{pmatrix}f\\p\\u\\1\end{pmatrix}.
$$

The one-mode bound in physical impulse units is $\sqrt{2}Q S_*T^2\eta\,p(1)\,mc$. Matrix exponentiation, direct time integration and an entrywise nonnegative Taylor sum with a norm bound on its remainder check this number. The calculation includes a small upward numerical allowance. With one wavelength per box dimension, $w_1$ and $\zeta$ increase with $k$, so the worst ceiling occurs at the high end of the band. In the common box, $w_1$ is the same for every mode. The number of nonzero integer lattice pairs $\pm(n_x,n_y)$ inside the disk $n_x^2+n_y^2\le2048^2$ is 6,588,364. Multiplying by the common-box one-mode ceiling gives the conservative simultaneous result without assuming phases of different modes cancel. The padded calculation increases $D_*$ and $\dot D_{\max}$ by 1% and reduces the lower axial acceleration by 1%. [Code and complete results](noise_bound.py) · [Machine-readable output](NOISE_BOUND.json)

## Arbitrary three-dimensional wavevectors and finite birth pitch

The time-domain simulations resolve the cold transverse sector. The comparison bound can also allow oblique and longitudinal waves, both charges' distinct gyro orbits, and every birth momentum in the supported ball $|\mathbf p|/m\le b_0$. The linearization here is around the prescribed **spatially homogeneous collisionless orbits** of one effective gauge sector. The perturbation fields include all six Maxwell components. This extension gives an upper bound without needing to identify which oblique modes grow.

For one complex Fourier mode take the convention $\delta f(\mathbf x)=\operatorname{Re}[\delta f_k\exp(i\mathbf k\cdot\mathbf x)]$. In a cohort's Eulerian perturbation equation its free streaming adds $-i(\mathbf k\cdot\mathbf v/a)\mathbf P$ and the same scalar times $d$. This is a pure phase. The vacuum Maxwell curl operator is skew-Hermitian for every direction of $\mathbf k$, including its zero-frequency longitudinal part. Both properties preserve the corresponding free propagation norms. The pair source is locally charge neutral at birth. In a longitudinal perturbation its opposite charges can later separate. Maxwell evolution then includes the charge response through the full current.

The velocity derivative is $\delta\mathbf v=(I-\mathbf v\mathbf v^{\mathsf T})\delta\mathbf p/(m\gamma)$, of operator norm $1/\gamma$. Thus both the velocity-current and density-current estimates continue to use $A_*$. The forcing for a tagged charge has norm at most $\sqrt{2}r\|(\widetilde{\mathbf E},\widetilde{\mathbf B})\|$ by $|\mathbf v|\le1$. Its signed guide term in scaled momentum is the following matrix action (the charge sign can change its overall sign).

$$
 \frac{QB_0}{\gamma}\{\mathbf P-\mathbf v(\mathbf v\cdot\mathbf P)\}\times\hat{\mathbf z}.
$$

The $\mathbf P\times\hat{\mathbf z}$ term is a rotation. The rate of norm growth from the remaining term is at most $Q|B_0|\,|\mathbf v_\perp|/\gamma\le Q|B_0|b_0/\gamma^2$. On the unperturbed orbit its perpendicular momentum can only rotate and redshift. Integrating this rate using the earlier lower bound on $|p_z|/m$ gives the following upper bound for the logarithm of momentum-propagator growth.

$$
 g=Q\max|B_0|b_0\left[\min(T,b_0/\alpha)
       +\arctan([\alpha T-b_0]_+)/\alpha\right].
$$

Consequently every cohort's momentum perturbation is bounded by $P(t)=\sqrt{2}\exp(g)\int^t r F\,dt'$ where $F$ now bounds the full complex electromagnetic norm. Replace the one in the **second row, first column** of the comparison matrix by $\sqrt{2}\exp(g)$. For the unnormalized complex mode, one marked Poisson pair has a current Fourier coefficient of norm at most $2w_1$ before division by $Q$. This follows by averaging the two signed velocity vectors, which may have different ballistic phase factors. The Poisson increments from different births are centered for nonzero modes even if the marks of the two particles in a pair are correlated. Use $\eta=(2\sqrt{w_1D_*}+2w_1)/Q$ with the tagged-source allowance. The extra force in any real spatial phase has norm at most $\sqrt{2}$ times the complex field norm, already included in the new $P$. The impulse RMS ceiling becomes $Q S_* T^2\eta\,p(1)\,mc$. The half-step histories give padded $g\le0.00340$ in both channels. The numerical suite checks the velocity Jacobian, the Hermitian part of its guide generator, the Lorentz force norm and the free Maxwell generator at 513 finite-pitch states. [Independent generator checks](MODE_CHECK.json)

There are **17,990,702,160** nonzero integer wavevector pairs modulo sign in $n_x^2+n_y^2+n_z^2\le2048^2$. For each pair use the same common-box upper bound evaluated at the high-frequency cutoff. Even adding all RMS amplitudes by the triangle inequality, including an additional 1% in the guide field upper bound, gives **less than $2.031\times10^{-4}\,mc$ and $4.518\times10^{-4}\,mc$** for the dedicated and shared histories. The preceding target exceeds these totals by more than **2.66 million** and **1.49 million**, respectively. The calculation makes no cancellation assumption about different modes. The homogeneous $k=0$ component is the supplied background and is held fixed. The bounds concern the nonzero spatial modes. The coherent box size, ultraviolet cutoff, and absence of incoming fields are explicit. The result is conditional on linearized local Maxwell and collisionless particle dynamics. It does not include perturbed birth rates or momentum laws driven by the fluctuating field, collisions, nonlinear mode coupling, correlations between separate births, or a responding bulk.

Both bounds assume local Maxwell propagation, homogeneous deterministic expansion and mean birth history, linear response and zero externally seeded perturbation fields. The acquired kicks predicted from discrete Poisson fluctuations are small enough to be consistent with that linearization. The calculation leaves open wavelengths below the stated cutoff length, different spatial statistics of pair creation, perturbations in the parent or bulk, and a nonlinear change in the background. It measures the impulse available by this early deadline. Modes that grow later would address a different question.

From the study root run `python -I -B modes/noise_bound.py --history history --output /tmp/NOISE_BOUND.json`.
