# What does each fermion carry away?

The quiet shell produces **0.769816 pairs per normalized comoving volume** with the factor $y^2 H_0^3$ removed. At the readout time $H_0\tau=2.10$, an individual particle carries mean proper energy **$14.6207 H_0$**. The pair-energy band $60<\omega<120$ contributes 2.30% of the total particle count and 4.14% of the outgoing energy for the full modeled interval $q<\omega<120$. The earlier study checked the energy in that upper band. The full interval supplies most of the possible scattering targets below $\omega=60$. [Previous ordered energy calculation](https://doi.org/10.5281/zenodo.24374531) · [Numerical moments and controls](DISTRIBUTIONS.json) · [Momentum spectra](SPECTRA.csv)

## From pairs to individual momenta

Use comoving momenta in units of $H_0$ and natural units $c=\hbar=1$. For a pair with total comoving energy $\omega$ and momentum magnitude $q$, the spin-summed matrix element for a scalar producing two massless Dirac fermions is independent of direction in the pair rest frame. If $u$ is the cosine of the rest-frame emission angle to the boost, then $u$ is uniform on $[-1,1]$. A boost to the shell frame gives the two comoving momentum magnitudes

$$
 k_1=(\omega+qu)/2,\qquad k_2=(\omega-qu)/2 .
$$

Thus either daughter has constant marginal density $1/q$ on $[(\omega-q)/2,(\omega+q)/2]$. The momenta in one pair are correlated. This formula integrates the one-particle distribution and retains the angle average exactly.

For the 24 native parent momentum weights $w_q$, define the differential pair count per unit comoving volume, with $y^2H_0^3$ removed, by

$$
 R_q(\omega)= w_q \frac{\omega^2-q^2}{8\pi^2}\langle Z_q(\omega)^\dagger Z_q(\omega)\rangle .
$$

The archived ordered source amplitude $Z$ includes the native gate and two free parent mode functions. The daughter field starts in its conformal vacuum. For either one sign separately, $\psi$ or $\bar\psi$, the comoving distribution $D(k)=dN/dk$ and proper momentum at the readout are

$$
 D(k)=\sum_q \frac{1}{q}
 \int_{\max(L_q,2k-q)}^{\min(U,2k+q)} R_q(\omega)\,d\omega ,
 \qquad p=k/a_f.
$$

An empty integration interval contributes zero. The complete modeled range has $L_q=q$ and $U=120$. The low and high contributions are split at $\omega=60$. This cutoff is on **total pair energy**. Either sign of daughter can occupy momenta on both sides of $k=30$.

Write $M_j=\int k^j D(k)\,dk$. The code integrates $M_{-1},M_0,\ldots,M_6$ using the conditional angular average. For nonnegative integer $j$ that average is

$$
 \langle k^j\rangle_{q,\omega}=(\omega/2)^j
 \sum_{\substack{0\le l\le j\\ l\,\mathrm{even}}}
 \binom{j}{l}\frac{(q/\omega)^l}{l+1},\qquad
 \langle k^{-1}\rangle_{q,\omega}=\frac{\ln[(\omega+q)/(\omega-q)]}{q}.
$$

The logarithm at threshold is integrable against $R_q\propto\omega^2-q^2$. The massless-pair number is $M_0$ and its conformal energy is $2M_1$ for both signs. At scale factor $a_f$, both-sign proper number and energy densities are $2y^2H_0^3 M_0/a_f^3$ and $2y^2H_0^4 M_1/a_f^4$. The proper energy per comoving volume is $2y^2H_0^4 M_1/a_f$, the normalization checked against the preceding study.

| Wall history | Pair count $M_0$ | Number density, both signs, divided by $y^2H_0^3$ | Energy density, both signs, divided by $y^2H_0^4$ | Mean proper daughter energy, in $H_0$ |
|---|---:|---:|---:|---:|
| Quiet | 0.769816 | 0.815965 | 11.93001 | 14.62074 |
| Near positive | 0.766784 | 0.812826 | 11.88747 | 14.62487 |
| Near negative | 0.782563 | 0.829593 | 12.12298 | 14.61318 |
| Far | 0.769900 | 0.816181 | 11.93541 | 14.62347 |

## An upper bound on occupation at every momentum

The unpolarized one-sign distribution has two spin states. Each has dimensionless occupation $f(k)=\pi^2 y^2 D(k)/k^2$ at this order. Define

$$
 B_q=\int ds\, |F(s)|\sqrt{\langle\chi_q(s)^2\rangle}.
$$

The ordered expectation defines a positive seminorm, $\|A\|^2=\langle A^\dagger A\rangle$. Applying its triangle inequality to $Z_q=\int ds\,F(s)e^{i\omega\eta}\chi_q(s)$ bounds the source power by $B_q^2$ for every frequency. The remaining integral of the positive phase-space factor has an elementary bound. Removing the energy cutoff can only increase it. For $k,q>0$ one has

$$
 \frac{1}{8q k^2}\int_{\max(q,2k-q)}^{2k+q}(\omega^2-q^2)\,d\omega
 =\begin{cases}\tfrac12+\tfrac{k}{3q},& k\le q,\\
                1-\tfrac{q^2}{6k^2},& k\ge q .
 \end{cases}
$$

Both branches are at most one. Therefore

$$
 f(k)\le y^2\sum_q w_q B_q^2 .
$$

The numerically evaluated coefficient is **0.0163959** on the quiet wall and **at most 0.0164889** across the four histories, on the retained parent momentum grid. The upper bound applies pointwise, including the limit at low $k$. The largest measured full-band occupancy averaged over the archived momentum bins is smaller, about $0.00229 y^2$ on the quiet wall. For an illustration at $y=0.1$, the uniform bound gives $f\le 0.000165$. In a factorized kinetic equation, two outgoing Pauli factors then obey $(1-f_3)(1-f_4)\ge 0.99967$ for this **initial** emitted distribution. This conclusion uses a single unpolarized species, an initially empty daughter field, and leading-order production. Later evolution, pre-existing populations and phase correlations require their own treatment.

## Numerical checks and replay

The calculation draws the archived canonical parent solutions and ordered operator sign from the pinned preceding study. It evaluates the pair spectra directly on Gauss–Legendre frequency nodes. For bands beginning at $\omega=q$ it maps the node $x\in[0,1]$ to $q+(U-q)x^4$. This resolves the logarithm in the inverse-momentum moment. A separate Chebyshev interpolation of $R_q$ yields piecewise antiderivatives for the exact number and energy in each of 128 bins over $0<k<64$. A fresh Fourier integral folded through the permitted angular interval checks eight bins, including the lowest and highest bins and the interval around the maximum quiet occupation.

The largest differences across the twelve wall–band comparisons were $1.5\times10^{-10}$ relative to archived outgoing energies, $3.6\times10^{-11}$ under joint frequency and native-time refinement for raw moments, and $2.6\times10^{-11}$ for the sum of absolute bin changes divided by total number or energy. Bin totals agree with directly integrated number and energy to $1.1\times10^{-14}$ relative. Five-point Gauss angular integration reproduces the analytic nonnegative moments within $5\times10^{-16}$ relative. Independent direct folding of the eight selected bins differs from the cumulative construction by at most $2.4\times10^{-14}$ of the total number or energy. On all four walls, the eight full-band moments agree with the sum of the independently integrated high and low bands to $1.0\times10^{-14}$ relative. These controls use the previous study's fixed native parent grid and wall histories.

Each CSV row reports the integral over a momentum bin for one particle sign, with the factors of $y^2$ and $H_0$ removed as above. The two **per-log** columns evaluate the proper densities at the bin midpoint, giving $kD(k)/a_f^3$ for number and $k^2D(k)/a_f^4$ for energy. The last column gives $f/y^2$ averaged with phase-space weight $k^2dk$ across the bin.

From the study root, the following regenerates both machine-readable files. The repository path points to the directory containing the pinned historical research capsules.

~~~sh
OPENBLAS_NUM_THREADS=1 python -IB daughters/daughters.py --repo ../.. --output /tmp/hdblast-daughters
~~~

[Calculation](daughters.py) · [Moment and check data](DISTRIBUTIONS.json) · [Physical and comoving momentum bins](SPECTRA.csv)
