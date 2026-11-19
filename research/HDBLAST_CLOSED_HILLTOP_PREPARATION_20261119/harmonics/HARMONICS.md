# Scalar harmonics during closed hilltop residence

A finite band of scalar fluctuations on the exact stationary solution reaches the scale of the nonlinear hilltop force well before each of the six archived, homogeneous, critical launches leaves its linear matching region. The result depends on which state is prepared and which spatial scales are measured. Even the nine real quadrupole components alone reach the force benchmark before the homogeneous matching time in every case. [Numerical results](HARMONICS.json) · [Transfer matrices](TRANSFER.csv)

## Field and geometry

Use Einstein-frame proper time $x=H_{\rm top}\tau$, dimensionless radius $a=H_{\rm top}a_{\rm phys}=\cosh(hx)/h$, and $h^2=V_*/V_{\rm top}$. The stationary solution has $\dot\Theta=V_{,\Theta}=0$. A first-order field perturbation therefore has vanishing matter density and momentum sources for the linear Einstein constraints, and $\delta\Theta$ is invariant under an infinitesimal change of time slicing. The term proportional to the stationary potential in the metric equations is already included in the de Sitter background. The scalar obeys the minimally coupled Klein–Gordon equation on this metric.

Expand $\delta\Theta=\sum q_{n\alpha}Y_{n\alpha}$ with real harmonics orthonormal on the unit $S^3$. The eigenvalue of $-\Delta$ is $k_n^2=n(n+2)$, the degeneracy is $(n+1)^2$, and their summed square at any point is $(n+1)^2/(2\pi^2)$. Each coefficient obeys

$$
q_n''+3h\tanh(hx)q_n'
 +\left[k_n^2h^2\operatorname{sech}^2(hx)-\mu^2\right]q_n=0,
\qquad \mu^2=-3(V_E/V_{\rm top})_{,\Theta\Theta}|_* .
$$

In order of $\phi_m=+0.5,-0.5$, $h$ equals one to input precision and $\mu^2$ equals 207.2548854 and 140.9655115. The highest shell with strictly negative frequency squared at the bounce is 13 and 10, respectively. A higher shell first attains zero frequency squared at $x=h^{-1}\operatorname{arcosh}(hk_n/\mu)$. The expanding solution eventually brings each fixed shell into the band. For $n=16$ the crossing times are 0.58944 and 0.89635. These crossing times refer to the sign of the restoring term for the field $q$, which has a friction term. A change of canonical dependent variable can shift the sign change of its effective frequency.

The primary band spans $2\le n\le16$. The computation includes $n=0$ as a homogeneous transfer control and $n=1$ as a separate four-component field multiplet. The $n=1$ scalar harmonics have vanishing traceless Hessian on $S^3$, so the metric's scalar decomposition has an exceptional dipole sector. At an exactly stationary field the decoupled matter dipole still satisfies the displayed equation. Its interpretation on a rolling, constrained background requires a separate gauge and constraint analysis. A finite Gaussian matter state on a fixed sphere specifies neither the second-order constraints nor the global linearization-stability conditions of compact de Sitter gravity.

## Transfer for both initial data

Let $s=\tanh(hx)$, $D=\sqrt{9+4\mu^2/h^2}$, $L=(D-3)/2$, $a_n=n/2+(3+D)/4$, $b_n=n/2+(3-D)/4$. The solutions with initial $(q,q')=(1,0)$ and $(0,1)$, respectively, are

$$
U_j(x)=\cosh^L(hx)(s/h)^j
 {}_2F_1\!\left(\tfrac12-a_n+\tfrac j2,
 b_n+\tfrac j2;\tfrac12+j;s^2\right),\qquad j=0,1.
$$

This follows by substituting $z=-\sinh^2(hx)$, factoring $(1-z)^{n/2}$, and solving the resulting Gauss equation about $z=0$. Its second solution carries the odd factor $\sinh(hx)/h$. A Pfaff transformation gives the displayed argument $s^2$ and avoids hypergeometric evaluation at large negative arguments. Differentiating $U_j$ gives both rows of the transfer matrix. Abel's identity gives $\det T(x,0)=\cosh^{-3}(hx)$, preserving the canonical commutator as the physical sphere expands.

The table contains all four entries for shells 0 through 16 at the tachyonic crossings, common sample times, and the three archived residence times on each branch. Late-time entries can be too nearly proportional for their rounded products to recover the determinant in double precision. The final column reports the analytic determinant. A direct integration of both columns agrees to $1.83\times10^{-11}$ in the largest absolute entry error divided by $\max(1,|T_{ij}|)$. The analytic-column determinants on $0\le x\le0.3$ agree with Abel's identity to $1.8\times10^{-12}$ relative error. ODE evolution from the identity for intervals of length 0.07 starting at every tabulated time gives maximum relative local-determinant error $1.65\times10^{-13}$. This latter check also tests friction near the late residence endpoints. [Transfer table](TRANSFER.csv) · [Validation records](HARMONICS.json)

## Finite Gaussian preparation

The dimensionless canonical field is $\Theta=\varphi/M_{\rm Pl}$. With $v=q'$ the quadratic action gives $p=(M_{\rm Pl}/H_{\rm top})^2a^3v$ per orthonormal real harmonic. Set $\hbar=1$, choose zero means and prepare an independent pure reference-oscillator Gaussian in every retained real component at the bounce. Its positive reference mass is $m_{\rm ref}=H_{\rm top}$. With $\omega_n^2=h^2 k_n^2+1$ and $\kappa=(H_{\rm top}/M_{\rm Pl})^2 h^3$, its symmetrized covariance is

$$
\Sigma_{n,0}=\frac{\kappa}{2}\operatorname{diag}(\omega_n^{-1},\omega_n),
\quad [q,v]_0=i\kappa,\qquad \Sigma_n(x)=T_n\Sigma_{n,0}T_n^\mathsf{T}.
$$

The preparation saturates the canonical uncertainty relation. The chosen positive oscillator specifies its width at the bounce even for shells whose actual squared frequency is negative. The $n=0$ harmonic is $1/\sqrt{2\pi^2}$, so the homogeneous field coordinate has bounce commutator $i\kappa/(2\pi^2)$. The comparisons below use only shells 2 through 16. The $n=1$ control in the data uses the same preparation on that shell separately.

For an angular filter $W_n$, the field variance at a point equals
$$
\sigma_W^2(x)=\sum_n\frac{(n+1)^2}{2\pi^2}W_n^2(\Sigma_n)_{qq}.
$$

Isotropy makes the same value the expectation of the squared field averaged over the whole sphere. We evaluate the full band, $2\le n\le8$, the single $n=2$ shell, and two angular heat filters $W_n=\exp[-\ell^2 k_n^2/2]$ with $\ell=0.25,0.5$ on the full band. These filters have fixed angular scale and physical smoothing length proportional to $\ell a_{\rm phys}(x)$. At $t=10^{-3}$ the initial RMS values for full band, $2\ldots8$, quadrupole, $\ell=0.25$, and $\ell=0.5$ are respectively 0.01718, 0.009091, 0.002430, 0.003369, and 0.001005. They scale as $\sqrt t$ at the bounce. Dynamic variances and shell power fractions are tabulated at the bounce, halfway to the first full-band diagnostic, and at that diagnostic. The preparation is a finite-mode model. Inclusion of arbitrarily high modes would require a field state and a prescription for the composite operators such as stress energy.

## Comparison with homogeneous departure

On each side of the reconstructed maximum we find the displacement at which the exact local force first differs from its linear approximation by 10%. In negative and positive field directions, respectively, the radii are (0.0198263, 0.0185957) on the positive-minimum branch and (0.209742, 0.0786871) on the negative-minimum branch. The force marker uses the smaller of these two radii. Independently we monitor the magnitudes of three expectation values of the quadratic energy of the unfiltered finite band, normalized by $V_*=3M_{\rm Pl}^2H_{\rm top}^2 h^2$,
$$
\frac{1}{6h^2}\left\{
\langle (\delta\Theta')^2\rangle,
\ a^{-2}\langle|\nabla_\Omega\delta\Theta|^2\rangle,
\ -\mu^2\langle\delta\Theta^2\rangle\right\}.
$$

The energy marker occurs when any term reaches 0.1 in magnitude. Keeping the signs separately displays cancellations between kinetic, gradient and potential corrections. These are auxiliary smallness checks. No stress energy was fed back into Einstein's equations. The archived homogeneous critical history reaches $|\Theta-\Theta_*|=10^{-5}$ at $x_{\rm match}$, its linear-to-nonlinear handoff. This provides a reproducible comparison time for the specified detuning and decay model.

| $\phi_m$ | $t$ | $x_{\rm match}$ | Full-band force | $n=2$ force | $\ell=0.5$ force | Full-band energy |
|---:|---:|---:|---:|---:|---:|---:|
| +0.5 | $10^{-3}$ | 3.382 | 0.0285 | 0.1927 | 0.2574 | 0.1560 |
| +0.5 | $10^{-10}$ | 8.755 | 0.7310 | 0.7872 | 0.8561 | 0.8134 |
| +0.5 | $10^{-20}$ | 16.430 | 1.6200 | 1.6636 | 1.7347 | 1.7019 |
| −0.5 | $10^{-3}$ | 3.028 | 0.2587 | 0.3649 | 0.4468 | 0.2359 |
| −0.5 | $10^{-10}$ | 8.401 | 1.0573 | 1.1101 | 1.1976 | 1.0393 |
| −0.5 | $10^{-20}$ | 16.077 | 2.1588 | 2.2035 | 2.2926 | 2.1408 |

Times after the earlier full-band diagnostic quantify only the extrapolated linear transfer for that preparation. A filtered measurement can have small variance while shorter scales already explore a nonlinear region. As a different preparation, the quadrupole results describe the field of a truncation containing that shell alone. Neither benchmark establishes a sharp validity time for a Gaussian with unbounded support. For example, in the positive-minimum case at $t=10^{-3}$ the bounce preparation already assigns weight 0.264 to a field measurement at one point outside the two force radii. The corresponding negative-minimum weight is $2.32\times10^{-6}$. The records give these one-point Gaussian tails in log space. Paths through the nonlinear region, correlations between field measurements at different times, the backreacted geometry, and probabilities of eventual radiation equality require additional dynamics.

## Longest quadrupole residence at the measured force radius

One can optimize the initial quadrupole state over all preparations with finite second moments. For the expectation of the spatially averaged squared field, the smallest possible maximum on $[0,X]$ is $9\kappa B_2(X)/(4\pi^2)$, where $B_2=T_{q\leftarrow v,n=2}$. Each real component has endpoint commutator $[q(0),q(X)]=i\kappa B_2(X)$. The uncertainty relation gives a lower bound $\kappa B_2(X)$ on the sum of the two endpoint variances per component. An independent pure Gaussian per component with equal endpoint variances and zero symmetrized covariance between endpoints attains the bound. Every intermediate field is a linear combination of the endpoint fields with nonnegative weights of sum at most one. This follows from the one-dimensional maximum principle because the quadrupole's restoring coefficient is strictly negative throughout. A state attaining the endpoint bound consequently satisfies the same variance budget throughout the interval. Spatial averaging supplies the unit-volume factor $1/(2\pi^2)$.

Setting this optimum equal to the square of the smaller measured 10% force radius gives limiting $X$ values **0.45590, 1.67556, 3.44807** on the positive branch and **0.81076, 2.33785, 4.53702** on the negative branch, in order $t=10^{-3},10^{-10},10^{-20}$. Each is earlier than the archived homogeneous matching time. This calculation allows unrestricted initial squeezing and measures RMS over the sphere and quantum state. Pointwise field bounds, trajectory probabilities and a physical preparation with a prescribed energy are separate questions. [Measured radii and optimized times](HARMONICS.json)

## Reproduction

Use Python with NumPy and SciPy from the study directory. The program reads the pinned prior hilltop data and reconstructs the same spline potential for the force comparison.

~~~sh
OPENBLAS_NUM_THREADS=1 python -IB harmonics/harmonics.py --prior prior --output-dir /tmp/hdblast-harmonics
~~~
