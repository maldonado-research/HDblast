# What does time averaging select on the round de Sitter bounce?

A finite time resolution selects a regular state of the nonzero massless scalar waves. Its physical energy density averaged through the bounce falls below that of the exact comparison state. Its density exceeds the comparison density later in the expansion. Its pressure at the bounce depends appreciably on how long energy was sampled.

We use one real, minimally coupled, massless scalar with constant kinetic normalization set to one, $\hbar=c=1$, $a(x)=\cosh(hx)/h$, spatial volume $2\pi^2a^3$, and unit-sphere harmonics normalized to one. The physical scale $h$ is constant. Let $z=hx$, $\theta=hT$, $\eta=\arctan(\sinh z)$ and $c=\cos\eta$. The nonnegative weight $w_\theta(z)$ is proportional to $\exp[-1/(1-(z/\theta)^2)]$ for $|z|<\theta$, vanishes outside, and integrates to one in $dz$. Olbermann's smooth sampling function is the square root of the weight, with the corresponding Jacobian for proper time. We average the local physical density along a comoving observer. Spatial homogeneity gives the same result everywhere.

| Proper-time halfwidth $hT$ | Sampled density difference | Bounce density difference | Bounce pressure difference | Density difference at $hx=0.5$ |
|---:|---:|---:|---:|---:|
| 0.1 | −0.01154834 | −0.00847271 | +2.23898458 | +0.00777988 |
| 0.2 | −0.00925346 | −0.00631094 | +0.56519138 | +0.00660532 |
| 0.3 | −0.00783363 | −0.00508158 | +0.25547561 | +0.00619572 |

Each entry gives the time-smeared state's value minus the exact reference state's value, in units of $h^4$. The zero mode is held in the same state on both sides. The negative entries are energy differences. Absolute renormalized densities also require the reference state's renormalized stress and a choice of gravitational counterterms.

## Minimization and short wavelengths

Each label $n\ge1$ has $L^2=(n+1)^2$ real oscillators with $\lambda=L^2-1$ and proper frequency $\omega=\sqrt{\lambda}/a$. The normalized reference solution and its proper-time derivative are

$$
u_L=\frac{h(Lc+i\sin\eta)e^{-iL\eta}}{\sqrt{2L\lambda}},\qquad
\dot u_L=-\frac{ih^2\lambda c^2e^{-iL\eta}}{\sqrt{2L\lambda}}.
$$

The conserved normalization is $ia^3(u_L^*\dot u_L-\dot u_L^*u_L)=1$. Write $\langle g\rangle_w=\int w_\theta(z)g\,dz$ and define the two quadratic forms

$$
A_L=\tfrac12\langle |\dot u_L|^2+\omega^2|u_L|^2\rangle_w,\qquad
B_L=\tfrac12\langle \dot u_L^2+\omega^2u_L^2\rangle_w .
$$

A normalized mode $v_L=\alpha_Lu_L+\beta_Lu_L^*$ has $\alpha_L>0$ and $\alpha_L^2-|\beta_L|^2=1$ after choosing its arbitrary phase. Its sampled density per oscillator, before the spatial factor, equals $(1+2|\beta_L|^2)A_L+2\alpha_L\operatorname{Re}(\beta_L^*B_L)$. Positive kinetic and gradient energy together with the Wronskian give $A_L>|B_L|$. The minimum is $D_L=\sqrt{A_L^2-|B_L|^2}$, attained at $\beta_L/\alpha_L=-B_L/(A_L+D_L)$. The gain $D_L-A_L=-|B_L|^2/(A_L+D_L)$ avoids catastrophic cancellation for high harmonics. In the minimizing basis the anomalous form vanishes. Any further normalized squeeze with negative-frequency coefficient $b$ raises the objective by $2D_L|b|^2$. Occupying its oscillator also raises the energy. For the even sampling used here $B_L$ and $\beta_L$ are real, $v_L(-x)=v_L(x)^*$ and the resulting diagonal stress is even.

Direct substitution gives

$$
A_L=\frac{h^4}{4L}\big\langle 2L^2c^4-c^2\cos(2\eta)\big\rangle_w,\qquad
B_L=\frac{h^4}{4L}\big\langle c^2[\cos(2\eta)+iL\sin(2\eta)]e^{-2iL\eta}\big\rangle_w .
$$

Set $F(\eta)=w_\theta(z)c$ and use $dz=d\eta/c$. Integration by parts gives

$$
B_L=\frac{h^4}{4L}\int
\left[2F\cos(2\eta)+\tfrac12 F' \sin(2\eta)\right]e^{-2iL\eta}\,d\eta .
$$

All boundary terms vanish. The bracket is smooth with compact support. Any further $k$ integrations by parts give $|B_L|\le h^4 C_{k,\theta}L^{-1-k}$, while $A_L/L\to h^4\langle c^4\rangle_w/2>0$. Thus $|\beta_L|\le C'_{k,\theta}L^{-2-k}$ for every $k$ at fixed nonzero $\theta$. The constants depend on the time resolution. The zero-width limit at fixed $L$ yields the instantaneous bounce ground mode, whose large-$L$ comparison with this reference has $|\beta_L|^2\sim 1/(16L^4)$. Exchanging these limits loses the short-distance control.

For completeness the reference rescaling $a u_L$ has exact conformal phase frequency $W=L\lambda/(L^2+\tan^2\eta)$ and squared amplitude $1/(2W)$. The phase frequency solves the adiabatic frequency equation for $L^2-2\sec^2\eta$. Its expansion about $L=\infty$ has positive leading term $L$ and can be differentiated term by term to every order on compact conformal-time intervals. It therefore has the all-orders positive-frequency expansion of a Hadamard state in this closed Robertson–Walker spacetime. The rapidly decreasing change of modes above adds a smooth kernel to the two-point function. Any regular Gaussian choice for the single remaining mode also changes the kernel smoothly. Together these facts give the same local Hadamard singularity. This construction uses the explicit round background. [Olbermann, *States of low energy on Robertson–Walker spacetimes*](https://arxiv.org/abs/0704.2986)

## The conserved stress difference

The sum over a complete unit-sphere shell gives $\sum_m |Y_{nm}|^2=L^2/(2\pi^2)$. For $d_k=|\dot v_L|^2-|\dot u_L|^2$ and $d_g=\omega^2(|v_L|^2-|u_L|^2)$ the isotropic differences are

$$
\Delta\rho=\sum_{L\ge2}\frac{L^2}{4\pi^2}(d_k+d_g),\qquad
\Delta p=\sum_{L\ge2}\frac{L^2}{4\pi^2}(d_k-d_g/3).
$$

The factors of $a$ already occur in the normalized modes, their proper derivatives and $\omega$. On a compact spacetime region the modes, harmonics and their derivatives of any fixed order grow at most polynomially in $L$. The bound on $\beta_L$ then yields absolute local convergence, including after any prescribed number of derivatives. The difference thus defines a smooth tensor and yields a distribution against smooth compactly supported spacetime tests. Each mode obeys the massless Klein–Gordon equation, giving $\dot{\Delta\rho}+3(\dot a/a)(\Delta\rho+\Delta p)=0$. A common local renormalization prescription cancels in the difference.

The program computes the oscillatory integral by Fourier transformation and independently checks the resulting coefficients by Gauss–Legendre quadrature through $n=4095$. Their largest absolute difference is $1.1\times10^{-15}$. In a third quadrature it minimizes the directly evaluated original kinetic and gradient energies for the first 24 shells and verifies stationarity in the full complex squeeze direction to $3.1\times10^{-16}$ in units of $h^4$. Other checks cover the time reflection, the Wronskian, and a finite-difference continuity residual below $4\times10^{-12}$ in units of $h^5$. On the reported time points, increasing the final shell from $8191$ to $16383$ changes either stress component by less than $10^{-12}h^4$. The ultraviolet conclusion follows from smooth smearing and integration by parts. Tiny individual high-harmonic coefficients can reach floating-point roundoff.

## What happens when the time window shrinks?

The width scan also admits two analytic limits. Write $g(y)=\theta w_\theta(\theta y)$ for the fixed unit-integral profile, $\widehat g(k)=\int g(y)e^{-iky}\,dy$, and $J(k)=\widehat g(k)-k\widehat g'(k)$. Then

$$
\frac{\Delta\rho(0)}{h^4}=-\frac{\log(1/\theta)}{32\pi^2}+O(1),
\qquad
\lim_{\theta\to0}\frac{\theta^2\Delta p(0)}{h^4}
=\frac{1}{16\pi^2}\int_0^\infty k\widehat g(k)\,dk
=0.0223264970150
$$

for the stated bump. Thus improving time resolution has very different effects on the bounce density and pressure.

For fixed $k=2L\theta>0$ as $\theta\to0$, substitute $z=\theta y$ in $B_L$. The factor $\cos(2\eta)+iL\sin(2\eta)$ tends to $1+iky$ and the oscillatory phase tends to $e^{-iky}$. Hence $4LB_L/h^4\to J(k)$ and $4L^2\beta_L\to-J(k)$. In the same scaled coordinate, $\eta(\theta y)/\theta$ has derivative $\operatorname{sech}(\theta y)$ bounded away from zero for $\theta\le0.3$. The smooth amplitude after the first integration by parts, including its measure, has bounded derivatives of every order uniformly in $\theta$ on $[-1,1]$. Subsequent integrations by parts in $y$ give the joint bound $|\beta_L|\le C_j L^{-2}(1+L\theta)^{-j}$ for each $j$.

At the bounce the exact pressure difference per shell divided by $h^4$ is $L[\beta_L^2(2L^2/3-1)-\alpha_L\beta_L(4L^2/3-1)]/(4\pi^2)$. Its leading term in the scaled sum is $LJ(2L\theta)/(12\pi^2)$. The joint bound controls the tails when taking the Riemann-sum limit, giving $(48\pi^2)^{-1}\int_0^\infty kJ(k)\,dk$. Integration by parts converts this to the pressure coefficient above. It can also be evaluated without an oscillatory integral using $\int_0^\infty k\widehat g(k)\,dk=-\int g'(y)/y\,dy$, with the value at zero defined by continuity.

The exact density difference per shell divided by $h^4$ is $L[\beta_L^2(2L^2-1)+\alpha_L\beta_L]/(4\pi^2)$. For $L\le 1/\theta$, expansion of the even integrand gives $\beta_L=-[1+O((L\theta)^2+L^{-2})]/(4L^2)$. This yields the shell term $-1/(32\pi^2 L)$ with an error bounded by a constant times $[(L\theta)^2+L^{-2}]/L$. The sum of those errors is bounded uniformly in $\theta$. The joint bound above also controls the sum for $L>1/\theta$ by a constant. This recovers the logarithm from the instantaneous launch state at the bounce, with the time resolution supplying the physical preparation scale.

For $\theta=0.003125$, summing through $n=255999$ gives $\theta^2\Delta p(0)/h^4=0.02232652331$. Halving the width from $0.00625$ to $0.003125$ recovers $0.9999853$ times the predicted logarithmic density slope. These divergences occur in the family of states as the preparation width tends to zero. Every finite width satisfies the smooth-state result above.

## The homogeneous degree of freedom

Write $\phi_0(x,\Omega)=q(x)/\sqrt{2\pi^2}$ and let $P=a^3\dot q$ with $[q,P]=i$. When the action has constant unit kinetic coefficient and zero mass, $P$ is conserved and the proper Hamiltonian is $P^2/(2a^3)$. The exact operator solution is

$$
q(x)=q(0)+I(x)P,\qquad
I(x)=\frac{h^2}{2}\left[\operatorname{sech}(hx)\tanh(hx)+\arctan(\sinh hx)\right].
$$

Its density and isotropic pressure coincide, $\rho_0=p_0=\langle P^2\rangle/(4\pi^2a^6)$. Under this evolution the covariance transforms by $q\mapsto q+IP$. In particular $\operatorname{Var}q(x)=\sigma_q^2+2I C+I^2\sigma_P^2$ for initial symmetrized covariance $C$. The determinant $\sigma_q^2\sigma_P^2-C^2\ge1/4$ is conserved. Every positive time weight minimizes the zero-mode energy by seeking $\langle P^2\rangle=0$. On the real noncompact field line this infimum is approached by broad normalizable packets and has no normalizable minimizer. Removing the constant mode by a shift quotient would specify a different algebra of observables.

At fixed initial standard deviation $\sigma_q$, a centered uncorrelated minimum-uncertainty Gaussian attains the smallest $\langle P^2\rangle=1/(4\sigma_q^2)$. It is time-reversal symmetric at the bounce. For this choice the separate peak energy density and its ratio to the round de Sitter source density are

$$
\rho_0(0)=\frac{h^6}{16\pi^2\sigma_q^2},\qquad
\frac{\rho_0(0)}{3h^2/(8\pi G)}=\frac{G h^4}{6\pi \sigma_q^2}.
$$

Here $h$, $q$ and $G$ must be expressed in mutually consistent physical units. The ratio decreases as $\operatorname{sech}^6(hx)$ along the fixed background. Giving this mode a fraction $\epsilon$ of the background density therefore requires $\sigma_q^2\ge G h^4/(6\pi\epsilon)$ for this Gaussian. For general packets the condition is $\langle P^2\rangle\le 3\pi\epsilon/(2G h^4)$, including the mean momentum. If a separate renormalized calculation bounds the absolute nonzero-mode energy density throughout the interval by $R$, a sufficient condition for $|\rho_{\rm total}|\le\epsilon\,3h^2/(8\pi G)$ is
$\sigma_q^2\ge h^6/[16\pi^2(\epsilon\,3h^2/(8\pi G)-R)]$ with a positive denominator and the specified Gaussian. The sampled differences in the table do not supply that absolute bound. A bound on energy density alone also leaves pressure and full semiclassical evolution to be tested.

From the study root, run

~~~sh
OPENBLAS_NUM_THREADS=1 python -I -B checks/state_check.py --output /tmp/STATE_CHECK.json
~~~

[Calculation](state_check.py) · [Numerical results and convergence tables](STATE_CHECK.json)
