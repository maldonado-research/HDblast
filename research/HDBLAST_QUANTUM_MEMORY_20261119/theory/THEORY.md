# What a finite positive bath can predict

A force that varies slowly can be represented accurately by a few oscillators even when their short-time force correlations differ considerably from those of the continuum. Here we can bound both errors and show how the answer depends on the time window.

## Static calibration and the imaginary-frequency test

Consider a single fixed parent momentum $0<q<\ell$ and a conformal-frequency band $[\ell,H]$. In this experiment $\ell=60$ and $H=120$. With $\hbar=1$, the massless vacuum pair spectrum used by the model gives the positive squared-coupling measure
\[
d\gamma(\omega)=\omega(\omega^2-q^2)d\omega/(4\pi^2).
\]
Push $d\gamma/\omega^2$ forward to $x=\omega^{-2}$, obtaining $d\mu=(1-q^2x)dx/(8\pi^2 x^2)$ on $[a,b]=[H^{-2},\ell^{-2}]$. Its mass $C=\int d\mu$ is the static response. Its first moment $D=\int x\,d\mu$ is the magnitude of the finite-band slope. Positive $N$-point Gaussian quadrature $Q_N$ in $x$ exactly integrates degrees $0,\ldots,2N-1$. Thus even one node matches both $C$ and $D$. Each node gives an oscillator with real frequency $\omega_i=x_i^{-1/2}$ and squared coupling $g_i^2=\mu_i/x_i>0$.

Let $z=\nu^2\ge0$ denote squared imaginary frequency. The finite-band response is $W(z)=\int d\mu/(1+zx)$ and $W_N(z)=Q_N[1/(1+zx)]$. Its first $2N$ Taylor coefficients at zero match. There is also a global result. Write $m=(a+b)/2$, $h=(b-a)/2$. For $z>0$,
\[
0\le\frac{W(z)-W_N(z)}{W(z)}
\le {\rm sech}^2\!\left[N\,{\rm arcosh}\frac{1/z+m}{h}\right]
\le 4\left(\frac{H-\ell}{H+\ell}\right)^{2N}.
\]
The last factor equals $4/3^{2N}$ for this band. In particular the uniform relative shortfall of the four-node high-band response is at most $0.0610\%$ on the entire positive imaginary-frequency axis. The bound concerns the high-band response before taking differences against calibrated contact forces. A relative error normalized by a small residual can be much larger. The sign for a parent inverse response with bath correction $-y^2W$ reverses.

To prove the bound, let $p_N$ be the monic orthogonal polynomial for $d\mu$. Its roots are the quadrature nodes. For any polynomial $v$ of degree at most $N-1$, minimize $\int[1-(1+zx)v(x)]^2d\mu/(1+zx)$. At the minimizer the residual is orthogonal to every polynomial of degree $<N$, and it equals one at $x=-1/z$. Hence it is $p_N(x)/p_N(-1/z)$. Exact polynomial quadrature, or completing the square in this minimization, shows its minimum is $W-W_N$. In particular the sequence $W_N$ is increasing and bounded above by $W$, strictly below it for any continuous positive band at $z>0$. As a trial residual use $T_N((x-m)/h)/T_N((-1/z-m)/h)$, where $T_N$ is the Chebyshev polynomial. Its magnitude is bounded by the reciprocal of the denominator on the band. This proves the stated relative bound. At $z=0$ the two responses coincide exactly.

## Vacuum fluctuations and the real-time error

For the stationary, initially free vacuum bath define
\[
E(t)=\int d\mu(x)\,x^{-1/2}e^{-it/\sqrt{x}},\qquad
E_N(t)=\sum_i\mu_i\omega_i e^{-i\omega_it},
\]
where $t$ is elapsed **conformal** time. For a coupling $-uB$ the retarded kernel is $-\theta(t)\operatorname{Im}E(t)$, the symmetric force noise is $\operatorname{Re}E(t)/2$, and the ordered vacuum two-point function $\langle B(t)B(0)\rangle$ is $E(t)/2$. For example, in the Fourier convention $\int dt\,e^{i\Omega t}$, the noise spectrum equals $\operatorname{sgn}(\Omega)\operatorname{Im}\chi_R(\Omega)$ in these units. Each positive oscillator separately satisfies this vacuum fluctuation–dissipation relation and supplies a positive ordered spectral measure. A stationary thermal occupation multiplies the symmetrized spectrum by $\coth(\beta|\Omega|/2)$. The numerical bounds below use the vacuum state. A deterministic calculation of mean retarded force alone omits its vacuum fluctuations. See the [independent ordered energy calculation](../quantum/SPECTRUM.md).

An explicit finite-time bound follows from analyticity. For the dimensionless variable $y=\ell^2 x\in[(\ell/H)^2,1]$, let $c,d$ be the midpoint and half-width, and parameterize a Bernstein ellipse by
\[
y(\theta)=c+\tfrac d2[(r+r^{-1})\cos\theta+i(r-r^{-1})\sin\theta].
\]
Choose $1<r<(H+\ell)/(H-\ell)$ so its interior has positive real part. Set
\[
M_r(T)=\max_{\theta}\frac{\ell}{H} |y(\theta)|^{-1/2}
               \exp\!\bigl(T\ell\,|\operatorname{Im}y(\theta)^{-1/2}|\bigr).
\]
Truncating the Chebyshev series at degree $2N-1$ gives a uniform polynomial error no larger than $2M_r(T)r^{-(2N-1)}/(r-1)$, since the degree-$k$ coefficient is bounded by $2M_r(T)r^{-k}$. The continuum and discrete integration functionals are positive with common mass $C$. Both integrate this polynomial exactly. Therefore
\[
\sup_{0\le t\le T}\frac{|E_N(t)-E(t)|}{CH}
\le\min\!\left(2,\,\frac{4M_r(T)r^{-(2N-1)}}{r-1}\right).
\]
This is an absolute error normalized by a common scale. It remains meaningful at zeros of the retarded or symmetric kernel. Taking the minimum over admissible ellipses is valid. The executable bounds the maximum on every ellipse arc. If $\operatorname{Re}y\ge u>0$ and $|\operatorname{Im}y|\le v$, then $|y|^{-1/2}\le u^{-1/2}$ and $|\operatorname{Im}y^{-1/2}|\le\min((8u)^{-1/2},v/(2u^{3/2}))$.

Sampling also gives a bound for the intervening times. With spacing $\delta$, linear interpolation of the complex difference differs by at most $\delta^2\sup |(E_N-E)''|/8$. The latter derivative is bounded by $2CH^3$, so adding $(H\delta)^2/4$ to the largest normalized sample error covers the grid gaps. This bound assumes exact sample values. The archived samples use double precision.

At the largest native momentum, the four-node relative error at $\nu=60$ is $4.3292\times10^{-8}$. Over $0\le t\le0.70$, its largest sampled normalized complex-kernel error is **0.68361**. The retarded component alone reaches **0.68347** on the same scale. The interpolation allowance is 0.00090. Sixteen inverse-frequency nodes reduce the maximum sampled complex error to **0.02590**. With 32 such nodes the sampled maximum is below $3.5\times10^{-11}$, and the analytic all-times-on-this-interval bound is $6.44\times10^{-7}$. The 32-node result here uses a Gaussian rule in inverse squared frequency. The full shell benchmark independently integrates its reference band in frequency. These figures describe pointwise kernels. A smooth forcing waveform can average many of their oscillations. [Kernel values and error bounds](CHECK.json)

More precisely, a normalized kernel bound $\epsilon$ gives $|\delta F(t)|\le CH\epsilon\,\|u\|_1$ for the induced force when the elapsed time since the drive began is at most $T$. For a smooth pulse with zero endpoint values, $|\delta W_{\rm drive}|\le CH\epsilon\,\|\dot u\|_1\|u\|_1$ follows from $W_{\rm drive}=-\int\dot uF\,dt$. A noise measurement smeared by a real test function $f$ has variance error at most $CH\epsilon\,\|f\|_1^2/2$ on the same interval. These are conservative upper bounds. None of the three uses a relative error near a spectral zero.

### Why a time horizon matters

The continuous high band has an integrable density, so $E(t)\to0$ by the Riemann–Lebesgue lemma. For every fixed finite rule, simultaneous phase recurrence gives arbitrarily late times $t_k$ with $\exp(-i\omega_it_k)\to1$ for all its nodes. Consequently $E_N(t_k)\to E_N(0)>0$. At $t_k+\tau$ the same argument returns $E_N(\tau)$, including its retarded part when the sine sum is nonzero. Long-time uniform convergence to this continuum is impossible for a fixed finite bath. The positive-frequency-axis bound and finite-lag error estimate have different physical content.

On the moving shell, elapsed conformal time is $\int ds/a(s)$ for proper time $s$. The factors from the scale factor and scalar gate belong in the source multiplying the conformal kernel. Energy, pressure and scalar force require the Hamiltonian and its variations, including the chosen compensation and the initial state. Equal static and kinetic moments do not independently fix those stresses, reactions after switching or the pair-energy reservoir. The finite Gaussian bath represents the specified vacuum two-point channel through order $y^2$; its higher correlation functions need not equal those of an interacting fermion theory. A test of spacetime dimension or cosmological reheating needs additional observables.

## Independent numerical construction

The check calculates moments from elementary powers and logarithms of $y$. It constructs the monic orthogonal polynomials with 180–320 digit decimal arithmetic. A Jacobi continued fraction tests the sign and the bound before converting the Gaussian nodes to machine precision. The continuum time reference is the closed Fourier transform of a quadratic polynomial on the band, evaluated with three spherical Bessel functions and separately checked by direct frequency integration. The actual source nodes and squared couplings are compared against the high-precision moment construction for every native momentum and every compressed order through 16. This procedure uses no shell-trajectory fit. [Reproduction code](check.py) · [Compressed rules](../source/RULES.json)
