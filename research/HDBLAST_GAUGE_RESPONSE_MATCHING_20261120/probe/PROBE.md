# Measuring the linear response of the added photon modes

When a scalar changes the electromagnetic coupling, short photon waves react even if they carry little energy after the coupling stops. To measure that response under controlled conditions, hold the radius of the closed universe fixed, prescribe a weak variation of the coupling and monitor the force from the additional waves. A comparison at two frequencies removes two different cutoff dependencies. It still leaves a predictable, frequency-dependent response and a finite memory of a short pulse.

## Retarded response on the static sphere

Use the conformal time of the unit three-sphere and define $s=(\ln I)'$. Each real transverse photon oscillator of label $L$ has conformal frequency $L$ and there are $2(L^2-1)$ such oscillators. Start all of them in the conformal vacuum. To linear order in the prescribed coupling variation the complex anomalous covariance $m_L=X_L+iY_L$ satisfies

$$ m'_L=-2iL m_L+s,\qquad X''_L+4L^2 X_L=s'. $$

The retarded solution is $X_L(\eta)=\int_{-\infty}^{\eta}s(v)\cos[2L(\eta-v)]dv$. The number in the mode begins at second order, $n_L=|m_L|^2+O(s^4)$ with the pulse shape held fixed while its amplitude varies. This follows either from the pure-state relation $n_L(n_L+1)=|m_L|^2$ or by integrating $n_L'=2s X_L$.

In the following sums the lower band, $L\le48$, is the reference. The added modes have $49\le L\le K$. Write $w_L=L(L^2-1)$, $d_L=L-1/L$ and $x=\omega^2/4$. With Fourier convention $e^{-i\omega\eta}$, define the dimensionless susceptibility of the added band by

$$ 4\sum_L w_L \widehat X_L=\chi_K(\omega)\widehat {s'},\qquad
\chi_K(\omega)=\sum_L \frac{4w_L}{4L^2-(\omega+i0)^2}
=\sum_L\frac{d_L}{1-x/L^2}. $$

The last expression applies below the first added pole, $|\omega|<98$. In this band a monochromatic linear response is real. At a pole the retarded prescription and the finite time of the drive matter. The poles occur at twice the individual mode frequency because pair creation couples the vacuum to two excitations. A finite pulse has a spread of frequencies and can leave real excitations.

In the rolling calculation's units, with $v=d\Theta/d(H_0t)$, $g=d\ln I/d\Theta$, $\epsilon=(H_0/M_{\rm Pl})^2$, the proper-time parent equation contains the force $-[2\epsilon/(\pi^2a^4)]g\sum_Lw_L X_L$. At fixed radius and locally constant $g$ one has $s'=a^2g\, dv/d(H_0t)$. Moving the quasistatic force to the inertial side adds $\Delta Z=\epsilon g^2\chi_K(0)/(2\pi^2 a^2)$ to the parent velocity's kinetic coefficient. The coefficient is positive. This response measurement calibrates the local, weak-field term at the declared radius.

## A second calibration and an independent held-out frequency

The first two terms of the expansion of $\chi_K$ are $A_K+B_Kx$, where $A_K=\sum d_L$ grows quadratically with the cutoff and $B_K=\sum d_L/L^2$ grows logarithmically. The remainder after fitting the value and slope at zero is

$$ R_K^{(0,0)}(x)=x^2\sum_L \frac{d_L}{L^2(L^2-x)}. $$

One can also fit the value at zero and at a nonzero probe frequency with $b=\omega_m^2/4$. Fit a line in $x$ to the two values. Subtracting that line leaves the exact residual

$$ R_K^{(0,b)}(x)=x(x-b)\sum_L \frac{d_L}{(L^2-x)(L^2-b)}. $$

For $0\le x,b<49^2$, the omitted tail above $K$ has magnitude at most

$$ \frac{|x(x-b)|}{2K^2[1-x/(K+1)^2][1-b/(K+1)^2]}. $$

This follows from $d_L\le L$ and $\sum_{L>K}L^{-3}\le1/(2K^2)$. The infinite sum can also be evaluated independently. With $a_*=49$, $R_\infty^{(0,0)}(x)=(x-1)[\psi(a_*)-(\psi(a_*-\sqrt{x})+\psi(a_*+\sqrt{x}))/2]+x\zeta(3,a_*)$. Here $\psi$ is the digamma function and $\zeta$ is the Hurwitz zeta function. Partial fractions of $L/(L^2-x)$ give the formula. Subtracting $(x/b)R_\infty^{(0,0)}(b)$ gives the limit for the two-frequency calibration.

Fit at $\omega=0,32$ and evaluate at $\omega=64$. The values below use the susceptibility normalization defined above. The last column bounds the additional absolute residual from all shells above the listed cutoff after the two-frequency calibration.

| Highest label $K$ | Response remaining after matching at zero | Response remaining after matching at 0 and 32 | Tail bound for the latter |
| ---: | ---: | ---: | ---: |
| 96 | 937.1574 | 188.1402 | 49.2164 |
| 192 | 1689.5306 | 222.7098 | 11.0441 |
| 384 | 2408.7601 | 230.8365 | 2.6899 |
| 1536 | 3830.6593 | 233.3412 | 0.1668 |
| Infinite band | Diverges logarithmically | 233.5078 | 0 |

The finite residual at the held-out frequency measures dispersion. The two-frequency fit extrapolates a weak, static-geometry response. A nonlinear implementation of an additional derivative term in the cosmology would require its own equations, matching prescription and stability analysis.

## Time-domain pulse and energy delivered

Independently impose $\ln I=\lambda\sin^6(\pi\eta/T)$ for $0<\eta<T=0.25$ and zero elsewhere. The extended coupling is $C^5$. It starts and ends at the same value. The pulse can be solved as an initial-value problem. Its cosine coefficients for $j=1,2,3$ are $\lambda(-15,6,-1)/32$, with fundamental $\nu=2\pi/T$. At unit amplitude, denote these three coefficients by $c_j$ and set $v_{Lj}=-c_j (j\nu)^2/[4L^2-(j\nu)^2]$. Inside the pulse the exact linear solution is

$$ X_L=\sum_j v_{Lj}[\cos(j\nu\eta)-\cos(2L\eta)]. $$

After switch-off, propagate the final $X_L,Y_L$ freely. The second cosine inside the switch enforces the initial conditions. It would be lost by treating the pulse as a sum of eternal monochromatic waves.

For $K=128$, the static term $A_Ks'$ alone has a **10.68%** relative root-mean-square force error during the pulse. Adding the slope term $-B_Ks'''/4$ reduces it to **2.55%**. These numbers compare the approximate and exact $4\sum w_LX_L$ on 601 evenly spaced points from start to switch-off. After switch-off the covariances continue to rotate, while the coupling derivative vanishes and energy transfer stops.

The outgoing coefficient of $a_{\rm physical}^4\rho_\gamma/\lambda^2$ to leading order is $D=\sum w_L |m_L(T)/\lambda|^2/\pi^2$. Integrating the exact retarded work, $\int s(4\sum w_L X_L)d\eta/(2\pi^2\lambda^2)$, gives the same value. Each of the two local derivative terms gives zero net work for this closed excursion by integration by parts. At $K=384$ the retarded solution gives $D=6.73442015$ for the added modes, with at most $2.23\times10^{-7}$ contributed by all higher modes at this order. To obtain this ceiling, write

$$ \frac{n_L(T)}{\lambda^2}=2\left[1-\cos(2LT)\right]\left(\frac{45\nu^6/4}
{[4L^2-\nu^2][4L^2-4\nu^2][4L^2-9\nu^2]}\right)^2 + O(\lambda^2). $$

Bounding the first factor by four, each denominator by $4L^2[1-(3\nu/(2(K+1)))^2]$, and $\sum_{L>K}L^{-9}$ by $1/(8K^8)$ gives the stated omitted-energy bound. It applies to this pulse and its weak-amplitude coefficient.

A second integration uses the full covariance equations in real canonical position and momentum, including the varying coefficients at finite amplitude. At $K=128$ and $\lambda=0.0002$ the relative root-mean-square discrepancy in the force is $1.99\times10^{-9}$ and the relative outgoing-energy discrepancy is $2.26\times10^{-7}$. The finite-radius probe is an external source experiment. Neither its prescribed amplitude nor its outgoing energy is a measurement of cosmological heating. [Executable probe](probe.py) · [Numerical results](PROBE.json)
