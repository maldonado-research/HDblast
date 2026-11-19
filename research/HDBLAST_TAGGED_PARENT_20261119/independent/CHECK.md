# Switched Dirac-pair response of one added zero-momentum parent

This private calculation evaluates a proposed leading-order tagged-parent probe on the published fine $K=8$, $y=1$, feedback-on wall histories. It takes the wall as an external input. The initial zero-momentum test mode is the unit-Wronskian physical oscillator ground at $T=0.70$. An added parent is compared with the vacuum of this same mode on the same wall. Here $y$ in the denominator denotes the perturbing Yukawa coupling of the probe. The original bath coupling used to generate the prescribed wall is a separate model input.

Use proper time $s=H_0\tau$, $a=1$ at the original quiet crossing, $d\eta=ds/a$, $m=100|\phi-0.5|$ and $v=a\chi$. For the declared switch $f(s)=\sin^4[\pi(s-s_l)/(s_r-s_l)]$ on its window and zero outside, numerically evaluate the specified coefficients
$$
\Delta N_{\mathrm{pair}}/y^2=(4\pi)^{-1}\int |(fv)'|^2 d\eta,\qquad
\Delta E_{\mathrm{pair}}/(y^2H_0)=(8\pi^2a_{\mathrm{cut}})^{-1}\int_{\mathbb R}|\omega|^3|F(fv)(\omega)|^2d\omega
$$
with $F(g)(\omega)=\int d\eta\,g(\eta)e^{-i\omega\eta}$ and the common physical readout $s_{\rm cut}=2.10$. These expressions include both frequency signs. In the constant-scale-factor, constant-mass limit without switching transients the number rate reduces to $y^2m/(8\pi)$ for the normalized positive-frequency rest mode.

| Proper-time switch | Quiet number coefficient | Near-positive number coefficient | Number change | Quiet energy coefficient | Near-positive energy coefficient | Energy change |
|---|---:|---:|---:|---:|---:|---:|
| 1.15–2.05 | 0.49414018 | 0.49380704 | −0.0674180% | 10.5285293 | 10.5214515 | −0.0672251% |
| 1.25–1.95 | 0.40835495 | 0.40806420 | −0.0712015% | 8.2529372 | 8.2451276 | −0.0946280% |

The shorter switch changes the absolute outcome and the fractional energy contrast. For the wider quiet run, the Fourier contributions to pair number are 0.33919269 at $\omega<0$ and 0.15494751 at $\omega>0$. For the narrower quiet case these contributions are 0.27767450 and 0.13068047. One frequency half alone omits a substantial component of the specified inclusive coefficient.

The solver propagates $(\chi,\Pi=a^3\dot\chi)$ with exact symplectic oscillators on 4,096 through 32,768 midpoint-held conformal-time intervals. It reconstructs the proper/conformal clock mapping separately on 32,769 proper-time nodes and interpolates the measured $A$ and $\phi$ with their saved derivatives. It obtains $(fv)'=a\dot f\,v+f(a^2h\chi+\Pi/a)$ without numerical waveform differentiation. Direct time integration is checked against FFT Parseval integration. The largest absolute number disagreement at the finest step is below $1.7\times10^{-8}$. Doubling time resolution changes every pulse contrast by less than $4.3\times10^{-8}$ percentage points. At fixed finest resolution increasing the zero-padding factor from 16 to 32 changes every energy by less than $2.4\times10^{-9}H_0$ per $y^2$. The calculation samples an unbounded daughter continuum; the energy fraction above $|\omega|=2a(s_l)M_5/H_0$ is under 0.00010 in the four finite calculations. This conservative threshold is for individual massless-daughter energies on the switch support. It does not bound ultraviolet completion errors in this process.

These numbers condition on a supplied tree-level formula, parent preparation, background, switching function and detector time. Applying them at strong coupling would require higher orders, and the switch can supply work. They describe an extra, tagged parent. They do not quantify the production, thermalization or cosmological radiation density of a full coupled gas.

~~~sh
OPENBLAS_NUM_THREADS=1 python -I -B tagged_transfer.py /path/to/HDblast-repository > replay.json
~~~
