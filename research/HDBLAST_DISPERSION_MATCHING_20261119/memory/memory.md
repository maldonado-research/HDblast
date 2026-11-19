# Testing the retarded force with two positive pulses

A response to a force can be measured after the force is switched off. At that time, local terms proportional to the drive and its derivatives vanish. The remaining response probes the spectrum of the bath. We consider one prescribed drive and evaluate the mean force in linear response. It is a diagnostic of the Gaussian bath used in HDblast. The simulated shell and its energy density are examined separately in the gravity calculation. [Test specified before evaluation](https://github.com/maldonado-research/HDblast/issues/157#issuecomment-7001918750)

For a native momentum $q$, a change to conformal time $d\eta=ds/a$ gives a bath mode with frequency $\omega$ and an applied source proportional to $aF\chi$. Set the proportionality and coupling factors to one and define the reported response with the positive sign of the retarded sine propagator. The reaction coordinate in a Hamiltonian with an interaction $+u\sum_j g_jx_j$ carries the opposite sign. All comparisons below use the same sign. In conformal variables the published spectral measure is $g^2dN=\omega(\omega^2-q^2)d\omega/(4\pi^2)$ for $\omega>q$. For a positive lag $t$ the causal continuum kernel is

$$
K_q(t)=\lim_{\epsilon\downarrow0}\int_q^\infty
 e^{-\epsilon\omega}\frac{\omega^2-q^2}{4\pi^2}\sin(\omega t)\,d\omega
=-\frac{qt\sin(qt)+\cos(qt)}{2\pi^2t^3}.
$$

Shifting the integral by $q$ reduces it to the integrals of $x$ and $x^2$ times $e^{(it-\epsilon)x}$. Polynomial contact terms are supported at zero lag. We only convolve the displayed kernel at strictly positive lag, after each pulse ends. The same response follows by integrating the centered pulse Fourier transform against the spectral weight. [Independent frequency calculation](pulses.py)

A pulse centered at zero is $u(t)=\cos^4(\pi t/T)$ for $|t|<T/2$, zero elsewhere. We selected the native momentum closest to two, $q=1.818314$, and the width $T=0.642001$ from the integrated conformal duration of the fine quiet switching interval. The external force is prescribed in conformal time. Only this duration and momentum are taken from the shell history. Its late readout times extend beyond that simulation. Let $A_q=\int u(t)\cos(qt)dt=0.237488$. The single-pulse late response has leading term $-qA_q\sin(qv)/(2\pi^2v^2)$.

We compare two positive pulses separated by $d=\pi/q$ with the same two pulses separated by $2d$. Together they supply equal integrated drive in both experiments. Their combined responses at time $v$ after the second center are $B(v)+B(v+d)$ and $B(v)+B(v+2d)$, with $B$ the single-pulse response. The half-period separation cancels the leading threshold oscillation. On the phase sequence $qv=2\pi n+\pi/2$, the half- and full-period responses have leading powers $v^{-3}$ and $v^{-2}$. Their ratio tends to $d/v$. Doubling the delay should then approximately halve the ratio.

| Phase index $n$ | Half-period response / full-period response |
|---:|---:|
| 4 | 0.120513 |
| 8 | 0.061837 |
| 16 | 0.031160 |
| 32 | 0.015613 |
| 64 | 0.007811 |

A sharp finite cutoff introduces its own endpoint oscillation. Integration by parts gives a contribution proportional to $-(\Lambda^2-q^2)\widetilde u(\Lambda)\cos(\Lambda v)/v$ for a single pulse, apart from the common spectral normalization. At fixed cutoff this can eventually exceed the continuum threshold signal. At cutoff 60 and $n=64$ the full-period response with a sharp edge is about 13.5 times the continuum value. A smooth cutoff that rolls to zero with two continuous derivatives over the top fifth of the band reproduces the continuum at the tested phases with increasing accuracy. At cutoff 480 the maximum relative discrepancy over the three pulse choices and five readouts is below $10^{-7}$. We also compared discrete quadratures to the sharp-band integral. The original 64-node rule performs poorly for these long oscillatory delays. The 4096-node rule agrees at both independent quadrature readouts to within $10^{-7}$. [Results for every cutoff and phase](PULSES.json)

This gives a concrete way to test the assumed continuum threshold. A residual measured after a double pulse would have to exhibit the predicted phase and decay with controlled band edges. The test could also succeed in an ordinary four-dimensional oscillator bath with the same spectral density. Its cosmological application requires a physical source for the probes and a measurable signal.
