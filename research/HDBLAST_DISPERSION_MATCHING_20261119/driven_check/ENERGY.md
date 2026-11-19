# Work left in the high-frequency bath

The measured quiet-shell drive deposits **2.57397203** units of conformal energy in the modeled 60–120 band. This number refers to the response of an initially stationary, externally forced oscillator bath, with the coupling coefficient factored out. The two local approximations give exactly zero net work for the same drive. The energy calculation sets a lower bound of **1.063%** on the relative force error of any approximation with zero net work on these waveforms. [Fourier calculation](energy.py) · [Numerical checks](ENERGY.json) · [Direct oscillator evolution](../driven/DRIVEN.json)

## Energy from the measured waveform

Use conformal time $d\eta=ds/a$ and the native parent momentum weights $w_q$. For a real parent trajectory, define $u=aF\chi$ and let each initially stationary oscillator solve $Q_j''+\omega_j^2Q_j=u$. Primes denote conformal-time derivatives. The reaction convention is $R=\sum_jg_j^2Q_j$. At the end of the drive, the sine and cosine solutions of the oscillator equation give

$$
\frac{1}{2}\sum_j g_j^2[(Q_j')^2+\omega_j^2Q_j^2]
=\frac{1}{2}\sum_j g_j^2\left|\int u(\eta)e^{i\omega_j\eta}\,d\eta\right|^2 .
$$

The right side only requires the external drive and the spectral weights. In the native model, $g_j^2$ represents $\omega(\omega^2-q^2)d\omega/(4\pi^2)$. At every $q$, the verifier evolves the free parent twice, from canonical coordinates $(\chi,\Pi)=(1,0)$ and $(0,1)$, on the tabulated fine quiet wall. It uses $F(s)=\sin^4[\pi(s-1.35)/0.70]$ inside the gate. The scalar gate is calibrated on this same monotone quiet history. The integrand satisfies $u\,d\eta=F\chi\,ds$, which retains the scale factor implied by the conformal clock. If the two Fourier amplitudes are $v_\chi$ and $v_\Pi$, the covariance average of their squared combination is $X|v_\chi|^2+P|v_\Pi|^2+2C\operatorname{Re}(v_\chi^*v_\Pi)$. This uses the tabulated fork covariance directly. The verifier then sums the positive spectral power over the high-frequency band and all 24 weighted parent modes.

| Measurement | Conformal energy or work |
|---|---:|
| Fourier power, 160 time nodes and 96 frequency nodes | 2.5739720310 |
| Independent causal convolution, $-\int u'R\,d\eta$ | 2.5739720387 |
| Terminal oscillator energy at switch-off | 2.5739720333 |
| Oscillator energy at final readout | 2.5739720333 |
| Direct oscillator work, $\int uR'\,d\eta$ | 2.5739720343 |

The last three rows come from the independent time-domain solver. The first two use distinct integrals over the same external waveform. The Fourier value changes by less than $7\times10^{-11}$ in relative terms across three quadrature resolutions. The causal convolution difference across 48 and 80 nodes provides a separate estimate of its integration accuracy.

## A lower bound on force error

Write $\langle u',R\rangle$ for the conformal-time integral on the gate, summed with native momentum weights and averaged with the parent covariance. Since $u=0$ at each endpoint, the energy balance reads $E=\langle u,R'\rangle=-\langle u',R\rangle$. Each proposed local force is of the form $R_{\rm loc}=C_qu-D_qu''$ with real, time-independent coefficients for each mode. Its integrated work reduces to

$$
\int uR_{\rm loc}'\,d\eta=
\left[\frac{C_q}{2}u^2-D_q\left(u u''-\frac{(u')^2}{2}\right)\right]_{\rm open}^{\rm close}=0
$$

because the gate and its first derivative vanish at each endpoint. The static approximation sets $D_q=0$. More generally, for any approximation whose integrated work on this weighted drive is zero, Cauchy–Schwarz gives

$$
\frac{\|R-R_{\rm loc}\|}{\|R\|}
\ge \frac{E}{\|u'\|\,\|R\|}.
$$

Independent evaluation of the retarded convolution gives $\|R\|=36.3869186$ and $\|u'\|=6.6534154$. The resulting floor of 1.063% is below the measured errors of 20.42% for the static force and 10.53% for the two-coefficient force. A model with a dissipative term can exchange net work, so this particular bound applies only to zero-work approximations on the specified drive. Converting the transfer into a statement about quantum particle production or the stress driving a cosmological expansion would require the corresponding fluctuations, state and gravitational dynamics.

## Replay

From the study root, run

~~~sh
OPENBLAS_NUM_THREADS=1 python -IB driven_check/energy.py \
  --native ../HDBLAST_COVARIANCE_FEEDBACK_20261119/native \
  --driven driven/DRIVEN.json --output /tmp/ENERGY.json
~~~
