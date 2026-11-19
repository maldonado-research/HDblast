# Finite retarded response of the massless Dirac pair channel

The outgoing spectrum fixes the part of the response that dissipates energy. Two local normalization conditions also fix a finite response during the switch. This calculation tests both statements on a single mode of the archived expanding wall.

## Dispersion relation

Let $q$ be the comoving wave number, $z$ the conformal frequency in the upper half-plane and $u=z^2-q^2$. The forward Fourier convention is $e^{+i z \eta}$. The Dirac pair weight is $K_q(\omega)=(\omega^2-q^2)\Theta(\omega-q)/(8\pi^2)$. If $F=f a\chi$ drives the bath with interaction Hamiltonian $+y F B$ in conformal time, its first-order induced amplitude is $-y\,\Sigma_{\rm ret}*F$. The sign follows by solving one forced oscillator. With a maximum pair energy $\Lambda$ and squared invariant mass $v=\omega^2-q^2$, the uncalibrated susceptibility is

$$
\Sigma_{\Lambda}(u)=\int_q^\Lambda
\frac{2\omega K_q(\omega)}{\omega^2-z^2}\,d\omega
=\frac{1}{8\pi^2}\int_0^{\Lambda^2-q^2}\frac{v}{v-u}\,dv .
$$

Choose a spacelike reference $u=-\mu^2$, $\mu>0$, and require that both the renormalized value and its derivative with respect to $u$ vanish there. At finite cutoff, take the value and derivative of the same regulated integral. This gives an absolutely convergent dispersion representation as the cutoff grows,

$$
\Sigma_R(u)=\frac{(u+\mu^2)^2}{8\pi^2}\int_0^\infty
\frac{v\,dv}{(v+\mu^2)^2(v-u)}
=\frac{u\log[\mu^2/(-u)]+u+\mu^2}{8\pi^2}.
$$

The logarithm is continued from the retarded half-plane. For real $\omega>q$, $\operatorname{Im}\Sigma_R=\pi K_q(\omega)$. The negative-frequency value is the complex conjugate of the positive-frequency value. The light-cone limit is $\mu^2/(8\pi^2)$ and the slope diverges logarithmically there. The infinite-cutoff spectral function depends on the Lorentz invariant $u$ with the retarded frequency prescription. The external trajectory, expanding wall and switch still specify a particular conformal frame. The main energy ledger's static normalization at $z=0$ is a different prescription, evaluated at $u=-q^2$ for each spatial mode.

Changing $\mu_1$ to $\mu_2$ adds the real polynomial $[u\log(\mu_2^2/\mu_1^2)+(\mu_2^2-\mu_1^2)]/(8\pi^2)$. Its action on a history is local, with $uF=-F''-q^2 F$. At separated times the retarded kernel follows by inserting $e^{-\epsilon\omega}$ in the sine transform and taking $\epsilon\downarrow0$. For a positive lag $t$,

$$
\Sigma_{\rm ret}(t)=2\int_q^\infty K_q(\omega)\sin(\omega t)\,d\omega
=-\frac{\cos(qt)+qt\sin(qt)}{2\pi^2t^3}.
$$

An ordinary time integral of this last expression is valid when the source and readout are separated by a positive gap. At contact its distributional extension depends on the local normalization conditions.

## One native mode and an independent energy check

We used the fine retained quiet wall, the upper of its two central Gauss momentum nodes ($q=4.256227571450422$), and the wide switch $f=\sin^4[\pi(s-1.35)/0.70]$ on $1.35<s<2.05$, zero elsewhere. The free complex mode has the instantaneous positive-frequency scalar normalization at the fork. We propagated its canonical variables in proper time $s$ on the archived wall. Its switched conformal history is $F=f a\chi$. A Fourier–Laplace inversion evaluates the kernel on this history using $\operatorname{Im} z>0$. The code also integrates the separated-time kernel directly over the source. The response below is a mode-function amplitude, per unit drive, with the sign convention $R=\Sigma_R*F$. Complex mode phases refer to the fork.

| Proper time | $R$ at $\mu=2$ | $R$ at $\mu=4$ | Direct separated-time integral |
|---|---:|---:|---:|
| 1.70 | $-17.42452-5.23949i$ | $-14.46826-3.73540i$ | Inside switch |
| 2.07 | $0.00791814-0.00085391i$ | $0.00791814-0.00085391i$ | $0.00791814-0.00085391i$ |
| 2.10 | $0.00562515-0.00092097i$ | $0.00562515-0.00092097i$ | $0.00562515-0.00092097i$ |

The finite response during the switch shifts by the predicted local polynomial. At the two post-switch readouts the direct integral and the frequency calculation agree within $1.4\times10^{-9}$ in amplitude. The post-switch values describe the remaining induced bath polarization per prescribed input. The scalar's backreaction force carries another factor of $f$ and vanishes after the switch.

For a real prescribed source $F_r=\operatorname{Re}F$ the external conformal work delivered to the bath over the full cycle, divided by $y^2$, is

$$
W=-\int F_r'(\eta)\, (\Sigma_R*F_r)(\eta)\,d\eta
=\int_q^\infty \omega K_q(\omega)
\,\bigg|\int F_r(\eta)e^{i\omega\eta}\,d\eta\bigg|^2 d\omega.
$$

This follows from Fourier transformation using the positive retarded cut and the conjugate negative-frequency cut. The local subtraction terms integrate to endpoint terms proportional to $F_r^2$ and $(F_r')^2$, both zero for this switch. We evaluated the two sides independently by time integration of the retarded response and by positive-frequency Gauss integration of the source spectrum.

| Maximum comoving pair energy | Work from retarded response at $\mu=2$ | Energy from pair spectrum |
|---:|---:|---:|
| 40 | 4.34106124 | 4.34106124 |
| 80 | 8.06541094 | 8.06541094 |
| 160 | 8.07160804 | 8.07160804 |
| 800 | 8.07160807 | 8.07160807 |
| Infinite | 8.07160807 | 8.07160807 to cutoff 800 |

Changing $\mu$ from 2 to 4 changes the infinite-cutoff integrated work by less than $5\times10^{-12}$ in this calculation. A hard cutoff of 80 already gives 99.92% of the converged work, yet the complex response at proper time 2.10 still differs by 0.426 from its infinite-cutoff value. Agreement in an integrated radiation observable alone can give a misleading indication of convergence for a separate, time-resolved response.

The work here is per $y^2$ for one real prescribed conformal mode amplitude, in the simulation's units. The quantum ensemble in the main ledger has its own state and phase-space normalization. This two-point response specifies one causal component needed for backreaction. Renormalizing the full stress tensor on the evolving wall would also require local covariant conditions for the geometry, the source and the time-dependent coupling.

## Replay

From the capsule run `python memory/retarded.py` or pass `--native /path/to/native` to the directory holding the preceding capsule's `FORK.json` and `WALLS.zip`. The script reads the two files, checks dispersion quadrature, retarded symmetry, the light-cone limit, pre-switch causality, time refinement, the switch-scale shift, separated-time agreement and the positive work sum rule. It writes the numerical results to standard output. [Computed results](MEMORY.json)
