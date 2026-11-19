# Where the pair energy can come from

The observable is the energy carried by massless Dirac pairs after the coupling has switched off. At order $y^2$ one can compute that energy from the previous pair spectrum and the parent response from the *same* vacuum two-point function. This identifies the terms that must balance on the prescribed expanding wall. Energies below are per unit comoving volume in $y^2 H_0^4$ after integrating over parent momentum. A physical energy density at readout has an additional factor $a_f^{-3}$. [Previous preparation and spectrum](https://doi.org/10.5281/zenodo.24357388)

## The finite Hamiltonian

Use dimensionless shell proper time $s=H_0\tau$, scale factor $a$, comoving parent momentum $q$, parent field $\chi$ and its canonical momentum $\Pi=a^3\dot\chi$. For one real canonical mode at $q$, the generator is

$$
H=\frac{\Pi^2}{2a^3}+\frac{a q^2+a^3 m^2}{2}\chi^2
 +\sum_j\frac{\pi_j^2+\omega_j^2 x_j^2}{2a}
 +yf\chi B+\frac{\kappa y^2}{2}a f^2 C_\Lambda\chi^2,
\qquad B=\sum_jg_jx_j.
$$

Here $\omega_j$ is the **total comoving energy of the pair**, $q<\omega_j<\Lambda$. The two calibrations use $\kappa=0$ and $1$. The bath starts in its conformal vacuum, independent of the parent Gaussian state. Its oscillator strengths obey

$$
g_j^2=2\omega_j K_q(\omega_j)\,\Delta\omega_j,\qquad
K_q(\omega)=\frac{\omega^2-q^2}{8\pi^2},\qquad C_\Lambda=\sum_j g_j^2/\omega_j^2.
$$

The weights $\Delta\omega_j$ stand for quadrature weights on this interval. The $q$ integral has measure $q^2dq/(2\pi^2)$ over $0<q<8$, already included in the native weights. With $d\eta=ds/a$, the positive-frequency vacuum correlation of $B$ is $\sum_j K_q(\omega_j)\Delta\omega_j\,e^{-i\omega_j(\eta-\eta')}$. It reproduces the Dirac pair spectral density. In the first-order transition amplitude $ds\,f\chi=d\eta\,f(a\chi)$, recovering the previous switched amplitude and hence the energy $\int (\omega/a_f)\,dN$ at order $y^2$. At this order the two-point function suffices. Higher orders require the corresponding fermionic correlation functions.

Completing the square in $x_j$ produces a static energy shift $-y^2 a f^2 C_\Lambda\chi^2/2$. Setting $\kappa=1$ cancels this shift at each instant. This is a choice of finite Hamiltonian calibration. The numerical cutoff selects comoving pair energy, and its proper energy $\Lambda/a$ changes with the expansion.

## Work identity in proper time

Let $X_0=\langle\chi^2\rangle_{y=0}$ on the same trajectory. Write $X_2,P_2$ for the $y^2$ coefficients of the parent variance and momentum variance, $T$ for the $y$ coefficient of $\langle\chi B\rangle$, and $D$ for the $y^2$ coefficient of the *proper* bath excitation energy above its vacuum. Quadratic parent energy at this order is

$$
E_{\chi,2}=P_2/(2a^3)+(a q^2+a^3m^2)X_2/2.
$$

At either endpoint $f=0$, so both the interaction energy and the explicit calibration term vanish. The factorized initial condition gives $E_{\chi,2}(s_i)=D(s_i)=0$. Applying $d\langle H\rangle/ds=\langle\partial H/\partial s\rangle$ at fixed canonical variables gives a per-mode work identity,

$$
E_{\chi,2}(s_f)+D(s_f)=W_f+W_m+W_a.
$$

The three integrated work rates, with $h=\dot a/a$, are

$$
\begin{aligned}
\dot W_f &= \dot f\,[T+\kappa a f C_\Lambda X_0],\\
\dot W_m &=\tfrac12 a^3\,\dot{(m^2)} X_2,\\
\dot W_a &=h\left[-\frac{3P_2}{2a^3}+\frac{a q^2+3 a^3 m^2}{2}X_2-D+\frac{\kappa}{2}a f^2 C_\Lambda X_0\right].
\end{aligned}
$$

The term $-hD$ redshifts the massless daughters. At an interior time the energy on the left also contains $fT+\kappa a f^2 C_\Lambda X_0/2$. The subtraction here compares coupled and uncoupled histories of the same prepared state. A separate subtraction compares two initial parent states. The two comparisons answer different questions. Equating the outgoing energy to parent energy lost requires the *integrated sum* of the three work terms to vanish. Closing the displayed identity with a prescribed wall accounts for the work done by that wall and by the switch. Evolving the wall and the field that controls the switch would be needed to identify their own energy sources.

## An independent cutoff test without evolving any daughters

At order $y$, the correlation $T$ is independent of $\kappa$. Switching between the two calibrations therefore shifts the switch work by an amount that can be calculated entirely from the free parent,

$$
W_f^{(1)}-W_f^{(0)}=\sum_q w_q C_\Lambda(q) I_q,
\qquad I_q=\int ds\,a f\dot f X_0
      =-\tfrac12\int ds\, f^2\frac{d(aX_0)}{ds}.
$$

The superscripts on $W_f$ denote $\kappa$. Propagating the three free canonical covariances gives $d(a X_0)/ds=\dot a X_0+2 C_{\chi\Pi,0}/a^2$. Both forms of the integral use the archived retained-state quiet trajectory, the fork covariance and the wide switch. They agree to $3\times10^{-12}$ or better in each momentum-weighted mode on either spatial grid. At fixed $q$ the continuous-frequency static calibration has the exact form

$$
C_\Lambda(q)=\frac{\Lambda^2-q^2-2q^2\ln(\Lambda/q)}{8\pi^2}.
$$

For $\Lambda>8$ the work shift on this wall is therefore $A\Lambda^2+B\ln\Lambda+C$, using the numerical cutoff in native $H_0$ units. The table uses continuous pair frequency and the archived 24-node parent momentum rule.

| Retained quiet wall | $A$ | $B$ | $C$ | Shift at $\Lambda=60$ | Shift at $\Lambda=120$ |
|---|---:|---:|---:|---:|---:|
| Coarse | 0.00238836 | −0.161980 | 0.216298 | +8.15118 | +33.83316 |
| Fine | 0.00239086 | −0.162149 | 0.216524 | +8.15973 | +33.86862 |

The outgoing energy coefficient at a given cutoff is the same in these two calibrations because their parent solutions first differ at $y^2$, two powers beyond the free parent entering the leading pair yield. The positive quadratic work difference shows how strongly an assignment of energy to the switch can depend on which static parent response was fixed. This statement concerns perturbative coefficients. An arbitrarily large cutoff at fixed finite $y$ can take the model outside its perturbative regime. [Numerical coefficients and independent checks](DISPERSION.json)

## Static dispersion and the remaining ultraviolet dependence

For a constant background choose $a=1$. The retarded bath response at complex frequency $z$ with $\mathrm{Im}\,z>0$ is

$$
\Sigma_\Lambda(z)=\int_q^\Lambda\frac{2\omega K_q(\omega)}{\omega^2-z^2}\,d\omega
=\frac{\Lambda^2-q^2+(z^2-q^2)[\ln(\Lambda^2-z^2)-\ln(q^2-z^2)]}{8\pi^2}.
$$

The logarithms are continued from the upper frequency half-plane. The scalar inverse response is $q^2+m^2-z^2-y^2 f^2[\Sigma_\Lambda(z)-\kappa\Sigma_\Lambda(0)]$. For $q<E=\sqrt{q^2+m^2}<\Lambda$, $\mathrm{Im}\,\Sigma(E+i0)=\pi K_q(E)$. At constant $f=1$ the leading occupation decay width is $\Gamma/y^2=\pi K_q(E)/E=m^2/(8\pi E)$, recovering the standard time-dilated massless-Dirac width. The direct complex quadrature and the spectral discontinuity check both pass. [Executable comparison](dispersion.py)

The subtraction $\Sigma_\Lambda(z)-\Sigma_\Lambda(0)$ leaves $z^2\ln(\Lambda^2)/(8\pi^2)$ at fixed $z,q$. Thus static matching alone leaves a frequency-dependent matching condition open. The $q$ dependence of the static calibration also matters for interpreting it as a spatial field theory. A local, covariant stress tensor and matching of the gravitational, scalar and switching sources are further physical inputs for a regulator-independent energy attribution on the wall. The width check constrains the on-shell dissipative normalization, while the finite balance constrains consistency of the chosen regulator and time-dependent Hamiltonian.

To replay the analytic and native calibration checks, pass the directory containing the authenticated parent study's FORK.json and WALLS.zip.

~~~sh
OPENBLAS_NUM_THREADS=1 python -IB dispersion.py /path/to/native
~~~
