# How much of the missing band's force is captured by two local terms?

On the fine quiet shell, a static approximation to the force from pair energies between 60 and 120 has a **20.42%** relative RMS error during the open gate. Including the kinetic coefficient reduces this to **10.53%**. The local terms capture about half of the departure from the static response in the chosen norm. The exact omitted-band force persists after the external drive ends. These results concern the prescribed waveforms and the driven retarded force of a finite bath. [Computed responses and numerical controls](DRIVEN.json)

## Drive and norm

Use conformal time $d\eta=ds/a$ and the preceding study's fine quiet trajectory. Each of its 24 parent momenta $q$ has retained fork covariance $\left(\begin{smallmatrix}X&C\\C&P\end{smallmatrix}\right)$. We initialize two real free-parent solutions from the columns of its lower Cholesky factor. Their outer products sum to this covariance. The drive in channel $l=1,2$ is

$$
u_{ql}(\eta)=a(s)F(\phi_b(s))\chi_{ql}(s).
$$

The gate is the same scalar-dependent coupling as in the previous native calculation, calibrated on fine quiet and active from $s=1.35$ to $2.05$. We propagate the free modes from the original fork through $s=2.10$. Summing over both columns makes the result invariant under an orthogonal rotation of the two real solutions. [Archived trajectories and coupling](https://doi.org/10.5281/zenodo.24360096)

For the eliminated comoving band $60<\omega<120$, integrate the retarded equation for oscillators initially at rest,

$$
Q_{\omega ql}''+\omega^2Q_{\omega ql}=u_{ql},\qquad
R_{ql}=\int_{60}^{120}g_q^2(\omega)Q_{\omega ql}\,d\omega,\qquad
g_q^2(\omega)=\frac{\omega(\omega^2-q^2)}{4\pi^2}.
$$

Primes denote conformal-time derivatives. Here $R$ is the positive retarded response for positive forcing in this equation. The local approximations obtained from the frozen susceptibility are $R_0=C_qu$ and $R_2=C_qu-D_qu''$, with $C_q=\int g_q^2/\omega^2$ and $D_q=\int g_q^2/\omega^4$. Both coefficients are fixed by the spectral measure before loading any force data. We calculate $u''$ by differentiating the free-parent equation, scale factor and scalar gate, including $F_{\phi\phi}$. The native momentum weights $w_q$ already contain the phase-space measure. The in-gate squared norm is

$$
\|e\|^2=\sum_{i=1}^{24} w_i\sum_{l=1}^2\int_{\eta(1.35)}^{\eta(2.05)} |e_{q_i l}(\eta)|^2\,d\eta.
$$

## Force error, late reaction and energy

| Quantity | Fine quiet result |
|---|---:|
| Static relative error $\|R-R_0\|/\|R\|$ | 0.2042065 |
| Static and kinetic relative error $\|R-R_2\|/\|R\|$ | 0.1052940 |
| Remaining fraction of the nonstatic reaction $\|R-R_2\|/\|R-R_0\|$ | 0.5156254 |
| Exact in-gate force norm | 36.38692 |
| Exact force norm over $2.05<s<2.10$ | 0.770056 |
| Postgate norm divided by in-gate norm | 0.0211630 |

All 24 individual momenta show improvement. Their static relative errors range from 19.57% to 21.15% and the corresponding two-term errors range from 10.00% to 10.97%. After closure $u$ and $u''$ vanish, so both local approximations give zero. The 2.12% postgate comparison is a ratio of **integrated norms over windows of different duration**. At the particular endpoint $s=2.10$, the instantaneous weighted norm of the exact force is 4.00684.

The bath's final classical conformal excitation energy, summed over both columns and native weights, is 2.573972033. Direct integration of the applied power $\sum w_qg_q^2 u_{ql}Q_{\omega ql}'$ gives 2.573972034. In closed-cycle work, each real local approximation gives exactly zero, since $-\int u'R_2\,d\eta$ consists of endpoint terms in $u^2$ and $(u')^2$, both zero at the endpoints. This drive tests deterministic polarization and its memory. Computing quantum particle production or changing the five-dimensional geometry would require the corresponding fluctuation, stress and scalar-current terms. [Independent Fourier energy and force checks](../driven_check/ENERGY.json)

## Numerical replay

The solver propagates the free real quadratures with DOP853 and evolves the forced oscillators using fourth-order Runge–Kutta. Each proper-time cell is split into two RK4 steps. Midpoint quadrature in $ds/a$ evaluates the force norms, and Simpson quadrature evaluates the work. Cell boundaries include every native background-interpolation knot and the gate endpoints. Midpoints select interior values for the piecewise second derivative of the archived interpolant. The primary calculation has 2,277 cells and 96 Gauss frequencies. At 1,330 cells the two-term error ratio changes by $2.95\times10^{-8}$ and the postgate norm ratio by $2.72\times10^{-7}$. Reducing the Gauss rule to 48 frequencies at fixed primary time resolution changes the two force-error ratios by less than $6\times10^{-16}$. These numerical tests resolve the much larger discrepancy between the driven reaction and its two local approximations. [Full time and frequency checks](DRIVEN.json)

~~~sh
OPENBLAS_NUM_THREADS=1 python -IB driven.py --repo /path/to/repo --output /tmp/DRIVEN.json
~~~
