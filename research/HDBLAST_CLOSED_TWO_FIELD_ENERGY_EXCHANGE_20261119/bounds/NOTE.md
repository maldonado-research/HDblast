# Bounds for any free massless scalar wave on the given spacetime

Changing the phase or wavelength of a daughter wave can change the sign of work done by the parent. Its possible size is nevertheless bounded on a given background over a finite time interval. This calculation uses local energy positivity on a whole closed three-sphere. It requires no cutoff in the daughter wavelengths. The assumption of a prescribed parent field and geometry enters when evaluating the bound on the archived first-passage solutions.

Use the unit-lapse proper time of the archive. Let $K=Z\dot\chi^2/2$, $G_r=Z(D_{\hat r}\chi)^2/2$, and $G_t=Z[(D_{\hat\theta}\chi)^2+(D_{\hat\varphi}\chi)^2]/2$. Hats denote spatial orthonormal components. All three contributions are nonnegative. Write $D=H_r+2H_t$ for the volume expansion rate, $s=\partial_t\ln Z$, and $E=\int(K+G_r+G_t)dV$. Stress conservation for a massless daughter with action $-Z(\Theta)(\partial\chi)^2/2$ gives

$$
\dot E=\int dV\left[(-D-s)K+(D-2H_r+s)G_r+(D-2H_t+s)G_t\right].
$$

The spatial flux integrates out over the closed slice. The terms proportional to $s$ sum to the direct coupling power $P_Z=\int s(G_r+G_t-K)dV$. The other terms come from the evolving metric. The identity also applies to a positive ensemble of such classical solutions. It supplies no positivity bound on an arbitrary vacuum-subtracted quantum energy difference.

Let $l(t)$ and $u(t)$ be the global minimum and maximum over the three rates multiplying $K,G_r,G_t$ and over the spatial slice. Energy positivity gives $lE\leq\dot E\leq uE$ and $|P_Z|\leq\|s\|_\infty E$. Integration yields

$$
\exp\!\left(\int_{t_0}^{t_1}l\,dt\right)
\leq \frac{E(t_1)}{E(t_0)}
\leq \exp\!\left(\int_{t_0}^{t_1}u\,dt\right),\qquad
\frac{|W_Z|}{E(t_0)}\leq\int_{t_0}^{t_1}\!dt\,\|s(t)\|_\infty
\exp\!\left(\int_{t_0}^t u\,dt'\right).
$$

The denominators require positive starting energy. Zero classical initial data give the identically zero solution. The upper bound also holds for the time integral of the absolute coupling power. The expansion rates retain their radial and tangential values even where expansion and contraction coexist.

| Archived first passage, time window $x_*\pm0.005$ | Lower and upper bounds on $E_1/E_0$ | Upper bound on absolute coupling work divided by initial energy |
|---|---:|---:|
| $A=+0.001$, poles | 0.98291 to 1.00971 | 0.007138 |
| $A=-0.001$, equator | 0.94824 to 1.08138 | 0.042212 |

The inequalities are analytic. The coefficients in the table are numerical evaluations on the stipulated evolving backgrounds, rounded to the precision shown. Extremum searches include the poles and equator and refine spatial extrema between sampled angles. Independent runs increase the spatial degree from 48 to 80, the number of candidate angles from 121 to 241, and the time quadrature from 33 to 65 nodes. Comparing them tests the numerical evaluation. The reported digits carry no rigorous interval-arithmetic certificate. The result bounds energy *relative to the daughter’s own initial energy in that window*. It concerns these two short intervals. Each new background and longer time interval requires its own coefficients. [Coefficients and convergence](BOUNDS.json) · [Calculation](envelope.py)
