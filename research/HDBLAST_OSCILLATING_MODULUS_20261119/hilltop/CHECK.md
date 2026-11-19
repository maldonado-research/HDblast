# Independent check of a very small closed-universe launch displacement

The reported displacements are relative to the maximum of one specified four-dimensional potential. The very small numbers arise from exponential growth near an unstable maximum. We checked that the early matching procedure retains the expansion and curvature of a time-symmetric closed launch through six thresholds. [Numerical results](CHECK.json) · [Independent calculation](check.py)

## Exact transfer at the maximum

On the closed de Sitter trajectory with $a(x)=\cosh(h_t x)/h_t$ and $h_t=\sqrt{u(\Theta_t)}$, the linear displacement gain $D(x)$ obeys

$$
D''+3h_t\tanh(h_t x)D'-\mu D=0,\qquad
D(0)=1,\quad D'(0)=0,\quad
\mu=-3u''(\Theta_t)>0.
$$

The substitution $z=-\sinh^2(h_t x)$ gives the solution regular and even at the bounce,

$$
D(x)={}_2F_1(A,B;1/2;z),\qquad
A,B=\frac{3\pm\sqrt{9+4\mu/h_t^2}}4 .
$$

Evaluation of this expression checks the primary differential-equation calculation of $\ln D$ and its derivative. At the six thresholds the largest absolute differences are below $3.2\times 10^{-13}$ and $8.8\times10^{-14}$, respectively.

## Keeping the nonlinear geometry from the bounce

We also integrated the unswitched nonlinear scalar and Raychaudhuri equations directly from the tiny launch, without replacing the early geometry by de Sitter. Define $d=\Theta-\Theta_t$, $L=\ln|d|$, $g=\Theta'/d$ and $N=\ln a$. Before the field reverses direction the equations are

$$
L'=g,\qquad g'=-3hg-g^2-3u'(\Theta_t+d)/d,\qquad
N'=h,\qquad h'=e^{-2N}-(dg)^2/2.
$$

The initial conditions are $g=h=0$ and $N=-\frac12\ln u(\Theta_t+d_0)$. The numerator of the potential derivative is a polynomial in the local spline coordinate. We impose the stationary-point condition on its constant term and evaluate the polynomial divided by $d$ directly at small displacements. This permits relative coordinates much smaller than a floating-point step in the original field. The first neglected machine-precision residual in the recentered derivative is recorded in the JSON. The calculation switches to the ordinary nonlinear force once $|d|>2\times10^{-5}$ in reduced Planck units.

At $|d|=0.03$, the primary matched trajectory and the nonlinear bounce agree in time and $\ln a$ to within $2.22\times10^{-6}$ and in velocity to $3.6\times10^{-10}$ relative across all six cases. Thus a transfer using the closed scale factor and its growing velocity carries the curvature information needed at the measured accuracy. Keeping only a time shift on a flat trajectory would discard that information.

At $t=10^{-3}$ we propagated the direct bounce states into the same subsequent scalar evolution and averaged decay model. At the central positive- and negative-minimum seeds, $\ln(g_{\rm decay}/g_{\rm crit})$ relative to zero was $+6.64\times10^{-6}$ and $+1.08\times10^{-6}$, respectively. Both narrow brackets retained opposite signs, with endpoint magnitudes greater than 0.00538. This independently checks the launch transfer. The late oscillator averaging and single decay width are shared assumptions of both tests. [Threshold and averaging checks](HILLTOP.md)

## What precision means for the physical interpretation

For $t=10^{-3}$ the central displacement magnitudes in this potential are approximately $6.41\times10^{-25}$ and $1.25\times10^{-19}$ in reduced Planck units. An unknown change in the dimensionless potential slope would move the maximum by $\Delta\Theta_t\simeq-\Delta u'/u''$. In an absolute field coordinate, a slope error comparable to $4.43\times10^{-23}$ or $5.89\times10^{-18}$, respectively, would already shift the maximum by a full central displacement. These are required slope precisions for locating those launches absolutely. This study provides no five-dimensional or quantum calculation fixing a physical maximum to either precision. A shift of the maximum can leave a well-defined classical sensitivity relative to the shifted maximum. Changes to its curvature, the decay width, or the effective potential away from the maximum would also change that sensitivity.

The comparison with $H_{\rm top}/(2\pi M_{\rm Pl})$ uses the archive's reference displacement. A distribution over initial states or a probability for obtaining the tuned launch would require physical input beyond this evolution.

Run from the study root after reproducing the primary thresholds.

~~~sh
OPENBLAS_NUM_THREADS=1 python -IB hilltop/check.py --archive archive --input /tmp/HILLTOP.json --output /tmp/HILLTOP_CHECK.json
~~~
