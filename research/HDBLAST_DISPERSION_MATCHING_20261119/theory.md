# Calibrating a restoring force and its frequency dependence

A static force measurement constrains one coefficient of a dynamical response. A second measurement can fix its slope. At frequencies away from both measurements, the positive spectral weight of the pair field puts a sign and a quantitative bound on what is still possible.

## Frozen finite-band inverse response

Use the frozen conformal normalization $a=1$, fix $q>0$, set $z=i\xi$ with $x=\xi^2\ge0$, and write $s=q^2$, $L=\Lambda^2$, $\lambda=y^2\ge0$. This is the response of a stationary oscillator. The full time-dependent shell contains additional sources. The parent mode has a restoring term $A$ and kinetic coefficient $B$. The massless Dirac pair weight $K_q(\omega)=(\omega^2-q^2)/(8\pi^2)$ for $\omega>q$ fixes the vacuum two-point response at order $y^2$. Eliminating a representative bath oscillator gives a minus sign for the correction to the parent *inverse* response. With $t=\omega^2$,

$$
D_L(x)=A_L+B_Lx-\lambda W_L(x),\qquad
W_L(x)=\int_s^L \frac{d\sigma(t)}{t+x},\qquad
d\sigma(t)=\frac{t-s}{8\pi^2}\,dt .
$$

The bath quadrature used in the shell satisfies $g_j^2=2\omega_jK_q(\omega_j)\Delta\omega_j$, giving the discrete analog $W_L(x)=\sum_jg_j^2/(\omega_j^2+x)$. Here $x$ is the square of Euclidean frequency. A slope with respect to real $z^2$ therefore carries the opposite sign. A susceptibility defined as the displacement per applied force is $1/D_L$ in this normalization. [Spectral density and previous continuum check](https://doi.org/10.5281/zenodo.24358059) · [Finite shell dynamics](https://doi.org/10.5281/zenodo.24360096)

For the Dirac density,

$$
W_L(x)=\frac{L-s-(s+x)\log[(L+x)/(s+x)]}{8\pi^2}.
$$

Matching $D_L(0)=k_R$ fixes $A_L=k_R+\lambda W_L(0)$. If the bare kinetic coefficient is held fixed while $L$ grows, the remaining frequency slope grows logarithmically. Matching also $D_L'(0)=Z_R$ fixes

$$
B_L=Z_R-\lambda\int_s^L\frac{d\sigma(t)}{t^2}
=Z_R-\frac{\lambda}{8\pi^2}[\log(L/s)-1+s/L].
$$

After both conditions the held-out inverse response is $D_L(x)=k_R+Z_R x-\lambda x^2\int_s^L d\sigma/[t^2(t+x)]$. It lies below its calibrated tangent for $x>0$ and a nonempty positive bath. These conditions are illustrative physical matching conditions on one mode at a time. Treating them as calibrated in nature would require actual measurements. A single local four-dimensional counterterm can also restrict how the fitted coefficients vary across $q$.

## Interpolation and extrapolation

More generally, choose two nonnegative reference points $a,b$, and let $P(x)$ be the affine line through $D(a),D(b)$. When $a=b$, use $D(a),D'(a)$ and its tangent. On any positive spectral band or quadrature,

$$
C(x)\equiv\frac{P(x)-D(x)}{(x-a)(x-b)}
=\lambda\int\frac{d\sigma(t)}{(t+a)(t+b)(t+x)}.
$$

At reference points the right side supplies the continuous limit. The identity follows by subtracting the line through the two values of $1/(t+x)$ inside the spectral integral. The constants $A,B$ cancel identically. In particular $D$ lies **above** the secant for $a<x<b$, and **below** its continuation outside. The normalized residual $C$ is positive and decreases with $x$. More generally, $(-1)^n C^{(n)}(x)\ge0$ for every $n\ge0$, by differentiating the integral. These sign tests concern Euclidean frequency, where every denominator is positive. Inside the physical real-frequency cut the principal-value integral does not have this ordering.

### Bounding a held-out frequency without choosing a high-energy tail

Suppose the spectrum, including its coupling normalization, is known up to $L$ and that its remaining weight is nonnegative on $[L,\infty)$. Assume the twice-subtracted integral converges and that the inverse response has only the affine contact terms $A+Bx$. Subtract the known contribution from $C$ to obtain its missing tail $U(x)$. The two reference conditions fix the sign of this missing curvature. They alone leave its normalization undetermined for an arbitrary positive tail.

Suppose one additional response value at $c$ gives $U(c)>0$. Then at a different, held-out $x$,

$$
\min\!\left(1,\frac{L+c}{L+x}\right)
\le \frac{U(x)}{U(c)}
\le \max\!\left(1,\frac{L+c}{L+x}\right).
$$

Indeed, the ratio is the average of $(t+c)/(t+x)$ with probability weights proportional to $d\sigma/[(t+a)(t+b)(t+c)]$ above $L$. A spectral atom at $L$ attains one endpoint. Moving an atom arbitrarily far into the tail approaches the other endpoint. If the tail has known finite upper support $H$, replace $1$ with $(H+c)/(H+x)$; atomic spectra attain both endpoints exactly. If $U(c)=0$, positivity forces the coupling-weighted tail measure to vanish and hence $U(x)=0$. A measured $U(c)<0$ fails these assumptions unless its uncertainty includes zero.

The interval width in units of $U(c)$ is $|x-c|/(L+x)$ for unbounded upper support. Multiplying by $(x-a)(x-b)U(c)$ converts it into a held-out *inverse-response* interval, reversing the endpoint order when necessary. For example, choose $q=4$, $a=0$, $b=4$, $c=9$, $x=25$ and $\Lambda=80$. The tail ratio must lie between $0.997510$ and $1$. Integrating the actual Dirac tail gives $0.998754$. This is a manufactured verification with a known spectral law. With observations, uncertainty in all the calibration values and the known band would need to be propagated, especially when their difference $U(c)$ is small. [Independent quadrature and endpoint checks](theory_check.json)

### Known massless Dirac ultraviolet continuation

When the stated Dirac spectral density holds all the way to infinity we can bound the tail absolutely. Using $t\ge L$ in each positive factor gives

$$
\frac{\lambda}{8\pi^2 L}
\frac{1-s/L}{(1+a/L)(1+b/L)(1+x/L)}
\le U(x)\le \frac{\lambda}{8\pi^2 L}.
$$

To see this write the tail integrand as $dt/t^2$ times $(1-s/t)/[(1+a/t)(1+b/t)(1+x/t)]$. Its integral gives
$U(x)=\lambda[1/L-(s+a+b+x)/(2L^2)+O(L^{-3})]/(8\pi^2)$. For a common calibrated value and slope at zero, $D_L(x)-D_\infty(x)=x^2 U(x)>0$ and the excess approaches zero as $\Lambda^{-2}$. The rate assumes the specified ultraviolet spectrum. A positivity-only extension with unknown normalization gives the weaker calibrated tail interval above.

## What a finite positive Hamiltonian permits

For a canonical finite bath with positive oscillator masses, a positive parent kinetic term requires $B_L>0$. The potential energy is positive definite precisely when $D_L(0)>0$, by completing the square in the bath coordinates. If $k_R,Z_R>0$ and $\lambda>0$ are both held fixed, then $B_L$ passes through zero at a unique finite $L_*>s$ determined by

$$
\log(L_*/s)-1+s/L_*=8\pi^2 Z_R/\lambda.
$$

The left side grows from zero to infinity and has positive derivative for $L_*>s$. Beyond this value the fixed local kinetic normalization would require a negative bare inertia. Mathematically, $D_L(0)=k_R>0$ and $D_L(x)\sim B_Lx<0$ at large $x$, so the finite inverse develops a zero at positive Euclidean $x$, corresponding to a growing or decaying exponential in real time. As a reproducibility example, $q=4$, $Z_R=0.03$, $k_R=1$, $\lambda=1$ give $\Lambda_*=21.17464$. At $\Lambda=\sqrt2\,\Lambda_*$ the inverse crosses zero at $x\simeq947.97$. These parameters illustrate an exact property of the Gaussian surrogate. Perturbation theory for a physical Yukawa field loses control when the logarithmic correction becomes comparable with the kinetic coefficient. Extrapolating this ceiling to a cosmological scale or to an interacting ultraviolet completion would be unjustified. At $q=0$ the static slope also has an infrared divergence for this continuum density. A nonzero spacelike reference or an infrared prescription is then needed.

## Contact terms and the scope of a calibration

The positivity test cancels two affine local terms. Adding a higher-derivative term $\gamma(x-a)(x-b)$ to $D$ preserves both reference conditions and shifts $C(x)$ by $-\gamma$. A cubic term $\gamma(x-a)(x-b)(x-c)$ can preserve the extra reading at $c$ while shifting the held-out value arbitrarily. The same examples work for a repeated reference. Positivity of a dissipative spectrum by itself does not constrain arbitrary higher-derivative contact terms. A stronger physical claim would need a specified operator content, high-energy behavior or more independent measurements.

Static and kinetic measurements with a constant coupling also leave curvature terms and operators involving derivatives of the switched coupling unconstrained. Two frozen propagators can agree while the stress driving an expanding wall differs. The held-out bounds apply to the stated linear matter response with its allowed contacts. Deducing a physical cosmological stress requires matching the scalar, metric and switching sources in a common action.

Numerical checks compare the raw logarithmic expression with the subtracted quadrature, test both sides of the calibration interval, test sharp atomic bounds in exact rational arithmetic, and locate the finite kinetic ceiling. Run `python -I -B theory_check.py` to reproduce the reported numbers. [Executable checks](theory_check.py)
