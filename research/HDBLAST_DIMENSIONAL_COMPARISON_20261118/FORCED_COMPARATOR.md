# A smooth four-dimensional pulse that matches the two-window statistic

A minimally coupled scalar on the four-dimensional background shell can reproduce every limiting HDBLAST two-window ratio in the registered band. Only the mass is varied. The example has zero initial field and velocity, a fixed source profile, a fixed onset and a time-independent gain. These are continuum solutions on the same background clock and sphere. The forcing mechanism is specified independently in four dimensions. Equality of the statistic provides no identification of the two preparations or their dynamics.

## Fixed pulse and local equation

Use the intrinsic metric $h=\rho_b^2[-dq^2+\cosh^2q\,\kappa]$, the past-directed clock $x=q/\ell$, and the transported harmonic $Y=n_1^2-n_2^2$, with $-\Delta_\kappa Y=8Y$. Set $\ell=10^{-6}$ and $H=25/16$. The scalar is $a(x)Y$. Let $c=129/128$, $r=1/128$ and let $p(x)$ be the positive unit-integral multiple of $\exp[-1/(1-((x-c)/r)^2)]$ for $|x-c|<r$, extended by zero. Its support is $[1,65/64]$. In $(\Box_h-m^2)\Psi=-J$ choose

$$
 J(q,n)=(\rho_b\ell)^{-2}\cosh^{-3/2}(\ell x)\,p(x)Y(n),\qquad \mu=m\rho_b\ell\in[6,7].
$$

The coefficient of $Y$ in the source is nonnegative. The source has quadrupolar spatial signs. An overall sign reversal is available by orienting the input oppositely. The field is exactly quiet through $x=1$ by zero-data uniqueness. With $u=\cosh^{3/2}(\ell x)a$, differentiating the metric wave operator gives

$$
 u''+[\mu^2+D(x)]u=p(x),\qquad D(x)=\ell^2[(35/4)\operatorname{sech}^2(\ell x)-9/4].
$$

Here primes denote derivatives in $x$. The nonzero mass scale and pulse duration are expressed in the same clock units as the HDBLAST packet. Their physical realizability is an additional assumption.

## Exact solution comparison and uniform error

Let $u_f$ solve the equation with $D=0$. After the pulse, the sine-kernel solution and the symmetry of $p$ give

$$
 u_f(x)=\frac{C(\mu)}{\mu}\sin(\mu(x-c)),\qquad
 C(\mu)=\int p(s)\cos(\mu(s-c))\,ds\ge1-(7r)^2/2>0.99.
$$

On $[1,H]$, the retarded kernel is bounded by $1/6$, $|D|\le7\ell^2$, and $\|u_f\|_\infty\le1/6$ from the unit source mass. With $k=7\ell^2(H-1)/6$, Duhamel's formula gives $\|u\|_\infty\le1/[6(1-k)]$ and $\|u-u_f\|_\infty\le k\|u\|_\infty$. Since $\log\cosh q\le q^2/2$, the conversion back to $a$ costs at most $(3/4)(\ell H)^2\|u\|_\infty$. Thus $\|a-u_f\|_\infty<10^{-12}$. After multiplication by $\mu/C(\mu)$, the absolute error in either window mean is below $10^{-11}$.

For either window $[l,b]$ after the pulse, its free mean times $\mu/C(\mu)$ equals

$$
 \sin[\mu((l+b)/2-c)]\,\operatorname{sinc}[\mu(b-l)/2],\qquad \operatorname{sinc}z=\sin(z)/z.
$$

Use $W_1=[5/4,11/8]$ and $W_2=[49/32,25/16]$, each normalized by its own length. Their center lags are respectively $39/128,69/128$. Their half-widths are $1/16,1/64$. The first normalized mean exceeds $4/5$ throughout $6\le\mu\le7$. Indeed its sine decreases across this interval of phases inside $(\pi/2,\pi)$, and $\operatorname{sinc}(\mu/16)\ge1-(7/16)^2/6$. Including the curvature error gives the stated lower bound. A degree-23 sine Taylor polynomial with absolute remainder $|z|^{24}/24!$ gives the following outward intervals for the exact curved-space ratio $R=\overline a_{W_2}/\overline a_{W_1}$.

| Scaled mass $\mu$ | Enclosure of $R$ |
| --- | --- |
| 6 | $[-0.098,-0.097]$ |
| 7 | $[-0.720,-0.719]$ |

The linear ODE solution depends continuously on $\mu$ on the compact clock interval, and its first window mean stays away from zero. The intermediate value theorem therefore covers the entire HDBLAST band $(-0.406,-0.330)$. For any selected limiting ratio in this band a single mass and a common gain match both window means. If the finite-amplitude target ratio varies, a mass matching that ratio may also vary. No equality of the full waveform, its energy budget, the generating source or the finite-amplitude field follows from this fit. The centered slope for the comparator evaluates to $a$ on these background cuts because the scalar is a multiple of $Y$.
