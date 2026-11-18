# The normal launch and shell signal give a signed mixed coefficient

Let $w$ denote the small positive inward width of a wedge of launch directions around a fixed tangent-sky direction. Let $d=v_1^2-v_2^2$, $\epsilon$ be the canonical shell amplitude, and $f(\lambda)$ be its first scalar trace on the ray for $d=1$. Use the affine frequency convention and the same segment as for the tangent scalar flux. Inward normal launch has background scalar gradient $p_{\rm in}=-p_b>0$. In the leading uniform average over normal launch component $s\in[0,w]$, the scalar-derivative term is $p_{\rm in}s+d\epsilon f'$. With Green weight $\lambda(L-\lambda)/(3L)$ this square gives the quadratic form

$$
P w^2+H w d\epsilon+F(d\epsilon)^2,\qquad
P=\frac{p_{\rm in}^2L^2}{54},\qquad
H=\frac{p_{\rm in}}{3L}\int_0^L\lambda(L-\lambda) f'(\lambda)\,d\lambda,
$$

where $F$ is the axis scalar-flux coefficient. The factors follow from the first two uniform moments $w/2$ and $w^2/3$, including the cross factor of two when squaring. The constants in the registered normalization obey

$$
1.17988\times10^{-21}<P<1.43502\times10^{-21},\qquad
5.35599\times10^{-23}<H<5.38988\times10^{-23}.
$$

Thus the scalar contribution has opposite mixed terms on interchanging the two quadrupole axes. These are coefficients of the joint small-width, small-amplitude expansion. A finite aperture held fixed while amplitude varies can sample higher orders. A determination of the complete optical signal must also account for other terms in the optical equations. [Fixed-affine rays and geometry](https://github.com/maldonado-research/HDblast/blob/a011dee2b490faa56ae81d3165b4b5fa5bbfe636/research/HDBLAST_FINITE_LIGHT_CONES_20261118/DERIVATION.md)

## Exact cancellation in the mixed term

Use $t=\lambda/L$, $y=(x-a)/\Delta$, $a=49/32$, $b=99/64$, $\Delta=b-a$ and $\ell=10^{-6}$. With $G=3/(\rho_b^3p_b)$ and the normalized past waveform $u$, write $f=G\ell\Delta^2 y^2g(y)u(a+\Delta y)$. The factor $g$ combines two powers of $\sin\theta/(q-q_0)$ with $\operatorname{sech}^{3/2}q$, so $C^3\leq g\leq1$ for $C=1-10^{-11}$. Define $J=L/(\rho_b\ell\Delta)\in[1,K]$, $K=1+10^{-11}$. Integration by parts uses the zero boundary values of $t(1-t)$ and gives $H=(p_{\rm in}L/3)\int_0^1(2t-1)f\,dt$. In particular the exact identity $p_{\rm in}G=-3/\rho_b^3$ retains the common shell calibration. Put

$$
m=\int_0^1(2t-1)y^2 gu\,dt,\qquad
H=-\frac{\ell^2\Delta^3 J}{\rho_b^2}m .
$$

The replay bounds $g(dt/dy)$ in $[C^3/K,K]$ and replaces $2t-1$ by $2y-1$ with absolute error at most $2(K-1)$. This follows from $1/K\leq dt/dy\leq K$ and $t(0)=0$. The polynomial weight $(2y-1)y^2$ has one sign on each of the 64 registered cells because its sign change is on a grid edge. Its exact cell integral is the difference of the endpoint values of $y^4/2-y^3/3$. Multiply this weight by the cell enclosure of $gu(dt/dy)$ in the correct order and add the integrated clock error. The result is $-0.238507<m<-0.237009$. The positive $P$ interval follows directly from the archived $p_b,\rho_b$ boxes and $L=\rho_b\ell\Delta J$. [Gain, units and transfer bounds](https://github.com/maldonado-research/HDblast/blob/a011dee2b490faa56ae81d3165b4b5fa5bbfe636/research/HDBLAST_DIMENSIONAL_COMPARISON_20261118/SINGLE_FIELD_TEST.md) · [Calculator](wedge_coefficients.py) · [Exact intervals](WEDGE_COEFFICIENTS.json)

## How much of the scalar flux can interference remove?

The two coefficients and the tangent flux share their physical scale. Denote the positive lower endpoint of the flux polynomial moment by $m_0$. The tangent calculation gives $F\geq(G\ell\Delta^2)^2m_0^2/[3 K^3(59/420)]$. Combining this inequality with the *exact* formulas for $P,H$ gives

$$
\frac{H^2}{4PF}\leq\frac{9}{2} K^3(59/420)\frac{m^2}{m_0^2}<0.589961 .
$$

Completing the square in $w$ now bounds the quadratic form from below by $[F-H^2/(4P)](d\epsilon)^2$. Its coefficient is greater than $0.410039 F$ and greater than $3.51758\times10^{-25}$. Thus even an unconstrained choice of signed width in this quadratic approximation cannot cancel the entire scalar term when $d\epsilon\ne0$. The cancellation estimate applies to the scalar contribution. A sign-resolved finite-amplitude mixed optical signal, in particular a difference of two positive quantities, requires control of the remaining optical expansion.
