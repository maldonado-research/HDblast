# A lower bound on the scalar flux along the tangent light rays

Take the same receiving observer, shell-tangent sky direction $v=(v_1,v_2,v_4)$, affine normalization, and background endpoint clocks $a=49/32$, $b=99/64$ in $x=q/\ell$ as in the finite-light-cone calculation. Write $d=v_1^2-v_2^2$, $\Delta=b-a=1/64$, and $L$ for the fixed physical affine length. For the registered linear scalar field restricted to that ray, define

$$
\mathcal F(v)=\frac1{3L}\int_0^L \lambda(L-\lambda)
               (\partial_\lambda\Phi_1)^2\,d\lambda.
$$

**The scalar contribution obeys $\mathcal F(v)\geq c_\phi d^2$, with $c_\phi>8.57864\times10^{-25}$ per squared canonical amplitude.** Either of the two tenth-sky caps has mean greater than $4.70475\times10^{-25}$. The full tangent-sky mean exceeds $2.28763\times10^{-25}$. The one-third refers to the three-dimensional optical screen in five spacetime dimensions. The result is stated for the coefficient multiplying amplitude squared. A second derivative in amplitude carries a factor of two. The inputs are a computed shell waveform and its certified error enclosures. [Null rays and screen conventions](https://github.com/maldonado-research/HDblast/blob/a011dee2b490faa56ae81d3165b4b5fa5bbfe636/research/HDBLAST_FINITE_LIGHT_CONES_20261118/DERIVATION.md) · [Physical scalar identification](https://github.com/maldonado-research/HDblast/blob/a011dee2b490faa56ae81d3165b4b5fa5bbfe636/research/HDBLAST_SHELL_RESPONSE_20261118/BRIDGE.md)

## Geometry and units

The unperturbed pullback of the scalar to the shell is constant. Along a tangent ray on the umbilical shell the linear pullback is therefore $\Phi_1=S(-\ell x)d\sin^2\theta$, with no first-order correction from displacing that ray tangentially. The bulk contribution $\operatorname{Ric}_5(k,k)=(k\phi)^2$ along the actual tangent ray is consequently quadratic to leading order. This calculation covers shell-tangent launches. [Bulk field equations and shell condition](https://github.com/maldonado-research/HDblast/blob/a011dee2b490faa56ae81d3165b4b5fa5bbfe636/research/HDBLAST_JOINT_CURVATURE_20261118/GRADIENT_POWERS.md)

Set $t=\lambda/L$ and $y=(x-a)/\Delta$. In the background $\lambda/\rho_b=(\sinh q-\sinh q_0)/\cosh q_0$ and $\sin\theta=\lambda/(\rho_b\cosh q)$. Hence, for $q_0=\ell a\leq q\leq\ell b$,

$$
(q-q_0)\operatorname{sech}q\leq\sin\theta\leq q-q_0 .
$$

The upper bound also follows by integrating $d\theta/dq=\operatorname{sech}q$. Since $dt/dy$ is $\cosh(\ell x)$ divided by its mean in $y$, the constant $K=1+10^{-11}$ gives $1/K\leq dt/dy\leq K$, $t\geq y/K$ and $1-t\geq(1-y)/K$. The same bound also follows from $\cosh q\leq1/(1-q^2)$ on this short interval. Put $C=1-10^{-11}$; both $\operatorname{sech}q$ and $\operatorname{sech}^{3/2}q$ exceed $C$ since $(1-(\ell b)^2)^2>C$.

Let $u(x)=\ell\cosh^{3/2}(\ell x)S(-\ell x)/G$ be the normalized past scalar coefficient. The published physical gain has $|G|\geq2493872451/250000000000$. The background geometry and transfer equations give

$$
\Phi_1=d\,G\ell\Delta^2 f(y),\qquad f(y)=y^2g(y)u(a+\Delta y),\qquad C^3\leq g(y)\leq1 .
$$

Here $g=\operatorname{sech}^{3/2}q\,[\sin\theta/(q-q_0)]^2$ with its continuous limit at the vertex. All factors of the shell radius cancel in $t$ and in this expression. In particular

$$
\mathcal F(v)=\tfrac13\int_0^1t(1-t)(\partial_t\Phi_1)^2\,dt
 \geq \frac{(d\,G\ell\Delta^2)^2}{3K^3}\int_0^1 y(1-y)f'(y)^2\,dy .
$$

The interval input for $u$ is the archived future master interval widened by $3\times10^{-6}$ on either side for past transfer. The replay verifies the common gain, wave-table hash, normalization, clock and transfer bound against the archived reports. [Registered waveform normalization](https://github.com/maldonado-research/HDblast/blob/a011dee2b490faa56ae81d3165b4b5fa5bbfe636/research/HDBLAST_DIMENSIONAL_COMPARISON_20261118/SINGLE_FIELD_TEST.md)

## Weighted integral without derivative samples

Choose $w(y)=y(1-y)$ and $p(y)=y+y^2$. The product $wp=y^2-y^4$ vanishes at both endpoints. Integration by parts and Cauchy–Schwarz give

$$
\int_0^1 w f'^2\,dy\geq
 \frac{\left[\int_0^1(2y^3-4y^5)g(y)u(a+\Delta y)\,dy\right]^2}{\int_0^1wp^2\,dy},
\qquad \int_0^1wp^2\,dy=59/420.
$$

On each of 64 continuous time cells, enclose $gu$ by interval multiplication, then integrate the polynomial exactly. On the sole cell containing $1/\sqrt2$, integrate its positive and negative monomials separately to maintain an enclosure. The bracketed moment lies strictly between $0.2468863$ and $0.2490624$. Its positive lower endpoint gives the bound quoted above. This uses a fixed polynomial test and finite-width data alone. The code prints exact rational endpoints. [Calculator](flux.py) · [Results](FLUX.json)

An independent point-difference bound uses the intersections of the neighboring closed intervals at the 63 interior grid points. Write $F=(x-a)^2gu$ and $w_x=(x-a)(b-x)$. The dimensionless integral is bounded below by $(dG\ell)^2/(3\Delta K^3)$ times $\int_a^b w_x(F')^2 dx$. On each of the 62 intervals between successive interior points, the contribution is at least the squared gap between their value enclosures divided by $\int 1/w_x\,dx$. Use zero when the enclosures overlap. Convexity of $1/w_x=[1/(x-a)+1/(b-x)]/\Delta$ allows its endpoint trapezoid as a rational upper bound on this denominator. Summing the disjoint intervals gives a second coefficient greater than $8.4783\times10^{-25}$. Neither calculation assumes a bound on a differentiated waveform or supplies an upper bound on the complete optical deficit.

## Angular averages

The observer’s *tangent* unit sky is $S^2$. Conditional on $v_1=z$, the mean of $d^2$ around the circle is $(19/8)z^4-(7/4)z^2+3/8$. On $v_1\geq4/5$, $z$ is uniform on $[4/5,1]$. Thus the mean of $d^2$ is $10283/18750$ on that cap and on its partner $v_2\geq4/5$. On the full sky, $E[v_1^4]=1/5$ and $E[v_1^2v_2^2]=1/15$, giving $E[d^2]=4/15$. Every individual direction in either cap satisfies $|d|\geq7/25$. The nodal directions $d=0$ have zero scalar flux at this order; additional optical contributions are separate. The averages are over shell-tangent launch directions and the underlying optical screen has three dimensions.
