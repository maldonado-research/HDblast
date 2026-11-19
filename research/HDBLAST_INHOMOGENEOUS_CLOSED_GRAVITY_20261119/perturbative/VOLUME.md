# Why the initially larger constrained slice expands more slowly

At proper time $x=0.12$, the quadratic change in the fractional proper-volume radius has coefficients $+15.8521$ and $+10.7131$ on the two archived hilltops. Relative to each slice's own initial radius, the corresponding changes are $-0.75250$ and $-0.36733$. These quantities have different signs because the constraint already requires the perturbed slice to start larger. [Numerical coefficients](VOLUME.json)

## The exact growth bound

Use unit lapse and zero shift, and follow the initially orthogonal geodesics until the first caustic or the end of the calculation. Let $J=\sqrt{\det\gamma}$ in comoving coordinates, define $r=(J/J_0)^{1/3}$ along each geodesic, and take the smooth limits at the coordinate poles. Write $\Pi=\partial_x q$, $\theta=3r'/r$, and $\sigma^2=\sigma_{ij}\sigma^{ij}$, with no factor of one half in this convention. The reduced canonical-scalar equations give $R_{nn}=\Pi^2-3u$. The irrotational Raychaudhuri equation then gives

$$
\frac{r''}{r}=u-\frac{\Pi^2+\sigma^2}{3}.
$$

The normalization follows from the homogeneous equation $a''/a=u-\Pi^2/3$. Scalar-gradient terms cancel in $R_{nn}$. Gradients still enter the spatial constraint, the scalar evolution and the shear. With $b=\cosh(hx)$ and $h^2=u(0)$, the equation has the exact comparison form

$$
\left[b^2\left(\frac{r}{b}\right)'\right]'
= -b r\left[h^2-u(q)+\frac{\Pi^2+\sigma^2}{3}\right].
$$

Time symmetry gives $r(0)=1$, $r'(0)=0$. While $u(q)\le h^2$ and $J>0$, $r/b$ is nonincreasing at each point. Therefore $V(x)/[V(0)b^3]$ is nonincreasing and $H_V=V'/(3V)\le h\tanh(hx)$. Equality through an interval requires the nonnegative source to vanish there. This is a finite-amplitude consequence of Einstein gravity for this choice of slicing and matter. Absolute volume at a given time also depends on the initial constraint solution.

The fitted stationary potentials are concave over the two field windows used in the previous study. On every intersecting polynomial segment we set $t=(\Theta-\Theta_i)/(\Theta_{i+1}-\Theta_i)\in[0,1]$. Write $u=N/(v_0 f^2)$, $N=(1-\phi/\phi_m)^2(1+b_\phi\phi)$. In the numerator of $d^2u/dt^2$ with denominator $v_0 f^4$, the polynomial is $N''f^2-4N'ff'-2Nff''+6N(f')^2$. For a polynomial $P$, $P(t)\le P_0+\sum_{j\ge1}\max(P_j,0)$. Applying this bound with exact rational representations of the stored spline coefficients certifies positive $f$ and negative curvature on every segment. The least restrictive curvature bounds, rounded towards zero, are $u_{,qq}\le-43.06$ over $|q|\le0.04$ on the positive-minimum branch and $u_{,qq}\le-34.81$ over $|q|\le0.08$ on the negative branch. The archived residual-slope correction enforces stationarity at $q=0$ and preserves these curvatures. The bounds apply to the reconstructed potential.

## Independent second-order solution

For $q=A T(x)U_2(\cos\chi)+O(A^2)$, use the round average $\langle U_2^2\rangle=1$, the eigenvalue $k^2=8$, and $\mu^2=-3u_{,qq}(0)>0$. At the stationary bounce $a(x)=\cosh(hx)/h$ and $H=a'/a$. The scalar obeys

$$
T''+3HT'+(8/a^2-\mu^2)T=0,\qquad T(0)=1,\quad T'(0)=0.
$$

The regular, reflection-even constrained branch has no order-$A$ metric source or independent metric perturbation. To fix the initial spatial scale, insert $\gamma_{ij}=h^{-2}(1+A^2p)^4\Omega_{ij}+O(|A|^3)$ in $^{(3)}R=|Dq|^2+6u$ at $K_{ij}=\Pi=0$. The angular mean of the order-$A^2$ equation gives $-24h^2\langle p\rangle=8h^2-\mu^2$. For $R_V=(V/2\pi^2)^{1/3}$ and $R_V/a=1+A^2 B+O(|A|^3)$, this fixes $B(0)=2\langle p\rangle=(\mu^2/h^2-8)/12$ and $B'(0)=0$.

At the same order, $\sigma_{ij}\sigma^{ij}$ and the variance of the local expansion are both $O(A^4)$. Spatial averaging the local acceleration, or expanding the exact comparison equation, gives

$$
B''+2HB'=F(x),\qquad F(x)=-\frac{\mu^2}{6}T^2-\frac{1}{3}(T')^2.
$$

The differences $H_V-H$ and $R_V''/R_V-h^2$ have coefficients $B'$ and $F$. The integrated Hamiltonian constraint provides a further identity. On a round closed background the integral of the metric-derivative terms in the first curvature variation vanishes, leaving the mean curvature shift $-12B/a^2$. Thus

$$
2HB'-2B/a^2=\frac{1}{6}\big[(T')^2+(8/a^2-\mu^2)T^2\big].
$$

We solve the two second-order ODEs together and independently differentiate the analytic hypergeometric $T$ to evaluate a Green-function integral,

$$
B(x)=B(0)+\int_0^x\frac{\sinh[h(x-s)]}{h\cosh(hx)\cosh(hs)}
                 \cosh^2(hs)F(s)\,ds .
$$

Across five times through $x=0.12$, the largest absolute difference between the two methods for $T,T',B,B'$ is below $3\times10^{-13}$. The averaged constraint residual is below $2\times10^{-12}$ in reduced units. [Implementation](volume.py) · [Numerical output](VOLUME.json)

| Coefficient per $A^2$ | $\phi_m=+0.5$ | $\phi_m=-0.5$ |
|---|---:|---:|
| Initial $R_V/a-1$ | 16.604574 | 11.080459 |
| $R_V/a-1$ at 0.12 | 15.852073 | 10.713133 |
| $R_V/[R_V(0)\cosh(hx)]-1$ at 0.12 | -0.752500 | -0.367326 |
| $H_V-H$ at 0.12 | -24.098431 | -10.239313 |
| $R_V''/R_V-h^2$ at 0.12 | -717.217143 | -255.578580 |

For $|A|=0.001$, multiply each coefficient by $10^{-6}$. In particular the quadratic result predicts positive proper-volume-radius acceleration on both branches through the comparison time.

## Signed nonlinear comparison

At $x=0.12$ we average each observable for equal-magnitude positive and negative launches before dividing by $A^2$. This retains the even amplitude response. The $|A|=0.001$ radius coefficients from the coupled Einstein–scalar evolution are 15.85198 and 10.71298. The corresponding expansion-rate coefficients are -24.11003 and -10.23844. Across the two branches, all four measured coefficients (radius, expansion rate, acceleration, and fractional growth since launch) agree with the second-order predictions within 0.093% at this amplitude and 0.024% at half amplitude. The signed pairs retain the absolute size increase and the slower growth. The decrease in coefficient error by a factor near four on halving supports an even correction of order $A^4$ to these observables over this range. [Coupled comparison](COMPARISON.json) · [Comparison script](compare.py)
