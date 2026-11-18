# Fixed observed frequency ratio on two sky caps

Use elapsed proper-clock labels $q=\ell x$, $\ell=10^{-6}$, receiver $x_0=49/32$, background emitter $X=99/64$ and $d=X-x_0=1/64$. Define $s_i=\sinh q_i$, $c_i=\cosh q_i$. The received source frequency sets the background redshift $z_*=c_0/c_1-1$, between $-2.40479\times10^{-14}$ and $-2.40478\times10^{-14}$. It is a blueshift. The initial condition for the geodesic source population is rest on the common quiet slice $q=0$. Sources held at constant angular coordinates in the conformal representative give a separate clock prescription.

Let $v$ be the receiving direction, $D(v)=v_1^2-v_2^2$ and $\theta=\int_{q_0}^{q_1}\operatorname{sech}q\,dq$. The source lies at $Y=D(v)\sin^2\theta$ on the background ray. With the conformal metric variation $2U(q)Y$, the static-coordinate clock has first fractional shift $-U_1Y$ in its emitted-to-received photon-energy ratio at fixed source coordinate. The geodesic correction to that shift, in elapsed-time velocity conventions, is $+2D(v)\sin\theta\cos\theta\operatorname{sech}q_1\int_0^{q_1}U(q)dq$. The shift in source coordinate at fixed observed ratio is this total first fractional shift divided by $\tanh q_1$. The background angular diameter distance is $\rho_b\cosh q_1\sin\theta$. Its logarithmic derivative in $q_1$ divided by $\tanh q_1$ is $1+\cot\theta/\sinh q_1$. Adding the source metric’s fractional conformal perturbation at fixed coordinate gives the following distance coefficients per unit canonical amplitude and per $D(v)$:
$$
\mathcal C_{\rm supported}=-\frac{\cos\theta}{s_1}U_1\sin\theta,\qquad
\mathcal C_{\rm geodesic}=\frac{\cos\theta}{s_1}
\left[\frac{2}{c_0}\int_0^{q_1}U(q)dq-U_1\sin\theta\right].
$$
These formulas use an observer at the unperturbed quadrupole critical zero. A common observer-time shift that affects both caps isotropically cancels from their half-difference. Directions and endpoint frequency refer to the observer frame, and the result is for angular diameter distance. Each cap’s angular average of $D(v)$ is $\pm18/25$.

## Certifying the sign and scale

Write $U=(\ell B_b/\rho_b)M$ and $J=\int_0^X M dx$. Use $e=\|M-z_0\|_\infty$ bounded by the parent’s `absolute_mean_error`, where $z_0$ is the double integral of the enclosed input starting at the front $x=1$. The proof of this uniform bound is in the parent metric reconstruction. [Metric reconstruction, including the uniform Volterra error](https://github.com/maldonado-research/HDblast/blob/main/research/HDBLAST_JOINT_CURVATURE_20261118/METRIC.md)

An input at $x=t$ contributes $X-t$ to $z_0(X)$ and $(X-t)^2/2$ to its integral. Compute their combination $2\int z_0-d z_0(X)$ with a single kernel $(X-t)(x_0-t)$, changing sign at the cell boundary $x_0$. All weights are integrated exactly on each cell and paired with its signed continuous input enclosure. With $C=1+10^{-11}$ one has $1\le c_i,\sinh q_1/q_1\le C$ and $d/C^2\le\sin\theta/\ell\le d C$. The latter follows from $\sin\theta=(s_1-s_0)/(c_0c_1)$ and integration of $\cosh q$. For $H\ge\|M\|_\infty$ given by the parent’s `z_bound`, the error in replacing $2J/c_0-M(X)\sin\theta/\ell$ by the single-kernel result is at most
$$
[2(X-1)+d]e+
\left[2(X-1)(1-1/C)+d(C^2-1)\right]H .
$$
Finally $\cos\theta\ge1-(\ell d C)^2$ and $\ell/\sinh q_1=1/[X(\sinh q_1/q_1)]$. These are comparison factors close to one. The physical scale $\ell B_b/\rho_b$ lies between $2.60778\times10^{-8}$ and $2.87606\times10^{-8}$ using the previous branch enclosure.

| Quantity | Outward enclosure |
| --- | --- |
| $M(X)$ | $(0.23776519,0.23868027)$ |
| $J$ | $(0.03184748,0.03197999)$ |
| Geodesic velocity term divided by the magnitude of the supported distance term | $(17.0792,17.2164)$ |
| Supported fractional cap-distance half-difference coefficient | $(-4.9925,-4.5093)\times10^{-11}$ |
| Geodesic fractional cap-distance half-difference coefficient | $(7.2804,8.0630)\times10^{-10}$ |

The last two rows are the derivative at zero amplitude of $(D_+-D_-)/(2D_{A,0})$. Changing the sign of the canonical perturbation reverses the small-amplitude half-difference. In the model the freely falling source’s velocity contributes more than seventeen times the magnitude of the supported clock’s term. The inverse redshift derivative is large on this near-turnaround branch. The above results concern derivatives at the fixed nonzero target blueshift and do not supply a numerical amplitude radius, timing precision, frequency resolution, or astronomical identification. [Exact calculation](redshift.py) · [Rational output](REDSHIFT.json)
