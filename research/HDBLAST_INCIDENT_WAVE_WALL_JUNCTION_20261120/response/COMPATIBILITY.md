# Can a gravitating wall follow the prescribed excursion?

In the reflecting-boundary calculation a chosen force law moves the boundary. The following test asks what surface density and pressure would produce such motion if the wall obeyed five-dimensional gravity and the two attached bulk regions kept their original static fields. In this restricted geometry, any nontrivial excursion that begins and ends at rest violates the surface null energy condition at some time. For a bounded excursion extending through all proper times, the energy averaged along a complete light ray on the wall is strictly negative. The original tuned scalar boundary condition also fixes the wall at rest.

## Junction equations

Take the five-dimensional metric $ds^2=dy^2+a(y)^2(-dt^2+d\mathbf{x}^2)$ with $a=e^A$. Join two identical copies of the increasing-$y$ region at the homogeneous wall. Use proper time $\tau$, proper radial speed $s=dy/d\tau$, $\gamma=\sqrt{1+s^2}$, and the normal pointing toward increasing $y$ in the retained copy. The $(t,y)$ components of the tangent and normal are $U^M=(\gamma/a,s)$ and $n^M=(s/a,\gamma)$. With $K_{ab}=e^M_a e^N_b\nabla_M n_N$, their mixed extrinsic curvatures obey

$$
K^i{}_i=A_y\gamma\quad\text{(each spatial direction)},\qquad
K^\tau{}_\tau=n\cdot\nabla_UU=A_y\gamma+\frac{\dot s}{\gamma}.
$$

Let $S^a{}_b=\mathrm{diag}(-\rho,p,p,p)$ include all surface stress, including tension. The mirror Israel equation is $K^a{}_b=-(\kappa_5^2/2)(S^a{}_b-S\delta^a_b/3)$. Its spatial component and the difference of its temporal and spatial components give

$$
\kappa_5^2\rho=-6 A_y\gamma,\qquad
\frac{\kappa_5^2}{2}(\rho+p)=\frac{\dot s}{\gamma}
=\frac{d}{d\tau}\operatorname{arsinh}(s).
$$

The last quantity is the derivative of the wall's rapidity. The null energy condition $\rho+p\geq0$ would make rapidity nondecreasing. Equal initial and final rapidity then imply constant rapidity throughout. This argument applies also to asymptotic rest at both ends. A pure tension contributes zero to $\rho+p$.

For the archived scalar $\phi=-\tanh y$, $W=1-\phi+\phi^3/3$, $\phi_y=W_\phi$ and $A_y=-W/3$. The original scalar jump condition requires the normal derivative to equal the static value $W_\phi(\phi)$ at the wall, in the normalization with tuned tension $2W$. For a moving cut with the same scalar on both sides the derivative is $\gamma\phi_y$. Since $\phi_y=-\mathrm{sech}^2 y$ is nonzero at finite $y$, this independent junction condition requires $\gamma=1$. A brane source which couples directly to the scalar would change the equation. The normal-derivative mismatch is quadratic in a small speed, while the required surface enthalpy is linear in the acceleration.

The bulk scalar can exchange energy with a moving cut. Differentiating the density above and using $H=A_y s$ gives $\kappa_5^2[\dot\rho+3H(\rho+p)]=-6A_{yy}s\gamma=2\phi_y^2s\gamma$ in the archived geometry. The rapidity argument already includes this effect through the junction equations.

## Energy averaged along light rays

Write $r=\operatorname{arsinh}(s)$ and consider an affinely parametrized light ray confined to the induced spatially flat metric $-d\tau^2+a(y(\tau))^2 d\mathbf x^2$. Its tangent $\ell^a$ obeys $d\tau/d\lambda=C/a$, with constant $C>0$ fixed by its conserved spatial momentum. The surface averaged null energy is

$$
\int S_{ab}\ell^a\ell^b\,d\lambda=C\int\frac{\rho+p}{a}\,d\tau
=\frac{2C}{\kappa_5^2}\left[\frac{r}{a}\right]_{-\infty}^{+\infty}
+\frac{2C}{\kappa_5^2}\int\frac{A_y}{a}\,s\operatorname{arsinh}(s)\,d\tau .
$$

For a smooth bounded rest-to-rest excursion the boundary term vanishes. The scale factor stays positive and bounded, so $\lambda$ spans the full real line when $\tau$ does. On this branch $A_y=-W/3<0$ everywhere and $s\operatorname{arsinh}(s)>0$ whenever $s\ne0$. Thus the averaged surface null energy is strictly negative for any such nonconstant trajectory with convergent integrals. This applies to complete null geodesics of the *induced* metric. Five-dimensional null rays obey a different geodesic equation.

## Gaussian check

For $y=y_0+d\exp[-\tau^2/(2T^2)]$, $T>0$, the speed extrema occur at $\tau=\pm T$. Write $h=\rho+p$, $h_+=\max(h,0)$ and $h_-=\max(-h,0)$. Integrating the rapidity changes gives

$$
\int h_+\,d\tau=\int h_-\,d\tau=
\frac{4}{\kappa_5^2}\operatorname{arsinh}
\left(\frac{|d| e^{-1/2}}{T}\right).
$$

The integrated negative magnitude is positive for every $d\ne0$. Its equality to the integrated positive magnitude is an integral over the wall's proper time. The table includes the affine integral with $C=a(y_0)$ so that $d\tau/d\lambda\to1$ at both ends.

| $d$ | $T$ | $\kappa_5^2 h(0)$ | $\kappa_5^2\int h_-d\tau$ | Peak scalar derivative ratio $\gamma$ | $\kappa_5^2\int S_{ab}\ell^a\ell^b d\lambda$ |
|---:|---:|---:|---:|---:|---:|
| 0.1 | 1 | −0.2 | 0.242464 | 1.001838 | −0.00802592 |
| −0.1 | 1 | 0.2 | 0.242464 | 1.001838 | −0.00713564 |
| 0.4 | 0.5 | −3.2 | 1.871829 | 1.111505 | −0.29276706 |
| 0 | 1 | 0 | 0 | 1 | 0 |

All examples use $y_0=0.3$ in wall units. The accompanying calculation integrates $h_+$ and $h_-$ separately, evaluates the curvature from numerical derivatives of the metric and tangent, and checks the surface energy identity. Across these cases both proper-time integrals agree with the analytic result to $10^{-11}$. Direct quadrature of the affine null energy and its integrated-by-parts form agree to $2\times10^{-16}$. The largest sampled absolute residuals in the coordinate junction and energy identity are below $2\times 10^{-12}$ and $7\times 10^{-12}$, respectively.

The calculation isolates the assumption that both mirrors keep the same *frozen* scalar and warp geometry throughout the motion. An incident pulse that gravitates can perturb the bulk fields, so the junction equations then contain additional terms. Establishing a physical response requires solving for those fields and matching the moving wall with its stated action. The mechanical boundary experiment measures a response under its selected force law. It supplies no evidence for a higher-dimensional origin of the Big Bang.
