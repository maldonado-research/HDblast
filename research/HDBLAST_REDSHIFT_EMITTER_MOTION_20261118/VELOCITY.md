# Which velocity errors can mimic a cap signal?

The distance at a specified redshift depends on which emitters provide that redshift. An additional smooth line-of-sight velocity can cancel the complete predicted distance quadrupole at first order. Sampling antipodal sky patches removes a common velocity dipole and leaves a test of spatially varying source motion. Here velocities refer to emitter velocities relative to the specified fiducial ensemble. They are derivatives in the canonical perturbation amplitude, with the speed of light equal to one.

## Redshift selection

At zero amplitude let $Z_0=1+z_0=\cosh q_0/\cosh q_1$ and $D_{A,0}(L)=L$, for receiver-normalized past affine length $L$. The logarithmic redshift slope is nonzero:
$$
G=L\left.\frac{d\log Z_0}{d\lambda}\right|_L
=-\frac{L}{\rho_b}Z_0\tanh q_1
=-\frac{(\sinh q_1-\sinh q_0)\sinh q_1}{\cosh^2q_1}.
$$
Here $q_0=(49/32)\ell$, $q_1=(99/64)\ell$, and $\ell=10^{-6}$. The magnitude of $G$ is about $2.417\times10^{-14}$. An additive error in log frequency is strongly amplified by selecting endpoints at fixed redshift in this model.

Write the emitting rest frame as $u_s$ and the past null tangent at emission as $k=\omega(-u_s+e_s)$ with $e_s$ a unit spatial vector. A first-order boost by $w\perp u_s$ changes $\log Z$ by $b=e_s\cdot w$ at the fixed event, since $g(u_s,k)=\omega>0$. With the vertex unchanged, angular-diameter distance at that event is unchanged by this boost. The source screen is one realization of the null screen quotient $k^\perp/\langle k\rangle$, and adding a multiple of $k$ to a connecting vector preserves its transverse inner products. Unknown intrinsic source size or orientation when inferring a distance is an additional measurement problem.

The change in the affine length of the emitter selected at the original redshift is $-Lb/G$. Since $D'_{A,0}=1$, its additional fractional distance change is $-b/G$. These first derivatives compare source ensembles that coincide at zero amplitude. For finite amplitudes a velocity perturbation must preserve a smooth congruence and the chosen nearby redshift branch. A sufficiently small smooth rapidity realizes every finite additional velocity derivative on this compact sky patch.

## Cancelling the signal

Put $Q(v)=v_1^2-v_2^2$ on the receiving sky. Let $A Q(v)$ be the fiducial fractional distance derivative at fixed redshift. The two caps $v_1\ge4/5$, $v_2\ge4/5$ have opposite mean quadrupoles, and their half-difference is $h=18 A/25$. An additional radial velocity derivative $b(v)=GAQ(v)$ makes the full first-order distance signal vanish. This is a smooth emitter field on the endpoint sphere, at positive radius from the vertex. It extends smoothly to a neighborhood of that sphere, for instance by a radial cutoff. Reversing the sign of a larger correction can reverse the signal.

To erase the two-cap half-difference alone requires $\|b\|_\infty\ge|Gh|$, from the difference of the two normalized cap means. This threshold is attainable with a smooth field equal to opposite constants on the two disjoint closed caps and smoothly interpolated outside them. It is the worst-case threshold for velocity errors of arbitrary angular form. For the normal geodesic emitters it is between $1.7596\times10^{-23}$ and $1.9489\times10^{-23}$ in speed per unit canonical amplitude. For the supported-coordinate ensemble its interval is $(1.0899,1.2067)\times10^{-24}$ in the same units. These values describe the ideal model’s sensitivity to an extra velocity field. They do not specify the accuracy of any instrument.

A common dipole can also change the half-difference between the two **one-sided** caps. For $b(v)=V\cdot v$, their mean radial velocities are $9V_1/10$ and $9V_2/10$, since the axial cosine of a uniform cap direction is uniform on $[4/5,1]$. The least common-dipole speed cancelling the contrast is $\sqrt{2}|Gh|/(9/10)$, attained by $V=(Gh/(9/10))(e_1-e_2)$. It cancels both signed cap means. For geodesic emitters the interval is $(2.7650,3.0623)\times10^{-23}$. Here a common dipole means a single vector in the receiver’s spatial sky frame after background spatial parallel transport; a general coherent flow need not have that angular structure.

## An antipodal control and a velocity gradient bound

For each chosen axis average its positive and negative caps, then difference these two axis averages with the same factor of one half. The fiducial contrast remains $h$ because $Q$ is even under $v\mapsto -v$. Any odd radial velocity pattern, including an arbitrary common dipole, cancels. A common monopole cancels between the two axes. A quadrupolar line-of-sight velocity still survives and can precisely reproduce the signal.

For an independent spatial smoothness bound, let $K_s$ bound the operator norm of the covariant spatial gradient of the extra velocity field on the **entire background emission ball**, using its physical three-metric at $q_1$. Transport velocities along the shorter spatial great-circle arcs through the shell point at angular label $n_*$ on this same slice. The spatial ray direction at emission is the parallel transport of $v$. The two antipodal endpoints are separated by $2a_1\theta$, with $a_1=\rho_b\cosh q_1$ and $\theta=\int_{q_0}^{q_1}\operatorname{sech}q\,dq$. Therefore $|b(v)+b(-v)|/2\le a_1\theta K_s$, and the absolute bias in the fractional antipodal contrast is at most $a_1\theta K_s/|G|$. For the model prediction to be erased at first order, a necessary condition is
$$ K_s\ge |Gh|/(a_1\theta). $$
For the geodesic baseline, the archived conservative $\rho_b$ interval makes this necessary threshold greater than $8.59\times10^{-18}$ in dimensionless speed per unit model length per canonical amplitude. An independently known bound on the physical velocity gradient, in these same conventions, would be required to use this as an exclusion. [Exact threshold replay](velocity.py) · [Rational output](VELOCITY.json)

## What a second distance or length calibration supplies

Write $d_z=L^{-1}\delta D_A|_z$ and $d_\lambda=L^{-1}\delta D_A|_L$. If the log redshift perturbation at the **independently identified affine endpoint** is known, then $d_\lambda=d_z+\delta\log Z|_L/G$. The boost terms cancel identically. Equivalently, if the actual affine length $\lambda_z$ of the redshift-selected endpoint can be independently inferred, the first variation of $[D_A(\lambda_z)-\lambda_z]/L$ is $d_\lambda$, since the derivative of the background residual $D_{A,0}(\lambda)-\lambda$ is zero everywhere. The same holds after angular averaging; for finite amplitude the raywise comparison or the comparison of the two root-mean-squared quantities uses the actual direction-dependent $\lambda_z(v)$. The redshift recorded on a surface already selected at fixed redshift supplies neither independent calibration. These identities leave the fixed-affine null-energy comparison intact and give no upper bound on distance from redshift alone. [Affine normalization and underlying optical comparison](https://github.com/maldonado-research/HDblast/blob/0a47ee1b2ed244a1a43a385a48ae814856fb6f7b/research/HDBLAST_FINITE_LIGHT_CONES_20261118/DERIVATION.md)
