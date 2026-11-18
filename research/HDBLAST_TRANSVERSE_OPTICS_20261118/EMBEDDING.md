# What an inclined five-dimensional null ray measures next to the shell

Fix one regular interior copy, its timelike shell with induced metric $h$, and its spacelike outward unit normal $n$. Write $g$ for the five-dimensional metric. Use the physical receiver and the fixed affine segment $0\le\lambda\le L$ of the registered shell calculation. The past-directed tangent $k$ satisfies $h(u,k)=1$ at the receiver, where $u$ is its future unit velocity. All curvature traces below are limits from this copy. The construction applies to a sufficiently small neighborhood of zero perturbation amplitude on this common ray domain. [Registered observer, clock and affine length](https://github.com/maldonado-research/HDblast/blob/a011dee2b490faa56ae81d3165b4b5fa5bbfe636/research/HDBLAST_FINITE_LIGHT_CONES_20261118/DERIVATION.md#L3-L19)

## The two tangential directions and the additional normal direction

The sign conventions are $R(X,Y)Z=\nabla_X\nabla_Y Z-\nabla_Y\nabla_X Z-\nabla_{[X,Y]}Z$, $R_{abcd}=g(e_a,R(e_c,e_d)e_b)$, and $K(X,Y)=g(\nabla_Xn,Y)$. On each copy the physical junction law gives $K=\chi h$, with $\chi=\sigma(\phi)/6$. It can vary along the shell. For tangent $X,Y$, the normal component of $\nabla^g_XY$ is $-K(X,Y)n$. In particular $K(k,k)=0$, so an affinely parametrized shell null geodesic also satisfies $\nabla^g_k k=0$ in each one-sided connection, with the same affine parameter. This holds on the actual perturbed shell. [Curvature conventions and Gauss equation](https://github.com/maldonado-research/HDblast/blob/a011dee2b490faa56ae81d3165b4b5fa5bbfe636/research/HDBLAST_INTRINSIC_SHELL_BUDGET_20261101/expansion/GAUSS_IDENTITY_STANDARD_APPLICATION.md#L9-L39) · [Actual junction and copy orientations](https://github.com/maldonado-research/HDblast/blob/a011dee2b490faa56ae81d3165b4b5fa5bbfe636/research/HDBLAST_SUBSIDIARY_CLOSURE_20261116/PHYSICAL_RECONSTRUCTION_AND_BOUNDARY_JOIN.md#L269-L297)

The ambient null screen is the positive three-dimensional space $k^\perp/\langle k\rangle$. It splits orthogonally into the intrinsic two-dimensional screen and $\operatorname{span}[n]$. Parallel-transport two orthonormal intrinsic screen representatives $E_1,E_2$ along the ray. These representatives are also ambient-parallel since $K(k,E_A)=0$. Although $\nabla^g_kn=\chi k$, the class $[n]$ is parallel. If an ambient-parallel vector is desired, $n-(\int_0^\lambda\chi\,ds)k$ represents that same class.

Define the optical matrix by $\mathcal T_{ij}=g(e_i,R(e_j,k)k)$ for the three screen directions. Its sign gives the Jacobi equation $D''=-\mathcal T D$. The Gauss identity gives $R^{(5)}_{E_A k E_B k}=R^{(4)}_{E_A k E_B k}$, because all the extra terms contain $K(k,k)$ or $K(k,E_A)$. For the mixed component the uncontracted Codazzi identity gives

$$
 g(n,R(E_A,k)k)=-(D_{E_A}K)(k,k)+(D_kK)(E_A,k)=0 .
$$

Here $D$ in the covariant derivatives is the intrinsic Levi–Civita connection. The two derivatives of $\chi$ multiply $h(k,k)$ and $h(E_A,k)$, respectively. Symmetry of the curvature tensor supplies the other mixed entry. A null-frame trace now gives an exact scalar in the last diagonal entry,

$$
 T_N=R^{(5)}_{nknk}=\operatorname{Ric}^{(5)}(k,k)-\operatorname{Ric}^{(4)}(k,k)
      =(k\phi)^2-\operatorname{Ric}^{(4)}(k,k).
$$

Indeed the two null-frame terms in the Ricci trace vanish by curvature antisymmetry. The last equality uses the five-dimensional Einstein–scalar equation and tangency of $k$. Equivalently $T_N=E_{kk}+(k\phi)^2/3$, with $E_{kk}=C^{(5)}_{nknk}$ in the repository's Weyl convention. Both normal signs give the same number. [Bulk equation and normal Weyl projection](https://github.com/maldonado-research/HDblast/blob/a011dee2b490faa56ae81d3165b4b5fa5bbfe636/research/HDBLAST_INTRINSIC_SHELL_BUDGET_20261101/expansion/GAUSS_IDENTITY_STANDARD_APPLICATION.md#L33-L83)

At the vertex choose orthonormal angular derivatives of the spatial launch direction at fixed received frequency. The transverse Jacobi matrix has zero value and identity derivative there. Uniqueness of its linear ODE gives, in the parallel screen,

$$
 D_5=\begin{pmatrix}D_4&0\\0&j\end{pmatrix},\qquad
 j''=-T_N j,\qquad j(0)=0,\quad j'(0)=1.
$$

Thus $j$ is the signed normal separation per infinitesimal inclination. At a grazing direction the one-sided three-dimensional beam Jacobian is $j\det D_4$. The usual shell angular-diameter distance is $D_A=\sqrt{\det D_4}$ before conjugate points. The extra Jacobi direction would require probes that can leave the shell.

## Realizing the one-sided inclination on the whole registered segment

In the background the shell has constant intrinsic sectional curvature and $\phi_0$ is constant along it. Both null Ricci contractions therefore vanish, so $j_0(\lambda)=\lambda$. Along the compact family of registered directions the actual curvatures converge uniformly to these background values as the amplitude tends to zero. This uses the strong lower-order metric and scalar regularity and continuous dependence of the geodesic equations. Set $\eta=L^2\sup|T_N|/6$. The integral equation $j(\lambda)=\lambda-\int_0^\lambda(\lambda-t)T_N(t)j(t)\,dt$ bounds $\sup|j/\lambda-1|$ by $\eta/(1-\eta)$ when $\eta<1$. Consequently $j/\lambda>1/2$ on the entire segment for sufficiently small amplitudes, with vertex value one. The two tangential Jacobi directions are already nonconjugate there. [Shell ray domain and Jacobi data](https://github.com/maldonado-research/HDblast/blob/a011dee2b490faa56ae81d3165b4b5fa5bbfe636/research/HDBLAST_FINITE_LIGHT_CONES_20261118/DERIVATION.md#L19-L29) · [One-sided field regularity](https://github.com/maldonado-research/HDblast/blob/a011dee2b490faa56ae81d3165b4b5fa5bbfe636/research/HDBLAST_SUBSIDIARY_CLOSURE_20261116/PHYSICAL_RECONSTRUCTION_AND_BOUNDARY_JOIN.md#L231-L241)

Here is how this Jacobi field is obtained from actual nearby bulk geodesics. Write $m=-n$ for the inward unit normal and $v$ for the shell observer's unit spatial sky direction. Launch with $k_\alpha=-u+v\cos\alpha+m\sin\alpha$ from the same vertex. The frequency normalization and null condition hold for all $\alpha$. Use the registered boundary coordinate $z$ with $z=0$ at the shell and $z>0$ in the interior. In a regular extension of the metric coefficients just across the boundary one may solve the ordinary geodesic IVP for all sufficiently small $|\alpha|$ on the common compact segment. Such a local differentiable extension can be made in collar charts, with Lorentz signature preserved by continuity. It supplies only an ODE neighborhood. For the resulting geodesics $\Gamma_\alpha(\lambda)$ put $F(\alpha,\lambda)=z(\Gamma_\alpha(\lambda))$. We have $F(0,\lambda)=F(\alpha,0)=0$ and

$$
 \partial_\alpha F(0,\lambda)=j(\lambda)\,dz(m)_{\Gamma_0(\lambda)} .
$$

The variation is orthogonal to $k$ and its screen class is $j[m]$. Its possible additional multiple of $k$ is killed by $dz$. The positive number $dz(m)$ is bounded away from zero on the compact shell segment. The quotient $F/(\alpha\lambda)$ extends continuously to both axes by integrating $\partial_\alpha\partial_\lambda F$ over the parameter rectangle. On the $\alpha=0$ edge it equals $(j/\lambda)dz(m)>0$, with its indicated vertex limit. Compactness in $\lambda$ and $v$ gives a common opening angle $\alpha_*>0$. For $0<\alpha<\alpha_*$, the ray lies at $z>0$ throughout $0<\lambda\le L$. Continuous dependence and the original margins keep it inside the physical collar and time slab. These restrictions are consequently solutions in the original interior, independent of the auxiliary extension. The same argument works uniformly after shrinking to a closed, sufficiently small amplitude interval and using collar extensions continuous in that parameter. The opening angle and amplitude threshold have no numerical calibration here.

The identical second copy has the same one-sided result with its own inward normal. The doubled metric can have a jump in its normal derivative. The calculation uses each regular curvature limit on the shell and requires no restriction of a distributional curvature tensor to a gliding generator. It supplies no propagation law for a probe crossing a junction, a global shortcut or a shell-to-bulk emission and detection mechanism. The three-dimensional determinant measures the local change in a pencil at the grazing direction, rather than a finite solid angle integrated into the bulk. [Orientation and regularity of the joined copies](https://github.com/maldonado-research/HDblast/blob/a011dee2b490faa56ae81d3165b4b5fa5bbfe636/research/HDBLAST_SUBSIDIARY_CLOSURE_20261116/PHYSICAL_RECONSTRUCTION_AND_BOUNDARY_JOIN.md#L394-L411)

## Fixed-affine first variation

Let $\delta=\partial_\epsilon|_0$ with the physical receiver normalization and $L$ held fixed. Since $k\phi_0=0$ on the background ray, $\delta\operatorname{Ric}^{(5)}(k,k)=0$ even for the varying ray and tangent. Thus $\delta T_N=-\delta\operatorname{Ric}^{(4)}(k,k)$. Combining the normal Jacobi equation with the two-dimensional Sachs variation yields

$$
 \delta j(L)=\int_0^L t(L-t)\,\delta\operatorname{Ric}^{(4)}(k,k)\,dt=-2\,\delta D_A(L).
$$

In particular the first variation of $j\,D_A^2$ vanishes at fixed affine length. Shell beaming can therefore have a linear response while the extra separation compensates in the full three-dimensional pencil. At a variable central endpoint $L_\epsilon$ with $L_0=L$, the variations of both $j$ and $D_A$ also contain $\delta L$ because both background slopes equal one. For a shell emitter selected at fixed redshift its motion determines this endpoint term. A redshift selected for an inclined bulk ray would additionally require a bulk emitter model. [Affine Sachs variation](https://github.com/maldonado-research/HDblast/blob/a011dee2b490faa56ae81d3165b4b5fa5bbfe636/research/HDBLAST_FINITE_LIGHT_CONES_20261118/DERIVATION.md#L21-L33) · [Fixed-redshift endpoint selection](https://github.com/maldonado-research/HDblast/blob/a011dee2b490faa56ae81d3165b4b5fa5bbfe636/research/HDBLAST_REDSHIFT_EMITTER_MOTION_20261118/FIXED_REDSHIFT.md)
