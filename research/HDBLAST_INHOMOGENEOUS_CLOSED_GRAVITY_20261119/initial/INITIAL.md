# A time-symmetric constrained quadrupole on the three-sphere

Use the stationary point of the archived canonical field and write its displacement as $q$. Set $8\pi G=1$ in the previous study's dimensionless units. The potential entering the stress tensor is $U(q)=3u(q)$, where $u(q)=u_{\rm fit}(\theta_*+q)-u'_{\rm fit}(\theta_*)q$. The fit's residual slope is removed consistently in the energy and the scalar equation. Write $u_*=u(0)$, $a_b=u_*^{-1/2}$, and $\mu^2=-3u''(0)>0$.

On a unit round three-sphere with metric $\Omega$ and radial angle $0\le\chi\le\pi$, prescribe $q=A U_2(\cos\chi)+b$, $\Pi=0$, $K_{ij}=0$, and $h_{ij}=a_b^2\psi^4\Omega_{ij}$. Here $U_n$ denotes a Chebyshev polynomial of the second kind, with $U_2=4\cos^2\chi-1$. The radial round average has weight $(2/\pi)\sin^2\chi\,d\chi$, so $\langle U_n U_m\rangle_\Omega=\delta_{nm}$. Both sides of the momentum constraint vanish. The Hamiltonian constraint $R[h]=|Dq|_h^2+6u(q)$ becomes

$$
-8\Delta_\Omega\psi+6\psi=\psi|Dq|_\Omega^2+6\psi^5u(q)/u_* .
$$

Reflecting through the equator preserves this initial family. On even radial harmonics the linearized operator at $\psi=1,A=b=0$ has eigenvalues $8n(n+2)-24$, all nonzero. This selects the locally unique even near-round branch by the implicit function theorem. The unrestricted linearized operator also has degree-one zero modes. The calculation addresses the near-round even solutions.

## Analytic geometry through second order

For $b=0$, expand $\psi=1+A^2\sum c_n U_n+O(A^3)$. The polynomial identity $U_2^2=U_0+U_2+U_4$ and $\Delta U_n=-\lambda_n U_n$, with $\lambda_n=n(n+2)$, give $|DU_2|^2=8U_2^2+\Delta(U_2^2)/2=8U_0+4U_2-4U_4$. The potential contributes $-\mu^2 U_2^2/u_*$ to the second-order source. Thus

$$
c_n=\frac{8-\lambda_n/2-\mu^2/u_*}{8\lambda_n-24},\qquad n=0,2,4.
$$

| Archived branch | $\mu^2/u_*$ | $c_0$ | $c_2$ | $c_4$ |
|---|---:|---:|---:|---:|
| $\phi_m=+0.5$ | 207.2549 | 8.302287 | −5.081372 | −1.257470 |
| $\phi_m=-0.5$ | 140.9655 | 5.540230 | −3.424138 | −0.862890 |

The proper volume satisfies $V/(2\pi^2a_b^3)=\langle\psi^6\rangle_\Omega=1+6c_0 A^2+O(A^3)$. Proper averages use $\langle f\rangle_V=\langle f\psi^6\rangle_\Omega/\langle\psi^6\rangle_\Omega$. The coordinate quadrupole has $\langle q\rangle_\Omega=0$, whereas its initial proper mean is $6c_2 A^3+O(A^4)$. To impose zero proper mean, solve the Hamiltonian constraint again for each trial constant $b$, then root-find on $\langle q\rangle_V$. The unique local offset begins at $-6c_2 A^3$. For $A=+0.001$ the offsets are $+3.08132\times10^{-8}$ and $+2.05160\times10^{-8}$ for the positive and negative branches. The largest residual proper mean across both centered families is below $7.1\times10^{-20}$.

## Initial proper-time acceleration

Take unit lapse and zero shift, so the normals are geodesic. Write $H_{\rm loc}=\theta/3$ for their expansion scalar divided by three. Initially the expansion and shear vanish. For a canonical scalar $R_{nn}=\Pi^2-3u$, hence Raychaudhuri gives $\dot H_{\rm loc}=u(q)$ on this slice. The scalar gradient contributes to the constraint and cancels from $R_{nn}$. For the whole-sphere volume radius $a_D=(V/2\pi^2)^{1/3}$, differentiation at this instant gives

$$
\frac{\ddot a_D}{a_Du_*}=\frac{\langle u\rangle_V}{u_*}=1-\frac{\mu^2}{6u_*}A^2+O(A^3).
$$

The field has zero normal velocity. For its proper mean, differentiating the volume measure twice and integrating the physical Laplacian over the compact slice gives $\ddot{\langle q\rangle}_V=\langle -3u'\rangle_V+3\operatorname{Cov}_V(u,q)$. For the proper-mean-centered $A=+0.001$ runs this is $+0.00108240$ and $-0.0000955102$, respectively.

| Branch, $A=+0.001,b=0$ | $V/(2\pi^2a_b^3)-1$ | $\min \dot H_{\rm loc}/u_*$ | $\ddot a_D/(a_Du_*)$ | $\langle q\rangle_V$ |
|---|---:|---:|---:|---:|
| $\phi_m=+0.5$ | +0.0000499950 | 0.999685863 | 0.999965338 | −3.08113e−8 |
| $\phi_m=-0.5$ | +0.0000332255 | 0.999788842 | 0.999976517 | −2.05152e−8 |

Both stationary controls have the round metric. All sampled perturbations start at a whole-volume minimum with positive expansion acceleration throughout the slice. This is an instantaneous result in the radial, reflection-even Einstein–scalar family. Radiation, decay and a tensor radiation sector are absent from these data.

## Numerical checks and replay

`constraint.py` solves the projected nonlinear residual with its analytic Newton matrix in even modes through degree 24, using 128 Gauss–Chebyshev nodes. A 1,537-point endpoint-inclusive angular grid checks the pointwise constraint. At $A=+0.001$ the positive-branch residual falls from $5.0\times10^{-6}$ at degree 4 to $8.5\times10^{-15}$ at degree 12. The negative branch goes from $4.9\times10^{-7}$ to $6.6\times10^{-15}$. Across all reported degree-24 data the scaled residual is below $1.1\times10^{-14}$. The global bound $\psi\ge\psi_0-\sum_{n\ge2}(n+1)|\psi_n|$ verifies positivity of the represented solution, including both poles. Increasing to degree 32 and 192 nodes changes the uncentered solution by less than $7\times10^{-17}$ in the global coefficient-based supremum bound. These figures report numerical residuals and differences between truncations.

Averaging the two amplitude signs isolates the leading geometry. Halving $|A|$ reduces the absolute difference between the measured $c_n$ and its second-order prediction by about four in all six comparisons. `INITIAL.json` includes both signs, both amplitudes, controls, angular refinements and proper-mean-centered initial states.

For a standalone replay, extract `HDBLAST_CLOSED_HILLTOP_PREPARATION_20261119.zip` from the previous release and set `PREVIOUS` to its `research/HDBLAST_CLOSED_HILLTOP_PREPARATION_20261119` directory. The ZIP SHA-256 is `6e86ae7026f07cbc7fe48f72d6b257d67b5cb82f16da229df8d749d1483e48e2`. The fit and force routines come from that release. [Previous archive](https://doi.org/10.5281/zenodo.24378057)

```sh
OPENBLAS_NUM_THREADS=1 python initial/constraint.py --previous "$PREVIOUS" --output /tmp/gravity-initial.json
```
