# Angular-diameter distance at fixed physical affine length

## Observer, rays and sky area

Let $h_\epsilon$ be the registered intrinsic shell metric and $\rho_b$ the background shell radius. The zero-amplitude metric is $h_0=\rho_b^2[-dq^2+\cosh^2 q\,\kappa]$, where $\kappa$ is the unit $S^3$ metric. Elapsed proper time at the material observer is $\rho_b q$ and $q$ increases into the registered past. The receiving observer has initial sphere label $n_*=(0,0,1,0)$, and the vertex is its event at $q_0=\ell a$, $a=49/32$, $\ell=10^{-6}$. Form an actual orthonormal spatial triad from the transported label axes $n_1,n_2,n_4$. Launch each past-directed null ray by its unit sky direction $v=(v_1,v_2,v_4)$ and normalize its affine tangent $k$ by $|h_\epsilon(u_\epsilon,k)|=1$ at reception. Affine parameter $\lambda$ increases into the past and vanishes at the vertex. [Shell metric and observer](../HDBLAST_JOINT_CURVATURE_20261118/METRIC.md)

Set $b=99/64$, $q_1=\ell b$ and fix the physical affine length
$$
 L=\rho_b\frac{\sinh q_1-\sinh q_0}{\cosh q_0}.
$$
This value is common to all rays and amplitudes. In the background the angular path traverses a great circle of the spatial sphere with angle $\theta$. Conservation of angular momentum and the null condition give
$$
 \lambda(q)=\rho_b\frac{\sinh q-\sinh q_0}{\cosh q_0},\quad
 \theta=\arctan(\sinh q)-\arctan(\sinh q_0),\quad
 \sin\theta=\frac{\lambda(q)}{\rho_b\cosh q}.
$$
The perturbed endpoints at $\lambda=L$ can move in time, position and frequency. Fixing $q_1$ in each perturbed geometry would instead introduce a first-order endpoint term.

The terminal background clock lies $1/64$ in scaled time inside the upper limit $25/16$. Compactness of the background rays and continuous dependence of their geodesic and parallel-frame equations therefore keep the perturbed rays in the common shell neighborhood for sufficiently small amplitudes. The screen Jacobi matrix $\mathcal D$ maps vertex angular separations to transverse lengths. It satisfies $\mathcal D(0)=0$, $\mathcal D'(0)=I$ in this physical normalization. The optical tidal matrix vanishes in the constant-curvature background. A double integration of the perturbed matrix equation gives a bound on $\|\mathcal D(\lambda)/\lambda-I\|$ uniform over the sky and the segment when $L^2\sup\|\mathcal T\|$ is small. Thus the actual rays remain free of conjugate points and their area depends differentiably on $\epsilon$ at zero. The area used below is the ray-labelled Jacobian integral. At the fixed nonzero $L$ the background endpoint sphere is also embedded, and compact $C^1$ stability gives embedded nearby endpoint spheres for sufficiently small amplitudes. The radius of that amplitude neighborhood is unquantified. [Field regularity](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_EVOLVED_PARAMETER_20261116/README.md) · [Common shell and observer domain](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_EFFECTIVE_OBSERVER_20261116/EFFECTIVE_OBSERVER_DOMAIN_AND_RESPONSE_BOUNDS.md)

## First variation of the light-pencil distance

Write $\delta=\partial_\epsilon|_{\epsilon=0}$ and $D_A=\sqrt{\det\mathcal D}$ before the first conjugate point, with the vertex limit fixed by $D'_A(0)=1$. A regular light cone has zero twist. In the background $D_A(\lambda)=\lambda$ and the shear vanishes, so the first variation of the Sachs focusing equation is
$$
 \delta D_A(v,L)=-\tfrac12\int_0^L\lambda(L-\lambda)\,\delta\operatorname{Ric}(k,k)\,d\lambda .
$$
The Ricci contraction uses the **actual null tangent** for each metric. It equals $(\operatorname{Ric}-Rh/4)(k,k)$, and that entire tensor vanishes at zero amplitude. Consequently variations of the evaluating event, frame and ray contribute no additional forcing. The variation of the physically normalized Jacobi initial data vanishes as well. [Sachs, gravitational optics (1961)](https://doi.org/10.1098/rspa.1961.0202) · [Kristian and Sachs, observational geometry (1966)](https://doi.org/10.1086/148522)

The shell's metric tangent is $2F h_0$ up to a coordinate Lie derivative, $F=U(q)Y(n)$ with $Y=n_1^2-n_2^2$. The resulting contraction is $-2\,d^2F/d\lambda^2$ along the unperturbed affine ray. Integrating twice by parts gives
$$
 \delta D_A(v,L)=L[F(0)+F(L)]-2\int_0^L F(\lambda)\,d\lambda .
$$
The vertex factor vanishes because $Y(n_*)=0$. Along the great circle, $Y=(v_1^2-v_2^2)\sin^2\theta$. An independent calculation checks the signs using a conformal metric. For $e^{2\epsilon F}h_0$ the ray path is unchanged, its observer-normalized affine measure picks up $e^{2\epsilon F-\epsilon F(0)}$, and its angular-diameter distance at the corresponding event picks up $e^{\epsilon F}$. Varying at fixed physical affine length therefore gives the same endpoint and integrated terms. This check uses only the first derivative of the actual metric. [General Jacobi transport and area distance, §III.4](https://arxiv.org/html/1410.8473#S3.E18) · [Intrinsic curvature tangent](../HDBLAST_JOINT_CURVATURE_20261118/METRIC.md)

Put $M(x)=\mu\rho_b^2 U(\ell x)/\ell$, $H(x)=M(x)\operatorname{sech}^2(\ell x)$, and let $x(t)$ be the reference ray clock corresponding to $\lambda=Lt$. Then
$$
 \frac{\mu\rho_b^4}{\ell L^3}\delta D_A(v,L)=(v_1^2-v_2^2)C,\qquad
 C=H(b)-2\int_0^1t^2H(x(t))\,dt.
$$
The uniform enclosure $H_-\le H\le H_+$ on the full second registered window gives $H_--2H_+/3\le C\le H_+-2H_-/3$. The first replay verifies $0.051613<C<0.107100$. This already fixes the sign independently of the subdivision used below. [Uniform input](../HDBLAST_JOINT_CURVATURE_20261118/METRIC.json) · [Uniform replay and output](OPTICS.json)

## Refining the bound from the waveform cells

The archived calculation supplies a normalized, conjugated scalar input $u(x)$ on continuous clock cells. Its tabulated future intervals are widened on both sides by $3\times10^{-6}$ to enclose the past input. Define the causal double integral $z_0(x)=\int_1^x (x-s)u(s)\,ds$ after the front and set $z_0=z'_0=0$ at $x=1$. On a cell starting at $x_i$, write $\tau=x-x_i$. For each endpoint choice $u_i^\pm$ propagate
$$
 z_i^\pm(\tau)=z_i^\pm+v_i^\pm\tau+\tfrac12u_i^\pm\tau^2,\qquad
 v_i^\pm(\tau)=v_i^\pm+u_i^\pm\tau .
$$
Positive integration kernels bound the actual $z_0$ between these two functions despite sign changes in $u$. Both endpoints of the optical window lie on the archived grid. Set $\Delta=b-a$ and $s=(x-a)/\Delta$. The replay integrates each polynomial exactly against $(x-a)^2 dx/\Delta^3$ on the window to obtain bounds $I_\pm$ on $\int_0^1s^2 z_0(a+s\Delta)\,ds$. It uses every input cell from the front to $b$ to enclose the initial states entering those integrals.

Two errors convert this result into an enclosure of $C$. First, the parent’s Volterra estimate and known time-dependent metric factor give a uniform error $e_M$ between $M$ and $z_0$, and a common bound $Z$ on the magnitude of $z_0$ and $M$ over the interval. The parent report calls the uniform estimate used for its means **absolute_mean_error**. With $X=25/16$ an error for $H-z_0$ is $e_H=e_M+(\ell X)^2Z$. Second, $t=(\int_a^x\cosh(\ell y)\,dy)/(\int_a^b\cosh(\ell y)\,dy)$. Since $1\le\cosh(\ell y)<K=1+10^{-11}$, both $t/s$ by its limiting value at zero and $dt/ds$ lie between $K^{-1}$ and $K$. Hence the relative weight $t^2dt/(s^2ds)$ lies between $K^{-3}$ and $K^3$. The absolute error bound used for the combined endpoint-minus-integral is
$$
 E=\tfrac53 e_H+\tfrac23(K^3-1)Z,\qquad
 z_0^-(b)-2I_+-E\le C\le z_0^+(b)-2I_-+E.
$$
The exact rational enclosure gives $0.0811747<C<0.0826894$. It includes the transfer from the future scalar inputs, the endpoint evolution, geometry and reparametrization of the affine window. These intervals enclose the model calculation. An experimental error budget would require a physical calibration. [Cellwise replay](beam_distance.py) · [Rational endpoints](BEAM.json)

A third calculation avoids enclosing the correlated endpoint and integral independently. By linearity the unperturbed double-integral expression $z_0(b)-2\int_0^1s^2z_0(a+s\Delta)\,ds$ is a single integral against the original acceleration $u(t)$. Its kernel on $1\le t\le b$ is
$$
 {\cal K}(t)=
 \begin{cases}
 \Delta/2+(a-t)/3,& t\le a,\\
 \Delta(1/2-y/3-y^4/6),\quad y=(t-a)/\Delta,&a<t\le b .
 \end{cases}
$$
This follows by integrating $(b-t)_+-2\int_0^1 s^2(a+s\Delta-t)_+\,ds$. The kernel is continuous, positive before $b$ and zero at $b$. On each cell it can therefore be integrated first and multiplied by the lower and upper acceleration endpoints in their original order. Including the transfer, geometry and clock corrections gives the independent central-pencil enclosure **$0.0817743<C<0.0820899$**. The cap and fractional response intervals reported below retain the bounds from the state-propagation replay. [Combined-kernel replay](kernel_check.py) · [Rational results](KERNEL.json)

## From pencils to two finite patches

The positive cap has $v_1\ge 4/5$, solid angle $\Omega=2\pi/5$, and the comparison cap has $v_2\ge 4/5$. Define
$$
 D_{\rm cap}(L)=\left[\frac1\Omega\int_{\rm cap}D_A(v,L)^2\,d\Omega\right]^{1/2}.
$$
The solid angle is fixed on the actual receiving observer’s unit sky. Since $D_A(v,L)=L$ throughout the background sky, $\delta D_{\rm cap}$ is the cap average of $\delta D_A$. On the positive cap $v_1$ is uniform on $[4/5,1]$, and the conditional mean of $v_2^2$ is $(1-v_1^2)/2$. Thus the average of $v_1^2-v_2^2$ is $18/25$ on that cap and $-18/25$ on the other. The full-sky mean vanishes. Multiplying the refined coefficient by $18/25$ gives the positive cap enclosure $0.0584458<(18/25)C<0.0595364$.

Using $\mu=1/(\rho_b B_b)$, the positive cap’s fractional distance response per unit canonical amplitude is
$$
 \frac{\delta D_{\rm cap}}{L}=\frac{18}{25}C\,\frac{\ell B_b}{\rho_b}\left(\frac{L}{\rho_b}\right)^2 .
$$
The scale $L/\rho_b$ lies between $\ell\Delta$ and $K\ell\Delta$ because the affine Jacobian is $\cosh(q)/\cosh(q_0)$. Inserting the archived enclosures for $\rho_b$ and $B_b$ gives the dimensionless interval $(3.7210,\;4.1805)\times10^{-25}$. This is the derivative of a fractional change at zero amplitude. It does not specify the amplitude of an actual source. [Input shell coefficients and checksums](../HDBLAST_JOINT_CURVATURE_20261118/SHA256SUMS) · [Output in both normalizations](BEAM.json)

## Four-dimensional comparison

For an actual regular vertex cone in four-dimensional Einstein gravity with positive Newton constant, matter satisfying the pointwise null energy condition gives $\operatorname{Ric}(k,k)\ge0$. Before any conjugate point on each ray, $D_A>0$ and
$$
 D_A''=-(|\sigma|^2+\operatorname{Ric}(k,k)/2)D_A\le0,\qquad
 D_A(0)=0,\quad D_A'(0)=1 .
$$
Thus $D_A(\lambda)\le\lambda$, and the same ceiling holds for $D_{\rm cap}$ at fixed length. The calculated excess in one cap at sufficiently small nonzero $\epsilon$ excludes a match within this comparison class. Which cap exceeds $L$ depends on the sign of $\epsilon$. The result concerns light constrained to the physical shell metric. It allows different four-dimensional descriptions if their matter violates the condition or if light couples to another metric.

An astronomical application would need identified endpoints, an independent affine-length calibration in the observer’s frequency convention, and source areas or other distance information for those same rays. Extended sources with unknown size, a catalogue at fixed redshift, and beams past caustics require different treatments. For example, suppose a specified emitter velocity field supplies a differentiable redshift $z(\lambda,\epsilon)$ along a chosen direction and $z_0'(L)\ne0$. Implicit differentiation at fixed observed redshift gives
$$
 \left.\delta D_A\right|_z
 =\left.\delta D_A\right|_L
 -\frac{D'_{A,0}(L)}{z_0'(L)}\,\left.\delta z\right|_L,
 \qquad D'_{A,0}(L)=1.
$$
The second term can change the sign. Bounding it requires a source and velocity model independent of the measured distance. This local shell calculation leaves the cosmological origin question open. [Redshift and optical propagation equations](https://arxiv.org/html/2604.05836#S0.E3)
