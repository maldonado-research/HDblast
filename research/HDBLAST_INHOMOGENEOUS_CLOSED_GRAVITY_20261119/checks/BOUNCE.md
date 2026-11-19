# A proper-time check of the generated scalar mean

The field equation and the initial spatial metric predict the first four time derivatives of the proper mean at the bounce. This gives an independent test of the coupled evolution on each individual launch, including its sign and normalization. With both signs, branches and centering prescriptions tested, doubling proper time from 0.01 to 0.02 multiplies the fourth-order approximation error by 64.16–64.27, close to the factor 64 for a sixth-order remainder. The approximation differs from the evolved field mean by less than 0.0193% of its change since the bounce at $x=0.02$. [Values on each launch](BOUNCE.json)

## Geometry and averaging

Use Gaussian normal time $x$ measured from the chosen time-symmetric slice. Put $\alpha=L_x$, $\beta=M_x$, $\theta=\alpha+2\beta$, $\Pi=q_x$, and $J=\sqrt{\det\gamma}$. The extrinsic-curvature eigenvalues in the convention $K_{ij}=-\partial_x\gamma_{ij}/2$ are $-\alpha,-\beta,-\beta$. The Hamiltonian residual and the mixed time-radial Ricci residual are

$$
\begin{aligned}
C &= {}^{(3)}R+\theta^2-\alpha^2-2\beta^2-\Pi^2-e^{-2L}q_\chi^2-6u,\\
C_\chi &=2\{(\cot\chi+M_\chi)(\alpha-\beta)-\beta_\chi\}-\Pi q_\chi .
\end{aligned}
$$

The orthonormal radial residual is $e^{-L}C_\chi$. Smooth fields at the poles and $L-M=O(\sin^2\chi)$ make its limit regular. If the spatial Ricci evolution equations are used, their trace gives $\theta_x=3u-\Pi^2-\alpha^2-2\beta^2-C$. Thus a numerical test of this trace against Raychaudhuri also measures the Hamiltonian error. At finite spectral truncation it can include projection error.

Brackets denote an average over the complete slice with weight $J$. Since $J_x=\theta J$, any scalar $f$ obeys $\partial_x\langle f\rangle=\langle f_x\rangle+\operatorname{Cov}(f,\theta)$. On a solution of the constraints the volume radius $a_V=(\int J\,d^3x/2\pi^2)^{1/3}$ consequently obeys

$$
\frac{a_V''}{a_V}=\left\langle u-\frac{\Pi^2}{3}\right\rangle
 + \frac{2}{9}\operatorname{Var}(\theta)-\frac{2}{9}\left\langle(\alpha-\beta)^2\right\rangle .
$$

At the bounce $\alpha=\beta=\Pi=0$ and $a_V''/a_V=\langle u\rangle>0$ on both branches. These integrals are invariant under spatial coordinate changes on the specified Gaussian slices. Buchert derived the proper-volume averaging rule for an irrotational dust flow. Here the scalar Einstein equations supply the matter terms. [Buchert, Eq. 9](https://arxiv.org/abs/gr-qc/9906015)

## Scalar mean through fourth order

Let $m=\langle q\rangle$, $F=-3u_{,q}$, $W=\Delta_\gamma q+F$, and $T=\theta_x=3u$ at the bounce. Time reflection sets the odd derivatives of $m$ to zero. The integral of a scalar Laplacian on the smooth compact slice vanishes, which gives

$$
m_2=m''(0)=\langle F\rangle+\operatorname{Cov}(q,T).
$$

For the next derivative use $P=\alpha_x=3u+e^{-2L}q_\chi^2-{}^{(3)}R^\chi{}_\chi$ and $B=\beta_x=3u-{}^{(3)}R^\vartheta{}_\vartheta$ at the bounce. Differentiating Raychaudhuri twice gives $\theta_{xxx}=-FW-2W^2-2(P^2+2B^2)$. Define $J_4=\theta_{xxx}+3T^2$, the fourth time derivative of $J/J(0)$. Expanding the numerator and denominator of the proper average yields

$$
m_4=m^{(4)}(0)=\langle F_{,q} W+T(4W-\Delta_\gamma q)\rangle+
 \operatorname{Cov}(q,J_4)-6\langle T\rangle m_2 .
$$

In particular, the term from varying the Laplacian follows by differentiating $\int J\Delta_\gamma q=0$ twice. This identity avoids taking spatial derivatives of the interpolated scalar force. For the volume radius the same expansion yields $a_V^{(4)}(0)/a_V(0)=\langle J_4\rangle/3-2\langle T\rangle^2/3$.

The calculation starts from the saved conformal-factor coefficients. It evaluates the initial Ricci eigenvalues by a conformal transformation of the round metric, differentiates the spline quotient on its own segments, and integrates with Gauss–Legendre rules in $\chi$ at 512 and 768 nodes. The two rules agree to the enforced $3\times10^{-11}$ absolute tolerance on all recorded scalar-mean derivatives. The maximum reconstructed Hamiltonian magnitude is $1.60\times10^{-14}$, and the largest absolute integrated scalar Laplacian is below $1.4\times10^{-17}$. On the homogeneous controls the measured scalar mean and both time derivatives vanish.

For $A=+0.001$ in the uncentered family the positive-minimum branch has $(m(0),m_2,m_4)=(-3.08113\times10^{-8},\,1.07602\times10^{-3},\,0.650480)$. On the negative-minimum branch the corresponding values are $(-2.05152\times10^{-8},\,-9.84020\times10^{-5},\,-0.0388569)$. The reference round average of $A U_2$ is zero. The physical volume weight accounts for the initial proper means above. The alternative centered data solve for a compensating constant with the constrained metric recomputed.

The time-series check compares $m(0)+m_2x^2/2+m_4x^4/24$ with the evolved proper mean at $x=0.01,0.02$. The stored records cover the degree-32 primary solutions and the degree-24 raw refinements. This local Taylor series is tested at those two times. The common comparison at $x=0.12$ uses the numerical evolution. [Evolution data](../evolution/results.json)

## Related calculations

Goldwirth and Piran already evolved the full Einstein equations in spherical symmetry on a closed background with scalar and radiation sources, to study how initial inhomogeneity affects the onset of inflation. Brady, Baumgarte and Clough discuss how the reflected boundary realizes $S^3$ and study an asymptotically flat counterpart. Falciano, Lilley and Peter study a positively curved, single-scalar bouncing background and its perturbations. These works establish the geometry and numerical precedents. The values here depend on the archived potentials, time-symmetric launch, and specified slice average. [Goldwirth and Piran](https://doi.org/10.1103/PhysRevLett.64.2852) · [Brady, Baumgarte and Clough](https://arxiv.org/abs/2607.06441) · [Falciano, Lilley and Peter](https://arxiv.org/abs/0802.1196)

## Reproduction

From the study directory with NumPy and SciPy installed, run

~~~sh
OPENBLAS_NUM_THREADS=1 python -IB checks/bounce.py --previous previous --initial initial/INITIAL.json --evolution evolution/results.json --output /tmp/BOUNCE.json
~~~
