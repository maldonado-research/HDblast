# Can an infalling observer cross the throat after a null wave passes?

A sufficiently short pulse can leave an open neighborhood of every finite point of the future Poincaré horizon free of the wave. For a wave that fades slowly, an infalling observer can instead encounter unbounded tidal curvature in finite proper time. The distinction depends on the tail measured along the observer's path. There is a second restriction on the earlier plane-wave example. A pulse of finite duration that remains a perfectly collinear source-free plane-front wave cannot also be localized across an infinite planar wall.

## The exact null sector

Use the tuned branch with wall-width unit one, warp factor $a=e^{A(y)}$, and $A(0)=0$. On the retained half-line $y\ge0$,

$$
\phi=-\tanh y,\quad A=-\frac13\left[y+\frac23\log\cosh y+\frac16\tanh^2 y\right],\quad
k=-A_y(\infty)=5/9.
$$

The previous study added a Maxwell sector and imposed a fixed warp and scalar during a separated null electromagnetic wave. That requirement selects $G=C a^2$. Here $C>0$ and the four-dimensional Einstein–Maxwell normalization $\lambda>0$ are constants. The exact plane-front example has

$$
ds^2=dy^2+a^2\left[-2\,du\,dv+dx_1^2+dx_2^2
       -\frac{\lambda}{2}E(u)^2(x_1^2+x_2^2)du^2\right],
\qquad F=E(u)\,du\wedge dx_1 .
$$

The four-dimensional null Ricci component is $\hat R_{uu}=\lambda E(u)^2$. The two transverse curvatures are $\hat R_{uiuj}=\lambda E^2\delta_{ij}/2$. Our sign convention gives $R_{TiTi}=+k^2$ for an orthonormal spacelike screen direction in pure AdS and a unit timelike tangent $T$. The shell stays fixed, satisfies the tuned matching conditions, and has no added electromagnetic surface sources. The disturbance travels *along* the wall. It serves as a regularity control for a proposed higher-dimensional signal. [Earlier solution and assumptions](https://github.com/maldonado-research/HDblast/blob/7fd4e3066499ed0e5264ec20a1bd97b634fa40ae/research/HDBLAST_BULK_PHOTON_PROPAGATOR_20261120/bounds/BOUNDS.md) · [Domain-wall plane waves](https://arxiv.org/abs/hep-th/9909130) · [Related horizon tidal tests](https://arxiv.org/abs/hep-th/0012076)

## The observer and the tail

On the axial timelike geodesic $x_1=x_2=0$, let $s$ denote proper time, $P=a^2du/ds>0$ and $B^2=2(a^2du/ds)(a^2dv/ds)$. Then

$$
\frac{dy}{ds}=\sqrt{B^2/a^2-1},\qquad
\frac{du}{dy}=\frac{P}{a\sqrt{B^2-a^2}},\qquad
\frac{ds}{dy}=\frac{a}{\sqrt{B^2-a^2}}.
$$

Here $B>1$ gives a geodesic entering the throat from $y=0$. The horizon lies at finite proper time $s_h$. Since $A_y\to-k$, the leading relations are $s_h-s\sim a/(kB)$ and $u\sim P/(kBa)$. More generally the projected timelike geodesic on the four-dimensional slice has an affine parameter satisfying $d\hat s/ds=a^{-2}$, a constant timelike norm, and constant $du/d\hat s$ because $du$ is parallel in the slice. Each future timelike infaller reaching the throat in this geometry samples $u\to+\infty$.

The screen vectors $e_i=a^{-1}\partial_{x_i}$ on the axial geodesic are parallel along $T$. Direct contraction gives

$$
R_{TiTi}=A_y^2-A_{yy}(dy/ds)^2 +\frac{\lambda P^2 E(u)^2}{2a^4},\qquad
R_{TT}=4(A_{yy}+A_y^2)-3A_{yy}B^2/a^2+\frac{\lambda P^2 E(u)^2}{a^4}.
$$

The wave contributes to these freely falling tidal measurements. The scalar invariants built by contracting powers of the curvature without derivatives can meanwhile retain their background values for this null plane wave. A fixed-position observer with a frame rescaled from the slice measures null stress proportional to $a^{-2}E^2$ at fixed $u$. Holding $u$ fixed when sending $y\to\infty$ follows a different path from any one of the timelike infallers above.

For $E(u)\sim E_0 u^{-p}$ with $E_0\ne0$, the added Ricci signal scales as $(s_h-s)^{2p-4}$. The ordinary Maxwell component $|F(T,e_1)|=P|E|/a^3$ scales as $(s_h-s)^{p-3}$. Thus $p<2$ gives unbounded freely falling curvature and obstructs a regular $C^2$ metric extension along this observer. At $p=2$ the wave's tidal signal tends to a finite, generally nonzero limit. At $p>2$ it vanishes. The ordinary two-form needs $p\ge3$ for a bounded component. The canonical stress includes $G\propto a^2$, which explains the two different thresholds. These are results for the nonoscillatory power asymptote stated here. General extendibility requires more than bounded curvature along one trajectory. In particular, a finite flux integral $\int E^2 du$ per unit transverse area permits $p>1/2$ and leaves open the tidal divergence.

In the tested family $E=(1+u^2)^{-p/2}$ the measured powers are as follows. The fits use proper throat depths $y=28,32$. Errors in the displayed powers are below $3\times10^{-7}$. The separate geodesic calculation evolves the full axial proper-time system to $y=4,8,12$, with elapsed-time and conserved-momentum checks.

| Amplitude tail $p$ | Ricci wave signal power | Measured ordinary field power |
| ---: | ---: | ---: |
| 0, eternal constant field | −4 | −3 |
| 1 | −2 | −2 |
| 1.5 | −1 | −1.5 |
| 2 | 0 | −1 |
| 2.5 | +1 | −0.5 |
| 3 | +2 | 0 |

For the compact pulse use $E(u)=\exp(1-1/[1-(u-2)^2])$ on $1<u<3$ and zero outside. With $B=1.4$, $P=0.9$, $\lambda=1$, the added Ricci signal is 26.4241 at the pulse center. Its support ends at $s=1.17321$, leaving 0.40211 of the observer's proper time before the horizon at $s_h=1.57532$.

## A local continuation across finite horizon points

In the vacuum region write $r=a$. The metric differs from a Poincaré patch of AdS of radius $\ell=1/k$ by a single term,

$$
g_{\rm vac}=g_{\rm AdS}+b(r)\,dr^2,\qquad
b(r)=\frac{1/A_y^2-1/k^2}{r^2}.
$$

On the actual wall, $A_{yy}=-\mathrm{sech}^4y/3$, $A_y+k=O(e^{-4y})$, and $r\sim c e^{-ky}$. Thus $b=O(r^{4/k-2})=O(r^{26/5})$ and $\phi+1=O(r^{2/k})=O(r^{18/5})$. Differentiating the convergent expansions in $e^{-2y}$ near the throat gives the corresponding derivative orders. Background boosted curvature anisotropy $|A_{yy}|/r^2$ tends to zero as $r^{26/5}$. The first and second derivatives needed for a $C^2$ continuation of the vacuum metric and scalar exist.

The function $r$ is a regular embedding coordinate on the extension of AdS. In light-cone coordinates set $U=r u$, $V=r v$, $X_i=r x_i$, $\ell r=X^{-1}+X^4$, $W=X^{-1}-X^4$. The hyperboloid identity gives $-2UV+X_1^2+X_2^2=-\ell^2+\ell r W$. At any finite horizon point $r=0$, $U$ is nonzero. Consequently a neighborhood of that point in $r>0$ has $u=U/r$ outside the support of any compact sandwich. Locally the wave vanishes identically. One may continue the vacuum metric and scalar by taking $b=0$ and $\phi=-1$ for $r<0$. The limiting scalar is an extremum of the archived potential and these fields match to at least second order. This constructs a local Einstein–scalar continuation at finite points of the AdS horizon. It imposes no global boundary condition at conformal infinity and establishes neither geodesic completeness nor uniqueness of the continuation.

One can do more for the particular algebraic tails used in the table. Write $X_\perp^2=X_1^2+X_2^2$. For $r>0$ their additional metric and field become

$$
\delta g=-\frac{\lambda X_\perp^2}{2}\frac{r^{2p-4}}{(r^2+U^2)^p}(r\,dU-U\,dr)^2,
$$

$$
F=\frac{r^{p-3}}{(r^2+U^2)^{p/2}}
       (r\,dU\wedge dX_1-X_1\,dU\wedge dr-U\,dr\wedge dX_1).
$$

At $p=2$ the additional metric has a smooth local expression through $r=0$ even though $F$ diverges. At $p=3$ both extra tensors have smooth expressions. Continuing the explicit signed powers to $r<0$ and taking $b$ and $\phi$ from the mirrored leaf realizes the same local metric and null-field solution on each open side. The numerical chart tests eliminate $V$ using the hyperboloid equation, check the Lorentzian signature and determinant, and verify the invariants when the field is finite. In all cases the gauge weight tends to zero at the horizon and $G_{,\phi}$ diverges at the limiting vacuum. Existence of a local tensor continuation in this restricted sector does not supply a regular interacting gauge theory there or determine how arbitrary charged perturbations cross.

## Why a sandwich field cannot also be a finite collinear beam

Keep the flat transverse plane and the common parallel null direction of the separated wave, but allow $F=du\wedge(e_1 dx_1+e_2 dx_2)$ with $e_i(u,x_1,x_2)$. At each $u$, the Bianchi identity requires $\partial_1e_2-\partial_2e_1=0$ and the source-free Maxwell equation requires $\partial_1e_1+\partial_2e_2=0$. This follows directly from $|\det \hat g|=1$ and $\hat g^{uv}=-1$, for arbitrary plane-front function $H(u,x_1,x_2)$. The radial weight does not change these tangential constraints. In transverse Fourier space the amplitude must be both perpendicular and parallel to each nonzero wave vector. A square-integrable transverse field supported in Fourier space only at zero is zero. The exact constant-screen wave has infinite transverse area. Compact duration resolves its horizon-tail problem, but it does not localize that beam.

There is a quantitative version. For a smooth transverse envelope of finite norm and second moment, center its transverse coordinates and define

$$
R^2=\frac{\int |x|^2|\mathbf e|^2}{\int |\mathbf e|^2},\qquad
D^2=\frac{\int (|\mathrm{div}\,\mathbf e|^2+|\mathrm{curl}\,\mathbf e|^2)}{\int |\mathbf e|^2}.
$$

Plancherel gives the numerator of $D^2$ as $\int|\nabla\mathbf e|^2$. Integration by parts in two dimensions gives $\int|\mathbf e|^2=-\mathrm{Re}\int \mathbf e^*\!\cdot (x\!\cdot\!\nabla)\mathbf e$. Cauchy–Schwarz yields **$DR\ge1$**. A circular Gaussian envelope attains equality. Doubling its RMS width halves the least RMS violation of the two perfectly aligned Maxwell equations. The code independently differentiates isotropic and elliptical envelopes on a transverse Fourier grid and recovers products $1,1,1.25,2.125$. This bound is an application of the two-dimensional uncertainty inequality in the stated normalization. It does not bound solutions with additional field components or with diffraction.

Simply replacing the transverse plane with a flat periodic two-torus also has a limitation in this ansatz. The Einstein equation for a globally periodic $H$ would read $\Delta_\perp H=-2\lambda|\mathbf e|^2$. The integral of a periodic Laplacian is zero and $\lambda>0$, so a nonzero null flux is excluded. Localized or compact higher-dimensional source models must relax a stated condition, for example by allowing angular spread and a responding spacetime. [Calculation](horizon.py) · [Numerical results](HORIZON.json)
