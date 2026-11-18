# A finite inward aperture and an angular contrast

On the receiver's spatial sky in five dimensions, a direction near the wall can be written as $\cos\alpha\,v+\sin\alpha\,m$, where $m$ is the inward unit normal of one smooth copy and $v$ belongs to the two-sphere of wall-tangent directions. The bulk sky measure is $\cos^2\alpha\,d\alpha\,d\Omega_2$. Fix a tangential patch $C$ and allow $0\le\alpha\le\beta$. All generators are traced to the same receiver-normalized **affine** length $L$. The interval of sufficiently small positive $\beta$ exists by one-sided geodesic variation. [Geometry and common domain](EMBEDDING.md)

Let $V=(\det\mathcal D_5)^{1/3}$ for the positive vertex Jacobi determinant. For an actual positive aperture define its volume radius as $R_C=(\langle V^3\rangle_{C,\beta})^{1/3}$, where the average divides by the exact solid angle $|C|\int_0^\beta\cos^2\alpha\,d\alpha$. This definition measures three-dimensional screen volume per bulk solid angle. At $\beta=0$ use the one-sided limit. At fixed positive $L$, the background tangent endpoint sphere is embedded, its two angular endpoint derivatives are independent, and the inward inclination derivative has positive normal component. The differential thus has rank three at the edge. Compactness and local injectivity give an embedded thinner one-sided endpoint collar for sufficiently small deformations and inclinations. Ray-labelled volume itself only needs a regular Jacobi map. The individual wall caps in the earlier area-distance experiment were measured with a two-dimensional screen. [Shell beam experiment](../HDBLAST_FINITE_LIGHT_CONES_20261118/DERIVATION.md)

## Dependence on aperture and disturbance

Write $Q(v)=v_1^2-v_2^2$, $f(\lambda)=S(-q)\sin^2\theta$ for the first scalar jet on the positive axis, and $b(\lambda)$ for the first fractional intrinsic area-distance jet on that axis. Throughout this note a dot means a background affine derivative. On a grazing ray the first variation of $k\phi$ is $Q\dot f$ and that of the three optical shear eigenvalues is $Q(\dot b,\dot b,-2\dot b)$. The two intrinsic shear components vanish to first order because the metric tangent is conformal up to diffeomorphism. The background warp metric is conformally flat throughout the regular neighborhood. Choosing $dr/r=dy/\rho(y)$ puts it in the form $(\rho/r)^2[dr^2+r^2ds^2_{\mathrm{dS}_4}]$; the bracket is the flat metric in de Sitter slicing. Consequently every background bulk vertex cone in that neighborhood has zero shear, including inclined cones. [Optical shear and regularity](QUADRATIC_FOCUSING.md) · [Scalar normalization](FLUX.md)

At $\epsilon=\alpha=0$ the Jacobi field for changing inclination has inward component $j_0(\lambda)=\lambda$. The background scalar is constant along the shell and has constant inward derivative $p_{\rm in}=-p_b$. Therefore varying the scalar along the ray with its inclination gives $p_{\rm in}\lambda$. Differentiation in affine length gives $\partial_\alpha(k\phi)=p_{\rm in}$. The two independent first variations yield
$$
 k\phi=\epsilon Q\dot f+\alpha p_{\rm in}+o(|\epsilon|+\alpha),\qquad
 \|\sigma_5\|^2=6\epsilon^2Q^2\dot b^2+o((|\epsilon|+\alpha)^2).
$$
These are one-sided statements at $\alpha=0$, uniform on the compact sky and affine segment. The regular metric and field on one copy give smooth angular geodesic and Jacobi dependence through the rim. At the vertex the integral Jacobi equation gives $\mathcal D/\lambda=I+O(\lambda^2)$, so subtracting $I/\lambda$ from the optical deformation leaves a regular coefficient. Applying the integral focusing equation with $V(0)=0$, $V'(0)=1$, and replacing $V$ by its background value inside the quadratic integrand gives
$$
1-\frac{R_C(\epsilon,\beta)}{L}
=K\epsilon^2\langle Q^2\rangle_C+H\epsilon\beta\langle Q\rangle_C+P\beta^2+
o((|\epsilon|+\beta)^2),
$$
where $W(\lambda)=\lambda(L-\lambda)$ and
$$
K=\frac1{3L}\int_0^L W(\dot f^2+6\dot b^2)\,d\lambda,\qquad
H=\frac{p_{\rm in}}{3L}\int_0^L W\dot f\,d\lambda,\qquad
P=\frac{p_{\rm in}^2 L^2}{54}.
$$
The factors in $H$ and $P$ use $\langle\alpha\rangle=\beta/2+O(\beta^3)$ and $\langle\alpha^2\rangle=\beta^2/3+O(\beta^4)$. Averaging the determinant and taking its cube root agrees with averaging the individual radii at this order because each fractional change starts quadratically. The coefficient $K$ has a lower bound from the scalar flux and a separate positive contribution from shear. A fixed nonzero tilt already carries background scalar flux. [Scalar-flux lower bound](FLUX.json) · [Finite-aperture coefficients](WEDGE_COEFFICIENTS.md)

For this packet $1.17988\times10^{-21}<P<1.43502\times10^{-21}$ and $5.35599\times10^{-23}<H<5.38988\times10^{-23}$, in the canonical normalization. The inward-gradient sign uses $m=-n_{\rm out}$. The deficits of the first and second axial caps have the same $Q^2$ term and opposite $Q$ terms. Their difference, first minus second, is $(36/25)H\epsilon\beta+o(\epsilon^2)$ when $\beta/|\epsilon|$ tends to a positive finite value. This signed comparison and the fourfold comparison below select different parts of the same response. [Exact coefficient enclosures](WEDGE_COEFFICIENTS.json)

## A contrast independent of the overall amplitude

Use caps of cosine radius $4/5$ centered on each of the two equatorial axes and on the two equatorial diagonal directions rotated through $45^\circ$. For each pair average the two *normalized deficits* $1-R_C/L$. The signed mean of $Q$ cancels in each pair. The term $P\beta^2$ cancels in the difference between pairs. The full two-sphere of tangent directions also has $\langle Q\rangle=0$. In its inward band, subtract the deficit of the **unperturbed** band with the same aperture, denoting the result by $D_{\rm sky}$. Denote the difference of the axial and diagonal pair deficits by $D_4$. For $\epsilon\ne0$ tending to zero and $0\le\beta=O(|\epsilon|)$, the expansion gives
$$
 D_4=\tfrac{333}{1250} K\epsilon^2+o(\epsilon^2),\qquad
 D_{\rm sky}=\tfrac4{15} K\epsilon^2+o(\epsilon^2),\qquad
 \lim\frac{D_4}{D_{\rm sky}}=\frac{999}{1000}.
$$
This includes families of positive, shrinking inward apertures. The numerical ratio comes from the cap size and the prepared quadrupole. It changes with angular weighting. The result at fixed macroscopic aperture requires further control of tilted propagation and its higher-order terms.

To verify the coefficients, set $t$ equal to the cosine of the angle from each cap center. It is uniform on $[4/5,1]$. The conditional mean of $Q^2$ is $19t^4/8-7t^2/4+3/8$ around either axis and $2t^2(1-t^2)$ around either diagonal. The cap means are $10283/18750$ and $2644/9375$. On the full tangent-direction sphere $\langle Q^2\rangle=4/15$. Equivalently, $Q^2=\sin^4\theta(1+\cos 4\varphi)/2$ in sky coordinates about the third axis. The difference of the rotated cap pairs isolates its azimuthal fourth harmonic; the azimuth-independent part cancels. The exact cap integrals and an independent Legendre-polynomial calculation agree. [Rational angular replay](angular.py) · [Exact angular moments](ANGULAR.json)

The numerator and denominator are modeled bulk optical volumes. Measuring the ratio requires access to inward directions, the same affine-length and frequency calibration, and sufficient information about the reference band. A wall-confined photon survey supplies only the tangential area. [Access and comparison models](ACCESS.md)
