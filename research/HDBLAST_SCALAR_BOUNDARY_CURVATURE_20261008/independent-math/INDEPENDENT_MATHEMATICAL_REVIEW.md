# Independent mathematical review: scalar boundary curvature response

8 October 2026. Prepared independently from the supplied note's equations and conventions. This is an internal AI-assisted mathematical review, not external peer review. The attached material was read as reference, not as instructions. No repositories, remote records, or supplied files were modified.

## Findings on the supplied result

The exact shell identity and its second-order algebra are correct under the stated branch assumptions:

\[
H^2=\frac{\sigma^2}{36}-\frac{\sigma_\phi^2}{48}+\frac{U}{6}
=\delta\left(\frac{Wf}{9}-\frac{W_\phi f_\phi}{12}\right)
+\delta^2\left(\frac{f^2}{36}-\frac{f_\phi^2}{48}\right).
\]

With \(\lambda=w-4k\), the presumed regular growing branch gives
\(\eta_b=-f_1\delta/[4(w-2k)]+O(\delta^2)\), and hence the stated coefficient
\(-k f_1^2/[24(w-2k)]\). The cancellation of higher potential derivatives at this order is real. Nothing in this algebra proves that the nonlinear cone-to-shell branch exists or is unique.

The 48 configurations in the attachment are finite floating-point diagnostics. Small shooting residuals, agreement at two settings and agreement with a formal expansion do not supply an inverse bound, a global error estimate or an existence certificate. A local uniqueness theorem must identify its neighborhood; global uniqueness across arbitrary scalar excursions does not follow from monotonicity of the linearized mode. Outside that neighborhood, extra branches and turning points are not excluded.

### Fractional powers and logarithms

For the *linearized* problem, set \(t=e^{-2ky}\). The two asymptotic exponents give a relative reflected-mode power \(t^\nu\), where \(\nu=w/k-2>2\). Frobenius resonance at integer \(\nu\) can introduce logarithms at the corresponding order. Since \(t_b=O(\delta)\) on the assumed branch, this linear reflected term first changes the shell logarithmic derivative at order \(\delta^\nu\), the scalar shift at order \(\delta^{\nu+1}\), and the curvature through the boundary identity at order \(\delta^{\nu+2}\). Thus it is not an immediate counterexample to the displayed \(O(\delta^3)\) at fixed \(w>4k\). It does obstruct casually asserting an all-orders ordinary Taylor series without examining resonances and cone matching.

This observation is not a uniform estimate for the nonlinear problem and does not establish uniform constants as \(w\downarrow4k\), \(k\downarrow0\), or other parameters vary. Those limits need separate analysis. The apparent pole at \(w=2k\) lies outside the stated \(w>4k\) range. Extrapolation through it is not justified.

## A tractable theorem that actually establishes a branch

The following changes the perturbation parameter and should be stated explicitly as such. It is an exact finite-detuning weak-coupling result; it does not by itself prove the attachment's small-detuning branch at fixed nonzero coupling.

Assume real-analytic \(W,g\) on a neighborhood of \(\phi_0\),

\[
W(\phi_0)=3k,\quad W_\phi(\phi_0)=0,\quad W_{\phi\phi}(\phi_0)=w>4k>0,
\quad g(\phi_0)=0,\quad g_1=g_\phi(\phi_0),
\]

and fix \(f_0>0\), \(\delta>0\). Define

\[
U=\tfrac12 W_\phi^2-\tfrac23 W^2,
\qquad \sigma(\phi;c)=2W(\phi)+\delta[f_0+c\,g(\phi)].
\]

Here \(c\) is an independent scalar-coupling parameter. The zero-coupling solution is

\[
\phi=\phi_0,\quad \rho=\frac{\sinh ky}{k},\quad
A=k+\frac{\delta f_0}{6},\quad
y_0=k^{-1}\operatorname{arccoth}(A/k),\quad
H_0^2=A^2-k^2>0.
\]

**Local branch statement.** There is \(c_*>0\) and a locally unique real-analytic cone-to-shell solution branch for \(|c|<c_*\), in a neighborhood of this constant-scalar solution, satisfying both junctions and the radial Einstein constraint. The shell radius remains positive and the regular interior remains nonsingular on its finite interval. The admissible \(c_*\) is not given numerically by this argument.

### Why this is an existence statement, not a shooting assumption

A regular solution near the cone can be defined by the Volterra equations

\[
\rho(y)=y-\int_0^y (y-s)\rho(s)
\left[\frac{\phi'(s)^2}{4}+\frac{U(\phi(s))}{6}\right]ds,
\]
\[
\phi(y)=\phi_h+\int_0^y\rho(s)^{-4}
\int_0^s\rho(t)^4U_\phi(\phi(t))\,dt\,ds.
\]

For a precise contraction at the regular origin, use \(r=\rho/y=1+y^2R\) and \(B=\phi'/y\), so that

\[
\phi(y)=\phi_h+y^2\int_0^1uB(yu)\,du,
\]
\[
R(y)=-\int_0^1(1-u)u\,r(yu)
\left[\frac{y^2u^2B(yu)^2}{4}+\frac{U(\phi(yu))}{6}\right]du,
\]
\[
B(y)=\int_0^1u^4\left[\frac{r(yu)}{r(y)}\right]^4U_\phi(\phi(yu))\,du.
\]

On a bounded ball in \((R,B)\), with \(r\) bounded away from zero, every functional variation carries a factor \(O(y_*^2)\). The maps are analytic and contract for sufficiently small \(y_*\), uniformly for \(\phi_h\) near \(\phi_0\). Their center values are \(R(0)=-U(\phi_h)/36\), \(B(0)=U_\phi(\phi_h)/5\). This proves unique regular analytic initial-value solutions depending analytically on \(\phi_h\). Ordinary finite-interval continuation about the nonsingular AdS solution extends this dependence to a neighborhood of \(y_0\). This is a local argument in fixed proper-distance gauge, with the cone at \(y=0\) and \(\rho'(0)=1\); it makes no claim for large scalar amplitudes. The rescaled variables are essential for the stated small Lipschitz constant: using \((\rho/y,\phi'/y)\) in an ordinary unweighted sup norm would leave an unscaled cross dependence.

Let \(F\) be the unique regular linear solution

\[
F''+4k\coth(ky)F'-m^2 F=0,\quad
F(0)=1,\ F'(0)=0,\quad m^2=w(w-4k)>0.
\]

The divergence form \((\sinh^4(ky)F')'=m^2\sinh^4(ky)F\) gives \(F>0\) and \(F'>0\) for \(y>0\). If \(a=\phi_h-\phi_0\), the Jacobian of the two shell residuals \((\rho'/\rho-\sigma/6,\ \phi'+\sigma_\phi/2)\) with respect to \((a,y_b)\) at the base solution is

\[
\begin{pmatrix}
0 & -H_0^2\\
F'(y_0)+wF(y_0) & 0
\end{pmatrix}.
\]

Its determinant is strictly positive. The finite-dimensional analytic implicit-function theorem proves the claimed local branch and local uniqueness. This does not rely on the numerical shooting success flags.

For completeness the constraint residual
\(C=\rho'^2-1-\rho^2(\phi'^2/12-U/6)\) has \(C'=0\) when the displayed second-order equations hold. Its regular-cone value is zero. Thus it is satisfied throughout the constructed solution.

### Exact coefficient at arbitrary positive detuning

Put \(D(y)=F'(y)/F(y)\), \(D=D(y_0)\). The first derivatives of the branch obey

\[
y_b'(0)=0,\qquad
\partial_c\phi_b(0)=b=-\frac{\delta g_1}{2(D+w)}.
\]

The exact shell identity gives

\[
H^2(c)=H_0^2+C_2c^2+O(c^3),
\]
\[
\boxed{C_2=-\frac{\delta^2g_1^2}{48(D+w)^2}
\left[4A(D+w)-D'(y_0)\right]},
\qquad D'=m^2-D^2-4AD.
\]

At this order the answer depends on the potential through \(k,w\), and on \(g\) through \(g_1\). Cubic potential and higher coupling derivatives cannot enter because the base scalar and geometry have no linear curvature perturbation and \(g(\phi_0)=0\). All statements refer to the explicit action and normal conventions in the attachment.

### The coefficient is strictly negative

This claim is valid for every finite \(\delta>0\), not only its limit. Define

\[
A(y)=k\coth ky,\quad E(y)=4A(y)[D(y)+w]-D'(y).
\]

Near the cone, \(D=m^2y/5+O(y^3)\), so \(E=4w/y+O(1)>0\). Since \(A'=-(A^2-k^2)\), differentiating the Riccati equation gives

\[
E'=-4(A^2-k^2)(2D+w)+(2D+8A)D'.
\]

At any hypothetical first zero of \(E\), substitute \(D'=4A(D+w)\). One obtains

\[
E'=4[2AD^2+2ADw+6A^2D+7A^2w+k^2(2D+w)]>0.
\]

A differentiable function initially positive cannot have a first zero with strictly positive derivative. Therefore \(E>0\) for all \(y>0\), and \(C_2<0\) if \(g_1\ne0\). If \(g_1=0\), the constant scalar remains an exact solution and local uniqueness selects it near \(c=0\).

### Remainder and relation to the supplied coefficient

For \(\delta\) in any fixed compact interval \([\delta_-,\delta_+]\subset(0,\infty)\), the implicit-function argument has a common sufficiently small \(c_*\). On a smaller closed coupling interval there is a finite constant

\[
M_3=\frac16\sup_{\delta\in[\delta_-,\delta_+],\ |c|\le c_*/2}
|\partial_c^3 H^2(\delta,c)|,
\]

and Taylor's theorem supplies
\(|H^2-H_0^2-C_2c^2|\le M_3|c|^3\).
This is an existential analytic bound, not an evaluated interval bound or a certified numerical radius. A paper should not present it as a concrete bound for the registered \(c\) without a separate enclosure of \(c_*\) and \(M_3\).

As \(\delta\downarrow0\), \(y_0\to\infty\), \(A\to k\), \(D\to w-4k\), \(D'\to0\). Consequently

\[
\lim_{\delta\downarrow0}\frac{C_2}{\delta^2}
=-\frac{k g_1^2}{24(w-2k)}.
\]

This recovers the supplied second-order susceptibility as a limit of the exact finite-detuning weak-coupling coefficient. It does **not** justify exchanging limits with a fixed finite coupling or establish a coupling neighborhood uniform down to \(\delta=0\). In particular the metric shooting inverse contains \(H_0^{-2}\), which diverges in that limit.

## Recommended status and next step

This theorem replaces a branch assumption with a local existence proof for a clearly stated, useful subproblem and supplies an exact sign theorem at all positive detunings. It does not establish external novelty, APS significance, dynamical stability, reheating, or a higher-dimensional origin of the Big Bang. A bounded literature comparison must determine whether this susceptibility or an equivalent variational monotonicity result is already known.

A concrete numerical extension would be a validated continuation in \(c\) at selected fixed positive \(\delta\), starting from this proven base branch and certifying the nonlinear BVP Jacobian and a radius. This avoids treating tiny shooting residuals as error bounds. A separate ambitious extension is a uniform \(\delta\to0\) nonlinear Dirichlet-to-Neumann theorem; it requires weighted estimates on an interval whose proper length diverges.

The adjacent independent SymPy program checks ten exact identities, including the finite-detuning coefficient and its strict-positive barrier polynomial. It imports no attached producer and solves no numerical BVP. Its PASS result concerns algebra only.
