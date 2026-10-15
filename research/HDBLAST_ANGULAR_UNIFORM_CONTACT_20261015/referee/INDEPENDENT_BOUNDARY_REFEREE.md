# Independent referee: angular-uniform mixed boundary evolution

15 October 2026. Internal mathematical review, not external peer review.

## Verdict

The frozen-background linear problem does admit a positive auxiliary energy with constants independent of a full tensorial angular cutoff. The nonsymmetric twelve-dimensional Robin matrix is not an obstruction: a bounded, block-diagonal change of field variables removes it. This conclusion requires a genuine spatial energy, a boundary-adapted splitting, and a justified norm for inhomogeneous boundary data. It does **not** justify an estimate with arbitrary boundary (L^2) Neumann forcing. Glancing persists on the constrained subspace, although it does not destroy homogeneous energy well-posedness.

The exact nonlinear boundary laws also admit a local nonlinear change of field variables that makes their twelve normal rows homogeneous Neumann. I find no algebraic obstruction to that reduction. An exact nonlinear contact theorem still requires a compatible quasilinear existence construction using the already established fixed-finite-regularity exact initial family; it cannot be inferred solely from the linear estimate or the established second coefficient.

## 1. Full direct-metric Robin matrix

Use the original proper-normal (y) component labels, even after replacing the independent radial coordinate by (z=-x=\log(b/r)). Thus the normal metric unknown below is (k_{yy}), not the transformed coordinate component (k_{zz}). Let

\[
 H=(k_{ab})\in\operatorname{Sym}^2(T^*\mathcal B),\qquad
 k=k_{yy},\qquad w=\delta\phi,\qquad
 \ell(H)=h_0^{ab}H_{ab},\quad J=(h_{0,ab}).
\]

There are four Dirichlet components (k_{ya}=d_a) and twelve normal components (U_N=(H,k,w)). Write (a=a_b), (s=\sigma'_0), (q=\sigma''_0), so that (s=-2p_b). The homogeneous part of the complete direct-metric boundary operator is

\[
 \partial_yU_N=B_yU_N,\qquad
 B_y=
 \begin{pmatrix}
 2aI_{10}&aJ&(s/3)J\\
 2a\ell&-4a&4s/3\\
 0&-s/4&-q/2
 \end{pmatrix}.
 \tag{R1}
\]

Here (\ell J=4). This is the linearization of the exact twelve normal laws and agrees with equation (5) of `GH_BOUNDARY_SOURCE_AND_CONSTRAINT_COMPLETION.md`. Since (\partial_z=-\rho\partial_y),

\[
 \partial_zU_N=MU_N+g,\qquad M=-\rho_bB_y,
 \quad g=-\rho_bf_y.
 \tag{R2}
\]

For the pointwise decomposition into trace-free (H) and a trace
\(\theta=\ell H\), with the boundary tensor (h_0|_b) extended constantly in the normal direction for this matrix calculation, the nine trace-free components have eigenvalue (2a) in (B_y). The remaining three-dimensional block, ordered as (\theta,k,w), is

\[
 B_{\rm tr}=
 \begin{pmatrix}
 2a&4a&4s/3\\
 2a&-4a&4s/3\\
 0&-s/4&-q/2
 \end{pmatrix}.
 \tag{R3}
\]

Extending the trace using (h_0(y)^{-1}) instead changes its derivative by (-2a\theta). Those two trace conventions must not be mixed. Neither the scalar Regge–Wheeler two-by-two matrix nor its weight \(\operatorname{diag}(3,8)\) represents (R1).

## 2. Angular principal symbol and the constrained glancing set

After harmless multiplication of the metric rows by (-2) and division by the common conformal factor, all sixteen unknowns have scalar wave principal part

\[
 P_0=\partial_t^2-\partial_z^2-\operatorname{sech}^2t\,\Delta_{S^3}.
 \tag{R4}
\]

Tensor connections, the scalar background gradient, the potential, and gauge-source coefficients give first or zeroth order terms. Angular tensor connections are accounted for with a fixed smooth product bundle connection; replacing a tensor wave operator by a scalar angular Laplacian on arbitrary raw components without those terms would be incorrect.

Freeze a boundary point and set (c=\operatorname{sech}t\). A decaying Laplace–Fourier mode has the form

\[
 e^{st+i\eta\cdot\omega-\lambda z}A,
 \qquad \lambda^2=s^2+c^2|\eta|^2,
 \qquad \Re\lambda>0 \quad(\Re s>0).
\]

The principal boundary symbol is (A_D=0), (-\lambda A_N=0). It is injective for \(\Re s>0\). At limiting glancing, \(s=i\tau\), \(\tau^2=c^2|\eta|^2\), \(\lambda=0\), its Neumann part vanishes. That fact is compatible with the ordinary reflecting-wave energy estimate and is not by itself an instability test.

The harmonic constraint principal symbol does not eliminate this glancing set. Let (\zeta_a\) be the nonzero null tangential covector. With (A_{ya}=0), at \(\lambda=0\), the tangential constraints are

\[
 \zeta^bH_{ba}-\frac12\zeta_a(k+\ell H)=0,
 \tag{R5}
\]

and the normal constraint is automatically zero. The four rows in (R5) have rank four. For example, in a null frame they independently control (H_{uu}), two (H_{ui}), and (k+H_{ii}). The constrained Neumann amplitude space has dimension (12-4=8); in particular the scalar amplitude remains free. Thus any proposed argument that constraints restore an unweighted uniform Neumann boundary determinant is false.

The full frozen normal symbol is (-\lambda I-M). A possible root has bounded \(|\lambda|\leq\|M\|\). Such roots cannot yield growth whose real part tends to infinity with angular frequency: if (s=\gamma+i\tau\), then \(\Re\lambda\geq\gamma\). Taking \(\gamma>\|M\|\) excludes all such roots. This is consistent with a finite exponential energy bound, not a claim of spectral stability.

## 3. Exact linear boundary removal and positive energy

Let (P_D,P_N) denote the fixed normal/tangential field splitting. Choose a smooth compactly supported scalar \(\theta(z)\) with \(\theta(0)=0\) and \(\theta'(0)=1\), and set

\[
 T=I_D\oplus\exp(\theta(z)M(t,\omega)),\qquad U=TV.
 \tag{R6}
\]

The bounded smooth bundle maps (T,T^{-1}) preserve the four Dirichlet rows; (T(0)=I) and (T_z(0)|_N=M). Consequently (R2) becomes exactly

\[
 V_D=d,\qquad \partial_zV_N=g.
 \tag{R7}
\]

One only needs the normal ODE at the boundary. Requiring (T_z=MT) on an infinite cylinder is unnecessary and can create avoidable growth in (T). The compactly supported choice (R6) is enough.

For any scalar-principal wave system, conjugation gives

\[
 T^{-1}P(TV)=P_0V+\widetilde A^\mu\mathcal D_\mu V+\widetilde B V.
 \tag{R8}
\]

The coefficients are bounded on a fixed collar and finite time slab. In particular, angular derivatives of (T) do not create a new second-order principal matrix.

Choose a positive fiber inner product built from the Riemannian product metric (dt^2+dz^2+\kappa_{S^3}), with (P_D,P_N) orthogonal, and a product connection \(\mathcal D\) preserving the splitting. This is an analytic control norm, not the Lorentzian metric contraction on tensor components and not a physical gravitational energy. For homogeneous boundary data,

\[
 E_0[V]=\frac12\int_{z\geq0,S^3}
 (|V_t|^2+|V_z|^2+c(t)^2|\mathcal D_\omega V|^2+|V|^2)
\]

satisfies

\[
 E_0'(t)\leq C E_0(t)+\|\widetilde F(t)\|_{L^2}\sqrt{2E_0(t)}.
 \tag{R9}
\]

The boundary flux vanishes because (V_{D,t}=0) and (V_{N,z}=0). No symmetry of (M) has been used. The constants depend on the fixed background slab, (T,T^{-1}), and finitely many coefficient norms, not on a harmonic label.

For high regularity, commute time and angular derivatives using this product connection. They preserve the homogeneous boundary rows. Their commutators contain only derivatives controlled at the same energy order. Recover two and more normal derivatives from (R8), whose normal coefficient is uniformly nonzero. Commuting the spacetime tensor connection or arbitrary normal derivatives first and simply declaring the resulting boundary conditions homogeneous would leave a gap.

All-angular estimates may be established directly. To express them as cutoff-uniform estimates, use full tensorial SO(4) isotypic projections: the rotational background, boundary splitting and transformation are equivariant, so these projections commute with the system. Componentwise scalar harmonic truncations do not supply this property.

## 4. A conservative, explicit source space

An (L^2) bulk energy does not control (V_t|_{z=0}). Hence the inhomogeneous flux \(\langle g,V_t\rangle_{z=0}\) cannot be bounded by an unproved elementary trace estimate. A safe, deliberately nonoptimal construction is a smooth normal lift

\[
 R_D=\chi(z)d(t,\omega),\qquad
 R_N=z\chi(z)g(t,\omega),\qquad \chi=1\text{ near }0.
 \tag{R10}
\]

Then (W=V-R) satisfies homogeneous mixed boundary conditions and is forced by \(\widetilde F-\widetilde P R\). For integer (m\geq2\), define

\[
 X_m(V;t)=\sum_{j=0}^{m+1}\|\partial_t^jV(t)\|_{H^{m+1-j}(\Sigma)},
\]

\[
 \mathcal F_m(F)=
 \sum_{j=0}^{m}\|\partial_t^jF\|_{L^1_tH^{m-j}}
 +\sum_{j=0}^{m-1}\|\partial_t^jF\|_{L^\infty_tH^{m-1-j}},
\]

\[
 \mathcal B_m(d,g)=\sup_tX_m(R;t)+\mathcal F_m(\widetilde P R).
 \tag{R11}
\]

These are concrete finite norms of the displayed lift, not an assertion of optimal boundary regularity. A simpler sufficient data class is bounded/continuous mixed time–angular derivatives of (d,g) through total order (m+2), with compact angular manifold and finite time interval. For example, (\partial_t^j(d,g)\in C_tH^{m+2-j}(S^3)) for (0\leq j\leq m+2\) suffices for (R11).

Tangential energy, normal recovery, and Gronwall give the conservative estimate

\[
 \sup_{0\leq t\leq T}X_m(V;t)
 \leq C_{m,T}\left[
 \sum_{j=0}^{m+1}\|V_j\|_{H^{m+1-j}}
 +\mathcal F_m(\widetilde F)+\mathcal B_m(d,g)\right],
 \tag{R12}
\]

where (V_j=\partial_t^jV|_{t=0}) are the recursively computed initial jets, not independently prescribed data. The normal-recovery part is the reason for retaining the displayed (L^\infty) lower-order source terms. Existence for smooth compatible data follows by the homogeneous mixed wave construction with bounded lower-order perturbations; the energy gives uniqueness and convergence of compatible approximations. Completion in the finite norm requires the usual compatibility trace conditions below.

Glancing confirms why the boundary norm qualification is material. In the constant scalar model with boundary forcing, the transformed solution is \(\widehat V=-\widehat g e^{-\lambda z}/\lambda\). Its bulk energy density integrated in (z) contains

\[
 \frac{|s|^2+c^2|\eta|^2+|\lambda|^2}
 {2\Re\lambda\,|\lambda|^2}|\widehat g|^2.
\]

For (s=\gamma+iN\), (c|\eta|=N\) and fixed positive \(\gamma\), this grows like \(N^{1/2}\gamma^{-3/2}|\widehat g|^2\). A uniformly bounded (L^2)-boundary-to-energy map is therefore not established by homogeneous flux cancellation; it fails in this model. The lift norm avoids claiming such a map.

## 5. Finite corner conditions and subsidiary uniqueness

Compute the jets (U_0,U_1,U_2,\ldots\) recursively from the wave equation, differentiating its time-dependent coefficients and source. To obtain (U\in\bigcap_{j=0}^{m+1}C^j_tH^{m+1-j}\), the classical trace compatibility conditions are

\[
 P_DU_j|_{z=0}=\partial_t^jd(0),\qquad 0\leq j\leq m,
\]

\[
 \partial_zP_NU_j|_{z=0}
 -\sum_{r=0}^j\binom jr(\partial_t^rM)(0)P_NU_{j-r}|_{z=0}
 =\partial_t^jg(0),\qquad 0\leq j\leq m-1.
 \tag{R13}
\]

These ranges follow from the value trace threshold (H^{1/2+}) and normal trace threshold (H^{3/2+}). Smooth solutions require the recursion at every order. Vanishing data and sources in an open initial shell neighborhood give all of them automatically. A statement of merely “compatible data” should identify the chosen regularity and recursion.

The linearized subsidiary covector obeys a scalar-principal wave equation with one Dirichlet row (C_y=0) and four Robin rows \(\partial_yC_a-2a_bC_a=0\). The same transformation and positive energy prove all-angular uniqueness for its zero-data problem. This eliminates an incoming constraint mode once (C|_{t=0}=C_t|_{t=0}=0) has actually been established from the complete initial gauge and physical constraint data. Codazzi and the scalar junction supply the tangential subsidiary rows; they are not extra primal conditions.

## 6. Independent assessment of the exact nonlinear state-flow proposal

Subtract the background field values and write the exact twelve normal laws at (z=0) as

\[
 \partial_zU_N=\mathfrak B(t,\omega,U_N),\qquad \mathfrak B(t,\omega,0)=0.
\]

The right side is smooth while the normal lapse is positive and the induced metric is invertible Lorentzian. It depends only on (U_N), not on normal derivatives or tangential derivatives of the four Dirichlet fields. Let \(\Phi(s,t,\omega,v)\) be the local finite-dimensional flow of \(\mathfrak B\), and define

\[
 U_N=\Phi(\theta(z),t,\omega,V_N),\qquad U_D=V_D,
 \quad \theta(0)=0,\quad\theta'(0)=1.
 \tag{R14}
\]

Then \(U_{N,z}=\mathfrak B(U_N)+D_v\Phi\,V_{N,z}\) at the boundary, so the exact nonlinear rows become \(V_{N,z}=0\). The flow fixes zero and its Jacobian is invertible for a sufficiently small field neighborhood and normal collar. A cutoff in \(\theta\) gives a localized smooth change of variables.

Every reduced Einstein–scalar field has the same scalar principal coefficient \(g^{\alpha\beta}(U)\). Applying the chain rule and multiplying by the inverse state Jacobian therefore preserves that scalar second-order principal part. All extra terms have at most first derivatives of (V), quadratically. The four Dirichlet rows still imply \(g^{za}=0\) at the boundary, for the solution and for any frozen iterate satisfying those rows. Thus an appropriate tangential time multiplier has zero boundary flux for the transformed homogeneous mixed conditions.

This removes the most obvious nonlinear derivative-loss mechanism arising from lifting a differentiated Robin source. It is a legitimate route to a conditional quasilinear theorem. The following points remain obligations of that theorem:

1. A regularity index, e.g. an integer (s\geq5) in four spatial dimensions, and a compatible mixed-boundary iteration/approximation. Naive iterates with the same initial displacement and velocity need not satisfy higher corner compatibility for their separately frozen equations. Tangential/time energies and normal recovery must be included.
2. Reuse the exact fixed-finite-regularity initial family established in `EXACT_LOCALIZED_INITIAL_DATA_AND_PRECONTACT.md`: for fixed integer (k\geq12), its metric and metric momentum have (H^{k+2}\times H^{k+1}) regularity, are parameter-smooth in that topology, and equal the background near the shell. This already resolves the corresponding bulk initial integrability question under its stated reflection/localization hypotheses. Complete its initial harmonic gauge and nonlinear mixed-boundary compatibility in the actual boundary chart. Background agreement on an open shell collar supplies all defined corner identities. A second-order polynomial is not a substitute for this exact family, and that family's second correction need not equal the previously selected frozen coefficient correction.
3. Preservation of Lorentzian signature, a timelike fixed boundary and the required coefficient bounds over the target finite slab. Small-data continuation around a smooth background is plausible but must be quantified to reach beyond the chosen first contact time.
4. Parameter differentiability in stated spaces. Quasilinear solution maps generally lose derivatives upon parameter differentiation. Two formal variation equations suggest two losses for the variation values, but that alone does not prove convergence and continuity of difference quotients at those orders. The established short-precontact corollary conservatively reserves four Sobolev losses for actual (C^2) parameter dependence; retain that margin or prove a sharper mixed-boundary comparison theorem. The fixed (k\geq12) initial family supplies substantial finite regularity. Do not promote this family to spatial (C^\infty) on one parameter interval, or import its short precontact lifespan as a lifespan beyond prescribed contact.
5. The nonlinear subsidiary system has solution-dependent Robin coefficients. Once the primal solution has the asserted regularity, the linear positive-energy argument applies to those fixed coefficients, but initial subsidiary vanishing still needs its separate constraint argument.

I recommend recording the frozen-background uniform theorem separately from these nonlinear obligations. The previously established finite-angular second coefficient remains unchanged in either case.
