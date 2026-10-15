# Angular-uniform control of the registered linear contact problem

Internal mathematical checkpoint. Artifact identifiers and machine receipts use UTC; this work began 14 October 2026 Pacific time.

## 1. Scope and statement

This note concerns the **full linearized reference-harmonic Einstein–scalar operator at the registered background** on a fixed smooth positive-cone collar. It is not a nonlinear contact existence theorem. The previous selected second-coefficient construction remains unchanged.

The model is
\[
 E_{AB}=R_{AB}-\phi_A\phi_B-\tfrac23U(\phi)g_{AB},\qquad S=\Box_g\phi-U'(\phi),
\]
with background
\[
 g_0=dy^2+\rho(y)^2[-dt^2+\cosh^2t\,\kappa_{S^3}],\qquad\phi_0=\phi_0(y).
\]
Use the relative wave-map vector \(C^A=g^{BC}(\Gamma[g]^A{}_{BC}-\Gamma[g_0]^A{}_{BC})\) and reduce by \(E_{AB}-\nabla_{(A}C_{B)}\), with the factor \(1/2\) in symmetrization. All fifteen metric and one scalar unknown share the scalar wave principal part after multiplication of their equations by nonzero background factors. Tensor connection and matter terms are at most first order in the unknowns.

Set \(z=\log(b/r)\), inward from the shell, so \(dy=-\rho\,dz\). Work on \(0\le t\le T\), with support inside a finite controlled collar \(0\le z\le Z\); extend background coefficients smoothly beyond that support only as a construction device. \(T,Z\) are fixed and the background is smooth with \(\rho>0\) there. If the initial/source radial support is bounded by \(z\le L+t\), choose \(Z>\max(L,\delta)+T\), where \(\delta\) contains the fixed transformation and lifting collars. This causal margin keeps the extended coefficients outside the domain of dependence. Constants below may depend on these intervals and this fixed background, but not on angular frequency or an angular cutoff.

**Linear control theorem.** With homogeneous coordinate Dirichlet data, homogeneous physical Robin data, and the normal relative-harmonic boundary condition, the complete reduced linear system has a positive auxiliary \(H^1\times L^2\) control estimate uniform over all angular modes. Smooth, fully compatible data and forcing give a unique smooth solution on the controlled collar. Inhomogeneous boundary data are permitted through the explicit smooth lifting norm in section 5. No estimate for arbitrary \(L^2\) Neumann boundary forcing is asserted. The same control mechanism gives uniqueness for the subsidiary vector system.

The norm is an auxiliary analytic control norm. It is not a physical gravitational energy, a ghost criterion, or a spectral stability statement.

## 2. The complete boundary matrix

Let \(a_b=\rho_y/\rho|_b\). The common physical copy-even laws are
\(K_{ab}=\sigma h_{ab}/6\), \(n\phi=-\sigma'/2\), with \(a_b=\sigma_0/6\). At the shell the exact reference-gauge completion is
\[
 g_{ya}=0,\qquad h_{ab,y}=\sigma N h_{ab}/3,\qquad
 \phi_y=-\sigma'N/2,\qquad
 N_y=2\sigma N^2/3-a_bN^3\operatorname{tr}_{h}h_0,
\]
where \(g_{yy}=N^2\). Linearize at \(N=1,h=h_0,\phi=\phi_0\), using direct metric unknown \(k=k_{yy}\), tangential metric \(k_{ab}\), and scalar \(w\). With \(\tau=\operatorname{tr}_{h_0}k_{ab}\), the twelve Robin rows are
\[
 \partial_y k+4a_bk-2a_b\tau-\tfrac43\sigma'_0w=f_k,
\]
\[
 \partial_y k_{ab}-2a_bk_{ab}
       -a_bh_{0,ab}k-\tfrac13\sigma'_0h_{0,ab}w=f_{ab},
\]
\[
 \partial_y w+\tfrac12\sigma''_0w+\tfrac14\sigma'_0k=f_w.
 \tag{1}
\]
The four \(k_{ya}\) rows are Dirichlet. The same linear left sides apply to the second Taylor coefficient; its right sides are the previously verified complete quadratic sources, including \(k_{yy}=2\nu_2+2\nu_1^2\). Do not substitute \(2\nu_2\) for the metric wave unknown.

Define \(B_y\) by the zeroth-order terms in (1). Inward coordinates give
\[
 -\partial_zU_N+M(t)U_N=g,\qquad M=\rho_bB_y,\quad g=\rho_b f.
 \tag{2}
\]
Here \(U_N=(k,k_{ab},w)\) has rank twelve. This full matrix generally is not symmetric for a naive positive component metric. The restricted RW two-channel symmetrizer is not needed and is not assumed to define a closed physical bulk system.

Use the global field bundle
\[
 E_D=dy\odot T^*(\mathbb R_t\times S^3),\qquad
 E_N=\mathbb R\,dy^2\oplus\operatorname{Sym}^2T^*(\mathbb R_t\times S^3)\oplus\mathbb R.
\]
It has ranks four and twelve. Give it the fixed positive product fiber metric obtained from \(dt^2+\kappa\), with the two subbundles orthogonal, and the product connection preserving the split. Replacing the physical tensor connection by this connection changes only bounded first- and zeroth-order bulk terms. The actual \(M\) is a smooth SO(4)-equivariant endomorphism of \(E_N\). It is also angular-parallel for this fixed product connection: \(h_0(t)\), its inverse trace and the scalar coefficients are parallel along the round sphere.

## 3. Removal of the nonsymmetric Robin matrix

Choose a smooth scalar cutoff \(\chi(z)\) supported in a fixed thin collar, with \(\chi(0)=0\), \(\chi'(0)=1\). Define
\[
 T_D=I,\qquad T_N(t,z)=\exp[\chi(z)M(t)],\qquad U=TV.
 \tag{3}
\]
Then \(T(0)=I\), \(\partial_zT_N(0)=M\), and \(T\) and \(T^{-1}\), together with each fixed finite number of derivatives, are bounded on the controlled interval. No commutativity of \(M\) with its time derivative, no diagonalizability, and no boundary energy symmetry are required. Away from the thin collar \(T=I\).

The exact boundary calculation is
\[
 -\partial_z(T_NV_N)+M T_NV_N=-\partial_zV_N\quad(z=0).
\]
Thus homogeneous (2) becomes \(\partial_zV_N=0\); Dirichlet conditions remain \(V_D=0\). Inhomogeneous data become \(-\partial_zV_N=g\), \(V_D=d\).

Because the second-order principal symbol is scalar on the full bundle, conjugation does not introduce a new matrix principal symbol. In product-connection form the transformed operator is
\[
 P_TV=P_0V+A^\alpha(t,z,\Omega)D_\alpha V+B(t,z,\Omega)V,
 \quad P_0=\partial_t^2-\partial_z^2-c(t)^2\Delta_\kappa,
 \quad c(t)=\operatorname{sech}t.
 \tag{4}
\]
The sign convention is \(\int\langle V,\Delta_\kappa V\rangle=-\int|D_\Omega V|^2\). The coefficients in (4) are bounded smooth matrices/tensorial maps on the fixed collar. In particular they do not contain an angular cutoff. Derivatives of \(T\) contribute only first- and zeroth-order terms.

## 4. The uniform energy estimate and existence

For homogeneous mixed boundary data define
\[
 \mathcal E[V](t)=\frac12\int_{[0,\infty)\times S^3}
 [|V_t|^2+|V_z|^2+c(t)^2|D_\Omega V|^2+|V|^2],dz\,d\mu_\kappa.
 \tag{5}
\]
Solutions have controlled compact radial support; a smooth extension of coefficients outside it has no effect by finite propagation. On the finite time interval, \(c(t)\) is bounded above and away from zero, so (5) is equivalent to the standard positive \(H^1\times L^2\) norm.

Multiply (4) by \(V_t\), integrate by parts in \(z\) and on the closed sphere, and differentiate \(c(t)^2\). The sole boundary flux is
\(-\langle V_t,V_z\rangle|_{z=0}\). It vanishes: \(V_{D,t}=0\) and \(V_{N,z}=0\), with the D/N subbundles orthogonal. Coefficient first-order terms are bounded by \(C\mathcal E\); the differentiated sphere weight contributes \(c c'\int|D_\Omega V|^2\), also bounded by \(C\mathcal E\). Hence, if \(P_TV=F_T\),
\[
 \mathcal E'(t)\le C\mathcal E(t)+\sqrt{2\mathcal E(t)}\,\|F_T(t)\|_{L^2}.
\]
A harmless regularization at zero energy gives
\[
 \sqrt{\mathcal E(t)}\le e^{CT}
 \left[\sqrt{\mathcal E(0)}+C\int_0^t\|F_T(s)\|_{L^2}\,ds\right].
 \tag{6}
\]
The factor \(C\) can be enlarged to absorb the equivalence constants of \(T\). Thus the original \(U\) is controlled in the same positive Sobolev norm. Every constant depends on a fixed finite number of derivatives of the background, bundle connection and \(T\), and not on harmonic degree.

For a cutoff statement use orthogonal **full-bundle SO(4) isotypic projections** onto an exhausting sequence of finite sets of complete irreducible types, not componentwise scalar cutoffs. Each type has finite multiplicity in these finite-rank homogeneous bundles over SO(4)/SO(3). The background operator, D/N split, \(M\) and \(T\) are equivariant and commute with these projections. Apply (6) to projected data and forcing: projection norm is at most one and the same constants apply at every cutoff. The previous finite-dimensional radial construction supplies each smooth projected solution. Its Volterra contraction uses bounded smooth coefficients and data; its regularity step uses the compatible corner recursion. Thus the construction extends directly from the previously sufficient vanishing-corner case to general fully compatible smooth projected sources. The uniform estimate for differences gives the energy-space limit for approximated data; higher estimates below give the smooth limit for smooth data. The limit satisfies the full system and boundary conditions. Equation (6) proves uniqueness. This establishes a genuine all-angular linear statement, rather than reinterpreting a finite-mode constant as uniform.

For higher regularity, commute tangential time derivatives and angular derivatives formed with the fixed split-preserving product connection (or the SO(4) action). They preserve homogeneous Dirichlet/Neumann conditions. Curvature and coefficient commutators have controlled order. Time differentiation of \(c^2\Delta_\kappa\) is controlled by the collection of commuted energies with one additional angular derivative and one fewer time derivative. Do not commute the physical spacetime connection and assume it preserves the boundary split.

The equation recovers second and higher normal derivatives from time/angular derivatives and forcing. These steps give, at each integer \(m\ge1\), a bound for \(\sum_{j=0}^{m}\|\partial_t^jV(t)\|_{H^{m-j}}\), with \(V_0\in H^m,V_1\in H^{m-1}\), compatible jets, and sufficient forcing norm
\[
 \sum_{j=0}^{m-1}\|\partial_t^jF_T\|_{L^1([0,T];H^{m-1-j})}
 +\sum_{j=0}^{m-2}\|\partial_t^jF_T(0)\|_{H^{m-2-j}}.
 \tag{7}
\]
The last sum is empty for \(m=1\). For \(0\le j\le m-2\), the fundamental theorem of calculus bounds \(\sup_t\|\partial_t^jF_T(t)\|_{H^{m-2-j}}\) by its displayed initial norm plus \(\|\partial_t^{j+1}F_T\|_{L^1_tH^{m-2-j}}\). These are the pointwise forcing norms used in normal recovery. For smooth data all these quantities are finite. This is a sufficient regularity class, not a sharp minimal trace theorem. At each fixed \(m\) the constants remain angular-cutoff independent.

## 5. Boundary data and corners

Equation (6) does not bound arbitrary Neumann forcing by its \(L^2\) boundary norm: \(V_t|_{z=0}\) has no such energy trace estimate. Glancing Neumann frequencies must not be dismissed or used as an automatic ill-posedness claim.

A concrete sufficient smooth lifting is
\[
 W_D(t,z,\Omega)=\eta(z)d(t,\Omega),\qquad
 W_N(t,z,\Omega)=-z\eta(z)g(t,\Omega),
 \quad\eta(0)=1,
 \tag{8}
\]
with \(\eta\) a fixed compact collar cutoff. Set \(R=V-W\). Then \(R_D=0\), \(R_{N,z}=0\), and \(P_TR=F_T-P_TW\). Apply (6) or (7) to \(R\), and add the explicit norms of \(W\). This gives an honest estimate for smooth boundary data, with its stated lifting regularity cost and constants independent of angular cutoff. It is not a lossless \(L^2\)-boundary-control theorem.

Compute initial time jets \(J_j=\partial_t^jV(0)\) recursively from (4). Smooth compatibility is exactly
\[
 (J_j)_D|_b=\partial_t^jd(0),\qquad
 -\partial_z(J_j)_N|_b=\partial_t^jg(0)
 \quad(j\ge0).
 \tag{9}
\]
For a finite Sobolev theorem require every identity in (9) whose traces are defined at that regularity; imposing the smooth identities is a conservative sufficient hypothesis here. The previously selected compactly supported initial data and Taylor sources vanish in an initial shell collar and meet all smooth corner identities after the stated gauge completion.

## 6. Constraints and what has advanced

The reduced control theorem permits arbitrary admissible forcing; such forcing does not automatically satisfy physical constraints. The homogeneous subsidiary statement in this section applies to homogeneous physical linearized equations, or to the actual differentiated Taylor sources satisfying the linearized/differentiated Noether and Codazzi identities.

For the actual model
\(\nabla^AE_{AB}-\tfrac12\nabla_B\operatorname{tr}E=-S\phi_B\). Under this compatible-forcing hypothesis, the reduced equations with \(S=0\) give the homogeneous subsidiary wave equation \((\Box_0+\operatorname{Ric}_0)C=0\) at the relevant coefficient order. Codazzi together with \(K_{ab}=\sigma h_{ab}/6\) and \(n\phi=-\sigma'/2\) gives \(E_{na}=0\), even for field-dependent tension. At the background its boundary conditions are
\[
 C_y=0,\qquad\partial_yC_a-2a_bC_a=0.
\]
Apply the same collar conjugation and uniform control to this rank-five vector bundle. Zero compatible subsidiary initial data imply \(C=0\). Initial physical constraints and the initial harmonic gauge completion provide these data at each already established linear/Taylor order. No extra primal tangential harmonic conditions are imposed.

The previous finite-angular proof is thereby supplemented by an all-angular **linear** estimate at the registered background, including the full twelve-row boundary matrix. This closes the specific frozen-linear control target. The companion nonlinear boundary-straightening note supplies a local field-space version of (3), as an algebraic reduction rather than an evolution theorem. An exact localized initial-data family is already available for the selected reflection-invariant seed at any one stipulated sufficiently high finite Sobolev class. The remaining construction must combine that existing input with a compatible quasilinear mixed-boundary existence proof, an initial harmonic-gauge completion, an explicit parameter-regularity budget and a lifespan extending through the prescribed contact time before asserting an exact nonlinear contact family. No Big Bang origin, external novelty, global stability or submission readiness follows from this note.

The independent boundary referee note supplies an explicit principal glancing calculation: four harmonic rows leave an eight-dimensional constrained Neumann amplitude space, including the scalar mode. This agrees with homogeneous energy control and does not establish lossless boundary-input control.

## 7. Primary literature and attribution

This argument is a model-specific implementation of familiar wave boundary methods, not a certified novelty claim. Fournodavlos–Smulevici, [arXiv:2104.08851v1](https://arxiv.org/abs/2104.08851v1), provide the vacuum umbilic reference-gauge mechanism and a scalar Robin estimate in Appendix A; their four-dimensional constant-umbilicity theorem is not transferred here to a five-dimensional scalar shell. The exponential normal change of variables also appears in the coupled Robin literature, including [arXiv:2202.00760](https://arxiv.org/abs/2202.00760). Its symmetric/symmetrizable matrix hypotheses and low-regularity boundary-input results must be kept separate from the directly proved nonsymmetric full-bundle homogeneous estimate here.

The [targeted literature synthesis](literature/ANGULAR_UNIFORM_CONTACT_LITERATURE.md) records primary linear Neumann existence and variable-coefficient boundary-flux references, the exact statements checked, and bounded search coverage. None alone establishes the pending mixed quasilinear Einstein–scalar contact theorem. The [exact nonlinear boundary-straightening companion](NONLINEAR_BOUNDARY_STRAIGHTENING.md) is a proved state-space reduction with an explicitly unproved evolution program. The [independent boundary referee](referee/INDEPENDENT_BOUNDARY_REFEREE.md) includes the constrained glancing analysis and source-norm qualification.
