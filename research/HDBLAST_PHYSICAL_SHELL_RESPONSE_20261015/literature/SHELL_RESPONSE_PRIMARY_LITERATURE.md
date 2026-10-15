# Primary literature for the physical-shell response target

Task date: 2026-10-15. This is a bounded search and source comparison, not an exhaustive review, a novelty certification, or independent mathematical peer review.

The directly useful result is the moving-hypersurface pullback identity: the linear physical shell scalar includes the embedding displacement. The retrieved literature does **not** supply a theorem that a nonzero compact bulk perturbation must produce a nonzero HDBLAST shell scalar. The source conditions make a direct calculation for a selected profile the strongest next route.

## 1. Define the physical response before observing it

Let \(j_\varepsilon:\Sigma\to M_\varepsilon\) describe the physical shell, and set
\[
 \psi_\varepsilon=j_\varepsilon^*\phi_\varepsilon .
\]
In an identification with the background spacetime, with embedding velocity \(Z_1\),
\[
 \Phi_1=\partial_\varepsilon\psi_\varepsilon|_0
       =j_0^*(\phi_1+\mathcal L_{Z_1}\phi_0).
\]
This is the scalar specialization of Marc Mars, *First and Second Order Perturbations of Hypersurfaces*, [arXiv:gr-qc/0507005](https://arxiv.org/abs/gr-qc/0507005), Lemma 1, equation (7). That lemma assumes a \(C^2\) family of fields and embeddings; it is a geometric differentiation identity.

Mars's Proposition 5, equation (43), gives \(Z_1\mapsto Z_1-S_1\) under a spacetime identification change. The scalar tangent gains \(\mathcal L_{S_1}\phi_0\), so the two terms cancel. A coordinate value \(\phi_1|_{y=y_b}\) represents the physical response only after the shell-fixing gauge has been specified.

Shinji Mukohyama, *Perturbation of junction condition and doubly gauge-invariant variables*, [arXiv:hep-th/0006146](https://arxiv.org/abs/hep-th/0006146), Sections II–III, independently explains why spacetime gauge transformations and intrinsic shell reparameterizations are independent. In HDBLAST the background shell scalar is constant, \(\psi_0=\phi_0(y_b)\). Consequently \(\Phi_1\) also has no first-order intrinsic reparameterization term: \(\mathcal L_X\psi_0=0\). This is the relevant simplification; background constancy of the **bulk** scalar is not required.

For a sufficiently regular family and a specified finite shell observation region, consider
\[
 Q(\varepsilon)=\int_\Sigma w_\varepsilon
       (\psi_\varepsilon-\psi_0)^2\,d\mu_\varepsilon,\qquad w_0\geq0.
\]
Ordinary differentiation at the vanishing background difference gives
\[
 Q'(0)=0,\qquad Q''(0)=2\int_\Sigma w_0\Phi_1^2\,d\mu_0 .
\]
Weight and measure derivatives disappear at this order. A geometrically specified observation domain, weight, and clock are needed if \(Q(\varepsilon)\) itself is to be advertised as a fully nonlinear invariant. The displayed Hessian is already invariant under first-order shell reparameterization when \(\psi_0\) is constant. This standard calculus identity neither supplies \(\Phi_1\neq0\) nor determines a detector signal.

## 2. Do not replace the bulk-scalar system with a vacuum master equation

Hideo Kodama, Akihiro Ishibashi and Osamu Seto, *Brane World Cosmology — Gauge-Invariant Formalism for Perturbation*, [arXiv:hep-th/0004160](https://arxiv.org/abs/hep-th/0004160), provide a general gauge-invariant framework allowing bulk matter and shell junction conditions. Their scalar reduction in Section IV.2.3, equations (151)–(153), uses a vacuum background with a constant-curvature two-dimensional orbit space. Their discussion states that the scalar master potential was introduced only under this orbit-space condition.

This framework motivates deriving gauge-invariant bulk-plus-shell variables. Its displayed vacuum master equation cannot simply be imported into the present nonconstant bulk Einstein–scalar background. Section III.3.3 also explains that the junction boundary condition depends on the intrinsic matter model. HDBLAST's actual potential-dependent scalar law and Israel law must remain in the calculation.

## 3. A nonzero bulk state need not be visible in one boundary component

Tatsien Li and Bopeng Rao, *Approximate boundary synchronization by groups for a coupled system of wave equations with coupled Robin boundary conditions*, [arXiv:2202.00746](https://arxiv.org/abs/2202.00746), study
\[
 U_{tt}-\Delta U+AU=0,\qquad
 \partial_\nu U+BU=0
\]
with an observation of selected linear combinations \(D^TU\), constant finite coupling matrices, and separated Dirichlet/Robin portions of a smooth bounded Euclidean boundary. Their applicable well-posedness assumptions include a symmetric, or suitably symmetrizable, boundary matrix.

Lemma 2.1 identifies the largest subspace invisible to \(D^T\) and invariant under both coupling matrices. Proposition 4.2 makes full composite rank a **necessary** condition for incomplete boundary observation. The introduction and Section 5 explicitly state that this rank condition is generally not a sufficient continuation theorem. It does not turn an incomplete boundary trace into full Cauchy data. Full observation and special two-equation structures admit stronger statements, with further conditions and adequate observation time.

This is a caution, not a counterexample specific to HDBLAST. HDBLAST has different time-dependent warped coefficients, constraints, derivative couplings, and a mixed metric/scalar boundary operator. Neither the nonzero bulk norm nor a finite matrix rank test establishes \(\Phi_1\neq0\).

Li–Rao's separate [2022 internal-observation article, DOI 10.5802/crmath.341](https://doi.org/10.5802/crmath.341), has an abstract claiming necessary and sufficient Kalman rank under its own setup. Only that publisher abstract was retrieved here. It concerns **internal** observation and must not be cited as sufficiency for a Robin boundary observable.

Eva Sincich and Sergio Vessella, *Wave equation with Robin condition, quantitative estimates of strong unique continuation at the boundary*, [arXiv:1603.08778](https://arxiv.org/abs/1603.08778), Theorem 2.1, is another useful model. Its scalar equation has time-independent, symmetric, uniformly elliptic Lipschitz spatial coefficients \(A(x)\), a Lipschitz Robin coefficient, and a local \(C^{1,1}\) boundary. The theorem quantitatively propagates smallness/flatness of the **boundary trace**. It is not a ready-made finite-time observability theorem for HDBLAST's coupled, time-dependent Einstein–scalar system.

## 4. Recent related papers: useful comparisons with explicit limits

Kosei Fujiki, Hiroki Kanda, Michitaka Kohara and Tadashi Takayanagi, *Brane Cosmology from AdS/BCFT*, [arXiv:2501.05036](https://arxiv.org/abs/2501.05036), [JHEP03(2025)135](https://doi.org/10.1007/JHEP03(2025)135), study EOW-brane dynamics, mainly in \(\mathrm{AdS}_3\), with a scalar **localized on the brane** and pure gravity in the bulk. Sections 5 and 7 describe big-bang-like toy cosmologies for suitable potentials. The scalar location, dimensions, action and boundary system differ from HDBLAST. These constructions are related mathematical models, not evidence that a higher-dimensional blast created our Universe or that our shell response is nonzero.

Ruikson S. O. Nunes, Miguel R. Nuñez-Chávez and Marcelo M. Cavalcanti, *Exact boundary controllability for one-dimensional wave equations systems with a time dependent potential in non-cylindrical domains*, [DOI 10.3934/eect.2025080](https://doi.org/10.3934/eect.2025080), published online 4 November 2025 and listed in the 2026 volume, provide a newer moving-boundary comparison. The retrieved publisher abstract/introduction restrict the problem to one space dimension, an endpoint moving below characteristic speed, a time-only smooth bounded potential that vanishes after a suitable time, and Dirichlet control at the moving endpoint. The later theorem body was not fully retrieved or audited. These conditions do not directly cover HDBLAST's fixed-shell coupled Robin problem.

No directly applicable newer theorem or observational discovery was identified by this bounded search. That statement does not establish that none exists.

## 5. Recommended calculation and certificate

Retain the selected aligned quadrupole and its exact linearized Einstein–scalar initial data. Derive the actual coupled radial-time evolution and shell operator, with the embedding term included or a documented shell-fixing gauge. For one explicit annular profile, compute a witness for the weighted trace norm after an admissible interaction time.

A reproducible numerical approximation becomes a nonvanishing proof only with an error bound in the same observation norm. If \(\widetilde\Phi_1\) is the computed trace and
\[
 \|\Phi_1-\widetilde\Phi_1\|_{L^2(w_0d\mu_0)}\leq\delta
 <\|\widetilde\Phi_1\|_{L^2(w_0d\mu_0)},
\]
the reverse triangle inequality gives a strictly positive response norm and hence \(Q''(0)>0\). Establishing the bound must address the coupled bulk residual, initial constraints, actual shell residual, and observation-trace stability. A converged plot or several grid resolutions alone do not provide that certificate.

The retained \(C^2\) nonlinear family would then transfer the positive first-response coefficient to sufficiently small nonzero amplitudes. Second-order choices of the initial completion do not enter this Hessian once the first physical tangent is fixed. Prior to such a witness or a fully checked analytic response theorem, the correct status is a well-defined response target with nonvanishing still open.

## Provenance

The source ledger records eight primary-source comparisons, access level, equation/section locations, and explicit applicability limits. The bounded search made five web searches, four paper searches, seven targeted paper reads, seven full or partial URL retrievals, and four canonical metadata inspections. One metadata search's returned identifiers were retained in the receipt while its full original response was not separately saved. An unsuccessful paper-index full-text lookup was retained and labeled.

Raw retrieved bodies in \`sources/\` are internal research records. Publish the original synthesis and source ledger if useful; do not automatically redistribute full publisher or arXiv bodies. Experimental HTML included duplicated mathematical notation and noncanonical date headings in older papers; canonical arXiv identifiers/metadata determine publication dates.

