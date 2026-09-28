# Independent review of the proposed matter extension

22 September 2026. Reviewed `MATTER_EXTENSION_AND_RESIDUAL_VACUUM.md` and `verify_matter_extension.py`. This is a separate derivation and interpretation check, not another execution of the 47 symbolic assertions. No producer file was edited.

**Verdict:** the action conventions, junction coefficients, energy-transfer signs, projected Einstein equations, Weyl balance and conditional radiation thresholds are internally consistent. They specify a useful extension but do not close a four-dimensional evolution or establish particle production. One production-energy interpretation requires special care, explained below.

## Action and signs

For two bulk copies, one shell action and an outward normal from each interior, the scalar boundary variation is

\[
-\frac{2n\phi}{\kappa_5^2}-\frac{\sigma'}{\kappa_5^2}-j=0.
\]

The reported scalar junction therefore has the correct minus sign and factor two. The Israel equation likewise gives `k_s=(σ+κ₅²ρ)/6` and `k_0=[σ−κ₅²(2ρ+3p)]/6`. As a separate sign check, positive pure tension on the chosen interior has positive k_s, consistent with the original registered boundary orientation.

The matter Ward identity has `D_μT^μ_ν=−jD_νφ`. Its time component is `ρ̇+3H(ρ+p)=jv`, since `D_μT^μ_τ=−[ρ̇+3H(ρ+p)]`. This sign agrees with an increasing χ mass doing positive work on its energy.

Contracted Codazzi gives `D_μS^μ_ν=2T_bulk,nν` for this orientation. At ν=τ this independently implies

\[
\dot\lambda+\dot\rho+3H(\rho+p)
=-2T_{n\tau}^{\rm bulk}
=-2(n\phi)v/\kappa_5^2.
\]

Thus the action, matter ledger and bulk flux all have consistent signs. Changing the normal requires changing all these conventions together.

## Independent projection and Weyl check

With `X=(∇φ)²=−v²+w²`, direct projection of the bulk scalar stress gives

\[
F_{\tau\tau}=v^2/4-w^2/4+U/2,
\quad F_{ij}=(5v^2/12+w^2/4-U/2)h_{ij}.
\]

The extrinsic contributions are `3k_s²` and `−k_s²−2k_sk_0`. Using the report's `𝒲=−E_ττ/3`, tracelessness gives a spatial Weyl contribution `+𝒲`. These expressions reproduce both reported H² and Ḣ, including the negative scalar-normal-gradient term in H² and the `−2𝒲` term in Ḣ.

Differentiating `k_s` and using exchange independently gives

\[
\dot k_s=-wv/3-\kappa_5^2H(\rho+p)/2.
\]

Substitution into the differentiated Friedmann identity reproduces the report's Weyl balance. It requires v̇ and ẇ, so it is a compatibility identity, not a missing evolution equation that closes the system. In particular, setting 𝒲=0 during production is not warranted by the local junctions. The report correctly states this limitation.

The field and coupling dimensions also agree: Φ has mass dimension 3/2, g_b has −1/2, ḡ has 1, and q has 2. The registered pure-tension calculation supplies none of the missing mass scales or matter couplings.

## Conditional radiation tests

Under fixed scalar/tension, negligible Weyl contribution and no continuing exchange, the radiation Friedmann contribution is

\[
F_r=\kappa_5^2\sigma_f\rho_r/18+\kappa_5^4\rho_r^2/36.
\]

Solving `F_r=H_vac²` gives the stated dominance threshold. The acceleration instead contains `−κ₅²σ_fρ_r/18−κ₅⁴ρ_r²/12`, which gives the distinct stated deceleration threshold. With freely redshifting radiation, `ρ_r∝a⁻⁴` yields the N_max formula. These are correct within the benchmark.

They are not automatically valid throughout the original five-dimensional transition: the bulk scalar, its normal derivative, the matter source and Weyl term can evolve. The full report explicitly retains this caveat. The 0.135-e-fold toy budget additionally assumes a fixed four-dimensional Planck mass; its own tension comparison shows why that assumption cannot simply be applied to the transition. It must not be presented as a five-dimensional no-go bound.

## Particle-production interpretation to preserve

The Gaussian number and energy integrals in the verifier are mathematically correct. However, the familiar linear-crossing occupation `n_k=exp[−π(k²+m₀²)/q]` is an **after-crossing** particle occupation in the idealized scattering problem. At the exactly massless crossing, the evolution is nonadiabatic and an instantaneous particle decomposition is not unique.

Consequently, weighting that occupation by massless energy k gives `q²/(4π⁴)` as a **formal energy proxy**, not an independently demonstrated instantaneous renormalized energy density. The report calls it an estimate and makes the associated prompt-transfer budget conditional; that distinction must remain explicit in any summary. A physical transfer acting during the crossing can also alter the production calculation itself. Quantitative energy deposition needs mode evolution, a consistent stress/source prescription, and the signed exchange ledger.

The crossing-eligibility scan is correspondingly a local approximation screen of archived matter-free trajectories. It does not choose a viable coupling, verify its cutoff, control backreaction, or establish a thermal bath. No action-level error was found that invalidates that limited use.

## Scope of the exact verifier

The verifier correctly includes wrong-sign/factor controls and exact scalar-junction and expansion checks. Several tests are identities within assumed projection formulas; they are not numerical tests of a coupled five-dimensional solution. Independent derivation above supports those formulas. The earlier independent static-branch review supplies a separate check of the proposed scalar-profile candidate. Neither verification establishes that the original unstable shell reaches that candidate or successfully reheats.
