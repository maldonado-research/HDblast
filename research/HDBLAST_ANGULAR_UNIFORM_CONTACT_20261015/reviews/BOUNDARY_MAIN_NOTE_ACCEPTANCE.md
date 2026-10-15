# Independent review acceptance of the frozen linear and exact boundary-flow notes

15 October 2026. Internal AI-assisted mathematical review; not external peer review. Final complete read-back finished at 2026-10-15 00:25:04 UTC; both file hashes were then checked against the requested frozen snapshots.

Reviewed file: `../ANGULAR_UNIFORM_LINEAR_CONTACT.md`.

Accepted SHA-256: `cc559ed4a16c920310fca4c97c3af869508414b6d0ffc9fd8b4612e82989dff1`.

Reviewed companion: astra/NONLINEAR_BOUNDARY_STRAIGHTENING.md, to be copied byte-for-byte as NONLINEAR_BOUNDARY_STRAIGHTENING.md in the publication capsule.

Accepted companion SHA-256: `da221f388e0397e4a3b06ce40be1feff1d170311f3b7eaf79f3644b667dd6882`.

I accept the mathematical claim of this exact snapshot: the full reduced linearized Einstein–scalar system at the registered smooth background has a positive auxiliary mixed-boundary energy estimate uniform over full tensorial angular cutoffs; compatible smooth data have a unique smooth solution; the stated lifted boundary-data class is sufficient; compatible homogeneous subsidiary data vanish by the same method.

The review checked the following material points:

1. The twelve direct-metric boundary rows, signs under `dy=-rho dz`, and distinction between the normal lapse coefficient and the metric coefficient are consistent with the exact laws.
2. The localized block-diagonal conjugation preserves the four Dirichlet rows and exactly removes the nonsymmetric Robin matrix. The scalar principal symbol leaves only bounded first- and zeroth-order new terms.
3. The positive product fiber metric and split-preserving connection justify the zero homogeneous boundary flux, all angular commutators, and normal recovery.
4. Equation (7) is a sufficient high-regularity source norm: its time-integrated derivatives plus initial values control the pointwise source terms needed for normal recovery by the fundamental theorem of calculus. The same initial source terms control the recursively defined initial wave jets. The smooth corner assumption is conservative and sufficient.
5. Complete SO(4) isotypic projections commute with the fixed-background operator, splitting and transformation, have finite angular multiplicity, and contract the invariant energy norm. The previous finite-angular construction extends to general fully compatible smooth sources by its Volterra construction and compatible corner recursion; the difference estimates then give the all-angular limit. A componentwise scalar cutoff would not justify this argument. The added angular-parallel property of the background Robin endomorphism is correct for the stated product connection.
6. The corrected section 6 explicitly restricts homogeneous subsidiary propagation to physical/Bianchi-compatible forcing. Arbitrary admissible reduced-system forcing does not acquire physical constraint propagation automatically.
7. The glancing qualification is correct: the four harmonic rows leave an eight-dimensional constrained Neumann amplitude space, including the scalar. The note asserts no arbitrary-boundary-`L2` energy theorem, spectral stability, or nonlinear existence.

The companion symbolic checker records five algebra groups with seven explicit checks. Ordinary Python and `python -O` both pass. Under `python -O`, each of four actual corruptions is rejected with exit status 1: wrong scalar Robin coefficient, wrong inverse-metric trace factor, an extra scalar glancing constraint, and an incorrect conjugation derivative. These checks audit displayed algebra; they are not a computer proof of the PDE theorem.

I accept the companion's exact boundary-flow proposition. Its inward-coordinate sign, subtraction of the actual background normal derivative, lapse-to-metric state change, localized pointwise flow, invertible state Jacobian and exact homogeneous D/N reduction are consistent. Common scalar wave principal coefficients survive both pointwise changes of dependent variables. The retained Dirichlet rows force vanishing boundary inverse-metric normal/tangent coefficients. Its energy paragraph correctly requires spatial energy positivity near the background in addition to timelikeness of the time covector. Its five-step quasilinear discussion is explicitly an unproved program, not an existence theorem.

An exact localized initial family already exists at fixed finite regularity under the hypotheses in the earlier exact-precontact note, with metric/momentum in \(H^{k+2}\times H^{k+1}\), \(k\geq12\). This acceptance does not extend that family through contact or assert a compatible quasilinear mixed-boundary existence lemma, parameter regularity without derivative bookkeeping, or a prescribed-contact lifespan. The established finite-angular second-coefficient theorem is unchanged. The main note's literature section distinguishes model-specific implementation from novelty and does not import the four-dimensional vacuum or symmetric Robin theorems as proofs of the pending coupled nonlinear result; primary-source bibliographic verification is separately recorded in the literature agent's artifact.
