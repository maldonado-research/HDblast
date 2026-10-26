# Independent bounded review of the new angular structure audit

26 October 2026 UTC. Internal AI-assisted mathematical/software review by a separate agent. This review read `exact_angular_structure.py`, `ANGULAR_STRUCTURE_PROOF.md`, the normal/optimized reports and controls, and the frozen definitions already identified in the subsidiary evidence note. It ran only the **new** structure checker against the frozen repository, first before and then after the exact active embedding/restriction was added, each with a fresh output in `/tmp`. It did not run any frozen producer or re-prove the upstream PDE/source certificates.

**Conclusion:** no blocking mathematical or implementation defect found within the stated scope. The exact selected 22-component reduction follows from equality of the two complete TT equations and all their data. It is not inferred merely from a symmetry name or a vanishing initial value.

## Checks and reasoning

1. **Block inference.** The checker tests the support of every primitive matrix, including the derivative, scalar, curvature and boundary generators. Since assembled matrices are scalar linear combinations/products already represented by these primitives, the resulting `5+8+8+1+1` invariant coordinate blocks hold for the formal coefficient family. Strong connectivity of their union establishes only coordinate indecomposability; the proof correctly disclaims absolute irreducibility under general simultaneous changes of basis.

2. **Angular parity.** The reflection is `R=diag(1,1,1,-1)` on ambient sphere coordinates. Tensor pullback multiplies both monomials and every reflected covector slot by the appropriate sign. Omitting the slot signs would be wrong; the checker retains them. An orientation-reversing orthogonal transformation exchanges self-dual and anti-self-dual Weyl parts. Here it fixes the seed `Y=n0^2-n1^2`, so the explicit exchange `R*Tplus=Tminus` has the correct plus sign. The selected scalar, gradient and trace-free Hessian templates are invariant. This is angular orientation reflection, not antipodal parity and not shell-normal reflection.

3. **Data and evolution equality.** Primitive TT off-block entries vanish, their diagonal coefficients agree, and their full linear Robin rows vanish. All eight stored geometric tensor maps have equal TT term lists, including both simplified and general reduced maps. The generic shell calculation yields equal TT forcing, and the chosen initial adapter gives zero TT displacement and velocity. Thus their difference solves a homogeneous scalar wave problem with zero initial data and homogeneous Neumann boundary. Inherited uniqueness, explicitly retained as a conditional premise, forces that difference to vanish. The argument remains valid even if the common TT source is nonzero. It would not apply to added independent homogeneous second data with unequal TT values.

4. **TT equation and normalization.** The displayed potential agrees with the primitive TT entries `Sa=-2`, `Sh=2`, `Rb=1`, `Rd=4`, `Delta=-6`, `Pmet=1`, with all TT connection/coupling entries absent. The shell factor follows from `QKab_T=-(4/3)(a*h_s+sigma'*v_s/3)d_s` and the registered forcing factor `2*ell*rho`, giving `-(8/3)*ell*rho*(a*h_s+sigma'*v_s/3)*d_s`. The common coefficient `tau=tau_plus=tau_minus` multiplies `Tplus+Tminus`; orthogonality and individual norm squared `1/6` give Gram weight `1/3`. Its source is each original TT source, not half that source. The proof explicitly preserves this distinction.

5. **Shell-normal parity.** The doubled copy-even metric has odd normal-tangential components and even tangential/normal-normal components; TT modes are even in that sense. Evenness does not imply a zero one-sided normal derivative when the doubled shell admits a cusp. The proof correctly describes the TT primal row as inhomogeneous Neumann and does not delete its forcing. The five Dirichlet components and 18 ambient Robin components reduce to five and 17 in the selected 22-component representation.

6. **Physical subsidiary scope.** The proof accurately states the inherited conditional formal completion, including `C2=DC[k]+QC`, zero compatible initial subsidiary value/velocity, normal component Dirichlet and induced tangential Robin conditions. It does not append additional primal boundary rows. The preceding independent evidence audit gives the retained primary proof locations. The new result neither numerically verifies all constraints nor gives nonlinear realization. Its requirement for quantitative C2 or a declared weak norm before claiming pointwise Einstein residual control is correct.

7. **Active representation.** The added embedding `E` copies the common active TT coefficient to both ambient TT slots; its restriction `R` averages those two slots and preserves all others. The definitions give `R E=I22`, and the check `M E=E R M E` for all 25 primitive matrices establishes invariance and the reduced operator `R M E`. The inherited diagonal Gram pulls back to the reported 22 diagonal entries, with final entry `1/3`. Averaging equal source slots returns each original source without an extra factor.

## Verification

The new checker succeeded in isolated Python (`-I -S -B`); its independent replay report is byte-identical to `ANGULAR_STRUCTURE_NORMAL.json` and `ANGULAR_STRUCTURE_OPTIMIZED.json`, with SHA256:

`3a0ac02eac30388f75ddbf71fb9194bc3a888f67466479c2bc978fa558f8018c`.

The reviewed final checker SHA256 matches the producer pin in `STRUCTURE_CONTROLS.json`:

`7c4b36532a2d4fd174028d8d8957b32d18377231ed5510b08470d3888cc26946`.

The controls report records the four intended rejections and absence of erroneous output artifacts. This review inspected those records; it did not independently rerun the four mutation cases. The full seven frozen closures were authenticated before project imports and after the successful replay.

No further change is required for the stated bounded structure claim. Full numerical angular evolution, a continuous error ledger, and the appropriate numerical subsidiary/residual checks remain subsequent work.
