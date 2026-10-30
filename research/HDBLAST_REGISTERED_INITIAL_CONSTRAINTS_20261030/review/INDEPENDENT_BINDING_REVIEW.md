# Independent review of the actual initial-data binding

30 October 2026. This is an internal mathematical and implementation review, not external specialist peer review. Frozen sources were read without modification.

## Result and contribution

The physical prefix and ADM assembly binding checkers examined here correctly separate the implementation side from an independent producer of the expected physical and ADM formulas. No mathematical mismatch was found. This new join concerns the current actual registered functions, not a new discovery of the Hamiltonian/momentum right inverse: the retained October 15 proof and its older `review/check_reduced_sources.py` and `review/check_gauge_adm_identities.py` already give all-point source reductions and ADM identities. Finite nonspecial point tests are supplementary.

The independent audit `check_independent_profile_reductions.py` imports neither producer nor project modules. It uses a separately written sparse Laurent polynomial engine over exact rational coefficients, with free actual profiles and their stated differential relations. All 22 collected residuals vanish identically. Eight deliberately wrong formulas are rejected in fresh subprocesses for their named mathematical residual, with no output report. Normal and optimized reports are byte-identical. Input SHA256 checks run before and after use.

The frozen independent checker SHA256 is `dc6284055446a46c7d36b56a56ae2f2859c4f7c4b860ea61e1729a334c77eea8`. Both independent report SHA256 values are `868064712b2458b8fbb4d13b5c7729e1f7d578e214604168f0ae08dfb340ab56`.

## Independent physical derivation and checks

For the first spatial tangent `gamma1_yy=-4uY`, `gamma1_ij=2uY gamma0_ij`, the sphere conformal-curvature contribution is `+20u/rho^2`, while the angular radial-lapse Laplacian contributes `-32u/rho^2`. Keeping both terms gives the first physical Hamiltonian coefficient

```
-6 u_yy -36 a u_y -(24 a_y+48 a^2+4 p^2+12/rho^2)u
-2p v_y -2 Vprime v.
```

It vanishes after `p v=-3(u_y+2au)`, `p_y=Vprime-4ap`, and `a_y=-rho^-2-p^2/3`. No additional potential-valued background Hamiltonian relation is required for this cancellation. The radial momentum coefficient is `[-3(u_yt+2a u_t)-p v_t]/rho`; it also vanishes. The first angular momentum vanishes structurally because the first mixed physical tensor `K^J_I-delta^J_I tr K` has only its radial diagonal entry nonzero on the unperturbed spherical background; it must not be treated as an unproved literal zero.

The independent exact calculations also verify the unsimplified physical-Cauchy-affine Hamiltonian source against registered `aa,dd`, `Fmom=3 K_registered`, and `Gmom-D_y Fmom=p Pi_registered`. They verify the three scalar Hamiltonian inverses including the monopole denominators `-6,10,42`, plus the generic radial and angular momentum inverse. Their scalar compensator is an ordinary unit-normal velocity, and `K2` is covariant physical extrinsic curvature. They are not canonical momentum densities. Omitting the physical radial cross term `-6u_y u_t/rho` changes the affine convention and is rejected.

## Assembly, harmonic values, and normal derivative

The independently expanded exact ADM prescription gives the registered coordinate lapse/shift velocities. The metric conversion must retain both `g2_tt=-2rho^2 U^2` and `dot g2_ty=dot beta2^y-8U dot beta1^y`. The current suffix checker compares all 23 displacements and 23 scaled velocities to an independently assembled frame table. The five explicit mean overrides, scalar field, zero initial Hessian trace-free/TT entries, and every-slot `ell` velocity factor are included. The authenticated adapter provides a third join, rather than the independent expected table.

The generic ADM identity enforces `C=0` for arbitrary angular functions in the exact initial lapse/shift prescription. Consequently its scalar degrees 0,2,4 for temporal/radial components and gradient degrees 2,4 for angular components vanish. Those eight generic modal consequences are distinguished from the additional direct source-map calculation now supplied. That direct calculation uses the corrected first harmonic one-jet, including `B_t=-6rho(u_y+4au)`, and verifies the first linear harmonic covector vanishes. The authenticated quadratic map gives

```
QC_y=-24u u_y Y^2,
QC_0=24u u_t Y^2/rho,
QC_a=24u^2 Y D_aY/rho.
```

I independently collected these expressions directly from `connection_contractions`: its three nonzero profiles are respectively `12u u_y`, `-12u u_t/rho`, and `-12u^2/rho` before the map's `-2` factor. The independently evaluated linear `DC` on the second registered one-jet gives their negatives in all eight modal coefficients. The physical first harmonic row makes the full-metric-lowering cross term zero. No `ell` factor enters this ordinary physical-frame covector calculation. The added controls detect the longitudinal first one-jet and a reversed radial quadratic source sign.

For the future timelike initial normal and `C=0` on the whole initial slice, let `X_A=(nabla_n C)_A`. The reduced equation gives `E_nn=X_n`, `E_nI=X_I/2`, `E_IJ=0`, hence `tr E=-X_n` and the physical trace reversal gives `mathscr H_nn=X_n/2`, `mathscr H_nI=X_I/2`. Thus the physical Hamiltonian is `2 mathscr H_nn=X_n`. Using `E_nn` alone doubles that Hamiltonian factor. This is an implication under the actual reduced acceleration equations; it is not a separately prescribed acceleration or a new computed `QC_x` enclosure. Future initial normal and outward shell normal must remain distinct.

## Exact scope and retained premises

The physical prefix test and suffix assembly test compose because the latter allows free physical `u,u_y,u_t,v_t,f_l,f_l_y,K,Pi`, while the former binds those quantities to the actual packet/background. The final physical checker explicitly places all three `f_l`, all three `f_l_y`, and the scalar compensators in the packet-jet ideal, so its result composes with the suffix ideal check.

A source binding can remove the need merely to assume that the actual registered initial functions equal the retained analytic physical/ADM formulas. It does not prove a global subsidiary identity by initial rows alone. In particular, the actual first harmonic wave must have the same initial spacetime one-jet used by the ADM prescription, including the nonzero first mixed time velocity. The differentiated Bianchi relation, actual shell Codazzi/source identities, subsidiary boundary conversion, regularity and uniqueness remain analytic premises unless separately audited. Coordinate source intervals that merely contain zero do not establish those identities.

The finite packet-ideal checks support the open initial shell collar statement. Smooth flat support and closure under finite differentiation give the all-orders structural statement only together with the retained homogeneous corner recursion; a finite list of jets is not itself an all-orders corner proof.

Nothing in this audit proves a finite-amplitude nonlinear family, physical relative accuracy of the numerical candidate, an empirical higher-dimensional origin of the Big Bang, novelty against all literature, journal acceptance, or external peer review.

## Final source review

The updated physical checker binds the actual polynomial potential and all four differentiated coefficient arrays to the authenticated source AST and actual-model JSON, replaces a literal first angular zero by the physical trace identity, and reauthenticates its independent arithmetic engine after use. The ADM checker supplies both the generic exact initial gauge identity and the separate direct eight-row initial quadratic-source join just described. No mathematical mismatch was found in either implementation. Their reported check counts include component joins and consequences; they count verification bookkeeping rather than independent scientific theorems.

The final public scientific README was reviewed. It correctly describes an implementation-to-formula join, the direct eight initial quadratic-source rows, and retained actual first-wave initial matching and global propagation inputs. It explicitly preserves the candidate/error ledger and excludes physical smallness, finite-amplitude existence and a cosmological origin claim. No scope correction is required.

Reviewed publication files and SHA256 values:

| File | SHA256 |
| --- | --- |
| `checks/check_independent_profile_reductions.py` | `dc6284055446a46c7d36b56a56ae2f2859c4f7c4b860ea61e1729a334c77eea8` |
| `checks/check_physical_initial_binding.py` | `3df453e4dd53378977815954cdbc111e87ab12ad2b248269d84323f67814b5d4` |
| `checks/laurent.py` | `02d2c85ec4188149f4a3c53c62a7781d17cf0a2f8e4b53cf3da5b697689e0b5b` |
| `checks/check_registered_adm_binding.py` | `6b5ec8e23a6dea34273b9ddac1fdd85d3a0461ddf66bea810f3042096dd83ee5` |
| `reference/INDEPENDENT.json` | `868064712b2458b8fbb4d13b5c7729e1f7d578e214604168f0ae08dfb340ab56` |
| `reference/PHYSICAL.json` | `48e546a90b09b9b4e4466575ea1b74b8b3082be3cc82a51e5956b7b4ef93d132` |
| `reference/ADM.json` | `de60e029ff4afe24619d0a6a9036064d867efe3aed0380a288c118eac8ddd8f6` |
| `README.md` | `1b75f5ab2ce3470bba639b7b347053d2d6103bd79c3cbf86abcc64a8aff23cc7` |
