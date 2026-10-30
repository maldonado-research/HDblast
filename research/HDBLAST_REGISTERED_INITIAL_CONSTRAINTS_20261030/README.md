# Exact binding of registered initial data to the physical constraints

30 October 2026. Separate continuation for the retained admissible aligned
quadrupole, preserving every earlier manuscript, source capsule and numerical
candidate.

The registered actual initial functions have been compared with an independent
assembly from the existing physical Cauchy and ADM harmonic-gauge derivations.
The physical source reductions and all **23 displacement plus 23 scaled velocity
entries** agree exactly, including the five explicit mean overrides. No
registered mismatch was found.

This closes an implementation-to-formula matching step. The underlying
all-point analytic initial-constraint and gauge proofs already existed; this
checkpoint does not claim to have discovered those proofs or a new nonlinear
existence theorem.

## What is checked

The calculations collect rational Laurent polynomials, using the exact retained
background differential equations and formal packet derivatives. A successful
identity has every collected residual coefficient equal to zero. No selected
point, interval midpoint, nominal background, cutoff or mollified candidate is
substituted for the actual analytic data.

The physical audit checks the authenticated background and potential, first
reconstruction identities, and the original physical-Cauchy-affine Hessians.
In particular it verifies

```text
aa = A_physical,       dd = D_physical,
3 K_registered = F_mom,
p Pi_registered = G_mom - D_y F_mom.
```

The local `K_registered` denotes the coefficient of the physical extrinsic
curvature correction, not the full23 metric slot `K:0`. Covariant extrinsic
curvature and unit-normal scalar velocity are the physical variables. A
spacetime-affine path or canonical momentum density has different nonlinear
conversion terms and cannot silently replace them.

Substitution verifies the three Hamiltonian, three radial-momentum and two
angular-gradient momentum modal rows. The Hamiltonian denominators remain
`-6, 10, 42`, including the monopole. The scalar compensator cancels the radial
derivatives of the metric correction using `p_y + 4 a p = V_prime`.

The separate ADM calculation expands the exact initial lapse/shift gauge
prescription before converting it to spacetime/frame data. It retains
`g2_tt = -2 rho^2 U^2`, the mixed velocity correction
`-8 U beta1dot^y`, and the factor `ell` in **every** velocity entry for
`x=t/ell`. Its independent label table agrees with the authenticated actual
registration suffix, the source adapter and all mean overrides. Initial
trace-free Hessian and TT slots are also checked.

The generic exact ADM identity gives the initial harmonic cancellation.
A separate calculation evaluates the authenticated quadratic `QC` source map
on the corrected first initial one-jet and checks all eight modal rows of
`DC[initial second jet] + QC_initial = 0`. The first harmonic covector is
checked to vanish before using the simplified source adapter. These are
initial-source identities, with ordinary physical-frame derivatives. The
actual first gauge-wave's matching to that prescribed initial one-jet remains
part of the retained analytic chain; no evolving `QC_x` estimate is inferred.

## Initial subsidiary derivative and shell corner

The universal trace-reversal calculation is also checked. For future initial
unit normal `n`, with `C=0` on the slice and reduced equations
`E_ab = nabla_(a C_b)`, put `X_a=(nabla_n C)_a`. Then

```text
E_nn = X_n, E_nI = X_I/2, E_IJ = 0, tr E = -X_n,
(E - g tr E/2)_nn = X_n/2,
(E - g tr E/2)_nI = X_I/2.
```

Thus the physical initial constraints imply the initial normal subsidiary
derivative vanishes, under the retained reduced-acceleration equations.
Using trace-reversed `E_nn` alone as the Hamiltonian constraint loses the
required factor. This is an exact conditional implication, not a new direct
`QC_x` enclosure or a prescription of extra accelerations.

The actual registered one-jet belongs to the packet-jet vanishing ideal.
Finite differentiation preserves this differential ideal. Smooth flat packet
support on an open initial shell collar and the retained homogeneous reduced
recursion supply the higher corner-order argument; a finite table by itself
is not an all-orders corner proof.

## Reproduce

Use a complete repository checkout and a fresh output directory outside this
capsule:

```sh
python -I -S -B research/HDBLAST_REGISTERED_INITIAL_CONSTRAINTS_20261030/checks/run_checkpoint.py \
  --repo-root . --output-dir /tmp/hdblast-initial-binding-fresh
```

Only the Python standard library is required. CI uses Python 3.12.14. The runner
authenticates the closed capsule, reproduces the complete two historical
physical/gauge rosters, and executes three new reports in both normal and
optimized modes. Each complete report must match its reference bytes.
Mathematical mutation controls run as separate processes and must fail at the
intended identity without producing a successful report. Source hashes are
checked before use and after readback.

The independent abstract-profile audit uses its own arithmetic implementation
and imports no producer. It checks the source reductions separately; the
physical-prefix and ADM-suffix audits provide the current registration join.
The review notes distinguish these roles and do not call repeated coefficient
rows independent discoveries.

## Limits and next priority

The result applies to the retained exact background and first reconstruction,
the stated physical constraint/ADM formulas and the fixed angular product
identities. Nonzero `rho` and `p` justify the Laurent denominators. It binds
the actual initial functions used in the later enclosures, without changing
`W1_quadratic`, its source, or any previously published error bound.

The initial matching calculation does not independently reprove the full
first harmonic identity on a neighbourhood, differentiated scalar/Bianchi
identity, physical-shell Codazzi source identity, subsidiary boundary
conversion, boundary regularity or uniqueness. Their existing analytic proofs
remain identified inputs. The next priority is a similarly explicit binding
of the registered bulk/shell source maps to that propagation chain.

No small physical residual, finite-amplitude nonlinear family, observational
prediction, higher-dimensional origin of the Big Bang, novelty priority or
external specialist approval is established. This checkpoint creates no
Zenodo version or journal submission. Private editorial history is excluded.
