# Natural bulk/shell compatibility and first gauge-wave binding

30 October 2026. Exact implementation checks for the retained aligned
quadrupole. All earlier manuscripts, actual data, source maps and numerical
candidates remain unchanged.

The five exported natural Einstein normal-tangential source rows agree with
the intrinsic shell-source divergence. A complementary calculation binds the
registered first gauge reconstruction to its three-amplitude wave operator
and the corrected initial one-jet. These are exact source/implementation
joins to existing analytic proofs, rather than new Codazzi or nonlinear
existence theorems.

## Five natural Codazzi rows and the scalar condition

Use the one-copy outward residuals and second parameter derivatives

```text
R_ab = K_out,ab - sigma(phi) gamma_ab/6,
P = n_out phi + sigma_prime(phi)/2,
g(eps) = g0 + eps q + eps^2 k/2.
```

The established off-shell identity is

```text
E_na = D^b R_ba - D_a tr_gamma R - P D_a phi.
```

The new check evaluates the authenticated **natural QE** map on the same
coherent shell first jet used by the physical `QR` adapter. There are three
scalar `QE_y0` and two angular-gradient `QE_ya` coefficients. With background
`R0=P0=0`, zero mixed shell traces, and first physical metric rows `R1=0`,
it verifies

```text
QE_na = D^b QR_ba - D_a tr_gamma QR - 2 P1 D_a eta,
p1 = s_y - p h/2 + sigma2 v/2,    P1 = p1 Y.
```

Here `eta=vY` is the first scalar perturbation, while `sigma2` is the second
derivative of the shell tension. This scalar-residual calculation retains
`R1=0`; it is not a source identity with arbitrary first metric residuals.
Imposing the physical scalar row `P1=0` cancels the matter remainder in all
five modal rows. Omitting it generically leaves an explicit nonzero term.

The natural compatibility calculation leaves `Ay,By,Cy` and unused normal
jets free. It needs no first harmonic row, first bulk PDE, background bulk
ODE or de Sitter `1/c^2+theta^2=1` substitution. Only tangential derivatives
of shell trace equations are taken. A relation valid at the shell is never
differentiated normally. The independent check derives the intrinsic
Lorentzian FRW divergence directly from coordinate tensors and Christoffel
symbols before comparing the five exported map rows.

The evaluated quantity is `QE`, not `QEred`, a unit-principal forcing, or an
interval independently selected from the actual jet enclosure. The factors
of `rho` use ordinary physical time/frame derivatives.

## Physical shell sources, lowering and boundary normalization

A separate exact calculation expands the affine unit-normal metric/scalar
residuals and binds their source adapter to the raw numeric route and its
actual-shell wrapper. It checks the same first-jet interface

```text
h=A, X=-D/2, V=E, L=L_frame/2, d=M/2, v=s,
sigma0=6a, sigma1=-2p, sigma2=4(1+eta_b), sigma3=4.
```

It verifies the physical source coefficients, the Lorentzian trace, the
normal gauge source and all 23 primal forcing entries. With the first normal
derivative left free, it distinguishes full-metric and background lowering:

```text
QC_actual,y - QC_background,y = 2 H C1_y.
```

After imposing the first normal harmonic row, it verifies the same actual
`QC_y` shell adapter and

```text
GN = -2 tr_gamma QR - 2 QC_y,
G_K = -ell rho GN,
G_N = 2 ell rho QR00, G_u = 2 ell rho QR0a,
G_ab = 2 ell rho QRab, G_w = ell rho QP.
```

The sign change from outward `y` to inward `r` uses
`partial_y=-(ell rho)^(-1) partial_r`. One-copy and two-copy residuals and
second parameter derivatives are kept distinct. The 23-slot system has five
Dirichlet and eighteen Robin modal slots.

The induced subsidiary boundary conversion is also checked. Under the reduced
bulk equation and `E2_na=0`, the normal condition `C2_y=0` yields the tangential
condition `partial_y C2_a-2a C2_a=0` for lower coordinate covector components.
Its radially parallel frame version has coefficient `a`. These are induced
subsidiary conditions, not additional primal incoming boundary data.

## First reconstruction, gauge operator and initial match

The registered first metric/scalar pullback is compared with an independent
coordinate Lie derivative and frame conversion. The variable-`rho` physical
derivative conversion is included. Direct contraction of the first harmonic
operator then agrees with the registered raw gauge matrices and source.

For raw gauge amplitudes `(f,b,d)`, scaled coordinates `x=t/ell,r=z/ell`, and

```text
Res = U_xx-U_rr-P U_x-Q U_r-Rmat U,
s_f = -6 ell^2 rho^2 (u_y+4a u),
```

the normalization identities are

```text
ell^2 rho^2 DC_y = (s_f - Res_f) Y,
ell^2 rho DC_0 = Res_b Y,
(ell^2 rho/c) DC_a = -Res_d bar_D_a Y.
```

Here `DC_0,DC_a` are physical-frame covector components and `bar_D_a Y` is
the unit-sphere gradient.

The same audited operator, specialized to zero generator Cauchy data,
recovers `f_tt=-6rho^2(u_y+4au)` and the corrected first mixed metric
velocity used in the earlier initial `DC+QC` audit. The match retains the
analytic premise that the actual gauge wave has those Cauchy data and solves
its registered PDE. A coefficient identity is not a new PDE solution.

## Reproduce and assess the result

Use a complete repository checkout and a fresh output directory outside the
capsule:

```sh
python -I -S -B research/HDBLAST_BULK_SHELL_CONSTRAINT_BINDING_20261030/checks/run_checkpoint.py \
  --repo-root . --output-dir /tmp/hdblast-bulk-shell-binding-fresh
```

The distributed replay uses only Python's standard library; CI fixes Python
3.12.14. Complete reports must reproduce byte for byte in normal and
optimized modes. Deliberately corrupted formulas run in separate processes
and must fail at the intended mathematical identity without writing a success
report. Source closures are authenticated before use and after readback.
The written reviews identify which computations have separate arithmetic,
which use a shared arithmetic engine, and which analytic propositions remain
inputs. A private supplementary SymPy review is separate from this distributed
standard-library replay. All reviews here are internal AI reviews, without an
external specialist's assessment or journal acceptance.

This checkpoint strengthens the implementation links needed by the existing
subsidiary argument. It does not independently replace the complete
differentiated bulk Bianchi identity, first natural field equations, smooth
boundary regularity or subsidiary uniqueness proofs. Numerical error
certificates and `W1_quadratic` are unchanged. No physical relative accuracy,
finite-amplitude nonlinear family or higher-dimensional Big Bang origin is
established. No Zenodo version or journal submission is created.
