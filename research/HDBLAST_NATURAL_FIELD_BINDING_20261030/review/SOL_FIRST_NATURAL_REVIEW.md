# Independent first natural Einstein–scalar reconstruction review

30 October 2026. Private independent scientific and implementation review.
The independent result is **PASS**. No sign, coefficient, temporal constant,
Hamiltonian, ordinary-radial conversion, or scalar residual-derivative defect
was found in the registered sufficient reconstruction or in the frozen
producer reviewed below. Existing analytic sufficient-reconstruction proofs
are explicitly retained; this review does not claim a new existence theorem.

The executable is `independent_natural_review.py`, with its adjacent
`SOURCE_PINS.json`. It imports only Python standard-library modules. It does
not import the producer, its algebra helper, any project module, or a saved
scientific report. Selected authenticated source ASTs are evaluated in the
independent exact Laurent-polynomial algebra. The tensor and separated
calculations were completed before the frozen producer was replayed.

## 1. Independent natural-tensor construction

The tensor calculation uses a different coordinate realization from the
producer. Write unit dS4 in flat slicing,

```
gamma = -dt² + s² (dx1²+dx2²+dx3²),  s_t=s,
g0 = dy² + rho² gamma.
```

This is a local coordinate chart on the same constant-curvature Lorentzian
manifold. Every scalar coordinate jet of an arbitrary function `Psi` through
order three is an independent formal symbol. Thus this calculation is a
coefficient identity for arbitrary smooth functions, rather than evaluation
at selected points or a restriction to a particular separated harmonic.

Starting from the diagonal metric and its ordinary coordinate derivatives,
the checker constructs all background Christoffels and then

```
C^a_bc = (1/2) g0^ad (nabla_b h_cd+nabla_c h_bd-nabla_d h_bc),
deltaRic_bc = nabla_a C^a_bc - nabla_c C^a_ba.
```

It includes the upper and both lower connection slots in the divergence.
The timelike inverse metric has negative sign. Substituting
`h_yy=-4Psi`, `h_yμ=0`, `h_μν=2rho² Psi gamma_μν` gives all 15 independent
Ricci rows with exactly zero collected coefficient differences:

```
deltaRic_yy = -4Psi_yy-16aPsi_y+2rho^-2 Box_gamma Psi,
deltaRic_yμ = -3D_μ(Psi_y+2aPsi),
deltaRic_μν = gamma_μν[-Box_gamma Psi
 -rho²(Psi_yy+10aPsi_y+6(a_y+4a²)Psi)].
```

The six off-diagonal tangential entries are checked directly. In particular
the trace-free tangential Hessian cancels without dividing by a tensor
harmonic. This avoids any identification of the original
`Y=n1²−n2²` with the distinct degree-two function `Z2=n1²+n2²−1/2`.

The natural convention is

```
E_ab=Ric_ab-phi_a phi_b-(2/3)V(phi)g_ab,
S=Box_g phi-V_prime(phi).
```

The independent derivation uses

```
rho_y=a rho,  p_y=V1−4ap,  a_y=−rho^-2−p²/3,
V_y=pV1,  (V1)_y=pV2,  (V2)_y=pV3,
a²=rho^-2+p²/12−V/6.
```

Set `Pmaster=rho²Psi`, `chi=−3(Pmaster)_y/(rho²p)` and

```
M=Pmaster_yy−(2a+2p_y/p)Pmaster_y
  +rho^-2(Box_gamma+4)Pmaster−(2/3)p²Pmaster.
```

The full matter and potential terms then yield

```
deltaE_yμ=0,
deltaE_yy=2M/rho²,
deltaE_μν=−gamma_μν M.
```

The Hamiltonian relation is necessary for this reduction. Before imposing
it, the checker independently retains the explicit radial witness

```
deltaE_yy−2M/rho²
 =16Psi [a²−rho^-2−p²/12+V/6].
```

Consequently a collection of interval-selected ODE values satisfying only
some radial recurrences cannot be silently used as an exact Einstein
background. The Hamiltonian corruption subprocess fails with this exact
nonzero polynomial remainder.

## 2. Direct scalar calculation

The scalar row is constructed independently of a Bianchi argument:

```
deltaS = Box0 chi − h^ab Hess0(phi0)_ab
         − (g0^bc C^a_bc)(phi0)_a − V2 chi.
```

Here `Hess(phi0)_yy=p_y` and
`Hess(phi0)_μν=ap g0_μν`. The longitudinal divergence contraction is
nonzero: `g0^bc C^y_bc=−6Psi_y−24aPsi`. Retaining it gives

```
deltaS=chi_yy+4a chi_y+rho^-2 Box_gamma chi−V2 chi
       +6pPsi_y+4V1Psi
      =−3/(rho²p) (partial_y+2a)M.
```

Every equality is checked by direct collection of exact coefficients. The
scalar reconstruction already contains one master derivative; its direct
natural equation requires master jets through order three. A master PDE
identity on a neighbourhood supplies its radial derivative. A vanishing
master defect at one point alone does not supply this scalar premise.

## 3. Binding the current actual separated field

The checker evaluates isolated authenticated ASTs from
`actual_first_jets.frame_amplitudes`, with independent exact jet objects.
It obtains the actual longitudinal amplitudes

```
u=alpha theta W/rho²,
v=3alpha theta/(rho³p) (W_r/ell+kW),
A=−4u, D=−2u, L=2u, Bmix=C=E=Mmetric=0.
```

Here `theta=cosh(t)^(-3/2)`, `x=t/ell`, `r=z/ell`, and inward depth satisfies
`partial_y=−(ell rho)^(-1)partial_r`. In source names, the scalar coefficient
`Bcoef=3alpha/(rho³p)` is distinct from the mixed metric amplitude Bmix.
The actual `interval_core.background` alpha RHS and final `k,C,B` ASTs give

```
alpha_y/alpha=V1/p−5a/2,
k=alpha_z/alpha=rho(5a/2−V1/p),
Rcoef=alpha/rho², Bcoef=3alpha/(rho³p).
```

These prove the same exact scalar identity
`vY=−3(Pmaster)_y/(rho²p)`, with
`Pmaster=alpha theta W Y`. They do not differentiate independently chosen
coefficient intervals or unknown Taylor remainders.

The actual `physical_derivatives` AST is also evaluated and checked in all
six ordinary derivative rows. Its essential second-radial formula is

```
F_yy=F_rr/(ell²rho²)+aF_r/(ell rho).
```

The independent engine obtains the extra term by differentiating variable
rho. Omitting it fails a dedicated corrupted subprocess.

The retained potential is authenticated as the same registered polynomial:

```
Wsuper=1/3+eta²+eta³/3,
V=(1/2)(Wsuper_eta)²−(2/3)Wsuper²,
V coefficients=[−2/27,0,14/9,50/27,−1/6,−4/9,−2/27].
```

The checker reconstructs this polynomial independently and compares it to
the `interval_core.U` AST. It separately binds both formal and numerical Q
ASTs from `check_actual_potential.py` to

```
h=V1/p,
Q=2h²−10ah+(35/4)a²−V2−p²/6,
mathcalW=rho²Q−(35/4)H²,  H=tanh(t).
```

Let `Lalpha=alpha_y/alpha=h−5a/2`. Before simplification, the radial
transformation gives

```
(Lalpha)_y+Lalpha²−(2a+2p_y/p)Lalpha−(2/3)p²
  =−Q+(5/2)rho^-2.
```

Also `2Lalpha−2a−2p_y/p=a`, so the radial derivative part is
`W_yy+aW_y=W_rr/(ell²rho²)`: the two first-radial terms cancel only after
the variable-rho term has been retained.

For the original l=2 harmonic, `Delta_S3 Y=−8Y`. Since
`theta_t/theta=−3H/2` and
`theta_tt/theta=−3/2+(15/4)H²`, direct temporal differentiation gives

```
(Box_gamma+4)(theta W Y)/(theta Y)
 =−W_xx/ell²+[−5/2+(35/4)H²]W.
```

The first-time derivative coefficient vanishes. The separate `+5/2` radial
and `−5/2` temporal terms cancel. Both are checked separately before the
complete master conversion:

```
Rmaster=W_xx−W_rr+ell²mathcalW W,
M[alpha theta W Y]=−alpha theta/(ell²rho²) Rmaster Y.
```

Therefore the natural residuals are exactly

```
deltaE_yμ=0,
deltaE_yy=−2alpha theta/(ell²rho⁴) Rmaster Y,
deltaE_μν=alpha theta/(ell²rho²) gamma_μν Rmaster Y,
deltaS=−3alpha theta/(ell²rho⁵p)
       [Rmaster_r/ell+k Rmaster]Y.
```

The scalar equality is checked with a separate explicit product-rule
expansion of `Rmaster_r`, leaving `W_xxr` and `W_rrr` independent and
retaining the radial derivative of the full potential. It is not checked
only with a fresh independent symbol for an unspecified derivative.

A free comparator satisfying `W_xx−W_rr=0` still has defect
`ell²mathcalW W`. Its exact nonzero polynomial is retained in the report;
claiming this comparator is the actual master is explicitly rejected.

## 4. Exact PDE versus numerical enclosures

The exact smooth master equation and its temporal `tanh²` coefficient are
retained defining analytic inputs in authenticated `HIGHER_JETS.md` and
`ACTUAL_POTENTIAL_WZZ_DERIVATION.md`. The interval source implementation
does not solve the PDE or select a correlated exact background point.
`interval_core.background` returns background and reconstruction jets; it
does not evaluate the complete time-dependent master potential.

The actual numerical-potential expression is separately AST-bound to
`rho²Q−[0,(35/4)T²]`. `check_higher_jets.py:317` uses the same static Q
enclosure and the corresponding `35(Xell)²/4` temporal majorant. These
enclosures are consistent with the exact coefficient definition because
`|tanh(t)|<=t`; they are not a replacement for the analytic equation.
No inference that arbitrary independent interval selections satisfy the
master or natural equations is made.

## 5. Current harmonic representative

The complete 22-member PR106 source capsule is authenticated before the
calculation and after readback. Its retained independent audit already
binds the current implemented full frame pair to

```
(h_long+Lie_xi g0, chi_long+Lie_xi phi0).
```

This review directly rechecks the actual plus scalar shift and tests a live
isolated-source AST mutation to its minus-sign alternative. For the full
metric pullback it uses that completed independent PR106 binding, rather
than claiming to reproduce the full first-wave audit here.

The transfer of natural residuals uses standard exact natural covariance:

```
DE[Lie_xi g0,Lie_xi phi0]=Lie_xi E0=0,
DS[Lie_xi g0,Lie_xi phi0]=Lie_xi S0=0.
```

No new generic pure-gauge tensor calculation is claimed in this review.
The first gauge-wave PDE is needed to obtain harmonic `DC1=0`; it is not
needed for this natural-covariance transfer. Tangential gauge values b,d
need not vanish at the shell. Dynamic shell trace equations, first shell
junction binding and quadratic Bianchi/source specialization are separate
joins and are not silently inferred from this bulk calculation.

## 6. Authentication, rejection controls and reproduction

The fixed pin file is itself hash-checked before use. Thirteen directly read
scientific files are pinned by size and SHA-256. Complete closed-tree source
closures are checked for all members, extra/missing files and symlinks,
before evaluation and again after readback:

| Capsule | Members excluding manifest | Manifest SHA-256 |
|---|---:|---|
| REGULAR_CONE_EVOLUTION_20261013 | 73 | `a873132efe11755217ddb865fd0156041b8c1ad21ae1c45310183da4fb6140ae` |
| REGULAR_SHELL_SPECTRA_20261013 | 129 | `d1d47c66cef1ed79979c22ae7145dc66d35fa35be69dba900f9912ac01f31e2f` |
| INTERIOR_MASTER_JETS_20261015 | 40 | `b043a13a228d11296fbbde292aa15f7a024deab9ff915bdce452c52c73745bbb` |
| BULK_SHELL_CONSTRAINT_BINDING_20261030 | 22 | `7eae74dedfef544199d410665deb5ebf013d4cf484745428aba37dac7ba0671e` |

No imported-scientific dependency remains outside these selected AST/source
closures. Canonical AST hashes for the actual source functions and isolated
expressions are included in the independent report.

Fresh normal and optimized runs both pass 72 exact checks and 22 actual
corrupted subprocess controls. The controls require their precise intended
failed identity, nonzero exit status, absence of an output report, and an
exact nonzero rational polynomial remainder. They cover chi/alpha/k signs,
both potential matter terms, the nonharmonic scalar divergence, Hamiltonian,
angular eigenvalue, master +4, theta differentiation, either cancelling
constant, three Q coefficients, inward radial sign, master and scalar ell
powers, variable rho, comparator substitution, point-only scalar reasoning,
and the live source scalar Lie-sign mutation. Validation uses no `assert`.

```
python -I -S -B independent_natural_review.py --repo-root ROOT \
  --output INDEPENDENT_NORMAL_FINAL.json --run-controls
python -I -S -B -O independent_natural_review.py --repo-root ROOT \
  --output INDEPENDENT_OPTIMIZED_FINAL.json --run-controls
```

The two complete independent reports are byte-identical with SHA-256
`22f131475d4b45afa1ba8a810a8667a54e5e6f6d688ec76e3e084ba9ee190c21`.
The independent checker SHA-256 is
`018cc9744597cb6d9cde238567f088d5530e0551e6dca2f15f9f2d36f12f1ac2`;
its fixed pin file SHA-256 is
`b67a008b2a6f8799190ed87e0239fdbf4e59b5ee93a054ae04eaf00482d62960`.

## 7. Frozen producer review and reviewer replay

The separate frozen producer was read and reviewed after the independent
derivation. Its slice-Riemann-normal-coordinate route constructs the same
connection tensor, includes the derivative slots and Lorentz signature,
and directly computes the scalar equation. Its use of slice normal
coordinates is valid here: the simplified connection variation is a tensor
built from the slice metric, arbitrary scalar jets and warp data; the
vanishing slice Christoffel at the evaluation point is not used to discard
a required derivative of a non-tensor connection term. The independent
flat-slicing derivation supplies a separate coordinate check.

No mathematical or implementation defect was found in the frozen sources:

```
check_natural_master.py
1c48117b341fd6c80c468b7facee921b21936d9352d339a0677ea2b3f2ed4669
natural_algebra.py
99231502fd782ebb5d4c4a5a684370dd61ad6b2b8dc1afa5806dc2bacaa5d829
```

The producer authenticates its helper before import and all five retained
source capsule closures before scientific use and after readback. Fresh
reviewer normal/-O executions pass its 72 checks and 22 intended corrupted
subprocess failures in each mode. The reviewer outputs are byte-identical,
match the producer's frozen report, and have SHA-256
`7ed5135d18a89c71460b25e8b78f43f4e821b8ac15ec7e0e8150d6b64a2facbf`.
The independent checker additionally supplies the explicit unreduced
Hamiltonian witness, original Wsuper-to-V polynomial join, numeric-potential
AST joins and expanded `Rmaster_r` check described above.

These are implementation-to-proof coefficient identities. Actual background
and master existence, regularity, the neighbourhood PDE identity and
subsidiary uniqueness remain identified analytic premises. Nothing here
proves a finite-amplitude family, improved physical error sizes, completeness,
a cosmological origin, literature novelty or external peer review.
