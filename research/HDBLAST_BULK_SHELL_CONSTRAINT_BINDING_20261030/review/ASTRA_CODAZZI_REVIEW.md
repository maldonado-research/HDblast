# Adversarial review of exact natural bulk/shell implementation binding

30 October 2026. Internal independent mathematical and source-code review.

**PASS for the stated implementation-binding scope.** No unresolved
mathematical mismatch was found in the final natural Codazzi checker or its
physical-source/gauge/normalization companion. The five actual natural map
rows agree exactly with the divergence of the same actual shell QR adapter.
This is a binding of implementations to retained analytic identities, not a
new Codazzi theorem, first-wave solve or finite-amplitude existence result.

## Reviewed artifacts and fresh replay

The source hashes reviewed and independently replayed are:

| Artifact | SHA256 |
|---|---|
| `check_natural_shell_codazzi.py` | `c3d260c951d54a3e47bf2b6a62d096d49e5dff2cb8811c0fed04d69546754b31` |
| `codazzi_laurent.py` | `178b0a4e58db63eae6c24112ed9d91bf430f70dca9f6f8126c5e352e429f8fc2` |
| `check_shell_source_bridge.py` | `bff18f81fba6ebb514be6605208b550714696add65c3ae1692780a6dfa787514` |
| Complete natural Codazzi report, either mode | `e1a44b5c3125b9f15f7dc708a97ac0940189c3be393af25c8e1f7a9f9e82514b` |
| Complete shell bridge report, either mode | `89100b8a4755d759832018c09037c1398e2d2ec571ec0f13ceec838e1fb8fd7f` |

Fresh review executions used Python 3.12.14 with `-I -S -B`, and again with
`-O`, including all corrupted subprocess controls. Each reproduced the
corresponding finalized complete report byte for byte. The main report has
35 named joins and 11 rejected corruptions per mode; the bridge has 157 joins
and 18 rejected corruptions per mode. These counts describe verification
bookkeeping. Source hashes were checked before and after these review runs.
The private replay receipt preserves the actual command arrays and hashes.

The final accompanying notes `NATURAL_CODAZZI_BINDING_NOTE.md` and
`SHELL_SOURCE_BRIDGE_NOTE.md` were also read. Their qualifications match the
executable premises. Replaying a producer is an independent process check,
not a second independent mathematical algorithm.

## Independent mathematical findings

I reconstructed the intrinsic frame contraction after reading the proposed
formulas. For the shell metric `rho^2[-dt^2+c^2 kappa]`, put
`Hubble=theta/rho`, `j=1/(rho c)`, and `T=3qa-8qb`. In the temporal
divergence the `q00_t` terms cancel between divergence and trace derivative.
The remaining spatial connection contributes
`-Hubble*(3q00+T)`. In the angular divergence, contraction of the connection
on its two tensor indices contributes `-4Hubble*q0`. The sphere identities
`Delta Y=-8Y` and `div Hess Y=-6 DY` then give

```text
temporal: (-T_t/rho-Hubble*(3q00+T)-8j*q0, j*q0),
angular:  (j*(2q00-4qa+10qb)-q0_t/rho-4Hubble*q0, j*qb).
```

These are the correct Lorentzian-frame primitives. The independent intrinsic
review uses a stronger separate implementation: explicit coordinate metric,
Christoffel symbols and tensor contraction, without importing the main
checker or its arithmetic. Its written derivation is consistent with these
contractions.

The actual map inputs are ordinary derivatives of frame amplitudes. Thus
`D_y=-beta`, `L_frame,y=beta`, `E_y=M_y=0` are correct, because the
coordinate tensor derivative includes the additional `2a*m` warp term.
Only tangential derivatives of these shell rows are taken. The final code
does not infer any `yy` relation from them. `Ay`, `By`, `Cy` and unused normal
jets remain formal and cancel. The natural five-row identity requires no
normal harmonic equation, reduced acceleration or background bulk ODE.

The stronger intermediate map check is specifically **off scalar shell**:
it retains `R1=0`, while leaving

```text
p1=sy-p*h/2+sigma2*v/2,    P1=p1*Y
```

free. Here `p1` is the scalar harmonic amplitude and `P1` is the full first
scalar residual. Its exact result is `QE_na=L(QR)_a-2P1 D_a eta`. It is not an arbitrary
`R1`-free off-shell theorem. A separately written private monomial decoder
and SymPy calculation verified all five intermediate rows and their five
`P1=0` specializations. It imported neither producer nor their arithmetic.
Its angular conventions and the frozen source map are shared inputs;
independence is computational and derivational, not source-blind. This
supplement used SymPy 1.14.0 and remains private; it is not a dependency or a
stdlib/`-S` claim for the reproducible capsule.

Deleting the scalar term before imposing `P1=0` gives the concrete temporal
mode-zero residual

```text
vt*(h*p-sigma2*v-2*sy)/(6*rho) = -p1*vt/(3*rho).
```

Thus this is an effective scalar cancellation control. Changing `QP` after
the scalar factor has vanished would not test this cancellation.

## Quadratic normal, lowering and connection terms

With zero mixed first shell traces, the affine unit normal is rescaled by
`1-epsilon*H/2+3epsilon^2*H^2/8`. The second pulled-back natural mixed
residual is therefore

```text
QE_na-H*E1_na+(3/4)*H^2*E0_na.
```

The extra terms vanish because `E0_na=0` and the first Codazzi identity gives
`E1_na=0` when `R1=0`; this remains true with `P1` free because
`D_a phi0=0`. Variation of the shell divergence/raising multiplies `R0` or
`R1` and also vanishes. This explains why these terms are absent in the
restricted source join, rather than silently dropping them off shell.

The companion independently expands the physical normal factor and the full
relative-gauge expression. Before imposing the normal harmonic row, the
correct lowered normal Hessian is

```text
QC_actual_y=-H*Hy-a*H*tau+sigma1*eta*tau/3+2a*tr(A^2),
QC_actual_y-QC_background_y=2H*C1_y,
C1_y=Hy/2+2aH-a*tau-2sigma1*eta/3.
```

The final code verifies both map lowerings and their difference with `Hy`
free, then specializes the retained harmonic row for the shell adapter.
The reference-connection term is included. Its targeted corruption preserves
background and first order, removing only the quadratic contribution, and
fails at `GENERIC_QC_ACTUAL_0`. The Lorentzian contraction retains both
`-2V^2|DY|^2` and `4d^2|Hess Y|^2`.

The same sources satisfy `GN=-2 tr_gamma(QR)-2QC_y`. The companion verifies
the numeric amplitude-name factors, fixed tension derivatives
`(-2p,4(1+eta_background),4)`, one-copy/two-copy normalization, all 23 packed
entries, and the inward conversion `partial_y=-(ell*rho)^-1 partial_r`.
The induced covector Robin coefficients are correctly `2a` in coordinate
components and `a` in the radially parallel frame. They are not additional
primal incoming data.

## Resolved control issues and remaining scope

Two control refinements were requested during review and are present in the
final hashes. The temporal trace control now actually replaces
`tr_gamma QR=-q00+T` by `+q00+T`, changing both temporal and angular
derivatives. The reference-connection control now removes its quadratic
term with lower orders preserved. The companion's generic lowering check
also prevents assuming `QC_actual=QC_background` before harmonicity.
These were improvements to what the tests diagnose; neither revealed a
baseline source-map defect.

The numeric shell wrapper is exercised with a symbolic substitute for the
interface that supplies its actual first trace. This verifies the wrapper's
argument, tension and normalization conversions. It does not execute the
underlying interval jet producer or establish that an evolved field satisfies
the physical trace equations. The existing analytic first-wave shell rows,
regularity, Bianchi/subsidiary propagation and uniqueness remain identified
inputs. The separate first-wave binding must be assessed on its own evidence.

No interval containing zero substitutes for a coefficient identity here.
The accepted result is exact source/proof consistency on the declared
coherent shell jet class. It establishes no new physical accuracy bound,
nonlinear family, novelty priority or cosmological-origin claim.
