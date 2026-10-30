# Exact binding of registered actual initial profiles to the physical constraints

30 October 2026. This audit fills an implementation-to-proof join. The
Hamiltonian/momentum right inverses and their all-point analytic reductions
already appear in the retained 15–16 October construction. No new nonlinear
existence or cosmological result is claimed.

## Exact objects and assumptions

The implementation side is the authenticated AST prefix of
`HDBLAST_FULL_ANGULAR_RESPONSE_20261026/initial/full23_initial.py`, method
`Context._registered`, ending after its `Pi` assignment. Its code is evaluated
over an independent sparse Laurent algebra, with no import of its interval,
candidate or `Context` modules. The independent physical side transcribes
equations (3), (6), (7), (9), (10), (12) of the authenticated
`HDBLAST_COMMON_GAUGE_FORCING_20261016/inputs/EXPLICIT_COMPACT_SECOND_INITIAL_JET.md`.
This retains that note's physical constraint formulas; it does not rederive
their full coordinate Ricci Hessian.

The coefficient field is the exact rationals. Formal generators are
`rho,p,a,eta,alpha,E,Q0,...,Q8`; `Qj` represents the actual packet's ordinary
z derivative of order j, and `E=exp(-z)`. Negative exponents in the audited
expressions occur only for rho and p. The retained regular actual branch has
rho positive and p bounded away from zero; these analytic branch premises are
retained, rather than inferred from independent choices of interval points.
No nominal Taylor polynomial, rational packet, cutoff or mollifier is used.

The actual derivation is

```
Dz rho = -rho^2 a             Dz eta = -rho p
Dz p   = -rho(Vp - 4ap)       Dz a = 1/rho + rho p^2/3
Dz alpha = rho(5a/2 - Vp/p)alpha
Dz E = -E                    Dz Qj = Q(j+1)
Dy = -(1/rho) Dz.
```

The script extracts the authenticated numerical background source's `rhs`
dictionary and compares all five entries exactly with this derivation. It
also extracts that source's potential polynomial and its four differentiated
coefficient lists, checking equality with the independent polynomial and the
retained actual-potential JSON. These checks precede its use to evaluate the
registered prefix. The derivative of an unknown numerical Taylor remainder
never enters this calculation.

## Exact physical reductions

The independent producer begins with `u=alpha Q0/rho^2`, `ut=E u`, and the
coupled reconstruction

```
v  = -3(Dy u + 2au)/p,
vt = -3(Dy ut + 2aut)/p.
```

These agree exactly with the actual registered `v` and `vt`. In particular,
putting `h=Vp/p`, `ir=1/rho`, and `d=3/(rho p)` gives
`vt=E(v-du)`. The retained background equations give
`Dy p=Vp-4ap` and `Dy a=-1/rho^2-p^2/3`.

The first Hamiltonian coefficient is evaluated before cancellation as

```
-6u_yy -36a u_y -(24a_y+48a^2+4p^2)u
+20u/rho^2 -32u/rho^2 -2p v_y -2Vp v.
```

The distinct angular conformal-curvature contribution `+20u/rho^2` and
radial-lapse angular Laplacian contribution `-32u/rho^2` are both retained.
Substitution gives zero. The radial first momentum coefficient is
`[-3(u_yt+2a ut)-p vt]/rho=0`. For the angular first momentum equation,
`tr K1=(-2+3)ut/rho` and therefore the angular block of
`P1=K1-gamma0 tr K1` is exactly zero. The background scalar has no angular
gradient and Pi0 vanishes. This is the structural angular cancellation.

The independent physical Hamiltonian Hessian uses the complete unsimplified
`A,D` of equation (6). Collecting after the preceding substitutions proves

```
A = -2vy^2 -8p u vy +8p(h+2a)u v
    +(16p^2-192ir^2)u^2 -(2Vpp+8p^2/3)v^2
    -ir^2(12ut^2+2vt^2),
D = ir^2(36u^2-2v^2).
```

These are exactly the registered `aa,dd`. The unsimplified physical affine
momentum Hessian of equation (7) is

```
Fmom = (18u ut - v vt)/rho,
Gmom = (12u uyt -48a u ut -6uy ut -2vt vy)/rho.
```

The script proves `Fmom=3K_registered` and
`Gmom-Dy Fmom=p Pi_registered`. Thus it binds the actual simplified source
expressions, including the physical radial cross term, to the conventions
of the original analytic right inverse. This uses covariant physical K and
unit-normal scalar velocity Pi, not canonical momentum densities. Replacing
the radial physical source by the spacetime-affine Hessian removes
`-6uy ut/rho` and is explicitly rejected as a wrong-formula control.

For the scalar modes l=0,2,4, the actual registered metric coefficients are
checked through the complete linear constraint, not merely their stored
formula. With lambda=(0,8,24), `zeta=-3fy/(2p)`, the sum

```
deltaR = -3fyy -12a fy +2(lambda-3)f/rho^2,
deltaScalar = -2p zeta_y -2Vp zeta
```

equals the negative physical source sector. The monopole is included and
the denominators are `-6,10,42`. For `K2=K_registered gamma0` in four spatial
dimensions, the linear momentum operator is `-3dK-p Pi dy`. It cancels the
physical source in all three radial and both angular-gradient sectors.
These eight modal rows represent the five physical constraint components
for this fixed scalar/gradient angular family.

Every accepted difference has exactly zero collected Laurent coefficients.
The identity holds at all points where the displayed actual background and
packet identities hold, including the smooth zero extension by continuity;
the computation is not a finite collection of numerical point checks.

## Flat support and composition with the ADM audit

The algebra verifies membership in the packet-jet ideal for all first
profiles and all second sources/right inverses, including all three `f`,
`fy`, and `zeta` coefficients. First profiles have packet degree at least
one; second expressions have degree at least two. The ideal generated by
all packet derivatives is closed under Dz, and the only denominators are
smooth nonzero background factors. Therefore flat actual packet support is
preserved under any finite differentiation. Applying this structural fact
to all corner orders still retains the original smooth flat-support and
reduced-recursion argument; a finite Q0–Q8 table is not itself an all-orders
evolution or corner proof.

The separate ADM audit must bind the completed physical data to every actual
full23 value/velocity, including five mean overrides and ell time scaling.
This prefix audit does not by itself prove that later assembly, subsidiary
normal derivatives, shell Codazzi rows, Bianchi evolution or uniqueness.

## Reproduction and controls

Run `check_physical_initial_binding.py --repo-root ROOT --output REPORT
--run-controls` under `python -I -S -B`, then repeat with `-O`. It authenticates
every input and the shared arithmetic engine before use and after readback.
The normal and optimized report bytes agree. Each control is a separate
corrupted subprocess; success requires a nonzero exit, no output report,
and its specific expected nonzero identity. Controls change the
Hamiltonian denominator, scalar compensator sign, Pi cross term, scalar
velocity sign, Hamiltonian angular gradient, physical/spacetime-affine radial
momentum convention, spatial trace factor and y/z derivative scale.

This work establishes an independently reproducible registration join to
existing analytic physical-constraint formulas. It supplies no external
specialist assessment, finite-amplitude nonlinear solution, physical residual
smallness, new empirical prediction or higher-dimensional Big Bang origin.
