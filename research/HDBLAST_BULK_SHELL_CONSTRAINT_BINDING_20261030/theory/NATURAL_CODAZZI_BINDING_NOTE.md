# Exact natural bulk/shell Codazzi source-map binding

30 October 2026. This is a finite exact implementation audit of the unchanged
registered geometric sources. It binds them to the existing analytic physical
shell/Codazzi identity. It is not a new Codazzi theorem, a new evolution solve,
a nonlinear existence theorem or a cosmological-origin result.

## Objects, variables and retained physical rows

Use the static warped reference frame

```
g0=dy²+rho²[-dt²+c(t)² kappa_S3],
e_y=partial_y, e_0=partial_t/rho, e_a=bar_e_a/(rho c),
a=rho_y/rho, p=phi0_y, theta=c_t/c.
```

All source quantities are second parameter derivatives for
`g=g0+epsilon q+epsilon² k/2`, `phi=phi0+epsilon eta+epsilon² w/2`.
They are not the coefficients of epsilon². The normal points outward along
`+partial_y`. The natural equation is
`E_AB=Ric_AB-phi_A phi_B-(2/3)V g_AB`; it is not the reduced equation.

The affine first shell trace has

```
q_yy=hY, q_y0=q_ya=0,
q_00=-2XY, q_0a=V D_aY,
q_ab=2L Y kappa_ab+2d Hess_ab(Y), eta=vY.
```

Here the profile named `V` is a metric mixed amplitude; the scalar potential
is a different object. The actual adapter's seven amplitude slots are
`A=h,B=C=0,D=-2X,E=V,L=2L,M=2d,s=v`. Ordinary `y,t` derivatives of these
frame amplitudes are used throughout, not coordinate tensor derivatives and
not scaled `r,x` derivatives.

The one-copy outward physical residuals are

```
R_ab=K_ab-sigma(phi) gamma_ab/6,
P=n_out(phi)+sigma_prime(phi)/2.
```

The background rows are `R0=P0=0`, hence `sigma0=6a` and
`sigma1=-2p`. Let `beta=a h+sigma1 v/3`. The first physical metric rows give

```
D_y=-beta, L_y=beta, E_y=M_y=0.
```

Their ordinary tangential time derivatives supply the `yt` entries. The zero
mixed traces `B=C=0` imply their `t,tt` derivatives vanish. Their normal
`By,Cy` and mixed `Byt,Cyt` derivatives remain free. No shell relation is
normally differentiated. In particular `Ay` and all unused `yy` derivatives
remain independent formal variables. The harmonic normal equation for `Ay`
is not imposed.

For the main physical result also impose

```
s_y=p h/2-sigma2 v/2.
```

A stronger intermediate identity leaves that first scalar row unsolved:
the amplitude `p1=s_y-p h/2+sigma2 v/2` is free while `R1=0` is retained.
The full first scalar residual is `P1=p1 Y`, and `eta=vY`. The term
"off-shell" in checker row labels refers only to this **off-scalar-shell**
class. It does not claim that arbitrary nonzero `R1` has been omitted from a
fully off-shell expansion.

## The same actual shell Hessian

The existing affine unit-normal Hessian and its first physical-row
elimination give

```
A0=a h²/4-sigma1 h v/6-sigma2 v²/6,
QR_00=q00 Y²,             q00=-A0+2 beta X,
QR_0a=q0 Y D_aY,          q0=-beta V,
QR_ab=qa Y² kappa_ab+qb Y Hess_ab(Y),
qa=A0-2 beta L,           qb=-2 beta d.
```

The checker constructs this expected table independently, then verifies every
scalar/vector/tensor entry against the authenticated actual `shell_sources`
adapter. It also checks the Lorentzian trace, including the trace of the
angular Hessian:

```
tr_gamma QR=(-q00+T)Y², T=3qa-8qb.
```

Only the adapter's `QK00,QK0a,QKab` are consumed in this minimal-premise audit.
Its returned `QC_y,GN` are not used to assume a harmonic `Ay` restriction;
they are separate companion joins.

## Independent intrinsic divergence

Define `j=1/(rho c)` and `Hubble=theta/rho`. The shell connection and the
sphere identities

```
Delta Y=-8Y,
div Hess(Y)=-6 DY,
div(Y DY)=|DY|²-8Y²
```

give the physical-frame covector
`L(QR)_a=D^b QR_ba-D_a tr_gamma QR`. Its temporal primitive coefficients in
`(Y²,|DY|²)` are

```
(-T_t/rho-Hubble*(3q00+T)-8j q0, j q0).
```

Its angular coefficients in `(Y DY,Hess(Y).DY)` are

```
(j*(2q00-4qa+10qb)-q0_t/rho-4Hubble*q0, j qb).
```

The factor `rho^-1` converts ordinary `t` derivatives into the physical frame.
The `3Hubble` and `4Hubble` terms come from the Lorentzian shell FRW connection.
The companion independent review constructs the full four-dimensional
coordinate shell tensor, connection, divergence and trace explicitly, rather
than importing these primitive formulas or the producer arithmetic.

The fixed angular product projections are

```
Y²=(1/6,1/2,1), |DY|²=(4/3,2,-4),
Y DY=(1/4,1/2), Hess(Y).DY=(1,-2).
```

Thus there are three temporal scalar modal rows and two angular-gradient rows.

## Exact natural source-map equality and the scalar term

The existing exact identity is

```
E_na=D_b R^b_a-D_a tr R-P D_a phi.
```

It follows from geometric Codazzi after substituting
`K=R+sigma gamma/6` and `n phi=-sigma1/2+P`. The derivatives of the single
nonconstant shell tension cancel the background scalar normal derivative.
At the reference shell `D_a phi0=0`. With `R0=P0=R1=0`, its quadratic
coefficient is

```
QE_na=L(QR)_a-2 P1 D_a eta.
```

Normal and raising/connection variations multiply zero background or first
normal-tangential residuals. The fixed mixed shell traces make the perturbed
normal a pure normal rescaling; the first normal-tangential equation vanishes
from the same first Codazzi identity. No arbitrary first bulk PDE is needed.

The main calculation independently decodes the authenticated exported natural
`QE/y0` and `QE/ya` coefficient maps and checks all five coefficients of this
identity with `s_y` free. In the temporal scalar projection the matter term
has amplitude `-2 p1 v_t/rho` multiplying `Y²`; in the angular-gradient
projection it has amplitude `-2 p1 v/(rho c)` multiplying `Y DY`.
It then imposes the actual scalar row `p1=0`, equivalently `P1=0`, and
checks all five coefficients of

```
QE_y0=L(QR)_0, QE_ya=L(QR)_a.
```

The independently decoded coefficients also agree with the authenticated
registered `source_coefficients.evaluate_sources` implementation on both
coherent jet classes. Thus the natural bulk map and the shell adapter are
bound to the same QR source, rather than to unrelated upper bounds or interval
endpoint selections.

Every collected Laurent residual coefficient is exactly zero. The final
normal-tangential sources contain no free normal jets. `Ay,By,Cy` and unused
second normal derivatives cancel without being assigned values. Neither
background `a_y,p_y,Vprime,Vsecond,Vthird` nor a bulk ODE, harmonic `C1=0`,
reduced acceleration or the hyperbolic identity `c^-2+theta²=1` is needed.
Nonzero `rho,c` and static shell background parameters suffice for the
algebraic denominators/derivatives, together with the stated shell rows.

## Authentication and independent checks

Before project imports or evaluation, the checker authenticates the exported
map, actual evaluator, shell adapter and both relevant derivation notes.
Both imported project modules have standard-library-only dependencies.
The separately written sparse Laurent arithmetic helper has its own fixed
SHA256 pin and is authenticated again after use. No source file is edited.

The checker supplies 35 named collected identity rows. Their counts are
bookkeeping, not independent theorem counts: 12 shell-adapter entries,
3 Lorentzian traces, 5 off-scalar-shell source identities, 5 physical source
identities and 10 registered-evaluator joins.

There are 11 fresh corrupted executions in each mode. They test using reduced
QE or principal forcing in place of natural QE, losing the second-variation
factor, omitting the scalar Codazzi term **before** imposing P1=0, using a
positive rather than Lorentzian QR time trace, omitting either shell FRW
connection contribution, treating coordinate normal derivatives as frame
amplitude derivatives, reversing the scalar physical row or the background
scalar-tension sign, and dropping the rho frame factor in time derivatives.
Each must fail at its named nonzero polynomial residual and produce no success
report. Altering QP only after setting P1=0 would be an ineffective scalar-term
control and is not counted.

Normal and optimized mathematical reports must match complete bytes:

```sh
python -I -S -B checks/check_natural_shell_codazzi.py \
  --repo-root /path/to/repo --output CODAZZI.json --run-controls
python -I -S -B -O checks/check_natural_shell_codazzi.py \
  --repo-root /path/to/repo --output CODAZZI_OPTIMIZED.json --run-controls
```

## Scope

This strengthens a specific implementation-to-analytic-proof join for the
unchanged actual model. The existing all-point Codazzi theorem and physical
shell proof were already available. Application to the actual evolved first
field retains its proven physical shell trace/regularity identities. The
scalar/Bianchi, first-wave matching, subsidiary boundary conversion, smooth
regularity and uniqueness arguments remain separate identified inputs.

No point interval containing zero is used as an identity proof. No new
physical error bound, finite-amplitude nonlinear family, empirical Big Bang
origin, novelty priority or external specialist assessment is established.
