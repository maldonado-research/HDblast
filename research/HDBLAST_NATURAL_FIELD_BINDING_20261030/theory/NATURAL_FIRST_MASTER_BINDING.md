# Actual first master and the natural Einstein–scalar residuals

30 October 2026. This is an exact implementation-to-proof binding on the
retained regular collar. The sufficient reconstruction proof already exists
in `HDBLAST_REGULAR_CONE_EVOLUTION_20261013/reconstruction/`
`SMOOTH_AMBIENT_RECONSTRUCTION_AND_LOCAL_FLUX.md`, equations (3), (6)–(9).
The present calculation independently reconstructs the natural linear
operators and binds their use to the current actual first-field expressions.
It does not supply a new background, actual wave existence proof, nonlinear
solution or physical interpretation of the higher-dimensional hypothesis.

## Conventions and exact background inputs

Use outward proper coordinate `y` on one regular interior and

```text
g0=dy²+rho² gamma_dS4,    signature(gamma)=(-,+,+,+),
E_AB=Ric_AB-phi_A phi_B-(2/3)V(phi)g_AB,
S=Box_g phi-Vprime(phi),  p=phi0_y, a=rho_y/rho.
```

The retained actual smooth background solves, as local function identities,

```text
rho_y=a rho, p_y=Vprime-4ap, a_y=-rho^-2-p²/3,
a²=rho^-2+p²/12-V/6.
```

The Hamiltonian relation is essential. The checker substitutes
`V=6rho^-2+p²/2-6a²` and also verifies its derivative is `Vprime*p` under
the same differential recurrences. It is not enough to choose unrelated
interval endpoints satisfying selected differential bounds. Here `rho` and
`p` are nonzero on the retained collar; no claim at a zero of `p` is made.

`V` denotes the scalar potential, `Wsuper` its generating polynomial,
`mathcal W` the master-wave potential and `W(x,r)` the current actual master
amplitude. The current reconstruction coefficient
`alpha=(rho/rho_b)^(3/2)*p/p_b` is **not** the Sturm–Liouville weight called
`alpha=(rho²p²)^-1` in the old ambient reconstruction note.

## Independent connection and natural operator calculation

For an arbitrary smooth spacetime function `Psi`, reconstruct

```text
h_yy=-4Psi, h_y,mu=0, h_mu,nu=2rho² Psi gamma_mu,nu,
Pmaster=rho² Psi, chi=-3(Pmaster)_y/(rho²p).
```

The executable first constructs the background warped-product connection,
forms every component of `nabla h`, and computes

```text
P^A_ab = (g0^-1)^AB (nabla_a h_bB+nabla_b h_aB-nabla_B h_ab)/2.
```

It checks all 125 ordered entries against the universal tensor expressions

```text
P^y_yy=-2Psi_y,                 P^y_y,mu=-2D_mu Psi,
P^y_mu,nu=-rho²(Psi_y+6aPsi)gamma_mu,nu,
P^mu_yy=2rho^-2 D^mu Psi,       P^mu_y,nu=Psi_y delta^mu_nu,
P^mu_nu,lambda=delta^mu_nu D_lambda Psi+delta^mu_lambda D_nu Psi
               -gamma_nu,lambda D^mu Psi.
```

These tensor expressions are obtained **before** differentiating `P`.
At an arbitrary slice Riemann-normal-coordinate point, `gamma=diag(-1,1,1,1)`
and `partial gamma=0`. The reduced `P` expression has no slice Christoffel
left to differentiate. Differentiating it there and retaining the upper and
both lower tensor-slot connection terms therefore computes its full
covariant derivative; this is not a deletion of derivatives of a
coordinate connection. The calculation applies at every point. All ordinary
`Psi` jets through the required orders remain formal independent variables.
In particular it does not sample a solution or divide by a harmonic tensor.

Contracting
`deltaRic_ab=nabla_c P^c_ab-nabla_b P^c_ac` gives all 15 independent
entries, including

```text
deltaRic_yy=-4Psi_yy-16aPsi_y+2rho^-2 Box_gamma Psi,
deltaRic_y,mu=-3D_mu(Psi_y+2aPsi),
deltaRic_mu,nu=gamma_mu,nu[-Box_gamma Psi
  -rho²(Psi_yy+10aPsi_y+6(a_y+4a²)Psi)].
```

The tangential trace-free Hessian pieces cancel without a spectral
restriction. The same generated connection yields the nonzero longitudinal
harmonic covector `DC_y=-6Psi_y-24aPsi`, `DC_mu=0`.

The scalar variation is evaluated independently, directly as

```text
DS=Box0 chi-h^ab (nabla² phi0)_ab-DC[h]^a (dphi0)_a-Vsecond chi.
```

The `-DC.dphi0` term is retained. It is needed because this longitudinal
field is not harmonic. The result is

```text
DS=chi_yy+4a chi_y+rho^-2 Box_gamma chi-Vsecond chi
   +6pPsi_y+4Vprime Psi.
```

Define the already retained master operator

```text
M[Pmaster]=Pmaster_yy-(2a+2p_y/p)Pmaster_y
  +rho^-2(Box_gamma+4)Pmaster-(2/3)p²Pmaster.
```

All 15 natural Einstein entries and the direct scalar expression then give
the exact collected-coefficient identities

```text
DE_y,mu=0, DE_yy=2M/rho², DE_mu,nu=-gamma_mu,nu M,
DS=-3/(rho²p)*(partial_y+2a)M.
```

The scalar identity was not inferred solely from Bianchi. A value `M=0`
at one point leaves `-3M_y/(rho²p)` in `DS`; the needed master equation is a
neighbourhood identity, with its derivative, rather than a pointwise interval
containing zero.

## Binding the current actual reconstruction

The checker authenticates complete closed manifests before and after use.
It rejects symlink members/directories, extra or missing files and changed
bytes. It extracts only the needed authenticated source function/expressions;
it does not import any historical numerical producer. Its own Laurent helper
is hash-authenticated before import.

From `HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/actual_first_jets.py`
it executes the isolated `frame_amplitudes` function on formal exact inputs.
With the raw gauge vector zero this binds

```text
u=R*vartheta*W, v=B*vartheta*(W_r/ell+kW),
A=-4u, D=-2u, L=2u, Bmix=C=E=Mframe=0,
R=alpha/rho², B=3alpha/(rho³p), k=rho*(5a/2-Vprime/p).
```

It separately checks the actual scalar pullback `v+p*f`, the source temporal
log-derivative `vartheta_x/vartheta=-(3/2)ell*H`, the actual background
`rho_z,p_z,a_z,alpha_z` expressions and both source uses of `h=Vprime/p`.
The source alpha recurrence gives
`alpha_y/alpha=Vprime/p-5a/2`.
Consequently the implemented `v` is exactly
`-3(partial_y(alpha*vartheta*W))/(rho²p)`.

The same isolated-source execution binds all six ordinary physical-jet
adapter coefficients. In particular

```text
F_y=-F_r/(ell rho),
F_yy=F_rr/(ell rho)²+a F_r/(ell rho),
F_yt=-F_xr/(ell²rho), F_tt=F_xx/ell².
```

These are exact-function chain rules, not derivatives of numerical Taylor
remainders or independently selected interval values.

The current potential expression is extracted from the authenticated
`HDBLAST_INTERIOR_MASTER_JETS_20261015/potential/check_actual_potential.py`:

```text
Q=2(Vprime/p)²-10a(Vprime/p)+(35/4)a²-Vsecond-p²/6,
mathcal W=rho²Q-(35/4)H², H=tanh(t), vartheta=cosh(t)^(-3/2).
```

The **analytic actual master PDE input** is the displayed exact equation in
`HDBLAST_INTERIOR_MASTER_JETS_20261015/jets/HIGHER_JETS.md`, under “The exact
equation and initial data”: `U_xx-U_rr+ell² mathcal W U=0`, with
`U(x,r)=w(ell*x,ell*r)`. We call that same scaled actual amplitude `W` here.
Its smooth solution/existence and dynamic shell trace remain retained
analytic inputs. They are not proved by authenticating a file or by the
coefficient calculation. The helper's choice of the same PDE is explicit.

For the quadrupole `Delta_S3 Y=-8Y`, the independently collected radial
coefficient is `-Q+(5/2)rho^-2`; the temporal/spherical coefficient in
`(Box_gamma+4)(vartheta W Y)/vartheta` is
`-5/2+(35/4)H²` in addition to `-W_tt`. The two `5/2` terms cancel.
The executable uses the exact de Sitter quotient `cosh(t)^-2=1-H²` and
finds

```text
M[alpha*vartheta*W*Y]
   =-alpha*vartheta/(ell²rho²) Rmaster *Y,
Rmaster=W_xx-W_rr+ell² mathcal W W.
```

The direct scalar factor therefore becomes

```text
DS=3alpha*vartheta/(ell²rho^4 p)
   *[(Vprime/p-5a/2)Rmaster-Rmaster_r/(ell rho)]*Y.
```

Master derivatives through order three suffice. A free/comparator field
solving `W_xx-W_rr=0` leaves the generally nonzero defect
`ell² mathcal W W`; replacing the actual master by it cannot establish the
natural field equation.

## Transfer and limits

The preceding bulk-shell checkpoint already binds the current implemented
plus-sign Lie pullback from longitudinal to harmonic fields. Natural
covariance, `DE[Lie_xi g0,Lie_xi phi0]=Lie_xi E0=0` and similarly for `S0`,
is retained as the exact tensor identity that transfers the natural equation
to that first pair. The gauge wave's PDE supplies harmonicity separately;
it is not the reason pure-gauge natural residuals vanish. This audit does
not claim a new independent pure-gauge operator cancellation.

The original quadrupole `Y=n1²-n2²` is not identified literally with the
different degree-two quadratic response basis element `Z2`. This calculation
uses the original `Y` and its eigenvalue, with no full23 coefficient alias.
Dynamic shell trace/physical junction matching, quadratic bulk Bianchi map
binding, subsidiary uniqueness and regularity remain separately identified
tasks/analytic inputs. Existing proofs are preserved. No finite-amplitude
nonlinear spacetime, improved numerical error or cosmological observation
follows from this audit.

## Reproduction and effective controls

Run `check_natural_master.py --repo /path/to/HDblast --output report.json
--run-controls` under `python -I -S -B`, then under `python -I -S -B -O`.
The canonical reports must be byte-identical. Each wrong-formula control
runs in its own fresh subprocess, must fail its specifically named
mathematical identity and must leave no output report. Controls retain
nonzero witnesses for premature point-only scalar or comparator claims.
Their number is execution bookkeeping, not a count of new theorems.
