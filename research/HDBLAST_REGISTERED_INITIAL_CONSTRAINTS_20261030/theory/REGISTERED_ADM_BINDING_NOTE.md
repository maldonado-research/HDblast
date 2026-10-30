# Exact ADM binding of the registered actual initial assembly

30 October 2026. This is a source-binding calculation for the unchanged actual
second initial datum. It preserves the frozen data, all numerical candidates,
and the earlier analytic proofs. It is not a new existence theorem or a
finite-amplitude nonlinear solution.

The existing all-point lapse/shift argument in
`HDBLAST_COMMON_GAUGE_FORCING_20261016/inputs/EXPLICIT_INITIAL_HARMONIC_COMPLETION.md`
already proves the local initial harmonic prescription. The new calculation
checks the implementation join to the current
`HDBLAST_FULL_ANGULAR_RESPONSE_20261026/initial/full23_initial.py`, including its
five mean-slot overrides. It does not infer an identity from interval samples.

## Exact arithmetic and independence

The checker has its own standard-library collected Laurent-polynomial algebra
over the rational numbers and a parameter-series algebra through degree two.
Every test reduces the complete difference to the zero polynomial. Laurent
denominators use the retained nonzero `rho` and `p`; these are assumptions on
the actual regular background, not a new numerical enclosure here.

Before any project module import, the checker authenticates the full bytes of
both imported project modules (`source_adapters.py` and `linear_operator.py`),
the actual registration source, the two physical/ADM proof notes, and the
retained linear harmonic formula note/report. Those imported modules have only standard-library imports, so this is their complete
project-code dependency closure. The same files are reauthenticated after the
calculation.

The independent producer starts from the exact rational ADM formulas. The
implementation side compiles the authenticated `_registered` assembly suffix,
starting at `zero = packet*0`. Its physical inputs `u,u_y,u_t,v_t,f_l,f_{l,y},K,Pi`
are independent indeterminates in this calculation. No implementation
expressions are simplified or replaced in that suffix. This modular choice
makes the independence boundary explicit: a separate physical-initial audit
must bind the prefix that computes these profiles. In particular, this audit
alone does not prove `3K=Fmom` or `p Pi=Gmom-D_y Fmom`.

## Lapse, shift and metric expansion

Use ordinary parameter coefficients with the convention
`X(eps)=X0+eps X1+eps² X2/2`. Set

```
A = 1-4 eps U,
B = 1+2 eps U+eps² f/2,
P = 1+eps U,
N = rho P,
K_yy = -2 eps T/rho+eps² q/2,
K_ab/gamma0_ab = eps T/rho+eps² q/2.
```

The independent series expansion applies the exact equations (5)-(7) of the
retained ADM note. It recovers

```
dotN1 = rho T,
dotN2 = 4 rho² q-24 rho U T,
dotbeta1_y = -6 rho²(U_y+4aU),
dotbeta2_y = rho²[-74 U U_y-(3/2) f_y-194 a U²-3af],
dotbeta2_i = 22 U D_i U-(1/2)D_i f.
```

Substitution in the contracted-Christoffel ADM formulas gives zero at orders
zero, one and two for the time, radial and angular harmonic components. At
second order the generic zero identity implies the scalar degrees 0,2,4 of time and radial rows and
the gradient degrees 2,4 of the angular rows: eight modal rows. Those generic modal rows are consequences
of the complete zero-function identity, rather than eight additional
independent source-map evaluations.

Expanding `g_tt=-N²`, `dotg_tt=-2N dotN`,
`dotg_tI=gamma_IJ dotbeta^J`, `dotgamma=2NK` and `dotphi=N Pi` independently
recovers the complete spacetime second one-jet. In particular:

- `g2_tt=-2rho²U²` although `N2=0`;
- `dotg2_ty=dotbeta2_y-8U dotbeta1_y`, including the metric-shift cross term;
- `dotg2_tt` includes `-4N1 dotN1`;
- the scalar velocity includes `2U v_t`.

The coordinate-to-frame factors at the time-symmetric initial slice are
included explicitly. The checker then uses independent scalar product
coefficients `Y²=(1/6,1/2,1)` and gradient coefficients
`Y DY=(1/4,1/2)` and builds the canonical full23 order directly. It verifies
all 23 values and all 23 scaled velocities against `_registered`, and checks
the adapter and packer as a third join. Every ordinary `t` velocity is
multiplied by `ell` for `x=t/ell`. The five overwritten mean values and
velocities match the independent table. Both TT entries and both
trace-free-Hessian entries vanish initially.

## Direct initial harmonic source-map join

Separately from the generic ADM identity, the checker substitutes the corrected
first harmonic frame one-jet into the authenticated actual quadratic
`source_adapters.gauge_source` map. Its nonzero first metric coefficients are

```
A=-4u, D=-2u, L=2u,
A_y=-4u_y, D_y=-2u_y, L_y=2u_y,
A_t=-4u_t, D_t=-2u_t, L_t=2u_t,
B_t=-6rho(u_y+4au).
```

These are ordinary physical-frame `y,0,a` covector coefficients at `t=0`,
where `e_0=rho^(-1) partial_t` and `e_a=rho^(-1) bar_e_a`. The source-map
join has no `ell` factor; the separate initial PDE-array comparison checks
every ordinary `t` velocity is converted by `ell`.

The first harmonic rows vanish exactly for this jet. The complete actual
quadratic covector source map then gives

```
QC_y=-24u u_y Y²,
QC_0=24u u_t Y²/rho,
QC_a=24u² Y D_aY/rho.
```

The retained physical-frame linear `DC` formula is independently evaluated on
the generic registered second one-jet. It gives the negatives of these three
profiles, with the `f_l`, `f_{l,y}`, `q` and `Pi` dependencies cancelling.
This checks all eight actual scalar/gradient source-map coefficients, in
addition to the generic exact-ADM proof. It uses the full-metric-lowered source
convention of the authenticated map. Because the first harmonic covector is
zero, its extra lowering cross term vanishes. Replacing `B_t` by the
longitudinal value zero or reversing the radial `QC` sign is rejected.

The analytic match between the evolving actual first gauge wave and this
initial corrected one-jet is retained: it follows from the zero generator
Cauchy data and its wave-equation acceleration as described below. The direct
map calculation does not independently evolve that first wave.

## Initial subsidiary normal derivative

The checker independently verifies the five-dimensional trace-reversal
normal-projection algebra in a frame with timelike future initial normal `n`:

```
E_nn=X_n, E_nI=X_I/2, E_IJ=0, tr(E)=-X_n,
H_nn=E_nn+tr(E)/2=X_n/2, H_nI=X_I/2,
physical Hamiltonian=2 H_nn=X_n.
```

Thus zero physical constraints and `C=0` imply `nabla_n C=0` **provided the
reduced Einstein equations determine the accelerations**. At the initial
slice `partial_t=rho n`; terms from lower-order frame/normal variations
multiply already zero lower-order subsidiary or constraint jets. This is the
retained conditional algebraic implication, not a direct calculation of a new
`QC_x` enclosure or an independently prescribed acceleration.

The exact ADM prescription has nonzero first mixed metric velocity
`dotg1_ty=-6rho²(U_y+4aU)`. A longitudinal first one-jet with that velocity set
to zero is not interchangeable with it. Application to the registered evolving
first wave retains its analytic initial matching argument: zero gauge-generator
values/velocities and its wave equation give
`xi_tt^y=rho² C_long^y=dotbeta1_y`. This checker does not independently solve
that gauge wave or replace the existing neighbourhood harmonic identity.
With that match, second-order full-metric lowering has no extra first-gauge
cross term because the first harmonic covector already vanishes.

## Structural collar and limits

All finite registered initial coefficients lie in the ideal generated by the
first packet variables and the physical correction profiles. After the
separate physical prefix proof places those profiles in the packet-jet ideal,
all initial rows vanish on an open shell collar where the actual smooth packet
and all its derivatives vanish. Closure of a differential ideal under
finite differentiation explains the support mechanism. A finite one-jet
check is **not** by itself an all-orders corner proof; smooth flat support and
the retained homogeneous reduced recursion still supply that analytic step.

Global subsidiary evolution, physical scalar/source Bianchi relations, shell
Codazzi and Robin conversion, subsidiary uniqueness, finite-amplitude
nonlinear realizability, physical relative accuracy and a higher-dimensional
cosmological origin are not established by these initial algebraic identities.

## Verification

The checker runs in fresh `python -I -S -B` and `python -I -S -B -O` processes.
Both mathematical reports are byte-identical. There are 212 named collected
identities/joins, and each mode rejects 14 genuinely corrupted executions:
missing second metric lapse value, missing mixed metric conversion, halved
lapse cross term, omitted `f_y`, halved angular term, omitted lapse-velocity
cross term, missing `ell` scaling, changed mean override, wrong gradient
projection, wrong trace reversal, use of `E_nn` alone, wrong initial future
normal factor, use of the longitudinal first one-jet, and reversed radial
quadratic gauge-source sign. Counts are verification bookkeeping, not
independent theorems.

```sh
python -I -S -B checks/check_registered_adm_binding.py \
  --repo-root /path/to/repo --run-controls --output ADM.json
python -I -S -B -O checks/check_registered_adm_binding.py \
  --repo-root /path/to/repo --run-controls --output ADM_OPTIMIZED.json
```
