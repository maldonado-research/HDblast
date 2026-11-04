# A physical initial-profile class with conditional coefficient robustness

4 November 2026. Conditional follow-on; portable finite-rational component. The constants,
closed class radius, norm guard and final outward gates were fixed in
`FROZEN_DESIGN.json` before the new rational evaluation. Its SHA-256 is
`18b6defb16a5f79548666ceb539cc6152c6cfeaf36b68f4fc26c0cd31081900a`.
No gate or class radius was changed after evaluation.

This result retains the exact background, first packet, canonical second
data, second homogeneous preparation freedom, weighted constraint inverse,
conditional exact-family C2 dependence, coupled five-mean operator,
uniqueness, physical boundary rows, shell anchoring and target of the
preserved 4 November preparation-dependence theorem. The finite calculation
below does not establish these analytic premises.

## Geometric class and its reconstruction

The geometric identification is proved separately in
`../theory/PROOF.md`. On the specified
physical initial Cauchy section, continue the selected annular sphere as the
connected scalar-level component with fixed background scalar label
`v=phi0(y)`. Let `D` be the difference between its actual second area
derivative and the base family's second area derivative, divided by its
background area. Every quadratic term in the shared, generally nonspherical
first fields cancels in this difference. With raw homogeneous datum
`H_h=h gamma0+Lie_(zeta partial_y)gamma0`, raw scalar difference zero and
conditioned scalar datum `s=-p zeta`, the exact formula is

```
D=(3/2)(h+2a zeta)=(3/2)h-3a s/p,
s=k(h-2D/3), k=p/(2a).
```

There is no additional factor `1/2`: these are actual second parameter
derivatives. The area is invariant under diffeomorphisms transporting the
initial section, metric and scalar. Its label derivative is calibrated by the
fixed background,

```
D_r=(-ell rho p) partial_v D, ell=10^-6.
```

Thus this derivative is not an arbitrary perturbed radial coordinate choice.
The section and background calibration are part of the observable. A change
of physical section is not covered by an invariance claim.

Restrict to the existing admissible homogeneous directions: smooth `h`
supported in the packet, flat at both ends, with exact zero moment
`integral rho^2 h dy=0`, and with the same retained canonical second data.
This restriction does not assert that arbitrary C1 area profiles are
admissible. Define for the **complete added preparation**, on the causal
initial segment only,

```
R=max(sup_(1<=r<=5/4)|D|, sup_(1<=r<=5/4)|D_r|).
```

Let `chi_y=rho^2 h`, with the compact primitive supplied by the zero moment,
and put `z=chi/(ell rho^3)`, `b=ell rho a`, `g=ell rho p^2/(3a)`.
Using `rho_r=-ell a rho^2`, the local preparation constraint gives exactly

```
z_r=3bz-h,
h_r=3ell^2 z+(-b+g)h-(2/3)gD,
z(1)=h(1)=0.                                                    (1)
```

The first condition uses the exact zero moment, the second flatness. This
linear initial-value problem reconstructs `h` and `chi` on the causal
segment. It gives no bound on preparations outside that segment, which do
not affect the retained target by the same finite-propagation premise.

## Pinned raw coefficient checks

`produce.py` reads the following byte-pinned, preserved inputs directly:

| Input | SHA-256 |
| --- | --- |
| `HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/inputs/ACTUAL_ANNULUS_STATIC_JETS.json` | `3aa2c8d2e9ad652b72a34b4d8a6e6a6f14f72c32b8dcd1fbdc67556093265c22` |
| `HDBLAST_PREPARED_SHELL_VOLUME_20261101/response/references/FINAL_NORMAL.json` | `1a9ebffb6bee9332bc0fd0df0f5d1d9f7f2a796825c2a563e4c644974ba4d8a3` |
| `HDBLAST_SHELL_PREPARATION_DEPENDENCE_20261104/numeric/references/ORIGINAL_REPORT.json` | `78d63abd8caa6c183bd1d7ee7ba897be393bde14bcca1a37ee25073ea8d8b880` |
| local `inputs/FOLLOW_ON_REFERENCE.json` | `d9fd349165af7f7fd495b4979f1edaae49aaa75a3cadc0aaebbf8d1e07bdbb4c` |

All old capsule paths in the table are relative to the supplied repository's
`research/` directory. The local transfer reference contains only the two
operator constants and exact threshold. It is a cross-check, not an analytic
premise: the producer recomputes all three values from the pinned old raw
coefficient caps and exact base intervals.
The source's continuous-domain declaration covers every point of the packet
annulus, not sampled points. The script checks its initial slice, scale,
coverage and Taylor-jet normalization; first jet entries are actual first
derivatives since `1!=1`. Direct interval inclusions give

```
129<=rho,rho_b<=131, 3/25000<=-p<=1/7000,
111/1000<=a<=14/125, |p_y|<=1/4000, |a_y|<=1/10000,
|F'_b|<=1/1000, q=1201/3200000<1.
```

The selected background has not been replaced by a rounded point. Its exact
shell radius and force derivative remain enclosed quantities. The original
direction report is used only to cross-check the retained scale, target and
coupled contraction constant; no new direction or bump is evaluated.

## Rational ODE and conditioned-data bounds

Write `bM=ell*131*(14/125)` and
`gM=ell*131*(1/7000)^2/(3*(111/1000))`. Thus

```
bM=917/62500000, gM=131/16317000000000.
```

For the matrix in (1), its infinity norm is bounded by

```
max(1+3bM,3ell^2+bM+gM)<=mu=20001/20000.
```

For `x=mu/4=20001/80000`, the exact rational Taylor calculation checks

```
exp(x)<=sum_(j=0)^6 x^j/j! + (x^7/7!)/(1-x/8)<4/3.             (2)
```

The tail bound follows because ratios after the first omitted term are at
most `x/8<1`. Variation of constants, (1), (2), and `Delta=1/4` give the
continuous, whole-segment bounds

```
|z|,|h|<=H R, H=Delta*(4/3)*(2/3)*gM=131/73426500000000,
|h_r|<=Hr R,
Hr=3ell^2 H+(bM+gM)H+(2/3)gM.
```

There is no mesh, numerical ODE solution or exponential oracle in this
argument. The remaining exact identities are

```
k_r=-ell rho[p_y/(2a)-p a_y/(2a^2)],
s_r=k_r(h-2D/3)+k(h_r-2D_r/3), V_L=h_r+b h.
```

Use the rational absolute caps

```
kM=pM/(2aL), krM=ell rhoM[pyM/(2aL)+pM ayM/(2aL^2)],
S=kM(H+2/3), Sr=krM(H+2/3)+kM(Hr+2/3), V=Hr+bM H.
```

They bound `|s|/R,|s_r|/R,|V_L|/R`. The retained conditioned five-mean
initial difference, in `(K,L,N,J,w)` order, is exactly

```
W0=(h,0,0,h,s), W0_r=(h_r,0,0,h_r,s_r),
W0_x=(0,V_L,0,0,0).
```

For the positive Gram matrix `diag(1,2,1,3,1)`, triangle inequalities give

```
D0/R<=2H+S, Dr/R<=2Hr+Sr, Dv/R<=2V,
N=max((D0+T Dv)/3,Dr+Dv), T=5/4.
```

The script reconstructs the exact resulting constant

```
N/R<=798916710779559466175931203077/
     1861847711577000000000000000000000
   <0.000429099<1/2000.                                       (3)
```

The decimal in (3) is only an outward display; the reported rational and
the `1/2000` comparison use exact `Fraction` arithmetic. All five channels
and the mixed metric velocity are retained.

## Transfer to the two complete second observables

The retained homogeneous mean estimate is `X(W)<=N/(1-q)`. The initial
collar and finite propagation give the exact zero wedge `x+r<=1`, hence the
boundary scalar value obeys

```
|m(x)|<=(x-1)_+ N/(1-q).
```

The exact retained kernels satisfy `0<=k(T,x)<=T-x` and `0<=k_T(T,x)<=1`.
Integrating, with the exact response forcing factor `F'_b/2`, gives

```
|Delta v2|<=Cv N, Cv=(131ell)^2*(1/1000)*Delta^3/[12(1-q)]
                    =17161/767711760000000000,
|Delta Theta2|<=Ct N, Ct=131ell*(1/1000)*Delta^2/[4(1-q)]
                       =131/63975980000.                       (4)
```

For example `integral_1^T (T-x)(x-1)dx=Delta^3/6`; the scalar forcing's
`1/2` produces the denominator 12. These are complete coupled differences:
the identical actual first forcing, proper-clock and measure corrections and
first-quadratic Gauss terms cancel when subtracting families. Every such
term remains in the preserved absolute base intervals.

From the exact old negative upper endpoints and (4), the exact sufficient
threshold is their minimum ratio, reconstructed and checked against the
local transfer reference. Its outward display is

```
Ncrit>1.18646887e30.
```

For the frozen **closed** class `R<=2e33`, (3) yields `N<=1e30<Ncrit`.
The actual rational bound is stronger, below `8.5820e29`, but the final
gates use the originally fixed `1e30` guard. This gives absolute difference
caps below `2.23535e16` for volume and `2.04765e21` for expansion.
Adding both signs of the exact difference bounds to the exact preserved
base intervals gives the predeclared outward enclosures

| Whole class at the fixed proper target | Outward enclosure |
| --- | --- |
| `v2` | `[-7.6900e16,-4.1600e15]` |
| `Theta2` | `[-1.1070e22,-2.3980e21]` |

Both upper endpoints are strictly negative. The complete exact interval
endpoints are in `references/NORMAL.json`; the finite gates include both
outward containments and both strict signs.

## Replay and scope

The standard-library producer preserves **38 explicit finite checks** and no
`assert`-based validation. Isolated normal and optimized runs use
`-I -S -B` and optionally `-O`; their reports are byte-identical, with SHA-256
`16efc39b531ea3aea14aa6157f3f804df5d335a34294aa7cae7fb3fb9a6d00d5`.
`checks.py` passes **16** controls: **2 successful subprocesses, 12 actual
negative subprocesses, and 2 in-process checks**. Every negative checks exact
single-object JSON stdout, empty stderr, exit code 2, an exact rejection
witness, and absent output or preservation of preexisting bytes. All children
explicitly use the same isolated flags. All new oracle, PDE and bump-query
counts are zero.

Only top-level `scope`, `frozen_design`, `authenticated_inputs`, and the new
`provenance` differ from the source frozen report. Every other dictionary
field, including all numerical scientific fields and 38 check names, agrees
exactly. `references/FROZEN_SCIENTIFIC.json` is their canonical JSON
projection, SHA-256
`e940b780cb5a020371575f42deef7e080a69c0e605807d44f152c30893ea69ce`.
The source full-report SHA-256 was
`5b1b693b96087db54f419d232941587925eec1b68d5c41d9550526a25bbaa605`.
The producer pins and verifies this scientific projection hash at runtime.
No private input path is needed or emitted.

Replay into a fresh absolute output path with an existing parent directory, outside both the component and repository. Symlinked paths, parent escapes and overwrites are rejected:

```
python -I -S -B produce.py --repo-root /absolute/repository \
  --output /absolute/fresh/REPORT.json
python -I -S -B -O produce.py --repo-root /absolute/repository \
  --output /absolute/fresh/OPTIMIZED.json
python -I -S -B checks.py --repo-root /absolute/repository \
  --output /absolute/fresh/CONTROLS.json
```

This is a background-calibrated physical class of initial area **coefficient
differences**, in the inherited model and parameter units. It is not a
physical mechanism selecting this class or the base preparation. For a
direction multiplied by `lambda`, the actual class norm is
`|lambda|R(direction)`; the unscaled direction norm never controls arbitrary
`lambda`.

For each admissible fixed family, the retained C2 dependence and vanishing
first observable variations give the usual sufficiently-small-parameter
qualitative decreases from the background at the target. The class does not
control higher derivatives, nonlinear remainder sizes, preparations outside
the causal segment, or a uniform nonlinear parameter interval. No uniform
or effective finite-amplitude window, contraction, physical radiation or
cosmological mechanism follows from these coefficient inequalities.
