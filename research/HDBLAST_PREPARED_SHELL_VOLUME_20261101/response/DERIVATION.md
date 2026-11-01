# New explicit-family mean evolution and intrinsic shell observables

1 November 2026. Private conditional computation. The registered target is
`T=x=5/4`, with `ell=10^-6`. Coefficients are parameter derivatives,
`g_e=g0+e g1+e² g2/2+…`; the same convention applies to every observable.
The source worktree has HEAD `c1a07319e790cda389d7649d686e24e2f317fbe0`,
containing the merged main result `27dedaddff80a6539ba18fbdc6c200e2884cab4e`.
No retained capsule is modified.

## Preparation and conditioning

The family is the deliberately selected explicit bump preparation of
`../preparation/EXPLICIT_MEAN_DESIGN.md`, with
its corrected canonical momentum `tau_angular=-QM/(6a)`. It is not the
old weighted-only mean family. The `N1024` numerical campaign fixes the
same background shell boxes, first seed, bump, moments and compact primitive.
It supplies whole-cell interval bounds, including flat tails; no pointwise
sampling is used as a replacement for a continuous bound.

The original second spatial metric is `f gamma0+Lie_(xi d_y) gamma0`.
Its `A=f+2xi_y` can reach approximately `10^39`. Direct evolution of
`B=f+2a xi` has severe shell cancellation: the first Robin lift contains
`ell rho_b a_b integral A dx`, whose primitive cancels the leading `2a xi`.
The old pilot's initial data, total cap and accuracy allocation are therefore
irrelevant to this new preparation.

Use the checked second-only harmonic gauge adapter
`../gauge/CONDITIONED_GAUGE.md`. Initially
its vector has `(Z^y,Z^t,Z^y_t,Z^t_t)=(-xi,0,0,a xi)`. Its homogeneous
monopole gauge wave has normal Dirichlet and tangential Neumann shell rows.
The gauge vector and all initial jets vanish on the same shell collar.
The physical shell, initial anchor, Cauchy geometry and observables are
transported together. The first solution and its quadratic sources are
unchanged. The scalar acquires `-p xi`; it must not be set to zero.

In mean order `(K0,L0,N0,J0,w0)`, Gram diagonal `(1,2,1,3,1)`, the adapter is

```
U0=(f,0,-u²/3,f,-p xi),
V0/ell=(2rho QM/(3a)+32u ut/3,
        rho[-fy+a f-13u uy/3-a u²/3],
        -2rho QM/(3a)+4u ut/3,-4u ut/3,0),
U0_r=-ell rho (fy,0,-2u uy/3,fy,-py xi-p xiy).
```

The response producer independently evaluates these formulas with exact
rational endpoint arithmetic and checks inclusion in every separately
produced outward-rounded conditioned cell. It uses the authenticated
actual global `py` interval. Only the new-family data enter the calculation.

## Finite causal problem and complete error budget

The five means are an exact invariant SO(4) projection. The authenticated
23-amplitude primitive matrices have no entries to or from the five slots;
the producer verifies this and the Gram diagonal. Internal mean metric,
scalar, connection and Robin couplings are retained through the resolvent.
No whole angular solve is required to integrate these mean observables.

Restrict to the backward-closed cone

```
D_T={0<=x<=T, 0<=r<=T-x}, T=5/4, delta=T-1=1/4.
X(U)=max(sup |U|_H/3, sup |U_x|_H, sup |U_r|_H).
```

Both new initial data and the unchanged quadratic first-pair sources vanish
in the initial flat wedge `x+r<=1`. Finite propagation of the retained
compatible problem therefore gives the same exact wedge for the total
response and for differences from its free solution. This is a statement
about the linear coefficient problem on the retained background cone; it
is not an assertion about a finite nonlinear parameter cone.

For any compatible wedge-zero `z`, integration from the wedge gives
`|z|_H<=delta X(z)` on `D_T` and
`|z(x,0)|_H<=(x-1)_+ X(z)`. The reflected retarded kernel has
`X(Kf)<=T sup|f|_H`; the shell lift has `X(Lg)<=sup|g|_H`.
The retained coefficient caps imply

```
q=T(2*10^-4+delta*10^-6)+delta/2000
 =1201/3200000=0.0003753125.
```

The original actual bulk and all physical/normal-gauge shell forcing
histories remain included. Their authenticated caps are
`FH=17387513156573`, `GH=3851400`; hence
`S_X<=T FH+GH`. They are charged once. Their signs are not assumed, and
no comparison forcing replaces the actual forcing.

Let `D0,Dr,Dv` be uniform H upper bounds for the new initial displacement,
radial slope and scaled velocity. Reflection is an H isometry and the
flat collar removes corner defects. The actual, unevaluated free wave `B`
satisfies

```
B_X<=max((D0+T Dv)/3,Dr+Dv).
```

The exact total response satisfies `U=B+TU+S`, so

```
E_X=(q B_X+S_X)/(1-q) >= X(U-B).
```

This encloses every feedback iteration, every mean coupling and the complete
source history. No old-family `Ucap`, initial reconstruction error, P1 data,
two-term tail or omitted forcing allocation is reused. No numerical field
solver is claimed. The free scalar trace is exactly `B_w(x,0)=w_initial(x)`,
because its initial velocity is zero and its reflection is even.

## Proper scalar mean and full first-order correction

`M2` is the second parameter derivative of the shell scalar's normalized
proper-volume average on the same anchored intrinsic normal-geodesic slice.
For the retained quadrupole, with ordinary proper time `s=rho_b ell x`,

```
M2=w0-Cclock psi_s/3+(E-8d/3)psi,
Cclock_s=X, E=Z-H Cclock.
```

Here `psi` is the first shell scalar, `m00=-2X Y`,
`mab=2(ZY kappa+d Hess Y)`, `m0a=V D_aY`, and
`<Y²>=1/6`. The last term includes the perturbed volume measure. It is
essential even though its magnitude is tiny in this calculation. The
`d` term is not dropped. Second angular transport cancels between scalar
pullback and proper measure; this does not remove the displayed measure
correction.

The actual shell 2-jet oracle is evaluated on the whole rectangle
`x in [1,5/4], r=0`. Its raw metric/scalar/background intervals are retained
in the report. Write `Cmax=rho_hi ell delta sup|X|`, using the exact
pre-contact zero region. Then

```
Q_M=Cmax sup|psi_s|/3
       +(sup|E|+8sup|d|/3)sup|psi|
```

bounds the complete first correction throughout this history. Therefore

```
M2(T) in w_initial(T)+[-delta E_X-Q_M,delta E_X+Q_M].
```

For histories, the stronger varying bound
`|M2(x)-w_initial(x)|<=(x-1)E_X+Q_M` is integrated. Thus the weighted and
unweighted history errors are respectively

```
delta³ E_X/6+delta² Q_M/2,
delta² E_X/2+delta Q_M.
```

The free histories are enclosed by exact positive cell weights multiplying
the continuous new initial scalar boxes. There is no assumption that M2
is positive at every arbitrarily early contact time.

## Intrinsic volume and expansion without a poorly conditioned J trace

The exact timelike-shell Gauss identity with the retained pure-tension
junction rows is

```
R(gamma)=(D phi)²+F(phi),
F=2V-sigma'²/4+sigma²/3.
```

Put `v2=V2/V0` and use the intrinsic anchored proper-time slicing. The first
volume coefficient vanishes. The volume-weighted second Gauss identity is

```
v2_ss+4H v2_s-2v2/A² = F' M2/2+Q1,
v2(0)=v2_s(0)=0,
Q1=F''psi²/12-psi_s²/6+4psi²/(3A²)-5E²/(3A²)
     +E_s²-16E_s D_s/3+8D_s²/3.
```

`D` denotes the synchronous angular Hessian coefficient here. The harmonic
inputs suffice through `D_s=d_s-V/A-Cclock/A²`. Equivalently in longitudinal
variables, `P_s=u`, `E=u-HP`, `D_s=-P/A²`. We use the full harmonic expression
and absolute caps, so no unrecorded gauge cancellation is needed numerically.

One derivation uses the volume-weighted spatial scalar curvature. At a round
three-sphere, the functional `integral R3 dvol-(2/A²) integral dvol` is
critical. Its Hessian on first metric
`2(EY I+D HessY)` is `10E²/(3A²)`; the Hessian is zero on the pure spatial
diffeomorphism part. Thus the weighted second spatial scalar curvature is
`-4v2/A²+10E²/(3A²)`. The temporal term gives
`2v2_ss+8H v2_s+2<tr K1²-theta1²>`, with
`<tr K1²-theta1²>=-E_s²+16E_sD_s/3-8D_s²/3`. Combining this with
`<(Dphi)²>_2=-psi_s²/3+8psi²/(3A²)` gives the stated equation.
The full proper mean M2, rather than a coordinate scalar coefficient,
is required in the right side.

All coefficients are evaluated from authenticated shell boxes:

```
F'=2V'-sigma'sigma''/2+2sigma sigma'/3,
F''=2V''-(sigma''²+4sigma')/2
        +2(sigma'²+sigma sigma'')/3,
```

using the exact `sigma'''=4`. The absolute sum of every Q1 term is charged.

In scaled time, the exact positive Green kernel is

```
k(T,x)=sinh(2ell(T-x))/(2ell)
           *cosh²(ell x)/cosh²(ell T).
```

It follows from homogeneous solutions `tanh(ell x)` and
`1+tanh²(ell x)`. For `u=ell(T-x)`, `b=ell x`,
`cosh(u+b)>=cosh u cosh b`, hence
`0<=k/(T-x)<=tanh(u)/u<=1`. Therefore the first quadratic source history
has absolute Green integral at most `delta² sup|Q1|/2`.
Also `|k/(T-x)-1|<10^-11`: the lower ratio is at least
`1/cosh²(ell T)>1-2(ell T)²`, using `cosh t<=1+t²` for this tiny interval.
The reported history allowance charges this full relative error.

For expansion, `Theta2=v2_s` and

```
k_T=cosh²(ell x)/cosh²(ell T)
      [cosh(2u)-tanh(ell T)sinh(2u)].
```

On `1<=x<=T`, `|k_T-1|<10^-11` follows from the same cosh estimate,
`tanh(ell T)<=ell T` and `sinh(2u)<=2u cosh(2u)`.
The quadratic-source derivative history is charged with
`delta (1+10^-11) sup|Q1|`. Both volume and expansion calculations keep
the factor `(rho ell)^2`, and expansion then divides by `rho ell`.

## Scope and portable replay

The numerical method, grid, target and scientific endpoints are unchanged
from the separately frozen private feasibility campaign. Portable packaging
removes execution paths and elapsed timings from canonical artifacts and
registers the path-free source/code bindings before the portable final runs.
`CAMPAIGN_FREEZE.json` binds the N1024 inputs, T=5/4 cone, complete method and
three observables. `FROZEN_SCIENTIFIC_REFERENCE.json` fixes every original
scientific field. New response reports must match those fields exactly.

Use `producer.py --repo-root REPO --output-dir FRESH_EXTERNAL_DIRECTORY`.
The parent checkpoint authenticates the complete dependency closure before
scientific or third-party imports and replays normal and optimized modes.
Local response authentication separately checks all G/C/T/previous shell
budget source members and new response/initial/gauge/preparation inputs
before scientific imports and after calculation. The final reports retain
all actual first shell metric/scalar/background 2-jet intervals.

All signs are conditional formal second changes for the selected prepared
family at the same anchored proper clock. They do not establish negative
total volume, contraction, a radiation mechanism or a cosmological origin.
No effective nonzero parameter radius or full angular waveform is supplied.
