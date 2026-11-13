# An exact nonlinear quiet window and delayed observer budgets

5 November 2026. Conditional application of retained finite propagation and
mixed-boundary uniqueness to the selected exact HDBLAST family. This note
evaluates an early zero interval and rational improvements of sufficient
observer bounds. It does not evaluate the nonlinear field norms, admissible
amplitude radius, or a finite-amplitude response after that interval.

## 1. Fixed source interface and hypotheses

Use the selected registered exact canonical curve, complete relative-harmonic
gauge and actual shell inclusion `z=0` in the common background-labelled
chart. Let `ell=10^-6`, `d=ell/2`. The first seed is supported in
`ell<=z<=2ell`; the **whole** exact canonical correction is supported in
the closure of `ell/2<z<5ell/2`. Thus the full exact canonical Cauchy data
agree with the background on the fixed open collar `0<=z<d`, for every
amplitude in the retained common admissible interval. The registered lapse,
shift and natural harmonic completion preserve that full open zero collar.
Canonical constraint equality alone would not fix the lapse and shift.

Retain the exact quasilinear shell theorem and its local domain-of-dependence
and uniqueness statements. Specifically require common actual radial
propagation speed bounded by `c_*=2`, positive auxiliary wave forms, regular
state-flow straightening and homogeneous mixed boundary conditions, and
existence on a common slab containing the closed interval `|t|<=ell/4`.
For an application to the later central target retain the longer
`T_coord=13ell/8` slab and all its original analytic premises. Neither that
slab nor `c_*=2` is certified by rational support arithmetic alone. The
amplitude interval on which they apply remains unevaluated.

For the past problem apply the theorem to the fully reflected canonical
data **and corresponding gauge**. Reflection is tensor pullback, including
the mixed metric components and gauge jets. It is not the rule of negating
every metric time derivative. Write `sigma=-t>=0` for the reflected
coordinate time, distinct from elapsed proper time.

## 2. Exact quiet wedge, including the physical boundary

The retained straightened reduced system has the form

```
g^{mu nu}(x,V) D_mu D_nu V = F(x,V,DV),   F(x,0,0)=0,
V_D|_{z=0}=0,                 D_z V_N|_{z=0}=0.
```

The perturbation and its complete Cauchy velocity vanish at `sigma=0` on
`0<=z<d`. Compare the actual solution with `V=0` on the shrinking collar

```
W = { sigma>=0, z>=0, z+2sigma<d } x S3.
```

For any smaller compact truncated collar before `sigma=d/2`, the difference
equation has bounded local coefficients. Since `F(x,0,0)=0`, its lower-order
source is bounded by a constant times the difference and its first jets.
The homogeneous Dirichlet/Neumann shell rows make the fixed-shell energy
flux vanish. The sloping boundary moves inward at the common permitted
radial propagation speed; the retained local energy/domain-of-dependence
estimate prevents incoming difference data through it. Starting from zero
Cauchy difference, the localized energy inequality and Gronwall give zero
difference on the truncated collar. Exhaustion yields

```
V_e=0 on W, for every admitted e.
```

This uses the retained actual-speed bound and mixed-boundary uniqueness,
not only vanishing formal time jets or the linear packet's propagation.
The shell is not an independently driven boundary: its physical laws are
already encoded by the homogeneous difference rows.

The boundary-straightening cutoff width `h` is not a shorter physical
quiet-window distance. The change of unknown is an invertible point map
`U=Psi(x,V)` fixing zero. Its boundary state vector field vanishes at
`U=0`, so `Psi(x,0)=0` even in the cutoff transition. All x-derivatives
of the identically zero function `x -> Psi(x,0)` vanish. This says nothing
about V-derivatives of Psi, whose invertibility is retained. Consequently
the transformation does not source the background; the retained transformed
equation explicitly has `F(x,0,0)=0`. An artificial outer construction
does not enter this local wedge under the retained propagation margin.

Returning through this zero-fixing map and the actual fixed-shell inclusion
therefore gives the **full physical fields**, not merely a selected
projection or first Taylor coefficient:

```
h_e(t,n)=rho_b^2[-dt^2+cosh(t)^2 kappa_S3],
varphi_e(t,n)=phi_b,                       |t|<ell/4.
```

For each admitted amplitude the full fields are continuous, with the
retained stronger C4 regularity. Since the common solution domain includes
the endpoint, their shell restrictions also equal the background at
`|t|=ell/4` by a one-sided limit. This is a closed-endpoint equality; it
does not supply an open background spacetime neighborhood beyond the
causal apex. No claim of first arrival **at** `ell/4` follows. The signal
may remain zero longer.

## 3. Proper observers, cut area and moment in the window

Use the actual intrinsic normal geodesics anchored to the fixed common
initial physical section `(t,n)=(0,n)`, with its actual future unit normal
and transported round labels. On this interval the shell metric is exactly
the background, so geodesic uniqueness gives

```
q=-tau/rho_b,   x_e(q,n)=(-q,n),
u_e=partial_tau x_e=rho_b^-1 partial_t,    0<=q<=q_D=ell/4.
```

The constant-angle curves are geodesics of the displayed de Sitter shell
metric and have unit future proper velocity. The endpoint statement follows
from the continuous solution at the endpoint, without requiring an open
quiet neighborhood there. The normal flow is noncrossing and the cuts are
spacelike on this early interval by their exact background form. This
does not prove noncrossing after the interval or ambient physical shell
embedding outside the supplied fixed inclusion.

The entire augmented state `(x,u,P=D_n x,Q=D_n u)` is amplitude independent
on the interval, in any fixed amplitude-independent representation. Thus
at the **closed** reset time

```
S_D=rho_b ell/4,   X_e(S_D)=X_0(S_D),
v_disp,D=0,             V_D=partial_e X_e(S_D)=0,
W_D=partial_e^2 X_e(S_D)=0.
```

The continuation retains the original physical anchor, transported labels,
template and proper-time origin. Writing an active integration clock
`s-S_D` is only an internal translation of the same initial-value problem;
it does not redefine the observable or move its original central proper-time
bracket. At the reset cut the inherited velocity equals its actual future
unit normal because both the cut and flow are exactly background there.

In intrinsic product coordinates Q is zero. In an auxiliary embedding its
baseline may be nonzero, but still has zero parameter derivatives. A
parameter-dependent physical identification, anchor, representation or a
substituted covariant variable requires its own costs. The later algebraic
slicing normal is not being equated to transported geodesic observers.

The actual cut metric and measure are exactly

```
k_e(q)=rho_b^2 cosh(q)^2 kappa_S3,
dA_e(q)=rho_b^3 cosh(q)^3 dOmega.
```

With the registered transported template `Y=n1^2-n2^2`,
`integral_S3 Y^2 dOmega=pi^2/3`. Hence the actual-area denominator is
`rho_b^3 cosh(q)^3 pi^2/3>0`, the background-subtracted scalar numerator
vanishes, and the specified actual-area quadrupole is exactly

```
M_e(-rho_b q)=0,                   0<=q<=ell/4.
```

All its admitted parameter derivatives through order two vanish there.
This is the protected early zero interval; it is not the later central
sign-changing nonlinear zero sought by the remainder campaign.

## 4. Shorter active ODE duration with the same coefficient inputs

The reset state is known exactly. Start the retained coefficient-neighborhood
continuation/variation argument there, rather than propagating zero initial
variations through an unnecessarily long interval. For the target
`T_q=25ell/16`, the active durations are

```
Delta_q=T_q-q_D=21ell/16,
Delta_S=rho_b Delta_q,    Delta_S/(rho_b T_q)=21/25.
```

Require the same full-ball regular coefficient hypotheses and state scales
as the retained bootstrap, now on the post-reset background trajectories.
The quiet endpoint does not make `F_e`, `F_Xe` or `F_ee` zero on arbitrary
coefficient balls: their position components can extend beyond the quiet
window. No such off-trajectory zero is needed for the reset argument.
The actual post-reset coefficient radius and bounds remain to be supplied.

For uniform nonnegative `L,P1,H,C,P2` on that neighborhood, the sufficient
nonexit and parameter bounds become

```
a phi_L(Delta_S) P1 < r,
V_* = phi_L(Delta_S) P1,
W_* = phi_L(Delta_S)[H V_*^2+2 C V_*+P2],
phi_L(S)=(exp(LS)-1)/L,    phi_0(S)=S.
```

These use proper duration `Delta_S`; coordinate-time and proper-time
coefficient normalizations must not be mixed. In a finite chart chain all
transition derivatives and Hessians must still be charged. In a fixed
auxiliary embedding the derivative bounds must include its extension costs.

For a comparison using **identical** coefficient bounds before and after
the reset, put `S=rho_b T_q`, `k=21/25`. The integral representation gives

```
phi_L(kS)=k integral_0^S exp(L k r) dr
        <=k integral_0^S exp(L r) dr=k phi_L(S),  L>=0.
```

The earlier zero-data V majorant improves by a factor at most `21/25`.
Writing W as the sum
`H P1^2 phi^3 +2 C P1 phi^2 +P2 phi`, the respective termwise factors
are at most `9261/15625`, `441/625`, and `21/25`. These are inequalities
between nonnegative majorants, including their zero cases; no division
by a zero old bound is intended. The H-only term shrinks by at least
40.7296 percent, but this is a bookkeeping improvement under the same
inputs, not a measured nonlinear error or an evaluated amplitude increase.

## 5. A less restrictive sufficient clock-drift guard

Let `U1` be the uniform first-amplitude bound on the projected future
proper-velocity component in the same admitted active tube. For
`|e|<=a`, set `c=rho_b a U1`. The fundamental theorem in amplitude and
the background velocity give `|t_q+1|<=c`. Because `t+q=0` exactly up
to `q_D`, for `q_D<=q<=T_q` one now has

```
|t_e(q,n)+q| <= c(q-q_D) <= c*21ell/16.
```

Thus `c<1/42`, instead of the previous whole-duration `c<1/50`, is
sufficient for drift `<ell/32`. With `T_coord=13ell/8`,

```
T_coord-T_q=ell/16,
T_coord - sup(-t_e) > ell/16-ell/32=ell/32,
-t_q >=1-c >41/42 on the active interval,
-t_q=1 on the quiet interval.
```

The upper past-slab margin remains strict. Time monotonicity also keeps
the flow in the past of the anchor. This clock estimate is used only after
the flow/velocity bounds have been justified on the common domain; it must
not be used circularly to assume that domain. The threshold is `25/21`
times the old sufficient c threshold. No corresponding amplitude increase
is proved because `U1`, the field bounds, admissible neighborhood and radius
are unevaluated and may depend on a.

## 6. What can and cannot now be evaluated

This note evaluates the full-family quiet duration and exact early observer,
area, denominator and moment values, and improves the duration and clock
thresholds in existing conditional bounds. It does **not** evaluate the
post-reset `L,P1,H,C,P2` or r. The exact-family theorem supplies C2 parameter
dependence into a lower shell C4 topology, not numerical mixed-derivative
suprema or an explicit admissible amplitude. Its energy majorant constants,
state neighborhood, weighted inverse norms and boundary state-flow interval
remain qualitative. The selected actual first field and selected second
projections are derivatives at e=0; they cannot substitute for uniform
bounds over `|e|<=a` and every off-trajectory coefficient ball. Full scalar,
metric, physical embedding, area and denominator remainder estimates after
the quiet interval are still required.

This is an application of established local domain-of-dependence,
uniqueness, geodesic transport, continuity and Gronwall methods. It neither
proves the retained nonlinear PDE premises anew nor establishes a
higher-dimensional origin of the Big Bang.

## Sources

The immutable source identities and passage ranges are recorded in
`SOURCE_PINS.json`. All links below refer to the actual PR125 merge.

- [Registered exact-family binding](https://github.com/maldonado-research/HDblast/blob/eec64b96c751086d8c914ce66ddf613a5e718968/research/HDBLAST_NONLINEAR_OBSERVABLE_BUDGETS_20261105/inputs/REGISTERED_EXACT_FAMILY_BINDING.md), sections 3–4: whole correction support, actual-speed choice, fixed-shell chart and complete initial gauge.
- [Exact quasilinear shell family](https://github.com/maldonado-research/HDblast/blob/eec64b96c751086d8c914ce66ddf613a5e718968/research/HDBLAST_NONLINEAR_OBSERVABLE_BUDGETS_20261105/inputs/EXACT_QUASILINEAR_FAMILY.md), sections 1–4: common family and C2/C4 regularity, zero-fixing mixed-boundary reduction, localized energy propagation and parameter dependence.
- [Nonlinear bridge](https://github.com/maldonado-research/HDblast/blob/eec64b96c751086d8c914ce66ddf613a5e718968/research/HDBLAST_NONLINEAR_OBSERVABLE_BUDGETS_20261105/inputs/NONLINEAR_BRIDGE.md), sections 2–3: reflected gauge, longer coordinate slab, physical normal observers and actual-area quadrupole.
- [Observer-flow budget](https://github.com/maldonado-research/HDblast/blob/eec64b96c751086d8c914ce66ddf613a5e718968/research/HDBLAST_NONLINEAR_OBSERVABLE_BUDGETS_20261105/observer/OBSERVER_FLOW_BUDGET.md), sections 1–4: full physical pullback, augmented state, proper-clock scaling, chart costs and the existing majorants.
- [Coefficient bootstrap](https://github.com/maldonado-research/HDblast/blob/eec64b96c751086d8c914ce66ddf613a5e718968/research/HDBLAST_OBSERVER_DOMAIN_20261105/bootstrap/FLOW_BOOTSTRAP.md), sections 1–5: full coefficient-ball assumptions, first exit, delayed initial-value reset and clock/slab bookkeeping.
