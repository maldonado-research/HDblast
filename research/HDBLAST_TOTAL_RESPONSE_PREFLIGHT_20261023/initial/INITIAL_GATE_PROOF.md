# Bounded registered five-mean initial/free-wave gate

23 October 2026. **FAILED PROVISIONAL ACCURACY GATE; VALID COARSE ENCLOSURE.**
This work constructs one explicit differentiable approximation to the registered
physical initial data and certifies its free-wave error on the entire retained
cone. It stops before bulk iteration because the resulting certified bound does
not meet the provisional initial allocation `10^28`. This is not a proof that
an approximation of the actual solution cannot attain that accuracy.

## Authentication and scope

Before any numerical capsule import, `authenticate_inputs.py` checked the full
manifest closures P (107 files), N (111), G (168), and C (48). The P/N manifest
hashes were independently compared with the existing publication pin files;
G/C hashes were compared with the frozen forward-plan pins. Every declared
file's length and SHA256 and the absence of undeclared non-cache files were
checked. The identical closure check ran after the evaluation. `AUTH_BEFORE.json`
records every verified entry. All old capsules and reports remained read-only.

Only the initial footprint `[0,7/4]` was numerically covered: 96 whole cells of
width `1/128` on `[1,7/4]`, together with the exact zero interval `[0,1]`.
There was no n16/n32 source integration, bulk iteration, or remote write.

## Registered initial prescription

Mean order is `(K0,L0,N0,J0,w0)` with Gram diagonal `(1,2,1,3,1)`.
Let `f0` be the mean Hamiltonian right inverse, `K=Fmom/3`,
`Pi=(Gmom-Fmom_y)/p`, and `zeta0=-3 f0_y/(2p)`. Then

```
U0 = (0, 0, -u^2/3, f0, zeta0),
V0 = ell * (
 (2 rho K - 8 u ut)/6,
 rho*((-26 u uy - 2 a u^2)/6 - 3 f0_y/2 - 3 a f0),
 (-8 rho K + 44 u ut)/6,
 (2 rho K + 4 u ut)/6,
 (rho Pi + 2 u vt)/6
).
```

Here `ell=10^-6`, `dy=-rho^-1 dz`, and `partial_r=ell partial_z`.
The lapse-square initial value, mixed `y0` velocity, second scalar value,
and scalar velocity are all retained. `initial_frame_data` and `pack_initial`
produce these data. The two momentum combinations are evaluated in their
algebraically reduced forms: the producer supplies `Fmom=3K` and evaluates
`Gmom-Fmom_y=p Pi` directly. The artificial arguments `Fmom_y=0,Gmom=pPi`
are only a way to evaluate that difference through the retained adapter;
they do not alter the individually prescribed physical momentum profiles.
The independent supplementary exact checker proves the reduced `K`, `Pi`,
and Hamiltonian expressions equal the physical prescription.

## A single rationally specified differentiable representation

The nominal background functions `rhohat,etahat,phat,ahat,alphahat` are fixed
rational polynomials of degree four in `z=ell*r`. Their coefficients are
selected once from shell midpoints and the static ODE jet recurrence; their
complete rational coefficient lists appear in `INITIAL_GATE_NORMAL.json`.
These polynomials define functions, not separately chosen cell values.
They need not solve the background equations exactly; their approximation
error is included in the interval comparison.

Set `n=2^40`, `d=(r-1)(2-r)`, and

```
Qbase(r) = [1 - (1-4d)/(1+(n-4)d)]^n,
S(t) = 126 t^5 - 420 t^6 + 540 t^7 - 315 t^8 + 70 t^9.
```

Use `Qhat=0` for `r<=65/64`, `Qhat=Qbase*S(64r-65)` on
`[65/64,66/64]`, and `Qhat=Qbase` above `66/64` on this footprint.
Use `Ehat=(1+z/n)^(-n)`. All these formulas have rational coefficients.
All required denominators were enclosed away from zero during evaluation.
The power `n` is evaluated by forty squarings, without expanding its
polynomial coefficients.

`S(0)=0`, `S(1)=1`, and derivatives 1 through 4 vanish at both endpoints.
Consequently `Qhat` is C4 across both joins. Using the same rational background
functions, packet, and their differentiated jets in every expression gives
one C1 `Uhat0` and one continuous `Vhat0`; the scalar value involves at most
three packet derivatives and its radial derivative at most four. Their exact
zero extension near `r=0` makes the mixed odd/even reflection compatible.
`Uhat0_L=0` identically. The actual packet on `[1,65/64]` is not discarded:
its actual analytic flat-tail enclosures are included in the error ledger.

Define the candidate free wave by the retained d'Alembert formula using
`Uhat0,Vhat0`, with odd extension for L0 and even extension for the other
four means. This integral of explicitly specified rational functions defines
one continuous C1 candidate on D+. The candidate and its first derivatives
vanish on the exact flat wedge `x+r<=1`. The Dirichlet L0 boundary value
is exactly zero. No Robin completion or claim that this free wave equals
the total response is made.

## Whole-cell and whole-cone error proof

The actual background enclosures are N's retained shell-plus-ODE integral
bounds with differentiated ODE recurrences through order four. Actual packet
jets use the retained derivative polynomials and rational exponential
bounds, including the analytic endpoint tail. The nominal derivatives are
obtained by differentiating the single rational representation above.
On every closed cell, interval arithmetic encloses both actual and nominal
jets. Subtracting those enclosing intervals bounds their difference.
Independent interval subtraction can overestimate error; it remains a
valid bound. The exact cell bounds, component maxima, and outward square-root
bounds appear in the report. No interval midpoint differentiations are used.

The 112-bit outward rational grid is inherited from N's arithmetic. For
component caps `c_i`, the H error is bounded by
`sqrt(c_K^2+2 c_L^2+c_N^2+3 c_J^2+c_w^2)`, using a rational upper square-root
bracket. The resulting exact rational maxima have display values

| Error | Certified upper bound (display only) |
|---|---:|
| delta0, initial H value error | 5.482e31 |
| deltar, scaled radial H error | 4.561e33 |
| deltav, scaled time H error | 1.313e24 |

For `0<=x<=3/2`, reflection is an H isometry. The free-wave value error is
at most `delta0+x*deltav`; either first derivative error is at most
`deltar+deltav`. Therefore on every point and characteristic segment of D+,

```
||U_free - Uhat_free||_X
 <= max((delta0+(3/2)*deltav)/3, deltar+deltav)
 < 4.561e33                            (safe outward display).
```

This one bound includes inherited background uncertainty, background
polynomial and packet/cutoff approximation, velocity and derivative
conversion, whole-cell variation, and outward arithmetic exactly once.
There is no additional quadrature defect because the candidate is defined
by the exact d'Alembert integral of the rational representation. Numerical
evaluation of that integral would require its own integration bound.

## The failed gate and what it establishes

The initial allocation was `10^28`; the certified bound above is much larger.
Even at zero spatial cell width, the retained actual scalar-mean oracle at
`r=5/4` returns intervals with approximate endpoints

```
zeta0 in [-7.66420312853039e30, -5.711425298940109e30],
partial_r zeta0 in [7.312541126620708e31, 1.1170638762444665e32].
```

The corresponding interval radii are about `9.76388914795141e29` and
`1.9290488179119788e31`. No radius-`10^28` interval around any single value
can contain either entire retained interval. This proves the current
whole-enclosure verification route fails that pointwise allocation even
before spatial discretization error is considered. These radii are **not**
lower bounds on the actual approximation error: the true fixed background
may occupy a much smaller subset. This report neither excludes a tighter
correlated background proof nor proves mathematical impossibility. Finer
spatial cells alone do not shrink the zero-width input/oracle intervals.

The two-term total-response route requires a certified `e0<=2e28` and
`e1<=2e28`. This pilot has not certified the first allocation, and it has not
constructed B1. The retained iteration-tail estimate cannot repair that
missing premise. No evaluated total-response enclosure is claimed here.
The concrete next requirement is a tighter correlated initial/background
certificate or a redesigned initial error proof, followed by a new bounded
preflight; a larger source grid does not address this failed gate.

## Reproduction and limits

Run the producer under Python standard library only:

```
python -I -S -B certify_initial_representation.py --output fresh.json
python -I -S -B -O certify_initial_representation.py --output fresh-o.json
python -I -S -B check_initial_gate_ledger.py --report fresh.json --output ledger.json
```

For public replay, pass `--dependency-map deps.json`, where `deps.json` is
a JSON object mapping `P`, `N`, `G`, `C` to the four extracted capsule-root
paths. This mode authenticates the complete closures against four expected
manifest SHA256 values embedded outside those capsules in the new authenticator.
It does not require private planning or publication pin files. The default
mode uses the recorded research-work paths and additionally cross-checks their
publication/forward-plan external pin files. Numerical formulas are identical.
Dependency paths and pin-source descriptions in the public-mode report differ;
the exact numerical ledger must match the retained report. `check_initial_gate_ledger.py`
is portable and checks the exact report arithmetic, coverage, H weights,
free-wave estimate and honest failed-gate scope; it does not replay the
numerical enclosure calculation or reauthenticate absent dependency trees.

This is a five-mean initial/free-wave certificate. It is not an evaluated
full-angular field, exact numerical constraint identity, second-derivative
Einstein residual bound, invariant shell observable, finite-amplitude family,
or physical sign inference. The retained exact conditional subsidiary
completion applies to the exact total solution, not to this approximation.
