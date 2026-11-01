# Method, proof obligations, and scope

## Geometric quantities and parameter convention

Write `t=ell*x`, `ell=1/1000000`, `tau=rho_b*t`,
`A=rho_b*cosh(t)`, and `H=tanh(t)/rho_b`. A prime in the formulas below is a
proper-time derivative `d/dtau`. Every suffix `2` denotes an actual second
parameter derivative: `F_e=F_0+e*F_1+e^2*F_2/2+o(e^2)`.
The probability-normalized round `S3` average is used in the first-only
quadratic formulas. The geometric average in the curvature defect uses the
perturbed volume measure; its second monopole contribution cancels exactly
against the volume-radius term.

For the actual longitudinal shell master `u(t,0)`, define

```
P(tau) = integral_0^tau u(v) dv
E = u-H*P
E_prime = u_prime-H*u-H_prime*P
H_prime = 1/A^2
L = P/A^2
theta_1 = (3*E_prime+8*L)*Y
shear_1 = -L*(Hess(Y)+(8/3)*Y*kappa)
```

Here `Y=n1^2-n2^2`, with `<Y^2>=1/6` and `Delta(Y)=-8Y`. For the anchored
normal proper-time foliation, the exact first-only second coefficients are

```
C2 = (<R3> - 6/R_volume^2)'' = (10/3)*E^2/A^2
S2 = (<|shear|^2>)'' = (80/9)*L^2
Var2 = (Var(theta))'' = (3*E_prime+8*L)^2/3
Q_B2 = ((2/3)*Var(theta)-<|shear|^2>)''
     = (2/9)*[(3*E_prime+8*L)^2-40*L^2]
```

The separate geometry component proves the anchored Gaussian extraction,
angular moments, curvature cancellation, and intrinsic shell trace equation.
The numerical aggregation and clock conversion were independently checked
against those identities before this portable component was assembled.

## Actual trace recovery

The linearized intrinsic four-dimensional shell trace equation is

```
u_xx + a(x)*u_x + b(x)*u = f(x)
a(x) = 3*ell*tanh(ell*x)
b(x) = ell^2*[8*sech(ell*x)^2-4]
f(x) = ell^2*rho_b^2*F_prime*psi(x,0)/6
u(0)=u_x(0)=0
F_prime = 2*p_b*B_cal
        = 2*V_prime-8*a_b*p_b+sigma_second*p_b < 0
```

The exact first initial shell collar supplies the zero initial master data.
The physical fixed shell has `f_gauge|_b=0`; since its background scalar
pullback is constant, the public actual scalar jet `s` is exactly the
longitudinal shell scalar amplitude `psi`. Tangential gauge labels do not
alter it. The two expressions for `F_prime` are enclosed using authenticated
actual shell boxes; their intersection retains every actual coefficient.
The coefficient is certified strictly negative before history integration.

Every history source value is the public `actual_jet(...,mode='actual')`
scalar enclosure, with the complete inherited actual-minus-comparison C3
errors. Error is charged on the whole interval from zero, including `[0,1]`.
No exact precontact zero or additional finite-propagation premise is used.
No free comparison is substituted for the actual source.

Let `T=63/50`, and let `a_cap,b_cap` bound the absolute values of the two
lower-order coefficients on all 64 whole history cells. Set

```
F_int = sum(cell_length*sup_cell|f|)
den = 1-a_cap*T-b_cap*T^2/2 > 0
V_cap = F_int/den
U_cap = T*V_cap
G = a_cap*V_cap+b_cap*U_cap
```

From `|u(s)|<=s*sup|u_x|`, the integrated equation gives
`sup|u_x|<=F_int+(a_cap*T+b_cap*T^2/2)*sup|u_x|`.
Consequently `G` bounds the omitted feedback acceleration, and

```
u_x(x) = integral_0^x f(s) ds                    +/- x*G
u(x)   = integral_0^x (x-s)*f(s) ds              +/- x^2*G/2
integral_0^x u(v) dv = integral_0^x (x-s)^2*f(s)/2 ds +/- x^3*G/6
P = rho_b*ell*integral_0^x u(v) dv
u_prime = u_x/(rho_b*ell)
```

These are global rational bounds, not time-stepping approximations. For
kernel powers `p=0,1,2`, each source-cell weight is computed exactly as
`integral_cell_intersection[0,x] (x-s)^p ds`. It is nonnegative and monotone
in `x`; its values at the two target-cell endpoints enclose the whole target
cell. Multiplying its interval by the full source-cell enclosure remains
valid if the source changes sign or varies within the cell. No point
quadrature or endpoint sampling is used.

All decisions use rational outward interval arithmetic on the inherited
`2^-112` grid. Very small squared quantities can have a zero rounded lower
endpoint. The full-interval `P` enclosure also contains zero. Accordingly the
recorded mean squared shear lower bound remains zero; no stricter numerical
claim is made for it. Recovered `u`, `E`, curvature defect, expansion variance,
and Buchert coefficient exclude zero throughout the target interval.

## Source closure and portability

`SOURCE_PINS.json` pins the complete manifests and every selected member of
`HDBLAST_GEOMETRIC_BULK_SOURCES_20261016` and
`HDBLAST_EXACT_FAMILY_TT_TRANSFER_20261031`. The geometric capsule supplies
all executed numerical dependencies; the exact-family capsule retains the
applicability premise. Full closed inventories, byte counts, hashes,
non-symlink files, and safe unique member paths authenticate before the first
scientific import and after the last call. Essential local pins and preserved
protocol files also authenticate before and after. The running producer hash
is checked before and after the computation. No installed project package,
private producer, absolute origin path, mutable upstream producer, or network
resource is a dependency.

The frozen source campaign is 64 history cells covering `[0,63/50]`, eight
cells covering `[31/25,63/50]`, and the point `5/4`, all at physical `r=0`.
The complete metric, scalar, background, and global-master jets are retained
for all 73 calls. The original global-master enclosure is checked for
compatibility with each recovered actual master value and derivative; it is
not intersected to manufacture a sign. The scientific projection digest locks
the already frozen deterministic witness and supports strict mutation
controls. It is a reproducibility gate; the geometry proof and independent
raw-row reconstruction establish the mathematical interpretation.

## Preserved protocol and numerical limits

The protocol records both limitations encountered before the successful
trace campaign. The original method used exact shell cancellation of the
comparison master plus the full inherited global master error. Its symmetric
actual master boxes contained zero, so its curvature lower bound was zero
and its Buchert enclosure straddled zero. Its frozen 64-cell final report is
retained as a method limitation reference. A separate eight-cell trace pilot
also contained zero because whole-cell source bounds were too broad. A
separate 64-cell trace assessment reused the preserved raw actual source rows
only to choose the new method's freeze. A fresh 73-call final campaign then
recomputed every actual source row on the frozen 64/8 partition. Portability
introduced no refinement, new target, changed error budget, or changed
mathematical endpoint. Historical origin strings were removed from portable
protocol copies; original and normalized SHA256 receipts are preserved.

`Q_B2` is the Buchert kinematic coefficient. It is distinct from the older
Lorentzian shell scalar-gradient mean, often also called `Q2`. Positivity of
this kinematic coefficient and of the curvature defect does not imply a
finite-amplitude backreaction effect, a positive acceleration, a force,
radiation, an SI magnitude, hot matter production, or a cosmological origin.
The actual-family and first-field premises remain conditional; no effective
nonzero parameter radius or second-order remainder bound is supplied.
