# Actual-background derivative bounds for the wider candidate

15 October 2026 UTC. New rational interval calculation confined to this
archive directory. The background remains the retained near-`phi=+1` regular
branch at exact `delta=1/1000` and
`c=96714321425692/161839258930731`. The new candidate has
`ell=1e-6`, `Z=8ell=8e-6`, and `T=3ell=3e-6`. This is a newly specified
perturbation width; no response coefficient from the previous `ell=1e-10`
state is transferred.

**Safe strict constants, now proved over the wider collar:**

\[
\boxed{|\mathcal W_z|<3\,10^7},\qquad
\boxed{|\mathcal W_t|<\frac{21}{400000}},\qquad
\boxed{|\beta_t|<\frac{21}{160000000}}.
\]

The looser proposed `|W_z|<1e8` is also valid. More useful near the frozen
constant `W(0,0)` is the pointwise variation bound

\[
\boxed{|\mathcal W(t,z)-\mathcal W(0,0)|
\le3\,10^7z+\frac{35}{4}t^2}.
\]

On the whole rectangle its absolute variation is strictly below
`192000000000063/800000000000 = 240.00000000007875`.
At `z=0`, the temporal variation alone is at most
`63/800000000000 = 7.875e-11`. These bounds are conditional on the retained
analytic branch theorem. No numerical trajectory or wave response is
validated by this background audit.

## 1. Wider-collar first-exit proof

The same actual background equations in inward `z` are

\[
\rho_z=-\rho^2a,\quad\eta_z=-\rho p,\quad
p_z=-\rho(U_\eta-4ap),\quad a_z=\rho^{-1}+\rho p^2/3.
\]

The authenticated prior shell interval lies strictly in

\[
\rho\in[129,131],\quad\eta\in[-9/10^5,-8/10^5],\quad
p\in[-14/10^5,-12/10^5],\quad a\in[111/1000,112/1000].
\]

The new program independently evaluates the vector field over that rectangle
and verifies all four strict inclusions
`I_shell + [0,Z] F(I_outer) subset interior(I_outer)` at the new
`Z=1/125000`. Thus no component can have a first exit before `Z`. The
vector field is smooth near the compact rectangle with positive `rho`, so
ODE continuation supplies the full interval. Exact lower and upper margins
are recorded in `BACKGROUND_DERIVATIVE_BOUNDS.json`.

The proper-distance change is bounded by `131 Z = 131/125000`, much smaller
than the retained proper shell endpoint lower bound `T_b>=9`. The collar
therefore stays in the actual regular branch; it is not an arbitrary
extension of coefficient functions. This background result supplies neither
the nonlinear constraint inverse nor its separate full support margin.

## 2. Derivatives of the actual potential

To distinguish the superpotential from the master potential, write its field
derivatives as `W_super''`, and put

\[
h=U_\eta/p,\qquad
Q=2h^2-10ah+35a^2/4-U_{\eta\eta}-p^2/6,
\]
\[
\mathcal W(t,z)=-35\tanh^2t/4+\rho^2Q.
\]

The division is valid throughout the first-exit rectangle. Direct
differentiation along its exact ODE gives

\[
h_z=-\rho(U_{\eta\eta}-h^2+4ah),
\]
\[
Q_z=(4h-10a)h_z+(-10h+35a/2)a_z
+\rho pU_{\eta\eta\eta}-pp_z/3,
\]
\[
\boxed{\mathcal W_z=-2\rho^3aQ+\rho^2Q_z},\qquad
\mathcal W_t=-35\tanh t\,\operatorname{sech}^2t/2.
\]

Exact interval evaluation of the factored spatial expression gives the
display-only range
`[-21392210.89480382,28190745.284941513]`; the authoritative exact rational
endpoints are in the JSON. This proves the strict `3e7` bound. The inequalities
`|tanh t|<=t` and `sech²t<=1` on `0<=t<=T` give the temporal bounds and
the variation estimate without numerical differentiation.

A separately obtained expanded expression agrees algebraically:

\[
\begin{split}
\mathcal W_z={}&\rho^3[4h^3-30ah^2+60a^2h-35a^3/2
+(12a-4h)U_{\eta\eta}\\
&\qquad +(29a/6-3h)p^2+pU_{\eta\eta\eta}]
+\rho(35a/2-10h).
\end{split}
\]

`check_background_derivative_identities.py` verifies this against the total
derivative and factored expression, as well as the ratio and shell identities.
It checks six formal identities and rejects four incorrect formulas. Those
checks are algebra, not a formal certification of the analytic branch/PDE.

## 3. Correlated shell inertia and tail gap

At the shell, all parameters are functions of the same `x=eta_b`:

\[
p=-2x-x^2-\delta c/2,\quad a=\sigma(x)/6,\quad
\rho^{-2}=\delta A(x)+\delta^2B(x),
\]
\[
A=W_{\rm super}f/9-W_{\rm super}'c/12,
\quad B=f^2/36-c^2/48,\quad f=1+c+cx.
\]

The new program derives and checks the exact Hamiltonian cancellation and
certificate coefficients. It partitions the retained `eta_b` interval into
128 rational subintervals, encloses each complete interval with rational
arithmetic and verified square-root brackets, then takes the hull. This is
interval coverage, not a set of sampled trajectories. No old producer is run.

Because the tension detuning is affine,
`sigma''/2=W_super''=2+2x`. The exact relations are

\[
B_b=h_b-4a_b+W_{\rm super}''(x),\qquad
\gamma=\rho_b(h_b-5a_b/2),\qquad
\mu^{-1}=\rho_bB_b,\qquad \beta_0=\gamma+5\mu/2.
\]

The important cancellation is

\[
\boxed{\mu^{-1}-\gamma
=\rho_b[W_{\rm super}''(x)-3a_b/2]},
\]
\[
\boxed{\mu^{-1}-\beta_0
=\rho_b[W_{\rm super}''(x)-3a_b/2]-5\mu/2}.
\]

The `h_b` dependence cancels exactly. Subtracting independent broad ranges
for `mu^-1` and `gamma` would discard this useful correlation. Compact strict
outward rational bounds are:

| Quantity | Strict lower bound | Strict upper bound |
|---|---:|---:|
| `mu` | `2059814/10^9` | `2271616/10^9` |
| `gamma` | `202092285/10^6` | `247360282/10^6` |
| `beta0` | `202097960/10^6` | `247365436/10^6` |
| `mu^-1-gamma` | `2381203/10^4` | `2381231/10^4` |
| `mu^-1-beta0` | `2381151/10^4` | `2381174/10^4` |
| `W(0,0)` | `283488885/10^4` | `691450116/10^4` |

Both tail gaps are positive. The time-dependent boundary coefficient obeys

\[
\beta(t)=\beta_0-35\mu\tanh^2t/4,\qquad
\beta_t=\mu\mathcal W_t,
\]

so `|beta(t)-beta0| <= 63/320000000000000` on the observation window,
using the independently verified `mu<1/400`. These identities provide
parameters for a future refined inertial-tail comparison; this audit makes
no assertion about its resulting response or postpulse sign.

## 4. Provenance and reproduction

Six public reference files were extracted byte for byte from immutable main
commit `c490f17139f1266426995595f5ba0eaefb14c636`, with one separately labeled
local planning reference. `DERIVATIVE_INPUT_PINS.json` records original Git
paths, sizes and SHA-256 hashes. The new arithmetic authenticates the retained
analytic proof, certificate and prior shell JSON before using them. No
historical producer, reference-only program or external model is executed.

Run only the new files in this directory:

```sh
python -I -B bound_background_derivatives.py
python -I -O -B bound_background_derivatives.py
python -I -B check_background_derivative_identities.py
```

The first two executions must produce identical rational receipt bytes.
The symbolic checker uses SymPy 1.14.0. Analytic branch, regularity and any
future exact-family support assumptions remain mathematical inputs; no
numerical waveform, effective nonlinear epsilon range or observational
cosmological scale is inferred here.
