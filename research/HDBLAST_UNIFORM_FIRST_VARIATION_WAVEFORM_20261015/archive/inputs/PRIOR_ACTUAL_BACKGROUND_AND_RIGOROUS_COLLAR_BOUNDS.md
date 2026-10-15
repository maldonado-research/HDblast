# Actual regular-branch input and a rational collar enclosure

15 October 2026 UTC. Focused source audit and new exact-arithmetic deduction.
This note treats historical documents as research references. No archived
program was imported or executed, and no repository was modified.

**Result.** The retained analytic branch theorem supplies an actual background
at the exact detuning `delta=1/1000`. Its existing error bounds are sufficient
to enclose the shell and a very short inward collar without validating a
historical floating trajectory. The new rational calculation proves, for
`0 <= z <= 1/1250000000 = 8e-10`,

\[
\frac1{1000}<\mu<\frac1{400},\qquad
|\mathcal W(t,z)|<10^6,\qquad
|\beta_{\rm bdry}(t)|<1000
\]

for every real `t`, together with

\[
|C_b/2|<250,\quad |C(z)/2|<500,\quad
16641<\rho_b^2<17161,\quad
\frac9{625000000}<p_b^2<\frac{49}{2500000000}.
\]

These are deliberately coarse, conditional consequences of the retained
analytic theorem. They are not a computed background trajectory, a waveform,
an effective nonlinear amplitude radius, or a formal proof certification.

## 1. Fixed model and actual background being enclosed

Use the exact-decimal model and the near-`phi=+1` regular branch. With
`eta=phi-1`,

\[
W=\frac13+\eta^2+\frac{\eta^3}{3},\quad
U=\frac12W_\eta^2-\frac23W^2,\quad
\sigma=2W+\delta(1+c+c\eta),\quad
c=\frac{96714321425692}{161839258930731},\quad\delta=\frac1{1000}.
\]

The exact rational `c` equals `2/1.0357712571566784-4/3` when the decimal
is interpreted exactly. No equivalence to its binary64 evaluation is needed.
The full background equations, Hamiltonian constraint and shell laws are

\[
\rho_y=a\rho,\quad \eta_y=p,\quad
p_y=U_\eta-4ap,\quad a_y=-\rho^{-2}-p^2/3,
\quad a^2=\rho^{-2}+p^2/12-U/6,
\]
\[
\rho(0)=0,\ \rho_y(0)=1,\ p(0)=0,\quad
a_b=\sigma_b/6,\quad p_b=-\sigma_{\eta,b}/2.
\]

The October 11 theorem and exact rational certificate establish one regular
solution for every `0<delta<=1/500`, unique within their specified weighted
Dirichlet tube and with proper radial endpoint `T_b>=9`. They do not claim
global uniqueness among other branches. The branch has `|eta|<=1/100`,
`a>0`, `p<0` on the open interior. We use that actual solution at fixed
`delta=1/1000`, defined by the theorem, rather than defining a new background
by rounded historical endpoint values.

The theorem's sharper outward bounds (its Eq. (29)) give, with
`b_eta=eta_b`,

\[
\left|b_\eta+\frac{9c}{64}\delta\right|
<\frac{799}{250}\delta^2,
\]
\[
\left|\rho_b^{-2}-\frac{1+c}{27}\delta-
\left(\frac{(1+c)^2}{36}-\frac{c^2}{384}\right)\delta^2\right|
<\frac{127}{500}\delta^3.
\]

The separately retained certificate gives even smaller exact rational
constants. The new program checks those are below the outward constants
used here. It then uses rational interval operations; 80-bit dyadic square
roots are bracketed by explicitly checking their squared endpoints.

Resulting decimal displays, not the decision arithmetic, are

| Quantity | Outward exact-rational enclosure, displayed approximately |
|---|---|
| `eta_b` | `[-8.72327877384e-5, -8.08407877383e-5]` |
| `rho_b` | `[129.9244817666, 129.9250388364]` |
| `p_b` | `[-1.3712350160e-4, -1.2433842727e-4]` |
| `a_b` | `[.111377370423, .111377371419]` |
| `B_b` | `[3.38817111197, 3.73669937891]` |
| `mu=(rho_b B_b)^-1` | `[.00205977121593, .00227166178147]` |

The exact numerator/denominator endpoints are authoritative in
`ANALYTIC_SHELL_BOX_DELTA_1E3.json`. The table is only a short display.

In particular, because the detuning is affine,

\[
\boxed{\sigma''_b=4+4\eta_b},\qquad
B_b=\frac{U_\eta(\eta_b)}{p_b}-4a_b+2+2\eta_b>0.
\]

The small `sigma''` in the old stable S85 certificate does not belong to
this branch and boundary theory.

## 2. A strict first-exit enclosure of the inward collar

The current response uses inward `z=log(b/r)` and `dy=-rho dz`. The exact
background ODE becomes

\[
\rho_z=-\rho^2a,\quad \eta_z=-\rho p,\quad
p_z=-\rho(U_\eta-4ap),\quad a_z=\rho^{-1}+\rho p^2/3.
\tag{A}
\]

Choose the rational outer rectangle

\[
\rho\in[129,131],\quad
\eta\in[-9/10^5,-8/10^5],\quad
p\in[-14/10^5,-12/10^5],\quad
a\in[111/1000,112/1000].
\]

It strictly contains the shell box. Evaluate (A) by rational interval
arithmetic over that outer rectangle. For every component `i`, the program
checks the strict inclusion

\[
I_{i,\rm shell}+[0,L]F_i(I_{\rm outer})
\subset\operatorname{int} I_{i,\rm outer},\qquad
L=\frac1{1250000000}.
\tag{B}
\]

All four lower and upper margins are positive exact rationals recorded in
the receipt. If the actual solution had a first exit before `L`, integrating
its derivative up to that exit would put it strictly inside the same box by
(B), a contradiction. The vector field is smooth near this compact rectangle
because `rho>0`, so the usual finite-dimensional ODE continuation theorem
excludes a failure of existence while it remains there. This proves (A)
stays in the rectangle throughout the stated collar. It uses neither an
ODE integrator nor unverified sampled background values.

## 3. Coefficient bounds without numerical differentiation

Put `h=U_eta/p`. This division is safe throughout the outer rectangle.
Then the actual coefficients satisfy

\[
C=\rho(3a+2p_y/p)=\rho(2h-5a),
\]
\[
C_z=\rho^2[2h^2-10ah+5a^2-2U_{\eta\eta}-5p^2/3]-5.
\tag{C}
\]

The conjugated physical response master has

\[
\mathcal W=\frac{35}{4}\operatorname{sech}^2t-\frac{25}{4}
+\frac{C_z}{2}+\frac{C^2}{4}+\frac23\rho^2p^2.
\]

Substituting (C) eliminates all differentiated background coefficients:

\[
\boxed{\mathcal W(t,z)=-\frac{35}{4}\tanh^2t+
\rho^2[2h^2-10ah+\tfrac{35}{4}a^2-U_{\eta\eta}-p^2/6].}
\tag{D}
\]

The physical even boundary has `-w_z+mu w_tt+beta_bdry(t)w=0`, where

\[
\mu=(\rho_bB_b)^{-1},\qquad
\beta_{\rm bdry}(t)=\mu[\tfrac{35}{4}\operatorname{sech}^2t
-\tfrac{25}{4}]+C_b/2.
\]

The ranges `0<=tanh²t<=1` and
`-25/4<=35 sech²t/4-25/4<=5/2` imply the displayed simple bounds by exact
rational interval evaluation. The JSON includes the full coefficient boxes,
not just the coarse constants. A separate newly executed formal SymPy check
verified (C) by total differentiation along (A), and then verified (D);
three incorrect coefficient/sign variants were nonzero. That formal check
is algebra only, not the analytic proof of (A)–(B).

## 4. What historical numerical data exist

The September 22 `static_branch/solve_plus_branch.py`, result JSON and
`PLUS_BRANCH_PROFILES.npz` concern the same near-`phi=+1` branch and affine
detuning. At `delta=.001`, their refined floating row reports

\[
T_b\approx30.27668039588533,\quad \rho_b\approx129.92476284967807,
\quad\eta_b\approx-8.405268308664358\times10^{-5},
\quad p_b\approx-1.306991661942445\times10^{-4}.
\]

Those data are guidance only. Small numerical junction residuals and
agreement under tolerance changes do not supply an interval enclosure.
No archived producer was executed in this audit. They are not inputs to
the new rational shell/collar bounds.

No closed-form actual `Omega(s),Phi(s)` was found in the focused sources.
They are defined by the analytic regular ODE solution and exact coordinate
conversion `dr/r=dy/rho`, `r/y->1`, `rho=r Omega`, `phi=Phi(r²)`.
The shell-normalized `z` problem and bounds above do not require the
numerical value of conformal radius `b`.

The accessible `D-Blast 3` audit was limited to exact-model signature
searches and selected relevant references. S85 uses a quadratic tension
correction and a cone near `phi=-1`; M465/M468/M471 concern the different
reference/tail domains already excluded by the prior source map. None
supplies the present branch's shell transfer. This was not an exhaustive
reading of the entire archive or a claim of access to unmounted folders.

## 5. Provenance and reproducibility

`BACKGROUND_SOURCE_EXTRACTION.json` records 10 exact Git blobs and three
local reference files, with bytes, SHA-256 hashes and the observed Git
`origin/main` commit `4b85e53307d3f2531e84c45a39cd4acf0aad788e`.
Principal input pins are:

| Frozen source | SHA-256 |
|---|---|
| October 11 Jacobi endpoint theorem | `f1e2571a551ec122dd7755a13145bae73bd92bd7f7426bdc3002e2f52f8dbfc8` |
| Exact rescaled certificate | `22c6cd4e505ee2c7e549420355b94a5b3c079ebb6945314f5b285c13bb67835d` |
| Retained October 8 baseline | `544c6dbced51fbc68dccab3bbbe84210545849c6e7447e1b6bdae0fc98e76f51` |
| Current physical response proof | `a78f2239893fdd0e5ff653d4d5ce5cd9d3244c57c794ed2d3bb43bd4f5e1484b` |

From this directory run the newly written standard-library program:

```sh
python -I -B derive_shell_box.py
python -I -O -B derive_shell_box.py
```

Both modes passed. The program authenticates the two analytic inputs,
verifies that its outward constants dominate the retained exact constants,
derives the shell box, checks all strict first-exit inclusions, and checks
the stated simple coefficient bounds. It writes
`ANALYTIC_SHELL_BOX_DELTA_1E3.json`. No numeric ODE root or historical
binary profile enters any pass/fail decision.

Remaining work for a validated waveform is a separate constructive response
enclosure and physical-output error estimate using these actual coefficient
bounds. The present note supplies a valid small-collar input for that task;
it does not itself prove the response magnitude.
