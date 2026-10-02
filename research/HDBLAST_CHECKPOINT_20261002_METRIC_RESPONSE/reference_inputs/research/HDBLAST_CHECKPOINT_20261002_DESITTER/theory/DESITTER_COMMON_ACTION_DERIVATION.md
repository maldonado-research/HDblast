# Massive Euclidean de Sitter sources from one fixed-reference action

Prepared 2 October 2026. This is a mathematical benchmark for a free, minimally coupled real scalar. It does not choose physical HDBLAST masses, gravitational couplings, or EFT matching data. The state is the regular Euclidean state on S4 continued to Lorentzian de Sitter. Throughout, `z=H²>0`, `x>0`, and the reference `r>0` is held fixed in metric and scalar variations. The massless endpoint is excluded.

## 1. A finite action, rather than a trace definition of stress

Put `h=x-r`, `C=1/(16π²)`, and

```
K(t) = Σ[l≥0] d_l exp[-l(l+3)t],
d_l = (l+1)(l+2)(2l+3)/6,
A(t) = 6t² K(t),
A(t) ~ 1+2t+(29/15)t²+(74/63)t³+(149/315)t⁴+...,
c = h²/2-2zh+(29/15)z².
```

The heat density on S4 is `(4πs)^-2 A(zs)`. Define the complete finite Euclidean action density by

```
W(z,x,r) = -C/2 ∫[0,∞] ds/s³ {
  e^(-xs) A(zs)
  -e^(-rs)[1+(2z-h)s+c s²] }.
```

No arbitrary determinant normalization remains in this definition. In particular a field-independent but curvature-dependent term cannot subsequently be discarded when computing stress. The action is `Gamma_E=Vol(S4) W`, with `Vol(S4)=8π²/(3z²)`.

The subtraction agrees with the physical heat kernel through order `s²`. The bracket is `O(s³)` at zero, so the action integrand has a finite limit. Positive `x,r` ensure exponential infrared convergence. Differentiation in `x,z,r` is justified on compact subsets bounded away from zero by the differentiated heat expansion and exponentially bounded spectral sum. In these statements convergence is a mathematical property; any particular numerical quadrature needs its own error control.

Define sources by independent variations of this same action:

```
Q = 2 W_x,
rho = W - (z/2) W_z,
p = -rho,
J_phi = x_phi Q/2.
```

The metric formula follows by varying the S4 radius at fixed `x,r`, including its volume. Direct differentiation gives separately convergent integrals

```
Q = C ∫ ds/s² {e^(-xs)A(zs)
     -e^(-rs)[1+(2z-h)s]},

rho = -C/2 ∫ ds/s³ {e^(-xs)[A(t)-t A'(t)/2]
       -e^(-rs)[1+(z-h)s+(h²/2-zh)s²]},  t=zs.
```

These are the primary current and metric-stress definitions. Stress has not been inferred by imposing a trace identity. The `t²` coefficient cancels from `A-tA'/2`; this does not eliminate the trace anomaly.

## 2. Exact connection to the frozen finite matching

Splitting off the first three heat coefficients gives the exact identity

```
W = V(x,r)-12z F(x,r)+(29z²/15)/(32π²) ln(x/r)+W_rem,

V = [x² ln(x/r)-3x²/2+2rx-r²/2]/(64π²),
F = [x ln(x/r)-x+r]/(192π²),
W_rem = -C/2 ∫ ds/s³ e^(-xs)
        [A(zs)-1-2zs-(29/15)z²s²].
```

`W_rem` is retained exactly; the split does not approximate the full source by a heavy-mass expansion. It follows directly that

```
V(r)=V_x(r)=V_xx(r)=0,
F(r)=F_x(r)=0.
```

The geometric logarithm is the S4 restriction of

```
a2_geo = R²/72 + C_mnrs C^mnrs/120 - E4/360.
```

Its Euclidean coefficient is `+ln(x/r)/(32π²)`. Equivalently the Lorentzian loop coefficients multiplying `R²,C²,E4` are respectively

```
-ln(x/r)/(2304π²),
-ln(x/r)/(3840π²),
+ln(x/r)/(11520π²).
```

They vanish at `x=r`. These are loop derivative-expansion functions, distinct from the fixed local EFT counterterm constants selected by matching. The scalar derivative of a loop logarithm must be retained. Thus the action realizes exactly the frozen SMOOTH_FRW finite convention for `V,F,alpha,beta,gamma` on this constant-curvature, constant-mass background, with total derivatives absent on S4. The frozen FRW state and its heavy-regulator argument are not used to define the Euclidean state. The agreement established here is the covariant local finite-action convention.

For `x/z` large, `W_rem ~ -C(74/63)z³/(2x)+...`; this expansion is a check, not a replacement for the exact action at the finite benchmark masses. At `x=r` the local coefficients specified by the matching vanish, but the full curved-space effective action and its sources generally do not.

## 3. The reference derivative fixes the exact trace remainder

The reference dependence is particularly simple. Since `∂r h=-1`,

```
∂r {e^(-rs)[1+(2z-h)s+c s²]} = -c s³ e^(-rs),
W_r = -C c/(2r).
```

Dimensional homogeneity of the full finite action gives

```
z W_z+x W_x+r W_r=2W.
```

Combining this independent fact with the metric and scalar variations yields

```
-4rho = -xQ-2r W_r
       = -xQ + [h²/2-2zh+(29/15)z²]/(16π²).
```

Therefore the prior static bridge's explicitly conditional trace formula is correct for this newly constructed action. There is no need to amend its conditional statement. The reference contribution cannot simply be omitted: if `r` is allowed to track `x`, `2dW/dx=Q+2W_r`, which changes the scalar source. Another finite action changes the local remainder according to its actual variations.

A nonvacuous static consistency identity is

```
rho_x = Q/2-z Q_z/4
       = Q/2-H Q_H/8.
```

This is checked independently of the static Ward identity, which has no time derivatives to test in an invariant state.

## 4. Independent spectral continuation and exact digamma sources

Introduce dimensionless variables and the spectral zeta function

```
u=x/z, v=r/z, nu²=9/4-u,
Z(s;u)=Σ[l≥0] d_l [l(l+3)+u]^(-s),
D(u)=u²/2-2u+29/15.
```

The sum defines `Z` for `Re s>2`; its standard Mellin continuation has

```
Z(0;u)=D(u)/6.
```

Analytically continuing the same reference-subtracted Mellin integral, including all polynomial finite terms, gives

```
W = z²/(32π²) {
  -6 Z'(0;u)-D(u)ln(v)+uv-2v-v²/4 }.
```

This identity fixes the normalization linking the independent spectral route to the convergent proper-time action. It is not sufficient to compute a determinant only up to an `x`-independent term when stress is required.

For the finite part of the resolvent, put `a=3/2`, `n=l+a`, and

```
Psi(u)=psi(a+nu)+psi(a-nu).
```

For imaginary `nu`, the two terms are conjugates and their sum is real. At `u=9/4`, use the analytic coincident-argument limit. For every `u>0`, the arguments avoid nonpositive-integer poles.

To obtain the finite part without a determinant ansatz, subtract the first two large-`n` terms from

```
Z(s;u) = (1/3) Σ n(n²-1/4)(n²-nu²)^(-s).
```

The pole at `s=1` has residue `(2-u)/6`. The remaining convergent resolvent sum uses

```
Σ [n/(n²-nu²)-1/n] = psi(a)-Psi(u)/2.
```

Together with `zeta_H(-1,3/2)=-11/24`, this gives

```
FP Z(1;u) = [(u-2)Psi(u)-u+4/3]/6,
∂u Z'(0;u) = -FP Z(1;u).
```

The second identity follows by differentiating `∂u Z(s)=-s Z(s+1)` through its meromorphic continuation, retaining the pole before taking `s=0`.

The resulting exact current is

```
Q = z/(16π²) {
  (u-2)[Psi(u)-ln(v)]-u+v+4/3 }.
```

Independently applying `rho=W-zW_z/2` to the spectral action gives

```
rho = z²/(64π²) {
  u(u-2)[Psi(u)-ln(v)]
  -u²+4u/3+2uv-2v-v²/2-D(u) }.
```

These formulas are an independent continuation route for both physical sources. The trace identity is subsequently a check of these results. The expression for `rho` is not defined by the trace equation.

For stable source derivatives define `L=Psi-ln(v)`. Differentiating the same closed current gives

```
Psi_u = [psi_1(3/2-nu)-psi_1(3/2+nu)]/(2nu),
Q_x = C[L+(u-2)Psi_u-1],
Q_z = C[-2L-u(u-2)Psi_u+u-2/3].
```

At `nu=0`, the removable limit is `Psi_u=-psi_2(3/2)`. The derivatives hold `r` fixed; differentiating along the special curve `x=r=2z` would be a different operation.

At `u=v=2`, the factor multiplying `Psi-ln(v)` vanishes and the exact source reduces to

```
Q=z/(12pi²),  rho=11z²/(960pi²),  p=-rho.
```

The independent review obtains these same values by separately integrating the exact physical plane-wave modes against the frozen de Sitter adiabatic current, energy, and pressure subtractions. That mode calculation tests the finite convention in curved space; pressure is independently integrated there.

As `x→0+`, the homogeneous mode produces `Q~3z²/(8π²x)`. Removing that mode would change the chosen state and determinant. Although `rho` has a finite massive limit, this does not supply a de Sitter invariant massless Euclidean vacuum; the zero-mode obstruction remains.

## 5. A convergent determinant series with an explicit tail

A useful independent determinant computation needs neither differentiated quadrature nor fitted tail extrapolation. Let `q=u-9/4`, `a=3/2`. For `|q|<a²`, binomial expansion gives

```
Z'(0;u) = (2/3)[zeta_H'(-3,a)-zeta_H'(-1,a)/4]
 -q/3[zeta_H(-1,a)+psi(a)/4]
 +q²/6[1/2-psi(a)-zeta_H(3,a)/4]
 +(1/3)Σ[k≥3] (-q)^k/k
       [zeta_H(2k-3,a)-zeta_H(2k-1,a)/4].
```

The explicit `k=1,2` terms include the zeta poles; treating them as ordinary zero terms loses finite contributions. Put `B_k=zeta_H(2k-3,a)-zeta_H(2k-1,a)/4`, which is positive for `k≥3`. For a sum retained through `N≥2`, `theta=|q|/a²<1` gives

```
|R_N| <= |q|^(N+1) B_(N+1) / [3(N+1)(1-theta)],
|∂u R_N| <= |q|^N B_(N+1) / [3(1-theta)].
```

These follow by bounding the term ratio of the positive spectral series by `a^-2` and summing a geometric majorant. They are truncation bounds; special-function evaluation and finite arithmetic errors remain separate. Any `u>0` can be covered by splitting out finitely many low harmonics and using a larger starting `a`, but the bounded registered benchmark can stay in the displayed disk.

## 6. Proper-time endpoint and harmonic tail bounds

For `t≥T>0`, positivity and the nonzero spectral gap `lambda_1=4` imply

```
K(t)-1 <= [K(T)-1] exp[-4(t-T)],
-K'(t) <= [-K'(T)] exp[-4(t-T)].
```

Consequently the physical long-time tails obey

```
∫[T,∞] e^(-ut)K(t)dt
 <= e^(-uT){1/u+[K(T)-1]/(u+4)},

∫[T,∞] e^(-ut)K(t)dt/t
 <= E1(uT)+[K(T)-1]e^(4T)E1((u+4)T),

∫[T,∞] e^(-ut)[-K'(t)]dt
 <= [-K'(T)]e^(-uT)/(u+4).
```

The reference tails are incomplete-gamma integrals of explicit polynomials. The same positivity permits rigorous harmonic tails: for first omitted shifted harmonic `n=N≥3/2`, use `d_l≤n³/3`. If `tN²≥3/2`,

```
Σ[omitted] d_l e^[-l(l+3)t]
 <= e^(9t/4-tN²)/3 [N³+N²/(2t)+1/(2t²)].
```

For `-K'`, replace the envelope `n³` by `n⁵`, with monotonicity threshold `tN²≥5/2` and integral `e^(-tN²)[N⁴/(2t)+N²/t²+1/t³]`.

At small `t`, an asymptotic polynomial by itself is not a certified remainder. A constructive alternative is Euler–Maclaurin for `f_t(n)=(n³-n/4)e^(-tn²)`:

```
Σ[j≥0] f_t(a+j)
 = ∫[a,∞]f_t(n)dn + f_t(a)/2
   -Σ[k=1..m] B_(2k)/(2k)! f_t^(2k-1)(a) + R_m,

|R_m| <= 2 zeta(2m)/(2π)^(2m)
          ∫[a,∞]|f_t^(2m)(n)|dn.
```

Every derivative is a polynomial times a Gaussian, so a coefficientwise absolute bound is an explicit sum of upper incomplete gamma functions. Its small-`t` scaling is `O(t^(m-2))` for `K`, hence `O(t^m)` for `A`. Choosing `m≥3` gives an integrable action error bound at zero. This is a mathematical construction; implementations must report the actual `m`, endpoint, evaluated bound, and other numerical errors. Step-refinement agreement alone is not called a certified enclosure.

## 7. Meaningful wrong-formula controls and scope

The companion symbolic verifier checks finite matching, spectral residues and polynomial constants, reference differentiation, source pairing, the exact trace, heat coefficients, and general local-action shifts. It deliberately tests failures produced by treating `W` as `rho`, allowing `r=x` during differentiation, omitting the `F_x` current, doubling a finite local term, and replacing the `4/3` finite constant in the resolvent by the incompatible `-2/3`. All 37 exact assertions pass and all six deliberate wrong-formula controls are detected.

This completes an actual analytic quantum source construction for one explicit Euclidean state and one explicit benchmark convention. Numerical evaluation and error evidence are separate deliverables. No supplied physical mass law, no gravitational normalization, no dynamic shell solution, no particle-production population, no radiation transfer, and no thermalization result are inferred. A prospective first-order response must still retain both `rho` and `J_phi`, use independently resolved classical endpoint sensitivities, and state its response scale.
