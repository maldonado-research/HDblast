# Later retained states: exact endpoint defects and continuous reconstruction bounds

Prepared 8 October 2026 UTC. This document is a static proof and manufactured-control contribution. It reads no retained array and evaluates no physical source or trajectory. Numerical application requires a separately frozen implementation, authenticated inputs, independent review and immutable public byte GO under the inherited next-calculation contract.

## 1. Target and the strongest available claim

Fix the inherited scalar linear model, represented positive momentum `k`, represented positive weight `w`, exact represented `Pi`, and `L(t)=-1/t` on `[a,b]=[-9/2,-7/2]`. Set `U=u_1/epsilon`, `W=w_1/epsilon`, with division by the exact represented epsilon once. The exact comparison evolution is

```
U'=W,             W'=2 i k W-g(t),
U(a)=U_saved(a),  W(a)=W_saved(a).
```

Here `g` is the fixed prescribed real source; it is not a rounded source history silently promoted to the exact model. Let `t0=a<t1=-4<t2=b` be the retained interior observation times. The primary achievable numerical result is a certified rectangle for

```
D_U(tj)=U_saved(tj)-U_exact_from_saved_a(tj),
D_W(tj)=W_saved(tj)-W_exact_from_saved_a(tj),       j=1,2.
```

These defects include all deviations of the saved endpoint from the declared exact evolution, regardless of whether their origin was quadrature, implementation, source evaluation, phase arithmetic or roundoff. The defect does not identify a cause. It is a well-defined forward-error target even if the historical internal stages were never retained.

The initial defect for this target is exactly zero. The earlier BD preparation error is a different contribution: by linearity, the later saved-minus-BD-exact error equals this later defect plus homogeneous transport of the already enclosed saved-incoming-minus-BD-incoming error. Reuse that earlier enclosure once. Reinitializing the trajectory at the BD midpoint would change the target.

The producer `research/HDBLAST_CHECKPOINT_20261002_METRIC_RESPONSE_MEMORY_FOLLOWUP/independent/forced_metric.py` retains mode `u,w` at observation steps, while its full time histories contain aggregate stresses/source/geometry. Static code inspection does not establish possession of the unsaved per-mode path. A theorem for a newly specified interpolant through retained points is valid for that interpolant, not automatically the historical internal solver path.

## 2. Exact unitary propagator and regular zero-momentum limit

For real `k` and nonnegative lag `d`, define

```
E_k(d)=exp(2 i k d),
Phi_k(d)=integral_0^d E_k(v) dv
        =d exp(i k d) sinc(k d),                 sinc(0)=1.
```

Thus `Phi_0(d)=d`, `Phi_k'=E_k`, `E_k'=2ik E_k`, and

```
|E_k(d)|=1,
|Phi_k(d)|=|sin(kd)|/|k| <= min(d,1/|k|)        (k nonzero).
```

The integral and sinc definitions, or their entire series, must be used near `k=0`; subtracting two ordinary floating-point exponentials and dividing by a small `k` is not a validated implementation. The propagator is

```
S_k(d)=[[1,Phi_k(d)],[0,E_k(d)]].
```

Its Euclidean operator norm is exactly `(sqrt(4+|Phi|^2)+|Phi|)/2`. This follows by computing the two eigenvalues of `S* S`, whose trace is `2+|Phi|^2` and determinant is one. The maximum grows at most linearly in `d`; an `exp(2|k|d)` stability factor is unnecessary. The `W` equation is unitary, although the complete triangular `U,W` system is not unitary.

These statements use complex modulus. Multiplication by `E` is not an isometry for the Cartesian complex L1 norm. A rectangle `x in [xl,xu], y in [yl,yu]` supplies the modulus upper `sqrt(max|x|^2+max|y|^2)`, or the rational upper `max|x|+max|y|`. Rational square-root uppers can be certified by squaring; no unvalidated floating square root is required.

## 3. Piecewise a posteriori theorem, including exact jumps

Let `Uhat,What` be a specified piecewise absolutely continuous reconstruction on a finite partition. Let `f` be a represented source. On each open panel define

```
r_u=Uhat'-What,
r_w=What'-2 i khat What+f,
F=r_w-(f-g)+2 i (khat-k) What.
```

Then the exact error `e=hat-exact`, relative to the fixed-`k`, fixed-`g` target, obeys

```
e_U'=e_W+r_u,
e_W'=2 i k e_W+F.
```

If the reconstruction and target have jump differences `J_U(c),J_W(c)` at internal interfaces `c`, variation of constants gives, using right-hand values at the final time,

```
e_W(t)=E(t-a)e_W(a)+integral_a^t E(t-s)F(s)ds
       +sum_(a<c<=t) E(t-c)J_W(c),
e_U(t)=e_U(a)+Phi(t-a)e_W(a)+integral_a^t r_u(s)ds
       +integral_a^t Phi(t-s)F(s)ds
       +sum_(a<c<=t) [J_U(c)+Phi(t-c)J_W(c)].             (1)
```

For the bounded piecewise source in the inherited problem the exact target is absolutely continuous and has no state jumps. A jump in `f` creates no impulse by itself. In this case the `J` are just reconstruction jumps. If the target has explicitly prescribed delta-function impulses, their known jumps must first be subtracted; (1) must not count them as numerical defects.

Given nonnegative modulus envelopes `q_u>=|r_u|`, `q_F>=|F|`, and initial/jump modulus bounds, put `p_k(d)=min(d,1/|k|)`, with `p_0(d)=d`. Then

```
B_W(t)=B_W(a)+integral_a^t q_F(s)ds+sum B_JW(c),
B_U(t)=B_U(a)+p_k(t-a)B_W(a)+integral_a^t q_u(s)ds
       +integral_a^t p_k(t-s)q_F(s)ds
       +sum [B_JU(c)+p_k(t-c)B_JW(c)].                  (2)
```

The triangle inequality proves `|e_X(t)|<=B_X(t)`. Replacing `p_k` by the exact `|Phi_k|` gives the sharper same proof. For any fixed terminal time, the integral operator norm `F -> e_W(t)` is exactly the integral of its nonnegative pointwise forcing radius, because a complex forcing can align with the conjugate phase. The analogous `F -> e_U(t)` norm is exactly `integral |Phi|q_F`. No uniformly smaller general bound follows from these disk envelopes alone. Correlated real forcing or known polynomial residuals permit stronger bounds below.

For constant envelopes on source panels `[l,h]`, every integral in (2) is explicit. For an already completed panel and terminal `t>=h`,

```
integral_l^h q ds=q(h-l),
integral_l^h (t-s)q ds=q[(h-l)(t-h)+(h-l)^2/2].
```

The integral of `min(t-s,1/k)` is obtained by splitting at `s=t-1/k`, if inside the panel. All breakpoints and results are rational for rational `k,t,l,h,q`. Evaluating these positive bounds at the rightmost final time gives an all-time majorant. For nonzero initial data, `min(T,1/k)B_W(a)` is likewise valid on the entire interval. No sampling of a later trajectory is part of this all-time proof.

## 4. Direct stress kernels that preserve oscillatory cancellation

Write the direct mode operators as

```
X_m(U,W;t)=Re[A_X(t) U+B(t) W],
A_R=k+3L^2/(2k),            A_P=k/3-L^2/(2k),
B=-L/(2k)+i/2.
```

Let `(beta_R,alpha_R)=(1,3)` and `(beta_P,alpha_P)=(1/3,-1)`, so `A_X=beta_X k+alpha_X L^2/(2k)`. Equation (1) supplies the exact `W`-impulse kernel

```
H_X(t,d)=A_X(t)Phi_k(d)+B(t)E_k(d).                    (3)
```

Consequently the initial/state-jump `U` component and `r_u` act through `A_X(t)`, while the `W` component, `W` jump and `F` act through `H_X(t,d)`. Direct algebra also gives

```
H_R=i/2+[-L E+3L^2 Phi]/(2k),
H_P=i/6+i E/3+[-L E-L^2 Phi]/(2k).                    (4)
```

The leading density oscillations cancel. Bounding `U,W` separately before applying the density operator discards this useful cancellation. The sharp support for a disk uncertainty of radius `q` in a complex impulse is `|H_X|q`. For a real scalar forcing it is `|Re H_X|q`. For independent real/imaginary component radii `q_R,q_I` it is `|Re H_X|q_R+|Im H_X|q_I`. These are exact support functions, not empirical estimates.

Apply the retained finite measure before numerical divisions:

```
mu=w k^2/(2 Pi^2),           Q=w/(4 Pi^2),
mu A_X=Q k(2 beta_X k^2+alpha_X L^2),
mu B=Q k(-L+i k),
mu H_X=Q k[(2 beta_X k^2+alpha_X L^2)Phi+(-L+i k)E].  (5)
```

Equation (5) contains no inverse momentum and is continuous at `k=0`, with value zero for fixed finite `w` and bounded state impulses. There is no claim that the normalized single zero mode exists; this is a stable extension of the weighted coefficient. For fixed positive retained nodes it avoids artificial infrared ill-conditioning and large common denominators.

A real source-error impulse has the particularly simple weighted kernels, with `theta=2kd`,

```
Re(mu H_R)=Q[-L k cos(theta)+(3L^2/2)sin(theta)],
Re(mu H_P)=Q[-L k cos(theta)-(2k^2/3+L^2/2)sin(theta)]. (6)
```

For example, rational, globally valid pointwise uppers are

```
|Re(mu H_R)| <= Q[L k+(3L^2/2)min(1,2kd)],
|Re(mu H_P)| <= Q[L k+(2k^2/3+L^2/2)min(1,2kd)],       (7)
```

on the inherited positive-`L`, positive-`k` interval. Alternatively, the global phase-amplitude bounds are `Q sqrt((Lk)^2+(3L^2/2)^2)` and `Q sqrt((Lk)^2+(2k^2/3+L^2/2)^2)`. Take the smaller of separately proved bounds, or enclose the actual trigonometric kernel on a registered phase interval. Lag integration of (7) is rational for rational panel data and can be materially tighter near `k=0` than a phase-independent envelope.

For arbitrary complex state error, a simple rational finite-weight bound from (2) is

```
|mu e_R| <= Q k[(2k^2+3L^2) B_U+(k+L) B_W],
|mu e_P| <= Q k[|2k^2/3-L^2| B_U+(k+L) B_W].          (8)
```

Here `sqrt(k^2+L^2)<=k+L` bounds the complex coefficient multiplying `W`. The exact square root or component rectangle support can replace it. For a uniform pressure bound one must maximize the absolute coefficient over the declared `L` interval: `max(|2k^2/3-La^2|,|2k^2/3-Lb^2|)`, or use `2k^2/3+Lb^2`. Substituting only `Lb` into an absolute expression is not always valid.

Sum any per-node bound with the actual positive `w_j` and the original twelve prefixes. Exact signed interval sums or direct common-source kernel sums can retain more cancellation; positive triangles never assert unknown cancellation. These are finite-measure statements. Nodes, weights and continuum quadrature errors are separate targets.

## 5. Certified endpoint evolution by entire polynomial moments

On a panel starting at `c`, represent `f(c+s)=sum_(m=0)^M a_m s^m`, with exact rational real coefficients, `0<=s<=h`. Define the entire moments

```
M_m(k,h)=integral_0^h E_k(h-s)s^m ds
        =sum_(n>=0) (2ik)^n m! h^(n+m+1)/(n+m+1)!,
N_m(k,h)=integral_0^h Phi_k(h-s)s^m ds
        =sum_(n>=0) (2ik)^n m! h^(n+m+2)/(n+m+2)!.
```

Integrating the exponential series term by term is justified by uniform absolute convergence on the compact panel; the beta integral yields the displayed coefficients. Thus

```
W_f(c+h)=E_k(h)W_f(c)-sum a_m M_m(k,h),
U_f(c+h)=U_f(c)+Phi_k(h)W_f(c)-sum a_m N_m(k,h).        (9)
```

At zero momentum, `M_m(0,h)=h^(m+1)/(m+1)` and `N_m(0,h)=h^(m+2)/[(m+1)(m+2)]`. No recurrence dividing repeatedly by `k` is necessary.

For truncation after `n=N`, the absolute first omitted terms are

```
T_M=(2|k|)^(N+1) m! h^(N+m+2)/(N+m+2)!,
T_N=(2|k|)^(N+1) m! h^(N+m+3)/(N+m+3)!.
```

Every subsequent ratio is at most `rho=2|k|h/(N+2)`. If `rho<1`, `T_M/(1-rho)` and `T_N/(1-rho)` are rigorous tails. The deliberately loose common ratio also covers the ordinary exponential series; an implementation may tighten denominators with a separately proved formula. Rational coefficients and rational real/imaginary series arithmetic plus rational tails give exact rational rectangles. Panel subdivision with `2|k|h<=1` avoids expensive long series at high momentum. An exact represented phase multiplied using rational interval arithmetic has a bounded arithmetic enclosure; an ordinary complex library phase is not treated as exact.

If `|f-g|<=q_g`, exact same-initial flow differences satisfy

```
|W_f(t)-W_g(t)|<=integral_a^t q_g(s)ds,
|U_f(t)-U_g(t)|<=integral_a^t p_k(t-s)q_g(s)ds.        (10)
```

Use (6)-(7) for sharper source-induced direct stress uncertainty. Since the forcing discrepancy is real, the disk bound need not replace the real kernel. Add these model/source-representation uncertainties once to the arithmetic/series target enclosure, then subtract that complete target rectangle from the exact later saved value. This yields endpoint defect intervals. Component intervals that exclude zero certify an actual endpoint discrepancy; a positive upper alone does not do so.

For a sharper phase-sensitive integration of a known residual polynomial, apply the same moments to `F` and ordinary polynomial integrals to `r_u` in (1), instead of bounding its modulus first. This is still an a posteriori error computation; it does not require the historical solver's local-error estimator.

## 6. Representation, arithmetic and geometry are explicit hypotheses

1. **Stored-value semantics.** A finite stored floating-point bit pattern, decoded correctly, is an exact rational. A chosen rational interpolant through those values has no unknown interpolation-arithmetic roundoff if its coefficients and derivatives are computed exactly. The endpoint defects already contain the producer's deviation from the declared target. Do not add an invented epsilon-times-operation-count error to those exact inputs.
2. **Represented source coefficients.** If `a_m` is chosen exactly inside a certified coefficient rectangle of radius `e_m`, and an analytic remainder is bounded by `q_tail`, then `q_g(s)=q_tail(s)+sum e_m |s|^m` is valid. If a previously exported complete `q_g` already includes both terms, reuse that complete bound instead of adding them again. If the coefficients are uncertain and interval arithmetic directly propagates their uncertainty, do not add the same coefficient radii afterwards.
3. **Unknown reconstruction coefficients.** If an intended curve is `Uhat+delta U`, `What+delta W`, with separately proved envelopes for the functions and their derivatives, then its residual differs by at most `|delta U'|+|delta W|` in the first equation and `|delta W'|+2|k||delta W|` in the second (plus declared coefficient/source differences). The coefficient basis supplies these derivative bounds. Endpoint coefficient accuracy alone supplies no derivative bound for an unspecified curve. Choosing the exact center curve instead needs no such term, but changes the object certified to that center curve.
4. **Momentum coefficient.** The fixed represented `k` is exact by definition. If a different `khat` was used in the reconstructed operator, the correction `2i(khat-k)What` in (1) is mandatory. Even a real frequency uncertainty creates phase error, despite every individual exact frequency having unit modulus. This correction uses the known reconstructed `What`, so it requires no exponential Gronwall factor. A genuinely non-imaginary unknown growth coefficient would require a separately justified stability bound.
5. **Geometry.** Equations (3)-(8) use exact `L=-1/t`. If stress evaluation uses `Lhat`, its extra weighted discrepancy at a fixed reconstructed state is exactly

   `Q k[alpha_X(Lhat^2-L^2)Re Uhat-(Lhat-L)Re What]`.

   For `|Lhat-L|<=e_L`, `|Lhat^2-L^2|<=e_L(|Lhat|+|L|)` provides a product bound. Evaluate with the same fixed measure. An uncertain measure/Pi or uncertain epsilon normalization is another declared multiplicative uncertainty, not silently absorbed by a state residual.
6. **Export.** Rounded exported intervals must be enlarged outward. Serialized midpoint/tolerance alone is not a complete enclosure. Exact rational exports need no extra decimal-display allowance; display rounding never determines signs or acceptance.

An analytic global bound such as the inherited `|g|<=64` may be used with `f=0` if its pinned proof covers both registered source definitions and the full interval. This proves a conservative enclosure with no physical source callback, but can be far too wide to resolve an endpoint error. It is a fallback analytic bound, not a claimed accurate trajectory calculation. A useful small endpoint defect generally requires a tight represented source and validated phase/moment arithmetic.

## 7. Piecewise contacts, work and what the Ward check does not prove

For full normalized stress `X=X_m+C_X`, add independent `|Chat_X-C_X|` envelopes. Identical complete exact contacts cancel in the specific two-evolution state comparison. They do not cancel merely because a later stored full-stress value has been loaded. A comparison involving stored full stress requires independent direct-operator/geometry/contact/accumulation error bounds as applicable.

On open panels, for exact fixed `k,L`, reconstructed source `f`, and the inherited matching baseline, define

```
zeta_C=Chat_R'-L(Chat_R-3Chat_P)+3 h' Bhat_0+L g/(2k).
```

The full off-shell identity is

```
Rhat'-[L(Rhat-3Phat)-3h'Bhat_0]
 = [(2k^2+3L^2)Re r_u-k Im r_w-L Re r_w]/(2k)
   +L(f-g)/(2k)+zeta_C,                              (11)
```

where here `r_w=What'-2ikWhat+f` (same fixed `k`). If `khat` differs, first replace `r_w` by `r_w+2i(khat-k)What`. The canonical invariant has drift `c'=Re r_u-Im r_w/(2k)` for real `f`; dropping it by checking only the noncanonical density piece is incorrect.

For piecewise absolutely continuous stress, the integrated version of (11) adds the signed jumps of the **full** density `Rhat`, including its contacts, at every interface. A jump in a density contact can give a nonzero integrated defect while every open-panel derivative defect vanishes. Pointwise contact error bounds and interface jump bounds are different inputs. The historical full contact inventory remains a prerequisite.

Given actual full-stress error bounds `e_R,e_P`, direct work error obeys `integral L(e_R+3e_P)` plus any independently changed baseline/source-jet terms and the ledger quadrature error. Direct pressure time-integral error obeys `integral e_P` plus its own quadrature error. The identity `(-R/(3L))'=P` used for homogeneous incoming-error transport does **not** apply unchanged to a forced/residual/contact error. Pressure is independently evaluated, never repaired from a Ward residual.

Conservation alone cannot bound state error: every homogeneous perturbation is on shell. For real constant `A`, `delta U=A exp(2ik(t-a))`, `delta W=2ik delta U` has zero canonical defect and zero mode Ward defect, while its initial density is `3L(a)^2 A/(2k)` and can be arbitrarily large. This is an independent negative control against promoting a Ward PASS to trajectory accuracy.

## 8. Retained samples do not determine an unsaved trajectory

Let the finitely many sample times be `tj`, and define `P(t)=product_j(t-tj)^2`. For arbitrary real `A`, the perturbation

```
delta U(t)=A P(t),          delta W(t)=A P'(t)
```

vanishes in both stored coordinates at every `tj` and obeys the first kinematic equation exactly. In every open interval where `P` is nonzero, its `U` magnitude is unbounded as `|A|` grows. Its second residual is `A[P''-2ikP']`. Thus finitely many exact state values, even together with the kinematic equation, do not imply any finite universal between-sample error/residual bound. Repeated roots of larger multiplicity similarly hide any finite number of retained derivatives. A smooth compactly supported bump between samples gives the same obstruction for arbitrary prescribed finite endpoint jets.

This obstruction does not apply after the full exact ODE, forcing and initial state have been imposed; that IVP has a unique solution. It explains why samples cannot show that an unspecified historical interpolation obeyed that ODE to a bounded residual. A declared exact interpolant supplies the missing function and derivative information and can then be certified by (1). Its certificate answers a different, explicit question.

Finite momentum samples also cannot control an unspecified continuum state/forcing between nodes: a smooth bump supported in a gap vanishes at all nodes yet has arbitrarily large integral after scaling. Positive retained weights certify the finite sum, not the continuum integral or ultraviolet tail. Full stored stress aggregates at additional time points do not recover the missing per-node complex path: finitely many linear aggregates have a nontrivial kernel unless additional dynamics or regularity constraints are supplied.

## 9. Complete minimal study and honest stopping conditions

A minimal complete numerical continuation can certify the original later mode values at both interior observation times for all 49,152 node occurrences and twelve retained prefixes. It needs exact decoding of the selected later states; the authenticated incoming saved states; a certified representation of the same source on the interval; validated evaluation of (9)-(10); exact direct finite stress discrepancy intervals at both times; export/readback and independent checks. It can separately add the previous incoming-state contribution to give saved-minus-prescribed-BD endpoint errors. This requires no claim to possess the unsaved solver trajectory and no new physical model.

Before numerical execution, the complete proof and exact manufactured controls already establish the stable conditional endpoint/trajectory theorem, the removable zero-momentum limit, the correct source and contact accounting, and a constructive obstruction to recovering a historical continuous trajectory from snapshots. Those are completed mathematical results. No numerical enclosure for later retained states has been computed by this contribution.

A narrow but rigorous endpoint result can be successful even if its intervals are too wide to resolve a nonzero discrepancy; report such precision failure separately. A continuous certificate for the historical solver curve remains blocked without an authenticated dense output definition/stage record or a separately validated reconstruction of the original algorithm and its rounding errors. Momentum/time quadrature, complete contact evaluation and UV errors remain separate missing contributions to the historical gate.

Unchanged statuses: incoming finite-node errors **ENCLOSED**; later endpoint theorem **CONDITIONAL, NOT YET APPLIED TO RETAINED STATES**; historical continuous pressure/contact/time/momentum/UV certificate **UNRESOLVED**; metric calibration **FAIL**; higher-dimensional origin **NOT_ESTABLISHED**; external novelty **NOT_ASSESSED**. The mathematics here is established variation of constants, polynomial enclosures and linear support-function analysis; no novelty claim is made.
