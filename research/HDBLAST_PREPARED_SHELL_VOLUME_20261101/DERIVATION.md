# A conditional prepared shell-volume and expansion response

1 November 2026. Retain the regular near-`phi=+1` Einstein–scalar branch,
potential/tension, derivative conventions and analytic premises of
[the exact-family binding](../HDBLAST_EXACT_FAMILY_TT_TRANSFER_20261031/DERIVATION.md).
The new mean preparation below changes the nonlinear completion; it does not
identify all second data with the earlier weighted-only family.

## 1. Observables and the conditional statement

Anchor intrinsic future unit-normal geodesics at the initial physical shell
section. Their elapsed time is proper time, and their slices have metric
`h_e(tau)`. Define

```
V_e=integral_S3 dvol_h_e, Theta_e=partial_tau log V_e,
<phi>_e=V_e^-1 integral_S3 phi_e dvol_h_e,
M2=partial_e^2 <phi>_e|0,
v2=(partial_e^2 V_e|0)/V0, Theta2=partial_e^2 Theta_e|0.
```

All coefficients are second derivatives, with the parameter convention
`g_e=g0+e q+e^2 k/2+o(e^2)`. The first averaged variations vanish for the
registered quadrupole. Fix `T=5/4`, `ell=10^-6`, and
`tau_star=rho_b ell T`, using the exact selected background value `rho_b`.
Its certified interval encloses that value; no rounded clock is substituted.
The target is on the shell `r=0`, not the earlier interior sphere.

Under the retained conditional existence/regularity and actual-field premises,
including C2 parameter dependence in the retained lower C4 spacetime topology
and a common intrinsic normal tube through the target after family-dependent
parameter shrinkage,
the completed rational calculation proves

```
M2(tau_star) in [1.1296e30,1.9779e30],
v2(tau_star) in [-5.4520e16,-2.6521e16],
Theta2(tau_star) in [-9.0208e21,-4.4458e21].
```

The background shell has `A=rho_b cosh(tau/rho_b)` and
`Theta0=3 tanh(tau/rho_b)/rho_b>0` at the target. C2 dependence, the strict
coefficient signs and positivity of the spatial metric imply that for each
retained prepared family there is an uncomputed `delta>0` such that

```
0<|e|<delta implies 0<V_e(tau_star)/V0(tau_star)<1
                         and 0<Theta_e(tau_star)<Theta0(tau_star).
```

This is a smaller volume and slower positive expansion relative to the same
background at the specified proper time. It is not negative total volume,
total contraction or a volume-acceleration sign. No usable numerical amplitude
interval, SI calibration, signed finite time history or cosmological origin is
inferred.

## 2. A specified compact mean preparation

Use the retained exact first seed and outward proper-radial coordinate y,
with `dy=-rho dz`, `p=phi0_y` and `a=rho_y/rho`. Let `QH,QM` be its
normalized canonical constraint Hessian means. Choose the same flat packet
bump `b=Q(z(y))`, supported on `ell<z<2ell`, and define exact quadratures

```
I=integral rho^4 QH dy, J=integral rho^2 b dy>0,
f=(2I/(9J))b,
C[f]=-3(f_yy+3a f_y)/2-Vbulk f,
xi=-(rho^4 p^2)^-1 integral_(yL)^y rho^4(QH+C[f])ds,
yL=y(2ell)<yU=y(ell).
```

These are definitions using the true fixed background, not midpoint choices of
their numerical enclosures. The earlier mean identity gives
`integral rho^4 C[f]dy=-(9/2)integral rho^2 fdy`, so the primitive has exactly
zero total integral. Flat endpoint data and the nonzero p bound make xi extend
smoothly by zero outside the packet. Choose mean canonical second data

```
h2=f gamma0+Lie_(xi partial_y)gamma0,
(pi2/varrho0)_rr=0, (pi2/varrho0)_angular=-QM/(6a),
phi2=varpi2=0.
```

Here `varrho0` is the fixed background spatial volume density. The positive
a bound legitimates the momentum division. The exact radial
operators give `Pg h2=-QH` and `Pm tau2=-QM`; the momentum normalization is
`Pm(C,D)=-2[C_y+3a(C-D)]`.

Add this prescribed monopole as the `e^2/2` trial term in the retained weighted
implicit-function construction. Its quadratic residual has zero rotational
mean. Haar averaging is bounded on the full radial weighted spaces and maps
the reflection-invariant domain back into itself. Equivariance and injectivity
of the retained inverse imply that its correction has zero second monopole.
Thus the exact family's second mean is the displayed preparation. The unknown
nonmean second data and higher orders can still depend on the weights.

All supports lie inside the earlier correction annulus; the initial shell collar
is exactly background. The same corner, harmonic-completion and prescribed-slab
arguments apply after qualitative parameter shrinkage. The first tangent and
second TT initial pair are retained. None of this supplies a physical mechanism
that selects this mean preparation among other admissible choices.

## 3. Condition the second jet without changing the anchored response

In mean order `(K,L,N,J,w)`, the original displacement is
`(f+2xi_y,0,-u^2/3,f+2a xi,0)`. Its large radial Lie terms make direct
componentwise volume bounds poorly conditioned.

Use a pure second generator with initial values

```
Z^y=-xi, Z^t=0, Z^y_t=0, Z^t_t=a xi,
(Box_g0+Ric_g0)Z=0,
Z^y|shell=0, (Z^t)_y|shell=0.
```

It yields the conditioned displacement and scaled velocity

```
U0=(f,0,-u^2/3,f,-p xi),
U_x0=ell*(2rho QM/(3a)+32u u_t/3,
         rho[-f_y+a f-13u u_y/3-a u^2/3],
         -2rho QM/(3a)+4u u_t/3,
         -4u u_t/3,0).
```

The APIs and direct initial connection/gauge checks derive every entry.
The scalar term `-p xi` is essential; removing the metric Lie term while leaving
the scalar second datum zero would change the physical preparation.

The original prepared exact family already exists conditionally. A pure-second
pullback gives its transformed jet `k+Lie_Z g0`, `phi2+Z phi0`; Z is tangent
to the shell and preserves the initial anchor. The background vector-wave
equation makes this jet satisfy the same forced linear reduced system. This
does not claim that the full finite-parameter pullback is in nonlinear harmonic
gauge. Jet naturality and the original C2 family's geometric observables are
sufficient for the coefficient calculation.

The mean projection is autonomous because the background linear operator is
rotationally equivariant. Unknown nonmean second fields do not enter it, while
all quadratic forcing is fixed by the unchanged first field. This justifies
evolving all five means with the new data, without inserting an old five-mean
solution or substituting an exceptional zero harmonic into an unrelated master
formula.

## 4. A fresh total mean bound

On the retained causal triangle up to T, use Gram weights `(1,2,1,3,1)` and

```
X(U)=max(sup|U|_H/3,sup|U_x|_H,sup|U_r|_H).
```

The exact prepared data and forcing vanish in the linear-coefficient zero wedge
`x+r<=1`. This statement concerns the background-characteristic coefficient
problem; it does not identify all finite-parameter nonlinear cones.

Write `U=B+mathscrT U+S`, where B is the exact free reflected wave of the
conditioned data. The calculation retains the full coupled matrices and
physical shell lift, reflecting the complete source. Its new bound is

```
q=1201/3200000<1,
X(S)<=T*17387513156573+3851400,
X(U)<=(X(B)+X(S))/(1-q),
X(U-B)<=(q X(B)+X(S))/(1-q)=Eerr.
```

The new initial profiles give `X(B)<1.748e31` and `Eerr<6.563e27` with safe
outward rounding. No older total U cap or initial-approximation allocation is
reused. All supplied bulk and shell source histories are charged even though
their contribution is small relative to the prepared mean.

The scalar free boundary value is exactly `B_w(x,0)=-p(x)xi(x)` because its
initial velocity is zero and its reflection parity is even. The zero wedge
improves its boundary error to `(x-1)_+ Eerr`; the corresponding weighted
history error is `delta_x^3 Eerr/6`, with `delta_x=T-1`.

The proper scalar mean also needs the first clock and measure terms:

```
M2=w0-Cclock psi'/3+(Efirst-8d/3)psi,
Cclock=P+rho_b b_gauge, P'=u, Efirst=u-HP.
```

Here `psi'` is a proper-time derivative and `b_gauge,d` are the unchanged first
harmonic gauge amplitudes. The entire correction is bounded from actual first
metric/scalar/background jets. The background scalar is constant intrinsically,
so its normalization by the actual volume introduces no unknown second metric
mean. These bounds give the positive M2 interval and a signed weighted history,
not just a sign at a sampled point.

## 5. Recover volume and expansion from the intrinsic trace

Define `mathcalF=2Vbulk-sigma_phi^2/4+sigma^2/3`. Standard Gauss and intrinsic
averaging give `R4=(D phi)^2+mathcalF(phi)` and

```
v2''+4H v2'-2v2/A^2=(mathcalF'_b/2)M2+Q1,
Q1=mathcalF''_b psi^2/12-psi'^2/6+4psi^2/(3A^2)
   -5Efirst^2/(3A^2)+Q_D,2/2.
```

All primes on time amplitudes are proper-time derivatives. The scalar-field
derivative is evaluated directly as
`mathcalF'=2Vbulk'-sigma_phi sigma_phiphi/2+(2/3)sigma sigma_phi`;
the authenticated shell boxes give `mathcalF'_b<0`. Define the kinematical term by

```
Q_D,2/2=Efirst'^2-(16/3)Efirst'D'+(8/3)D'^2,
D'=d'-Vshift/A-Cclock/A^2.
```

Here D is the synchronous angular Hessian coefficient and Vshift is the retained
first harmonic shift amplitude. In the longitudinal representation the same
identity is `D'=-P/A^2`. The full harmonic estimate is charged. The complete first quadratic
term Q1 is charged using actual inputs; it is not set to zero. The initial shell
volume and its first time derivative have zero second changes by the open
background collar.

For dimensionless t, homogeneous solutions are `tanh(t)` and
`1+tanh(t)^2`. The exact retarded kernel is

```
G(t,s)=[tanh(t)-tanh(s)] [1-tanh(t)tanh(s)] cosh(s)^4.
```

The scaled kernel and proper-time factors are

```
k(x,s)=G(ell x,ell s)/ell
      =sinh(2ell(x-s))/(2ell)*cosh(ell s)^2/cosh(ell x)^2,
v2(x)=(rho_b ell)^2 integral_0^x k(x,s)[(mathcalF'_b/2)M2(s)+Q1(s)]ds,
Theta2(x)=rho_b ell integral_0^x partial_x k(x,s)[(mathcalF'_b/2)M2(s)+Q1(s)]ds.
```

They have the checked positive bounds
`0<=k(x,s)<=x-s`, `0<=partial_x k<=1`. The retained kernel-error bounds,
whole-cell M2 moments, source error history and all Q1 terms yield the signed
v2 and Theta2 enclosures. A target-point scalar sign alone would not establish
this conclusion; the time integrals and their errors are essential.

## 6. Evidence and limits

The initial campaign froze its N1024 cover, flat-tail treatment and targets
before evaluation. Gauge conditioning is separately frozen postprocessing of
the same raw cells; no grid or target was retuned. Independent rational code
reaggregates the I/J contributions, primitive, conditioned components, and
proper-response error ledger without importing numerical producers. Shared
background/actual-field boxes remain shared premises, not new observations.

The entire statement is conditional on retained analytic theorems and their
application to the specified preparation. Passing finite checks, source
authentication and CI is not external peer review, a proof of novelty or a
proof-assistant formalization. The general constraint, gauge, contraction and
averaging machinery is established mathematics. The contribution is this
specified preparation and the bounded model-specific intrinsic response.

The preparation is not uniquely selected by the blast hypothesis. No numerical
nonzero amplitude interval is certified, and the large coefficients supply no SI-scale prediction.
The full angular waveform, nonlinear stability, ordinary matter production and
higher-dimensional cosmological origin remain unresolved.
