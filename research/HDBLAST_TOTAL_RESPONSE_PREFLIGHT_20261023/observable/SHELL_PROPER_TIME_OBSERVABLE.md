# Anchored proper-time shell volume: bounded observable specification

23 October 2026. **ANALYTIC_SPECIFICATION_ONLY; NOT_NUMERICALLY_EVALUATED.**
All coefficient conventions below are parameter derivatives:
`g_e = g_0 + e q + e^2 k/2 + O(e^3)`.
This note supplies a kinematical second-order observable and its gauge law.
It supplies neither a finite-amplitude physical solution family nor a new
physical sign, force, backreaction or cosmological conclusion.

## 1. Registered geometry and the anchor

Use the physical, copy-even, timelike shell as an abstract four-manifold
`S`, with embedding `E_e:S -> M`, induced Lorentzian metric
`gamma_e=E_e^*g_e`, a specified compact initial spacelike section
`j_e:S^3 -> S`, and a chosen future time orientation. The initial section is
the registered initial shell section. Its identity is part of the input:
it is transported with the physical configuration under a gauge change.
Keeping a coordinate time surface fixed under a gauge transformation that
moves this physical section would define a different observable.

On the retained background, with `rho=rho_b` fixed,

```
s = rho t,       A(s) = rho cosh(s/rho),       H(s) = A'(s)/A(s),
gamma_0 = -ds^2 + A(s)^2 kappa,               j_0(S^3)={s=0}.
```

Here `kappa` is the unit round three-sphere metric with volume `2 pi^2`.
The retained model uses a fixed-shell chart and shell-preserving gauge.
The initial physical and gauge jets vanish on an open shell collar, so
the registered initial anchor agrees there with the background section.
Initial flatness alone would not justify fixing a later moving shell.
This is a choice of representation, not a rule that the shell may be held
fixed in every gauge.

Let `u_e` be the future unit vector intrinsically normal to `j_e(S^3)`.
Define the **intrinsic** normal exponential map

```
F_e(tau,z) = exp^(gamma_e)_(j_e(z))(tau u_e(z)).
```

The worldlines are geodesics of the shell's induced metric. They need not be
ambient spacetime geodesics or material-flow worldlines. This observer
convention is part of the observable. Whenever `F_e` is a diffeomorphism
onto a shell neighborhood, put

```
Sigma_(e,tau)=F_e({tau} x S^3),
h_e(tau)=F_e(tau,.)^*gamma_e,
V_e(tau)=integral_(S^3) dvol_(h_e),
R_e(tau)=(V_e(tau)/(2 pi^2))^(1/3),
Theta_e(tau)=partial_tau log V_e(tau).
```

`R` is the volume radius; it does not assert pointwise spherical symmetry.
`Theta` is the volume-averaged expansion: exactly
`Theta = V^(-1) integral theta dvol_h`, where
`theta=(1/2)h^{ij} partial_tau h_ij`. It is not the unweighted average
of the local second expansion coefficient.

## 2. Conditional local foliation and its limits

Sufficient assumptions are: a smooth parameter family of smooth timelike
shell embeddings and induced Lorentzian metrics on a common neighborhood;
a smooth compact embedded spacelike initial section without boundary;
and the common future orientation. The metric and section must remain
nondegenerate for parameters in a neighborhood of zero. Smoothness is a
convenient sufficient assumption, not a sharp regularity claim.

Then there exist `e_*>0` and `tau_*>0` such that the maps `F_e` give a
smooth foliation for `|e|<e_*`, `|tau|<tau_*`, with no focal point or
intersection of normal geodesics in this neighborhood.

Proof: smooth geodesic ODE dependence supplies `F_e` on a uniform short
interval over the compact initial section. At `tau=0`, its differential
maps the time direction to `u_e` and the three section directions to
`dj_e`; these four vectors form a basis. Apply the inverse function
theorem to `(e,tau,z) -> (e,F_e(tau,z))`. Compactness gives a uniform
locally invertible interval. Global injectivity follows after shortening:
otherwise a sequence of colliding pairs with times tending to zero has a
subsequence whose limiting section points agree because `j_e` is an
embedding. The pairs then lie in one inverse-function neighborhood, a
contradiction. Gauss' lemma gives
`F_e^*gamma_e=-d tau^2+h_e(tau)` with positive `h_e` after shortening.
Thus the proper-time clock has no caustic on this local interval.

This proves existence of a sufficiently short interval, not that the whole
retained computational shell history satisfies an evaluated no-caustic
bound. A useful sufficient no-focal-point estimate is as follows. In a
parallel orthonormal screen along each normal geodesic, let `J(0)=I`,
`J'(0)=B_0`, and `J''+R_u J=0`, with
`||B_0||<=B`, `||R_u||<=K` uniformly. Integral Gronwall in its elementary
supremum form gives

```
sup_[0,T] ||J-I|| <= (BT+KT^2/2)/(1-KT^2/2).
```

Therefore `BT+KT^2<1` ensures no focal point on `[0,T]`. Global
injectivity must still be supplied by a normal-tube/injectivity argument;
absence of focal points alone does not rule out distinct geodesics crossing.
No `B,K,T` certificate is evaluated here.

The retained `q,k` coefficient solution does not establish an exact
nonlinear Einstein/scalar/shell family. For this specification the second
variation can equivalently be read as a jet functional: choose any smooth
kinematical realization with the prescribed jets and anchor, form the
observable locally, and take its second parameter derivative at zero.
The answer depends only on those jets, not on the third and higher
parameter terms of the realization. Kinematical realizability is not a
nonlinear physical existence theorem.

## 3. Full embedding and slicing variation

Let `Psi_e(tau,.)=E_e o F_e(tau,.)` be the complete ambient embedding of
the proper-time spatial section. Extend its variation locally by maps
`D_e` with `Psi_e=D_e o Psi_0`, and define the generators by the pullback
expansion

```
D_e^* = I + e L_(Z1) + e^2/2 (L_(Z2)+L_(Z1)^2) + O(e^3).
```

These `Z1,Z2` include the physical shell displacement, movement of the
initial section, normal-clock correction and angular relabeling. `Z2`
is this pullback generator, not the uncorrected coordinate acceleration
of an embedding. Restricting the resulting tensors to spatial slots,

```
h1 = Psi_0^*[q + L_(Z1) g0],
h2 = Psi_0^*[k + 2 L_(Z1) q + (L_(Z2)+L_(Z1)^2) g0].       (1)
```

All bending and slicing terms are present in (1). The expression is
independent of the auxiliary extensions when they describe the same
embedding family. In particular no choice of extension can remove the
physical displacement terms while retaining the same observable.

For `B=h0^(-1)h1`, `C=h0^(-1)h2`, differentiation of the determinant gives

```
dvol_h = dvol_h0 [1 + e tr(B)/2
                 + e^2/2 {tr(C)/2+(tr B)^2/4-tr(B^2)/2}] + O(e^3).
                                                               (2)
```

Write `v1=V1/V0`, `v2=V2/V0`, and let `<.>_0` be normalized background
spatial averaging at the specified `tau`. Then

```
v1 = <tr B/2>_0,
v2 = <tr C/2+(tr B)^2/4-tr(B^2)/2>_0,
R1/R0 = v1/3,
R2/R0 = v2/3 - 2 v1^2/9,
Theta1 = partial_tau v1,
Theta2 = partial_tau (v2-v1^2).                              (3)
```

The determinant contractions in (2) use the positive spatial metric
`h0`, not the four-dimensional Lorentzian trace of the shell perturbation.
The square of a mean cannot replace the mean of these quadratic products.

## 4. Explicit clock/lapse/shift terms in fixed-shell coordinates

In a representation with the physical shell and anchor fixed at `y=y_b`
and `s=0`, write the induced metric as

```
gamma_e = gamma0 + e m + e^2 n/2 + O(e^3),
m = E0^*q,       n = E0^*k.
```

For a moving shell, first use the shell-only embedding generators `e1,e2`
in the same pullback convention as (1):

```
m=E0^*[q+L_(e1)g0],
n=E0^*[k+2L_(e1)q+(L_(e2)+L_(e1)^2)g0].
```

These are four-dimensional induced shell tensors before clock slicing.
Apply the clock construction once to them; using the complete already
clock-sliced tensors from (1) as inputs would double count the clock.
Let the proper-time clock in the fixed shell chart be

```
T_e(s,z)=s+e alpha(s,z)+e^2 beta(s,z)/2+O(e^3),
gamma_e^(-1)(dT_e,dT_e)=-1,
T_e(0,z)=0.
```

The future normal is `u_e=-grad_(gamma_e) T_e`. Differentiating the
eikonal equation and its initial value yields, with all spatial norms
using `h0=A^2 kappa`,

```
alpha_s = -m_ss/2,                  alpha(0,z)=0,
beta_s  = -n_ss/2 - m_ss^2/4 + |m_s. + D alpha|^2_(h0),
beta(0,z)=0.                                               (4)
```

This explicitly includes the time-time perturbation and the time-angular
shift. In the shell ADM convention
`gamma=-N_lapse^2 ds^2+sigma_ij(dz^i+beta_shift^i ds)(dz^j+beta_shift^j ds)`,
`(N_lapse)_1=-m_ss/2` and
`(N_lapse)_2=-n_ss/2+|m_s.|^2-m_ss^2/4`.
The remaining terms in (4) implement the tilted geodesic clock and may
not be omitted by identifying coordinate time with proper time.

Instead of transporting angular labels with the geodesics, compute each
level set as a graph over the fixed angular chart:

```
s=S_e(tau,z)=tau+e a(tau,z)+e^2 b(tau,z)/2+O(e^3),
a=-alpha,                    b=-beta+2 alpha alpha_s,
```

where the right sides are evaluated at `s=tau`. The graph and geodesic
parameterizations have the same volume; their angular labels differ by
a diffeomorphism of the compact sphere. The `h1,h2` in (5) are in graph
labels and can differ from (1)'s geodesic-label tensors by this spatial
pullback. Direct graph pullback gives

```
h1_ij = m_ij + a partial_s h0_ij,
h2_ij = n_ij + 2a partial_s m_ij + b partial_s h0_ij
       + a^2 partial_s^2 h0_ij
       + 2(m_si D_j a + m_sj D_i a) - 2 D_i a D_j a.         (5)
```

The last two terms come from the shift and the timelike part of the
graph metric. The term with `partial_s m` is an evaluation correction;
the `b` and `a^2` terms are the second clock and background-evaluation
corrections. Formula (5) plus (2) accounts for every term through order
two in this fixed-shell representation.

## 5. Exact mean contraction for the retained quadrupole

The retained shell producer supplies, in its background orthonormal frame,

```
m_00=-2X Y,   m_0a=V D_aY,
m_ab=2[Z Y delta_ab + d D_aD_bY],
X=u+b_gauge,t,   Z=u+tanh(t)b_gauge,
V=cosh(t)d_t-b_gauge/cosh(t),
Y=n1^2-n2^2,             Delta_kappa Y=-8Y.
```

Here `d` and `b_gauge` are the retained first gauge amplitudes; `b` in
(5) is the second graph-clock coefficient and is unrelated. In the
`(s,z)` chart these components become

```
m_ss=-2X Y,     m_si=A V D_iY,
m_ij=2 A^2 [Z Y kappa_ij+d D_iD_jY].
```

Use probability measure `dmu=dvol_kappa/(2 pi^2)` and define

```
Cclock(s)=integral_0^s X(v) dv,     E(s)=Z(s)-H(s) Cclock(s),
Bclock(s)=<b(s,.)>.
```

All primes below are `partial_s=rho_b^(-1) partial_t`
(`=(rho_b ell)^(-1) partial_x`, with `ell=10^-6`).
Let `N0(s)` and `J0(s)` be the **total**, not auxiliary, mean amplitudes:
`<n_ss>=N0` and `<tr_(h0)n_ij>=3 J0`. The second-order operator's `J`
slot multiplies the angular metric, so its mean trace is `3 J0`, not
`J0`. The retained five-mean Gram matrix is `diag(1,2,1,3,1)` in
`(K0,L0,N0,J0,w0)`.

Equation (4) becomes

```
alpha=Cclock Y,
<beta>' = -N0/2 - X^2/6 + (4/3)(V+Cclock/A)^2,
Bclock = -<beta> + Cclock X/3,          <beta>(0)=0.          (6)
```

Equivalently,
`Bclock'=N0/2+X^2/2+Cclock X'/3-(4/3)(V+Cclock/A)^2`.
The first spatial matrix is
`B=2[E Y I+d Hess_kappa(Y)]`. Exact sphere contractions are

```
<Y>=0,  <Y^2>=1/6,  <|DY|^2>=4/3,  <|Hess Y|^2>=8,
tr(Hess Y)=-8Y,
< (tr B)^2/4-tr(B^2)/2 > = E^2/2-(8/3)E d-(16/3)d^2.       (7)
```

Therefore `v1=0` and the full second volume coefficient is

```
v2 = (3/2) J0 + 3H Bclock
   - (Cclock/3)[3Z'+6HZ-8d'-16Hd]
   + (Cclock^2/2)(H'+2H^2)
   - 8 Cclock V/(3A) - 4 Cclock^2/(3A^2)
   + E^2/2 - (8/3) E d - (16/3) d^2.                       (8)
```

In this retained quadrupole specialization,

```
V0=2 pi^2 A^3,   R0=A,   Theta0=3H,
V1=R1=Theta1=0,
V2=V0 v2,       R2=A v2/3,       Theta2=v2'.                 (9)
```

These are formal identities, not evaluated size or sign statements.
The nonmean second angular modes integrate to zero at this order because
the second coefficient enters linearly; every quadratic first mode is
retained. Second-order angular transport also does not need to be
numerically solved just to integrate the volume. `N0` is required through
the clock history and `J0` directly; their total evolution remains the
coupled five-mean problem. `K0=W_yy` cannot replace them, and `W` cannot
replace the registered total `U=U_hom+W`.

## 6. Gauge law, including the physical shell and anchor

For a gauge diffeomorphism `phi_e` with pullback generators `xi1,xi2`,

```
q' = q + L_(xi1) g0,
k' = k + 2 L_(xi1) q + (L_(xi2)+L_(xi1)^2)g0.             (10)
```

Transport the same physical shell and initial section with the inverse
map: `E'_e=phi_e^(-1) o E_e`. Then its intrinsic metric is unchanged,
so its normal geodesics and proper-time labels are the same physical
ones. The complete section embedding transforms as
`Psi'_e=phi_e^(-1) o Psi_e`, and exactly

```
(Psi'_e)^*(phi_e^*g_e)=Psi_e^*g_e.                          (11)
```

Consequently `h_e,V_e,R_e,Theta_e` and their parameter derivatives at
fixed anchored proper time are invariant. For the chosen pullback
generator convention the combined embedding generators obey

```
Z1' = Z1-xi1,
Z2' = Z2-xi2+[xi1,Z1].                                    (12)
```

Substituting (10),(12) into (1) cancels the first- and second-order gauge
terms, including the commutator. Shell relabelings change `h_e` by a
spatial pullback and leave the integral unchanged. A residual tangential
gauge shift need not fix the anchor's coordinate equation; the anchor must
be transformed as physical data. A genuine change to the initial section
or observer prescription is a change of observable and need not leave
the result unchanged.

The standard second-order gauge law is verified in Bruni, Matarrese,
Mollerach and Sonego, *Perturbations of spacetime: gauge transformations
and gauge invariance at second order and beyond*, CQG 14 (1997) 2585,
[arXiv:gr-qc/9609040](https://arxiv.org/abs/gr-qc/9609040), eqs. (1.2),(1.3).
The varying-embedding differentiation rule is verified in Mars,
*First and second order perturbations of hypersurfaces*, CQG 22 (2005)
3325, [arXiv:gr-qc/0507005](https://arxiv.org/abs/gr-qc/0507005), Lemma 1,
eqs. (7),(8); its section 7 distinguishes ambient and hypersurface gauge.
The PDFs/texts were independently retrieved by the literature task at
`../literature/sources/`; these citations support standard machinery,
not a novelty claim for this observable specification.

## 7. What a numerical observable certificate would still require

This note does not evaluate (6)--(9), propagate input intervals into them,
prove a positive or negative invariant response, or evaluate a no-caustic
time. A follow-on certificate must supply:

* Total `N0,J0` shell histories and the registered first shell traces
  with joint correlations/continuous error bounds sufficient for the
  products and clock integrals. Isolated point intervals do not suffice.
* A bounded normal-foliation interval and, for finite parameter claims,
  a justified family and its parameter range. The local existence proof
  above does not furnish a numerical range.
* An error budget for (8) and, for `Theta2`, its proper-time derivative.
  The retained first C3 information can support the required first-trace
  time derivatives; total `N0,J0` need C1 shell control. All factors
  `1/(rho_b ell)` must be retained. No differentiated source or fourth
  first-gauge derivative is introduced by (4)--(9).
* The same physical initial second prescription and the conditional
  subsidiary bridge for the exact total solution. Replacing total data
  by zero changes the observable being evaluated.

The total mean pilot's shell dependency history is
`0<=x<=5/4` for the proposed cone. At the background level this is
`0<=s<=rho_b ell(5/4)`. Restrict the observable to an interior compact
proper-time interval whose needed clock history lies in that shell
history and for which the local foliation is justified; small finite
parameter excursions also require a surrounding solution neighborhood.
No claimed local theorem automatically extends to the full cone or a
larger proper-time interval.

## 8. Retained source join and check scope

Read-only sources used for conventions are the forward plan
`boundary-response-20261023/followon/NEXT_SCIENTIFIC_TARGET.md`, the common
gauge capsule's `COMMON_HARMONIC_PHYSICAL_SHELL_FORCING.md` sections 2,4,5,
`BOUNDARY_COMPATIBILITY_AUDIT.md` sections 1,2, and
`EXPLICIT_INITIAL_HARMONIC_COMPLETION.md`, together with
`FULL_SECOND_RESPONSE_CONSTRUCTION.md` and the exact projection
`linear_operator.py`/`source_adapters.py`.

`check_observable_algebra.py` checks finite formal identities for the
eikonal clock, graph pullback, determinant, gauge cancellation, sphere
contractions, mean formula, and radius/log-volume coefficients with exact
symbolic arithmetic. Such checks supplement the derivation. They do not
authenticate numerical capsule closures, evaluate the actual response,
establish a nonlinear family, or validate a physical sign.
