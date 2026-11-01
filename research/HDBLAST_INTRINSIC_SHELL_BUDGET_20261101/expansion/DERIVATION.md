# Anchored shell expansion and intrinsic curvature defect

1 November 2026. Portable source-bound geometry component, based on commit
`b51908d7e62e7ec80cef82438edd00923ff94723`. This component supplies an analytic derivation and finite exact checks.
All parameter coefficients are derivatives, not half-coefficients.

This applies standard Gaussian foliation, averaged Raychaudhuri and
normalized total scalar-curvature identities to the **selected exact family**
of `HDBLAST_EXACT_FAMILY_TT_TRANSFER_20261031/DERIVATION.md`, conditional on
that family's stated analytic premises. Its full first field agrees with the
registered first field. Its scalar second data need not agree with an older
formal completion. Nothing here replaces those second data with the older
completion or derives a sign of the complete volume response.

## 1. Clock, domain and notation

Use the physical initial shell section as anchor and its intrinsic future
unit-normal geodesics, as in the registered shell observable. On an interval
where their exponential map is injective,

```
gamma_e = -dτ² + h_e(τ),   K_e = (1/2)∂τ h_e,
θ_e = tr_h K_e,           s_e = K_e − θ_e h_e/3,
S_e = |s_e|_h²,           V_e = ∫S³ dvol_h,
Rv_e = (V_e/(2π²))^(1/3), H_D,e = V_e'/(3V_e),
<f>_e = V_e^−1 ∫S³ f dvol_h.
```

Primes on time-dependent amplitudes below mean `∂τ`. A subscript `2`
means the second parameter derivative at zero. `s_e` is shear; it is not
the tension `sigma(phi)`. There is **no factor 1/2** in `S_e`.

The background has `A=ρ_b cosh(τ/ρ_b)`, `H=A'/A`, `h0=A²κ`,
`R3,0=6/A²` and `R4,0=12/ρ_b²`. With `t=τ/ρ_b=ell*x`,
`∂τ=(ρ_b ell)^−1∂x`, `H=tanh(t)/ρ_b` and
`H'=sech²(t)/ρ_b²`. The volume radius is a volume statistic, not an
assumption that `h_e` is round or that `gamma_e` is FRW.

The old observable proved a sufficiently short local normal tube. For the
conditional exact family a fixed compact background interval inside its
common shell solution slab also admits this foliation after qualitative
parameter shrinkage: the background normal map is an embedding with no
focal point there, its derivative stays nonsingular under C1-small changes,
and compactness plus local invertibility prevent new coincident points.
Smooth geodesic dependence and a strict time margin keep the perturbed
maps inside the solution slab. This supplies no numerical parameter radius
or caustic estimate. In particular any numerical first-field history is
used only on a compact proper-time interval with such a surrounding margin.

## 2. Exact intrinsic averaged balances

In this curvature and future-second-form convention, direct Gaussian
connection/Gauss contraction gives

```
G4(u,u) = R3/2 + θ²/3 − S/2,
Ric4(u,u) = −θ' − θ²/3 − S,
R4 = R3 + 2θ' + 4θ²/3 + S.
```

For any scalar `f`, differentiation of the normalized spatial integral is

```
<f>' = <f'> + <f θ> − <f><θ>.
```

Consequently, with `Q_D = (2/3)(<θ²>−<θ>²)−<S>`,

```
3 H_D² = <G4(u,u)> − <R3>/2 − Q_D/2,
3(H_D' + H_D²) = −<Ric4(u,u)> + Q_D,
6(H_D' + 2H_D²) = <R4> − <R3> + Q_D.             (1)
```

These are the established kinematical averaging identities often used in
Buchert's formalism. No dust flow or four-dimensional matter equation of
state is assumed: these observers are shell-intrinsic geodesics, not a
specified matter flow. The last equation is the trace combination of the
first two, not an independent closing evolution law.

Set

```
F(phi) = −sigma_phi(phi)²/4 + 2Vbulk(phi) + sigma(phi)²/3.
```

The self-contained standard Gauss calculation in `GAUSS_IDENTITY_STANDARD_APPLICATION.md` supplies

```
R4 = −(u phi)² + |grad_h phi|² + F(phi),
G4(u,u) = (u phi)²/4 + 5|grad_h phi|²/12 + F(phi)/4 − E(u,u),
Ric4(u,u) = 3(u phi)²/4 − |grad_h phi|²/12 − F(phi)/4 − E(u,u). (2)
```

Here `E_ab=C5(n,e_a,n,e_b)` is the electric Weyl projection at the
shell. It cancels from the trace equation but is not determined by the
published bulk magnetic screen coefficient. Equations (1)-(2) are induced
geometry of the existing bulk scalar/tension model, not identification of
the right sides with actual four-dimensional fluid matter.

## 3. The first Gaussian field, including angular transport

Let `u(τ)` be the physical longitudinal metric amplitude at the shell and
`v(τ)` the first shell scalar amplitude, so `phi1=vY`,
`Y=n0²−n1²`, `Delta_κ Y=−8Y`. Use different letters for the unit vector
`u_e` and scalar amplitude `u`; all formulas below use the scalar amplitude.
The registered harmonic shell metric has first coefficients

```
m_ττ = −2(u+b_t)Y,
m_τi = A(cosh(t)d_t−b/cosh(t))D_iY,
m_ij = 2A²[(u+tanh(t)b)Yκ_ij + d D_iD_jY].
```

The gauge functions `b,d` and physical data vanish initially on the shell
collar. Define the histories

```
P(τ)=∫₀^τ u(s)ds,    J(τ)=∫₀^τ P(s)/A(s)² ds,
E=u−HP,              L=P/A²,    T=E'=u'−Hu−H'P.
```

The first eikonal correction is `alpha=(P+ρ_b b)Y`. The first normal
angular velocity in the original angular chart is
`−(d'+P/A²)grad_κ Y`. Integrating it gives angular displacement
`−(d+J)grad_κ Y`. Combining this displacement with the graph-slice metric
is essential: graph labels alone are not geodesic labels. The complete
Gaussian first spatial metric is therefore

```
h1 = 2A²[E Y κ − J Hess_κ Y].                     (3)
```

All harmonic `b,d` terms cancel. `phi1=vY` is unchanged because the
background scalar is constant on the shell. From `δ(h^−1 K)=∂τ(h0^−1h1)/2`,

```
θ1 = (3T+8L)Y,
(s1)^i_j = −L[(Hess_κ Y)^i_j + (8/3)Y δ^i_j].     (4)
```

The exact sphere moments are
`<Y²>=1/6`, `<|DY|²>=4/3`, `<|Hess Y|²>=8`, and
`<|Hess Y+(8/3)Yκ|²>=40/9`. Thus

```
(Var θ)_2 = (3T+8L)²/3,
<S>_2 = 80L²/9 >= 0,
Q_D,2 = (2/9)[(3T+8L)²−40L²].                    (5)
```

These coefficients use only the full first field and its clock history;
any second angular TT metric coefficient enters shear linearly at order
two and therefore cannot enter its squared norm until higher order. The
scalar-gradient/Q algebra was already present in October 16 work. The
application to the present selected exact family and anchored clock is the
relevant join here, not a new universal Q formula. `Q_D,2` can have either
sign; positivity of a shear square does not prove positive acceleration.

## 4. Curvature defect: cancellation of every second metric component

Define the intrinsic curvature defect at each fixed proper time

```
C(h) = <R3>_h − 6/Rv(h)².
```

For `h0=A²κ`, let `B=h0^−1 δh`. The first variation of the
Einstein-Hilbert integral `I(h)=∫R3 dvol_h` on the closed sphere is
`δI=−∫<G3,δh>dvol = A^−2∫tr(B)dvol`, because
`G3,0=−h0/A²`. Since `δV=(1/2)∫tr(B)dvol`,

```
δ<R3> = −(2/A²)<tr B>,
δ(6/Rv²) = −(2/A²)<tr B>,   hence dC|h0 = 0.       (6)
```

For an arbitrary C2 curve, the chain rule is
`C''=D²C[h1,h1]+DC[h2]`. Equation (6) kills **all** of `h2`, including
the exact family's still unknown scalar monopole. Diffeomorphism invariance
and criticality also imply `D²C[L_Xh0,k]=0` for every `X,k`: differentiate
`DC_h[L_Xh]=0` in direction `k` and use `DC_h0=0`.
The Hessian part of (3) is `L_{−J grad Y}h0`, so it contributes nothing.

It suffices to calculate along `h_e=exp(2eEY)h0`, which has the remaining
first tangent. The dimension-three conformal formula and integration by
parts on the closed sphere give

```
V_e/V0 = <exp(3eEY)>_0,
I_e/V0 = A^−2 <exp(eEY)[6+2e²E²|DY|²]>_0.
```

Using the exact moments above,

```
V2/V0=3E²/2,     I2/V0=19E²/(3A²),
<R3>_2=−8E²/(3A²),   (6/Rv²)_2=−6E²/A².
```

Thus for the actual selected family, regardless of its second data,

```
C0=C1=0,      C2 = 10 E²/(3A²) >= 0,              (7)
<R3>_2 = −4 v2/A² + 10E²/(3A²),  v2=V2/V0.       (8)
```

The conformal family is an independent kinematical counterexample to exact
FRW substitution: at `E!=0` it has positive second curvature defect. It is
not asserted to solve the shell/bulk field equations. For the current exact
family, a certified nonzero `E` at an observation proves `C2>0`; with C2
parameter regularity it also gives qualitative `C_e>0` for both sufficiently
small nonzero parameter signs. It gives no effective parameter radius.
This uses the well-known criticality/Hessian of normalized total scalar
curvature at a round Einstein metric; no new Yamabe theorem is claimed.

The checker also supplies an exact finite counterexample, rather than only
a nonzero Taylor coefficient: take `h=(1+Y/10)^4 h0`. Its conformal factor
is at least `9/10`. Set `D=<(1+Y/10)^6>` and
`N=<6(1+Y/10)^2+8|DY|²/100>`. The three-dimensional conformal formula gives
`C=(N−6D^(1/3))/(A²D)`. Exact sphere quadrature verifies the rational
inequality `N³−216D>0`, certifying a strictly positive defect for this
kinematical metric. Its exact fractions are recorded in both reports.

## 5. Exact second volume balance and the remaining physical unknown

The registered first quadrupole has `V1=0` and hence
`H_D,2=v2'/3`. Put `M2=(<phi>_e)_2`, the **proper-time spatial averaged**
scalar second derivative. It includes measure/clock correlations:

```
M2=<phi2>_0 + 2<vY w1>_0,
w1=(1/2)tr(h0^−1h1)=(3E+8J)Y,
M2=<phi2>_0 + v(3E+8J)/3.
```

Since `D phi0=0` on the shell, (2) gives

```
<R4>_2 = F'_0 M2 + F''_0 v²/6 − v'²/3 + 8v²/(3A²). (9)
```

The prime on `F` here denotes a scalar-field derivative at `phi_b`, not a
proper-time derivative. Inserting (8)-(9) into (1) yields

```
v2'' + 4H v2' − 2v2/A²
 = (F'_0/2) M2 + F''_0 v²/12 − v'²/6 + 4v²/(3A²)
   − 5E²/(3A²) + Q_D,2/2.                        (10)
```

This is a volume/trace consistency equation, not a closed scalar
evolution: the second scalar mean `M2` remains. It cannot be set to zero
because the initial scalar canonical curve is affine. That initial
condition does not impose zero scalar second history after interaction.
Nor can an older five-mean solution replace the selected exact family's
second solution without a full data match.

The unknown is not multiplied by a vanishing coefficient in this model.
At the background shell use `p=phi0_y=−sigma_phi/2`, `sigma=6a` and the
retained dynamic coefficient `Bcal=Vbulk'/p−4a+sigma_phiphi/2`. Directly,

```
F'_0 = −sigma_phi sigma_phiphi/2 +2Vbulk'+2sigma sigma_phi/3
     = 2p Bcal = 2p/(ρ_b mu).
```

The current branch has `p<0`, `rho_b>0` and the retained dynamic boundary
certificate has `mu>0`; thus `F'_0<0`, in particular `F'_0!=0`. A complete
volume acceleration sign still needs `M2` (or an equivalent coupled
second-order physical solution). Equations (7) and (5), plus certified first
histories, are the available progress without that solve.

## 6. Checks, source scope and limits

`check_shell_expansion.py` independently contracts an algebraic Gaussian
4-curvature, integrates ambient-polynomial sphere moments, computes the
normalized conformal functional along two different second jets, checks
gauge/clock cancellation, and checks the scalar volume equation. Its
deliberate wrong-sign/coefficient controls must be rejected. These checks
are finite exact algebra, not PDE solutions or specialist proof review.
`SOURCE_PINS.json` authenticates repository-relative mathematical inputs and the local proof components. No private-directory input is required.
Normal and optimized reports are intended to agree exactly.

The exact-family existence/regularity premises remain conditional. No new
four-dimensional fluid, radiation energy, material clock, nonlinear
cosmological origin, finite-amplitude quantitative sign or public theorem
is claimed.

## 7. First intrinsic trace ODE for a separate numerical pilot

The longitudinal representation has the full induced Lorentzian first
metric `m=2uY gamma0`; its harmonic representative differs by a
shell-preserving tangential Lie derivative. Since the background intrinsic
scalar curvature is constant, its first scalar-curvature variation is
unchanged by that tangential gauge. The standard conformal first-variation
formula in four dimensions is

```
R4,1 = −2uY R4,0 − 6 Box_gamma0(uY)
     = 6[u''+3Hu'+(8/A²−4/ρ_b²)u]Y.
```

The shell Gauss trace (2) has `R4,1=F'_0 vY`; the gradient-square term has
zero first variation. Thus the **actual first shell amplitude**, for the
same exact-family tangent, satisfies

```
u''+3Hu'+(8/A²−4/ρ_b²)u = F'_0 v/6,
u(0)=u'(0)=0.                                      (11)
```

With `t=τ/ρ_b`, `x=t/ell`, this is exactly

```
u_xx + 3ell tanh(ell*x)u_x
 + ell²[8sech²(ell*x)−4]u = ell²ρ_b² F'_0 v/6.       (12)
```

No replacement of physical proper time by perturbed coordinate time is
being made: `u` is a first coefficient on the background, so conversion
uses background `τ=ρ_b t`; equations (3)-(4) still implement the full
anchored clock when constructing the observable. This linear equation is
a consequence of the actual first bulk/shell equations, not a new boundary
condition added to the master problem. The retained master equation plus
dynamic scalar shell trace prove those equations in the pinned
`FIRST_SHELL_MASTER_BINDING_NOTE.md`.

The source `v` is precisely the first scalar returned at `r=0` by the
current actual-jet reconstruction: the shell-preserving normal gauge
amplitude vanishes there, leaving `chi_harmonic=vY`. Strict `F'_0<0` is
proved in section 5. It is therefore possible to recover an asymmetric
actual `u` enclosure by a validated causal ODE/Volterra calculation from
the scalar source, even when the old Dirichlet comparison for the master
has identically zero shell value. Strict signs still require a completed
source integration and a bound for the homogeneous coefficient feedback;
they do not follow from the displayed ODE alone. The separate numerical
pilot must retain its own frozen history mesh, source errors and controls.
