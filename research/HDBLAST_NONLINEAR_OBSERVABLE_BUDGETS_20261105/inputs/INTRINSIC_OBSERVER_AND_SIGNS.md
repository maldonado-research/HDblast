# Conditional intrinsic shell budget and remaining mean data

1 November 2026. Retain the Einstein–scalar model, branch, complete actual first
field, and conditional exact-family construction in
[`HDBLAST_EXACT_FAMILY_TT_TRANSFER_20261031`](../HDBLAST_EXACT_FAMILY_TT_TRANSFER_20261031/DERIVATION.md).
All coefficients below are parameter derivatives: `g_e=g0+e q+e^2 k/2+o(e^2)`.
The analytic localized-constraint, shell-evolution, regularity and actual-field
premises remain explicit. Numerical authentication and finite algebra do not
replace them.

## 1. Physical shell, anchors and common domain

The shell has spacelike outward normal `n` and induced Lorentzian metric `gamma`.
Use the retained initial physical shell section and its intrinsic future unit
normal geodesics as anchors. Before their first crossing, their proper-time
pullback is `gamma_e=-d tau^2+h_e(tau)`. These are geodesics within the shell;
they differ from the ambient bulk geodesics in the preceding checkpoint.

Define normalized actual spatial averages, volume radius and expansion by

```
<f>_e = (integral f dvol_h_e)/(integral dvol_h_e),
V_e = integral dvol_h_e, Rv_e = (V_e/(2pi^2))^(1/3),
H_D,e = V_e'/(3V_e), theta_e=(1/2)tr(h_e^-1 h_e'),
s_e=(1/2)h_e'-theta_e h_e/3,
Q_D,e=(2/3)(<theta_e^2>_e-<theta_e>_e^2)-<|s_e|^2>_e,
C_e=<R3(h_e)>_e-6/Rv_e^2.
```

A prime on time amplitudes means `partial_tau`. The background is
`h0=A^2 kappa`, `A=rho_b cosh(t)`, `H=tanh(t)/rho_b`,
`tau=rho_b t=rho_b ell x`. The observation interval
`31/25<=x<=63/50` lies strictly inside the retained `x<=3/2` shell slab.
The whole clock history begins at zero and ends by `63/50`; no larger exact-family
slab is inferred from an older numerical domain. Compactness, background
injectivity and stability of the intrinsic geodesic ODE supply a common normal
tube after qualitative parameter shrinkage. No effective parameter radius is
computed. The family C2 regularity in the retained lower C4 topology supplies
the differentiations, curvature evaluation and normalized integrals used here.

Under joint transport of the configuration and anchors, these constructions
are covariant. A different physical initial section, background observer or
background averaging convention can define different quantities.

## 2. Completion-independent first-field geometry

Write the longitudinal shell metric first amplitude as `u(t)` and the shell
scalar first amplitude as `psi(t)Y`, where `Y=n0^2-n1^2` and `Delta Y=-8Y`.
The shell-preserving normal gauge amplitude vanishes at the boundary. The
identified harmonic first field therefore has this same scalar pullback.
The complete intrinsic clock and angular transport give

```
P(tau)=integral_0^tau u(s) ds, J'=P/A^2,
E=u-HP, L=P/A^2, T=E',
h1=2A^2[E Y kappa-J Hess(Y)],
theta1=(3T+8L)Y,
(s1)^i_j=-L[(Hess(Y))^i_j+(8/3)Y delta^i_j].
```

The tangential harmonic gauge terms cancel after both clock and angular
transport are retained. Exact sphere averages give

```
(Var theta)_2=(3T+8L)^2/3,
<|s|^2>_2=80L^2/9,
Q_D,2=(2/9)[(3T+8L)^2-40L^2].
```

For `C(h)=<R3>-6/Rv^2`, the first derivative at a round three-sphere is zero
for every metric variation. Thus `C2=D^2C[h1,h1]` contains no second metric
component. Diffeomorphism invariance and that criticality remove the Hessian
term in `h1`. Computing the remaining Hessian along `exp(2eEY)h0` gives

```
C0=C1=0, C2=10E^2/(3A^2),
<R3>_2=-4v2/A^2+10E^2/(3A^2), v2=V2/V0.
```

The generic identities, clock cancellation, exact conformal counterexample and
independent checks are in [the expansion proof](expansion/DERIVATION.md).
These use established averaging and normalized total scalar-curvature machinery;
no new general Buchert or Yamabe theorem is claimed.

## 3. Actual first shell trace recovered from its scalar source

The retained outward shell laws are `K_out=sigma gamma/6` and
`n(phi)=-sigma_phi/2`, with `Ric5=dphi tensor dphi+(2/3)Vbulk g`.
Contracting the standard hypersurface Gauss equation gives

```
R4=(D_gamma phi)^2+F(phi),
F=-sigma_phi^2/4+2Vbulk+sigma^2/3.
```

The background scalar is intrinsically constant. The first induced metric in
longitudinal representation is `2uY gamma0`; a tangential gauge changes no
first scalar curvature because `R4,0` is constant. Its conformal variation
therefore yields the actual first trace equation

```
u_xx+3ell tanh(ell x)u_x+ell^2[8sech^2(ell x)-4]u
 = ell^2 rho_b^2 F'_b psi(x)/6,
u(0)=u_x(0)=0.
```

Here `F'_b` is a scalar-field derivative. The retained boundary coefficient
gives `F'_b=2p Bcal=2p/(rho_b mu)<0`. This ODE is a consequence of the
existing first field equations and shell laws, not a new boundary prescription.

The numerical recovery uses frozen whole history cells, including all inherited
actual scalar errors even on `0<=x<=1`. For damping/mass bounds `a,b`, horizon
`X`, and `B=integral |forcing|`, the Volterra estimate is

```
sup|u_x| <= B/(1-aX-bX^2/2), sup|u|<=X sup|u_x|.
```

The report verifies the denominator is positive. Bounding the remaining lower
order feedback by `G` charges errors `xG`, `x^2G/2`, `x^3G/6` in the
first, second and third forcing primitives. Whole-cell nonnegative kernel
integrals are exact; no point quadrature or finite-iterate substitution is used.
The actual first-field boxes are retained and independently reconstructed.

The completed continuous enclosures imply, with deliberately loose outward
decimal endpoints,

```
3e-24 <= C2 <= 2e-22,
1e-9 <= Q_D,2 <= 2e-8.
```

The lower bound on the shear-square coefficient remains zero. No strict shear
claim is inferred from an interval containing zero. The older symmetric
global-error method is preserved as an unsuccessful sign test; the new trace
method has its own frozen campaign.

## 4. Scalar-gradient budget and angular causal character

Put `Q=(D_gamma phi)^2`, `K=(U phi)^2`, `S=|grad_h phi|^2`, with `U` the
intrinsic unit normal. Since `D_gamma phi0=0`, all three start at order two.
Their second derivatives are independent of second metric/scalar data, first
observer changes, and first/second event or measure variations at the fixed
background anchors. In proper-time notation,

```
<K>_2=psi'^2/3,
<S>_2=8psi^2/(3A^2),
<Q>_2=<S>_2-<K>_2.
```

This identity was already in the October 16 intrinsic-scalar calculation.
The new actual-field campaign gives continuous `<Q>_2<-4e17` and
`<S>_2/<K>_2<1e-12` on the stated interval, within the current exact-family
domain. It does not import the older second zero bracket lying beyond that
domain. Details and error inputs are in [the scalar component](gradient/README.md).

For uniform round angular measure, set `Z=n0^2+n1^2` and
`Y=Z cos(2alpha)`. Then `Z` is uniform on `[0,1]`, the angle is uniform and
independent, and `|DY|^2=4(Z-Y^2)`. If
`R=<S>_2/<K>_2`, `r=R/(R+2)`, the negative coefficient fraction is

```
(2/pi)[arccos(sqrt(r))-sqrt(r(1-r))].
```

Its complementary fraction is at most `4sqrt(r)/pi<4sqrt(r)/3`.
The rational inequalities `R<1e-12` and
`r<1/(2e12)<[3/(4e6)]^2` prove that fraction is strictly less than `1e-6`.
The geometry component checks the angular distribution and degenerate cases.

For every fixed admitted exact family, uniform C2 Taylor convergence on the
compact time/angular domain preserves a compact negative-coefficient set of
mass greater than `1-1e-6` for all sufficiently small nonzero parameters.
Normalized actual proper volume converges uniformly to the background measure,
so the same strict fraction holds for actual proper volume. There is no effective
or completion-uniform amplitude interval. A temporal scalar gradient is not a
radiation fluid, a particle-production rate or a force.

## 5. The remaining second scalar mean is essential

The exact geometric trace average is

```
6(H_D'+2H_D^2)=<R4>-<R3>+Q_D.
```

Let `M2=(<phi>_e)_2`, with scalar, clock and measure all pulled back to the
same intrinsic proper-time slices. Then the full second volume equation is

```
v2''+4H v2'-2v2/A^2
 = (F'_b/2)M2 + F''_b psi^2/12-psi'^2/6+4psi^2/(3A^2)
   -5E^2/(3A^2)+Q_D,2/2.
```

The coefficient `F'_b` is nonzero. The new strict first-only bounds therefore
do not close this equation or establish a sign of volume acceleration.
An affine initial canonical scalar does not set its second postcontact mean
history to zero. The selected family need not share all second scalar data
with an older formal five-mean solution.

The [initial mean calculation](weighted-mean/README.md) spells out the radial
weighted equations and the canonical conversion. In its fixed reconstruction
gauge, write the averaged second initial spatial metric as
`diag(A2,B2,B2,B2)`. It proves

```
integral rho^2(A2+2B2)dy
 = (2/3) integral rho^4 QH,0 dy
 <= -(49/150) integral rho^4 v_initial^2 dy < 0.
```

All proper-radial integrals have increasing outward `y`; `dy=-rho dz` is
converted with positive support length at most `131ell`. The proof completes
squares, uses `Vbulk''>3`, and bounds the only positive term by the compact
first-constraint primitive. It holds for every admissible weight choice.
The defining function and weight parameters remain uninstantiated, so it
supplies no unique numerical initial mean profile. Its homogeneous kernel
and the freedom in second data are made explicit. This initial fixed-gauge
moment is not shell energy, a force or a cosmological observable.

## 6. Qualitative exact-family consequence and scope

For each retained C2 family, compactness and the strict continuous coefficient
bounds give an uncomputed `delta>0` such that throughout the stated proper-time
interval, for `0<|e|<delta`,

```
C_e >= (3e-24/4)e^2 > 0,
Q_D,e >= (1e-9/4)e^2 > 0,
<Q>_e <= -(4e17/4)e^2 < 0.
```

These follow by the integral Taylor formula; all three zeroth and first
derivatives vanish. The same interval may be reduced for the proper-volume
angular fraction. A shared first tangent supplies no uniform nonlinear radius
over arbitrary completions, and none is claimed. No sign at a specified
nonzero amplitude or conversion to an empirical physical scale is established.

The finite exact checks and independently written rational reconstruction
support the calculation under shared source premises. They are not external
peer review or a proof-assistant formalization. The generic geometry is
standard; the claimed progress is the model-specific trace recovery, continuous
intrinsic enclosures, exact-family applicability and initial-mean constraint.
The quantitative second-order mean response for this selected construction,
created hot matter and higher-dimensional Big Bang origin remain unresolved.
