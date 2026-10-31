# Independent proper-time screen formula review

The reviewed formula is conditional on the exact-family binding and on a
common invertible Gaussian tube. The floating coordinate fixtures below are
independent diagnostics; the all-point algebra is exact.

**Result:** the radial quadratic formula, all three first-Weyl combinations,
the complete relative Lorentz tilt, its reduction to the two coordinate
position variations, and the triangular geodesic-variation ODE agree with an
independent derivation. The independent coordinate curvature fixtures give
nonzero checks of the signs, second-parameter factors, and normalization.
This review supplies no evolved numerical sign.

## Gaussian screen and the two-position reduction

Let the first geodesic variation have off-screen part
`xi_off=Y(alpha partial_t+beta partial_y)`. At a fixed initial radial label,
the angular connecting vector has off-screen variation
`D_A Y(alpha partial_t+beta partial_y)`. Angular trajectory displacement and
its derivatives contribute only within the angular screen. Resolve these
vectors in the boosted background normal plane after dividing by R:

```
p=(rho cosh(chi) alpha-sinh(chi) beta)/R,
q=(-rho sinh(chi) alpha+cosh(chi) beta)/R.
```

Thus the relative orthonormal screen variation is
`delta e_a=v_a(p Uchi+q nchi)` plus an in-plane rotation. Gauss' lemma gives
`g(U,E_A)=0` for the initial normal geodesic congruence. Differentiating this
and `g(n,E_A)=0` in the common metric gives
`delta U=p v^a e_a`, `delta n=-q v^a e_a`, plus a radial boost. This fixes
both tilt signs without using a geodesic angular-velocity equation.

Indeed, orthonormality gives
`g(delta U,e_a)+g(U,p v_a U)=delta U_a-p v_a=0` and
`g(delta n,e_a)+g(n,q v_a n)=delta n_a+q v_a=0`.
The positive Gram inverse square root yields orthonormal axes. Its relative
in-plane map is orthogonal; its first rotation multiplies the complete zero
first magnetic tensor. Consequently no angular geodesic variation, angular
connecting-field variation, or radial-label derivative is required to
evaluate this selected second coefficient. They may still be used to prove
uniform tube invertibility.

For scalar first fields the complete first screen tensor vanishes for both
radial basis normals, by the trace/Hessian symmetric-curl identity. This
holds at every event and for every radial rapidity. Hence the complete
moving-event first tensor, first radial boost correction, and first in-plane
rotation vanish. The cancellation does not apply separately to `xi d C1`.
The metric-dependent baseline epsilon and screen-frame corrections remain
inside the complete time and radial quadratic terms.

## All-point curvature and Lorentz algebra

Set `H=theta/rho`, `m=M_t/(2rho)-E/R`, `j=M_y/2-C/R`. The radial Codazzi
calculation must use a timelike radial hypersurface. Its screen Gauss product
has the minus sign arising from its timelike normal. With

```
l_r=L_y/2-aA/2+HB,
beta_r=j+aM, beta_t=m+HM,
eta=R(C_t-theta C)/(2rho)-R E_y/2-B/2,
eta_r=-eta-B,
```

the four quadratic coefficients of
`Aax_ab=epsilon_(a^{cd}v_c Hess(Y)_b)d` are

```
-A beta_r/2+B beta_t-(C/R)(L-2M),
-(L/2-M)j,
M(-l_r+j),
m eta_r.
```

Their exact symbolic sum is

```
F_y=((-A+4M-L)M_y-2M L_y)/4
    +C(A-L)/(2R)+(E/R-M_t/(2rho))eta.
```

The spatial registered radial normal differs at first order from the full
normal to the radial hypersurface by `-B Y e0`. Its change in the second
screen tensor is proportional to the complete zero first time tensor, so
the displayed radial quadratic answer applies to the registered normal.

For the first Riemann tensor the Hessian coefficients are

```
g_R=-L/(2R^2)+H m-a j,
e_tR=-m_t/rho-2H m-D/(2R^2)+a j,
e_yR=-j_y-2a j-A/(2R^2),
k_R=eta/R^2-m_y-a m-H j.
```

Schouten cancellation gives the Weyl combinations
`At=g+et=g_R+e_tR`, `Ay=ey-g=e_yR-g_R`, `K=k=k_R`. Their expansions agree
with the producer's `At,Ay,K` expressions. Independently extracting these
three combinations from the full coordinate Weyl tensor also agrees.

Use the conventions
`C_abcd=g(delta_ac h_bd+delta_bd h_ac-delta_ad h_bc-delta_bc h_ad)`,
`C_aU bU=e h_ab`, `C_aU b n=k h_ab`, with isotropic terms omitted only
after symmetric contraction. Varying all four Weyl slots under the relative
Lorentz frame change yields

```
delta H=-[p(g+e)+q k] Aax.
```

The angular observer slot contributes `-p g Aax`; the two final angular
slots contribute `-p e Aax-q k Aax`. The first slot contributes zero:
its U part repeats an antisymmetric normal and its n part is `C_nUcd=0`.
The exact tensor calculation uses arbitrary symmetric h and arbitrary v,
so this is an all-point identity, independent of a numerical fit.

With `c=cosh(chi),s=sinh(chi)`, boosted combinations are
`Echi=c^2 At+s^2 Ay+2cs K`, `Kchi=cs(At+Ay)+(c^2+s^2)K`.
Substitution of p,q gives

```
p Echi+q Kchi=
 [rho alpha(c At+s K)+beta(s Ay+c K)]/R.
```

The second parameter derivative contributes twice the first-frame times
first-Weyl term. Since `Aax=2(Tplus-Tminus)`, the scalar tilt coefficient
is therefore `-4[p Echi+q Kchi]`. The old fixed probability normalization
extracts this coefficient exactly because each TT norm is 1/6 and their
cross product vanishes.

The resulting full formula is

```
h_geo=3(c U_Tx-s U_Tr)/(2ell rho^2 cosh(t))
      +(4/R)[c F_t+s F_y
              -rho alpha(c At+s K)-beta(s Ay+c K)].
```

The linear radial sign is independently checked using the full coordinate
curvature of an angular metric `R^2 u Y Hess(Y)`, whose curl is Aax. Its
first tensor is exactly
`[c u_t/rho+s u_y] Aax/(2R)`, with no `a u` or `H u` term. Thus the TT
radial derivative becomes `-s U_Tr/(ell rho)` under `dy=-rho ell dr`.
The use of the exact metric second derivative k introduces no extra half.

## Triangular geodesic ODE and a valid uniform bound

Write `X=rho xi^t`, `Y=xi^y`, `P=rho delta V^t`, `Q=delta V^y`, with
proper time as independent variable and c,s the background rapidity factors.
The background connection gives exactly

```
X'=a s X+P,                 Y'=Q,
P'=-a s P-2a c Q-2a_y c s Y-SP,
Q'=-2a c P-(2a^2+a_y)c^2 Y-SQ,

SP=aB-(c^2 D_t-s^2 A_tmetric)/(2rho)-c s D_y-s^2 B_y,
SQ=c^2[B_t/rho-a(A+D)-D_y/2]
   +c s[A_tmetric/rho-2aB]+s^2 A_ymetric/2.
```

Here `A_tmetric` is the first metric-amplitude derivative, not the invariant
At. Direct coordinate connection variations are

```
delta Gamma^t_tt=-D_t/2+a rho B,
delta Gamma^t_ty=-D_y/2,
delta Gamma^t_yy=-B_y/rho-aB/rho+A_tmetric/(2rho^2),
delta Gamma^y_tt=rho B_t-rho^2[a(A+D)+D_y/2],
delta Gamma^y_ty=A_tmetric/2-a rho B,
delta Gamma^y_yy=A_ymetric/2.
```

They prove SP,SQ algebraically. Angular variation cannot force this subsystem:
the background radial/time connection is angle independent and its mixed
radial/angular coefficients in these equations vanish. At a fixed label the
scalar forcing is proportional to Y. Initial values must come from the
authenticated initial field; they are zero only if the chosen launch lies
in the exact unperturbed initial collar.

For zero initial values and uniform caps on `[0,T]`, take

```
S>=max(|SP|,|SQ|),
M>=max(|as|+2|ac|,2|ac|),
K>=max(2|a_y c s|,|(2a^2+a_y)c^2|).
```

For `W=sup max(|P|,|Q|)`, integration of the ODE and
`|Y(t)|<=integral_0^t W` gives
`W<=TS+MTW+(KT^2/2)W`. Therefore

```
W<=TS/(1-MT-KT^2/2),
|Y|<=TW,
|X|<=TW/(1-sup|as| T),
```

provided both denominators are positive. This proves the proposed bound
without integrating a new four-state interval ODE. Each input cap must still
come from the complete authenticated actual-jet API over the declared tube,
with physical/scaled derivative factors and inherited errors retained.

## Reproduction, independence, and limits

`symbolic_audit.py` verifies the radial four-term identity, all first-Weyl
combinations, the complete Lorentz tensor contraction, hyperbolic Gaussian
reduction, and geodesic source formulas symbolically. It uses explicit
exceptions rather than Python assertions.

`coordinate_weyl_check.py` constructs the complete five-dimensional metric
in a nonnormal stereographic S3 chart. Angular jets are computed by exact
rational differentiation before conversion to binary64. Curvature uses
the complete metric two-jet, full inverse metric, first-kind Christoffels,
Riemann, Ricci, Schouten, and Weyl tensor, with degree-two parameter
convolution. It imports no production scientific module. It independently
checks eight general nonzero fixtures at two rational sphere points, with
every amplitude value and first/second physical time/radial jet populated,
nonzero a, nonzero theta, and both angular tilt types. It also checks six
linear normal-derivative fixtures. The comparison tolerance is 2e-10 times
the relevant scale; these are consistency witnesses, not interval bounds.

Named wrong formulas omit the radial quadratic term, reverse the radial
lapse sign, reverse the B sign in eta, omit the old epsilon/volume term,
drop mixed Weyl curvature, reverse radial screen tilt, or halve the tilt
second derivative. Each must have a substantial nonzero coordinate residual.
Normal and optimized runs must yield identical payloads.

The Gaussian congruence and its inverse are C2 in parameter on any common
regular compact tube when the metric family has the stated parameter C2 and
piecewise spacetime C4 regularity. This follows from finite-dimensional ODE
differentiation and the inverse function theorem after recording a uniform
Jacobian margin. It does not prove that an unspecified tube is covered or
free of caustics. The background shooting, proper-clock value, tube, exact
first-jet bounds, exact-family TT binding and new combined TT enclosure are
separate required inputs. No old sign is inherited by this review.
