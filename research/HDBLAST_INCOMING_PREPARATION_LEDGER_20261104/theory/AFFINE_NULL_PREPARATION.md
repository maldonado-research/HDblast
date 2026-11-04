# A constraint-forced second expansion on an initial null cut

4 November 2026. Conditional partial characteristic construction. This solves
one exact null Einstein constraint through second order after declaring free
state data. It does not construct the full past Einstein--scalar solution,
select a physical state from the existing action, or predict the current shell
volume/expansion coefficients.

## 1. Declared past cut, generator and fields

For orientation, the actual radial background is

```
g0=rho(z)^2[-dt^2+dz^2+cosh(t)^2 kappa],
phi0=phi0(z), dy=-rho dz.
```

One possible past initial characteristic cut is the outgoing surface
`Cminus={t-z=tminus}`, with `tminus<-2ell`, restricted to the actual positive
collar and a segment preceding the present t=0 slice. Its corner at the shell
is `(tminus,z=0)`. To encompass the current packet, a strictly larger causal
margin is needed, not merely equality tminus=-2ell. An incoming pulse crossing
this outgoing cut can subsequently propagate toward the shell. Prescribing
it is a new preparation-state datum, not already supplied by current Cauchy
values. The other radial null orientation can be treated similarly.

On this background a future affine generator has the form
`k0=c rho(z)^(-2)(partial_t+partial_z)`, c>0, with c fixed by the corner clock
normalization. The spherical cross-section radius is
`R0(lambda)=rho(z(lambda))cosh(t(lambda))`, not merely rho. In particular

```
vartheta0=3R0'/R0=3c rho^-2[tanh(t)-rho a],
phi0'=-c p/rho,  p=phi0_y, a=rho_y/rho.
```

Primes in this note mean affine lambda derivatives, not y,r or proper-shell
clock derivatives. These formulas satisfy the background Raychaudhuri identity exactly. The
retained static equations give `rho_z=-a rho^2` and
`a_y=-1/rho^2-p^2/3`, hence `a_z=1/rho+rho p^2/3`. Put
`A=rho a`, `T0=tanh t`, `B2=rho^2 p^2`, and `Dplus=partial_t+partial_z`.
Then

```
Dplus(T0-A)=A^2-T0^2-B2/3,
Dplus(rho^-2)=2A rho^-2,
vartheta0'=3c^2 rho^-4[2A(T0-A)+A^2-T0^2-B2/3],
vartheta0^2/3=3c^2 rho^-4(T0-A)^2.
```

Adding the last two expressions cancels all A,T0 terms and leaves
`-c^2 rho^-4 B2=-(Phi0')^2`, exactly as required. These are algebraic
identities, not new numerical past-domain enclosures.

For a smooth one-parameter family, choose physical null cuts C_epsilon and
pull them back to fixed `(lambda,omega)` labels so that k=partial_lambda is
affine for every family member and its corner normalization is specified.
The cut embedding and this normalization are part of the proposed state
prescription. Let q_epsilon,AB be the positive three-metric of its sections,
Phi_epsilon the **pulled-back scalar**, and

```
B_AB=(1/2)partial_lambda q_AB,
vartheta=q^AB B_AB,
sigma_AB=B_AB-(vartheta/3)q_AB.
```

The null cut is hypersurface orthogonal, so twist vanishes. At the spherical
background, `q0=R0^2 kappa`, `sigma0=0`. Write jth actual parameter derivatives
as vartheta_j, sigma_j, Phi_j (epsilon^2/2 convention). Phi_j includes changes
of the cut embedding and generator through the pullback; it must not be
replaced by the old harmonic-coordinate scalar coefficient without that
transport calculation.

## 2. Exact second Raychaudhuri constraint

For the retained bulk action `R/2-|dphi|^2/2-Vbulk(phi)` in five dimensions,

```
Ric_AB=partial_A phi partial_B phi+(2/3)Vbulk g_AB,
Ric(k,k)=(Phi')^2.
```

The potential term drops only because k is null. The exact affine optical
constraint is therefore

```
vartheta'+vartheta^2/3+|sigma|_q^2+(Phi')^2=0.             (1)
```

Differentiating at the background gives

```
vartheta1'+(2vartheta0/3)vartheta1=-2Phi0'Phi1',
vartheta2'+(2vartheta0/3)vartheta2
 =-(2/3)vartheta1^2-2|sigma1|_q0^2
   -2(Phi1')^2-2Phi0'Phi2'.                              (2)
```

No shear-second term survives because sigma0=0. Varying the inverse section
metric in the shear norm produces no additional second term for the same
reason. The factor2 in each square follows from actual second parameter
variation, not a coefficient lacking the Taylor factorial. All generator
variation terms are already incorporated by the fixed affine pullback; with
an arbitrary changing nonaffine generator, equation(2) would require extra
terms.

Here is a concrete **candidate state condition**, not a derived law:

```
Phi2=0 on Cminus,
vartheta2(lambda0,omega)=0.                              (3)
```

The first condition fixes free second scalar data on the declared null cut.
It is not equivalent to zero scalar second data on the present Cauchy slice.
A pulse flat near the corner can be paired with background corner data for
this equation, but the other characteristic constraints and physical shell
corner conditions still have to be solved. A separately chosen second free
conformal section metric may affect those equations; it does not enter (2)
through sigma2 at this order.

## 3. Explicit forced solution and its sign

The integrating factor for (2) is exactly `R0^2`, since
`2vartheta0/3=2R0'/R0`. In the general case,

```
R0(lambda)^2 vartheta2(lambda)
 =R0(lambda0)^2 vartheta2(lambda0)
  -2 integral_lambda0^lambda R0(s)^2
    [vartheta1(s)^2/3+|sigma1(s)|^2+(Phi1'(s))^2
                                  +Phi0'(s)Phi2'(s)] ds.  (4)
```

Every angular label is retained in this identity. Under (3), put
`S1=vartheta1^2/3+|sigma1|^2+(Phi1')^2>=0`; then

```
vartheta2(lambda,omega)
 =-2 R0(lambda)^(-2) integral_lambda0^lambda R0(s)^2 S1(s,omega) ds
 <=0.                                                    (5)
```

The inequality is strict downstream of any nonzero integral of S1. Angular
averaging also preserves it, because R0 is independent of omega and the
background sphere measure is fixed. If Phi2 is not fixed as in (3), the
last term of (4) has no prescribed sign and this conclusion cannot be made.
Even with (3), first variation vartheta1 is constrained by the first equation
of (2); it is not an independently assignable function unrelated to Phi1.

Define the transported log area-density difference
`L_epsilon=log(dmu_qepsilon/dmu_q0)`. Then
`partial_lambda L_j=vartheta_j`. If additionally `L2(lambda0)=0`,

```
L2(lambda,omega)
 =-2 integral_lambda0^lambda R0(s)^2 S1(s,omega)
                    [integral_s^lambda R0(v)^(-2)dv] ds <=0. (6)
```

Thus a vanishing homogeneous second expansion at the corner does not imply
vanishing second null geometry: the first pulse forces the double-integral
term. Equation(6) is an explicit source functional for this one constraint.

Actual area has a different second variation. At a fixed affine section,

```
A2(lambda)/A0(lambda)=average_kappa[L2(lambda,omega)+L1(lambda,omega)^2]. (7)
```

Consequently a negative second optical expansion or log area does not by
itself prove negative actual area. Nor is vartheta2 the current proper-time
shell expansion Theta2. These null sections are not the scalar-level sections
inside the physical Cauchy slice used to define D. Their clock, embedding,
scalar-label displacement and full canonical evolution must be transported
before any comparison with the H1 preparation bound.

## 4. What remains unsolved

The actual first pulse is presently specified as mixed Cauchy data, not as a
past-cut free state. Its compatible Phi1 and sigma1 on Cminus require an
actual past development and cut transport or an independently supplied
incoming state proven to develop into that packet. Conditions(3) are new
second state assumptions. Equations(4)--(6) give their necessary and unique
Raychaudhuri response after the corner constants are set; they do not solve
all null Einstein constraints, the scalar wave equation, the physical shell
junction law or a nonlinear characteristic initial/boundary problem.

A physical source is not invented here. With background past data, no incoming
state and no added source, uniqueness cannot create the nonzero pulse. A
model must specify a nontrivial state or generating process, its coupling and
constraint/conservation laws. Only after that past problem is well posed and
computed can the present second state be tested against the retained h-only
sector and its complete-increment H1 norm. The affine focusing sign alone
provides no H1 radius, finite nonlinear parameter window, radiation production
or cosmological-origin evidence.


## Retained source links

- [Actual first packet and initial master](../../HDBLAST_COMMON_GAUGE_FORCING_20261016/inputs/HIGHER_JETS.md).
- [Physical first Cauchy data and nonzero second constraint witness](../../HDBLAST_COMMON_GAUGE_FORCING_20261016/inputs/EXPLICIT_COMPACT_SECOND_INITIAL_JET.md).
- [Exact harmonic lapse/shift velocities](../../HDBLAST_COMMON_GAUGE_FORCING_20261016/inputs/EXPLICIT_INITIAL_HARMONIC_COMPLETION.md).
- [Canonical action/momentum normalization](../../HDBLAST_SHELL_NOETHER_20261013/manuscript/SHELL_NOETHER.tex).
- [Five mean physical boundary rows](../../HDBLAST_PREPARED_SHELL_VOLUME_20261101/gauge/MEAN_INITIAL_GAUGE.md).
- [Conditional continuous operator premises](../../HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023/analytic/EXPLICIT_PREMISES.json).
- [Retained coupled preparation comparison and proper scalar meaning](../../HDBLAST_SHELL_PREPARATION_DEPENDENCE_20261104/DERIVATION.md).
- [Earlier direct H1 energy theorem](../../HDBLAST_H1_PREPARATION_RESPONSE_20261104/energy/PROOF.md).
