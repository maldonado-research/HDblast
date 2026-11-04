# Incoming-state constraints and a future quadratic-source ledger

4 November 2026. This is a conditional continuation of the retained five-dimensional
Einstein–scalar shell model. It does not add a physical blast source, solve a past
characteristic problem, or change the earlier source capsules.

## 1. Fixed model and normalization

Keep the selected background, first quadrupole packet, actual canonical data,
physical shell, harmonic conditioning and anchored observer of the
[prepared-family checkpoint](../HDBLAST_PREPARED_SHELL_VOLUME_20261101/).
The bulk action has normalization `R/2-|dphi|^2/2-V(phi)`. Two identical interior
copies have one-copy outward laws `K_ab=sigma*gamma_ab/6` and
`n_out(phi)=-sigma_phi/2`. Parameter coefficients below are actual derivatives:
`X_eps=X0+eps*X1+(eps^2/2)*X2+o(eps^2)`.

Write `p=phi0_y`, `a=rho_y/rho`, with `dy/dz=-rho` in the static background
`g0=dy^2+rho^2[-dt^2+cosh(t)^2*kappa]`.

The conditional future operator assumptions, mean projection, source enclosures,
finite propagation and uniqueness are retained. Authenticating their bytes does
not reprove them. The current source-state join checks the full physical and
canonical first data and the full harmonic one-jet. The October 31 native-Cauchy
velocity update changed a second-response approximation; it did not replace the
registered actual first packet or its quadratic sources.

## 2. Registered first packet and characteristic content

The normalized scalar master obeys, on the current initial section,

```
w(0,z)=q_ell(z), w_t(0,z)=exp(-z)q_ell(z),
ell=10^-6, supp(q_ell) subset [ell,2ell].
```

The nonzero smooth packet is flat at both endpoints. In normalized coordinates
`x=t/ell,r=z/ell`, this is `U0=Q(r)`, `U0_x=ell*exp(-ell*r)*Q(r)`.
Its time derivative is nonzero; a time-symmetric background does not make the
perturbation time symmetric.

Define the normalized-master principal characteristics `Jplus=w_t+w_z` and
`Jminus=w_t-w_z`. The first travels toward the shell at z=0, the second away.
Either pure-one-way condition would impose a homogeneous first-order ODE on
q with a zero endpoint. Uniqueness would give q identically zero. At the unit
packet's center q=1,q_z=0, both characteristics are positive.

For their conventional quadratic weights `Eplus/minus=(1/4)integral Jplus/minus^2 dz`,
compact endpoint integration by parts gives

```
Eplus-Eminus=(1/2)integral exp(-z)q_ell^2 dz >0,
Eplus+Eminus=(1/2)integral[(q_ell,z)^2+exp(-2z)q_ell^2] dz.
```

The endpoint Cauchy estimate `integral q^2<=ell^2*integral q_z^2/4` therefore gives

```
0 < (Eplus-Eminus)/(Eplus+Eminus) <= ell^2/4 =1/(4*10^12).
```

Each share differs from one half by at most `1/(8*10^12)`. The independent
calculation also gives the stronger elementary constant 1/6. These are diagnostics
of the normalized master, not physical Noether energies or radiation fluxes.
A spatially varying conjugation does not preserve these numerical weights.
Qualitative one-way tests on the compact physical mode must retain its own
conjugation. No result here excludes a past incoming state that develops mixed
present data through potential scattering or earlier shell interaction.

The [registered-packet audit](sources/REGISTERED_PACKET_AUDIT.md) records the
source and full first-data chain. A present-slice condition on the homogeneous
addition can select the declared reference, but changing the particular solution
changes that declaration. This is not an independently derived retarded rule.

## 3. Actual second mean, forcing and corner compatibility

In order `(K,L,N,J,w)` and `H=diag(1,2,1,3,1)`, the actual second mean satisfies

```
U_xx-U_rr=A U_x+B U_r+Dmat U+S0,
L(x,0)=0, PN U_r(x,0)=M PN U(x,0)+G0.
```

Here PN removes the orthogonal Dirichlet L slot. The actual adapter is

```
S0=P0(-2ell^2*rho^2*Q_E_reduced, +ell^2*rho^2*Q_phi),
G0=P0(2ell*rho*Q_R_tangent, ell*rho*Q_P, -ell*rho*G_N),
G_N=-2trace_gamma(Q_R)-2(Q_C)_y.
```

The Q terms are Hessians in the first solution, with no hidden Taylor factor 1/2.
These are actual complete quadratic sources in the retained gauge and frame.
Changing the second preparation does not change those first-quadratic sources.
The open initial shell collar gives the retained zero corner values needed by
the reduced problem. It does not make the bulk quadratic constraint zero:
a retained nonzero Hamiltonian-Hessian witness excludes setting every absolute
second Cauchy coefficient to zero.

With `P=PN(U_x+U_r)` and `u=PN U(x,0)`, the inhomogeneous boundary identities are

```
PN(U_x-U_r)=P-2(Mu+G0), u_x=P-Mu-G0.
```

Thus the incoming characteristic combination at the physical shell is constrained
by its junction law. It is not an additional freely prescribed source while
that law is held fixed. The [full proof](theory/ACTUAL_PACKET_AND_FORCED_MEAN.md)
retains actual canonical momentum normalization and all mixed metric velocities.

## 4. Forced tilted-energy estimate

Use the same constant-principal cone `D_T={0<=x<=T,0<=r<=T-x}`, `T=5/4`, and
future conditional matrix bounds as the
[H1-response result](../HDBLAST_H1_PREPARATION_RESPONSE_20261104/):
`||A||,||B||<=10^-4`, `||Dmat||<=10^-6`, `||M||<=1/2000` in H.
These remain explicit conditional analytic assertions, not new matrix proofs.

For `V=U_x,Q=U_r,k=alpha=1/2`, define

```
e=(|V|_H^2+|Q|_H^2+2k<V,Q>_H+|U|_H^2)/2,
E(x)=integral_0^(T-x)e dr+alpha|PN U(x,0)|_H^2/2.
```

The moving outer flux is `-(1-k)|V-Q|_H^2/2-|U|_H^2/2<=0`. The Dirichlet row
also contributes nonpositively. The homogeneous nonsymmetric Robin estimate
retains a residual `-|P|_H^2/8`. The source term satisfies the exact identity
`2e-|V+kQ|_H^2=(1-k^2)|Q|_H^2+|U|_H^2>=0`; Young's inequality therefore charges
bulk forcing by `e+|S0|_H^2/2`.

Absorbing the extra boundary G0 terms into the residual and boundary state gives

```
E' <= (301/100)E+(1/2)||S0(x,.)||_L2_r,H^2
                         +(5002001/4000000)|G0(x)|_H^2.
```

No symmetry or derivative of M, and no derivative of G0, is used. Rational Taylor
bounds give `exp((301/100)T)=exp(301/80)<44`. Hence

```
E(x)<=44[E(0)+B_src],
B_src=(1/2)integral_D_T |S0|_H^2 dx dr
       +(5002001/4000000)integral_0^T |G0|_H^2 dx.
```

This is a future reduced-system estimate on an existing compatible solution.
It constructs neither its initial completion nor any past retarded inverse.
The energy and B_src are analytic norms, not physical injected energy.

## 5. Continuous cell-clipped source budget

The [source ledger](sources/SOURCE_LEDGER.md) authenticates the retained continuous
source cells and current-state adapters. Intersecting their half-unit rectangles
with `D_T` and `x+r>=1` uses exact rational polygon areas. Five positive-area cells
partition active area `9/32`; boundaries have zero area. The shell source is active
on an interval of length `1/4`. Direct projection to all five H-weighted means gives

```
integral_D_T |S0|_H^2 <=8526418129275428520152526,
integral_0^T |G0|_H^2 <=3708320490000,
B_src<=4263209064642351515776089.
```

A separate full-angular Gram route supplies a conservative cross-check. These
are exact postprocessed ceilings from saved continuous enclosures, not source
samples or a fresh proof of those enclosures.

The boundary state gives `|w(x,0)|^2<=4E(x)<=176[E(0)+B_src]`. For the auxiliary
zero-Cauchy forced reduced component alone, the resulting trace ceiling is
`27392057158547`. That component generally fails the inhomogeneous Einstein
constraints by itself; the complete compatible initial correction remains needed.
Its ceiling is also looser than the retained Duhamel trace bound. No numerical
improvement of an observable or new absolute proper-volume/expansion sign follows.
Absolute proper observables retain their first-quadratic and clock terms, which
cannot be dropped or assigned to an unphysical reduced component.

## 6. A partial constraint-compatible past-null construction

On the actual background, extend the past null cut `t-z=tminus<-2ell` from its
shell corner `z=0` to `z=2ell`. Its packet-overlap segment is `ell<=z<=2ell`.
The lower limit called corner below is that shell intersection, with additional
corner data declared here; their full junction compatibility remains to be solved.
Use the common affine generator

```
k=c*rho^-2(partial_t+partial_z), c>0,
R0=rho*cosh(t), theta0=3c*rho^-2(tanh(t)-rho*a),
Phi0'=-c*p/rho.
```

Future k moves toward increasing z, away from the shell. Data on this outgoing
null segment can encode an incoming wave crossing it. This convention must not
be reversed. The background ODE `a_y=-rho^-2-p^2/3` verifies exactly
`theta0'+theta0^2/3=-(Phi0')^2`; the bulk potential drops out of null Ricci.
No future operator bound is being asserted on this past segment.

Primes denote affine derivatives. Phi_i and theta_i are actual parameter derivatives
of the pulled-back scalar and null expansion on the cut. For a twist-free
three-dimensional screen, fixed affine normalization across families and zero
background shear, the second trace constraint is

```
theta2'+(2theta0/3)theta2
 =-(2/3)theta1^2-2|shear1|_0^2-2(Phi1')^2-2Phi0'Phi2'.
```

The explicit new free-data choice `Phi2=0` on the cut, with corner `theta2=0`,
leaves the unique Raychaudhuri solution

```
theta2(lambda)=-2/R0(lambda)^2 * integral_corner^lambda
 R0(s)^2[theta1(s)^2/3+|shear1(s)|_0^2+(Phi1'(s))^2] ds <=0.
```

The logarithmic area-density second variation has a corresponding nonpositive
double integral if its corner value is declared zero. Actual area instead has
`A2/A0=average(log_area2+log_area1^2)`, so its sign is not fixed by that observation.
Neither null expansion is the present shell's proper-time expansion.

This solves one trace constraint under declared data. Remaining null constraints,
corner and Israel/scalar compatibility, a nonlinear characteristic development,
and transport to the present preparation fiber remain open. The
[full null note](theory/AFFINE_NULL_PREPARATION.md) states these hypotheses; this
is not a new existence theorem or a physical preparation law.

## 7. Stable causal identification on the retained preparation image

For a homogeneous preparation difference with the other canonical sectors fixed,
let h be its compact metric-mean increment, w its conditioned initial scalar mean,
and define `chi_y=rho^2 h`, `Z=chi/(ell*rho^3)`, `b=ell*rho*a`. The auxiliary Z is
distinct from the physical radial coordinate z. The retained constraint gives

```
Z_r=3b Z-h,
h_r=3ell^2 Z-bh+(2ell*rho*p/3)w, Z(1)=h(1)=0.
```

The free scalar Neumann wave has boundary value w(x), since its initial scalar
velocity is zero. In the max C1 norms using normalized x and r, the coupled
resolvent and zero wedge give

```
||m-w||_C1 <=kappa||w||_C1,
kappa=24020001201/63975980000000 <47/125000=0.000376.
```

Thus `(1-kappa)||w||_C1<=||m||_C1<=(1+kappa)||w||_C1`: a complete scalar history
causally identifies this restricted preparation, stably on its admissible image.
It does not identify exterior data, other canonical sectors or arbitrary nonlinear
families, and does not show that an arbitrary prescribed history is in the range.

The physical history is the second-family difference of the normalized scalar
mean `(integral phi_Sigma dV)/V`, at fixed anchored proper clock
`tau=rho_b*ell*x`. The background scalar is constant on the shell. Its volume
second variation cancels against normalization; shared first measure, pullback
and clock terms cancel in the difference, and the pure second normal gauge
vanishes at the shell. An unnormalized scalar integral would have an extra volume
term. See the [inverse proof](theory/SCALAR_HISTORY_INVERSE.md).

## 8. Vacuum exchange is not a radiation mechanism

The present shell has `S_ab=-sigma*gamma_ab`, `rho_shell=sigma`, `p_shell=-sigma`.
The action-derived scalar law and Codazzi identity give
`D_a S^a_b=2n_out(phi)D_b phi=-sigma_phi D_b phi`, hence the local proper-flow
balance `dot(rho_shell)+Theta(rho_shell+p_shell)=sigma_phi dot(phi)`.
This is exchange with vacuum tension, not a populated thermal particle sector.
The [conservation note](literature/MODEL_CONSERVATION_NOTE.md) derives the formal
matter extension with an explicit stress/matter-metric convention. That extension
is not implemented here. Strictly traceless classical radiation has no trace-only
conformal source; particle production would require specified interactions, state,
pump history, backreaction and thermalization.

The [primary-source comparison](REFERENCES.md) limits every transfer of a cited
result. No physical scale, observational prediction, uniform nonlinear amplitude
window or higher-dimensional origin of the Big Bang is established.
