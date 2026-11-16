# Uniform nonlinear relative-harmonic Cauchy data on the admitted initial interval

13 November 2026 UTC. Derivation with an independent mathematical review recorded in AUDIT.md. The result concerns
the complete initial field values and first time velocities. It does not
assert a uniform nonlinear evolution, a proper-observer interval or an
observable remainder. All constants below come from the evaluated whole-ball
composition and boundary-chart interfaces, with their stated analytic premises.

The main conclusion is the uniform completed-data bound

```
sup_(|e|<=10^-6402)
 [||partial_e^2 V_e(0)||H14(E)
       +||partial_e^2 partial_t V_e(0)||H13(E)] <= B_C,
B_C=10^10400.                                           (0)
```

The original admitted initial interval is retained. The Cauchy map is C2
there, with an open parameter neighborhood, in the precise target spaces
defined below. This is not C2 dependence of an evolved family in X14.

## 1. Fixed canonical curve and the actual completion problem

Retain the selected initial curve from the nonlinear constraint-ball
checkpoint on `|e|<=a_init=10^-6402`, with an open parameter neighborhood.
The calculation first concerns one bulk copy, on the physical collar
`Sigma=[0,6ell] x S3`, `ell=10^-6`. The correction annulus is
`ell/2<z<5ell/2`; put

```
r0=(z-ell/2)(5ell/2-z)/ell^2, kappa0=10^-8,
q=kappa0^2 r0^2, w=exp(-1/r0),
X=(dot H16_(q,w) scalar + dot H14_(q,w) one-form)^Gamma,
Gamma=(Z2)^4.
```

Gamma is the retained group of four coordinate reflections of the sphere.
The direct sum is Hilbert, and each dotted weighted space is the closure
of compact smooth tests in
`||T||Hm(q,w)^2=sum_(j=0)^m integral w^2 q^(2j)|D^jT|^2 dmu_gamma0`.
Write `u(e)=(N(e),Z(e))` for the selected potentials in X. The accepted
bounds on one copy are

```
||u(e)||_X<=K e^2, K=10^6400,
||u'(e)||_X<=.0202, ||u''(e)||_X<=8K,
u(0)=u'(0)=0.
```

The background spatial metric is
`gamma0=dy^2+rho^2 kappa_S3`, with `dy=-rho dz`, `129<=rho<=131`.
Use `D=nabla_gamma0` and the fixed canonical density
`rho0=sqrt(det gamma0)`. The canonical family is

```
gamma=gamma0+h,      h=e h_seed+A_H N(e),
phi=phi0+psi,        psi=e psi_seed,
pi_act=rho0 tau,     tau=e tau_seed+A_M Z(e),
varpi_act=rho0 omega, omega=e omega_seed.                 (1)
```

All four seed sectors are retained. The lapse and shift are the exact
registered choices, with no replacement by background lapse:

```
f=u_profile(y)Y(Omega), Y=n_1^2-n_2^2,
ell_lapse=e f, n=1+ell_lapse,
alpha=rho n, beta^i=0.                                 (2)
```

The lapse alpha is the temporal ADM lapse. It is distinct from the radial
component root `sqrt(g_yy)` used later in the boundary chart.

For d=0,1,2 use the accepted literal graded norms

```
E_d(T)=sum_(j=0)^(14-d)||w^-1 q^(j+d-2) D^jT||_L2(gamma0),
q=kappa0^2 r0^2, w=exp(-1/r0), kappa0=10^-8.
```

The accepted four-sector seed sum is below `10^200`. The same bound
holds for E0(f) directly: the registered seed satisfies
`h_seed=f(-4dy^2+2gamma_angular)`, hence
`f=(1/2)tr_gamma0 h_seed`. The trace has Hilbert norm two in four
dimensions and commutes with D, so `E0(f)<=E0(h_seed)<=10^200`.
This is a bound on the prescribed lapse profile, not on an unrelated
completed wave-chart representative.

The general correction maps give the literal-sum bounds
`E0(A_H N)<=10^30||N||` and `E1(A_M Z)<=10^27||Z||`.
Therefore the enlarged physical state
`z=(h,psi,tau,omega,ell_lapse)`, with norm
`E0(h)+E0(psi)+E1(tau)+E1(omega)+E0(ell_lapse)`, obeys

```
sup ||z|| <=10^-6201,
sup ||z'|| <=10^201,
sup ||z''|| <=10^6432.                                 (3)
```

The second line includes both the actual canonical seed and f. In the
third line only the two correction maps contribute: psi, omega and the
prescribed lapse perturbation are affine in e. This does not assert that
the physical normal scalar momentum or any completed wave velocity has
zero second derivative.

## 2. Canonical-to-physical conversion with the retained half-action

Set `J=det(gamma)/det(gamma0)`, `lambda=J^-1/2`, and
`T=tr_gamma tau`. The temporal normal is
`u_normal=alpha^-1(partial_t-beta)`. The retained action is
`R/2-|dphi|^2/2-U`, so its actual canonical definitions are

```
pi_act=sqrt(det gamma)(K^ij-gamma^ij tr_gamma K)/2,
varpi_act=sqrt(det gamma) Pi,
K=(1/2)Lie_(u_normal) gamma, Pi=u_normal(phi).
```

In four spatial dimensions the exact inverse map is consequently

```
K_ij=2lambda[gamma_ik gamma_jl tau^kl-T gamma_ij/3],
Pi=lambda omega,
tr_gamma K=-(2/3)lambda T.                              (4)
```

The older generic gauge companion displays the doubled-action canonical
densities. Its physical gauge identities apply unchanged, but its
canonical inverse formula cannot be substituted for (4).

## 3. All sixteen wave fields and their exact first time velocities

Use the fixed background identification `(t,y,Omega)` on the initial
section t=0. The metric values are

```
g_tt=-alpha^2, g_ti=0, g_ij=gamma_ij, phi=phi0+psi.       (5)
```

The relative-harmonic gauge is exactly

```
C^A=g^BC(Gamma[g]^A_BC-Gamma[g0]^A_BC),
R^A=g^BC Gamma[g0]^A_BC.
```

At t=0, the reference temporal coefficients
`GammaBar^0_ij=GammaBar^0_00=0`, while the cross inverse components
`g^0i=0` for every amplitude. Hence `R^0=0` identically on this entire
initial family, even though the other reference contractions depend on
the perturbed metric. The exact gauge completion gives

```
gamma_t=2alpha K, phi_t=alpha Pi,
alpha_t=alpha^2 tr_gamma K,
beta_t^i=alpha^2(Gamma_gamma^i-R^i)-alpha gamma^ij d_j alpha.
```

Substitution of (4) gives the complete diagonal/spatial and scalar
velocities:

```
g_ij,t=4rho n lambda[gamma_ik gamma_jl tau^kl-T gamma_ij/3],
g_tt,t=(4/3)rho^3 n^3 lambda T,
phi_t=rho n lambda omega.                              (6)
```

For the four mixed velocities define the actual covector

```
Ccov_i(gamma)=gamma^jk(D_j h_ki-(1/2)D_i h_jk).
```

The reference contribution is
`GammaBar^i_tt=rho D^i rho`, with the index raised by gamma0. Lowering
the shift completion with gamma gives the exact identity

```
g_ti,t=alpha^2 Ccov_i(gamma)+gamma_ij rho D^j rho-alpha d_i alpha
 =rho^2[n^2 Ccov_i(gamma)+a_warp h_iy
          -n D_i n+a_warp(1-n^2)theta_i],              (7)
```

where `theta=dy`, `a_warp=rho_y/rho`, and
`h_iy=h_ij(D^j y)`. Equivalently the background-subtracted last terms
are
`rho h_i^j rho_j-rho(2ell_lapse+ell_lapse^2)d_i rho
 -rho^2(1+ell_lapse)D_i ell_lapse`.

Formula (7) retains the metric dependence of the reference contraction.
At the registered first seed it gives
`g_ty,t,1=-6rho^2(u_profile,y+4a_warp u_profile)Y` and zero first
angular mixed velocity. The full nonlinear mixed covector need not
vanish. Specifying only initial mixed metric values would omit four
required Cauchy velocities.

Equations (5)--(7) solve the five relative-harmonic equations C=0
initially, using the original lapse and zero shift. No independent second
time accelerations are added. The subsidiary normal derivative uses the
reduced equations and the actual constraints, as in the retained source;
this note does not assume an evolved solution to define its Cauchy pair.

## 4. Uniform graded composition and physical Sobolev conversion

The companion `GRADED_CAUCHY_COMPOSITION.md` works on the explicit open
five-field ball

```
Q delta(z)<1/1024, Q=10^100,
delta(z)=E0(h)+E0(psi)+E1(tau)+E1(omega)+E0(ell_lapse).    (8)
```

The curve (3) lies strictly inside this ball. The accepted graded product
algebra concerns physical correction fields, not raw weighted normal
potentials. Constant identities and fixed background fields act in a
separately bounded multiplier algebra.

In particular, the positive powers of rho and the covector `rho D rho`
are quantitatively controlled, rather than folded into an unknown constant.
Their scaled covariant jets through order fourteen obey

```
sup |q^j D^j(rho^k)| <=131^3*5^j*j!, |k|<=3,
sup |q^j D^j(rho D rho)| <=131^2*6^j*j!,
max_(j<=14) of these bounds <10^28.                     (9)
```

The recurrences include all derivative slots; their respective next-stage
factors are bounded by `5(j+1)` and `6(j+1)`. The finite covariant Leibniz
sum then gives a common multiplier norm `10^33` on every required E_d.
The fixed density rho0 remains parallel and fixed throughout.

For `A=(I+H)^-1`, `H=gamma0^-1 h`, the full matrix Neumann/logarithm
algebra on (8) gives

```
E0(A-I)<=2E0(h),
E0(DA[h_i])<=4E0(h_i),
E0(D^2A[h_1,h_2])<=16Q E0(h_1)E0(h_2).
```

For the positive fractional determinant `lambda=J^-1/2=1+s`, use
`L=log det(I+H)` and `lambda=exp(-L/2)`. The resulting bounds are

```
E0(s)<=4E0(h),
E0(Dlambda[h_i])<=4E0(h_i),
E0(D^2lambda[h_1,h_2])<=16Q E0(h_1)E0(h_2).             (10)
```

These are uniform estimates for the full inverse and determinant on the
positive ball. They do not restrict the metric to a scalar-conformal form.
Their convergent series give C2 maps on the whole open ball, including
spatial derivatives falling on arbitrary state factors.

Let W(z) be the background-subtracted raw Cauchy pair given by (5)--(7).
Its literal graded block sum uses E0 for positions, E1 for velocities,
and a factor two for the mixed metric covector to cover its two symmetric
slots. The finite formula contains at most 32 monomials, at most six
nonbackground factors per monomial, and coefficient/contraction norm at
most sixteen per monomial. Every term contains at least one physical
perturbation factor. The remaining factors are the unital maps
`I+H,A,1+ell_lapse,lambda`, whose module norms are below two.

In the rescaled ideal norms `e_d=Q E_d`, each factor's first direction
costs at most `4 sigma_i` and its second direction at most
`16 sigma_1 sigma_2`, with `sigma_i=Q delta(z_i)`.
The product rule gives respective six-factor ceilings

```
state: 32 Q delta,
first direction: 768 sigma_i,
second direction: 10752 sigma_1 sigma_2.
```

Combining these counts with the background multiplier `10^33`, and
undoing the Q rescaling, yields the convenient uniform roundings

```
||W(z)||graded <=10^40 delta(z),
||DW(z)[z_i]||graded <=10^40 delta(z_i),
||D^2W(z)[z_1,z_2]||graded <=10^140 delta(z_1)delta(z_2). (11)
```

The determinant Hessian, inverse-metric Hessian, nonzero mixed velocity
and both parameter allocations in every product are included in (11).
No pointwise high-derivative assumption is made on the actual first
parameter derivative.

Define the raw physical norm using the positive bundle metric
`dt^2+gamma0`, the connection extended by `D dt=0`, and the full measure
`dvol_gamma0`. Set

```
X=(g-g0,phi-phi0)|t=0, Xdot=partial_t(g-g0,phi-phi0)|t=0,
R_j=||partial_e^j X||H14(gamma0)
                  +||partial_e^j Xdot||H13(gamma0), j=0,1,2. (12)
```

The metric/scalar Hilbert sum is dominated by the literal block sum in
(11). For either grade d=0 or 1 the exact scalar factorization

```
D^j T=[w q^(2-j-d)] [w^-1 q^(j+d-2)D^j T]
```

and `sup w q^-12<10^226` give
`||T||H^(14-d)(gamma0)<=10^226 E_d(T)`. No derivatives of a hidden
norm equivalence are used here. Consequently the raw map and its first
and second Fréchet derivatives have respective bounds
`10^266 delta`, `10^266 D_i`, and `10^366 D_1D_2`.

Apply the actual curve chain rule, including its acceleration z'':

```
R_0<=10^266*10^-6201=10^-5935,
R_1<=10^266*10^201=10^467,
R_2<=10^266*10^6432+10^366*(10^201)^2
       =10^6698+10^768<=10^6700.                       (13)
```

The second term is the raw nonlinear composition Hessian. It has not
been removed because some canonical sectors are affine.

## 5. Boundary chart and final completed norm

The target Cauchy-data norm is

```
||V_e(0)||_H14(E,dz^2+kappa)+||partial_t V_e(0)||_H13(E,dz^2+kappa),
E=R_v+Sym^2(R dt+T*S3)_H+R_w+(R dt+T*S3)_D,
```

with full Frobenius norm and fixed product connection and measure. The
companion `UNIFORM_BOUNDARY_CAUCHY_TRANSPORT.md` proves its uniform chart
estimate on the sufficient raw ball `R_0<=10^-3228`. Equation (13)
places the actual whole curve strictly inside that ball directly. This
admission is not inferred by applying Taylor's theorem with the Hessian
bound being proved.

The fixed raw-to-product fiber identification I sends
`(dt,dy,rho theta)` to `(dt,dz,theta)` on every tensor and derivative
slot. It is not the coordinate pullback `dy=-rho dz`. The retained
linear operator bound is `10^3000` through order fourteen. Its actual
normalization L(z) is

```
x=delta g_yy,
H_tt=delta g_tt/rho_b^2,
H_tA=rho delta g_(t,phys A)/rho_b^2,
H_AB=rho^2 delta g_(phys A,phys B)/rho_b^2,
w=delta phi,
D_t=delta g_yt/rho_b,
D_A=rho delta g_(y,phys A)/rho_b.
```

With its evaluated background derivatives, the combined linear map
obeys `||L I||Hm->Hm<=A_lin=10^3216`, m<=14. Put

```
P=L I X, Qvel=L I Xdot,
r(P)=(sqrt(1+x)-1,H,w,D), T(z,P)=F_-z(r(P)).             (14)
```

Here F_-z is the actual boundary-frozen state flow and rho_b is the
fixed shell radius. The root is the proper-y component root, not the
temporal lapse and not a Schur-complement radial lapse.

On the sufficient raw ball, `||P||H14<=10^-12` and the evaluated H3
embedding gives `sup|P|<=10^-5`. The root and both signs of the actual
flow remain below `1.01*10^-5<10^-3` in the star state norm. Thus the
entire curve uses the proved inner state patch, with strict room for
function-space differentiation. The flow's total mixed z/state derivatives
through eighteen, together with the root derivatives, give a point
coefficient envelope `10^420` for T.

The companion's all-mode Sobolev product proof expands at most 18!
spatial chain terms. At most one field factor needs its highest derivatives
in L2; all other factors have three derivatives available for the tensor
H3 embedding. Its uniform state derivative bounds in H14 and H13 are

```
C_1=10^437, C_2=10^444, C_3=10^451.                     (15)
```

Velocity arguments need only H13. Base-state factors are controlled by
`10^7||P||H14<=10^-5`, which is how the higher spatial chain terms are
uniformly summed. These are whole-state estimates, not the older A/B
maps evaluated only at zero amplitude.

At t=0 the explicit first time derivative of the de Sitter coframe
normalization vanishes for every state: its factors are powers of
`sech(t)`, and the subtracted background is even. The root and autonomous
flow have no remaining explicit time dependence. Therefore the exact
completed pair is

```
V=T(P), V_t=T_1(P)[Qvel].                              (16)
```

For primes denoting the original amplitude derivatives,

```
V''=T_1[P'']+T_2[P',P'],
(V_t)''=T_1[Qvel'']+T_2[P'',Qvel]
           +2T_2[P',Qvel']+T_3[P',P',Qvel].             (17)
```

All coefficients are evaluated at the same nonzero P. The factor two
and the third state derivative term are compulsory. The largest flow
order needed is 14+2 or 13+3, within its eighteen-jet certificate.
A third state derivative of the chart is not a third amplitude derivative
of the family. Only its existing C2 parameter regularity is used.

Combining (15)--(17), `||Qvel||H13<=1`, and A_lin, the companion proves
on the whole raw ball

```
||Cauchy(e)|| <=10^3653 R_0,
||Cauchy'(e)|| <=10^3661 R_1,
||Cauchy''(e)|| <=10^3661 R_2+10^6884 R_1^2.             (18)
```

Insert the actual raw bounds (13):

```
sup ||Cauchy(e)|| <=10^-2282 <=10^-2200,
sup ||Cauchy'(e)|| <=10^4128 <=10^4200,
sup ||Cauchy''(e)|| <=10^10361+10^7818
                         <10^10362<10^10400.           (19)
```

This proves (0) on the full original initial interval. The value and first
derivative estimates in (19) are separate direct bounds; they are not
assumptions used to manufacture a whole-family radius from zero jets.

Equations (3),(13) were computed on one bulk copy. For the retained union
of two identical reflected copies, the disjoint-union Hilbert norm gains
exactly sqrt(2), conservatively at most two. The duplicated raw data still
satisfy the chart admission criterion. The completed value and first
derivative bounds can gain at most two, and the sufficient Hessian join
can be bounded by
`2*10^10361+4*10^7818<10^10400`. Thus all three rounded final ceilings
`10^-2200`, `10^4200`, and `B_C=10^10400` remain valid for that specified
identical union. No bound for an arbitrary number of copies is asserted.

## 6. C2 regularity, support and scope

The admitted canonical branch is C2 on an open neighborhood of the closed
interval. The five-field state curve remains strictly inside (8).
The accepted pointwise estimate gives
`|ell_lapse|<=10^12 E0(ell_lapse)<=10^-6189<1/2`, so the
prescribed lapse is positive. The raw completed curve remains strictly
inside the chart ball. The graded
composition, fixed linear transport and nonlinear chart are C2 on open
neighborhoods of these sets. Their composition is therefore C2 into
the retained H14-position/H13-velocity Cauchy space on an open parameter
neighborhood of the stated closed interval. All endpoint derivatives in
(0),(19) have this ordinary two-sided interpretation.

The canonical seed and corrections vanish outside their original fixed
annulus. The prescribed lapse perturbation has the smaller seed support.
The algebraic completion and its spatial derivatives are local, and every
background-subtracted expression (5)--(7) vanishes when those inputs do.
Thus positions and velocities, and their two parameter derivatives,
vanish on the same fixed open shell and outer collars. The physical
zero-extension results prevent annulus-boundary distributional layers.
Finally `F_-z(0)=0` and (16) preserve these collars after straightening.

The construction retains the original physical initial section, amplitude,
normal orientation, full nonzero-velocity seed, half-action normalization,
and relative-harmonic representative. It does not replace the completed
mixed velocities by zero or choose a new lapse. The canonical identities
`phi''=0` and `varpi_act''=0` do not become such identities for Pi,
the wave velocities, the radial root or the boundary-chart fields.

No independent higher time Cauchy derivatives have been estimated here.
A quantitative nonlinear initial time-jet recursion and the variable-state
PDE estimates on a prescribed slab remain separate steps. In particular
(0) is not an evolved X14 bound, an observer-domain estimate or a
finite-amplitude observable remainder.

## 7. Source ledger

The two independently accepted complete interfaces are the reviewed originals
listed below; public editorial changes and hashes are recorded in AUDIT.md:

| Companion | Accepted SHA256 |
|---|---|
| `GRADED_CAUCHY_COMPOSITION.md` | `b7d1fdb558f985dcd5054db7bb9893f41cc700e33b2b537bdc0728fcf316d812` |
| `UNIFORM_BOUNDARY_CAUCHY_TRANSPORT.md` | `a303689363f59115d071c5af7badf24a11d52103cbe3e02163cc80a89a8b2c1e` |

Retained paths below are relative to
the repository `research/` directory:

| Source | SHA256 |
|---|---|
| `HDBLAST_NONLINEAR_CONSTRAINT_BALL_20261113/NONLINEAR_CONSTRAINT_BALL.md` | `04f35658bd66d2c3cfb9d1b11cbf942d349595e0231dd11b3a8d847241ab1f54` |
| `HDBLAST_NONLINEAR_CONSTRAINT_BALL_20261113/PHYSICAL_CORRECTION_MAPS.md` | `6c8b365adb37b7ca15cc59c071a1cad44606d7c0200b66225de1f06fe49cc2d0` |
| `HDBLAST_NONLINEAR_CONSTRAINT_BALL_20261113/SEED_ENVELOPE.md` | `17bbd5f7da7266cd2aa91fd76e74e13e5f0df6ada9cfeecb6462a9a80603456e` |
| `HDBLAST_INITIAL_TRANSPORT_20261106/INITIAL_TRANSPORT.md` | `8a0cfa469037427f8636344c7635fcdd8249039564d16242d3f050ce912ed0bc` |
| `HDBLAST_INITIAL_TRANSPORT_20261106/inputs/GAUGE_CAUCHY_COMPLETION_06G.md` | `392c4ee59957797b9812062fc4691cf006cc95c806ba7a5004b8fd43509dfd07` |
| `HDBLAST_INITIAL_TRANSPORT_20261106/inputs/REGISTERED_EXACT_FAMILY_BINDING.md` | `69bfe846cc4fda8ec6a981969fbdd12fff7862b5ca7cf31a19d696b0a839165a` |
| `HDBLAST_QUASILINEAR_SHELL_EVOLUTION_20261015/gauge-parameter/INITIAL_GAUGE_SUBSIDIARY_AND_PARAMETER.md` | `74bce09080d4015b2eccf8c183b33f47b80516216a9547837e3e0e34251d7b03` |

The retained analytic background and selected weighted-family premises
remain visible in these sources. The new constants are evaluated analytic
ceilings in the same model normalization; no physical-unit interpretation
or external confirmation is inferred.
