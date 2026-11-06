# Known quadratic proper-area sources and the remaining curvature projection

6 November 2026 UTC. Conditional mathematics for the registered affine
amplitude e, at e=0. A subscript m means the m-th amplitude derivative,
not the coefficient of e^m; the latter is divided by m!. This note evaluates
known quadratic sources. It does not evaluate the full second field,
uniform nonlinear remainder or admissible amplitude radius.

## 1. Fixed physical observation and retained inputs

Keep the registered Einstein–scalar branch, ell=10^-6 packet, nonzero
initial velocity exp(-z)Q(z/ell), original physical initial anchor,
transported round labels n on S3, and proper q=-tau/rho_b. The exact family
and its higher completion are fixed. Retain the smooth common linear and
nonlinear developments, reconstruction, gauge, uniqueness and observer
domain premises of the linked sources. In particular this is an intrinsic
shell geodesic congruence, initially normal to the physical anchor. Before
its cut/flow fails it is hypersurface orthogonal, with pulled-back metric

```
X_e^*h_e = -rho_b^2 dq^2 + k_e(q,n),
k0=rho_b^2 cosh(q)^2 kappa,   G_e=k0^-1 k_e,
j_e=sqrt(det G_e),   j0=1.
```

The background scale rho_b is fixed, independent of e, and 129<=rho_b<=131.
All quantities below use the model units; no SI or observational calibration
is made. The fields, observers, initial anchor and labels must be transported
together under a gauge change. Coordinate clock pieces alone are not used.

The pinned [first-cut theorem](inputs/FIRST_CUT_BOUND.md) supplies the actual
past first field u(q), |u|<U=23/(2*10^9), and

```
I(q)=integral_0^q u(s) ds,
J(q)=integral_0^q sech(s)^2 I(s) ds,
G1=2[u-tanh(q)I]Y Id-2J Hess(Y),  Y=n1^2-n2^2.
```

The first tangent vanishes through q=ell. Set T=25ell/16, D=T-ell=9ell/16;
then |I|<=U(q-ell)_+ and |J|<=U(q-ell)_+^2/2. For higher full-family
quantities only the retained nonlinear quiet interval q<=qD=ell/4 is
available. Do not extend it to ell. Its proof is
[the exact quiet-window note](https://github.com/maldonado-research/HDblast/blob/a715e8ee2963467b8988b7984a1cd0d402468f3d/research/HDBLAST_QUIET_SHELL_INTERVAL_20261105/EXACT_QUIET_WINDOW.md).

## 2. Actual first time derivative: complete endpoint energy

The pinned [higher-jet derivation](inputs/HIGHER_JETS.md) and
[report](inputs/HIGHER_JETS_NORMAL.json) bound the actual future shell
trace y_+(x), x=t/ell, and its x derivative on 0<=x<=3. The report's
exact.boundary.yprime is the positive rational B_yprime used below.
No new master/PDE enclosure is run here.

The [same-state past/future proof](inputs/PAST_FUTURE_PROOF.md), sections
2–3, retains H=L2(0,infinity) direct-sum R with norm squared
integral v^2+Lambda*b^2, Lambda>1000, and the linked positive operator
A0>=5500ell^2. Its complete dynamic endpoint is part of this Hilbert
evolution, not an extra static Robin condition. For
Delta=Z_minus-Z_plus, the initial data are (0,-2g), ||g||<=ell. The temporal
potential omega=-(35/4)ell^2 tanh(ell*x)^2 is decreasing on 0<=x<=3,
and A0+omega is positive there. The complete energy therefore obeys

```
E_Delta(x)<=2ell^2,  ||Delta_x||H^2<=4ell^2,
|y_-'(x)-y_+'(x)|^2 <=4ell^2/Lambda <ell^2/250.
```

Since (253/4000)^2>1/250, |delta y_x|<(253/4000)ell. This uses the
same reflected nonzero velocity, full endpoint mass, common smooth extension
and local support as the prior proof. It is not an independent master solve.
Writing u=rho_b^-2 cosh(q)^(-3/2)y_-(q/ell) gives

```
|u_q| <= B_uq :=129^-2 [B_yprime/ell+253/4000+(3/2)T*(190ell)]
             <293/10000=0.0293,
F(q):=3u_q-3tanh(q)u+5sech(q)^2 I,
|F|<=Fmax:=3B_uq+3TU+5UD <11/125=0.088.
```

These derivatives are with respect to dimensionless q, not physical proper
time. Smooth prefront flatness gives F=0 for q<=ell. The original report's
actual input joins, scaled domain and pins are retained.

## 3. Exact angular and determinant decomposition

Let A=diag(1,-1,0,0), s=n1^2+n2^2, and A_T its compression to T_n S3.
For unit tangent v, the sphere Hessian is
HessY=2A_T-2Y Id. Define

```
H=HessY+(8/3)Y Id=2A_T+(2/3)Y Id,   trH=0,
K(n)=tr(H^2)=8(1-s)+(8/3)Y^2,       0<=K<=8,
B(q)=u-tanh(q)I+(8/3)J,
G1=2BY Id-2JH,                    j1=3BY.
```

Here tr(A_T)=-Y and tr(A_T^2)=2-2s+Y^2 prove the formula for K.
As |Y|<=s<=1, K<=8(1-s)+(8/3)s^2<=8. The operator bound
||H||op<=2 is sharp. For unit v perpendicular to n,
v1^2+n1^2<=1 and v2^2+n2^2<=1, which imply
|v^T A v+Y/3|<=1. Equality is attained at n=e3, v=e1 or e2.

For the relative area density, differentiating the determinant gives

```
j2=(1/2)trG2+F_det,
F_det=(1/4)(trG1)^2-(1/2)tr(G1^2)=3B^2Y^2-2J^2K.
```

This is an exact separation of the known quadratic determinant term from
the unknown true second-cut metric. It does not set G2 to zero. The
round-label mean <f>=(2pi^2)^-1 integral f dOmega satisfies

```
<Y^2>=1/6, <Y^4>=3/40, <Y^2s>=1/8,
<K>=40/9, <Y^2K>=8/15, <Y F_det>=0,
<F_det>=B^2/2-(80/9)J^2.
```

These follow from the exact S3 moments
<product n_i^(2a_i)>=product(2a_i-1)!!/[4*6*...*(4+2sum(a_i)-2)].
Odd coordinate powers vanish. The quarter-turn exchanging n1,n2 reverses
Y and H, while K and F_det are even; its zero quadrupole projection is
consistent with the retained parity, not a new independent prediction.

For the registered normalized denominator, use the fixed weight

```
<f>_w=<Y^2 f>/<Y^2>,
D_e/D0=<j_e>_w,  D1=0,
(D2/D0)_det=<F_det>_w=(27/20)B^2-(32/5)J^2.
```

Only this known determinant contribution has been evaluated; the other
piece is (1/2)<trG2>_w. Since
|B|<=Bmax:=U[1+TD+(4/3)D^2], exact arithmetic gives

```
-2.95e-41 < <F_det> < 6.613e-17,
-2.12e-41 < (D2/D0)_det < 1.786e-16.
```

The inequalities enclose the entire interval 0<=q<=T, in the derivative
convention stated at the top. At Y=0,n=e3, F_det=-16J^2 although j1=0.
Unknown G2 can still change the total sign or magnitude of j2.

## 4. First shear and the exact curvature transfer identity

Let C_e=(1/2)k_e^-1 partial_q k_e, theta_e=tr C_e and
Sigma_e=C_e-(theta_e/3)Id. The normalized background is
C0=tanh(q)Id, theta0=3tanh(q), Sigma0=0. Since k0 is a scalar warp,
C1=(1/2)partial_q G1. Its trace and trace-free parts give

```
theta1=F(q)Y,    Sigma1=-sech(q)^2 I(q) H,
tr(Sigma1^2)=sech(q)^4 I(q)^2 K(n).
```

The future unit expansion and shear are -theta/rho_b and -Sigma/rho_b.
For example ||sigma1||op<=2UD/129<1.01e-16 in inverse model-proper-time
units. The round-label mean contraction is (40/9)sech^4(q)I^2.

With the curvature convention Ric_h0=3h0/rho_b^2, define the *full*
pulled-back intrinsic shell curvature scalar
R_e=Ric_(h_e)(partial_q X_e,partial_q X_e). It includes field,
embedding, clock and observer variations, and R0=-3. The intrinsic
geodesic Raychaudhuri identity is

```
theta_q=-theta^2/3-tr(Sigma^2)-R_e.
```

It reproduces theta0_q=3sech^2(q). Differentiating twice in e gives

```
theta2_q+2tanh(q)theta2=-(2/3)theta1^2-2tr(Sigma1^2)-R2.
```

Put L_e=log j_e; then L2_q=theta2 and j2=L2+j1^2. All initial
amplitude variations vanish. Solving the linear equation and integrating
once more supplies the exact decomposition

```
L2 = L2_exp + L2_shear + L2_R,
L2_exp(q,n)=-(2/3)Y^2 integral_0^q cosh(s)^2 F(s)^2 [tanh(q)-tanh(s)] ds,
L2_shear(q,n)=-2K(n) integral_0^q sech(s)^2 I(s)^2 [tanh(q)-tanh(s)] ds,
L2_R(q,n)=-integral_0^q cosh(s)^2 R2(s,n)[tanh(q)-tanh(s)] ds.
```

These are zero-data particular components. They do not solve for R2 or
the full area. The negative kernels of the two known components do not
imply net contraction, net focusing, a force or a singularity.

## 5. Quantified known terms for the observable denominator

Let Cosh=(1-T^2/2)^-2<1+3e-12. The bound follows from
cosh(T)<=sum_m(T^2/2)^m; (2m)!>=2^m. Also
0<=tanh(q)-tanh(s)<=q-s. Using the exact weighted moments gives

```
<L2_exp>_w=-(3/10)integral cosh(s)^2 F(s)^2[tanh(q)-tanh(s)] ds,
<L2_shear>_w=-(32/5)integral sech(s)^2 I(s)^2[tanh(q)-tanh(s)] ds,
<j1^2>_w=(81/20)B^2.
```

For q<=T the first two integrals start at ell, and

```
0<=-<L2_exp>_w <= (3/20)Cosh Fmax^2 D^2 <3.65e-16,
0<=-<L2_shear>_w <= (8/15)U^2 D^4 <7.062e-42,
0<=<j1^2>_w <=(81/20)Bmax^2 <5.357e-16.
```

The D^4/12 factor in the shear integral comes from integrating
(s-ell)^2(T-s). Define the known sum
K_D=<L2_exp+L2_shear+j1^2>_w. Then

```
D2/D0 = K_D + <L2_R>_w,
-3.65e-16 < K_D <5.357e-16.
```

This is an alternative exact closure to the determinant decomposition in
section 3; do not add both decompositions and count first-field terms
twice. Mean pointwise log density <L2> is different from the second
derivative of log total area; the latter includes the variance of j1.
Likewise the fixed weighted denominator uses the weights stated above,
not an e-dependent averaging measure.

The remaining second-denominator input has therefore narrowed to a
specific intrinsic-curvature projection. If a future computation supplies
Cr>=sup_(0<=q<=T)|<R2(q)>_w| on the same physical cylinder, the retained
full nonlinear quiet interval gives

```
|<L2_R>_w| <= C_R Cr,
C_R=(1/2)Cosh*(T-ell/4)^2 <8.614e-13.
```

Cr is not evaluated here. In particular the ell/4 lower integration limit
cannot be replaced by ell. This closes neither uniform parameter
derivatives away from e=0 nor the scalar, odd remainder, common-domain
and higher embedding estimates needed for the full observable.

## 6. Scope and next step

Bulk energy conditions alone do not determine this intrinsic shell R2:
the Gauss equation contains projected bulk curvature and shell geometry.
Even knowledge of G1 and its shear leaves a second intrinsic field input.
As a purely geometric illustration, G_e=exp(eG1+e^2 c(q)Id), with c smooth
and zero through ell/4, has the same G1 but L2=3c with either sign. This
is not asserted to satisfy the registered Einstein–scalar equations or
constraints; no freedom to adjust their fixed exact curve is inferred.

The practical next target is an actual bound on <R2>_w, or an equivalent
second-field evolution estimate, followed by compatible uniform mixed
field/observer bounds on one admitted amplitude interval. Determinant and
Raychaudhuri identities are standard. The new evaluated application and
explicit unresolved projection are conditional research progress, without
an originality theorem, external peer review, journal acceptance, calibrated
prediction or confirmation of a higher-dimensional Big Bang origin.
