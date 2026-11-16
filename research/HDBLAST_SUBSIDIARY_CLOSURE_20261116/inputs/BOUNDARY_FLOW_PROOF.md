# Quantified retained boundary straightening: four Dirichlet and twelve Neumann fields

Conditional mathematical checkpoint, 6 November 2026. This quantifies the actual retained
background-subtracted boundary flow. It proves a local chart, its inverse,
and exact primal equivalence. It supplies no evolved norm, nonlinear
amplitude interval, observer flow, or clock-domain conclusion. Background
existence/enclosures retain their source hypotheses. Inputs are read-only.

## 1. Source convention, exact state and frozen generator

`SOURCE_PINS.json` pins the retained public inputs by immutable source path,
byte count, SHA256 and Git blob identity. The exact source locations establishing the convention are:

* `ANGULAR_UNIFORM_LINEAR_CONTACT.md` §1: `gbar=dy²+rho(y)²q_dS`,
  `q_dS=−dt²+cosh(t)²kappa_S3`, `z=log(b/r)`, `dy=−rho dz`.
  Section 2 gives the exact proper-y boundary completion.
* `NONLINEAR_BOUNDARY_STRAIGHTENING.md` §1 equations (1)–(4):
  proper-y radial component square root `N=sqrt(g_yy)`, full-collar
  background subtraction, and the autonomous boundary-frozen generator
  `B(u)=−rho_b[f(q0,b+u)−f(q0,b)]`.
* `REGISTERED_EXACT_FAMILY_BINDING.md` §1: the same fixed background chart,
  exact current c, delta, near-phi=+1 branch and wider-collar bounds.
* `MODEL_AND_ARCHIVE_CONTEXT.md` §§1–2: the tube `|phi−1|<=.01` and the
  shell relations `a_b=sigma_b/6`, `p_b=−sigma'_b/2`.
* `INITIAL_GAUGE_SUBSIDIARY_AND_PARAMETER.md` (2), §§2–3: relative harmonic
  gauge and its temporal ADM completion. Its temporal lapse alpha is a
  different variable from the radial component square root here.
* `EXACT_FINITE_REGULARITY_SHELL_EVOLUTION.md` §3: sixteen-field scalar
  principal reduction and fixed Dirichlet/Neumann split.
* The authenticated `ACTUAL_POTENTIAL_JET_BOUNDS.json` records the actual
  wider rectangle and exact enclosures
  `rho in [129,131]`, `eta in [−9/100000,−1/12500]`,
  `a in [111/1000,14/125]`, `p in [−7/50000,−3/25000]`.
  It states its conditional retained-background scope explicitly.

The coframe `(dt,cosh(t)theta^i)` is q_dS-orthonormal and independent of z;
write `J=diag(−1,1,1,1)`. Define the sixteen **perturbation** variables

```
u=(v,Htilde,w,D),  n=1+v=sqrt(g_yy),
Htilde=(h−hbar(z))/rho_b²,  w=phi−phibar(z),
D_a=g_ya/rho_b.                                           (1)
```

There are ten symmetric Htilde entries, one v, one w and four D entries.
All rho_b factors are constant shell normalizations. At the boundary put
`H=J+Htilde`, `eta=eta_b+w`, `sigma_b=sigma(eta_b)` and
`sigma'_b=sigma'(eta_b)`. Then the literal retained generator, expressed
in these constant-scaled perturbations, is

```
B_v = −rho_b{(2/3)[sigma(eta)n²−sigma_b]
                   −a_b[n³ tr(H^{-1}J)−4]},
B_H = −rho_b[ sigma(eta)n H−sigma_b J ]/3,
B_w =  rho_b[ n sigma'(eta)−sigma'_b ]/2,
B_D = 0.                                                 (2)
```

Here `eta=phi−1`,

```
sigma(eta)=2/3+delta(1+c)+delta*c*eta+2eta²+(2/3)eta³,
delta=1/1000,  c=96714321425692/161839258930731.
```

Using the exact background identity `sigma_b=6a_b`, the first row
is equivalently `B_v=−rho_b[(2/3)sigma*n²−a_b*n³tr(H^{-1}J)]`.
The subtracted form makes zero preservation manifest. The boundary rows
are exactly `u_D=0`, `u_N,z=B_N(u_N)`. To verify all
signs directly, the proper-y outward normal on `g_ya=0` is
`n_out=n^{-1}partial_y`. The physical laws give
`h_y=sigma*n*h/3` and `phi_y=−sigma'*n/2`. The relative harmonic vector
satisfies on the boundary

```
g^{BC}Gamma[g]^y_BC = n_y/n³−2sigma/(3n),
g^{BC}Gammabar^y_BC = −a_b tr_h(hbar_b).
```

Tangential derivatives of the vanishing Dirichlet trace are zero, so they
contribute nothing here. Consequently `C^y=0` gives
`n_y=2sigma*n²/3−a_b*n³ tr_h(hbar_b)`. Subtract the actual background
normal derivative, apply `partial_z=−rho_b partial_y` at the shell and
the constant scales in (1), obtaining (2). In particular `B(0)=0` exactly:
`(2/3)sigma_b−4a_b=0`, and both displayed subtractions vanish at zero.

Equation (2) is autonomous in the auxiliary flow parameter and independent
of D. Its coefficients are the shell values, as in the retained source.
It is not the physical bulk radial ODE and its parameter is not time.
The nearby raw conformal component `sqrt(g_zz)` obeys a different raw
boundary formula; see [coordinate dictionary](COORDINATE_DICTIONARY.md). No raw conformal flow
is substituted for (2) in this primary proof.

## 2. Outer domain, invertibility and exact derivative bounds

For complex perturbations set

```
||u||_* = max(|v|, ||Htilde||_op, |w|, ||D||_2),
```

using the positive Hermitian Euclidean operator norm on matrices. Let
`O={||u||_*<.01}`. On its closure, `|n|<=1.01`,
`||H^{-1}||<=100/99`, and `|tr(H^{-1}J)|<=400/99` by Neumann series.
Also `|eta|<=.0101`, since `|eta_b|<.0001`. This last scalar-center
allowance is necessary: the outer perturbation disk is not centered at
eta=0. Exact rational bounds give throughout this complex domain

```
|sigma|<s=.67, |sigma'|<p=.042, |sigma''|<q=4.1.           (3)
```

The weaker eta_b bound follows even from the retained tube plus
`−2p_b in [.00024,.00028]`: the polynomial sigma' is increasing on the
real tube, and its values at ±.0001 bracket that interval. The checker
verifies this without selecting a wrong scalar root or dividing by sigma'.

Set `r=.01`, `nmax=1+r`, `I=(1−r)^{-1}`, `T=4I`, `rho_max=131`,
`a_max=14/125`, and `p_b,max=.00028`. Component bounds for B are

```
Rv=(2/3)*131*s*nmax²+131*a_max*nmax³*T,
RH=131*[s*nmax*(1+r)+s]/3,
Rw=131*[nmax*p+p_b,max]/2.
```

For a unit *-norm variation use
`d(H^{-1})=−H^{-1}(dH)H^{-1}` and sum the contributions of the three
nontrivial input blocks. The corresponding operator-row bounds are

```
Lv=(4/3)*131*s*nmax +3*131*a_max*nmax²*T
     +(2/3)*131*p*nmax² +4*131*a_max*nmax³*I²,
LH=131*[s*(1+r)+nmax*p*(1+r)+nmax*s]/3,
Lw=131*[p+nmax*q]/2.
```

D has zero derivative row. Exact arithmetic gives

```
||B||_* <120.767<140,       ||DB||_{*->*}<365.050<400.     (4)
```

Thus these are complex operator bounds in the declared block norm, not
coordinate bounds silently interpreted as operator norms. The positivity
of sigma' is unnecessary and its zeros are nonsingular points of B.

For real states, H is Lorentzian and n>.99. The original metric-state
conversion is locally invertible: `g_yy=n²`, `g_ya=rho_b D_a`,
`h=hbar+rho_b²Htilde`. Its derivative in the n variable is 2n, bounded
away from zero. On the actual collar, `rho_b/rho<=131/129`; the much
smaller inner images below lie within the explicit metric/energy patch.
The entire outer flow domain is only a chart/analyticity domain; it is
not silently identified with the tighter energy patch.

## 3. Complex existence and real forward/inverse bounds

Let `ell=10^{-6}` and `F_z` be the flow of (2). Fix any real state center
`u_c` with `||u_c||_*<=.003`. Allow independent complex increments at
most `.001` in each of the sixteen scalar coordinates (including the
ten independent symmetric matrix entries). Such a matrix increment has
operator norm at most `.004`, and the D-vector increment has norm at
most `.002`. Every such initial state has * norm at most `.007`.

For complex `|z|<=8ell`, Picard integration along the segment from zero
to z maps the ball of analytic trajectories of displacement at most
`140*8ell=.00112` into O, since

```
.003+.004+.00112=.00812<.01.
```

Its Lipschitz constant is at most `400*8ell=.0032<1`. The contraction
therefore constructs the holomorphic solution on the whole complex time
disk, uniformly on every stated initial polydisk. Uniform analytic Picard
iteration proves holomorphic dependence on those initial coordinates.
This establishes both signs of real flow, with no numerical trajectory
or unproved continuation interval invoked. D stays constant.

Since B(0)=0 and O is convex, `||B(u)||_*<=400||u||_*`
(in fact at most 4 on the entire outer domain, including complex states). For real states
and either sign of z, Gronwall and the variational equation give

```
||F_{±z}(u)||_* <= exp(400|z|)||u||_*,
||D F_{±z}(u)||_{*->*} <= exp(400|z|)
                         <=exp(.0032)<1.01,             (5)
```

on the verified trajectories. In particular the real inner patch
`||u||_*<=.001` has images within `||F_{±z}(u)||_*<.00101<.003`
for every `|z|<=8ell`, hence also for `|z|<=6ell`. Those images lie in
the independently constructed inverse domain. Uniqueness gives
`F_{−z}(F_z(u))=u`, so F_z is a diffeomorphism onto its image with inverse
F_{−z}. The derivative of that inverse satisfies the same (5).

For two inner states the straight segment between their images lies in
the .003 ball. Integrating the forward derivative along the original
segment and the reverse derivative along the image segment gives

```
exp(−.0032)||u−u'||_* <= ||F_z(u)−F_z(u')||_*
                       <= exp(.0032)||u−u'||_*.          (6)
```

The inverse statement concerns the image, and does not assert that the
.001 ball maps onto itself. B(0)=0 also implies F_z(0)=0 exactly, so this
is the retained zero-background chart and preserves the locus on which
physical perturbations vanish.

## 4. State and normal-coordinate jets through total order eighteen

Fix real `|z|<=6ell` and real `||u||_*<=.003`. A complex z circle of
radius ell and the initial coordinate polydisk of radius `.001` lie
strictly within the domain just constructed: `6ell+ell<8ell`.
On that complex domain `||F_z(u)||_*<.01`, so even the loose bound 2
is valid. Banach-valued Cauchy estimates give, for every coordinate
multi-index alpha,

```
||partial_z^m partial_u^alpha F_{±z}||_*
 <=2*m!*alpha!*10^{6m}*1000^{|alpha|}.
```

Every scalar coordinate of a unit *-norm direction has modulus at most 1.
Expanding k Fréchet directional derivatives has at most `16^k` terms,
and `alpha!<=k!`. Therefore, for both the chart and its inverse,

```
||partial_z^m D_u^k F_{±z}||_*
 <=8*m!*k!*10^{6m}*16000^k <10^126,     m+k<=18.          (7)
```

The extra factor four is harmless rounding. All 190 final integer
inequalities are verified exactly by the checker. The derivatives in (7)
are normal-coordinate z jets at fixed state and Fréchet state jets; they
are not eighteen solution derivatives, nor eighteen iterates of an
unquantified solution-dependent physical unit normal. The sign change
F_{−z} only changes signs in z derivatives, so the same Cauchy estimate
applies. No differentiation of a numerical flow approximation is used.

## 5. Exact Neumann rows, sixteen-field primal equivalence, and reconstruction

On the audited collar choose the local retained flow parameter `theta(z)=z`
and define `u(z)=F_z(V(z))`, keeping D unchanged. At z=0,
`F_0=Id`, `D F_0=Id`, and `partial_z F_z|_0=B`. The chain rule gives
exactly

```
u_D|_0=0,  u_N,z|_0=B_N(u_N|_0)
    if and only if V_D|_0=0,  V_N,z|_0=0.                (8)
```

There is no discarded nonlinear boundary remainder. These are the ten
Israel rows, one scalar row and one normal relative-harmonic row, with
four coordinate Dirichlet rows. Appending four tangential harmonic rows
would overprescribe the primal system. Normalization of the scalar row
never divides by sigma', so its possible zero inside the patch is harmless.

For every allowed physical or frozen iterate with D trace zero, the
inverse metric has `g^{za}=0` at the boundary, for all four tangent indices
including time. A fixed positive fiber metric and orthogonal D/N split
therefore give zero full principal flux: the D time trace is zero and the
N normal trace is zero. This states boundary cancellation, not a complete
energy/coercivity or commutator estimate.

The reduced metric/scalar equations have one common scalar metric wave
principal symbol on fifteen metric components and one scalar. The
metric-to-u map and u-to-V flow are invertible pointwise maps. Their
Jacobians multiply the second derivatives after substitution; multiplying
by the inverse Jacobians leaves the actual metric as the common principal
coefficient on all sixteen fields. Chain-rule remainder terms contain at
most first derivatives of V and explicit coordinate/chart derivatives.
Thus (8) is exactly equivalent to the original reduced primal problem
on the chart. Gauge propagation still uses the retained subsidiary proof
and its physical initial constraints; no additional primal gauge rows
are introduced here.

To remove any lapse or normalization ambiguity, reconstruction in the
actual fixed conformal-z chart is

```
g_zz=rho(z)² n²,
g_za=−rho(z)rho_b D_a,
h_ab=rho(z)² J_ab+rho_b² Htilde_ab,
phi=phibar(z)+w.                                         (9)
```

The minus sign comes from dy=−rho(z)dz. Dividing (9) by moving rho(z)²
in the coframe `(dz,dt,cosh(t)theta)` gives

```
Gamma_zz=n²,
Gamma_za=−(rho_b/rho(z))D_a,
Gamma_ab=J_ab+(rho_b/rho(z))²Htilde_ab.                    (10)
```

Here `n=sqrt(g_yy)` is a coordinate radial square root, not temporal ADM
lapse. Off D=0 it is not the Schur-complement lapse of radial hypersurfaces.
Equations (9)–(10), rather than a mistaken moving normalization of H in
(2), connect the chart to the separate metric energy certificate.
For inner images the background-normalized matrix perturbation is bounded
by `2(.00101)+(.00101)² +(131/129)²(.00101)
 +(131/129)(.00101)`, which is below 1/40 even with this redundant sum;
the energy author's sharper block estimate may be used instead.

The old source permits a compact collar parameter theta. This proof
chooses an admissible local extension with `theta=z` on the audited collar,
as specified in the coordinate dictionary. Extending
that choice to a global cutoff requires charging the cutoff derivatives;
no claim about such uncomputed derivatives is made. Since F_z(0)=0,
physical background and vanishing perturbation support are preserved
throughout the present collar. No new physical family, time evolution,
or observer congruence is introduced by the auxiliary flow.
