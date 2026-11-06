# Actual scalar-curvature source and a conditional compensation test

6 November 2026 UTC. Conditional mathematics for the registered ell=10^-6
first field, physical shell, original anchor and transported round labels.
Amplitude subscripts denote derivatives at e=0; the coefficient of e² is
one half the second derivative. These are model-normalized coefficients,
without SI calibration, physical energy or force interpretation.

## 1. Observation and domain hypotheses

Let q=-tau/rho_b, T=25ell/16, Y=n1²-n2² on the unit S3. Fix the original
proper-time origin and normal intrinsic shell geodesics, with labels n
transported from the same physical anchor. Assume a C² Einstein–scalar
shell family and its observer map exist smoothly on one common target
cylinder 0<=q<=T before failure of the normal flow. Retain the full physical
matching, reconstruction, constraint, gauge and uniqueness assumptions of
the pinned sources. This is an additional prescribed-cylinder hypothesis,
not an evaluated nonlinear existence time or amplitude radius.

The registered binding's shorter 3ell/2 coordinate slab by itself does not
establish this longer proper-time cylinder. The retained first *linear*
master and its actual source estimates cover the larger range through 3ell.
These two domain statements are distinct. The common nonlinear cylinder,
clock drift, embedding and uniform parameter estimates remain open inputs.

In physical Gaussian observer coordinates write

```
X_e^*h_e=-rho_b² dq²+k_e(q,n),
k0=rho_b²cosh²(q)kappa, G_e=k0^-1 k_e, j_e=sqrt(det G_e),
<f>=(2pi²)^-1 integral_S3 f dOmega, <f>_w=<Y²f>/<Y²>,
D_e/D0=<j_e>_w, D0>0.
```

The background rho_b is fixed with respect to e. Retain the exact quiet
shell/observer collar through qD=ell/4, including the original anchor, so
all initial amplitude variations in the area equation vanish. The first
linear field has the longer quiet interval q<=ell. A completion comparison
must retain the same full first physical jet and this exact quiet collar,
not just one coordinate scalar component. No freedom to modify the fixed
registered completion is asserted.

## 2. Full intrinsic curvature from the bulk equations

The pinned physical model and one-copy outward shell laws are

```
Ric5=dphi tensor dphi+(2/3)V(phi)g,
K_out=(sigma(phi)/6)h, N_out(phi)=-sigma'(phi)/2,
W=1-phi+phi³/3,
V=(W')²/2-(2/3)W²,
sigma=2W+delta(1+c_model phi),
delta=1/1000, c_model=96714321425692/161839258930731.
```

These are regular one-sided bulk traces, not distributional thin-shell
curvature. Use the curvature convention R_ABCD=k(g_ACg_BD-g_ADg_BC)
for constant sectional curvature k. Let X=partial_q X_e have squared
length -rho_b² and let E_ab=C5(N,a,N,b) be the electric bulk Weyl tensor
on the shell. Gauss and the five-dimensional Weyl decomposition give

```
Ric4(X,X)=(2/3)Ric5(X,X)+(rho_b²/3)Ric5(N,N)
          -(rho_b²/12)R5-E(X,X)-3(sigma/6)²rho_b².
```

Put P(phi)=sigma'²/16-V/2-sigma²/12. The full physical scalar pullback
is psi_e=phi_e evaluated on the moving shell and actual observer event.
The full Weyl contraction is mathcal E_e=C[g_e](N_e,X_e,N_e,X_e) at that
same event. Substituting the scalar field equation and shell laws yields
the exact intrinsic identity

```
R_e:=Ric_(h_e)(X_e,X_e)
   =(3/4)(partial_q psi_e)²-(rho_b²/12)|grad_(k_e) psi_e|²
     +rho_b² P(psi_e)-mathcal E_e.
```

In the intermediate bulk formula the potential term is -rho_b²V/2,
not -rho_b²V/3. The coefficient follows from -4/9+2/9-10/36=-1/2.
As a check, a constant-curvature bulk with Ric5=lambda g gives
Ric4(X,X)=-3rho_b²(lambda/4+(sigma/6)²). The actual warped background
is conformally flat, mathcal E0=0, psi0=phi_b is constant intrinsically,
and the background constraint gives P_b=-3/rho_b² and R0=-3.

Write the actual first physical scalar as psi1=S(q)Y. Its first pullback
has no clock/label term because psi0 is constant on the shell. Differentiating
twice, including the full normal, tangent, event and shell identification,
and using the exact sphere moments gives

```
<R2>_w=K_R(q)+Z(q),
K_R=(27/40)(S_q)²-(1/5)sech²(q)S²+(9rho_b²/20)P''_b S²,
Z=rho_b²P'_b <psi2>_w-<mathcal E2>_w.
```

Here <Y²>_w=9/20 and <|grad_kappa Y|²>_w=6/5. The squared first
gradients use only k0 because all intrinsic background scalar derivatives
vanish. In contrast, psi2 and mathcal E2 are completed physical pullbacks.
Bare coordinate second scalar or Weyl components cannot replace them.

The same actual background sources give P'_b=-p_b/(2rho_b mu_b), with
unscaled mu_b=(rho_b B_b)^-1 and p_b<0. Thus 3<rho_b²P'_b<5.
The unknown second scalar projection cannot be dropped merely because
the initial coordinate scalar was prescribed affinely. It can cancel
against the unknown second electric-Weyl projection; neither is bounded here.

Only three fixed angular test functions enter. With z=n1²+n2² set
H2=z/2-1/4 and H4=Y²-z/2+1/12. Then

```
Y²=1/6+H2+H4,
Delta_S3 H2=-8H2, Delta_S3 H4=-24H4,
<H2²>=1/48, <H4²>=19/720,
<f>_w=<f>+6<H2 f>+6<H4 f>.
```

This selects specific scalar projections in degrees 0, 2 and 4. It does
not truncate the nonlinear solution, all modes in those multiplets or raw
Weyl tensor slots before the physical pullback is formed.

## 3. Evaluated actual first-scalar source

The [scalar-source derivation](SCALAR_SOURCE_BOUNDS.md) transfers the
retained actual future waveform to this same nonzero-velocity past state.
It uses the complete positive endpoint-mass operator and its graph energy,
retains all endpoint terms and the original affine normalization, and
does not rerun the older PDE enclosure. On 0<=q<=T it obtains

```
|S|<51500, |S_q|<3107000000000,
-766143/10^9 < P''_b < -766128/10^9,
S(q)=0 for q<=ell,
S(5ell/4)<-37376.
```

The derivative is with respect to dimensionless q. These bounds are
continuous, conditional on the retained analytic source estimates and
linear reconstruction. Pointwise K_R then has a coarse strict enclosure
(-16000000000,6520000000000000000000000). This upper bound is not an
estimate of its typical value, an observed effect or a full curvature bound.

For the target T define the positive area Green kernel

```
W_T(s)=cosh²(s)[tanh(T)-tanh(s)],
Cosh=(1-T²/2)^-2 > cosh²(T), D=T-ell=9ell/16.
```

Let q*=5ell/4. Since S(ell)=0, Cauchy–Schwarz and the actual witness imply

```
integral_ell^q* S_q² ds > 37376²/(q*-ell),
W_T(s)>=(T-q*)/Cosh on ell<=s<=q*,
(27/40)integral_0^T W_T S_q² ds > (27/32)*37376²/Cosh.
```

The negative S² terms are bounded by a positive Cmax satisfying
Cmax>=1/5+(9rho_b²/20)|P''_b|. The first support gap gives

```
integral_0^T W_T Cmax S² ds <= Cosh*Cmax*sup|S|²*D²/2.
```

The final exact arithmetic uses the sharper retained rational bounds:
the positive term exceeds 1178689535.997, while the negative allowance is
below 0.002524. Consequently the **known scalar-curvature source integral**
satisfies the strict bound

```
I_R:=integral_0^T W_T(s)K_R(s) ds > 1170000000.
```

Equivalently its zero-data particular contribution to the weighted second
relative log-area derivative is less than -1170000000. This statement
concerns one exact source decomposition, not the full second area or its
sign. Affine amplitude rescaling changes these coefficients; the registered
amplitude is essential. Large derivative coefficients are compatible with
small first metric values and do not by themselves invalidate a small-
amplitude expansion or provide a numerical amplitude interval.

## 4. Conditional compensation criterion

The pinned previous proper-area theorem gives, at q=T,

```
D2/D0=K_D-I_R-I_Z,
I_Z=integral_0^T W_T(s)Z(s) ds,
-3.65e-16<K_D<5.357e-16.
```

These are the direct Gauss and Raychaudhuri decompositions; do not add an
alternative determinant correction or an independently shifted Bianchi
source decomposition to them. The known scalar term and unknown Z together
constitute the full curvature contribution once each.

Let Zmax=sup_(0<=q<=T)|Z(q)|, still unevaluated. Only the retained full
nonlinear quiet interval through ell/4 applies to Z. Therefore

```
|I_Z|<=C_R Zmax, C_R=(1/2)Cosh*(T-ell/4)²<8.614e-13.
```

It follows that the **additional sign hypothesis** D2(T)/D0(T)>=0 would
require

```
-I_Z>1170000000-5.357e-16,
Zmax>1.35e21.
```

The last inequality uses a conservative strict outer threshold. This is a
necessary budget conditional on that sign hypothesis, not a lower bound
on the actual unknown Z without it. As a useful sufficient conditional
test, a future actual bound Zmax<=10^21 would instead imply D2(T)/D0(T)<0.
Neither input nor resulting total sign is established in this checkpoint.

Positive physical area D_e>0 does not require D2>=0. A negative second
derivative is consistent with positive area for sufficiently small e.
Thus compensation is not a requirement for physical positivity, a
singularity theorem or a claim of net focusing. A second derivative at
zero cannot replace a uniform finite-amplitude remainder bound.

## 5. Remaining data and nearest known methods

The contracted Bianchi identity for this scalar model gives

```
nabla^a E_ab=[(2/3)Box_h psi+P'(psi)]psi_b-(1/12)partial_b Q,
Q=h^-1(dpsi,dpsi).
```

This is a conservation consequence, not an intrinsic constitutive law for
the missing bulk anisotropic Weyl stress. The known general nonclosure
of projected brane equations is discussed in the primary literature; the
scalar coupling and shell laws here must be derived in the present model.
No vacuum, Randall–Sundrum or observational constants are imported.

The earlier explicit second jets, prepared five-mean responses and
robustness classes use different prescriptions or domains and do not
supply the current psi2 or electric-Weyl bound. The weighted constraint
inverse and its numerical choices/constants must be instantiated or kept
as explicit inputs, not silently replaced by an older preparation. The
initial zero coordinate scalar second pair does not close its evolved
physical pullback.

The next targets are the three second-scalar projections and the two
nonconstant scalar anisotropic-Weyl stress projections, together with
the actual initial inverse and common development domain. Uniform mixed
metric, scalar, observer/clock, embedding, area and denominator estimates
on an admitted amplitude interval remain necessary for the original
nonlinear remainder and central-zero question. No new Big Bang origin,
physical force, journal acceptance, external peer review or observational
confirmation follows from this conditional calculation.
