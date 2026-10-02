# Actual-root and metric-response prerequisite audit

2 October 2026. This is an independently derived analytic audit of the saved
stationary and matched-stress checkpoints. It performs no new root search,
physical mode evolution, quantum source evaluation, or numerical response
experiment. The accompanying exact symbolic checks are algebra and synthetic
wrong-formula witnesses. They do not establish numerical accuracy, full metric
response, coupled stability, heating, or a physical cosmological history.

## 1. The response background must be the selected root

The saved stationary common-action model selects fixed parameters
`r=2z0`, `phi_b0`, `b`, and `gamma`, and solves a regular bulk branch with
endpoint `(phi*, z*=H*²)`. Its quadratic mass law is

```
x*=r [1+b(phi*-phi_b0)/2]²,
x_phi*=b r [1+b(phi*-phi_b0)/2],
x_phiphi*=b²r/2.
```

These parameters remain fixed under perturbation. In particular, `r` does not
become `2H*²` or `x*`. The physical background is not the H=1 calibration
without an explicit and consistently dimensional coordinate rescaling.

The following values are transcribed from existing saved producer rows
(`b=+1`, resolution `refined`), without recomputing any source or root:

| gamma | saved H*² | saved x* | saved u*=x*/H*² |
|---:|---:|---:|---:|
| 0 | 5.924014794328956e-5 | 1.1848029588657912e-4 | 2 |
| 0.01 | 5.924014794480220e-5 | 1.1848029588657418e-4 | 1.9999999999488487 |
| 1e6 | 5.939257863457309e-5 | 1.1847979917229192e-4 | 1.9948586489444573 |

The gamma=0.01 row passes the project's sampled hierarchy diagnostic. The
gamma=1e6 row fails it and is only a mathematical stress case. The tiny
gamma=0.01 departure from 2 does not establish an error bound on replacing
the actual-root propagator by plane waves: that would require bounds on the
propagator, subtraction, history, and momentum integration as a function of
u. The large-amplitude case provides a clear reason not to identify the two
propagators at all.

## 2. The arbitrary-root BD mode and its causal Green function

For a minimal scalar on a prescribed flat de Sitter chart,

```
a(eta)=-1/(H* eta), eta<0,
v_k=a chi_k,
Omega_k²=k²+a²x*-a''/a=k²+(u*-2)/eta²,
nu*²=9/4-u*,
v_k(eta)=sqrt(-pi eta)/2 exp[i pi(2nu*+1)/4]
                         H_nu*^(1)(-k eta).
```

Choose the usual real positive nu for `0<u*<=9/4` and positive imaginary nu
for larger u. For imaginary nu the displayed complex phase has modulus
`exp[-pi Im(nu*)/2]`; dropping this factor loses canonical normalization.
The positive-frequency past asymptotic and Wronskian are

```
v_k~exp(-ik eta)/sqrt(2k),
v_k v_k*'-v_k' v_k*=i.
```

At u*=2, nu*=1/2, the phase cancels the `-i` in the half-order Hankel function
and gives exactly the inherited plane wave. Positive x* defines the
Euclidean/BD state; resetting an instantaneous vacuum at a finite eta_i
does not define that state.

For eta>eta', the retarded Green function is

```
G_R,k(eta,eta')=i[v_k(eta)v_k*(eta')-v_k*(eta)v_k(eta')],
G_R,k(eta',eta')=0,
partial_eta G_R,k|eta=eta'+=1,
delta v_k(eta)=-integral_eta_i^eta G_R,k(eta,eta')
                                     delta Omega_k²(eta')v_k(eta') deta'.
```

The fundamental Green function is fixed by the mode equation and Wronskian;
the composite variance/stress susceptibility also depends on the state.
Gaussian occupation or coherence changes composite response even when the
fundamental retarded Green function is unchanged. The initial state must be
held fixed, with source and metric perturbations zero near eta_i.

At fixed geometry, the independent Kubo and forced-mode integrands coincide:

```
delta |v_k(eta)|²
  =-2 integral_eta_i^eta delta Omega_k²(eta')
                  Im[v_k*(eta)² v_k(eta')²] deta'.
```

For plane waves this is `-integral delta Omega² sin(2kDelta)/(2k²)`, the
already validated sign and normalization. For metric sources the canonical
frequency can depend on k; the forcing formula still applies per mode, but
it is not a single mass-only composite kernel.

## 3. Lapse and scale-factor channels derived before gauge fixing

Let the prescribed perturbed homogeneous geometry be

```
ds²=-N(eta)² deta²+A(eta)² dX²,
N=a(1+n), A=a(1+h), x=x*+delta x,
F=sqrt(A³/N)=a(1+c), c=(3h-n)/2,
v=F chi.
```

The exact canonical frequency is

```
Omega_k²=N²(k²/A²+x)-F''/F.
```

Its first variation is

```
delta Omega_k²=2k²(n-h)+2a²x* n+a²delta x-c''-2L c',
L=a'/a.
```

This follows directly from the kinetic coefficient `A³/N` of the physical
minimal action. Replacing c by h before imposing a conformal gauge drops
lapse and kinetic contacts.

In conformal gauge n=h, the source reduces to

```
delta Omega_k²=a²delta x+2a²x* h-h''-2L h'.
```

At the exact calibration x*=2H*² and L=-1/eta,

```
delta Omega_k²=a²delta x+4h/eta²-h''+2h'/eta.
```

This is a useful bounded next calibration source. It is not yet the complete
renormalized metric Hessian: explicit variations of observables and all
subtractions remain necessary.

## 4. An exact diffeomorphism null test

Use the infinitesimal time-relabeling convention

```
n -> n-T'-L T,
h -> h-L T,
delta x -> delta x-x*' T=delta x.
```

A pure time coordinate perturbation has `n=-T'-L T`, `h=-L T`. For arbitrary
a(eta) and constant x*, the canonical source identity is

```
delta Omega_k²=-T (Omega_k²)'-2T' Omega_k²-T'''/2.
```

The forced-mode solution, when T and its derivatives vanish before the source
window, is exactly

```
delta v_k=-T v_k'+T'v_k/2.
```

The `+T'v/2` term is a canonical density contact. Treating v as a scalar
would omit it. With physical variance `Q=F^-2 integral dmu |v|²`, the explicit
`-2c Q` factor then gives

```
delta Q=-T Q'.
```

For the static covariantly renormalized root, Q0 is constant and the complete
pure-gauge response must vanish. The renormalized rho0 and p0 are also
constant, giving zero pure-gauge scalar density/pressure response.

**Important finite-band qualification:** a fixed comoving momentum cutoff
generally makes the subtracted de Sitter Q_K and rho_K time dependent.
Therefore the finite-K covariance test is `delta Q_K=-T Q_K'`, not zero.
Declaring zero at finite K would misidentify genuine cutoff behavior as a
gauge failure. Covariant continuum subtraction and finite-band remainder
tests must be distinguished.

Equivalently in cosmic time, if lapse is `1+n_t` and spatial fractional
perturbation is h, the proper Hubble variation is

```
delta H=dot h-H* n_t,
delta R=6 delta dot H+24H* delta H.
```

Both are gauge invariant on a constant-H background. Constant shell phi*
makes delta phi itself invariant under tangential time reparameterization.
Gauge fixing the lapse before variation discards the Hamiltonian constraint
in a coupled solve; the lapse equation must first be retained.

A nontrivial compact time relabeling does not stay in conformal gauge:
its pure-gauge `n-h=-T'`. Therefore a numerical pure-gauge fixture must retain
the general lapse channel or use an explicit coordinate pullback. Setting
`n=h` and `h=-L T` for a varying compact T describes a physical perturbation,
not a null test.

## 5. Explicit observable and subtraction contacts

In conformal gauge define `D_k=v_k'-L v_k`. The unrenormalized minimal
bilinears are

```
e_k=[|D_k|²+(k²+a²x*)|v_k|²]/(2a⁴),
p_k=[|D_k|²-(k²/3+a²x*)|v_k|²]/(2a⁴),
Q_k=|v_k|²/a²,
delta D_k=delta v_k'-L delta v_k-h'v_k.
```

Their direct first variations are

```
delta e_k=-4h e_k+{Re[D_k* delta D_k]
  +(k²+a²x*)Re[v_k*delta v_k]
  +(a²delta x+2a²x*h)|v_k|²/2}/a⁴,

delta p_k=-4h p_k+{Re[D_k* delta D_k]
  -(k²/3+a²x*)Re[v_k*delta v_k]
  -(a²delta x+2a²x*h)|v_k|²/2}/a⁴,

delta Q_k=[2Re(v_k*delta v_k)-2h|v_k|²]/a².
```

These formulas contain the explicit metric prefactors, kinetic connection
contact and mass contacts. Feeding the effective source into the validated
mass-only closed stress formulas does not retain them.

Fixed positive r still varies inside the metric-dependent reference
frequency. In conformal gauge put

```
w=sqrt(k²+a²r), Delta=a²(x*-r)-a''/a,
U2=Delta/(2w)-w''/(4w²)+3w'²/(8w³),
Q_sub,k=[1/(2w)-U2/(2w²)]/a²,
delta w=ra²h/w,
delta Delta=a²delta x+2a²(x*-r)h-h''-2Lh'.
```

Then

```
delta U2=delta Delta/(2w)-Delta delta w/(2w²)
  -delta w''/(4w²)+w''delta w/(2w³)
  +3w'delta w'/(4w³)-9w'²delta w/(8w⁴),

delta Q_sub,k=-2h Q_sub,k+a^-2[-delta w/(2w²)
                        -delta U2/(2w²)+U2 delta w/w³].
```

When h=0, these reduce exactly to the inherited fixed-geometry mass
subtraction `delta Q_sub,k=-delta x/(4w³)`. Setting delta w=0 with a metric
perturbation keeps a fixed canonical frequency rather than a fixed physical
reference mass and changes the model. The stress W2/W4 expressions require
their full analogous variation before combining bare and subtraction
integrands; separate divergent pieces do not define renormalized contacts.

The physical current is

```
delta j=(x_phi*/2)delta Q+(x_phiphi* Q0/2)delta phi,
delta x=x_phi* delta phi.
```

The second term is the quadratic mass-law contact. For `b=0` the current
channel vanishes, while the metric-metric quantum response remains.
Generating functional derivatives of densitized sources additionally have
measure contacts; physical j and physical rho should not be confused with
those densitized derivatives.

## 6. Ward and trace targets, including the actual-root anomaly

On a constant-x, constant-H de Sitter root with p0=-rho0, linearized energy
exchange in conformal time is

```
delta rho'+3L(delta rho+delta p)=Q0 delta x'/2.
```

The lapse and expansion variations multiply vanishing background rho0+p0
or derivatives of the constant background. They therefore do not contribute
to this scalar equation at first order. Passing this equation does not by
itself validate the metric contacts, pressure, or the off-diagonal kernels.

The inherited common-action trace identity gives

```
-delta rho+3delta p=(delta Q''+2L delta Q')/(2a²)
                        -x*delta Q-Q0delta x+delta A_r.
```

Let `M=x*-r-2H*²`. In proper cosmic time the matched anomaly variation is

```
delta A_r=1/(16pi²){M delta x-(M/6+H*²/90)delta R
 -(delta ddot R+3H*delta dot R)/30
 +(delta ddot x+3H*delta dot x)/6}.
```

This follows from direct variation of

```
a2(r)=[x-r-R/6]²/2+(Riem²-Ric²)/180+box R/30-box x/6,
Riem²-Ric²=-12H²(dot H+H²).
```

At the matched exact reference x*=r=2H*², the algebraic terms become
`-2H*²delta x+(29H*²/90)delta R`. The derivative terms remain. A metric
contact cannot be inferred from the static trace remainder alone because
derivatives of curvature enter the response. At finite K use the explicitly
varied finite-K trace remainder, not the continuum anomaly as an exact
finite-band identity.

## 7. The matrix and the bulk/shell boundary problem still missing

The homogeneous physical response contains mass/scalar, lapse and spatial
metric channels. In a compact notation it is a retarded matrix relating
`(delta rho,delta p,delta j)` to `(n,h,delta phi)`, including local
observable, measure, reference-subtraction and mass-law contacts.
One mass-only column cannot determine its metric columns. The appropriate
object is the in-in/closed-time-path response, not the sphere-restricted
Euclidean stationary Hessian.

The static derivatives `rho_z,rho_phi,j_z,j_phi` constrain specified
equilibrium source derivatives of the same action. They do not equal a
finite-time retarded kernel or prove its zero-frequency limit for a specified
switching history. Likewise the retarded matrix is not symmetric under
exchanging both operator labels and times: exchanging times maps retarded
to advanced response. Symmetry of the doubled contour Hessian must not be
misapplied to a single retarded kernel.

For the declared doubled bulk, one shell, positive-normal convention, the
dimensionless generic junction signs are consistent with

```
2(K_mn-K h_mn)=-sigma h_mn+gamma T_mn,
2 n^A partial_A phi=-sigma_phi-gamma j.
```

On an invariant background `K_mn=q h_mn`, `T_mn=-rho h_mn`, these reduce
to the pinned model's `q=(sigma+gamma rho)/6` and
`w=-(sigma_phi+gamma j)/2`. Their first variations require

```
2(delta K_mn-h_mn delta K-K delta h_mn)
  =-sigma delta h_mn-sigma_phi delta phi_b h_mn+gamma delta T_mn,

2 delta(n^A partial_A phi)
  =-sigma_phiphi delta phi_b-gamma delta j,

delta K=h^mn delta K_mn-K^mn delta h_mn.
```

The trace/contact contraction in delta K is necessary. These covariant
linearized junction equations still need bulk scalar/gravitational
perturbation equations, constraints, gauge choices, regularity/ingoing
conditions, and initial data. They are not solved by the stationary shooting
Jacobian in `(ell,y_b)`.

For a normal shell displacement zeta, the actual induced quantities include

```
delta phi_b=delta phi_bulk(y_b)+w_b zeta,
h_b=h_bulk(y_b)+q_b zeta.
```

A radial gauge displacement Y changes bulk values by `(-w_bY,-q_bY)` and
zeta by `+Y`, leaving these combinations invariant. Setting zeta=0 without
retaining the corresponding bulk gauge and normal derivative equations
omits physical moving-endpoint terms. Static endpoint sensitivity identities
provide useful checks, but not the dynamical brane-bending operator.

Homogeneous FRW response also cannot test the nonconformal/tensor parts of
the covariant metric-metric kernel or choose a Weyl-squared finite constant
that vanishes on the homogeneous backgrounds. Full stability must specify
the perturbation sector and the gravitational EFT operators it retains.

## 8. Fastest bounded contribution and stopping boundary

The fastest substantial next contribution is the exact homogeneous
canonical operator, explicit metric observable/subtraction contacts, and
pure-diffeomorphism test above, with independent algebraic review. They form
an actionable prerequisite, not a completed numerical metric response.

A separately public frozen experiment can then calibrate the homogeneous
metric column at the existing exact plane-wave reference using one prescribed
smooth metric pulse. Freeze lapse/gauge, pulse, source translation, output
times, cutoff/refinement ladders, direct stress/subtraction implementation,
independent route, pure-gauge diagnostic, finite-K trace/Ward treatment,
precision, gates and stopping rule before any physical computation. The
existing mass-channel implementation is a normalization check. It is not
the independent metric route if both implementations import its contact
algebra.

The separately required actual-root propagator can next be checked by Hankel
mode identities and independent evolution on the prescribed selected root,
with explicit fixed r and state. Only after those numerical metric and
shifted-root prerequisites and bulk/shell constraint matching are validated
can a coupled stability or initial-data task be considered. Completing this
document alone must leave the automation's metric prerequisite blocked.

## Reproduction and provenance

```
/workspace/hdblast-cloud-setup/venv-frw/bin/python verify_prerequisites.py --output CHECKS.json
/workspace/hdblast-cloud-setup/venv-frw/bin/python -O verify_prerequisites.py --output CHECKS_OPTIMIZED.json
```

`INPUT_PINS.json` records the exact source documents and saved producer rows
used in this audit. The script uses explicit exceptions, including under
Python optimization. Its checks cover the action-derived canonical
frequency, conformal/exact-reference limits, time-gauge frequency and mode
contacts, variance covariance, retarded Green jump, direct trace variation,
full order-two variance subtraction variation, Ward linearization, and radial
brane-bending invariants. Eight synthetic wrong-formula mutations must be
detected. No interpretation of these algebraic checks as numerical response
or coupled-gravity validation is allowed.
