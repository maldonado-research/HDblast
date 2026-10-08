# A conditional bulk/brane model with an exact response and a four-dimensional replica

8 October 2026. This is a new **conditional analytic testbed**, not a modification
of the existing HDBLAST scalar calculation. Its geometry is flat, its fields are
free apart from a quadratic boundary coupling, and gravity is nondynamical. It
does not derive a blast, the prescribed B or zB profiles, a cosmological
background, reheating, or a Big Bang. The elementary results below carry no
external novelty claim.

The useful advance is an action-derived joint restriction: a positive continuum
spectral density fixes both brane dissipation and Gaussian equilibrium noise.
The same construction exhibits an exact four-dimensional continuum-bath rival.
Consequently the restriction can test this action against named smaller model
classes; it cannot identify the number of spacetime dimensions by itself.

## 1. Model, conventions and hypotheses

Use natural units hbar=c_light=1, signature (+----), coordinates
X=(t,x1,x2,x3,y), and the fixed half-space y>=0. The brane is y=0 and the outward
normal from the bulk points toward negative y. There is exactly one spacelike
extra coordinate. This is a one-sided bulk, with no second Z2 copy or full-line
factor of two. The real bulk field is Phi and the real brane field is q.
The complete scalar action is

```text
S = 1/2 int_(y>0) dt d^3x dy
        [Phi_dot^2 - |grad_x Phi|^2 - Phi_y^2 - M^2 Phi^2]
  + int_(y=0) dt d^3x
        [1/2 q_dot^2 - 1/2 |grad_x q|^2 - 1/2 m^2 q^2
         - c/2 Phi_0^2 + g q Phi_0 + J q].                 (1)
```

Here c is a Robin parameter, not the speed of light; Phi_0 is the trace of Phi.
Assume M>0, m^2>0, c>0, g real and nonzero, and

```text
g^2 < m^2 (c+M).                                          (2)
```

There are no hidden boundary kinetic terms, counterterms or omitted scalar
contacts in (1). J is an independently declared infinitesimal probe used to
define susceptibility. Its specification is not a dynamical bulk source. For
classical statements take finite-energy initial data and smooth compactly
supported probes; use smooth solutions for the displayed integrations by parts
and their energy-norm closure otherwise. Spatial plane waves and sinusoidal
forcing denote Fourier components/distributions, not finite total-energy states.

Mass dimensions are [Phi]=3/2, [q]=1, [c]=[M]=[m]=1,
[g]=3/2 and [J]=3. Fourier conventions are exp(-i omega t+i k.x), with retarded
transforms analytic for Im omega>0. The classical bulk response has no incoming
homogeneous wave from infinity; explicitly prepared incoming waves/noise are
added separately. Smooth spacetime-smeared field observables are used in quantum
statements. Pointwise noise variance and renormalized vacuum stress are not
defined by the formulas below.

The quantum model is the Gaussian quantization of this positive quadratic
system. An equilibrium assertion means the stationary KMS state of the **full
coupled** system at inverse temperature beta>0, or its ground-state beta=infinity
limit. It does not mean an abruptly factorized bath and brane preparation. An
incoming-scattering description must populate any bound normal mode separately;
it agrees with full equilibrium only when that mode has the same temperature.

## 2. Theorem and proof: stability, retarded reduction and positive spectrum

**Theorem.** Under (1)-(2):

1. The classical finite-energy initial-value problem admits a positive
   self-adjoint spatial operator and a unique energy-space wave evolution.
   There is no negative kinetic-energy sector or tachyonic scalar normal mode.
2. With vanishing incoming mean bath field, the exact brane retarded response is

   ```text
   q(omega,k) = chi_R(omega,k) J(omega,k),
   chi_R^-1 = k^2+m^2-(omega+i0)^2 - g^2 G_R(omega,k),
   G_R = 1/(c+s_R),
   s_R = sqrt(k^2+M^2-(omega+i0)^2).                       (3)
   ```

   The square root is the branch analytic for Im omega>0 with positive real
   part there. At omega>sqrt(k^2+M^2), s_R=-i p with
   p=sqrt(omega^2-k^2-M^2)>0. This sign is fixed by waves traveling into y>0.
3. The induced kernel has the positive spectral representation

   ```text
   G_R(omega,k) = int_(M^2)^infinity
       rho_b(u) du / [k^2+u-(omega+i0)^2],
   rho_b(u) = sqrt(u-M^2) / [pi (c^2+u-M^2)].              (4)
   ```

   There is no bath bound state for c>0. For omega>threshold,

   ```text
   Im G_R = p/(c^2+p^2) > 0,
   Re G_R = c/(c^2+p^2).                                 (5)
   ```

4. The retarded kernel is causal and has memory. It is not replaced by an
   arbitrarily prescribed pulse or a local damping coefficient. Equation (3)
   is exact for this declared quadratic model, without a weak-coupling or
   derivative approximation.

**Proof, boundary variation.** Varying the bulk action yields the boundary term
+Phi_y(0) delta Phi_0 because the integration domain is y>=0. Combined with the
boundary action this gives

```text
(d_t^2 - Laplacian_x - d_y^2 + M^2) Phi = 0,
Phi_y(0) = c Phi_0 - g q,
(d_t^2 - Laplacian_x + m^2) q - g Phi_0 = J.               (6)
```

The Robin condition is an action-derived scalar boundary condition, not an
Israel gravitational junction condition.

**Proof, stability and well-posedness.** For a decaying smooth bulk profile,

```text
int_0^infinity (Phi_y^2+M^2 Phi^2) dy - M Phi_0^2
    = int_0^infinity (Phi_y+M Phi)^2 dy >= 0.              (7)
```

Let B be the integrated bulk spatial quadratic energy plus c int Phi_0^2,
without the factor 1/2. Then B>=(c+M)||Phi_0||^2. If
r=|g|/[m sqrt(c+M)]<1, Cauchy-Schwarz and 2ab<=a^2+b^2 imply

```text
|2g int q Phi_0| <= r [m^2||q||^2+B].                    (8)
```

The coupled spatial quadratic form is therefore bounded below by (1-r) times
the uncoupled positive form, including the brane gradient term. The reverse
upper bound follows by the same inequality. The trace map on the half-space
H^1 space is continuous; thus the coupled and uncoupled form norms are
equivalent, the form is closed and positive on
H^1(R^3) direct-sum H^1(R^3 x R_+), and its self-adjoint operator generates the
unique energy-space wave evolution. Smooth data in its operator domain satisfy
(6). Strict positivity excludes exponentially growing normal modes and implies
retarded analyticity in the upper frequency half-plane. This is scalar-model
stability; it does not analyze gravitational backreaction or another signature.

Condition (2) is the sharp zero-momentum positivity boundary. At equality a
zero-frequency k=0 generalized mode appears. If g^2>m^2(c+M), choosing
Phi proportional to exp(-My) and very slowly varying spatial profiles makes
the spatial form negative. The resulting negative spectral sector is unstable.
The strict stable domain, rather than the equality/unstable continuation, is
used for all claimed predictions.

**Proof, retarded boundary elimination.** For Im omega>0 the decaying solution
of the transformed bulk equation is Phi(y)=Phi_0 exp(-s_R y). Substitution in
(6) gives (c+s_R)Phi_0=gq. A nonzero incoming/homogeneous bath field instead adds
Phi_h,0, so that Phi_0=Phi_h,0+g G_R q. Substitution into the q equation gives

```text
[d_t^2+k^2+m^2]q(t,k)
  - g^2 int_(-infinity)^t G_R(t-t',k) q(t',k) dt'
    = J(t,k) + eta(t,k),
eta = g Phi_h,0.                                         (9)
```

Initial transients, state correlations and any surviving bound mode must be
included if the lower limit or preparation is changed. This equation does not
assert that arbitrary product initial data have stationary equilibrium noise.

**Proof, spectrum and causal kernel.** The normalized generalized Robin modes
for p>0 are

```text
u_p(y) = sqrt(2/pi) [p cos(py)+c sin(py)]/sqrt(p^2+c^2),
u_p'(0)=c u_p(0),
u_p(0)^2=(2/pi)p^2/(p^2+c^2).                             (10)
```

They are the continuous spectral transform of the nonnegative half-line Robin
operator; the putative bound profile exp(-a y) would require -a=c and therefore
does not exist for c>0. The identity underlying (4) is, for Re s>0,

```text
(2/pi) int_0^infinity p^2 dp/[(p^2+c^2)(p^2+s^2)]
   = (s-c)/(s^2-c^2) = 1/(c+s).                          (11)
```

It follows by partial fractions and int_0^infinity dp/(p^2+a^2)=pi/(2a).
At s=c the continuous limit is 1/(2c). Setting u=M^2+p^2 yields (4).
Equivalently,

```text
G_R(t,k)=theta(t) int_0^infinity dp u_p(0)^2
  sin(sqrt(k^2+M^2+p^2)t)/sqrt(k^2+M^2+p^2).              (12)
```

The integral is an oscillatory distribution, with smooth smearing when needed.
Every integrand is a retarded four-dimensional Klein-Gordon kernel, so the
brane kernel also vanishes outside the brane causal cone. The square-root
threshold proves nonlocal time dependence for M>0. No instantaneous influence
or superluminal propagation follows from eliminating the bulk. Equations
(5) follow directly from s_R=-ip. This completes the proof.

## 3. Bound modes and why their state cannot be omitted

The bath alone has no bound mode, but the **coupled** q/Phi system can have one.
Write z=omega^2-k^2 and, below the bulk threshold,

```text
D(z)=m^2-z-g^2/[c+sqrt(M^2-z)].
D'(z)=-1-g^2/[2 sqrt(M^2-z)(c+sqrt(M^2-z))^2] < 0.
```

Strict stability gives D(0)>0. There is exactly one positive bound mass squared
z_b in (0,M^2) iff m^2-M^2-g^2/c<0. If that quantity is positive there is none;
at zero the mode reaches threshold and is not an L^2 bulk bound state. Its
brane residue is

```text
Z_b = [1+g^2/(2 s_b (c+s_b)^2)]^-1 > 0,
s_b = sqrt(M^2-z_b).                                     (13)
```

For example, the manufactured parameters M=2, c=3, g=1, m^2=13/4 have
z_b=3, s_b=1 and Z_b=32/33. This stable mode does not decay into the continuum.
An incoming thermal bath alone cannot erase an independently prepared bound
mode. For a simple dissipative-only test design one can declare the stronger
domain m^2>M^2+g^2/c, which excludes it. The manufactured no-bound-mode example
M=2,c=3,g=1,m^2=9 lies in that domain.

Above threshold the q spectral density is nonnegative:

```text
rho_q(u) = g^2 rho_b(u) / |m^2-u-g^2 G_R(u+i0)|^2.
```

When a bound mode exists add Z_b delta(u-z_b). These positive spectral sectors
justify Gaussian ground/KMS covariances after smearing. No assumption about a
hot universe, thermalization mechanism or observed temperature is introduced.

## 4. Joint noise, response and signed energy flow

Define the symmetrized equilibrium spectral covariance by

```text
S_q(omega,k) = int dt exp(i omega t)
                 (1/2)<{q(t,k),q(0,-k)}>_connected,
```

with the momentum-conserving delta distribution stripped off. For the full
coupled KMS state, the spectral theorem for its positive oscillators gives

```text
S_q(omega,k) = coth(beta omega/2) Im chi_R(omega,k).         (14)
```

This includes bound-mode delta functions when present; at zero temperature
coth(beta omega/2) becomes sign(omega). For incoming Gaussian bath modes at
that temperature, the force noise in (9) has

```text
N_eta(omega,k) = g^2 coth(beta omega/2) Im G_R(omega,k).     (15)
```

N_eta is even in omega and nonnegative. On the continuum, away from isolated
bound poles, Im chi_R=g^2 Im G_R |chi_R|^2, hence
S_q=|chi_R|^2 N_eta. This continuum identity is **not** a prescription to drop
the bound-mode covariance. Equations (14)-(15) use the one-half
anticommutator convention; using an unsymmetrized covariance changes the Bose
factors. Gaussianity makes the mean and two-point functions determine higher
correlations by Wick's rule. A non-Gaussian bath is outside this contract.

For fixed k, |Im G_R|~1/|omega| at large frequency. The raw equal-time force
variance consequently has a logarithmic ultraviolet divergence even before
spatial coincidence is considered. Smooth spacetime test functions (including
time smearing) or a separately declared physical regulator are necessary;
spatial smearing alone does not cure this fixed-k divergence. This does not spoil the convergent
retarded integral (11), positivity of smeared covariances, or the classical
energy theorem. It prohibits presenting an unsmeared finite noise amplitude,
vacuum energy or stress tensor without additional work. No theorem requiring
compact spectrum or exponentially decaying bath density is used here.

For finite-energy classical fields define

```text
E_bulk = 1/2 int d^3x dy [Phi_dot^2+|grad_x Phi|^2+Phi_y^2+M^2 Phi^2],
E_q    = 1/2 int d^3x [q_dot^2+|grad_x q|^2+m^2 q^2],
E_R    = c/2 int d^3x Phi_0^2,
E_int  = -g int d^3x q Phi_0.
```

When no flux reaches an omitted outer boundary, (6) gives

```text
dE_bulk/dt = -int Phi_dot_0 Phi_y(0),
d(E_bulk+E_R)/dt = g int q Phi_dot_0,
dE_q/dt = int [J q_dot+g Phi_0 q_dot],
d(E_bulk+E_R+E_q+E_int)/dt = int J q_dot.                 (16)
```

For a finite bulk slab retain its outward flux term. This ledger includes the
interaction energy, so it does not assign an ambiguous instantaneous exchange
term to a uniquely separate brane reservoir. The positive-y energy flux is
S_y=-Phi_dot Phi_y. For a positive-frequency outgoing wave with q amplitude Q,
Phi(y)=gQ exp(ip y)/(c-ip), the time-averaged flux per brane area is

```text
<S_y> = omega g^2 p |Q|^2/[2(c^2+p^2)] >= 0.              (17)
```

It equals minus the average work of the self-force on q. Below threshold the
evanescent field has zero averaged propagating flux. A no-incoming retarded
experiment describes energy sent **from the driven brane into the bulk**.
Bulk-to-brane energy supply requires an independently prepared incoming state
or a dynamical source and its initial energy. Nothing in (17) derives such a
source or a cosmological event. Equilibrium noise is not a net blast.

## 5. A falsifiable restriction, and its exact limits

Assume an ideal calibrated probe of J and q, the canonical normalization in
(1), known measurement response, and the same parameters across frequencies
and momenta. These are prospective measurement premises, not an existing
HDBLAST data product. Define the continuum dissipative coefficient

```text
Gamma(omega,k) = -Im chi_R^-1(omega,k) = g^2 Im G_R(omega,k).
```

For omega>0 the model jointly requires

```text
Gamma = 0                                      if omega^2<k^2+M^2,
Gamma = g^2 p/(c^2+p^2)                         if p>0,
p/Gamma = (c^2+p^2)/g^2,
omega_threshold(k)^2-k^2 = M^2,
Re[g^2 G_R]/Gamma = c/p                         if p>0,
S_q/Im chi_R = coth(beta omega/2)                on regular equilibrium bands.
                                                               (18)
```

Thus p/Gamma is affine in p^2 with shared slope 1/g^2 and intercept c^2/g^2;
Gamma peaks at p=c with value g^2/(2c). The threshold, dispersive response,
absorption and thermal fluctuation factor cannot be chosen independently.
The response identifies g^2, not the sign of g, which changes under Phi->-Phi.
A failed restriction can reject this action or its calibration/state premises.
Agreement would support only the declared family relative to tested rivals.

**Restricted rival theorem.** A finite number of four-dimensional free fields
with a finite-dimensional local quadratic action, finite derivative order,
translation-invariant coefficients, the same probe/measurement coordinate and
no unmodeled continuum produces a rational retarded q response in frequency
at each fixed k. For the ordinary Lorentz-invariant second-order version the
response is rational in z=omega^2-k^2. Our nonzero-g, M>0 response cannot equal
such a rational response on an open regular frequency interval.

To see this, invert the putatively equal response where it is nonzero. Equation
(3) would make 1/(c+sqrt(k^2+M^2-omega^2)) rational. Solving for the square root
would make sqrt(k^2+M^2-omega^2) rational. Its square has simple roots at
omega=+/-sqrt(k^2+M^2), whereas every zero/pole of a squared rational function
has even order. Local equality extends analytically, giving a contradiction.
For intervals on the continuum, use the local analytic continuation of the
specified boundary branch away from the threshold; the same algebra applies.
The argument permits rational frequency-dependent readouts, but not arbitrary
nonlocal measurement filters that themselves supply the missing square root.

This is an exact ideal-function separation. It supplies **no proved positive
distance for a specified finite noisy experiment**, no likelihood, and no
exclusion of interacting four-dimensional theories, whose multiparticle spectra
can be continuous. Positive spectral quadrature can approximate the integral
on compact upper-half-plane domains or after suitable finite-time/smooth
measurement filtering; its error, UV tail and experimental resolution need
their own bound. An exactly closed finite bath has discrete real-axis spectral
delta functions: approximation is not uniform convergence of its pointwise
infinite-time absorption to the continuous density. A massless M=0,k=0 special case has
G_R=1/(c-i omega), so observing only that sector would even lose this
square-root-versus-rational distinction. The declared M>0 avoids that exception.

**Exact continuum rival.** Introduce four-dimensional scalar fields chi_p(x),
p in (0,infinity), with mass squared M^2+p^2 and the action

```text
S_4 = int d^4x [1/2 (partial q)^2 - 1/2 m^2 q^2 + Jq]
    + 1/2 int_0^infinity dp int d^4x
        [(partial chi_p)^2 - (M^2+p^2)chi_p^2]
    + g int d^4x q int_0^infinity dp u_p(0) chi_p.         (19)
```

This has infinitely many species labeled by p, with no geometric extra
coordinate. It is local in the four coordinates x, and the coupling is defined
as a quadratic form. Its stability follows from
int dp u_p(0)^2/(M^2+p^2)=1/(c+M) and (2). The unitary Robin spectral transform
Phi(x,y)=int dp u_p(y)chi_p(x) maps (1) into (19), including the Robin energy.
Thus matching the full Gaussian preparation (including bound modes), J,
measurement map and contacts gives **identical smeared q quantum states and
correlation functions, and identical outcome probabilities under the same
declared measurement protocol**. Noncommuting q operators at different times
are not assigned an ordinary joint path probability. Classical Gaussian
versions have identical laws as smeared Gaussian random distributions. Both
(14) and (18) are inherited by
this stable four-dimensional rival. A freely selectable temperature or noise
law alone is unnecessary: this microscopic replica realizes them together.

The model therefore advances the earlier identifiability argument: the
encompassing rival is now exhibited by an explicit positive quadratic action.
It does not assert that a finite-species four-dimensional theory has the same
exact response or that every arbitrary effective kernel has such a completion.

## 6. Scope of completion and next useful calculation

This prototype supplies an action, signature, exact flat embedding, scalar
matching, a sharp stability domain, a causal induced kernel, positive spectral
measure, a declared Gaussian state class, a signed energy ledger and a joint
response/noise restriction. It also gives a precise finite-free-sector rival
that cannot match the ideal response and a continuum rival that matches it
exactly. These are conditional analytic results; all numerical checks bundled
here are manufactured symbolic or diagnostic controls.

A prospective physical test still needs an actual probe/readout and sensor
calibration; chosen parameter units/supports and preparation; smooth measurement
filters; finite-window/UV error bounds; declared finite, interacting and
continuum rivals; an independent likelihood/holdout contract; and evidence that
the prototype is relevant to a physical system. None exists in this subtask.
For a cosmological origin model, gravity, bulk source dynamics, an evolving
background, backreaction, energy supply and thermalization remain additional
unsolved requirements. Appending the old B and zB pulse shapes would not solve
them. The old de Sitter scalar mass squared two is not M, m or a limit silently
taken in this model.

The next narrowly defined mathematical task is to predeclare smooth probe
filters and compare the finite-resolution smeared covariance/response to a
finite set of four-dimensional masses with a rigorous quadrature and tail
bound. That would quantify what precision and rival restrictions make (18)
testable, while retaining the exact continuum degeneracy. No physical source,
target or likelihood evaluation is authorized by this note. The separate
trajectory gate still requires its root freeze, independent review and immutable
public GO before such protected execution.

Current project statuses remain: historical metric calibration **FAIL**;
full continuous certificate **UNRESOLVED**; higher-dimensional Big Bang cause
**NOT_ESTABLISHED**; external mathematical novelty **NOT_ASSESSED**.
