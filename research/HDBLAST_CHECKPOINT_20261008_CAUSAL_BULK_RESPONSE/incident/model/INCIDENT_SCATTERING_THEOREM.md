# Incoming bulk energy gives transient brane motion and unit reflection

8 October 2026. Additive pure analytic continuation of the accepted flat scalar
half-space model. The accepted `physical-model` source files and manifest are
unchanged. This note specifies an incoming coherent excitation and follows its
energy. It does not derive the initial excitation, gravity, a cosmology,
thermalization, or any observed signal.

## 1. Same action, new preparation

Use the one-sided Minkowski half-space y>=0, signature (+----), with the scalar
action and conventions in `../physical-model/HALF_SPACE_RESPONSE_THEOREM.md`.
In particular, the boundary terms are -c Phi_0^2/2+g q Phi_0 and

```text
(d_t^2-Laplacian_x-d_y^2+M^2)Phi=0,
Phi_y(0)=c Phi_0-gq,
(d_t^2-Laplacian_x+m^2)q-g Phi_0=0.                       (1)
```

There is no brane probe: J=0. Assume M>0, c>0, m^2>0, real g and
0<g^2<m^2(c+M). The action is autonomous and quadratic. The positive total
energy and self-adjoint spatial operator established for that action continue
to apply. No additional brane channel, dissipative term, time-dependent
coupling or other interaction is introduced.

For a spatial Fourier coordinate k and p>0 define

```text
omega(p,k)=sqrt(k^2+M^2+p^2)>0,
D=k^2+m^2-omega^2=m^2-M^2-p^2.
```

The incoming wave travels toward the brane from large positive y. Its harmonic
coefficient is normalized to one:

```text
Phi=Re[A (exp(-ipy)+R exp(ipy)) exp(-i omega t+i k.x)],
q=Re[A Q exp(-i omega t+i k.x)].                          (2)
```

This infinite-duration plane wave is a Fourier component. Its fluxes are
time-averaged flux densities per brane area. It is not a finite total-energy
state. Finite-energy packets are specified in section 4. A coherent quantum
excitation means a Weyl displacement of the full coupled Gaussian ground state
in these continuum modes; its excess mean field follows (1). Any bound-mode
coherent occupation is zero. Ground-state fluctuations are not set to zero,
and no unsmeared vacuum energy or noise variance is used.

## 2. Exact reflection and brane amplitudes, including D=0

**Proposition 1.** For every real p>0 and g nonzero, (1)-(2) have the unique
stationary scattering amplitudes

```text
F=(ip-c)D+g^2,
R=[(ip+c)D-g^2]/F,
Q=2ipg/F,
1+R=2ipD/F.                                             (3)
```

They are finite on this open continuum and satisfy |R|^2=1.

**Proof.** The boundary and brane equations form the two-by-two linear system

```text
(ip-c)R+gQ=ip+c,
-gR+DQ=g.                                               (4)
```

Its determinant is F. For real parameters,

```text
|F|^2=(g^2-cD)^2+p^2D^2>0,                             (5)
```

because a vanishing second term with p>0 requires D=0, at which the first term
is g^4>0. Cramer's rule gives (3), without dividing by D. The numerator of R
is -F*, so R=-F*/F and |R|^2=1. This also proves that no real pole is hidden by
the algebra.

For D nonzero one may eliminate q first and write
Phi_y(0)=(c-g^2/D)Phi_0. That quotient is not the definition at D=0. If m^2>M^2,
the bare-brane frequency lies in the continuum at p=sqrt(m^2-M^2), and (3) gives

```text
D=0:       R=-1,       Q=2ip/g,       Phi_0/A=0.           (6)
```

The derivative boundary condition is still satisfied: ip(R-1)=-gQ=-2ip.
The brane amplitude is finite for fixed nonzero g. A zero trace Phi_0 does not
mean zero stored brane energy. This is the bare-brane frequency, not a claim
that the maximum of every susceptibility or packet observable occurs there.
When m^2<=M^2 there is no D=0 point with p>0.

For fixed D nonzero, g->0 gives Q->0 and
R->R_0=(ip+c)/(ip-c), the uncoupled Robin reflection. This limit is not uniform
at D=0: (6) has Q proportional to 1/g. At exactly g=0,D=0, (4) fixes R_0 but
does not determine the independent free brane oscillator amplitude. Choosing
zero initial brane excitation sets that amplitude to zero. Substituting
g=D=0 in (3) produces 0/0 and is invalid. Infinite-duration harmonic excitation
must not be used to infer a finite-energy packet divergence in this limit.

## 3. Response compatibility and complete energy accounting

The same incoming wave in the uncoupled Robin bath has

```text
R_0=(ip+c)/(ip-c),
Phi_free,0/A=1+R_0=2ip/(ip-c),
G_R=1/(c-ip),
chi_R=[D-g^2/(c-ip)]^-1.
```

The driven brane and full outgoing amplitude obey the exact identities

```text
Q=chi_R g (Phi_free,0/A),
R=R_0+gG_R Q.                                           (7)
```

Thus the incident calculation uses the accepted retarded response with the
proper Robin incoming field and its phase. Replacing Phi_free,0/A by one would
give the wrong normalization and phase. The direct reflected wave and the
brane-radiated wave interfere coherently.

The positive-y flux of the real harmonic field is S_y=-Phi_dot Phi_y. Its
incoming, outgoing and net time averages are

```text
<S_in>  = -omega p |A|^2/2,
<S_out> =  omega p |A R|^2/2,
<S_net> =  omega p |A|^2 (|R|^2-1)/2=0.                 (8)
```

Incoming/outgoing interference makes no contribution to this signed current.
There is only one open bulk channel at fixed k, and all incident stationary
flux returns in that channel. The stationary work on q is zero:

```text
<g Phi_0 q_dot>
   = -omega g |A|^2 Im[(1+R)Q*]/2=0,                  (9)
```

since DQ=g(1+R) with real D. Nevertheless the time-averaged bare brane energy
density is generally positive,

```text
<E_q>/brane volume = (omega^2+k^2+m^2)|A Q|^2/4.          (10)
```

Equation (10) is a plane-wave mean density, not global energy. During a finite
pulse, energy can build up in q and later return to the bulk.

For finite-energy fields the complete classical energy remains

```text
E=E_bulk+E_q+E_R+E_int,
E_R=(c/2)int d^3x Phi_0^2,
E_int=-g int d^3x q Phi_0,
dE/dt=0,                                               (11)
```

where E_bulk and E_q are the positive kinetic, gradient and mass energies of
the accepted model. For a slab 0<=y<=Y, the corresponding signed ledger is
dE_slab/dt=-int d^3x S_y(Y), including E_R and E_int and assuming no flux through
the omitted spatial faces. The bulk boundary sign gives

```text
d(E_bulk+E_R)/dt=g int q Phi_dot_0,
dE_q/dt=g int Phi_0 q_dot,
dE_int/dt=-g int(q_dot Phi_0+q Phi_dot_0).                (12)
```

The separate outgoing radiation gG_RQ carries positive flux. Counting that
flux as irreversible absorption of the incident wave would omit its
interference with R_0. Equations (7)-(9), or the full ledger (11), prohibit that
interpretation in this closed one-channel model. The incoming preparation is
the energy supply; the action does not create that supply from nothing.

## 4. Finite continuum packets leave no persistent brane excitation

**Proposition 2.** Fix k and choose a(p) in C_c^infinity((p_0,p_1)) for finite
0<p_0<p_1. Define the complexified solution

```text
Phi_k(t,y)=(1/sqrt(2pi)) int dp a(p)
        [exp(-ipy)+R(p)exp(ipy)]exp(-i omega(p,k)t),
q_k(t)=(1/sqrt(2pi)) int dp a(p)Q(p)exp(-i omega(p,k)t).
                                                               (13)
```

Real and imaginary parts separately solve the real equations. At each finite
time the solution has finite energy in the half-line Fourier fiber. Its bound
spectral projection is zero. For every positive integer N and finite Y,
q_k, q_dot_k, Phi_k, Phi_dot_k and Phi_y,k are O(|t|^-N), uniformly for
0<=y<=Y as t tends to either plus or minus infinity. Every quadratic local
energy component, including the boundary and interaction terms in absolute
value, tends to zero. The full half-line energy is conserved and is carried
by an outgoing continuum packet at late times.

**Proof of decay and finite energy.** On this compact support, (5) makes R,Q
and all their derivatives bounded. At finite t the Fourier integrals and their
y derivatives are rapidly decreasing as y tends to infinity, giving finite
fiber energy. Also

```text
omega'(p)=p/omega(p,k)>=v_min>0.
```

For any compactly supported smooth amplitude f, integration by parts has no
endpoint term and yields

```text
int f(p)exp(-it omega(p))dp
   =(1/(it))int Lf(p)exp(-it omega(p))dp,
Lf=d_p[f/omega'],
|int f exp(-it omega)dp|<=|t|^-N ||L^N f||_1.             (14)
```

Use f=aQ, -i omega aQ and a[exp(-ipy)+Rexp(ipy)], with the corresponding
time/y derivative factors. For y in a fixed finite interval all resulting
L1 norms are bounded uniformly. This proves the stated rates. Equivalently,
a change from p to omega gives a smooth compact frequency density and the
Riemann-Lebesgue conclusion; (14) supplies the stronger rapid-decay estimate.
No such rate is claimed for packets touching p=0, nonvanishing endpoint
amplitudes, nonsmooth profiles, or constants uniform as g tends to zero.

For clarity a nonnegative local majorant is

```text
E_loc,abs=1/2[|q_dot|^2+(k^2+m^2)|q|^2
  +int_0^Y (|Phi_dot|^2+|Phi_y|^2+(k^2+M^2)|Phi|^2)dy
  +c|Phi_0|^2+2|g||q||Phi_0|].                          (15)
```

It is O(|t|^-2N). The actual signed local interaction energy is bounded in
absolute value by this expression. This avoids assuming positivity of an
arbitrarily truncated interaction-energy assignment.

**Asymptotic outgoing energy.** In the complex fiber norm, use

```text
E_k=1/2[|q_dot|^2+(k^2+m^2)|q|^2
 +int_0^infinity(|Phi_dot|^2+|Phi_y|^2+(k^2+M^2)|Phi|^2)dy
 +c|Phi_0|^2-2g Re(q*Phi_0)].                           (16)
```

This is the sum of the real-part and imaginary-part energies. The incoming
component of (13) lies asymptotically in y>0 as t->-infinity; the outgoing
component lies there as t->+infinity. On the wrong half-line their phases have
derivative of magnitude at least |y|+v_min|t|. Repeated nonstationary-phase
integration gives vanishing half-line energy norms for those components.
On the correct half-line their energies approach the full-line Plancherel
energies. The cross term vanishes by Cauchy-Schwarz, and the boundary terms
vanish by (14). Consequently

```text
E_k(-infinity)=int dp omega^2 |a(p)|^2,
E_k(+infinity)=int dp omega^2 |a(p)R(p)|^2
             =E_k(-infinity).                          (17)
```

These conventions fix the Fourier normalization in (13). Real physical
components obey their own identical conservation law; (17) uses the declared
complex fiber norm. No global bulk energy disappears when local energy decays.

**Global brane energy qualification.** A fixed k labels a Fourier fiber, not a
normalizable plane wave in the three brane coordinates. To claim a finite
global brane energy, take a(k,p) smooth and compactly supported in both k and
p, with |k|<=K and a common 0<p_0<p_1. Then
omega'>=p_0/sqrt(K^2+M^2+p_1^2)>0 uniformly. Applying (14) and Plancherel in k
proves that the integrated brane energy and the entire finite-y-layer energy
majorant decay as O(|t|^-2N). The total bulk-plus-brane energy remains constant,
with the same incoming/outgoing equality integrated over k. A real physical
packet is obtained by taking real parts or the usual conjugate Fourier data.

## 5. Why a continuum packet cannot permanently populate a bound mode

The coupled model can have a bound mode with mass squared z_b in (0,M^2).
Write s_b=sqrt(M^2-z_b)>0 and choose its brane coordinate q_b=1. Its spatial
wavefunction and mass relation are

```text
Phi_b(y)=g exp(-s_b y)/(c+s_b),
m^2-z_b=g^2/(c+s_b),
||b||^2=1+g^2/[2s_b(c+s_b)^2].                          (18)
```

Its conserved Hilbert-space pairing with a continuum mode in (3) is zero:

```text
<b,psi_p>=Q+g/(c+s_b)[1/(s_b+ip)+R/(s_b-ip)]=0,         (19)
```

after inserting D=g^2/(c+s_b)-s_b^2-p^2. This can also be obtained from Green's
identity for the common self-adjoint operator. Equation (19) proves that (13)
has exactly zero bound-mode projection for both its position and velocity
data, not merely a small one. Spectral
projections commute with autonomous linear evolution, so no subsequent
continuum-to-bound transfer occurs.

If a bound component is initially present, it persists as an oscillation with
frequency sqrt(k^2+z_b), and the corresponding conclusion q->0 is false.
A packet placed far from the brane at a finite starting time is not automatically
orthogonal to the bound mode; both position and velocity projections must be
checked or the exact incoming scattering preparation (13) used. Setting only
q and q_dot to zero at that starting time does not remove a bulk overlap with
the bound tail. Similarly, a freely specified
initial q excitation can populate the bound sector. Neither case is capture
from an asymptotically pure continuum input. At a threshold resonance the
constants may deteriorate as the support approaches p=0; that limit is excluded
by p_0>0. If no bound mode exists, the bound-occupation condition is vacuous.

For a coherent quantum packet on the coupled ground state, the displacement
changes these classical mean fields and their excess energy, while connected
Gaussian two-point functions remain those of the background state. It produces
no new thermal noise spectrum or derived equilibrium temperature. Any bound
occupation above the ground-state baseline remains zero. These statements
concern the specified free Gaussian preparation, not arbitrary quantum states.

## 6. Scientific consequence and scope

This specified incoming bulk pulse can produce nonzero transient brane motion
and temporarily stored q energy. Its coherent energy is fully reflected into
the bulk continuum asymptotically, with no residual q excitation for the packet
class above. A narrow resonance or a long dwell time does not change the unit
reflection or establish permanent absorption. Finite observation times can see
transients, and the decay constants are not uniform in all parameter limits.

Thermodynamic capture, reheating and a hot Big Bang do not follow in this
autonomous free one-channel testbed. Additional accessible brane channels,
interactions, time-dependent backgrounds/couplings or cosmological dynamics
would define different models requiring a fresh energy/state analysis. This
result does not exclude such models or HDBLAST in general. The initial bulk
packet and its energy are declared inputs, not a derived primordial source.
The accepted four-dimensional continuum-bath representation also reproduces
this scattering preparation after the same spectral transformation; this new
channel supplies no independent dimensional-origin label.

All controls here are manufactured algebra and hypothesis checks. No original
retained arrays, physical source or trajectory callbacks, observed data,
likelihoods or network sources are used. Historical metric calibration stays
**FAIL**; the full continuous certificate stays **UNRESOLVED**;
higher-dimensional Big Bang cause stays **NOT_ESTABLISHED**; external novelty
stays **NOT_ASSESSED**.
