# Independent review of the incident-continuum extension

8 October 2026. **ACCEPT_CONDITIONAL_ANALYTIC_INCIDENT_EXTENSION** for the six
source files pinned by incident manifest SHA-256
`e7183b5c1b3e9a5c56e5a1542bfaec533c6f228699c6af805cf69113fbe13658`.
The final theorem SHA-256 is
`2914c9f287f60398bbe36846d30f04a6c7743aab535750534d13b3697b3ae6e5`.
No analytic blocker remains within the declared domain and preparation. This
is not physical-source generation, thermalization, a Big Bang result, a
general exclusion of HDBLAST, or authorization for protected numerical work.

The exact reviewed bytes are preserved under `reviewed-source/` and listed in
`REVIEWED_SOURCE_HASHES.json`. The accepted original model's six files and its
manifest `d82d09d1786b1317dbd5fb8852e618a4ecd63b26bb70b729e9be4f041a5bc15c`
were independently rehashed unchanged. Later changes to the incident source
require a fresh comparison and receipt.

## Scattering amplitudes, current and response

The uneliminated boundary system is

```
[(ip-c)  g] [R] = [ip+c]
[  -g    D] [Q]   [  g ] .
```

Its determinant is F=(ip-c)D+g^2. Solving the matrix gives the stated R and Q
without dividing by D. For real p>0 and g nonzero,
`|F|^2=(g^2-cD)^2+p^2D^2` cannot vanish: its imaginary component forces D=0,
then its real component is g^2. The numerator of R is -F*, which proves unit
reflection. Reflection is even in g and Q is odd, as required by field-sign
conventions.

At D=0 the matrix is nonsingular, R=-1, the bulk trace vanishes and Q=2ip/g.
The nonzero bulk derivative still matches -gQ. The source correctly describes
this as the bare-brane frequency, not a universal location of a response
maximum. At exactly g=D=0 the matrix has rank one and leaves the decoupled
brane oscillator undetermined. The g->0 limit at fixed nonzero D is regular;
the on-resonance harmonic limit is different and cannot define a packet-energy
divergence.

Independently evaluating -Phi_dot Phi_y for an arbitrary complex reflection
coefficient gives `omega*p*(|R|^2-1)/2`. The incoming/outgoing cross term has
zero real contribution to this signed current. The result is zero for the
accepted solution. A separate negative control changes R's phase while
retaining unit modulus: the boundary equation then fails. Unit modulus alone
therefore would not suffice to validate the matching calculation.

The uncoupled Robin field has boundary value 1+R0=2ip/(ip-c). Consequently
Q=chi_R*g*(1+R0) and R=R0+g*G_R*Q with precisely the earlier retarded kernel.
The free incoming trace is not one. Radiation from q interferes with direct
Robin reflection; its positive outward flux is not irreversible loss of the
incident energy. Positive temporary bare brane energy is consistent with zero
stationary mean work and unit total reflected flux.

The total energy includes the Robin and interaction terms with the accepted
signs. No additional dissipation, gravitational energy source or channel has
been introduced. As a useful independent consequence of the accepted
coercivity inequality, `E_q(t)<=E_total/(1-r)` with
`r=|g|/[m*sqrt(c+M)]<1`. Thus a fixed finite incoming energy bounds temporary
brane storage. This follows from the original full energy theorem; no change
to the frozen incident source is needed.

## Packet decay and spectral preparation

For a smooth compact a(k,p), common p support bounded away from zero, and
bounded k support, the group velocity p/omega has a strictly positive uniform
lower bound. F is nonzero on that compact support, so R,Q and all required
derivatives are bounded. Full-space Fourier extension followed by restriction
to y>=0 proves finite energy at each finite time. Fixed k alone proves fiber
energy, not global energy on the three-dimensional brane.

Repeated integration by parts with `L f=partial_p(f/partial_p omega)` gives
the stated all-N inverse-time bound without endpoint terms. Its constants are
uniform for y in a fixed finite slab and k in the compact support. Plancherel
in k then controls integrated brane energy and the finite-layer energy
majorant, including the absolute interaction term. The estimate does not
assert uniformity as p support reaches zero or g reaches zero.

The incoming component occupies the positive half-line asymptotically at
early time; the outgoing component occupies it at late time. On the opposite
half-line the phase derivative is bounded below by |y|+v_min|t|. Repeated
nonstationary integration bounds its energy norm by an integrable tail, and
Cauchy-Schwarz removes the cross-energy contribution. The source's
1/sqrt(2pi) normalization then gives
`E_in=int omega^2|a|^2 dp=E_out=int omega^2|aR|^2 dp` in its declared complex
fiber norm. That complex energy is the sum of real and imaginary solution
energies. Global bulk energy is conserved while energy leaves every fixed
layer; local decay is not disappearance of total energy.

The full bound-continuum inner product contains both q and the bulk integral:

```
Q + g/(c+s) [1/(s+ip)+R/(s-ip)] = 0,
D=g^2/(c+s)-s^2-p^2.
```

This identity holds exactly, and multiplying each mode by -i omega proves the
velocity projection is also zero. Omitting the q term gives a nonzero overlap.
The spectral projectors commute with the autonomous wave evolution, so zero
bound position and velocity remain zero. The result is no continuum-to-bound
transfer for this preparation, not dynamical erasure of an arbitrary bound
component.

For an explicit counterexample to a weaker initial-data premise, use the
manufactured M=2,c=3,g=1,m^2=13/4 model and initial q=q_dot=0,
Phi=(1+4y)exp(-y), Phi_dot=0. The derivative boundary condition is satisfied.
The normalized bound vector has q component sqrt(32/33), and the initial
overlap is `3*sqrt(32/33)/8`, not zero. Its persistent q contribution is
`(4/11)cos(sqrt(3)t)`. Such a component is present initially in the spectral
sense; it is not capture from a pure incoming continuum state. The final
source explicitly excludes this inference.

The delegated `packet-audit/PACKET_ANALYTIC_REVIEW.md` independently supplies
the finite/global energy proof, local decay argument, full asymptotic energy
normalization and bound orthogonality. It reread the final theorem hash above
and found no blocker.

## Exact controls and their scope

The independent script imports no producer code or project data. Its 35 checks
passed in normal and optimized Python with byte-identical receipts: 29 exact
symbolic checks, four negative controls, one exact matrix-rank check and one
finite endpoint-expansion bound. Failure checks remain active under `python -O`.

One additional manufactured packet makes the decay mechanism concrete without
sampling a trajectory. Choose a frequency window
`B(omega)=(omega-3)^3(4-omega)^3` on [3,4], zero elsewhere, and
`a(p)=sqrt(2pi)*(p/omega)*B(omega)/Q(p)` in the manufactured bound-mode model.
Then the complex q packet equals `int_3^4 B(omega)exp(-i omega t)domega`.
It is nonzero at t=0, where q=1/140. Exact endpoint integration gives
`|q|<=2316|t|^-4`, `|q_dot|<=18186|t|^-4` and a finite bare-energy bound
`174081564|t|^-8` for |t|>=1. These are loose exact bounds in manufactured
units, not fitted physical predictions. This window is C2 after zero
extension, so the control establishes a finite-order decay result; it is not
substituted for the theorem's C-infinity, all-N argument. No quadrature or
interval certificate is claimed.

The frozen producer script was inspected and replayed in both modes. Its 25
exact controls and eight negative controls passed, and both independently
replayed receipts match the frozen producer receipts byte for byte. There are
therefore 68 manufactured checks per mode, plus independent analytic review.

## Bounded scientific conclusion

This continuum preparation produces transient q motion and local storage,
followed by complete reflected continuum energy and no residual coherent
brane excitation in the stated packet class. In a coherent displacement of
the coupled ground state, connected covariance remains unchanged; vacuum
fluctuations are not residual coherent capture or newly generated thermal
noise. No observed temperature or thermalization mechanism follows.

The incoming packet and its energy are declared inputs. Extra accessible
channels, nonlinear interactions, time-dependent couplings/backgrounds or
gravity define different models and require a new analysis. The accepted
four-dimensional continuum representation reproduces the matched preparation,
so this extension does not identify a geometric extra dimension. These results
do not exclude HDBLAST generally or establish a primordial source.

No retained array, observed datum, protected physical source/target or
likelihood was evaluated; no network call or remote/publication mutation was
made. Historical metric calibration remains **FAIL**, the full continuous
certificate **UNRESOLVED**, higher-dimensional cause **NOT_ESTABLISHED**, and
external novelty **NOT_ASSESSED**.
