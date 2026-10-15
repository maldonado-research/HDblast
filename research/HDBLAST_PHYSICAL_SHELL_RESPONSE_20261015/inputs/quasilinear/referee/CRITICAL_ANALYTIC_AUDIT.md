# Critical analytic audit of the exact-contact program

Internal independent mathematical audit, 15 October 2026 UTC. This preliminary
note audits the established reduction and records requirements for a new proof;
it is not an approval of a yet-unwritten quasilinear existence theorem.

Reviewed inputs:

- `NONLINEAR_BOUNDARY_STRAIGHTENING.md`, SHA256
  `da221f388e0397e4a3b06ce40be1feff1d170311f3b7eaf79f3644b667dd6882`.
- `ANGULAR_UNIFORM_LINEAR_CONTACT.md`, SHA256
  `cc559ed4a16c920310fca4c97c3af869508414b6d0ffc9fd8b4612e82989dff1`.
- `NEXT_QUASILINEAR_CONTACT_TARGET.md`, SHA256
  `52918865bb4e63cd3ff1bdb92ddd70ee4006e38f34092fb93772e3ec0e6a330c`.

## 1. The elementary variable-coefficient energy is valid

Write the principal operator as

\[
 a^{00}u_{tt}+2a^{0i}D_i u_t+a^{ij}D_iD_j u,
 \qquad A=-a^{00}>0,\quad (a^{ij})>0,
\]

with a fixed positive fiber metric and a split-preserving fixed connection.
On a fixed compact time/collar domain the near-background assumptions supply
uniform lower bounds for both displayed positive forms. Define

\[
 E=\frac12\int[A|u_t|^2+a^{ij}\langle D_i u,D_j u\rangle+|u|^2].
\]

An integration by parts gives the principal boundary flux

\[
 \nu_i a^{0i}|u_t|^2+\nu_j a^{ij}\langle D_i u,u_t\rangle.
\]

At the shell, all normal–tangential inverse-metric entries vanish, including
`a^{zt}`. Thus the flux is `-a^{zz}<u_z,u_t>` for the inward coordinate.
The orthogonal four-dimensional Dirichlet and twelve-dimensional Neumann
blocks make it vanish. Every remaining term is bounded by coefficient first
spacetime derivatives times `E`, plus `sqrt(E)||F||_2`. This derivation does
not require selfadjoint lower-order matrices or a positive physical energy.
It does require the positive *spatial inverse block*: timelikeness of `dt`
alone would not justify this particular multiplier outside the chosen patch.

Tangential/time derivatives preserve the boundary conditions and differentiate
an identically zero `a^{za}|_b`. Hence there is no hidden principal boundary
forcing when these derivatives commute with the equation. Globally one may
use the sphere Killing action on the full tensor bundle; componentwise scalar
harmonic truncation is not a substitute for this invariant split.

## 2. What a high-order closure actually needs

For `m>=6` in four spatial dimensions, let

\[
 J_m(u)=\sum_{j=0}^{m}\|\partial_t^j u\|_{H^{m-j}}.
\]

Time/angular energy estimates through tangential order `m-1` control all
terms having at most one normal derivative. Recover a term with `r>=2`
normal derivatives by differentiating

\[
 a^{zz}u_{zz}=-a^{tt}u_{tt}-2a^{zt}u_{zt}
 -2a^{zA}u_{zA}-2a^{tA}u_{tA}-a^{AB}u_{AB}+F
\]

only `r-2` times normally. The mixed terms have at most `r-1` normal
unknown derivatives; all coefficient commutators and differentiated
first-derivative nonlinearities also have strictly fewer unknown normal
orders. This is a viable induction because `a^{zz}` stays positive.

At the top tangential order, a term `D^alpha a(v) D^2u` must be estimated
by a tame product bound, for example coefficient high Sobolev norm times
`D^2u` in `L^infinity`, plus coefficient low `W^{1,infinity}` norms times
`J_m(u)`. Bounding all coefficient derivatives in `W^{m,infinity}` would
not close at the advertised finite data regularity. One must also retain
the pointwise differentiated source norms needed in normal recovery; an
`L^1_t` source energy bound alone does not provide them unless accompanied
by a time derivative/initial-value bound.

For a frozen system, the weak variational space is the common mixed
`H_D^1` space. The natural Neumann conormal equals `a^{zz}u_z` because all
`a^{za}` vanish at the shell. This gives a plausible Galerkin existence
construction in addition to the energy estimate; an estimate alone does
not establish existence. Arbitrary angular projections need not preserve
higher corner compatibility for a variable projected operator. The actual
selected data avoid this trap by vanishing on an open initial shell collar,
provided frozen forcing is flat on a common initial corner wedge.

## 3. An iteration trap and a concrete repair

Homogeneous boundary traces alone do not make Picard trials compatible.
There is an additional global initial-jet issue: high energies contain
`partial_t^j u(0)` for `j>=2`. If these depend arbitrarily on the previous
trial, a uniform high-energy estimate can inherit a term without a small
iteration-time factor. Restrict trials to the same exact nonlinear
initial time jets through the chosen finite order. The frozen equation
then reproduces those jets recursively.

A naive time Taylor polynomial `sum t^j J_j/j!` is not generally an
`C_t H^m` extension because `J_j` has only `H^{m-j}` spatial regularity.
A concrete repair is to use the already available short, protected local
reduced Cauchy evolution with the desired lapse/shift, then multiply it by
a smooth time cutoff after a common short interval. Its regularity is the
required `J_m` class, it realizes the true jets, and finite propagation
makes it zero in a common initial shell-corner wedge. Frozen iterates with
these jets and flat corner forcing retain all required corner identities.

For contraction, the difference forcing contains
`[a(v)-a(w)]D^2u`. With high solutions uniformly in `J_m`, its natural
bound is in the source norm for `J_{m-1}`, since `D^2u` is only `H^{m-2}`.
Contraction in `J_m` without extra data derivatives is therefore unjustified.
A lower-order contraction is consistent with the stated target.

## 4. Gauge, parameter, lifespan and scope checks

The old exact-precontact theorem chooses initial lapse and shift equal to
the background. That gives the prescribed physical canonical tangent but
need not give the particular longitudinal or harmonic spacetime metric
first tangent. If full tangent identity is claimed, prescribe initial
lapse/shift first parameter derivatives from the selected harmonic tangent;
then solve the exact harmonic conditions for their time velocities. The
choice remains background on the open shell collar.

For the exact subsidiary vector, the invariant tangential boundary row is
`(nabla_n C)_a-K_a^b C_b=0`, with `C_n=0`. In proper-y coordinate components
this is

\[
 \partial_y C_a-\frac{N\sigma}{3}C_a=0.
\]

The factor two relative to `K=σh/6` comes from the connection term in
`nabla_n C_a`. A missing connection term would give the wrong Robin factor.
This is again a scalar-principal variable-coefficient mixed problem, with
lower-order terms controlled by the actual sufficiently regular metric.

The natural formal variation values lose one Sobolev order for the first
variation and two for the second. Difference-quotient convergence needs
additional room. Reserving four total orders is conservative and plausible;
it must be proved in the through-boundary norm, not imported merely from
the old precontact Cauchy statement. No top-norm smooth solution map follows.

A local time from zero does not suffice. A perturbation energy inequality
on the *fixed background interval*, with constants controlled under a small
state/high-norm bootstrap, can give `J_m(V)<=C_T J_m(V(0))` and hence any
prescribed finite background interval for sufficiently small initial data.
The parameter radius can depend on that interval. Cone support must use a
uniform bound for the actual metric's coordinate speeds, for example
`1+eta` in z after sufficiently small near-background perturbation, rather
than the exactly background value one.

Finally, existence on an interval extending beyond the reference contact
window does not by itself prove a nonzero gauge-invariant shell response,
or identify an exact first-arrival time. Such an observable nonvanishing
claim requires its own calculation. Finite regularity and finite time do
not imply spatial C-infinity, global stability, novelty, or cosmological
origin evidence.
