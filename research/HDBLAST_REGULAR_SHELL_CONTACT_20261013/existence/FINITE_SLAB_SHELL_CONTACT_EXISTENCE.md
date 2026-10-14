# Finite-slab regular-master solutions through actual shell contact

13 October 2026 (America/Los_Angeles). Original proof and internal check, confined to the research
workspace. No repository is modified. The primary theorem below concerns a
scalar master equation. Its conversion to a coupled shell perturbation needs
the independent junction-sufficiency calculation; it is not inferred from a
general scalar Wentzell theorem.

## 1. Precise result and actual background domain

Write `b=r_b>0`, `s=R²−T²`, and assume the **actual** registered background
extends smoothly on `−s_-<s≤b²`, with `Omega>0` and `f!=0`. Analyticity is one
way the inherited work establishes that extension; arbitrary continuation of
coefficient functions is not substituted for the background equations. On
this domain the canonical master obeys

    (Box_eta−V(s)) psi=0,
    V=3L−2s L_s+s L²+(8/3)s f².

The shell `R²−T²=b²` has outward Minkowski unit normal `n_eta`. Fix either

    n_eta psi−c Box_Sigma psi+d psi=0,       c>0,       (even)
    n_eta psi+d psi=0.                                (odd)

Here `c=(Omega_b B_b)^−1`, and `d` is the appropriate shifted real Robin
coefficient from the inherited derivation. No positivity of `d` or `V` is
needed.

**Constructive theorem.** Choose

    0<tau_*<min{pi/4, 2 arctan(sqrt(s_-)/b)},
    0<a<d_0<tau_*/2.

For each spherical harmonic `H_l` on `S³`, with `lambda=l(l+2)`, take arbitrary
smooth radial initial profiles supported in the annulus described below by
`a≤x≤d_0`. Either boundary law has a unique smooth solution on the corresponding
closed time slab, smooth through and after its first causal shell contact.
Finite sums of such harmonics are allowed. The radial profiles are arbitrary
smooth functions, so this is an infinite-dimensional class. In particular
`l≥2` can be fixed throughout.

For time-symmetric data with a nonzero displacement profile whose support has
infimum depth `a`, the shell trace vanishes before contact and is nonzero at
some time in **every** interval immediately after contact. Thus the theorem
gives actual boundary interaction, not only a causal upper bound on support.

The proof establishes neither a completed Sobolev phase space for the coupled
gravitational problem nor a theorem covering every infinite angular expansion.
The energy estimate in section 4 applies more broadly to smooth scalar
solutions, but that broader estimate is not used as an unproved smoothness
theorem.

## 2. Two explicit fixed-domain coordinate systems

### 2.1 Ambient-time fixed ball

Put `A(T)=sqrt(b²+T²)`, `y=X/A(T)`, `h=T/A²`, and
`u(T,y)=psi(T,A(T)y)`. The shell becomes `|y|=1`, and

    D=partial_T−h y·grad_y,
    D²u−A^−2 Delta_y u+V(A²|y|²−T²)u=0.

The principal quadratic form controlling spatial derivatives is

    A² delta_ij−T² y_i y_j,

whose radial minimum on the unit ball is `b²>0`. In four spatial dimensions
the principal Hamiltonian density is

    (1/2)[A^4 u_T²+(A² delta_ij−T² y_i y_j)u_i u_j].

The boundary normal and wave operator become

    n_eta psi=(T/b)u_T+(b/A²)u_r,
    Box_Sigma z=−(A²/b²)z_TT−(4T/b²)z_T+A^−2 Delta_S3 z.

These formulas provide a direct fixed-ball formulation on every finite
ambient slab within the actual background domain. The following conformal
coordinates remove the mixed time-space derivatives and make the constructive
proof shorter.

### 2.2 Fixed round hemisphere

Use `0≤chi≤pi/2`, `|tau|<pi/2`, and

    q=b/(cos(tau)+cos(chi)),
    T=q sin(tau),       R=q sin(chi).

Then exactly

    eta=q²[−d tau²+d chi²+sin²(chi)dOmega_3²],
    s=b²[cos(tau)−cos(chi)]/[cos(tau)+cos(chi)].

The shell is `chi=pi/2`; write `q_b=b sec(tau)`. The compact rectangle
`|tau|≤tau_*`, `0≤chi≤pi/2` stays in the actual background because

    −b² tan²(tau_*/2) ≤ s ≤ b².

This is the relevant domain bound. The shell can have a large ambient time
while remaining at `s=b²`; bounding the whole rectangle by minus the square
of the largest shell time is unnecessarily restrictive.

Set `psi=q^−3/2 u` on the hemisphere `M=S⁴_+`. The bulk equation is

    u_tautau−Delta_M u+U(tau,chi)u=0,
    U=9/4+q²V(s).                                           (2.1)

The curvature term `9/4` follows either by direct differentiation or from
the five-dimensional conformal wave operator on the Einstein cylinder,
whose scalar curvature is 12 and conformal coupling is 3/16.

The shell law becomes

    partial_nu u+partial_tau(mu u_tau)−mu Delta_Gamma u+beta u=0, (2.2)
    mu=c cos(tau)/b>0,
    beta=(d−3/(2b))q_b
         −mu[(3/2)sec²(tau)+(3/4)tan²(tau)],

or, for odd parity,

    partial_nu u+alpha u=0,
    alpha=q_b[d−3/(2b)].                                    (2.3)

Here `Gamma=S³`, and `partial_nu=partial_chi`. In particular the exact
coefficient of `u_tau` in (2.2) is `mu'`, not an independently chosen damping.
The identity used in this conversion is

    Box_Sigma(q_b^−3/2 z)=q_b^−7/2[
      −z_tautau+tan(tau)z_tau+Delta_Gamma z
      +((3/2)sec²(tau)+(3/4)tan²(tau))z].

All `U,mu,beta,alpha` are even in `tau`.

The rectangle contains **every point** of the full moving-ball ambient slab

    |T|≤T_0,       R≤sqrt(b²+T²),
    T_0=b tan(tau_*/2).

Indeed the inverse time coordinate is
`tau=arctan((T+R)/b)+arctan((T−R)/b)` and has, for positive `T`, its maximum
`2 arctan(T/b)` at the center. The rectangle does not cover arbitrary large
ambient slabs. This coverage restriction is built into the theorem.

## 3. Corner compatibility and admissible data

For smooth general hemisphere data let `u_0=u|tau=0`, `u_1=u_tau|tau=0`.
Define successive formal bulk time jets recursively by

    u_(j+2)=Delta_M u_j−sum_(k=0)^j binom(j,k) U^(k)(0)u_(j−k).

The full even compatibility hierarchy is, for every `j≥0`,

    partial_nu u_j
      +sum_(k=0)^j binom(j,k){
         mu^(k)(0)u_(j−k+2)+mu^(k+1)(0)u_(j−k+1)
         −mu^(k)(0)Delta_Gamma u_(j−k)+beta^(k)(0)u_(j−k)}=0

on `Gamma`. The odd hierarchy is

    partial_nu u_j+sum_(k=0)^j binom(j,k)alpha^(k)(0)u_(j−k)=0.

For classical solutions the boundary displacement and velocity are the traces
of `u_0,u_1`. They are not extra freely specifiable fields. An auxiliary energy
completion may carry separate `L²` boundary components, but that observation
does not remove the trace compatibility required for smooth data.

Now put `x=pi/2−chi`, and choose `u_0=h_0(x)H_l`,
`u_1=h_1(x)H_l` with both profiles supported in `[a,d_0]`. They vanish in an
open shell collar, so every formal jet also vanishes there and every order of
both compatibility hierarchies holds. They also vanish near the origin, so
the harmonic notation introduces no coordinate singularity in the data.

At `tau=0`, `R=b tan(chi/2)` and `q=b/(1+cos(chi))`; the profiles therefore
give smooth compact ambient Cauchy data. More explicitly

    psi|T=0=q^−3/2 u_0,
    psi_T|T=0=q^−5/2 u_1,

because `q_tau=R_tau=0` and `T_tau=q` on this slice. In particular `u_1=0`
is exactly time-symmetric ambient canonical data.

For a completely concrete nonzero choice, set
`tau_*=min{pi/8,arctan(sqrt(s_-)/b)}`, `a=tau_*/4`, `d_0=tau_*/3`, and

    h_0(x)=exp[−1/((x−a)(d_0−x))] for a<x<d_0, and zero otherwise,
    h_1=0.

Multiply `h_0` by any nonzero real amplitude and choose any real `H_l` with
`l≥2`. These choices satisfy every strict inequality used in the proof,
including the optional positive-`s` support bound in section 6. The formula
is an existence choice; its small numerical magnitude is irrelevant.

## 4. Positive auxiliary control energy

For real functions the even weak equation has mass and stiffness forms

    m_tau(v,w)=int_M vw+int_Gamma mu vw,
    a_tau(v,w)=int_M(grad v·grad w+Uvw)
                  +int_Gamma(mu grad_Gamma v·grad_Gamma w+beta vw),

and reads `d_tau m_tau(u_tau,w)+a_tau(u,w)=0` for time-independent tests.
The natural form domain is `H¹(M)` with trace in `H¹(Gamma)`; the trace is
identified with the boundary component. For odd parity use

    m_tau(v,w)=int_M vw,
    a_tau(v,w)=int_M(grad v·grad w+Uvw)+int_Gamma alpha vw,

on `H¹(M)`.

On the compact rectangle `mu_min>0` and all coefficient derivatives are
bounded. In even parity adding `K m_tau(u,u)` for sufficiently large constant
`K` makes `a_tau+K m_tau` uniformly coercive. In odd parity use the trace
inequality `||u||²_L²(Gamma)≤epsilon||grad u||²+C_epsilon||u||²` first and
then choose `K`. Thus

    E_K=(1/2)[m_tau(u_tau,u_tau)+a_tau(u,u)+K m_tau(u,u)]>0

controls precisely the corresponding bulk and boundary first-derivative
norms. Differentiating and using the common boundary flux gives the exact
identity

    E_K'=−(1/2)m_tau'(u_tau,u_tau)+(1/2)a_tau'(u,u)
            +K m_tau(u,u_tau)+(K/2)m_tau'(u,u).

Consequently `|E_K'|≤C E_K` and
`E_K(tau)≤exp(C|tau−tau_0|) E_K(tau_0)`. The same estimate with source terms
gives scalar uniqueness and continuous dependence on each such slab.
The constant depends on the compact slab and coefficients; no uniform
infinite-time estimate is claimed.

This is an **auxiliary positive scalar control norm**. It is not the
Einstein–scalar action, a gravitational Hamiltonian, a positivity theorem
for physical shell kinetic energy, or a ghost criterion.

The forms permit the usual finite-dimensional Galerkin approximation:
choose a dense smooth form-domain basis, solve its positive-mass linear ODE
system, apply the displayed bound uniformly, and pass weakly on compact
time intervals. A general theorem by this route would also supply the
weak-form compactness, time regularization, and compatibility-preserving
approximation details. The displayed energy estimate alone does not justify
all orders of classical regularity. The next section supplies an explicit smooth
construction for the theorem's class and avoids that logical gap.

## 5. Smooth construction by characteristics and a boundary ODE

Fix one `H_l`; rotational invariance preserves its angular factor. Write
`u=u_l(tau,x)H_l` and conjugate the radial operator by

    v=cos(x)^(3/2)u_l.

On the actual collar the equation is exactly

    v_tautau−v_xx+W_l(tau,x)v=0,
    W_l=q²V(s)+(lambda+3/4)sec²(x).                         (5.1)

At `x=0`, the derivative of `cos(x)^(3/2)` vanishes. The boundary law becomes

    −v_x+mu v_tautau+mu'v_tau+(mu lambda+beta)v=0,            (5.2)

or `−v_x+alpha v=0` in odd parity. The support estimate below gives
`x≤d_0+|tau|≤d_0+tau_*<pi/2`. Choose a number between this upper bound and
`pi/2`, and extend `W_l` smoothly and boundedly beyond that collar solely
to solve a half-line problem. No resulting physical claim uses the
extension, which is causally inaccessible to the solution.

Let

    p=v_tau+v_x,      m=v_tau−v_x,
    z(tau)=v(tau,0),  w=z'.

The outgoing and incoming characteristic equations are

    p_tau−p_x=−W_l v,       m_tau+m_x=−W_l v.                 (5.3)

In even parity the boundary relation is the positive-mass ODE

    mu w'+(mu'+1)w+(mu lambda+beta)z=p(tau,0),
    z'=w,                    m(tau,0)=2w−p(tau,0).           (5.4)

In odd parity it is

    z'=p(tau,0)−alpha z,
    m(tau,0)=p(tau,0)−2alpha z.                             (5.5)

These relations show directly which characteristic is supplied by the
boundary. There is no prescription of two independent incoming wave data.

For `tau=t≥0` the Volterra representation is

    p(t,x)=p_0(x+t)−int_0^t W_l(s,x+t−s)v(s,x+t−s)ds,

    m(t,x)=m_0(x−t)−int_0^t W_l(s,x−t+s)v(s,x−t+s)ds, x≥t,

    m(t,x)=m(t−x,0)
              −int_(t−x)^t W_l(s,x−t+s)v(s,x−t+s)ds, x<t,

    v(t,x)=v_0(x)+(1/2)int_0^t[p(s,x)+m(s,x)]ds,

together with the integrated ODE (5.4) or (5.5).

Here is an explicit contraction argument, rather than a claim that the
characteristic formulas automatically solve the problem. Substitute the
formula for `p` into those for `m` and the boundary ODE, and regard `(v,z,w)`
as unknown in even parity. On an interval of length `delta≤1`, differences
of the `p` map are bounded by `delta ||W_l||_infty ||Delta v||`. Differences
of `m` are bounded by `2||Delta w||+2delta||W_l||_infty||Delta v||`.
After the final time integrals, differences of the new `v,z,w` are bounded
by a constant times `delta` times their joint uniform norm. The constant is
finite and can be bounded using

    1+||W_l||_infty
      +mu_min^−1[1+||mu'||_infty+lambda||mu||_infty+||beta||_infty].

Taking `delta` less than the reciprocal of ten times the sum of these
positive bounds is sufficient after using the maximum norm on the three
components. The odd map uses `(v,z)` and the analogous finite bound
`1+||W_l||_infty+||alpha||_infty`.

Banach iteration therefore gives a unique continuous solution on a short
interval. Differentiation of the characteristic integrals and boundary
ODE gives all derivatives for smooth data. The only possible matching
problem is the characteristic leaving the initial boundary corner; its
matching conditions are exactly section 3's hierarchy. Our profiles vanish
on a whole neighborhood of that corner, so the iterates and their limit
are identically zero there, satisfying every matching condition. One may
equivalently prove `C^k` regularity successively on compact characteristic
triangles; the equations for the highest derivatives are linear Volterra
equations, with the already constructed lower derivatives as smooth sources.

The auxiliary relation `v_x=(p−m)/2` follows from its initial validity and
the zero time derivative of its defect, using (5.3). The trace identity
`v(t,0)=z(t)` follows from its initial validity and `v_tau(t,0)=w` in (5.4),
or `v_tau(t,0)=p−alpha z` in (5.5). Thus no constraint was lost by the
first-order construction.

The contraction interval depends only on bounded coefficients, `mu_min`,
and fixed `lambda`, not the amplitude of the solution. Iterating finitely
many intervals covers `0≤tau≤tau_*`; the linear Volterra bounds or the
positive energy bound exclude a finite-slab blowup. The same construction
backward in time covers `−tau_*≤tau≤0`. Restarted data inherit compatibility
from the smooth solution. Finite sums of modes are then obtained by linear
superposition. Uniform estimates as `l→infty` are not required or claimed
for this constructive theorem.

## 6. Finite propagation, genuine activation, and ambient contact time

The characteristic formulas and uniqueness imply

    supp v(tau,.) subset [max(0,a−|tau|),d_0+|tau|].

The homogeneous boundary ODE stays zero until the outgoing characteristic
can bring a nonzero signal to it. Thus `z(tau)=0` when `|tau|<a`. The same
observation proves that the arbitrary coefficient extension used above has
no effect anywhere on the physical solution's support. This is the radial
finite-propagation statement for the constructed finite-harmonic family.
No general localized energy theorem for arbitrary angular data, particularly
with a possibly negative Robin potential, is inferred from the global
control estimate in section 4.

To show that contact really occurs, choose `h_1=0`, `h_0!=0`, with
`inf supp h_0=a`. All coefficients are time even, so uniqueness gives
`v(−tau,x)=v(tau,x)`. Suppose, for a number `tau_1>a`, that the shell trace
vanished throughout `0<tau<tau_1`. It would vanish on `−tau_1<tau<tau_1`,
and either boundary equation would imply `v_x(tau,0)=0` on that interval.

Now view (5.1) as a one-spatial-dimensional wave equation with `x` as time:

    v_xx−v_tautau−W_l(tau,x)v=0.

Zero displacement and normal derivative at `x=0` have the unique zero
solution in `x+|tau|<tau_1`, by the ordinary one-dimensional Cauchy energy
estimate, equivalently the characteristic Volterra formula. At `tau=0`
this would force `h_0(x)=0` for `0<x<tau_1`, contradicting the infimum of
its support. Therefore the boundary trace is nonzero at some point of
every `(a,a+epsilon)` contained in the constructed slab. The trace is
smooth and may be flat at the first-contact time itself; no nonzero jet
at exactly the support edge is asserted.

At the shell `T=b tan(tau)`. The first causal contact is consequently

    T_c=b tan(a).

The outer edge of the ambient initial annulus is
`R_0=b tan(pi/4−a/2)`, and the same value follows directly from
`R_0+T_c=sqrt(b²+T_c²)`:

    T_c=(b²−R_0²)/(2R_0).

Since `a<tau_*/2`, `T_c<T_0=b tan(tau_*/2)`. We have therefore constructed
actual smooth solutions on a full ambient moving-ball slab containing a
nonempty time interval after contact. These shell-collar data are a new
family. This proof does not claim that the earlier near-cone flux pair
reaches the shell before leaving its known actual-background domain.

One may impose the additional bound `d_0+2tau_*<pi/2`; for example choosing
`tau_*<pi/5` together with `d_0<tau_*/2` suffices. The support estimate then
implies `chi>|tau|` everywhere the solution can be nonzero, hence `s>0` on
its full causal support. These interacting packets can therefore remain
entirely within the original positive-`s` background collar. The actual
negative-`s` extension is still the stated basis for speaking about the
whole hemisphere rectangle and the zero perturbation around its center.

### 6.1 Time-odd activation and a displacement–velocity master pair

There is a second, equally bounded activation statement. Choose instead
`h_0=0`, `h_1!=0`, with `inf supp h_1=a` and the same compact collar and
fixed harmonic assumptions. Time reversal of the bulk and either boundary
equation preserves the system: `U,mu,beta,alpha` are even and `mu'` is odd.
Uniqueness now gives

    v(−tau,x)=−v(tau,x).

Suppose its trace vanished throughout `0<tau<tau_1` for some `tau_1>a`
inside the constructed slab. Oddness extends that vanishing to
`−tau_1<tau<tau_1`; either boundary equation again supplies `v_x=0` there.
The same sideways one-dimensional Cauchy uniqueness argument then forces
`v=0` in `x+|tau|<tau_1`. Differentiating at `tau=0` gives

    0=v_tau(0,x)=cos(x)^(3/2)h_1(x),       0<x<tau_1,

contradicting `inf supp h_1=a`. Thus the time-odd solution also has nonzero
shell trace somewhere in every interval immediately after `tau=a`, with
the same ambient contact infimum `T_c=b tan(a)`. This does not assert a
nonzero trace at exactly `a` or at every later instant.

For clarity the canonical transformation does not give the same numerical
hemisphere profile to displacement and velocity data. Define

    w(s)=(Omega(s)f(s)²)^(1/2)>0,
    Q_0(chi)=b/(1+cos(chi)),       R=b tan(chi/2).

Since `F=w psi`, the exact Cauchy-data conversion on `T=tau=0` is

    u_0=Q_0^(3/2)w(R²)^−1 F|T=0,
    u_1=Q_0^(5/2)w(R²)^−1 F_T|T=0.                         (6.1)

Indeed `w_tau=Q_tau=R_tau=0` and `T_tau=Q_0` on that slice; equivalently
`w_T=0` at fixed `R` when `T=0`. All weights in (6.1) are smooth and
strictly positive, so multiplication by them preserves support and its
infimum depth. Consequently, for any nonzero real smooth annular profile
`q_0(R)H_l` in this collar, **both** master Cauchy states

    (F,F_T)=(q_0 H_l,0),       (F,F_T)=(0,q_0 H_l)

undergo the proved boundary activation. Their hemisphere profiles differ
by one factor of `Q_0`, as required by (6.1). This is a statement about two
actual interacting solutions. Whether their physical Cauchy symplectic
pairing is nonzero, and how the total physical current is conserved with
a moving shell, requires the separate action-derived calculation. No
nonzero local shell-current claim follows from activation alone.

## 7. Reconstruction and remaining scope

On the physical domain return to

    F=(Omega f²)^(1/2)psi,
    h_AB=2F(s eta_AB−3X_A X_B),
    chi=−3(E+4)F/(2Omega²f).

All coefficients are smooth and `f!=0`, so the constructed canonical
solutions reconstruct smoothly through the cone wherever their domain
includes it. The inherited direct residual identity proves all bulk
linearized Einstein–scalar equations. At the shell this does **not by
itself** prove the complete coupled junction system. The independent
junction derivation must supply, and separately checks, the sufficient
embeddings `zeta_e=0` and `zeta_o=−P_o/(a_b rho_b²)`, with `P=r⁴F`.

The present theorem is scalar uniqueness within its specified compatible
solution class; it is not uniqueness modulo every gravitational gauge
transformation. It establishes no completeness across exceptional harmonic
or embedding sectors, no global or all-time boundary theory, no equality
of completed physical phase spaces, no nonlinear stability, and no sign
claim for the physical gravitational Hamiltonian.

## 8. Independent checks

`check_shell_contact_coordinates.py` checks the coordinate metric,
conformal bulk and boundary conjugations, radial characteristic reduction,
ambient contact-time identity, and strict negative controls. These are
algebraic differential-operator checks. The contraction, support,
time-reversal activation and actual-background-domain arguments above are
proofs in this note; no computer output is presented as a certified PDE
existence proof.
