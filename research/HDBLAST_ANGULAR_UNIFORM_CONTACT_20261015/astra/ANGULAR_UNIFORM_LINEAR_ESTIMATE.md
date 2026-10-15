# Angular-uniform control of the actual linearized harmonic shell problem

15 October 2026. Internal mathematical derivation. This is a linear theorem
on the fixed smooth positive-cone collar, not an assertion that the published
second coefficient belongs to an exact nonlinear family. A separate note
examines the nonlinear reduction and identifies the remaining data issue.

## 1. Result and its precise scope

The full nonsymmetric Robin matrix is not an obstruction to an angular-uniform
linear energy estimate. A bounded, invertible, block-diagonal change of
dependent variables turns it into homogeneous Neumann conditions. Because
all sixteen reduced equations have the same scalar wave principal symbol,
this change introduces only first- and zeroth-order bulk terms. A positive
component control norm then has exactly zero homogeneous boundary flux.

This proves estimates uniform in the complete angular spectrum, including
the monopole and all tensor, vector and scalar sectors. It does not require
a scalar-master truncation, a symmetrizer of the Robin matrix, a favorable
Robin sign, or a division by an angular eigenvalue. It is a control estimate;
it is not a positive gravitational Hamiltonian or a stability theorem.

The high-regularity estimate below uses explicit time derivatives of the
source. The elementary boundary lifting also has explicitly stated sufficient
regularity. No sharp glancing estimate for rough boundary data is claimed.

## 2. Actual principal operator and actual Robin matrix

Use the model and reference wave gauge of the established checkpoint:

\[
 E_{AB}=R_{AB}-\phi_A\phi_B-\tfrac23U(\phi)g_{AB},\qquad
 S=\Box_g\phi-U'(\phi),\qquad
 E^{\rm red}_{AB}=E_{AB}-\nabla_{(A}C_{B)}.
\]

At the background
\(g_0=dy^2+\rho^2(-dt^2+\cosh^2t\,\kappa)\),
the metric principal part is \(-\Box_0/2\), and the scalar principal part is
\(\Box_0\). Multiplying the metric and scalar rows by their respective
nonzero normalizations gives the same scalar principal operator.
In \(z=-\log(r/b)\), \(\partial_y=-\rho^{-1}\partial_z\),

\[
 \rho^2\Box_0=-\partial_t^2+\partial_z^2+
 \operatorname{sech}^2t\,\Delta_\kappa
 -3\tanh t\,\partial_t-3a\rho\,\partial_z,
 \qquad a=\rho_y/\rho.
 \tag{1}
\]

Tensor connection and scalar-metric couplings add at most first derivatives.
Thus on a fixed controlled collar the full system is

\[
 P U:=U_{tt}-U_{zz}-c(t)^2\Delta_\kappa U
       +A^0U_t+A^zU_z+A^iD_iU+BU=F,
 \quad c(t)=\operatorname{sech}t.
 \tag{2}
\]

There are four Dirichlet components \(U_D=k_{ya}\). Put the twelve other
components in \(U_N=(\ell,H_{ab},w)\), where
\(\ell=k_{yy}\), \(H_{ab}=k_{ab}\). Their actual homogeneous boundary
operator is

\[
 -U_{N,z}+M U_N=0,\qquad M=\rho_b\mathcal B,
 \tag{3}
\]

\[
 \mathcal B(\ell,H,w)=
 \left(
 4a_b\ell-2a_b\operatorname{tr}_0H-\tfrac{4\sigma'_0}{3}w,
 -2a_bH-a_b\ell h_0-\tfrac{\sigma'_0}{3}wh_0,
 \tfrac{\sigma'_0}{4}\ell+\tfrac{\sigma''_0}{2}w
 \right).
 \tag{4}
\]

The inhomogeneous second-order data are exactly the right sides of equation
(5) in `boundary/GH_BOUNDARY_SOURCE_AND_CONSTRAINT_COMPLETION.md`, multiplied
by \(\rho_b\). In particular \(\ell=2\nu_2+2\nu_1^2\) must be used;
dropping \(2\nu_1^2\) changes those sources. The present estimate changes
none of the previously derived source formulas or gauge transformations.

Equation (4) is generally nonsymmetric in a positive component norm. The
trace uses the Lorentzian background tangential metric; this does not make
the *control norm* indefinite.

One clean bundle implementation fixes the positive product metric
\(dy^2+dt^2+\kappa\), its induced positive metric on the metric/scalar
bundle, and a product connection preserving

\[
 E_D=dy\mathbin{\odot}T^*(\mathbb R_t\times S^3),\qquad
 E_N=\operatorname{span}(dy^2)\oplus
 \operatorname{Sym}^2T^*(\mathbb R_t\times S^3)\oplus\mathbb R.
\]

The two summands are orthogonal. The difference from the spacetime tensor
connection only changes the first- and zeroth-order coefficients of (2).
Alternatively use a global left-invariant orthonormal frame on
\(S^3\simeq SU(2)\), together with \(dy,dt\); the splitting then has fixed
component projectors. This avoids differentiating coordinate components
through singular angular charts. Normal derivatives in the principal energy
are the product derivatives, so no connection-generated boundary term is
silently omitted.

## 3. Boundary conjugation

Let \(\theta(z)\) be fixed and smooth, equal to \(z\) near zero and zero
outside a fixed collar. Extend the boundary endomorphism \(M(t,\Omega)\)
in the product trivialization, and set

\[
 T_D=I,
 \qquad T_N(z,t,\Omega)=\exp\bigl(\theta(z)M(t,\Omega)\bigr),
 \qquad U=TV.
 \tag{5}
\]

Then \(T|_{z=0}=I\), \(T_z|_{z=0}=0\oplus M\), and

\[
 (-U_{N,z}+MU_N)|_{z=0}=-V_{N,z}|_{z=0}.
 \tag{6}
\]

Consequently homogeneous boundary data become
\(V_D|_0=0\), \(V_{N,z}|_0=0\). The exponential is invertible regardless
of the sign, spectrum, diagonalizability or symmetry of \(M\). For example,
\(\|T^{\pm1}\|\le \exp(\|\theta\|_\infty\|M\|_\infty)\).
Derivatives of any fixed order are bounded by functions of the corresponding
derivative bounds of \(M\), the fixed collar and its cutoff. A replacement
\(I+\theta M\) would not supply general invertibility.

Since the principal coefficient is a scalar times the identity,

\[
 \widetilde P=T^{-1}PT
 =\partial_t^2-\partial_z^2-c^2\Delta_\kappa
   +\widetilde A^\mu D_\mu+\widetilde B.
 \tag{7}
\]

For instance, the new first-order coefficients include
\(2T^{-1}T_t\), \(-2T^{-1}T_z\), and
\(-2c^2T^{-1}D_iT\). No second-order matrix coefficient is introduced.
For the actual background, \(M\) and \(T\) are rotationally equivariant
bundle endomorphisms, including their trace/insertion blocks.

## 4. Basic energy and angular uniformity

For \(V_D=0\), \(V_{N,z}=0\), use

\[
 \mathcal E[V](t)=\frac12\int_{\Sigma}
 \bigl(|V_t|^2+|V_z|^2+c(t)^2|DV|^2+|V|^2\bigr)\,dz\,d\omega,
 \quad \Sigma=\mathbb R_+\times S^3.
 \tag{8}
\]

The added \(L^2\) term avoids a Poincare constant depending on radial support.
All norms in (8) use the fixed positive product fiber metric. On a finite
time interval \([0,T_*]\), \(c\ge\operatorname{sech}T_*>0\).

Integration by parts gives boundary flux
\(-\langle V_t,V_z\rangle|_{z=0}\). Its Dirichlet part vanishes because
\(V_{D,t}=0\); its Neumann part vanishes because \(V_{N,z}=0\). Cross terms
are absent because the chosen positive fiber metric makes the splitting
orthogonal. First-order couplings need not preserve the splitting: they are
estimated as bulk terms without integration by parts. Therefore

\[
 \frac{d}{dt}\mathcal E[V]
 \le C\mathcal E[V]+\|\widetilde F(t)\|_{L^2}
                            \sqrt{2\mathcal E[V]},
\]

\[
 \mathcal E[V](t)^{1/2}
 \le C e^{Ct}\left(\mathcal E[V](0)^{1/2}
                 +\int_0^t\|\widetilde F(s)\|_{L^2}\,ds\right).
 \tag{9}
\]

Here and below \(C\) depends on the controlled background, \(T_*\), a finite
set of coefficient bounds and the *fixed* conjugation collar; it does not
depend on any angular eigenvalue or angular cutoff. Equivalence under (5)
gives the same bound for \(U\) in
\(H^1(\Sigma)\times L^2(\Sigma)\).

In particular high angular frequency is part of the positive gradient term,
not a frequency-dependent zeroth-order coefficient estimated by its size.
This is precisely the estimate unavailable from the finite-mode
characteristic proof. At glancing, no division by the normal characteristic
root occurs. A zero of the unweighted Neumann boundary symbol does not
invalidate (9).

## 5. High Sobolev estimate with explicit sufficient forcing norms

For an integer \(r\ge0\), define

\[
 J_r(W;t)=\sum_{j=0}^r
  \|\partial_t^jW(t)\|_{H^{r-j}(\Sigma)},
 \qquad J_r=0\quad(r<0).
 \tag{10}
\]

**Theorem.** Fix \(m\ge1\). Assume smooth coefficients with bounded
derivatives through order \(m+2\) on a fixed controlled collar; smooth
compatible data may initially be used. For homogeneous Dirichlet/Robin
boundary conditions, the solution satisfies

\[
 \boxed{
 \sup_{0\le t\le T_*}J_m(U;t)
 \le C_m\left[
 \|U(0)\|_{H^m}+\|U_t(0)\|_{H^{m-1}}
 +J_{m-2}(F;0)+\int_0^{T_*}J_{m-1}(F;t)\,dt
 \right].}
 \tag{11}
\]

The bound is also valid by completion for limits of smooth compatible data
in the displayed norms. In particular when \(F=0\), there is no derivative
loss relative to \(H^m\times H^{m-1}\) Cauchy data. The forcing requirement
in (11) is sufficient; no claim that it is necessary or optimal is made.

Here is a proof that tracks the boundary and every normal derivative.
Let \(X_i\) be a fixed collection of sphere Killing fields spanning the
tangent spaces, or the global left-invariant fields of the bi-invariant
round sphere. Their iterates define the standard Sobolev norms, and
\([X_i,\Delta_\kappa]=0\) on component functions. In the bundle formulation
one can instead use the induced rotation Lie derivatives and the invariant
product connection, with bounded lower-order differences.

Apply (9) to \(W=\partial_t^jX^\alpha V\), for
\(j+|\alpha|\le m-1\). Every such \(W\) has homogeneous Dirichlet/Neumann
conditions. Commutators with first-order coefficients contain at most
\(m\) total derivatives of \(V\), with at most one normal derivative.
Commuting time through \(c^2\Delta\) produces terms such as
\((\partial_t c^2)\Delta\partial_t^{j-1}X^\alpha V\), also of total order
at most \(m\), all tangential. These are controlled by the same collection
of energies. Summing and applying Gronwall controls every derivative of total
order at most \(m\) with zero or one normal derivative, by initial commuted
energies and \(\int J_{m-1}(\widetilde F)\).

Recover two or more normal derivatives from the equation:

\[
 V_{zz}=V_{tt}-c^2\Delta V+
       \widetilde A^\mu D_\mu V+\widetilde BV-\widetilde F.
 \tag{12}
\]

After \(q-2\) normal and the remaining time/angular differentiations, the
right side has at most \(q-1\) normal derivatives of \(V\); terms of the
same total order have fewer normal derivatives. Induction on \(q\) yields

\[
 J_m(V;t)\le C_m\left[
   \Big(\sum_{j+|\alpha|\le m-1}
               \mathcal E[\partial_t^jX^\alpha V](t)\Big)^{1/2}
          +J_{m-2}(\widetilde F;t)\right].
 \tag{13}
\]

The instantaneous source term is controlled by its initial value plus
\(\int J_{m-1}(\widetilde F)\), using the fundamental theorem of calculus
for each \(\partial_t^j\widetilde F\) in
\(H^{m-2-j}\). Initial time jets are recursively computed from the PDE;
their required norms are bounded by
\(\|V_0\|_{H^m}+\|V_1\|_{H^{m-1}}+J_{m-2}(\widetilde F;0)\).
Finally multiplication by smooth \(T^{\pm1}\) gives (11).

Only commuting angular derivatives would not prove this full statement.
The time commutations and normal recovery are essential. Likewise (11)
does not license replacing its source norm by just
\(\|F\|_{L^1_tH^{m-1}}\) without another boundary regularity argument.

For existence, conjugate first and solve the fixed homogeneous mixed problem
by the scalar mixed-boundary wave Galerkin construction or a weak energy
construction plus bounded first-order perturbation. The estimate (9) gives
uniqueness and convergence at the energy level. Tangential/time difference
quotients and (12) give (11) and regularity. Smooth compatible coefficients
and data give smooth solutions. This construction never diagonalizes \(M\).

## 6. Inhomogeneous boundary data: an explicit sufficient estimate

Write \(U_D=d\) and \(-U_{N,z}+MU_N=g\). A coupling to a prescribed
Dirichlet trace, if present in a more general problem, is first absorbed
in \(g\). Because \(T|_0=I\), the transformed data are
\(V_D=d\), \(-V_{N,z}=g\).

Choose a fixed smooth cutoff \(\eta\), equal to one near zero and supported
in the conjugation collar. The elementary lift

\[
 Q_D=\eta(z)d(t,\Omega),\qquad
 Q_N=-z\eta(z)g(t,\Omega)
 \tag{14}
\]

has exactly these data. Then \(W=V-Q\) has homogeneous boundary conditions
and forcing \(T^{-1}F-\widetilde P Q\). Apply (11) to this forcing and its
shifted initial data. In particular, with

\[
 K_r(b;t)=\sum_{j=0}^r
  \|\partial_t^jb(t)\|_{H^{r-j}(S^3)},\quad b=(d,g),
\]

a completely explicit sufficient addition to the right side of (11) is

\[
 C_m\left[K_m(b;0)+\int_0^{T_*}K_{m+1}(b;t)\,dt\right].
 \tag{15}
\]

This includes time derivatives through order \(m+1\), and counts \(d,g\)
with the same spatial regularity. It is intentionally not a sharp trace
estimate. It suffices for the smooth actual Taylor sources; it must not be
advertised as a lossless estimate for arbitrary rough Neumann forcing.

One optional improvement in spatial boundary regularity is explicit as well.
Put \(\Lambda=(1-\Delta^E_\kappa)^{1/2}\), where \(\Delta^E_\kappa\) is the
nonpositive connection Laplacian for the invariant, split-preserving positive
product connection on the boundary bundle, and use

\[
 Q_D=\eta e^{-z\Lambda}d,\qquad
 Q_N=\eta\Lambda^{-1}e^{-z\Lambda}g.
 \tag{16}
\]

The sign in (16) is correct since \(-\partial_zQ_N|_0=g\).
The spectral integral gives, for every integer \(r\ge0\),

\[
 \|Q_D\|_{H^r(\Sigma)}\le C_r\|d\|_{H^{r-1/2}(S^3)},\qquad
 \|Q_N\|_{H^r(\Sigma)}\le C_r\|g\|_{H^{r-3/2}(S^3)}.
 \tag{17}
\]

Consequently (15) can be replaced by the sufficient norm

\[
 \begin{split}
 &\sum_{j=0}^{m}\left(
 \|d^{(j)}(0)\|_{H^{m-j-1/2}}+
 \|g^{(j)}(0)\|_{H^{m-j-3/2}}\right)\\
 &+\int_0^{T_*}\sum_{j=0}^{m+1}\left(
 \|d^{(j)}(t)\|_{H^{m+1/2-j}}+
 \|g^{(j)}(t)\|_{H^{m-1/2-j}}\right)dt.
 \end{split}
 \tag{18}
\]

Negative spatial Sobolev indices in (18) are understood spectrally. The
statement first applies to smooth compatible data and then to their closure.
It is still a sufficient lifting-based estimate, not an optimal general
hyperbolic boundary-data theorem.

## 7. Corners, constraints and the actual angular cutoff

Compute the formal initial jets \(U^{[j]}\) from (2), \(U_0,U_1\) and
the source jets. At integer order \(m\), the corner identities whose traces
are defined are

\[
 \Pi_D U^{[j]}|_0=d^{(j)}(0),\quad 0\le j\le m-1,
\]

\[
 \left[-\partial_zU_N^{[j]}+
  \sum_{\ell=0}^j{j\choose\ell}M^{(\ell)}(0)U_N^{[j-\ell]}
 \right]_{z=0}=g^{(j)}(0),\quad 0\le j\le m-2.
 \tag{19}
\]

Smooth existence uses compatibility to every order. In the established
contact construction, initial perturbations and quadratic sources vanish
on an open initial shell neighborhood, so all these identities hold, after
the stipulated complete gauge and initial constraint construction.

The subsidiary covector has normal Dirichlet and tangential Robin data:

\[
 C_y|_b=0,\qquad
 \partial_y C_a-2a_b C_a=0.
 \tag{20}
\]

Its common scalar principal wave operator and smooth lower-order terms fit
the same theorem, with \(M_C=-2\rho_ba_b I\). Thus a zero-source subsidiary
field with zero initial value and time derivative is zero in the full
angular space. At Taylor order two, the exact Bianchi and scalar identities
make its bulk equation homogeneous; Codazzi and the scalar junction give
(20). The initial gauge jets and *trace-reversed physical initial
constraints* supply the two zero initial conditions. They are assumptions
to verify in any new application, not consequences of (9) alone.

The first gauge-vector problem has the same mixed structure and therefore
also admits the all-angular estimate. No extra tangential harmonic primal
rows are introduced: the primal count remains four Dirichlet plus twelve
normal derivative rows.

For the actual rotational background, complete SO(4) isotypic projections
on tensor/covector/scalar fields commute with the physical linear operator,
the boundary operator and the invariant conjugation. Applying (11) to such
finite projections therefore gives one constant for every cutoff. A
componentwise scalar harmonic cutoff is not being substituted for this
bundle projection. The proof already holds without any projection.

## 8. Domain of dependence and limits of the conclusion

Take initial support \(z\le L\), a fixed finite \(T_*\), and sources within
the corresponding radial causal region. The principal radial propagation
speed is one. Work on the actual smooth collar containing
\(0\le z\le L+T_*\), and extend coefficients smoothly farther inward only
to construct a global half-cylinder solution. Standard local versions of
the same energy identity imply radial finite propagation; that extension
cannot affect the physical solution. Neither the energy nor its constants
use a support-dependent Poincare inequality.

The constants can depend on \(T_*\), the lower bound for \(c\), the minimum
of \(\rho\) on the controlled collar, and the relevant background and shell
coefficient derivatives. They are not asserted uniform near the cone,
as \(T_*\to\infty\), or under a shrinking conjugation collar.

This removes the specific all-angular *frozen linear* gap for the actual
sixteen-field boundary formulation. It also controls its smooth
inhomogeneous Taylor sources and its subsidiary uniqueness argument. It
does not on its own construct exact nonlinear constrained initial data,
prove an exact perturbation family, or establish any novelty claim.
