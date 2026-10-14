# Full nonlinear shell Noether identity on a matched reference domain

14 October 2026 UTC. New authoring calculation, separate from the frozen
first-order capsule. The result is a classical variational identity and a
necessary quadratic condition for a hypothetical sufficiently regular exact
family. It is not a nonlinear existence, stability, sufficiency, or phase-space
completeness theorem. The accompanying review is internal, not external peer
review.

## 1. Result and scope

For two Einstein–scalar interiors with both outward GHY terms, exact induced
matching, and one common scalar-dependent shell tension, the missing internal
cut cancellation can be derived before restricting the perturbations to a
scalar master image. On a fixed reference worldvolume the exact ingredients are

\[
 \theta_\Sigma+\delta\ell=\alpha+d_\Sigma\beta,\qquad
 B_\xi:=\sum_i Q_{\xi,i}+i_\xi\ell-\beta(G_\xi)
       =-i_{I\cdot\xi}\epsilon_\Sigma,                 \tag{1}
\]

where \(I_{ab}\) is the full Israel residual, \(\alpha\) contains both full
junction residuals, and \(\beta\) is the unrestricted GHY corner. Consequently

\[
 k_\xi-(\delta\beta)(\delta,G_\xi)
 =\delta B_\xi-i_\xi\alpha+d_\Sigma(i_\xi\beta(\delta)). \tag{2}
\]

This identity includes the normal derivative of the gauge vector, which
cancels between its Komar and GHY-corner terms. The integral of the last term
vanishes on the closed three-dimensional cut. The remaining terms vanish on
the nonlinear junction equations and their linearization. Combining (2) with
the bulk Noether identity proves

\[
 \Omega^{\mathrm{rel}}_\Phi(\delta\Phi,G_\xi\Phi)=0       \tag{3}
\]

for every full linearized solution and every fixed pair of smooth reference
bulk vectors tangent to the shell and inducing the same vector there. No
Killing assumption is needed for (3). The ten compensated background Lorentz
generators are included. An exterior boundary, a cut edge, a cone excision, or
fixed physical shell data changes the conclusion by its retained charge.

The moving embedding is handled by exact pullbacks in Section 7; in particular
the \(i_\chi L\) and \(Q_\chi\) terms are displayed rather than discarded.
The proof of (3) is the actual relative Noether identity (1)–(2), not the fact
that an auxiliary pullback is invariant under simultaneous relabeling.

## 2. Classical domain and conventions

Take two compact spatial four-balls crossed with a sufficiently short time
interval, each having a smooth timelike boundary identified with one common
worldvolume \(W\). Their common spacelike cut \(S\) is a closed three-sphere.
The interiors include all smooth geometry inside the shell: neither a regular
characteristic cone nor its vertex is removed. There is no physical exterior
boundary. Work locally in time if needed; auxiliary initial/final surfaces are
the surfaces defining a variational potential, not extra physical walls.

A safe sufficient differentiability assumption is piecewise smooth spacetime
fields and smooth timelike embeddings, with the hypothetical family twice
continuously differentiable in its parameter in a \(C^k\) topology on each
closed reference slab for some \(k\geq4\). Assume uniform nondegeneracy of the
metrics and induced shell metric, a collar for each embedding, and exact
matching throughout the family. This controls the displayed traces, products,
pullbacks, and two parameter derivatives on compact sets. We do not assert
that this is an optimal regularity class or that a family in it exists.
Weaker energy or Sobolev spaces require a separate trace and approximation
argument. The same calculation is valid for individual smooth solutions and
smooth full linearized solutions without assuming a smooth solution manifold.

Use signature \((-++++)\), \(8\pi G_5=1\), outward spacelike unit normals
\(n_i\), and

\[
 K_{i,ab}=h_a{}^A h_b{}^B\nabla_A n_{i,B},\quad
 p_i=n_i^A\nabla_A\phi_i,
\]
\[
 S=\sum_{i=\pm}\int_{M_i}\epsilon_i
       [R_i/2-(\nabla\phi_i)^2/2-U(\phi_i)]
   +\int_W\epsilon_\Sigma(K_++K_--\sigma(\varphi)).       \tag{4}
\]

Exact matching is part of configuration space, not an equation omitted from
the argument:

\[
 X_+^*g_+=X_-^*g_-=h,\qquad
 X_+^*\phi_+=X_-^*\phi_-=\varphi.                        \tag{5}
\]

Thus every displayed \(\delta h\) and \(\delta\varphi\) is common to both
sides. Differentiating (5) twice includes all evaluation, bending, and mixed
terms; it is not replaced by equality of second-order bulk fields at the old
surface.

**Orientation convention.** On each side \(\epsilon_{\Sigma i}=i_{n_i}
\epsilon_i\). Boundary forms are represented in one common *outward
co-oriented density convention* on \(W\), denoted \(\epsilon_\Sigma\), and
the same conversion is applied to \(\theta_i,Q_i,\beta_i\), not only to GHY.
Equivalently choose sidewise orientations mapping their outward boundary forms
to that common representative. If one instead insists on one global oriented
doubled atlas, let \(X_i^*\epsilon_{\Sigma i}=s_i\epsilon_\Sigma\); insert
\(s_i\) into every converted boundary form and the corresponding cut integral.
Unsigned Israel sums must never be combined with raw oppositely oriented form
pullbacks. All formulae below use the co-oriented convention.

Orient a future Cauchy section \(C_i\) as a future boundary of a slab and its
cut as \(S_i=\partial C_i\). The future-end orientation of the timelike side
is \(-S_i\), because \(i_u i_n\epsilon_i=-i_n i_u\epsilon_i\).
This is the reason for the minus sign in the relative corner potential below.
For an orthogonal local chart \(\epsilon_i=dt\wedge dy\wedge d^3x\),
\(n_i=\partial_y\), one has \(\epsilon_\Sigma=-dt\wedge d^3x\) and
\(S_i\) has orientation \(+d^3x\). These explicit signs also fix the
constraint normalization in Section 9.

## 3. Unrestricted first variation and the actual GHY corner

Write \(q_{AB}=\delta g_{AB}\), \(q=g^{AB}q_{AB}\), and
\(q^{AB}=g^{AC}g^{BD}q_{CD}\) with no minus sign in this definition.
The fixed bulk Lee–Wald potential is exactly the prior capsule's convention:

\[
 \theta^A=\sqrt{-g}\left[
   \tfrac12(\nabla_Bq^{AB}-\nabla^Aq)-\nabla^A\phi\,\delta\phi
 \right],\qquad
 \omega(\delta_1,\delta_2)=\delta_1\theta(\delta_2)
                         -\delta_2\theta(\delta_1).      \tag{6}
\]

On a fixed timelike boundary put

\[
 v_i^a=h^{ab}n_i^Aq_{i,Ab},\qquad
 \beta_i(\delta)=-\tfrac12 i_{v_i}\epsilon_\Sigma.        \tag{7}
\]

No radial gauge condition is placed on \(q_{nn}\) or \(q_{na}\). Direct
normalization of \(n\), variation of \(K\), and tangential integration by
parts give the exact local identity

\[
 \theta_i|_W+\delta(K_i\epsilon_\Sigma)
 =-\epsilon_\Sigma\left[
     \tfrac12(K_i^{ab}-K_i h^{ab})\delta h_{ab}
      +p_i\delta\varphi\right]+d_W\beta_i.             \tag{8}
\]

One can verify all signs in Gaussian normal coordinates for the background,
while leaving the variation unrestricted. Set
\(g=dr^2+h_{ab}dx^adx^b\), \(\lambda=q_{rr}\),
\(v_a=q_{ra}\), and \(q_h=h^{ab}q_{ab}\). Then

\[
 \delta K=\tfrac12\partial_rq_h-\tfrac12\lambda K-D_av^a,
\]
\[
 \theta_g|_W=\tfrac12[D_av^a-K^{ab}q_{ab}
                           +K\lambda-\partial_rq_h]\epsilon_\Sigma,
\]
\[
 \delta(K\epsilon_\Sigma)=
 [\tfrac12\partial_rq_h-\tfrac12\lambda K-D_av^a
                          +\tfrac12Kq_h]\epsilon_\Sigma.
\]

Their sum is (8). In particular the unrestricted corner is present even
though it vanished on the earlier conformal zero-shift radial calculation.

Define the full nonlinear residuals

\[
 I_{ab}=\sum_i(K_{i,ab}-K_i h_{ab})+\sigma h_{ab},\qquad
 P=p_++p_-+\sigma'(\varphi),                            \tag{9}
\]
\[
 \ell=(K_++K_--\sigma)\epsilon_\Sigma,\qquad
 \beta=\beta_++\beta_-,\qquad
 \alpha(\delta)=-\epsilon_\Sigma
       [\tfrac12 I^{ab}\delta h_{ab}+P\delta\varphi].    \tag{10}
\]

Summing (8) and varying the **single** tension gives (1)'s first identity.
The scalar term \(\sigma'\delta\varphi\) is indispensable. The shell has
no separately independent intrinsic kinetic field in (4).

## 4. Nonlinear embedding consistency, including the normal equation

Let \(E_{i,AB}=G_{i,AB}-T_{i,AB}\), with
\(T_{AB}=\nabla_A\phi\nabla_B\phi-
g_{AB}[(\nabla\phi)^2/2+U]\). Common induced fields imply the exact
Gauss and Codazzi residual identities

\[
 E_{i,n_in_i}=\tfrac12[K_i^2-K_{i,ab}K_i^{ab}-R[h]
                         -p_i^2+(D\varphi)^2+2U],       \tag{11}
\]
\[
 \sum_i E_{i,n_i a}=D_bI^b{}_a-PD_a\varphi.             \tag{12}
\]

For \(\Delta K_{ab}=K_{+,ab}-K_{-,ab}\),
\(\Delta p=p_+-p_-\), and \(\Delta K=h^{ab}\Delta K_{ab}\),
subtracting (11) gives

\[
 \boxed{\sigma\Delta K+\sigma'\Delta p
 =2(E_{+,n_+n_+}-E_{-,n_-n_-})
       +I_{ab}\Delta K^{ab}+P\Delta p.}                \tag{13}
\]

To check the factors, in four shell dimensions
\(\sum K=(4\sigma-\mathrm{tr}I)/3\) and
\(\sum K_{ab}=I_{ab}+(\sigma-\mathrm{tr}I)h_{ab}/3\).
Also \(\sum p=P-\sigma'\). These identities inserted in the difference
of (11) give (13) without expansion about a background.

The usual mean normal force convention divides the left side of (13) by two.
Its vanishing is unaffected. Thus the full normal force and tangential shell
conservation follow from the full bulk constraints and the two junctions.
They are not additional equations silently dropped from the variational
domain, and their quadratic corrections have not been inferred from a
first-order radial identity. In a moving-embedding presentation, these are
the corresponding diffeomorphism relations among embedding and field
Euler–Lagrange equations.

## 5. Exact tangent Noether cut identity

Fix smooth vector fields \(\xi_i\) on the reference bulks, tangent to their
fixed shell boundaries, with the same restriction \(\xi\) on \(W\). They
may have arbitrary normal derivatives. Let

\[
 G_\xi g_i=\mathcal L_{\xi_i}g_i,\qquad
 G_\xi\phi_i=\mathcal L_{\xi_i}\phi_i,
\]
\[
 Q_{\xi,i}=-\tfrac12 *_i d\xi_i^\flat.                  \tag{14}
\]

Minimal scalar coupling adds no separate \(Q\). For a single side let
\(A_i^a=h^a{}_B n_i^A\nabla_A\xi_i^B\). Tangency gives
\(h^a{}_B n_{i,A}\nabla^B\xi_i^A=-(K_i\cdot\xi)^a\), so

\[
 Q_{\xi,i}|_W=-\tfrac12 i_{A_i+K_i\cdot\xi}\epsilon_\Sigma,
 \qquad
 \beta_i(G_\xi)=-\tfrac12 i_{A_i-K_i\cdot\xi}\epsilon_\Sigma.
\]

Subtracting the second expression from the first and adding
\(i_\xi K_i\epsilon_\Sigma\) cancels \(A_i\) exactly. Both copies and the
single \(-\sigma\epsilon_\Sigma\) therefore give

\[
 \boxed{B_\xi:=\sum_iQ_{\xi,i}|_W+i_\xi\ell-\beta(G_\xi)
             =-i_{I\cdot\xi}\epsilon_\Sigma.}           \tag{15}
\]

This is an explicit residual potential, not an assumed zero boundary charge.
In particular, dropping \(\beta(G_\xi)\) leaves a spurious dependence on the
normal extension of \(\xi\).

Set \(k_\xi(\delta)=\sum_i(\delta Q_{\xi,i}-i_{\xi_i}\theta_i)|_W\).
Because \(\xi\) is field independent and the potential is covariant,

\[
 (\delta\beta)(\delta,G_\xi)
       =\delta[\beta(G_\xi)]-\mathcal L_\xi\beta(\delta).
\]

Use Cartan's identity and (10) to obtain the fully off-shell algebraic result

\[
 \boxed{k_\xi-(\delta\beta)(\delta,G_\xi)
   =\delta B_\xi-i_\xi\alpha(\delta)
                   +d_W(i_\xi\beta(\delta)).}          \tag{16}
\]

The explicit residual expression is

\[
 \delta B_\xi-i_\xi\alpha
 =-i_{(\delta I^a{}_b)\xi^b}\epsilon_\Sigma
  -\tfrac12(h^{cd}\delta h_{cd})i_{I\cdot\xi}\epsilon_\Sigma
  +[\tfrac12I^{ab}\delta h_{ab}+P\delta\varphi]
                                      i_\xi\epsilon_\Sigma. \tag{17}
\]

On \(I=P=0\), and for a full linearized solution \(\delta I=\delta P=0\),
this is zero. The exact differential in (16) integrates to zero because
\(\partial S=\varnothing\); if the cut has an edge its integral is a real
edge term. No omitted corner or independent shell charge remains.

## 6. Relative potential, off-shell generator, and gauge degeneracy

The timelike boundary's future-end sign specified in Section 2 gives

\[
 \Theta_C^{\mathrm{rel}}(\delta)
   =\sum_i\int_{C_i}\theta_i(\delta)-\int_S\beta(\delta),
 \qquad \Omega_C^{\mathrm{rel}}=\delta\Theta_C^{\mathrm{rel}}. \tag{18}
\]

This follows by inserting (10) in the slab first variation and using
\(\int_Wd\beta=-\int_{S_f}\beta+\int_{S_i}\beta\).
It includes the GHY corner once, not twice. Curling (10) proves conservation
of (18) for two full linearized solutions on a nonlinear shell solution.

For clarity the complete bulk residual form and constraint current are

\[
 \boldsymbol E_i\cdot\delta\Phi_i
  =\epsilon_i[-\tfrac12 E_i^{AB}\delta g_{i,AB}
                   +(\Box\phi_i-U')\delta\phi_i],
 \qquad
 C_{\xi,i}=E_i{}^A{}_B\xi_i^B\epsilon_{i,A}.             \tag{19}
\]

With (6) and (14), \(J_\xi=\theta(G_\xi)-i_\xi L=C_\xi+dQ_\xi\).
The exact fixed-parameter bulk identity is

\[
 \omega_i(\delta,G_\xi)=\delta C_{\xi,i}
      +d(\delta Q_{\xi,i}-i_{\xi_i}\theta_i)
      +i_{\xi_i}(\boldsymbol E_i\cdot\delta\Phi_i).      \tag{20}
\]

Define the full smeared relative constraint functional by its derived form

\[
 \boxed{H_\xi=\sum_i\int_{C_i}C_{\xi,i}+\int_S B_\xi
 =\sum_i\int_{C_i}E_i{}^A{}_B\xi_i^B\epsilon_{i,A}
             -\int_S i_{I\cdot\xi}\epsilon_\Sigma.}     \tag{21}
\]

Combining (16), (18), and (20) yields

\[
 \boxed{\Omega^{\mathrm{rel}}(\delta,G_\xi)
 =\delta H_\xi
   +\sum_i\int_{C_i}i_{\xi_i}
                         (\boldsymbol E_i\cdot\delta\Phi_i)
   -\int_S i_\xi\alpha(\delta).}                       \tag{22}
\]

The identity is off shell on the exact matching configuration space. It is
not obtained by guessing an Israel term to append to a bulk Einstein
integral. Scalar matching and the scalar junction enter \(\alpha\), although
there is no separate scalar diffeomorphism constraint in (21). Normal force
is related to these residuals by (13).

On a full nonlinear solution, with \(\delta\) a full tangent solution,
\(\delta C_\xi=0\), \(\delta B_\xi=0\), and all residual terms vanish.
This proves (3) for arbitrary common shell-tangent reference generators.
There is no reliance on separation of variables, reflection parity, a master
equation, or the earlier restricted current balance.

## 7. Moving domains and embeddings are retained

Locally in parameter and time, a smooth family of timelike embedded
boundaries with a collar can be represented by diffeomorphisms
\(Y_i:\bar M_i\to M_i\), restricting to \(X_i\) on the common fixed
reference \(W\). They can be extended inside the two compact balls. This is
a local change of variables, not a nonlinear field-equation theorem.

Write \(\bar g_i=Y_i^*g_i\), \(\bar\phi_i=Y_i^*\phi_i\), and
\(\chi_i=\delta Y_i\circ Y_i^{-1}\). Then exactly

\[
 \delta\bar\Phi_i=Y_i^*(\delta\Phi_i+\mathcal L_{\chi_i}\Phi_i). \tag{23}
\]

The action (4) on moving domains is identically its barred version on fixed
domains, including both GHY terms and the common tension evaluated on the
pulled-back common fields. Apply Sections 3–6 to the barred fields. Their
unrestricted potential in original variables is

\[
 \theta_i(\bar\Phi_i;\delta\bar\Phi_i)
 =Y_i^*\{\theta_i(\Phi_i;\delta\Phi_i)
          +i_{\chi_i}L_i+C_{\chi_i,i}+dQ_{\chi_i,i}\}.   \tag{24}
\]

Thus

\[
 \boxed{\Theta_C^{\mathrm{rel}}
 =\sum_i\int_{Y_i(\bar C_i)}
       [\theta_i(\delta\Phi_i)+i_{\chi_i}L_i+C_{\chi_i,i}]
   +\int_S\sum_iQ_{\chi_i,i}
   -\int_{\bar S}\bar\beta(\delta\bar\Phi).}           \tag{25}
\]

Pullbacks and orientation conversions are understood in the two displayed
cut terms. Equation (25) is exact off shell; \(C_\chi\) can be removed only
on shell (and from its curl only when the linearized constraints hold).
The \(i_\chi L\) term records the moved integration domain; matter generally
makes it nonzero. The \(Q_\chi\) and unrestricted \(\bar\beta\) terms record
the moving gravitational cut. The common shell is already present in the
barred boundary identity; adding a second shell volume or another GHY corner
would double count it. When taking a field-space curl, vary the pullbacks,
the cut locations, and \(\chi\) consistently. Formula (18) in barred
variables is a simpler unambiguous definition of that curl.

An ambient gauge transformation has \(\delta_K\Phi_i=\mathcal L_{K_i}\Phi_i\)
and \(\delta_KY_i=-K_i\circ Y_i\). A common worldvolume reparametrization
\(k\), extended by fixed reference vectors \(\xi_i\) tangent to the
reference boundary, adds \((Y_i)_*\xi_i\). The combined transformation gives

\[
 G_{K,k}\bar\Phi_i=\mathcal L_{\xi_i}\bar\Phi_i.         \tag{26}
\]

Its shell restrictions agree precisely because \(k\) is common. Equation
(22), applied to (26), proves the required gauge degeneracy, including the
embedding action. Different choices of interior extension of \(Y_i\) or of
\(\xi_i\) differ by reference diffeomorphisms with zero shell restriction;
(22) shows that they are null directions on shell. No relative boundary
charge was assumed absent to establish that statement.

The physically needed assertion is for these common transformations of the
full matched system. We do not assert that an arbitrary vector normal to a
*fixed* reference wall, with the wall held fixed, is an admitted gauge vector.
Normal physical motion is represented by (23) and is not set to zero.

## 8. Reducibility and the necessary quadratic condition

At the regular nonconstant background

\[
 g_0=\Omega(s)^2\eta,\quad \phi_0=\Phi(s),\quad
 s=X\cdot X,\quad \Phi_s\ne0,\quad W:\ s=b^2,
\]

the common continuous metric–scalar isometries are the ten Lorentz
generators. Indeed \(K\phi_0=0\) implies \(Ks=0\), and then
\(\mathcal L_Kg_0=0\) implies \(\mathcal L_K\eta=0\); a flat Killing
translation has \(Ks=2t\cdot X\), so it vanishes on an open set only when
\(t=0\). The shell identification allows only the common Lorentz action.
This does not classify any additional adjoint kernel in an unspecified weak
constraint domain.

Fix one such \(K\), and fix the common background worldvolume vector
\(k_0=K|_W\). The full generator is

\[
 G_{K,k_0}g_i=\mathcal L_Kg_i,\quad
 G_{K,k_0}\phi_i=K\phi_i,\quad
 G_{K,k_0}X_i=-K\circ X_i+(X_i)_*k_0.                    \tag{27}
\]

Then \(G_{K,k_0}(\Phi_0)=0\), including the embeddings. Keep \(K,k_0\)
fixed in field space. In radial/de Sitter coordinates \(K^y=0\) and its
tangential components are independent of \(y\). For a normal graph
\(X_i(x)=(y_b+\zeta_i(x),x)\), (27) gives
\(G\zeta_i=K\zeta_i\); hence
\(\mathcal K v:=DG_{K,k_0}|_0[v]\) is the usual Lie derivative on the
reconstructed fields and their graph displacements. A field-dependent
compensator would require additional variation terms and is not used.

Suppose an exact parameter-\(C^2\) family in Section 2's class exists, with
\(v=\Phi'_0\), \(w=\Phi''_0\). Apply (3) along it and differentiate:

\[
 0=(D\Omega_0[v])(v,G_0)+\Omega_0(w,G_0)
                      +\Omega_0(v,DG_0[v])
   =\Omega^{\mathrm{rel}}_0(v,\mathcal K v).             \tag{28}
\]

This is now based on the unrestricted nonlinear cut calculation, not on
postulating its conclusion. Equivalently, (22) at a reducible solution gives
\(DH_\xi|_0[w]=0\) for every smooth exact-matching variation \(w\).
Differentiate (22) once in the direction of a full linearized solution \(v\).
The derivatives of its bulk and shell residual terms vanish because
\(DE[v]=0\), \(D(\Box\phi-U')[v]=0\), \(DI[v]=0\), and \(DP[v]=0\).
Consequently

\[
 \boxed{D^2H_\xi|_0[v,v]
       =\Omega^{\mathrm{rel}}_0(v,\mathcal K v).}         \tag{29}
\]

This Hessian is taken in a local chart on the *exact matching* configuration
space, using the barred action. Because \(DH_0=0\), it is independent of a
second-order chart correction. From exact constraints along the family,
\(0=DH_0[w]+D^2H_0[v,v]\), proving (28) again.
The conventionally halved quadratic Taub functional is
\(T_\xi(v)=\tfrac12D^2H_\xi[v,v]\). Its surface source is explicitly the
second variation of the surface integral in (21), with barred geometry;
all second-order matching, tension, normal-normalization, and bending terms
are thereby included. An expanded component formula for every unrestricted
second-order perturbation is not supplied by this note.

## 9. Relation to the previously evaluated physical pairing

The bulk potential (6), action normalization, and unrestricted GHY identity
are the same as in the frozen packet calculation. There is an overall
Cauchy-orientation conversion to state explicitly. In a time chart, the scalar
part of (18) is \(\delta\pi_\phi\wedge\delta\phi\), so it pairs two data as
\(p_1q_2-p_2q_1\). The prior future-Cauchy convention gives
\(q_1p_2-q_2p_1\). Define the nonlinear future-Cauchy convention by

\[
 \Omega_{\mathrm{future}}:=-\Omega^{\mathrm{rel}},\qquad
 T_\xi(v)=-\tfrac12\Omega_{\mathrm{future}}(v,\mathcal K v). \tag{30}
\]

Its equality to the frozen reconstructed \(\Omega_{\mathrm{phys}}\) is
established below for the compact initial-support states used here. Equation
(30) alone does not identify arbitrary after-contact component currents.

For an orthogonal slice, (19) gives
\(\int_C C_\xi=-\int_C\sqrt q\,u^AE_{AB}\xi^B\).
For a vector tangent to the shell cut,
\(\int_SB_\xi=-\int_S\sqrt c\,I(\tau,\xi)\).
These signs match the independently derived integrated ADM momentum
constraint. Changing both signs by a global convention changes no necessary
zero condition.

There is a particularly direct, ambiguity-free comparison for the states at
issue. Their initial master data, full reconstructed perturbations, and shell
displacements vanish on an open neighborhood of the initial shell cut. The
same is true after applying a background Lorentz generator: its action is
local and does not enlarge support. Because the shell displacement and its
needed time jet vanish there, choose the freely extendible first variation
of the straightening map to be zero on a neighborhood of the **entire initial
Cauchy section**, including the time jet needed for the current. Merely
requiring it to vanish at the cut would not by itself eliminate the bulk
\(i_\chi L\) terms. With this choice, every term in (25) beyond the ordinary
bulk pairing vanishes at quadratic order for \((v,\mathcal K v)\); commuting
second parameter variations of the extension cancel in the antisymmetric
curl. Thus \(\Omega_{\mathrm{future}}=\Omega_{\mathrm{phys}}\) for the
initial states compared here. This is the previously computed physical
pairing, not an auxiliary scalar energy.

For general perturbations touching the cut, (25) supplies the unrestricted
completion. The old odd-graph expression

\[
 \omega_\Sigma^\mu=18c_o\sqrt{-\gamma}
   (P_1D^\mu P_2-P_2D^\mu P_1)
\]

was derived from the same action on its restricted graph family, with the
same moved-domain contributions. This note does **not** claim a new
component-by-component rederivation of that restricted formula from (25).
That extra comparison is unnecessary for the initial support argument above;
it would be required before using a different cut-potential convention to
identify arbitrary after-contact component currents without further work.
The previous total-current conservation, with its retained odd cut term,
already transports its initial pairing through contact on that family.

Hence the single-real-angular-harmonic calculation now clears the ten
necessary quadratic symmetry tests, within the stated classical family
hypotheses. Conversely, any admitted linearized state with nonzero initial
\(\Omega_{\mathrm{phys}}(v,\mathcal K v)\) fails that necessary condition.
Neither ten zeros nor a nonzero mutual pairing between unrelated states
proves nonlinear integrability. Solving second-order bulk constraints,
junction compatibility, other adjoint range conditions, and higher orders
remains separate work; no nonlinear closure of the scalar master image is
asserted.

## 10. Verification and provenance

The derivation uses the action and Lee–Wald convention recorded in
`junctions/FULL_PACKET_JUNCTION_SUFFICIENCY_AND_CURRENT.md` and
`basis/prior-checkpoint/covariant/COVARIANT_CURRENT_AND_CONE_LIMITS.md` of the
frozen regular-shell-contact capsule, and directly addresses the missing
identities isolated in the two private second-order forward notes.

`cut-sign-review.md` independently derives the unrestricted GHY sign, the
normal-extension cancellation in (15), and the exact residual expression
(17). The accompanying exact-symbolic checker tests the four-dimensional
trace conversion, Lorentzian normal-force factorization, Codazzi residual
factorization, the GHY corner coefficient, and the two-copy/single-tension
charge cancellation with deliberately wrong alternatives. It is an algebra
check, not a PDE existence or machine-checked proof of the differential-form
argument. Every normal and optimized-Python run must report the same JSON.

The construction follows standard covariant Noether and relative-boundary
variational methods. No originality claim is made. A separate source note
records primary literature and the distinctions between fixed external
boundaries, auxiliary embedding variables, and this physical internal shell.
