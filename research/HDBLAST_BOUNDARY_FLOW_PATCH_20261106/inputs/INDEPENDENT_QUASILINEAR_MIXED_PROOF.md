# Independent proof of the quasilinear mixed shell evolution step

Internal derivation, 15 October 2026. This note was developed without reading
the parallel Astra derivation. It proves an analytic wave-system statement;
the identification of the resulting reduced solution with the physical
Einstein–scalar solution still requires the initial-gauge and subsidiary
arguments. Those are separate model-specific inputs. No computation below is
a substitute for the estimates.

## 1. Precise analytic problem and statement

Let \(X=[0,Z]\times S^3\), with a fixed smooth product measure and a fixed
positive bundle metric on \(E=E_D\oplus E_N\), of ranks four and twelve.
The orthogonal splitting and the product connection \(D\) are fixed and
parallel. They are not the physical Lorentzian tensor norm or connection.
Changing to this connection introduces only smooth first-order terms.

Consider

\[
 g^{\mu\nu}(t,x,u)D_\mu D_\nu u=F(t,x,u,Du),\qquad
 u_D|_{z=0}=0,\quad D_z u_N|_{z=0}=0.                 \tag{1}
\]

The principal coefficients are real scalars, symmetric in their indices,
smooth in ((t,x,u)), and do not depend on \(Du\). The function \(F\) is
smooth in all its displayed variables and \(F(t,x,0,0)=0\). The background
is the zero state. On the fixed background slab \(0\le t\le T\), suppose a
fixed state neighborhood has

\[
 -g^{00}\ge\lambda,
 \qquad g^{ij}\xi_i\xi_j\ge\lambda|\xi|^2,
 \qquad |g|\le\lambda^{-1},                            \tag{2}
\]

for a fixed positive \(\lambda\). Shrink that neighborhood as necessary.
In addition, at the shell,

\[
 u_D=0\quad\Longrightarrow\quad g^{za}(t,0,\Omega,u)=0
 \quad(a=t,\Omega).                                    \tag{3}
\]

This is an identity for every boundary state used in the iteration.
Timelikeness of \(dt\) alone would not replace the spatial positivity in
(2). The exact state-space shell flow established in the previous capsule
supplies (1), including precisely homogeneous boundary conditions, and
supplies (3) from the four coordinate Dirichlet rows.

For integer \(s\ge1\) define the spatial state norm

\[
 A_s(u;t)=\|u(t)\|_{H^s(X)}+\|u_t(t)\|_{H^{s-1}(X)}.
                                                               \tag{4}
\]

Fix \(m\ge6\), a fixed finite \(T\), and a common radial support region
strictly inside the outer wall. Suppose the initial perturbation data are
in \(H^m\times H^{m-1}\), vanish in an open initial shell collar, and have a
protected short Cauchy evolution with these data. Its support is separated
from the shell for a short interval. The retained precontact Cauchy theorem
is one way to supply this last input. Equivalently, it suffices to supply a
trial seed with the properties proved in section 4 below.

Then sufficiently small initial \(A_m\) norm, with smallness depending on
the fixed background slab, gives a unique solution of (1) on the whole
prescribed interval \(0\le t\le T\), with

\[
 \partial_t^j u\in C([0,T];H^{m-j}(X)),\qquad 0\le j\le m.
                                                               \tag{5}
\]

The state stays inside the neighborhood in (2), and the actual
characteristic support stays inside the outer causal margin. Data
dependence is continuous in (4) and Lipschitz one derivative lower on
bounded sets with a common shell collar and common protected seed interval.
The registered family has precisely these common supports and intervals.
No top-norm differentiability of the
solution map is asserted. In particular \(m=14\) applies to completed
\(H^{14}\times H^{13}\) reduced initial data. The fixed-\(m\) statement
does not produce a common spatially smooth parameter family.

The proof first uses smooth shell-flat data, and then passes to finite
Sobolev data by the compatible smoothing argument in section 8.

## 2. Variable-coefficient energy and genuine linear existence

Write \(a=-g^{00}>0\), \(b^i=g^{0i}\), \(c^{ij}=g^{ij}\). For the frozen
operator \(L=-aD_t^2+2b^iD_iD_t+c^{ij}D_iD_j\), with any bounded smooth
first-order matrix terms included if desired, use

\[
 E[w]=\frac12\int_X
   \{a|w_t|^2+c^{ij}\langle D_iw,D_jw\rangle+|w|^2\}.
                                                               \tag{6}
\]

Multiply the equation by \(-w_t\). The mixed time/space term contributes
a boundary multiple of \(b^z|w_t|^2\); the spatial term contributes a
multiple of \(\langle w_t,c^{zj}D_jw\rangle\). At the shell (3) reduces
their sum to \(c^{zz}\langle w_t,w_z\rangle\). It vanishes because
\(w_{D,t}=0\), \(w_{N,z}=0\), and the fixed splitting is orthogonal. Thus

\[
 \frac d{dt}\sqrt{E[w]}
 \le C_\lambda(1+\|D_{t,x}g\|_\infty+\|B\|_\infty)\sqrt{E[w]}
       +C_\lambda\|Lw\|_2.                            \tag{7}
\]

Here \(B\) denotes the displayed lower-order coefficients. The estimate
uses an ordinary positive fiber metric, without a putative symmetrizer of
the original twelve-row Robin matrix. At an artificial outer wall impose
full homogeneous Dirichlet conditions; \(w_t=0\) kills all its fluxes.

For clarity, smooth frozen linear existence is not being inferred just
from (7). The following construction supplies it. Let \(P_N\) be finite
orthogonal invariant angular projections, separately on \(E_D,E_N\),
exhausting complete finite-dimensional SO(4) representation subspaces.
They commute with time and normal differentiation and with the SO(4)
action, preserve radial support, and have \(L^2\) norm at most one. Solve

\[
 P_N Lw_N=P_N f,\qquad w_N=P_Nw_N.                    \tag{8}
\]

No variable coefficient is asserted to commute with \(P_N\). At each
fixed \(N\) this is a finite-component one-dimensional radial wave system.
Its time mass \(P_NaP_N\) and radial stiffness
\(G_N=P_Nc^{zz}P_N\) are uniformly coercive. The first-order radial reduction
has a positive time matrix and symmetric radial matrix; at the shell the
space \(w_{D,t}=0,w_{N,z}=0\) is a maximal zero-flux space. One may construct
the system directly by radial variational Galerkin, with Dirichlet rows
essential and Neumann rows natural.

For this radial construction at a fixed \(N\), first compute that finite
angular system's exact smooth initial jets and subtract a smooth
boundary-preserving Borel lift of them. Its remainder has zero initial
data and forcing flat in time at zero. Thus radial Galerkin does not have
to project nonzero initial data or their powers through a radial stiffness
operator. The lift is used only for smooth existence at fixed \(N\); its
norm need not be uniform in \(N\). All uniform estimates below are applied
directly to the resulting solution with its original projected data and
forcing. Initial jet bounds for those estimates follow by their finite
triangular recursion, using the coercive inverse of \(P_NaP_N\) and
commuting its differentiated equation, exactly as in (9).

Here is the regularity step in that construction, so that it does not
assume \(w_{zz}\in L^2\) before proving it. Put the radial equation in
divergence form and differentiate it \(p\) times in time inside the finite
radial Galerkin space. For \(v_p=\partial_t^pw\), the differentiated radial
principal form produces terms

\[
 -\sum_{k=1}^p{p\choose k}
 \int\langle \partial_t^kG_N\,\partial_z\partial_t^{p-k}w,
                    \partial_z v_{p,t}\rangle.
\]

Add the corresponding cross terms with \(\partial_zv_p\) to the energy.
Their time derivatives cancel precisely these highest terms. The remaining
terms contain at most \(\partial_zv_p\) and lower time-level energies.
Each cross term is bounded by
\(\eta E[v_p]+C_\eta\sum_{q<p}E[v_q]\). Recursively adding large positive
multiples of lower energies gives a positive equivalent energy satisfying
an ordinary Gronwall inequality. Mixed principal terms are handled by the
same integration by parts; the boundary contributions vanish by (3).
All remainder Galerkin time jets at zero are zero. Hence
the radial weak limit has every time derivative in \(H^1\times L^2\).
Its equation and \(G_N^{-1}\) then give \(w_{zz}\in L^2\), and induction
gives a smooth fixed-\(N\) solution. This closes the finite-dimensional
existence step before any claim of uniform angular regularity.

Commute (8) with time and SO(4) generators. The commutator is exactly
\(P_N[L,\mathcal T^I]w_N\); projection disappears from the energy pairing
because the test field belongs to its range. Normal recovery uses

\[
 G_N\partial_z^2w_N=P_N(\text{all other terms}),\qquad
 \|G_N^{-1}\|_{L^2_\Omega\to L^2_\Omega}\le\lambda^{-1}. \tag{9}
\]

Differentiating this identity only differentiates its displayed factors;
it does not require a uniform high Sobolev multiplier norm of the inverse.
Consequently all estimates in sections 3 and 5 hold uniformly in \(N\).
For smooth coefficients and data their higher-order versions and compactness
give a smooth all-angular limit. Testing against finite angular sums
passes the equation, and (7) proves uniqueness.

The support-flat application has no projection/corner defect: initial
data and every initial forcing jet vanish on an open radial shell collar,
and \(P_N\) preserves that fact. For arbitrary smooth compatible linear
data one may first subtract a boundary-preserving smooth Borel lift of the
exact initial time jets. The remainder has zero initial data and forcing
flat at zero, so (8) is again compatible. A bare projection of arbitrary
compatible data would not have this property. This general smooth
extension is not needed in the prescribed-interval construction below.

## 3. Tangential estimates and normal recovery

Let \(\mathcal T\) consist of time differentiation and the infinitesimal
SO(4) action on the fixed bundle. These operators preserve the D/N split
and traces; the angular generators span the tangent directions of \(S^3\).
Define

\[
 e_m(w)^2=\sum_{|I|\le m-1}E[\mathcal T^Iw].          \tag{10}
\]

For every word \(I\), the commuted forcing consists of \(\mathcal T^If\)
and the full commutator. Its principal terms include every

\[
 (\mathcal T^\beta g^{\mu\nu})
 D_\mu D_\nu\mathcal T^{I-\beta}w,\qquad
 1\le|\beta|\le|I|,                                  \tag{11}
\]

as well as bounded connection/action commutators. None of these is dropped.
All commuted fields have zero boundary flux.

Recover normal derivatives by solving the equation for \(g^{zz}w_{zz}\).
To estimate a derivative with \(r\ge2\) normal factors and total order at
most \(m\), apply \(D_z^{r-2}\mathcal T^I\) to that identity. Every other
top principal term has fewer normal factors. A differentiated coefficient
multiplies a derivative of \(w\) of strictly smaller total order. Induct
first in total order, then in normal count. The same proof with (9) works
at angular cutoff. This is an actual normal-recovery argument; it does not
give homogeneous boundary conditions to normal derivatives.

For arbitrary smooth fixed coefficients it gives the sufficient estimate

\[
 \sup_t\sum_{j=0}^m\|\partial_t^jw(t)\|_{H^{m-j}}
 \le C_m\Big[A_m(w;0)
 +\sum_{j=0}^{m-1}\|\partial_t^jf\|_{L^1_tH^{m-1-j}}
 +\sum_{j=0}^{m-2}\|\partial_t^jf(0)\|_{H^{m-2-j}}\Big]. \tag{12}
\]

The coefficient dependence of \(C_m\) must not be used as a quasilinear
argument. The explicitly tame version needed for that purpose is next.

## 4. Compatible seed and iteration

Use a smooth short precontact Cauchy solution \(u^{\mathrm{pre}}\), equal
to the background in a shell wedge, on a protected interval ([0,h]).
Choose a fixed smooth time cutoff equal to one near zero and vanishing
before \(h\), multiply the solution perturbation by it, and extend by zero.
Call this seed \(u^0\). It satisfies the homogeneous mixed traces, is zero
in the initial shell wedge, and has exactly the nonlinear initial jets
through every available finite order. It obeys

\[
 \sup_{[0,T]} A_s(u^0)\le C_{s,h,T} A_s(\text{initial data}),
                                                               \tag{13}
\]

with the usual tame high bound when the low norm is controlled. The seed
is not asserted to solve (1) after the time cutoff. A polynomial
\(u_0+tu_1\) is not used as an \(H^m\) extension of \(H^{m-1}\) velocity.

For \(n\ge0\) solve, on the whole fixed interval,

\[
 g^{\mu\nu}(x,u^n)D_\mu D_\nu u^{n+1}
       =F(x,u^n,Du^n),                                \tag{14}
\]

with the same exact initial data and homogeneous mixed traces. At each
stage the coefficients satisfy (3), because \(u_D^n=0\). Differentiating
(14) at \(t=0\) shows inductively that all initial jets through order \(m\)
are the true nonlinear jets: computing the next current jet uses only
earlier current jets and preceding-iterate jets already matched by the
seed. These jets are zero in an open shell collar. Thus the smooth linear
construction in section 2 applies at every stage. This is stronger than
merely observing that the two uncommuted boundary traces vanish.

At Sobolev order \(m\), the required corner conditions include
\((J_j)_D|_b=0\) for \(j\le m-1\), and
\(D_z(J_j)_N|_b=0\) for \(j\le m-2\); these are exactly the integer orders
at which the corresponding traces are defined. Here all finite corner
identities inherited from the open zero collar are imposed.

## 5. Actual tame estimates; avoiding false mixed-time interpolation

There is a useful reason to use the spatial state norm (4). At one fixed
time, mixed time jets of arbitrary trial functions need not obey a tame
interpolation inequality. Instead, use their equations.

Set \(Y_n=(u^n,u_t^n,D_xu^n)\). Solving (14) for \(u_{tt}^{n+1}\) gives a
first-order evolutionary identity

\[
 \partial_tY_{n+1}
 =A^i(x,u^n)D_iY_{n+1}
       +H(x,Y_n,Y_{n+1}),                             \tag{15}
\]

where \(H\) is smooth and is affine in the current field. Spatial/time
connection terms are harmless smooth terms of the indicated orders.
Repeatedly use (15) to remove time derivatives. A spacetime derivative
of order \(r\ge1\) of \(u^n\) becomes a finite sum of spatial differential
polynomials of total derivative order at most (r-1) in recent \(Y\)'s.
Only finitely many indices, between \(n-r\) and \(n\), occur. When this
recursion reaches the seed, retain the protected Cauchy field
\(Y_{\mathrm{pre}}\) as an additional **known** field on the support of the
fixed cutoff and its derivatives. Its spatial high norm is bounded by
\(C_sD_s\). Cutoff derivatives are fixed smooth coefficients. One must not
bound higher time jets of the seed by its instantaneous \(A_s\) norm:
ratios of cutoff derivatives can be unbounded near a flat cutoff endpoint.
Thus seed terms contribute known \(C_sD_s\) terms to the estimates below,
rather than another instantaneous iterate norm.

The ordinary four-dimensional spatial product inequality is sufficient:
if \(\sum k_j\le q\), generalized Hölder and Gagliardo–Nirenberg give

\[
 \Big\|\prod_jD^{k_j}Z_j\Big\|_2
 \le C_q\sum_j\|Z_j\|_{H^q}
                 \prod_{\ell\ne j}\|Z_\ell\|_\infty,
                                                               \tag{16}
\]

with the usual harmless low terms on the compact cylinder. Composition
has the analogous Moser bound. Formula (16) follows, for derivative orders
summing to \(q\), by taking Hölder exponents \(2q/k_j\) and interpolating
each factor between \(L^\infty\) and \(H^q\); Young's inequality makes the
result linear in the high norms. A fixed low bound \(A_6\le R\) supplies
all required \(L^\infty\) norms, including \(D^2u\), by Sobolev embedding
and (15).

The hypothesis that \(g=g(x,u)\) has no derivative dependence is used
here in a decisive way. In (11) a coefficient derivative of order
\(|\beta|\ge1\) contributes spatial derivative order at most
\(|\beta|-1\) in \(Y\). The remaining differentiated \(D^2u\) contributes
at most \(|I|-|\beta|+1\). Their sum is at most \(|I|\le m-1\).
Similarly \(\mathcal T^{m-1}F(Y_n)\) has total order (m-1) in \(Y\).
Consequently, for a set \(I_n\) of at most (m+3) recent indices,

\[
 e_m(u^{n+1};t)
 \le C_{m,R}D_m
       +C_{m,R}\int_0^t\sum_{\nu\in I_n}A_m(u^\nu;s)\,ds,
                                                               \tag{17}
\]

where \(D_m\) is the initial \(H^m\times H^{m-1}\) size. The known protected
field's time integral has been absorbed into \(C_{m,R}D_m\), with a constant
depending on the fixed slab and cutoff. This is linear in
the high norms. Constants depend on the fixed background and low bounds,
not on high norms or cutoff.

Pointwise normal recovery costs one less order of forcing. Every forcing
derivative there has order at most (m-2), and every coefficient-derivative
term has total order at most (m-2) in \(Y\). Therefore the actual estimate
is

\[
 A_m(u^{n+1};t)
 \le C_{m,R}e_m(u^{n+1};t)
         +C_{m,R}\sum_{\nu\in I_n}A_{m-1}(u^\nu;t)
         +C_{m,R}D_{m-1}.                            \tag{18}
\]

In particular the pointwise \(F(u^n,Du^n)\) terms have been accounted for;
they have not been put under a time integral without justification.
Use the purely spatial interpolation

\[
 A_{m-1}(v)\le\eta A_m(v)+C_{m,\eta}A_1(v).          \tag{19}
\]

Choose \(\eta\) small enough to absorb the current high norm in (18)
and make the sum of preceding high coefficients less than, for example,
(1/4). Thus (17)–(19) are a genuine closed tame hierarchy, rather than
an assertion that every high coefficient norm is bounded by a low one.

## 6. Uniform prescribed lifespan and convergence

The low energy has the independent estimate

\[
 A_1(u^{n+1};t)
 \le C D_1+C_R\int_0^t
       [A_1(u^{n+1};s)+A_1(u^n;s)]\,ds.              \tag{20}
\]

Here \(F(x,0,0)=0\) gives
\(\|F(u^n,Du^n)\|_2\le C_R A_1(u^n)\).
It follows by induction that \(A_1(u^n;t)\le K_1D_m e^{K_1t}\), with
\(K_1\) fixed. Combine (17)–(19) with (20). An envelope

\[
 A_m(u^n;t)\le K_mD_m e^{K_mt}                       \tag{21}
\]

closes: the algebraic previous-iterate term consumes at most one quarter
of the envelope, and choosing the exponential rate sufficiently large
makes its Volterra integral consume at most another quarter. Increase the
initial envelope constant to dominate the seed and initial jets.
For \(m=6\), use a first-exit bootstrap \(A_6\le R\) to justify the low
constants; choose \(D_m\) so that the right side of (21) is less than
(R/2) for all \(t\le T\). It then prevents the first exit. Higher finite
orders obey the same tame estimate with low \(A_6\) controlled; no extra
smallness of their high norms is necessary.

This argument runs the linear solves from the original initial section
over the whole prescribed slab. It requires no restart at a nonzero shell
corner and no unproved Sobolev jet-extension assertion there. It proves
the required fixed finite lifespan directly. A global or large-data
continuation theorem is not needed for this conclusion and is not claimed.

For differences \(w^{n+1}=u^{n+1}-u^n\), the right side of the frozen
equation is

\[
 F(u^n,Du^n)-F(u^{n-1},Du^{n-1})
 +[g(u^{n-1})-g(u^n)]D^2u^n.                         \tag{22}
\]

Its \(L^2\) norm is at most \(C_R A_1(w^n)\), since \(D^2u^n\in L^\infty\)
and the principal coefficients contain no \(Du\). The initial difference
is zero and the mixed traces are homogeneous. Equation (7) gives

\[
 A_1(w^{n+1};t)\le C_R\int_0^t A_1(w^n;s)\,ds.      \tag{23}
\]

After absorbing the ordinary energy exponential, an exponentially weighted
supremum norm makes (23) a contraction on the entire fixed interval.
Alternatively iterate (23) to obtain factorial convergence. The high bound
(21) and compact interpolation give a limit strongly below the top order
and weakly at the top. It solves (1), has the correct initial data and
traces, and is unique by the same estimate (22) for two actual solutions.

Localized energy gives finite propagation for the local limiting equation.
Choose a fixed radial speed bound exceeding every characteristic radial
speed allowed by (2) and the state neighborhood. If the initial support
ends at \(L\), choose \(Z>L+c_*T\), also enclosing the fixed transformation
collar. The outer artificial wall is outside the actual perturbed domain
of dependence. Angularly projected intermediate equations need not have
finite angular propagation; that property is used only after passage to
the local equation.

## 7. Difference estimates needed for top continuity

For two smooth actual solutions \(u,v\), uniformly bounded in \(A_m\),
with the same fixed support gap and protected seed interval as above,
their difference equation is again (22). The same spatial jet substitution,
commutation, and normal recovery give the two estimates

\[
 \sup_{t\le T}A_{m-1}(u-v;t)
 \le C_{T,M}A_{m-1}(u-v;0),                          \tag{24}
\]

\[
 \sup_{t\le T}A_m(u-v;t)
 \le C_{T,M}\{A_m(u-v;0)
    +[1+\sup_{t\le T}A_{m+1}(v;t)]A_{m-1}(u-v;0)\}.  \tag{25}
\]

To see the only derivative loss in (25), commute (m-1) times through
\([g(v)-g(u)]D^2v\). Its highest term contains \(D^{m+1}v\) multiplied
by (u-v). Bound that factor in \(L^\infty\) by \(A_{m-1}(u-v)\), while
all remaining terms are bounded by \(C_M A_m(u-v)\), using (16). At order
(m-1), only \(D^mv\) occurs, so (24) closes with the already controlled
high \(A_m\) bound. Terms from the smooth first-derivative nonlinearity
obey the same tame product estimate. In normal recovery the forcing is
differentiated one fewer time, so it introduces no additional loss.
Initial time-jet differences satisfy these same bounds by their triangular
recursion. This proves (24)–(25), including their uniform dependence on the
single displayed extra derivative of \(v\).

## 8. Finite regularity and strong top time continuity

Choose smoothing operators \(S_N\) by extending the shell-flat data to a
boundaryless manifold, smoothing there, then multiplying by a fixed cutoff
which equals one on their support and vanishes on a smaller common shell
collar. The two Sobolev components are smoothed at their own orders. They
satisfy

\[
 \|S_Nd\|_{A_{m+1}}\le C N\|d\|_{A_m},\qquad
 \|(I-S_N)d\|_{A_m}\to0,\qquad
 N\|(I-S_N)d\|_{A_{m-1}}\to0.                       \tag{26}
\]

The last assertion follows by splitting off a fixed smooth approximation
and estimating the remaining high-frequency tail. Multiplication by the
fixed support cutoff preserves these properties. No auxiliary mixed
Laplacian power domains, and hence no unintended higher boundary identities,
are introduced. For each smoothed datum regenerate its actual nonlinear
precontact seed; thus every approximation has its own exact nonlinear
initial jets and the same open shell collar.

Let \(u_N\) be the corresponding smooth nonlinear solution. The high tame
bound gives \(\sup_tA_{m+1}(u_N)\le C N\|d\|_{A_m}\), while the \(A_m\)
and low bounds are uniform. Apply (25) with \(v=u_N,u=u_M\), \(M\ge N\).
Equations (24) and (26) show that the right side tends to zero as
\(N\to\infty\), uniformly in \(M\ge N\). Thus \(u_N\) is Cauchy in
\(C_tH^m\cap C_t^1H^{m-1}\). This supplies strong top time continuity;
it is not inferred from weak compactness and interpolation. The actual
nonlinear jet recursion then gives convergence and continuity of all
time jets in (5). The same three-term smoothing comparison proves
continuous dependence in the top norm. Equation (24) supplies the explicit
one-derivative-lower Lipschitz estimate.

The smoother approximations need not satisfy the physical constraint
equations: they are an analytic device for the reduced system. The limit
has exactly the original completed data, so the exact constraint and gauge
properties are recovered from those data and subsidiary uniqueness, not
from any assertion about the approximants.

## 9. What this closes and what remains separate

The analytic argument closes the smooth all-angular frozen solver, the
cutoff-independent tame mixed estimates, inherited finite corner jets,
compatible whole-interval Picard iteration, low-norm contraction, prescribed
finite lifespan for sufficiently small data, strong top Sobolev continuity,
and continuous initial-data dependence. It deliberately avoids a restart
argument requiring new nonzero-corner extension theory, and avoids a claim
of top-norm smooth dependence on parameters.

To apply it to the selected aligned quadrupole \(A=B=\mathrm{diag}(1,-1,0,0)\),
one must supply the exact initial harmonic completion and prove the exact
subsidiary boundary system in the actual variables. The prior fixed-finite
\(H^{14}\times H^{13}\) canonical family, its common open shell collar, and
the exact nonlinear boundary flow are the retained inputs. Parameter
variations need their own explicitly reserved derivative budget. The
physical first signal arrival is governed by the actual perturbed metric;
it need not equal the background reference time. No Big Bang, novelty,
spectral stability, or publication claim follows from this lemma.
