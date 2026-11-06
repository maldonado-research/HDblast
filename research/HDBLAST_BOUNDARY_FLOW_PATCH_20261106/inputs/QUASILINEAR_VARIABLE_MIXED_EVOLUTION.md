# Variable mixed wave evolution after exact boundary straightening

15 October 2026. Independent analytic note. The subject is the sixteen-field
reduced system after the already established pointwise boundary flow. This
note does not derive the harmonic initial completion, the subsidiary system,
or parameter differentiability; those are separate parts of the physical
contact theorem. It proves the variable-coefficient energy and smooth frozen
linear construction explicitly, gives the compatible Picard estimate, and
records precisely the finite-regularity passage needed for the intended
initial family.

## 1. Geometric and analytic hypotheses

Let M=[0,Z]×S³ and E=E_D⊕E_N be the fixed orthogonal product bundle, with
ranks 4 and 12 and a fixed positive fiber metric. Its product connection
preserves the splitting. All norms below are these auxiliary positive norms.
Use homogeneous V_D=0 and D_z V_N=0 at z=0. At the artificial outer end one
may impose Dirichlet conditions on every component, and choose Z beyond the
maximum perturbed causal support on the prescribed interval. There is no
physical second boundary in the resulting supported solution.

After changing the overall equation sign, write the frozen principal operator

    P_a u = a u_tt − 2 b^i D_i u_t − c^{ij} D_iD_j u,

where a=−g^{tt}>0, b^i=g^{ti}, c^{ij}=g^{ij}. A fixed product-connection
rewrite adds only first/zeroth order terms; these can be included in the
source or in the estimates below. Assume, uniformly,

    a ≥ λ,       c^{ij} ξ_i ξ_j ≥ λ |ξ|²,
    |a|+|b|+|c| ≤ λ⁻¹,
    b^z=0, c^{zA}=0 at z=0,  A angular.                 (1)

These are the stated positive energy margins, not merely the assertion that
dt is timelike. The scalar coefficient c^{zz} multiplies the identity on E.
For the nonlinear system (1) follows from the fixed coordinate Dirichlet
trace and the prescribed state neighborhood. Every iterate must satisfy it.

Let X_r(u;t)=Σ_{j=0}^r ||∂_t^j u(t)||_{H^{r−j}(M)}. Let A_r(t) be the same
sum of coefficient jet norms for (a,b,c), with the fixed smooth background
included separately when useful. Let

    F_r(f;t)=Σ_{j=0}^r ||∂_t^j f(t)||_{H^{r−j}}.

The nonlinear construction below uses m=12 or m=14. The linear estimates
also hold at each integer r=10,11,12,13,14, and at the one additional
smoothing order r=15. A uniform safe low index is 10 for r≤14 and 11 for
r=15. More precisely use ℓ_r=floor((r+1)/2)+3, giving respectively
8,9,9,10,10,11 for r=10,...,15. This deliberately avoids an unjustified
fixed-time interpolation among independent time jets; no all-integer-order
tame theorem is inferred from this finite derivative count.

## 2. Positive variable energy, including its exact flux

Set

    E_a[u]=½∫_M [a|u_t|²+c^{ij}⟨D_i u,D_j u⟩+|u|²].

Its coercivity constants depend only on λ and the fixed geometry. Direct
integration by parts in P_a u=f gives

    E_a'[u] = ∫_M ⟨f,u_t⟩ + ∫_M ⟨u,u_t⟩
      + ∫_M (½a_t−D_i b^i)|u_t|²
      + ½∫_M c_t^{ij}⟨D_i u,D_j u⟩
      − ∫_M (D_i c^{ij})⟨D_j u,u_t⟩
      + ∫_∂M [b^i n_i |u_t|²+c^{ij}n_i⟨D_j u,u_t⟩].       (2)

The fixed connection has bounded curvature; a local-coordinate rewrite
produces bounded lower order terms and the same flux. At z=0 the flux is

    −c^{zz}⟨u_z,u_t⟩
    =−c^{zz}(⟨u_{D,z},u_{D,t}⟩+⟨u_{N,z},u_{N,t}⟩)=0.   (3)

The two equalities use b^z=c^{zA}=0 and then the orthogonal, preserved D/N
split. A matrix principal coefficient mixing D and N would invalidate this
calculation. At the outer all-Dirichlet boundary u_t=0. Thus

    E_a[u](t)^(1/2) ≤ C exp(C∫₀ᵗ K₁)
       [E_a[u](0)^(1/2)+∫₀ᵗ ||f(s)||₂ ds],               (4)

where K₁ controls coefficient first derivatives and bounded first/zeroth
order terms. This is homogeneous boundary control; no L² Neumann boundary
input estimate is asserted.

For tangential commutation use ∂_t and the six SO(4) action generators on the
full bundle. Their lifts preserve E_D,E_N, commute with ∂_z, and span angular
derivatives modulo bounded zeroth order maps. Their commutators with the
fixed connection and with P_a have controlled differential order. They do
not commute with variable scalar coefficients, and all resulting coefficient
commutators are retained.

## 3. A genuine frozen variable-coefficient construction

Here coefficients, data and source are smooth and satisfy every corner
compatibility identity. Smoothness is used for this construction, not as an
assumption on the final exact initial family.

### 3.1 Remove the initial-jet obstruction before Galerkin projection

Let J_j be the recursively determined time jets at t=0. Compatibility says
(J_j)_D=0 and ∂_z(J_j)_N=0 at z=0, and J_j=0 at the outer Dirichlet boundary.
A smooth Borel sum

    W(t)=Σ_{j≥0} χ(t/ε_j) t^j J_j/j!

with ε_j decreasing sufficiently rapidly gives a smooth boundary-preserving
lift with ∂_t^j W(0)=J_j. At any fixed derivative order the tail is made
summable by the choice of ε_j; this is a direct construction in the Fréchet
space of smooth sections satisfying these fixed boundary rows. Set v=u−W
and r=f−P_a W (including lower terms if present). Then v has zero initial
data and r is flat at t=0 to every order. This step matters: directly
projecting u₀,u₁ into an elliptic eigenbasis generally destroys higher
corner compatibility for the projected equation.

### 3.2 A fixed elliptic basis and only L² projection bounds

Freeze the spatial divergence-form operator at the initial time,

    A_* w = −D_i(c^{ij}(0)D_j w)+w,

with the mixed boundary rows above. It is the positive self-adjoint operator
of the closed coercive form on

    H_D¹={w∈H¹: w_D|₀=0, w|_Z=0}.

Its domain is H² with those Dirichlet rows and ∂_z w_N=0: the natural
Neumann row is c^{zj}D_j w_N=c^{zz}∂_z w_N because c^{zA}|₀=0. Standard
elliptic integration by parts and boundary tangential difference quotients
prove

    ||w||_{H²} ≤ C (||A_*w||₂+||w||_{H¹}).                (5)

For clarity, the boundary part of this elliptic estimate uses the scalar
positive c^{zz} to solve for w_zz after obtaining all tangential first
differences in H¹; D and N are separate component conditions, with no
spatial interface of different boundary type. The interior estimate follows
from ordinary difference quotients. Compact H¹→L² and the variational
resolvent give a complete smooth orthonormal eigenbasis of A_*.

Let P_N be its L² projection onto the first N eigenvectors. Solve the
finite-dimensional weak equation for v_N∈ran P_N, with zero initial data.
The time mass matrix (⟨a e_j,e_k⟩) is positive definite, so the ODE exists
on the entire coefficient interval. Testing by v_N,t gives (2) without any
commutation claim for P_N and a,b,c. All time jets of v_N at zero vanish,
because r is flat. The ODE equations imply this recursively.

The following *discrete elliptic estimate* supplies the missing regularity
step. For A(t)=−D_i(c^{ij}(t)D_j)+1 and w_N∈ran P_N,

    A_*w_N=P_N A(t)w_N+P_N(A_*−A(t))w_N.

Since P_N is an L² contraction and c(t)−c(0) is small in C⁰ on a sufficiently
short interval, (5) gives, after absorbing the principal difference,

    ||w_N||_{H²} ≤ C[||P_N A(t)w_N||₂+||w_N||_{H¹}].      (6)

The difference of the first coefficient derivatives is bounded, need not
be small, and multiplies only H¹ derivatives. Constants in (6) are independent
of N. No H^s boundedness of P_N, and no pointwise normal recovery from the
projected PDE, is used.

### 3.3 Close time regularity before passing to the limit

Put v_N^{(q)}=∂_t^q v_N. Differentiate the projected equation q times.
Equation (6) and that differentiated equation give

    ||v_N^{(q)}||_{H²} ≤ C_q [||v_N^{(q+2)}||₂
          +||v_N^{(q+1)}||_{H¹}+||r^{(q)}||₂
          +Σ_{j<q}||v_N^{(j)}||_{H²}
          +Σ_{j≤q}||v_N^{(j)}||_{H¹}].                 (7)

Terms containing derivatives of a and b are bounded by the displayed time
and H¹ norms. The differentiated spatial coefficients create H² norms only
at strictly smaller q. Hence induction in q eliminates the lower H² terms.

Now apply (2) to v_N^{(p)}. Its commutator source includes D²v_N^{(p−1)},
which (7) bounds by v_N^{(p+1)} in L² and v_N^{(p)} in H¹, both in the same
p-th time energy, together with lower energies and forcing. Consequently,
for each p,

    d/dt Σ_{j=0}^p E_a[v_N^{(j)}]
       ≤ C_p Σ_{j=0}^p E_a[v_N^{(j)}]
           +C_p Σ_{j=0}^p ||r^{(j)}||₂².               (8)

The constants may involve finitely many smooth coefficient derivatives;
tameness is not needed for this existence construction. Initial energies
are zero. Gronwall gives uniform time energies and then (7) gives uniform
H² norms. This simultaneous argument resolves the apparent circularity
between time differentiation and normal-second-derivative control.

### 3.4 Limit first; full spatial regularity afterwards

Compactness gives a weak solution v, with every time derivative in H¹ and
then H², satisfying the variational equation and the natural Neumann row.
The already established time regularity makes the spatial equation

    A(t)v = r−a v_tt+2b^i D_i v_t + first/zeroth terms

an elliptic equation with an L² right-hand side. Higher spatial elliptic
regularity then applies to the *limit equation*: if all time derivatives
are H^s, this equation makes them H^{s+1}, since the spatial first derivatives
of v_t are H^{s−1}; start at s=2 and iterate, differentiating the equation
in time when needed. More explicitly, to gain H^{s+1}, all terms on the
right are H^{s−1}, and the elliptic estimate gains two derivatives. For time
derivatives, first improve j=0, then j=1, and so on: derivatives of A(t)
multiply only lower time derivatives, which have already gained the spatial
derivative at that stage. All time derivatives are available, so this
lexicographic induction in spatial order and then time order is unrestricted
for smooth input. Time continuity follows first in H¹ from the next time
 derivative, and then in every H^s from the elliptic equation. Thus u=v+W
is smooth and satisfies the strong equation and boundary rows.

Estimate (4) gives uniqueness. Repeat this short-interval construction a
finite number of times, using a new fixed elliptic operator and the smooth
compatible terminal data each time. This proves smooth frozen variable
solvability on any compact interval satisfying (1). It does not invoke a
Sommerfeld boundary theorem, angular commutation, or a finite-mode constant
whose limit has not been controlled.

## 4. Tame high-order a priori estimates

Let Y_m(u) be the square root of the sum of E_a[T^αu] over all words in the
time/angular generators with |α|≤m−1. Thus Y_m controls all derivatives of
total order at most m having at most one normal derivative. Commuting P_a,
using (2), and retaining every coefficient commutator gives

    Y_m(t) ≤ C Y_m(0)+C∫₀ᵗ [F_{m−1}(f;s)+X_m(u;s)
                                  +A_{m−1}(s) X_10(u;s)] ds.             (9)

Here and below C depends on λ, the geometry, m, and bounded low coefficient
norms; fixed smooth background coefficient norms are included in C. For
m=12 or 14 a direct finite Leibniz proof suffices. A coefficient commutator
has factors of total derivative orders q and m+1−q, with 1≤q≤m−1. Their
sum is ≤15, so one factor has order ≤7. Its L∞ norm is controlled by X_10,
since it retains at least three spatial derivatives in dimension four.
The other factor is controlled in L² by A_{m−1} or X_m. For compositions
F(x,v,Dv), two factors with total derivative order ≥8 would use at least
14 of the at most m−1≤13 differentiations. Thus at most one factor is high,
and the same finite Leibniz estimate applies. No interpolation of arbitrary
high time jets into low time jets is assumed.

For normal recovery, solve the uncommuted equation for c^{zz}u_zz. Apply
D_z^{q−2}T^α with q+|α|≤m. Every principal term has fewer normal derivatives
than D_z^qT^αu; coefficient commutators have total derivative order on u
at most m−1. Induction on q therefore gives

    X_m(u;t) ≤ C [Y_m(u;t)+F_{m−2}(f;t)+X_{m−1}(u;t)
                                  +A_{m−2}(t) X_10(u;t)].             (10)

Terms with one normal/tangential derivative in the principal operator are
included in this induction; they are not assumed to vanish inside the
collar. Only c^{zz}≥λ is needed for the algebraic divisions. No normal
commutation of the boundary conditions is performed.

Equations (9)–(10), rather than a schematic statement that “Sobolev products
close,” are the useful tame estimate. General forcing must have the displayed
time jets. A sufficient source norm is

    Σ_{j=0}^{m−1} ||∂_t^j f||_{L¹_t H^{m−1−j}}
      +Σ_{j=0}^{m−2} ||∂_t^j f(0)||_{H^{m−2−j}}.                       (11)

The second sum and the time derivative terms in the first bound the
pointwise F_{m−2} in (10) by the fundamental theorem of calculus. This is
not arbitrary L² boundary forcing.

For the linear first/second-variation estimates at r=10,11,12,13,14, the
same formulas (9)–(11) hold with m replaced by r and X_10(u) replaced by
X_{ℓ_r}(u). Coefficient low norms through order 10 suffice uniformly over
these five orders. At r=15 use ℓ_15=11 and coefficient low norms through 11.
Indeed one factor in every principal commutator has total derivative order
at most floor((r+1)/2), leaving at least three spatial derivatives in X_{ℓ_r}.
Crucially ℓ_r≤r−1 throughout this finite range, so the normal-recovery
remainder still involves only lower unknown jets and the FTC step remains
valid. This includes the r=10,...,13 linear estimates needed for parameter
variations, without silently asserting arbitrary-r fixed-time Moser bounds.

## 5. Compatible quasilinear Picard scheme

Consider P_{g(V)}V=F(x,V,DV), with F(x,0,0)=0, coefficients depending on V
and not on DV, and boundary geometry (1) whenever V_D|₀=0. Fix m=12 or14.
Assume ||V₀||_{H^m}+||V₁||_{H^{m−1}}=ε, with the initial perturbation
identically zero on fixed shell and outer collars. The exact nonlinear
initial recursion defines J_j∈H^{m−j}, 0≤j≤m, with

    I_m:=Σ_{j=0}^m ||J_j||_{H^{m−j}} ≤ C ε.                            (12)

The recursion is local, so every J_j is zero on the initial collars. For
smooth data this holds at every order. Formula (12) follows inductively
from a≥λ and the finite Sobolev product estimates; the jet map vanishes
at zero because the background-subtracted equation has F(x,0,0)=0.

Choose a boundary-preserving seed V⁰ with precisely these initial jets.
The already available short precontact Cauchy solution, multiplied by a
smooth time cutoff equal to one near zero and supported strictly within
its precontact lifetime, is a particularly simple seed. It is zero in a
shell corner neighborhood, preserves every jet, and satisfies
sup_t X_m(V⁰;t)≤C_seed ε. This seed uses a fixed positive short time; it does
not claim that the desired contact solution has already been constructed.
For a smooth data approximation one can use its reduced Cauchy solution.

Define

    P_{g(Vⁿ)} Vⁿ⁺¹ = F(x,Vⁿ,DVⁿ),
    Vⁿ⁺¹(0)=V₀,    Vⁿ⁺¹_t(0)=V₁,
    Vⁿ⁺¹_D|₀=0,   ∂_zVⁿ⁺¹_N|₀=0.                                  (13)

The frozen coefficient trace satisfies (1), because Vⁿ_D|₀=0. If Vⁿ has
jets J_j, differentiating its frozen successor equation shows recursively
that Vⁿ⁺¹ has the same jets. To compute J_r one uses only earlier successor
jets and predecessor jets of order at most r−1. Consequently the smooth
frozen theorem applies at every step. Merely assigning homogeneous boundary
traces, without this jet argument, would not establish compatibility.

### 5.1 The normal-recovery term does not prevent closure

Suppose the iterates lie in a fixed X_m ball of radius R inside the state
and margin neighborhood. Composition estimates give

    F_{m−1}(F(Vⁿ,DVⁿ)) ≤ C_R X_m(Vⁿ),
    F_{m−2}(F(Vⁿ,DVⁿ)) ≤ C_R X_{m−1}(Vⁿ).

Equations (9)–(10) therefore give, writing u=Vⁿ⁺¹ and v=Vⁿ,

    X_m(u;t) ≤ C_R I_m + C_R[X_{m−1}(u;t)+X_{m−1}(v;t)]
                    + C_R∫₀ᵗ [X_m(u;s)+X_m(v;s)] ds.                 (14)

The instantaneous lower norms are handled exactly by

    X_{m−1}(w;t) ≤ X_{m−1}(w;0)+∫₀ᵗ X_m(w;s) ds.                     (15)

Each term ∂_t^j w in H^{m−1−j} differentiates into a term contained in X_m.
Thus (12), (14), and (15) imply

    X_m(Vⁿ⁺¹;t) ≤ C_R ε + C_R∫₀ᵗ [X_m(Vⁿ⁺¹;s)+X_m(Vⁿ;s)] ds.         (16)

This resolves the pointwise source cost in normal recovery without hiding
it in an alleged small time factor. After Gronwall, an exponential majorant
C₀ ε exp(C₁t) bounds all iterates on any fixed [0,T], provided
C₀ ε exp(C₁T)<R. Constants depend on the chosen background interval, collar,
state neighborhood and finite regularity, not on angular cutoff or iteration.

### 5.2 Lower-norm convergence and uniqueness

Let Wⁿ=Vⁿ⁺¹−Vⁿ. Its initial data and mixed boundary rows are zero. Its source is

    F(Vⁿ,DVⁿ)−F(Vⁿ⁻¹,DVⁿ⁻¹)
      +(g(Vⁿ⁻¹)−g(Vⁿ))^{μν}D_μD_ν Vⁿ,

up to the fixed sign convention. The uniform high bound implies

    ||source||₂ ≤ C_R X_1(Wⁿ⁻¹).

The first derivative part of F is Lipschitz in L² under the bounded low
norms; the second displayed term is bounded by the L² coefficient difference
times the L∞ second derivative of Vⁿ. Therefore (4) yields

    sup_{s≤t}X_1(Wⁿ;s) ≤ C_R∫₀ᵗ sup_{r≤s}X_1(Wⁿ⁻¹;r) ds.             (17)

Repeated integration gives factorial convergence on the entire prescribed
finite interval. One can instead use a short-interval contraction. Spatial
interpolation first gives convergence of V and V_t one derivative below the
uniform norm. The equation, differentiated recursively in time, then gives
convergence in X_{m−1}; no interpolation of unrelated time jets is needed.
The same energy difference estimate proves uniqueness in the high class.

## 6. Finite Sobolev completion: what is justified here

The preceding argument yields a solution in

    V∈L∞([0,T];H^m)∩C([0,T];H^{m−1}),
    V_t∈L∞([0,T];H^{m−1})∩C([0,T];H^{m−2}),
    ∂_t^jV∈L∞ H^{m−j}∩C H^{m−1−j}  (0≤j≤m−1),                    (18)

and ∂_t^mV∈L∞L², by smooth approximation of the initial data preserving the
fixed zero collars, (16), (17), and weak compactness. This finite regularity
is already ample for the classical equation, traces, and contact geometry
when m=12 or14. It does not assert spatial C∞ for the finite initial family.
The equation and all required finite corner identities hold. The exact
canonical constraints need not hold for the analytic approximants: they
are conditions on the limiting exact datum, used subsequently in the exact
subsidiary argument.

For completeness, strong continuity in the top Sobolev norm and top-norm
continuous dependence require the usual additional smoothing comparison;
they do not follow merely from a bounded sequence and weak compactness.
A sufficient comparison, with u^h a spatially smoothed-data solution, is

    sup_t X_m(u^δ−u^h;t)
      ≤ C_T [initial top-jet difference
                  +sup_t X_{m−1}(u^δ−u^h;t) · sup_t X_{m+1}(u^h;t)]. (19)

For the product in (19), the required lower difference norm is X_{m−1},
not merely X₁. One places the *coarser, smoother* solution u^h in the high
second-derivative factor. The tame difference estimate gives (19) with
that norm, and X_{m+1}(u^h)≤C h⁻¹ follows from the same high energy argument. At this one additional order use low index 11 (bounded by the uniform X_m norm), so the elementary two-factor derivative count still supplies three spatial derivatives for the L∞ factor. More explicitly, first close the base-order iteration in X_m on the common prescribed slab. At order m+1, its coefficient and source product bounds have constants depending only on that already controlled low X_m bound; the high-order analogue of (16) is linear in X_{m+1} of the current and preceding iterate. Gronwall therefore bounds X_{m+1} by its possibly large initial norm C h⁻¹ on the same interval. There is no additional X_{m+1}≤R first-exit assumption: the state margins and lifespan are protected at order m, so the smoother data need not be small in the extra high norm.
For an H^m datum, standard compactly supported smoothing has initial
H^{m−1} difference o(h) and high H^m approximation error o(1). The lower
energy-and-normal-recovery estimate for differences with matching nonlinear
initial jet maps is Lipschitz at X_{m−1}, so its size is o(h). The product
in (19) tends to zero. Each fixed smoothed-data solution is already C in X_m
by the lower-regularity convergence argument at level m+1, so the Cauchy
comparison at level m is not circular. This is the Bona–Smith completion of (18) to

    ∂_t^jV∈C([0,T];H^{m−j}), 0≤j≤m.                                  (20)

This paragraph identifies all ingredients rather than claiming top-norm
Lipschitz dependence: the solution map at top norm is continuous, while a
uniform Lipschitz estimate is asserted only below the top norm. Compactly
supported initial smoothing is legitimate here specifically because the
initial data vanish on fixed open boundary collars. Arbitrary nonzero
corner data would need a separate compatible approximation construction.

## 7. Prescribed lifespan, support, and the continuation qualification

The construction in (13) is performed on the *whole prescribed background
interval*. The smallness condition C₀εe^{C₁T}<R is the required lifespan
bootstrap. It is not just a local theorem before shell arrival. The state
bound and Sobolev embedding keep the actual inverse metric inside the
chosen positive margins throughout [0,T]. In particular the boundary is
timelike and the state-space straightening flow remains defined.

Local energy on moving truncated cones proves finite propagation from (2),
with a speed bounded uniformly by the coefficient neighborhood, for example

    v_max ≥ sup (|b|+sqrt(|b|²+a λ_max(c)))/a.

If the initial support is in z≤L and Z>L+v_max T, the limit vanishes near
the artificial boundary, so the construction there has no physical effect.
The same local estimate compares with the zero solution and proves this
support statement even if a smoothing device used in an intermediate
construction was not exactly support preserving.

The *a priori* continuation estimate is also explicit: if the state and
margins stay in a compact admissible set and X_10 stays bounded on a finite
time interval, (9)–(10), the composition estimates, and (15) bound X_m by
Gronwall. This rules out high-norm blowup while those quantities stay bounded.
A full general restart theorem at an arbitrary nonzero shell corner requires
a compatible finite-Sobolev trace lift or an equivalent general linear
regularity theorem. The zero-collar smoothing argument in §6 alone does
not prove that general restart theorem. One can avoid that issue entirely
for the intended prescribed contact interval, since the whole-interval
small-data Picard construction starts from the original zero-collar data.
This note therefore proves the needed prescribed-interval reduced evolution
and the a priori continuation bound; it does **not** silently upgrade the
latter to a general arbitrary-corner continuation existence theorem.

## 8. Scope of the result and remaining interfaces

The analytic advance is a variable-principal mixed-boundary construction,
a tame finite-order estimate, exact-jet Picard closure, lower-norm convergence,
uniqueness, and small-data lifespan on a prescribed background interval.
The frozen construction does not rely on angular symmetry of the coefficients.
The finite m=14 class is available from the prior fixed k≥12 data family
(metric H^{k+2}, momentum H^{k+1}), subject to the separate harmonic completion
and the pointwise transformations preserving those classes.

To assert the physical Einstein–scalar contact family one must additionally
supply: (i) the exact harmonic initial completion and its regularity;
(ii) the exact variable subsidiary boundary identities and uniqueness;
(iii) both-copy matching and reflection invariance; and (iv) the two-parameter
variation/difference-quotient budget, if C² dependence is claimed. No such
physical statement follows merely from the energy algebra or its symbolic
checks. The general arbitrary-corner restart theorem is the additional
analytic qualification recorded above, not an obstacle to the constructed
whole prescribed interval from the zero-collar initial data.
