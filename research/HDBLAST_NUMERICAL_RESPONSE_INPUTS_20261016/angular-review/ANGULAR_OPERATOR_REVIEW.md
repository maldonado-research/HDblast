# Independent homogeneous angular-operator review

16 October 2026. The complete homogeneous 23-amplitude bulk operator, all
boundary rows, primitive angular maps and physical Gram matrix are accepted
within the stated conditional model. The sharp selected-space evaluation
constant is reviewed separately in `ANGULAR_EVALUATION_REVIEW.md`.
**The complete quadratic Einstein source is not assembled by this result.**
The current producer expressly leaves its Ricci and reduced-gauge derivative
adapter open. No evaluated response or continuous approximation residual is
certified here.

## Intrinsic derivatives and exact angular closure

The producer represents sphere tensors in ambient tangent slots. Its rule

    (D_i T)_J = P_i^a ∂_a T_J + sum_s n_(j_s) T_(J[j_s→i])

is correct for tangent covariant tensors. Differentiating n^j T_j…=0 shows
that projection of each tensor slot contributes exactly its displayed plus
term. The derivative slot is also projected; applying the rule a second
time acts on that new slot as well. Thus a missing derivative-slot connection
has not been hidden in normal-frame notation.

Replacing an ambient polynomial by another representative modulo |n|²−1
does not change this derivative on the sphere: the extra derivative is
proportional to Pn=0. The executable reduction eliminating n4² is therefore
compatible with the projected derivative. Exact zero reconstruction in this
quotient is an all-point identity, not a comparison at angular nodes.

For λ=l(l+2), direct sphere commutation gives

    div Df=−λf,                    Δ Df=−(λ−2)Df,
    Hess f=HessTF f−(λ/3)fκ,       div(HessTF f)=(2−2λ/3)Df,
    Δ(HessTF f)=−(λ−6)HessTF f,
    div(fκ)=Df,                    tr(fκ)=3f.

The last trace-free Hessian eigenvalue follows either by commuting twice
on the constant-curvature sphere or by its spin-two homogeneous-bundle
Casimir. The two chiral TT features have trace/divergence zero and rough
Laplacian eigenvalue −6. These identities determine all nine primitive
matrices independently. The audit compares them against the candidate's
exact polynomial projections and zero reconstruction checks.

The coefficient order and Gram matrix agree with the selected features:
five scalar copies of 1,Z,F; two vector copies of DZ,DF; one tensor copy
HZ,HF,T+,T−. k_y0 and both mixed vector sectors have physical weight two,
and Jκ has weight three. Every TT slot remains present. The boundary split
is five Dirichlet amplitudes (L0,L2,L4,v2,v4) and eighteen Robin amplitudes.

## Independent frame and curvature derivation

Use K=k_yy,L=k_y0,N=k_00, m_ab=Jκ_ab+H_ab, v_a=k_ya,u_a=k_0a.
The covariant connection acting on k has a minus sign in each covariant
slot. The time action divided by a is consequently

    O(K,L,N,J,v,u,H)=(-2L,-K-N,-2L,0,-u,-v,0).

For the angular cross term, apply the connection to D_a k before summing a.
Its unprojected actions are

    Ca: K=−2div(v), L=−div(u), v=DK−div(m), u=DL,
        m=Dv+(Dv)^T;
    Ch: L=−div(v), N=−2div(u), v=−DL, u=−DN−div(m),
        m=−Du−(Du)^T.

Substitution of the sphere identities above yields every Ca and Ch table
entry, including the β=2−2λ/3 signs and the trace coefficients ±2λ/3.

The connection squares were derived separately before comparison with the
producer. For example, with general a,h,

    sum Ω_a² K = −6a²K+6a²J+6ahL,
    sum Ω_a² L = 3(h²−a²)L+3ah(−K+N+2J),
    sum Ω_a² v = (−6a²+h²)v+5ah u,
    sum Ω_a² u = (−a²+6h²)u−5ah v.

The angular trace output is
2a²K+2h²N+2(h²−a²)J−4ahL; the N output is
6h²N+6h²J−6ahL. Every angular trace-free tensor, including either TT
feature, has coefficient 2(h²−a²). These yield exactly Sa,Sah,Sh.
The rough frame wave operator contributes

    a²(Sa−O²)+ah(Sah−3O)+h²Sh.

No frame derivative of a occurs in this zeroth-order term: the radial frame
is parallel and the time derivative of a is zero. Curvature coefficients
are handled separately by the reduced Ricci operator.

With b=−(a_y+a²), d=ρ^−2−a², the warped Riemann coefficients give raised
Ricci eigenvalues 4b in the normal direction and b+3d in tangent directions.
For
L_R(k)=(Ric·k+k·Ric)/2−Riem(k), independent contraction gives

    L_R K = b(4K+N−3J),
    L_R L = (7b+3d)L/2,       L_R v = (7b+3d)v/2,
    L_R N = bK+(b+3d)N+3dJ,
    L_R J = −bK+dN+(b+d)J,
    L_R u = (b+4d)u,
    L_R H = (b+4d)H,          L_R T±=(b+4d)T±.

The constant-curvature check b=d gives L_R(k)=5bk−b tr(k)g, and
L_R(g)=0. These fix both the mixed curvature coefficient and its sign.
Solving the reduced wave equation for U_xx introduces the displayed **−2**
curvature coefficient in D; no convention reversal is being made.

## Scalar gauge coupling and boundary signs

The scalar equation contains the full dC[k] coupling. In this frame

    dC[k]_y=(K_y+N_y−3J_y)/2−L_t/ρ+div(v)/(ρc)
              +a(4K+N−3J)−3hL.

Moreover k:Hess(phi0)=p_y K−ap N+3ap J. Before solving for the unit wave
principal part, their sum is

    (p_y+4ap)K−3phL

plus the displayed derivatives. Thus N and J cancel only after retaining
dC[k]. Under y-derivative conversion ∂y=−(ellρ)^−1∂r and time conversion
∂t=ell^−1∂x, the three derivative matrices are exactly

    Aphi: L←2w, w←L;
    Bphi: K←−4w, w←(K+N−3J)/2;
    Ephi: v←2w, w←λv.

The first two have coefficient ellρp; the angular one has coefficient
ell²ρp/c. The remaining scalar value row is
−ell²ρ²(p_y+4ap)K+3ell²ρ²phL. This independently verifies the cancellation
in `linear_operator.assemble()`.

The parallel-frame boundary rows are

    N_y=−aK−sigma' w/3,        J_y=aK+sigma' w/3,
    w_y=pK/2−sigma''w/2,
    K_y=−4aK−2aN+6aJ+(4sigma'/3)w,

with u,H,T± homogeneous Neumann and L,v Dirichlet. The conversion to r
multiplies By by −ellρ. In particular a coordinate-component 2aJ term
does not belong to this parallel-frame matrix. The audited interface also
packs the boundary sources with metric coefficient +2ellρ, scalar coefficient
+ellρ, normal-gauge coefficient −ellρ, and zero Dirichlet sources.

## Exact comparison strategy and scope

`check_independent_operator.py` contains independently derived rational
tables, rather than importing the candidate checker. Candidate modules are
loaded only as the implementation under review. Every entry of all 25
primitive 23×23 matrices is compared with the independent tables. The nine
intrinsic sphere primitives are compared with their closed-form commutation
matrices, after the candidate has reconstructed each result exactly.

For full A,B,D,M,PD,PN comparison, the reviewer implements a separate sparse
Laurent-polynomial class in 15 independent background symbols. It constructs
the operator from the independent tables and, in the scalar row, adds the
uncancelled Hessian and gauge contributions separately. Every one of 3174
full matrix entries is compared identically as a Laurent polynomial. This is
not a collection of generic background evaluations. Treating h,Ht,c,b,d as
independent symbols makes the comparison stronger than checking only the
physical relations between them. The physical definitions are imposed when
using the matrices.

The field/source/initial packing interfaces are checked in every coordinate,
including both TT positions. The bulk forcing coefficients are −2ell²ρ²
for the metric and +ell²ρ² for the scalar, consistent with factorial second
variations. Initial ordinary t velocities receive exactly one factor ell.
An incomplete Einstein-source dictionary is rejected. These interface checks
do not calculate the missing full quadratic Einstein source or validate
the producer's source-adapter formulas themselves.

Normal and optimized reports agree. Seven wrong-formula controls per mode
reject a time-connection sign, tensor rough-Laplacian shift, coordinate/frame
Robin confusion, angular connection sign, mixed Frobenius weight, curvature
sign, or omitted scalar gauge coupling. Algebraic reconstruction and formal
coefficient identities have all-point scope; future numerical interpolation,
source error enclosures and continuous residual certification do not follow
from them.

No blocker was found for the homogeneous operator or selected-space
evaluation theorem. The provisional source inventory must remain explicit
about the unfinished quadratic Ricci/reduced-gauge adapter. All predecessor
publication files were left unchanged.

## Selected actual solution and the refined pointwise consequence

The separately pinned `SELECTED_SPACE_CLOSURE.md` is accepted as an all-point
logical statement. Sym²V is multiplicity-free and the five projected vectors
of the aligned v⊙v are nonzero. An equivariant quadratic map therefore sends
that state into one fixed feature direction in each irreducible factor, with
the allowed multiplicity amplitudes. The selected initial construction adds
no independent free angular datum. Every homogeneous bulk/boundary map and
the half-line integral operators preserve these directions. Hence the actual
fixed-point data lie in E_v and its operator preserves E_v. Closedness and
uniqueness place the actual accepted formal solution in E_v; its existing
scaled x/r derivatives remain there as well.

This proof needs the naturality and selected-data hypotheses. Stabilizer
invariance alone is insufficient. It does not need numerical values of the
unfinished Ricci source adapter, since naturality already controls the
source's angular representation. That distinction is logically sound.

Consequently the sharp selected-space constant is valid for this actual
formal response and its scaled derivatives. Conditional on the independently
reviewed response norm bound ||U||_X<1.03e33, rational arithmetic gives
sqrt(321/19)<4.111 and hence

    sup |U|_F < 3·4.111·1.03e33 = 1.270299e34,
    sup |U_x|_F, sup |U_r|_F < 4.111·1.03e33 = 4.23433e33.

The factor three is the value weight in X. These are upper bounds for the
selected formal coefficient, not evaluated waveform values. They control
frame components and their scaled partial x/r derivatives, not unspecified
higher covariant derivatives. The sharper factor applies to an approximation
error U−V only if its reconstructed candidate V is also in E_v. Arbitrary
full-space errors retain the factor17 certificate. No equality of the two
TT amplitudes, orientation-reversing reduction, or 22-space assumption is
used here.
