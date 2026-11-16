# Explicit retained Einstein–scalar coefficient envelope, run 06g

6 November 2026. This conditional mathematical checkpoint gives coefficient
inputs for the retained evolution problem; it is not an evolved estimate.
It gives an exact finite formula for every one of the sixteen reduced PDE
rows after the literal proper-y, background-subtracted flow and moving-frame
reconstruction. On the specified real inner patch the coefficients of its
quadratic first-jet polynomial, including mixed base/state partial derivatives
through total order sixteen, have the conservative numerical ceiling
**10^150000**. The size is deliberately uneconomical. It is conditional on the
retained actual smooth background and its authenticated all-point collar box;
it is not conditional on an unspecified higher background derivative budget.
A first-jet bound, Sobolev product constants, and an evolved bootstrap remain
open. The separate outer complex flow domain does not imply this physical patch.

## Source pin and normalization, before reduction

`SOURCE_PINS.json` records actual byte hashes and the immutable public origins
of all inputs used. Specifically the 06f copies of
`INITIAL_GAUGE_SUBSIDIARY_AND_PARAMETER.md`, lines 30 and 46–57, pin the action
normalization **R−|dφ|²−2U**, the coupling and the relative harmonic reduction:

    E_ab = Ric_ab − φ_a φ_b − (2/3) U(φ) g_ab = ∇_(a C_b),
    C^a = g^{bc}(Γ[g]^a_bc − Γ[g0]^a_bc),   □_g φ = U′(φ).

Parentheses include one half. The metric signature is (−++++), with the
negative slot t rather than the first coordinate in the ordering below.
`MODEL_AND_ARCHIVE_CONTEXT.md` pins

    η=φ−1, W=1/3+η²+η³/3, U=(W′)²/2−2W²/3
     =−2/27+(14/9)η²+(50/27)η³−η⁴/6−(4/9)η⁵−(2/27)η⁶.

The actual fixed member has δ=1/1000 and
c=96714321425692/161839258930731 in its shell tension. These shell parameters
enter the retained flow F_z; they do not replace the bulk potential.
The **full** reference is g0=ρ(z)²[dz²−dt²+cosh²(t)κ_S³], φ0=1+η0(z).
Its Christoffel symbols and their derivatives below are functions of position.
Only the auxiliary boundary generator freezes shell coefficients.

## Exact local reconstruction

Use the two usual stereographic charts of S³, restricted to |ξ|≤1. Their
closed hemispheres cover the sphere. Put s=2/(1+|ξ|²), θ^i=s dξ^i, and
E=diag(1,1,s cosh t,s cosh t,s cosh t), in coordinate order
(z,t,ξ1,ξ2,ξ3). Thus the coframe is (dz,dt,cosh(t)θ^i). Any original round
orthonormal frame is related by a position-dependent orthogonal rotation;
the retained generator commutes with these rotations, since its tangential
operations are matrix inverse, multiplication, and trace against J4.

Write u=F_z(V)=(v,Htilde,w,D), n=1+v. Reconstruct frame components by

    G_zz=ρ² n²,
    G_ab=ρ² J4_ab+ρ_b² Htilde_ab,
    G_za=−ρ ρ_b D_a,            φ=φ0(z)+w,
    g_ij=Σ_AB E^A_i E^B_j G_AB.

This defines the sixteen-component map X=Ψ(x,V)=(g_ij,φ), with the metric
stored in its fifteen independent symmetric entries. All x derivatives of
Ψ are taken at fixed V. In particular they differentiate ρ, φ0, E, and F_z.
The conformal coordinate satisfies dy=−ρ dz; its sign and factor are already
present in G_za. The n here is the proper-y component root, not a temporal
ADM lapse or an interior Schur-complement lapse.

The physical region is 0≤z≤6·10^-6, |t|≤13·10^-6/8, ||V||star≤10^-3,
with star=max(|v|,||Htilde||op,|w|,||D||2). The retained flow gives
||u||star<.00101 and its reconstructed frame metric obeys
||ρ^-2G−J5||op<1/40. The 06f inverse metric bound is <1/16000.
The outer .01 state box is used solely for analytic flow estimates.

## Exact lower-order formula without omitted Ricci terms

All indices in this paragraph range over the five coordinate directions.
Let q=g^-1, k_iab=∂_i g_ab, p_i=∂_i φ, and K^a_bc=Γ[g0]^a_bc. Define

    A^a_bc = (1/2) q^{ar}(k_brc+k_crb−k_rbc),
    q_i^{ab} = −q^{ar} k_irs q^{sb},
    T_i^a_bc = −q^{ar} k_irs A^s_bc,
    C^a = q^{bc}(A^a_bc−K^a_bc),
    C_i^a = q_i^{bc}(A^a_bc−K^a_bc)
              +q^{bc}(T_i^a_bc−∂_i K^a_bc),
    C_b = g_ba C^a,
    N_ib = k_iba C^a+g_ba C_i^a−A^a_ib C_a,
    R0_ab = T_c^c_ab−T_b^c_ac
                  +A^c_cr A^r_ab−A^c_br A^r_ac,
    L_ab = R0_ab−(N_ab+N_ba)/2−p_a p_b−(2/3)U(φ)g_ab,
    Q_ab = 2 L_ab,              Q_φ = q^{ij} A^a_ij p_a+U′(φ).

Repeated indices are summed. T is precisely the derivative of Γ evaluated
with all metric second derivatives set to zero. This is an explicit rational
finite expression, not a symbolic placeholder for an unspecified lower term.
With a general symmetric second jet h_ijab, differentiating Γ and C adds
exactly −(1/2)q^{ij}h_ijab to Ric_ab−∇_(a C_b). Indeed substitution of
(1/2)q^{ar}(h_ibrc+h_icrb−h_irbc) cancels all second derivatives except
that contraction. Hence the original reduced equations are exactly

    q^{ij} ∂_i∂_j X = Q(x,X,∂X).

The scalar sign is fixed by □φ=q^{ij}(φ_ij−Γ^a_ij φ_a). Both ∂K and K
are retained. The exact finite checker tests fifteen nontrivial rational
metric rows and rejects freezing ∂K, but the cancellation above is an
algebraic identity for arbitrary jets, not a theorem inferred from samples.

## Literal point-map chain rule and quadratic structure

Put J=Ψ_V, P_i^A=∂_i V^A, k/p=Ψ_i+J P_i in the preceding formula. Then

    q^{ij} ∂_i∂_j V = f(x,V,P),
    f = J^-1 { Q(x,Ψ,Ψ_x+JP)
        −q^{ij}[Ψ_ij+2 Ψ_iV P_j+Ψ_VV(P_i,P_j)] }.                 (1)

This determines all sixteen rows and all their base/state derivatives by
finite differentiation. The map is invertible on its image, as required;
no assertion that a ball is mapped onto itself is used.
Q is of degree at most two in the first jet, so uniquely

    f^A = c^A + b^{A,i}_B P_i^B
                  + d^{A,ij}_{BC} P_i^B P_j^C,                  (2)

where d is chosen symmetric in the combined indices (i,B),(j,C).
For an explicit extraction: c=f(P=0), b=∂_P f(0), d=(1/2)∂²_P f(0).
Equation (1) and the displayed finite formulas are an evaluation algorithm
for the actual c,b,d, not merely a general quasilinear structural statement.
No derivative of V of order two enters (2). Also f(x,0,0)=0, because
F_z(0)=0, Ψ(x,0)=(g0,φ0), and that actual background solves the reduced system.

If the principal operator is expressed with the fixed product bundle
connection D_i=∂_i+Ω_i and product base connection γ^k_ij, add to f the term

    q^{ij}{Ω_i P_j+Ω_j P_i−γ^k_ij P_k
       +(∂_iΩ_j+Ω_iΩ_j−γ^k_ijΩ_k)V}.                            (3)

Here Ω is the scalar/vector/symmetric-tensor induced round-sphere connection,
trivial in z,t; γ is the product coordinate connection for dz²+dt²+κ.
This follows by expanding D_iD_j including its base covector index. The
coefficients are explicit from s: ∂_i log s=−2ξ_i/(1+|ξ|²),
γ^k_ij=δ^k_j∂_i log s+δ^k_i∂_j log s−δ_ij∂^k log s on angular indices,
and ω_i^a_b=δ^a_i∂_b log s−δ_ib∂^a log s for the orthonormal frame.
The usual induced tensor representation gives Ω. Its choice of component
sign convention must be used consistently in (3). This changes only c,b.
It does not freeze or replace the physical bulk K.

## Numerical all-point derivative inputs

The 06e source-jet note gives the actual coupled background equations and
all-point box. Converting its ODE literally to z, introduce β=1/ρ:

    ρ_z=−aρ², β_z=a, a_z=β+ρp²/3,
    p_z=−ρU′(η)+4ρap, η_z=−ρp.

At any center use absolute zeroth bounds (131,1/129,.112,.00014,.00009).
For n=0,...,18, ordinary Cauchy products of the currently known nonnegative
Taylor coefficients bound the right sides; use absolute polynomial
coefficients, divide order n by n+1 and round upward to the next integer.
Induction on n bounds every normalized derivative at every center.
The closed replay described in [REPRODUCE.md](REPRODUCE.md) evaluates this recurrence with Fraction and
verifies |∂_z^n(ρ,β,a,p,η)|<10^100 for n≤19. The centers are arbitrary
throughout the retained collar: this is neither a Taylor-convergence claim
nor a numerical trajectory. The sphere/time frame and reference below need
no continuation outside that real collar.

For the angular/time frame, take the complex coordinate polydisc radius
1/100 at each real chart center. The change of 1+Σξ_i² is at most
6/100+3/10000<1/10. The cosh argument has modulus <.011;
|cosh t|<2 and |cosh t|>1/2 follow from its absolutely convergent series.
Every entry of E and E^-1 is therefore bounded by 8 there. Cauchy's formula
bounds their ordinary partial derivatives through order nineteen by
8·19!·100^19<10^60. These assertions also cover s without cosh.

The retained 06f flow supplies mixed z/state derivatives through total order
18 below 10^126 in star operator norm. Coordinate unit vectors have star
norm at most one, so each independent component partial has this same bound.
A reconstruction entry is a sum of products of at most six factors drawn
from ρ, ρ_b, E, 1+v, Htilde, D, and φ0; the worst product is
ρ·ρ·n·n·E·E. The total-order Leibniz sum is at most 6^18 times the
product of the displayed bounds; it is <10^700, even allowing the finitely
many component sums. Thus every Ψ partial through eighteen is <10^700.
For g0 and its explicit inverse β² E^-1 J5 E^-T, the four-factor bound
4^19·10^(2·100+2·60)<10^350 holds. Differentiating their Christoffel
formula through seventeen therefore gives every K partial <10^720:
125·3·2^17·10^700<10^720. This includes the genuine spatially varying
reference derivatives ∂K required in Q. Product connections in (3) have
strictly smaller bounds by the same rational-frame calculation.

## Norm conversion and finite derivative algebra

In the sixteen independent coordinates let |.|∞ be coordinate maximum.
Then |u|∞≤||u||star≤4|u|∞. Also ||u||star≤sqrt(2)|u|2 and
|u|2≤sqrt(7)||u||star; these are inequalities, not equal norms. The symmetric
ten-entry tensor norm differs from the Frobenius norm by off-diagonal
multiplicity. In particular the 06f bound ||DF_-z||star<1.01 gives
||DF_-z||∞<4.04, not 1.01 in this coordinate norm.

At real |ξ|≤1, E^-1 is diagonal with entries ≤1. The inverse reconstruction
has differential factors 1/(2ρ²n), 1/ρ_b², −1/(ρρ_b), and 1 in its
respective components, each of magnitude ≤1 on the inner patch.
Therefore ||J^-1||∞<4.04 and its sum of absolute entries is <65.
Likewise the coordinate q has sum of absolute entries <25/16000<1.
These zeroth inverse estimates use actual inner-patch membership.

For a finite jet of a scalar, vector or matrix on the 21 variables (x,V),
use the sum of absolute entries of all Taylor coefficients ∂^α/α! through
order r. This is a submultiplicative algebra norm; a contracted tensor
product is bounded by the product of these entry sums. It is equally valid
to majorize coefficient by coefficient before taking the sum. There are
binom(21+18,18)<10^12 multiindices. Consequently **B=10^1000** bounds all
needed raw building-block jets: Ψ through 18, Ψ_xx through16, Ψ_V through17,
K through17, and the product connection terms with their required derivatives.

For completeness an inverse matrix jet has recursion

    Q_0=M_0^-1,
    Q_α=−Q_0 Σ_(0<β≤α) M_β Q_(α−β).

If ||Q_0||entry≤65 and the aggregate nonconstant M jet is at most B, its
aggregate inverse jet through order16 is bounded by
65 Σ_(j=0)^16(65B)^j ≤65(1+65B)^16 <10^17000.
(The middle inequality uses binomial coefficients ≥1.) Apply this to g and
to J separately, using the zeroth bounds just proved. Set M=10^17000; it
therefore bounds every required input including q and J^-1 in the jet
algebra through order16.

The explicit Q expression is a polynomial of degree at most seven in these
building blocks; the outer J^-1 makes the maximum degree eight. Here is
an entry-sum bound avoiding a claim about an unexpanded component count.
Regard k/p=Ψ_x+JP as polynomials in the eighty P variables, and sum the
absolute jet norms of their coefficients. Symmetric-entry duplication and
five copies of J are covered by 16M. The defining formulas then give

    A ≤24M², q_i ≤16M³, T ≤384M⁴, C ≤25M³,
    C_i ≤785M⁵, C_b ≤25M⁴, N ≤1785M⁶, R0 ≤1920M⁴.

For the potential, ||φ−1||≤2M and the sum of absolute U coefficients is
less than5, so U≤320M⁶; the analogous U′ estimate is ≤640M⁵.
Thus the combined Q metric/scalar polynomial coefficient norm is <10^4 M⁷.
The chain-rule subtraction before J^-1 has coefficient norm ≤36M²:
the three terms have bounds M², 10M², 25M², respectively. Product-connection
corrections (3) have bound <10^4 M⁵, including the coordinate functions V
and their state derivatives. A deliberately loose **10^30 M^8** covers
every coefficient and polarization factor with ample margin. Restoring
ordinary partial derivatives multiplies by at most16!, and

    16! · 10^30 · (10^17000)^8 < 10^150000.

This proves the advertised bound on the **sum of absolute coefficient
entries** at each ordinary mixed partial of total order≤16. It applies to
either (1) or the fixed-product version (3). In particular on |P|∞≤R,

    |∂_(x,V)^α f|∞ ≤ C(1+R)²,
    |∂_(x,V)^α ∂_P^β f|∞ ≤ 2C(1+R)²   (|β|=1 or 2),
    ∂_P^β f = 0                                  (|β|≥3),
    C=10^150000, |α|≤16.

The displayed P bounds are intentionally loose; the r=1 and r=2 bounds
can instead be linear and constant in R. This supplies coefficient
composition inputs for quantitative energy estimates. It does **not**
supply R for an actual evolving family, commute a solution sixteen times,
or bound a Sobolev norm of the solution.

## Verification and limits

The exact checker verifies the potential identity, 19 coupled background
steps with upward rational rounding, all raw and inverse-jet exponent
inequalities, fifteen dense rational reduced-Ricci principal rows, all
sixteen moving-map chain-rule rows, and three wrong-formula distinctions
(wrong Ricci sign, omitted ∂K, omitted factor two). Its tests use nonzero
reference and first jets. The general identities and bounds are proved
above; finite checks supplement rather than replace those arguments.

The actual background theorem and 06f flow theorem remain retained inputs.
The first-jet radius R, high-order initial gauge completion constants,
Sobolev/product/commutator constants, boundary lifting, and the quantitative
Picard/energy bootstrap are not evaluated here. No finite-amplitude interval,
nonlinear observable remainder, clock or area estimate, physical energy,
or evolved common domain follows from this coefficient ceiling alone.
