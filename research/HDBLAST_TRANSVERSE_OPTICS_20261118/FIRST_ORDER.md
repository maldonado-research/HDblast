# Normal Jacobi displacement and opening along the shell-tangent rays

Independent derivation, 18 November 2026. The positive shell light-pencil distance response is accompanied by a normal response of exactly minus twice its first variation. The three-dimensional bulk screen has zero first-order trace. For the registered central ray, the normal fractional displacement coefficient lies between **−1.152807 and −1.041260 times 10⁻²⁴**, and its parallel-frame opening coefficient lies between **−3.564043 and −3.219147 times 10⁻²⁴**. These numbers are derivatives at zero canonical amplitude, not finite-amplitude observations.

The calculation uses the smooth metric up to the boundary of either single bulk copy. It does not assign a geodesic crossing law to the doubled distributional shell. Its normal Jacobi field is defined from the one-sided curvature jet. An actual interpretation as the derivative of inward-tilted rays requires the inward ray variation described below. A cap here is in the shell observer's tangent **S² sky**, not a positive-measure cap in the full bulk S³ sky.

## 1. Fixed physical normalization

Use the vertex and registration of `HDBLAST_FINITE_LIGHT_CONES_20261118/DERIVATION.md`. At the observer label n*=(0,0,1,0), the vertex clock is q₀=ell a, with

    ell=10^-6,  a=49/32,  b=99/64,  Delta=b-a=1/64,
    q₁=ell b,
    L=rho_b [sinh(q₁)-sinh(q₀)]/cosh(q₀).

The null tangent has unit observed frequency at reception. The physical affine parameter λ increases into the registered past, starts at zero, and ends at the same L for every amplitude ε. The unperturbed shell is

    h₀=rho_b²[-dq²+cosh(q)² kappa].

The actual receiving sky triad is built from the transported label directions n₁,n₂,n₄ and orthonormalized in the physical metric. A direction v=(v₁,v₂,v₄) has |v|=1. Put Q(v)=v₁²−v₂². The background ray obeys

    λ(q)=rho_b [sinh(q)-sinh(q₀)]/cosh(q₀),
    sin θ=λ/[rho_b cosh(q)],
    Y(n(λ))=Q(v) sin²θ.

The proper observer time q, the affine parameter λ, and the normal angle are different quantities. The result is at fixed λ=L. Fixing an endpoint clock or redshift instead would require the corresponding endpoint correction.

## 2. Exact one-sided three-screen geometry

The actual shell and bulk equations are

    K_ab=c h_ab,  c=sigma(Phi)/6,
    Ric⁵_AB=phi_A phi_B+(2/3)V(phi)g_AB,
    Phi=phi|shell.

Here the normal is spacelike, K_ab=g(∇_a n,e_b), and the curvature convention is the one in the retained Gauss note. Since K(k,k)=0, an affinely parametrized shell null geodesic is also an affine bulk geodesic along the boundary. This is the local umbilic null-geodesic property; c need not be constant.

Choose two orthonormal shell screen vectors E₁,E₂, parallel along the ray in the shell. Their bulk normal derivative is zero because K(k,E_A)=0. Also ∇_k n=c k. Thus n is parallel as a class in k-perpendicular modulo k. The representative

    N(λ)=n(λ)-[integral_0^λ c(s) ds] k(λ)

is genuinely parallel and agrees with n at the vertex. Reversing the normal changes the representative orientation, not any of the scalar optical coefficients below. Use the inward orientation when interpreting a positive launch tilt into a single copy.

Define T_AB=R(E_A,k,E_B,k) so the screen Jacobi equation is J″=−TJ and tr T=Ric⁵(k,k). Gauss removes the extrinsic curvature products from the two shell screen entries. Codazzi makes R(n,k,E_A,k)=0: derivatives of c multiply h(k,k) or h(k,E_A), both zero. Consequently the exact optical matrix, at each amplitude, is

    T⁵ = diag(T⁴, T_N),
    T_N=(k Phi)²−Ric⁴(k,k).                                  (1)

The normal entry is also T_N=E⁵(k,k)+(k Phi)²/3, where E⁵_ab=C⁵(n,e_a,n,e_b). This keeps the Weyl projection and the nonnegative scalar contribution separate.

All background screen entries vanish: de Sitter gives T⁴₀=0, and Phi₀ is constant on the shell, so k₀Phi₀=0. With regular, observer-normalized vertex data,

    J_ε(0)=0,  J_ε′(0)=I₃,  J₀(λ)=λ I₃.

These statements require only the one-sided smooth metric and curvature along the shell. An auxiliary smooth extension can be used to construct the ray variation; its boundary Jacobi equation is independent of that extension. For sufficiently small ε, j_N(λ)/λ remains positive uniformly down to the vertex. A sufficiently small inward launch angle then gives rays remaining in the chosen smooth copy over this fixed segment. This last retention is qualitative. No numerical finite-tilt radius is inferred, and no assertion about rays crossing the distributional junction is needed.

## 3. First variation and exact trace cancellation

Write δ=∂_ε at ε=0. The accepted intrinsic metric tangent is 2Fh₀ plus a tangential Lie derivative, with F=U(q)Y(n). Its intrinsic Weyl tangent is zero. Along a background affine null ray,

    δRic⁴(k,k)=−2 F″(λ),
    δT⁴=−F″ I₂.

These expressions vary the actual null vectors and evaluating events. The vanishing background optical tensor removes their extra first-order forcing. They are not obtained by contracting an uncorrected Ricci variation against a fixed vector that ceases to be null.

Since Ric⁵(k,k)=(kPhi)² starts at order ε², equation (1) yields

    δT⁵=F″ diag(−1,−1,2),
    δJ⁵=diag(d,d,−2d),                                      (2)
    d(v,L)=L[F(0)+F(L)]−2 integral_0^L F(λ) dλ.

The Jacobi initial variations vanish. At this vertex F(0)=0. The shell angular-diameter derivative δD_A=d is precisely the archived PR154 coefficient. The normal Jacobi coefficient is therefore

    δj_N(v,L)=−2d(v,L),  δj_N′(v,L)=−2∂_L d(v,L).            (3)

The two shell screen directions expand equally while the normal direction contracts twice as much, for ε>0 and Q>0. The opposite axis or a negative ε reverses this linear strain. Equation (2) gives

    δ tr J⁵=0,
    δ det J⁵(L)=0,
    δ[(det J⁵(L))^(1/3)/L]=0.

The corresponding shell area derivative is 2Ld. Thus intrinsic two-screen defocusing does not imply five-dimensional volume defocusing. At first order the complete bulk optical tidal matrix is Weyl: its Ricci trace is zero. At quadratic order the canonical scalar null flux and the shear square are nonnegative; that is a separate calculation. None of the individual three tidal eigenvalues is required to be nonnegative by the null energy condition.

## 4. Normal displacement and an explicit angle convention

Launch a neighboring ray by a small inward normal angle α from the same vertex and receiving frame, with initial past-directed null direction −u+cos(α)e(v)+sin(α)n_in, where u is the future-directed unit receiver velocity. The α derivative at zero is the inward normal Jacobi field. Relative to the central tangent ray its signed normal displacement is

    r_N(ε,α,L)=α j_N(ε,v,L)+o(α).

The mixed derivative ∂_ε∂_α r_N|₀ is δj_N. This is a response to tilt; it is not a spontaneous normal drift of the central ray. The exactly tangent ray remains on the umbilic shell at every amplitude.

For the opening, parallel transport the receiving observer and inward screen vector along the central bulk ray. The reference ray has frequency one in that parallel observer frame. The normal projection of the covariant relative tangent is αj_N′+o(α). Equivalently, j_N′ is the infinitesimal angular opening factor in this parallel-frame comparison. Its mixed amplitude/angle derivative is δj_N′. A chosen material emitter observer would in general have a different frequency and screen; no claim that the same number is an angle measured by an unspecified emitter is made.

Set

    H(x)=M(x) sech²(ell x),
    M=mu rho_b² U/ell,
    kappa_opt=ell/(mu rho_b⁴)=ell B_b/rho_b³.

Then F(λ)=Q kappa_opt λ²H(x(λ)). Define the two distinct coefficients

    C_pos=H(b)−2 integral_0^1 t²H(x(Lt)) dt,
    C_open=H(b)+r_b H_x(b),
    r_b=L/(dλ/dx at b)=integral_a^b cosh(ell x) dx / cosh(ell b).

The exact response formulas are

    δj_N/L=−2Q kappa_opt L² C_pos,
    δj_N′ =−2Q kappa_opt L² C_open.                           (4)

In particular a terminal displacement coefficient alone does not determine the terminal opening. This calculation separately controls H_x.

## 5. A finite positive cap and its averages

Take C₊={v₁≥4/5} on the actual receiving observer's tangent S². Its solid angle is 2π/5, one tenth of that sky. At the background v₁ is uniform on [4/5,1] and E[v₂²|v₁]=(1−v₁²)/2, whence

    Q≥7/25 throughout C₊,   average_C₊ Q=18/25.

Define the scalar cap statistics using the fixed sky labels and normalized solid angle:

    J_N,cap(ε,L)=(1/Omega) integral_C₊ j_N(ε,v,L) dOmega,
    A_N,cap(ε,L)=(1/Omega) integral_C₊ j_N′(ε,v,L) dOmega.

Their background values are L and 1. Both first derivatives are 18/25 times the central v=e₁ value in (4). An RMS normal-width statistic has the same first derivative because j_N=L at ε=0, but is a different nonlinear statistic. These are averages of scalars measured in each ray's own transported frame, not vector averages across different endpoints. The two-copy geometry does not double these per-ray means.

The comparison cap v₂≥4/5 has mean Q=−18/25 and opposite linear normal response. The entire tangent-sky mean is zero. No cap has been treated as a finite three-dimensional bulk-sky solid angle, and no area or volume of an unspecified transverse aperture is assigned to these means.

## 6. Exact cell arithmetic and numerical sizes

The standalone `firstorder.py` reads only four byte-pinned archived numerical inputs in the retained sibling capsules. It imports no ancestor producer. It independently integrates the positive position kernel from PR154 and a new positive opening kernel against the continuous input cells, widened by the retained 3×10⁻⁶ transfer bound.

The endpoint equation is z″+V(x)z=u(x), with 0<V≤(5/2)ell², zero data through x=1, and H=sech^(7/2)(ell x)z. Write X=25/16, d₀=X−1, Z=||z||_infty bound, U*=||u||_infty bound, and theta=(5/4)ell²d₀². For z₀″=u with zero data, the elementary Volterra and hyperbolic-factor estimates give

    e_H=theta Z+(7/4)(ell X)² Z,
    e_Hx=(5/2)ell²d₀ Z+(7/4)(ell X)² U* d₀+(7/2)ell² X Z.

These bound |H−z₀| and |H_x−z₀′| on the full admitted interval. They use 1−sech^a(q)≤a q²/2 and |tanh(q)|≤|q|. With K=1+10⁻¹¹, Delta/K≤r_b≤Delta. Therefore the opening coefficient differs from z₀(b)+Delta z₀′(b) by at most

    e_open=e_H+Delta e_Hx+Delta(K−1)(U* d₀+e_Hx).

The flat opening kernel is simply (b−s)+Delta on 1≤s≤b, and is positive. Its exact integral uses 2,240 retained cells, preserving the lower and upper input order despite sign changes. The position coefficient uses the already derived positive piecewise polynomial kernel with the independent, slightly sharper e_H above. Its error is (5/3)e_H+(2/3)(K³−1)Z. This reproduces and slightly tightens the archived combined-kernel interval without asserting new information about the input cells.

Outward rounded bounds are

    0.0817743170857 < C_pos < 0.0820898087395,
    0.2528124152261 < C_open < 0.2537906850232,
    1.5625e-8 <= L/rho_b <= 1.562500000015625e-8,
    2.03007002760e-6 < L < 2.03007873184e-6 model lengths.

The report stores exact rational interval endpoints. The table below rounds outward; entries are first derivatives per canonical ε. Normal length entries are additionally per radian of infinitesimal inward launch tilt.

| Statistic | Central v=e₁ | Positive tenth-sky mean |
| --- | --- | --- |
| Normal fractional Jacobi response δj_N/L | (−1.152807, −1.041260)×10⁻²⁴ | (−8.300208, −7.497076)×10⁻²⁵ |
| Normal Jacobi response in model length units | (−2.340289, −2.113832)×10⁻³⁰ | (−1.685008, −1.521959)×10⁻³⁰ |
| Parallel-frame opening response δj_N′ | (−3.564043, −3.219147)×10⁻²⁴ | (−2.566111, −2.317786)×10⁻²⁴ |

The opening-to-fractional-displacement coefficient ratio is between 3.07970 and 3.10356. These are unrescaled dimensionless optical response sizes, unlike the order-0.08 normalized position coefficient.

For any independently supplied physical calibration L_phys, the leading normal displacement difference is

    r_N(ε,α,L)−r_N(0,α,L)
       =ε α L_phys [δj_N/L]+o(ε α)

in the mixed-derivative sense just specified. For example, *if* L_phys were calibrated as one metre, the central displacement coefficient would be about 10⁻²⁴ metres per ε per infinitesimal radian; the opening coefficient is about 3.4×10⁻²⁴ radians per ε per infinitesimal radian. No assignment of the model length to metres is established here. The scale of an actual source amplitude, the finite-amplitude remainder, finite-tilt corrections, and detector noise are all separate. In particular ε=1 is not asserted to be an admissible physical amplitude.

## 7. Independent five-dimensional metric check

The retained longitudinal reconstruction in `FULL_PACKET_JUNCTION_SUFFICIENCY_AND_CURRENT.md` has

    δg_yy=−4Psi,  δg_yμ=0,  δg_μν=2rho²Psi gamma_μν,
    Psi|shell=F.

Thus the radial lapse is N=1−2εPsi+O(ε²). In this first-order representative the y slices are conformal to gamma, so all radial second-form terms in R(n,k,n,k) are proportional to h(k,k) and drop out. The zero-shift radial curvature identity leaves

    R(n,k,n,k)=−N^−1 Hess_h N(k,k) + null-isotropic terms,
    δR(n,k,n,k)=+2 Hess_h₀ Psi(k,k)=+2F″.

This recovers the normal entry and the factor −2 in (3) directly from the actual bulk metric tangent, without using the Ricci-trace argument. The copy-even material shell has zero graph displacement in this gauge. Gauge changes and frame variations do not alter the first optical response because its background tidal tensor vanishes. A fully conformal five-dimensional metric would instead have δN=+F, yield δRic⁵(k,k)=−3F″, and fail this zero-linear-flux check; it is not the registered bulk reconstruction.

A second auxiliary cross-check is the five-dimensional plane-wave metric

    ds²=−2 dλ dv+dX₁²+dX₂²+dX₃²
        +εF″(λ)(X₁²+X₂²−2X₃²)dλ².

Along X=0 its optical tidal matrix is exactly εF″diag(−1,−1,2), and its null Ricci trace is zero. Direct Jacobi integration gives (2). The script checks the two integration expressions independently for the nonconstant polynomial F=λ²+λ³, producing d(1)=5/6. This plane wave is an auxiliary local tidal model, not a replacement Einstein–scalar/Israel solution or a reconstruction of the entire actual bulk.

## 8. Reproduction and retained scope

Run `python -I -B firstorder.py` and compare its JSON with `FIRST_ORDER.json`. Running with `-O` gives the identical report because checks use explicit exceptions, not Python assertions. The sibling inputs and their hashes are in the report; exact source locations and proof hashes are in `SOURCES.json`.

The arithmetic certifies the displayed finite interval operations conditional on the retained continuum input enclosures. The analytic statements are the Gauss/Codazzi, actual longitudinal reconstruction, one-sided screen geometry and Jacobi arguments above. This note establishes no finite-amplitude remainder, numerical inward-tilt domain, shell-crossing prescription, global focusing theorem, astronomical calibration, or physical confirmation of an extra dimension.
