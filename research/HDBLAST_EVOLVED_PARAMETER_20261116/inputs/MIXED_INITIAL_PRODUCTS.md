# All-state covariant coefficient and mixed initial-jet products

14 November 2026. This is a conditional proof, conditional on the retained actual smooth background, the authenticated collar bounds, and the actual reduced Einstein–scalar point map. It does not change an earlier certificate. It supplies an evaluated coefficient and product interface for initial temporal recursion; it proves neither a time-evolved solution bound nor a nonlinear existence theorem.

The read-only source tree is the repository at base `2592e90650cb440cab2ce74f726e3dcf9b9d1247`. The primary inputs are:

* `research/HDBLAST_EVOLUTION_BOUNDS_20261106/COEFFICIENTS.md`: the exact sixteen-row coordinate equation, all-state mixed ordinary base/state derivatives through order sixteen, the separate sharper principal inverse-metric jet algebra, and the flow/reconstruction inputs.
* `research/HDBLAST_BOUNDARY_FLOW_PATCH_20261106/METRIC_ENERGY_CERTIFICATE.md`, equation (2): the actual all-state lower bound `a=-q^{tt} >= 1/18000`.
* `research/HDBLAST_ORDERED_EVOLUTION_20261106/QUADRATIC_FORCING.md`, sections 2 and 4: the exact first-jet change and the evaluated full-tensor `H^3 -> L^infinity` constant.

All coefficient bounds below are on the actual real patch

    0 <= z <= 6*10^-6,  |t| <= 13*10^-6/8,  |U|star <= 10^-3.       (1)

No extension of a background-only derivative bound to this patch is used. The actual background reference Christoffel symbols and their derivatives remain in the coordinate coefficients throughout.

## 1. Norms, connections, and exact all-state equation

The field bundle has sixteen independent components: one scalar lapse-root perturbation, a symmetric four-tensor, one scalar, and a four-vector. Its positive fiber metric is the direct-sum Euclidean/Frobenius metric. Add `dt` as a parallel slot and use the fixed product connection `D` on the bundle and on the base metric

    dt^2 + dz^2 + kappa_S3.

Spatial Sobolev norms use full ordered covariant derivative tensors and the full measure `dz dmu_kappa`, on `M=[0,6*10^-6] x S^3`. For coefficient maps, the convention is stronger than a separate bound on each unit-direction evaluation: all displayed derivative/output/base slots are full Frobenius arrays. The norm of a coefficient with `m` state inputs is the operator norm from the `m` Hilbert state factors to that full array. A sum of absolute entries of the entire array is a valid upper bound. No high-rank coefficient map is Sobolev-embedded.

In a stereographic component trivialization, let `D_i=partial_i+Omega_i` and let `gamma` be the product-base connection. The retained coordinate equation is

    q^{ij} partial_i partial_j U = f_c(x,U,P),
    f_c = c_c + b_c P + d_c(P,P),   P=partial U.

Here `d_c` is symmetric in its combined derivative/fiber inputs. Put `Q=DU` and `K(U)=(Omega_i U)_i`. Expanding the bundle Hessian and then substituting `P=Q-Omega U` in the entire expression gives

    q:H U = f_c(x,U,Q-Omega U)
       + q^{ij}{Omega_i Q_j + Omega_j Q_i - gamma^k_ij Q_k
                          +(partial_i Omega_j-Omega_j Omega_i)U}.       (2)

The `gamma Omega` contributions cancel. The order of the last two connection matrices is shown before contraction; symmetry of `q` permits reversing that order only after contraction. In particular, retaining the old zero-order correction after replacing `P` by `Q` would be incorrect.

Define

    S(x,U)Q = q^{ij}(Omega_i Q_j+Omega_j Q_i-gamma^k_ij Q_k),
    T(x,U) = q^{ij}(partial_i Omega_j-Omega_j Omega_i).

Then the actual full-state coefficients are

    cbar = c_c - b_c K + d_c(K,K) + T U,
    bbar Q = b_c Q - 2 d_c(K,Q) + S Q,
    dbar = d_c.                                                       (3)

These are identities for arbitrary states in (1). Since the actual background is the zero section, `c_c(x,0)=0`, hence `cbar(x,0)=0` identically in every base variable.

Use the fixed time-independent spatial frame `(partial_z,E_1,E_2,E_3)`, where the `E_A` are round orthonormal fields. Thus spatial Hessians remain covariant Hessians, and no moving-frame differentiation is suppressed. In this convention `H_ti U=D_i U_t`. With `a=-q^{tt}`, solve the equation once for the acceleration:

    U_tt = A^i(x,U) D_i U_t + B^{ij}(x,U) Hsp_ij U
          + C(x,U) + E^mu(x,U) D_mu U
          + G^{mu nu}(x,U)(D_mu U,D_nu U),                        (4)

where

    A=2 a^-1 q^{ti},  B=a^-1 q^{ij},
    C=-a^-1 cbar,     E=-a^-1 bbar,     G=-a^-1 dbar.

Spatial indices run over four directions and `mu,nu` over all five. The signs follow from `q^{tt}=-a`. Equation (4) is exact on the entire patch. In particular `C(x,0)=0`.

## 2. Evaluated ordinary-coordinate coefficient ledger

The following three distinctions matter: ordinary coordinate derivatives versus covariant derivatives; coordinate first jets versus `DU`; and the sharper principal-metric branch versus the much larger lower-order coefficient branch.

The retained source proves:

    ordinary mixed (x,U) derivatives of c_c,b_c,d_c, order <=16:
        full coefficient entry sum <= 10^150000;
    aggregate normalized Taylor jet of q, order <=16: <10^17000.

The second assertion is proved in that source by the inverse-matrix Taylor recurrence, with the actual zeroth inverse and reconstruction bounds. It follows that each ordinary mixed derivative of `q` through sixteen has entry sum `<16! 10^17000 <10^18000`. The flow derivatives through eighteen (`10^126`) and the varying physical reference connection have already been charged in these actual inputs; they are not new unknown constants here.

Here are explicit bounds for the fixed product connection. On either stereographic chart `|xi|<=1`, use the complex coordinate polydisc of radius `1/100`. The denominator `1+sum xi_i^2` changes by less than `1/10`, and

    u_i=partial_i log s = -2xi_i/(1+sum xi_j^2),  s=2/(1+|xi|^2).

There `|u_i|<3`. Each coordinate `gamma` entry has modulus at most `9`, and each orthonormal-frame connection entry at most `6`. The induced scalar/vector/symmetric-tensor representation on the sixteen components has entry sum, including all five direction indices, below `2*10^4`; this also covers the much smaller base-connection array. Cauchy's formula therefore gives, for each ordinary derivative through seventeen,

    full entry sum of connection derivatives
      <= 2*10^4 * 17! * 100^17 < 10^60.                         (5)

The same estimate includes `partial Omega`, with the appropriate derivative allowance. In the state patch the sum of the sixteen independent state-component magnitudes is at most `.016<1`; derivatives of the state coordinate function have entry sum at most one for each fixed coordinate derivative, and higher derivatives vanish.

In (3) the largest product is `d_c(Omega U,Omega U)`, with five factors and only two connection factors. Total-order Leibniz gives at most `5^16<10^12`; entry-sum contraction is submultiplicative. The finite displayed additions, numerical multipliers, and the smaller `q,T,S` terms fit a factor `10^2`. Thus the ordinary mixed derivatives of the combined `cbar,bbar,dbar` coefficient arrays through sixteen have the evaluated common bound

    10^150200.                                                (6)

For the reciprocal of `a`, do not use `10^150000`. There are fewer than `10^12` multiindices in the twenty-one variables through order sixteen. The preceding sharper `q` branch bounds the aggregate nonconstant Taylor jet of `a` by `10^18012`; its actual zeroth reciprocal is at most `18000<10^5`. The scalar reciprocal Taylor recurrence, truncated at sixteen, gives

    aggregate jet(a^-1)
      <= 10^5 sum_(j=0)^16 (10^18017)^j
      <= 10^5 (1+10^18017)^16
      < 10^288293.

Restoring ordinary derivatives costs at most `16!<10^14`. Therefore every ordinary mixed reciprocal derivative through sixteen is `<10^300000`. This calculation is at an arbitrary state of (1); it does not invert only the background `a`.

The products in (4) use at most two of these factors. Their Leibniz sum is at most `2^16<10^5`; a further `10^2` covers the factor two and the five displayed coefficient groups. Consequently the combined ordinary mixed derivative coefficient arrays for `A,B,C,E,G`, through total order sixteen, are bounded by

    10^450225.                                                (7)

All these calculations use entry sums, so every component contraction has already been charged.

## 3. Full-state horizontal and vertical derivative conversion

An ordinary fixed-component base derivative is not the horizontal derivative of a nonlinear bundle map. For such a map `L(x,U)`, with its output, base, and state-input slots understood, define

    nabla_i^hor L
      = partial_i L + induced connection actions on all tensor slots
                       - L_U[Omega_i U].                    (8)

The input-slot actions have their usual dual signs. At the zero state the final term vanishes, but it is essential away from zero. With (8), the exact nonlinear chain rule is

    D_i[L(x,U(x))] = (nabla_i^hor L)(x,U) + L_U(x,U) D_i U.

Vertical derivatives are the ordinary fiber derivatives with their new state-input slots. Further horizontal derivatives act also on all previously created derivative and state-input slots. We allow arbitrary ordered mixed words in horizontal base and vertical state derivatives, of total order at most sixteen. This is stronger than the particular horizontal-then-vertical orders used by the temporal recursion.

Here is a finite array ledger for this conversion. Initially a lower coefficient has at most eight elementary tensor slots: two output slots and at most two inputs with three slots each. A vertical derivative adds at most two elementary slots; a horizontal derivative adds one. Thus through total order sixteen there are at most forty elementary slots. Using whole sixteen-dimensional state-slot matrices instead of elementary slots only decreases the number of required connection actions.

At an arbitrary coordinate center, use the finite Taylor entry-sum norm in twenty-one variables, truncating at order `16-k` after `k` differentiations. Ordinary normalized coefficient derivatives (7) give initial aggregate norm at most `10^450237`. From (5), the connection Taylor aggregate is at most `10^72`. The affine state-coordinate jet has aggregate norm less than `17`, so the coefficient array `Omega_i U` has aggregate norm below `10^74`.

On a Taylor polynomial of degree at most sixteen, a coordinate derivative increases its entry-sum norm by at most sixteen. The first-order differential part `partial_i-(Omega_i U)^A partial_UA` therefore costs less than `16*10^75`; the connection action on at most forty slots costs at most `40*10^72`. Taking the full array over the five new base directions multiplies by at most five. A vertical derivative, including the full array of sixteen possible new state directions, costs at most `16*16`. Both operations are therefore bounded by `10^80` times the previous aggregate norm. This estimate is valid after each truncation; it never differentiates beyond the supplied ordinary order sixteen.

After sixteen such operations the entry-sum bound is thus

    10^450237 (10^80)^16 = 10^451517.

It remains to read these arrays in orthonormal rather than coordinate base slots. The stereographic scale lies in `[1,2]`. At most forty elementary slots occur, and each coordinate/orthonormal conversion factor is at most two. Symmetric-entry duplication in the field/state slots and the direct-sum normalization are also covered by a factor `2^64<10^20`. We allow the larger factor `10^30`. Hence all full covariant arrays in question are bounded by `10^451547<10^452000`. At each point the quaternion and stereographic orthonormal frames differ by an orthogonal rotation. The tensor norm is invariant under that rotation; derivatives of a rotation are already precisely the connection actions in (8), and are not discarded.

We record the simpler common result:

    Cstar = 10^500000 bounds every ordered mixed horizontal/vertical
    derivative of A,B,C,E,G through total order sixteen,
    in the full Frobenius base/output-array and Hilbert state-operator norm. (9)

In particular (9) includes the sum needed for all spatial derivative directions in Sobolev norms; it is not merely a bound with those directions separately fixed. It holds for all sixteen fields and every angular mode. The frame connection and the nonlinear state argument have both been treated on the whole patch (1).

## 4. Exact zero-state extraction

The zero-state identity is preserved in (4), not replaced by a bound on a constant forcing term. Since `C(x,0)=0`, every horizontal base derivative of this restriction is also zero. Fiberwise integration on the star-norm convex ball gives

    C(x,U) = Ctilde(x,U)[U],
    Ctilde(x,U) = integral_0^1 C_U(x,theta U) dtheta.            (10)

Horizontal transport commutes with radial fiber scaling, so the horizontal derivative of the integrand is the derivative described in (8) at `theta U`; vertical derivatives add factors `theta<=1`. Thus `Ctilde` has the same bound `Cstar` through total mixed order fifteen, using (9) through sixteen. In particular twelve initial-row derivatives, with a reserve of two parameter derivatives if needed, are available. This observation supplies a coefficient reserve; it does not by itself assert a bound for a parameter-differentiated jet hierarchy.

Every term after differentiating (4), with (10), contains at least one actual field factor. This remains true after base differentiation. It prevents a nonzero background-volume contribution from entering the recursive estimate.

## 5. Evaluated mixed-jet product interface

The retained tensor embedding supplies, without an endpoint boundary condition,

    ||T||_infinity <= 10^7 ||T||_H3(M)                         (11)

for full tensor/direct-sum sections with at most sixteen rotating/spatial slots. Its proof uses the complete scalar spherical spectrum with full volume `2pi^2`, giving the Hilbert-valued evaluation constant `1/sqrt(12)`, the full tensor connection estimate below `1000`, and the radial trace factor `Z^-1/2+sqrt(Z)<1000` for `Z=6*10^-6`. The constant is independent of the number of components and involves no angular cutoff or unproved extension/cover constant.

Let `J_j=partial_t^j U|t=0` be known for `0<=j<=r+1`, `0<=r<=12`, and put

    S = sum_(j=0)^(r+1) ||J_j||_H^(14-j),     S <= 10^-12.      (12)

Equation (11) and `|J_0|star<=sqrt(2)|J_0|Hilbert` imply `|J_0|star<=2*10^7 S<=2*10^-5<10^-3`, so the actual state remains in the coefficient patch. The smallness condition is a state-membership condition as well as a product convenience.

Differentiate the right-hand side of (4), with (10), `r` times in time and then `k` times spatially, with `0<=k<=12-r`, keeping the order of covariant derivatives. Each resulting monomial has total field derivative order at most `r+k+2<=14`, where a time derivative and a spatial derivative each count one. At most one field factor has order at least eight. Choose a maximal-order factor for `L2`; its norm is at most `S`. Every other factor has order at most seven and is bounded in `L^infinity` by `10^7 S`, because the necessary three further spatial derivatives are included in (12). Its tensor rank is at most `2+7`, with the three-derivative embedding reserve still below sixteen. High-rank coefficient maps are controlled by (9), not by (11).

The time/spatial selector maps in `D_mu U`, in the Hessian, and in the principal mixed term are contractions or orthogonal projections on full arrays; they add no dimension factor. Spatial and time derivatives commute here because the product connection is time-independent. No two spatial derivatives are commuted, and no homogeneous boundary condition is assigned to a normal derivative.

For completeness the finite differentiation count is as follows. There are five normalized groups in (4); after replacing `C` by (10), each begins with at most two field factors. After `k` differentiations there are at most `k+2` field factors. A further derivative has at most two coefficient-chain choices (horizontal or vertical) and one choice for each field factor, hence at most `k+4` branches. Through twelve derivatives the number of ordered leaves is at most

    5*(4*5*...*15) = 5*15!/3! < 10^13.

This counting includes Leibniz multiplicities by retaining repeated leaves. Covariant slot permutations are isometries. A monomial with `m>=1` field factors is bounded by

    Cstar * S * (10^7 S)^(m-1) <= Cstar*S.

The output `H^(12-r)` square sum adds at most `sqrt(13)`. The full arrays, including all state and base input slots, are already charged in (9). Thus the right-hand side row has bound at most `sqrt(13)*10^13*Cstar*S`, and the convenient rounded interface is

    ||J_(r+2)||_H^(12-r) <= M*S,       M=10^500200.            (13)

There is considerable explicit numerical reserve in (13); the displayed term count and full-array norm convention, not an unspecified Moser constant, justify it. All homogeneous estimates use weak inequalities, so the zero hierarchy is included.

The recursion is triangular: after `r` time differentiations the highest jet in a coefficient is `J_r`, in a first-jet term it is at most `J_(r+1)`, and in the spatial-Hessian term it is `J_r`. There is no `J_(r+2)` on the right after the initial normalization by the actual all-state `a^-1`.

These statements first follow for smooth jets. They extend to the indicated Sobolev jets by smooth density, the strict margin inside (1), and the same finite low/high product estimates applied to differences; the extra coefficient derivative needed for a difference is available in (9). No additional temporal solution is assumed for this purely initial construction.

## 6. Scope of the handoff

The concrete delivered constants are the all-state coefficient `Cstar=10^500000`, tensor embedding `10^7`, hierarchy guard `S<=10^-12`, and row factor `M=10^500200`. The companion `NONLINEAR_INITIAL_JETS.md` uses these to evaluate the thirteen-row recurrence and its amplitude restriction. The present proof preserves the registered nonlinear coefficients, the exact fixed bundle first jet, all spherical modes, all derivative slots, and the actual zero-state source identity. It makes no evolved Picard, normal recovery, or proper-observable assertion.
