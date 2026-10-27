# Same-event infinitesimal frame independence for the registered first metric

## Statement and scope

Let `g(λ) = g0 + λ g1 + (λ²/2) g2 + o(λ²)` be a sufficiently smooth formal metric family on a five-dimensional region, with background

`g0 = dy² + ρ(y)² h`,

where `ρ > 0` and `h` is a four-dimensional Lorentzian constant-curvature metric. The proof is local and assumes all fields used below are smooth enough for the indicated derivatives. Fix an event, a background unit timelike vector `U0` tangent to a `y = constant` leaf, and the unit normal `n0 = ∂y`. Choose the orientations consistently. Let the registered first metric be

`g1 = -4uY dy² + 2uY ρ²h + L_ξ g0`,                                   (1)

with smooth scalar `uY` and an arbitrary smooth first-order vector field `ξ`. In the application, `Y` is the registered scalar quadrupole and `ξ` is the registered scalar gauge vector; the proof does not need either specialization. The master amplitude and the gauge vector have their common functional dependence exactly as in (1); do not widen their correlated metric components as unrelated functions.

For each λ, take two smooth oriented orthonormal frames `U,n,e_a` and `U',n',e'_a` for the *same metric at the same event*. Assume they have the same zeroth-order frame (the angular frame may instead differ by a fixed proper rotation if components are first identified by that rotation). Their first-order relative Lorentz transformation is otherwise arbitrary. Define

`H_ab(λ) = (1/2) ε_(a^{ cd} C_{b) U c d}(λ)`,                         (2)

where `a,b,c,d = 1,2,3` label the screen, round brackets denote **normalized** symmetrization, and the components are in the oriented orthonormal frame. Then `H1 = 0` and the second parameter derivatives in these matched frames obey

`H2' = H2`.                                                               (3)

The selected local components of these two equal screen tensors are consequently equal. Equality of a sphere-wide basis coefficient additionally requires the comparison at every point of that sphere using the same background identification, basis, orientation, and measure. Equation (3) allows arbitrary smooth first-order relative frame changes, including all boosts, angular tilts and screen rotations. It does **not** compare frames with distinct background observers or distinct events, and it is not a theorem about a general gauge transformation, an invariantly selected event, a constructed geometric clock, or existence of `g2`. It applies to the exact correlated first sector (1), not to a general scalar metric. It also does not turn the selected three-screen projection into the whole five-dimensional magnetic Weyl tensor.

## Curvature lemma

Consider locally any metric

`g̃ = N(y,x)² dy² + γ_ij(y) dx^i dx^j`, `γ = ρ(y)² h`, `N > 0`.

The sign of the unit normal is fixed by increasing `y`. The `y` leaves are intrinsically constant-curvature and their shape tensor is umbilic:

`K_ij = (1/(2N)) ∂y γ_ij = k γ_ij`, `k = ρ'/(Nρ)`.

Choose either standard sign convention for the Riemann tensor consistently. The Codazzi equation says that `R_{nijk}` is a linear combination `γ_ij q_k - γ_ik q_j` for some covector `q` (the overall sign is convention-dependent). In five dimensions its Ricci contraction fixes `R_{ni}` to the corresponding multiple `-3 q_i`. Consequently the Ricci/Schouten subtraction in the Weyl formula removes this whole tensor:

`C_{nijk} = 0`.                                                           (4)

For the all-tangential components, Gauss gives intrinsic constant curvature plus an umbilic quadratic term, hence a multiple of `γ_ik γ_jl - γ_il γ_jk`. The ambient Weyl trace subtraction preserves a tangential Kulkarni–Nomizu form, so for some symmetric leaf tensor `S`,

`C_{ijkl} = (γ ∧ S)_{ijkl}`
`= γ_ik S_jl - γ_il S_jk - γ_jk S_il + γ_jl S_ik`.                     (5)

There is **no extra factor 1/2** in this definition. Weyl tracefreeness in five dimensions and leaf dimension four implies `tr_γ S = 0` and

`C_{n i n j} = -2 S_ij`.                                                  (6)

Indeed, the tangential contraction of (5) is `2S_jl + γ_jl tr S`; taking its trace and using the Weyl traces proves `6 tr S = 0` and then (6). No field equation is used.

At `N = 1`, the background Weyl vanishes. One way to see this is that the full construction has the Lorentz symmetry of the constant-curvature leaf, so the leaf tensor `S` is proportional to `γ`; (6) forces its trace to vanish and hence `S = 0`. With (4)–(6), all background Weyl components are zero. Equivalently, this is the standard conformal flatness of a one-dimensional warped product with a constant-curvature fibre.

## Reduction of the registered first metric

Under the conformal rescaling `ḡ(λ) = exp(-2λuY) g(λ)`, the non-gauge contribution to the first derivative becomes

`ḡ1 = -6uY dy² + L_ξ g0`.

The `uY` term is exactly the derivative at zero of a family `N(λ)²dy² + ρ²h` with `N(λ)² = 1 - 6λuY` locally. The pure infinitesimal diffeomorphism in (1) contributes `L_ξ C0 = 0` to the first Weyl tensor. The conformal rescaling also leaves the first covariant Weyl tensor unchanged because `C0 = 0` (`C_{abcd}` rescales by the square of the conformal factor). Thus the registered `C1`, expressed at the background event, obeys (4)–(6) with `γ = ρ²h`. Differentiating the adapted normals and frames in the curvature lemma introduces no extra terms, again because `C0 = 0`. This argument does not invoke the master PDE, numerical field values, or any second-order field equation.

## Vanishing and arbitrary first-order frame changes

Take any unit timelike `U` tangent to the background leaf and choose an oriented screen orthogonal to `U` within the leaf. From (5), using screen orthogonality,

`C1_{b U c d} = δ_bc S_{Ud} - δ_bd S_{Uc}`.

Contracting with `ε_a^{cd}` gives `2 ε_{ab}{}^d S_{Ud}`, which is antisymmetric in `a,b`. Normalized symmetrization therefore gives `H1_ab = 0`. This holds for every tangent unit `U`, not merely one chart's time direction. The conclusion does not require `S_{Ui}` to vanish.

At the fixed background frame, any infinitesimal Lorentz change decomposes into:

1. A boost of `U` toward the three-screen, with the accompanying change of the screen basis. This varies `U` *within the leaf* while leaving `n` fixed. The displayed vanishing holds for every tangent unit `U`, so its derivative is zero.
2. A boost in the `U,n` plane and/or a tilt of `n` toward the three-screen, with the accompanying frame changes required by orthonormality. Differentiating the tensor expression inserts precisely one normal index in one of the four Weyl slots; every such term vanishes by (4) and Weyl index symmetries. In the frame formula, the orientation tensor has constant components; there is no omitted epsilon variation.
3. A rotation within the three-screen. This acts on the already vanishing `H1` and therefore contributes zero.

Thus the linear derivative with respect to an arbitrary relative first-order Lorentz transformation vanishes. These are all 10 generators of the five-dimensional Lorentz group: 3 screen boosts, 1 `U,n` boost, 3 normal-screen rotations, and 3 screen rotations.

Finally, `C0 = 0`. Expanding (2) in λ, the difference between two second parameter derivatives is twice the first-frame derivative of the projection of `C1`; possible relative second-frame variations multiply `C0` and vanish. The first-frame derivative was shown to vanish, establishing (3). This includes the usual factor of two from the convention `g(λ)=g0+λg1+λ²g2/2` and does not assert that the second-order Weyl tensor itself is invariant under arbitrary gauge changes.

## Algebra control and limitations

`check_frame_algebra.py` independently assembles a five-dimensional algebraic Weyl tensor from a completely general tracefree symmetric leaf `S`, verifies its curvature symmetries and traces, and symbolically differentiates all nine screen components under all 10 infinitesimal Lorentz generators. A negative control adds a nonzero mixed Weyl block while preserving all Weyl symmetries and traces; an infinitesimal radial boost then gives a nonzero response. Both controls run with explicit checks in normal and optimized Python. The script requires SymPy and is an internal algebra control for this proof, not an independent specialist review or a substitute for checking (1) against the physical reconstruction. The result does not certify transport of observers or event labels, calculate first-field error bounds, solve the second-order bulk and shell boundary problem, establish a finite-amplitude exact spacetime, or provide evidence for a higher-dimensional Big Bang origin.
