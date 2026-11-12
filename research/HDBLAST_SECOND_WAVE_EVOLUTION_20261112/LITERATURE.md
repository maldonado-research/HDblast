# Forced second-wave literature and method check

Scientific cutoff: **6 November 2026**. Retrieval: **12 November 2026 UTC**. This is a focused supplement to the retained first-wave review, not a repeat of its broad literature search. No later publication is used.

The useful findings concern mixed time/spatial source regularity, forcing contributions to initial corner jets, and the angular support of quadratic products. None of the inspected external results directly supplies the numerical estimate for the retained sixteen-field mixed Dirichlet/Neumann shell operator.

## Primary sources and actual coverage

| Source/version | Passages inspected | What can be used |
| --- | --- | --- |
| Matti Lassas, Tony Liimatainen, Leyter Potenciano-Machado and Teemu Tyni, *Stability estimates for inverse problems for semi-linear wave equations on Lorentzian manifolds*, [arXiv:2106.12257v1](https://arxiv.org/html/2106.12257v1), 23 June 2021. | §2, energy spaces (14), product estimate (15), corner recurrence (17), Proposition 1 and (18)–(19); relevant proof discussion in Appendix B. | The scalar all-Dirichlet theorem assumes compatible data and F in a mixed energy space E^s, and concludes u in E^(s+1). Its energy-space algebra supports the derivative bookkeeping below. Its boundary problem, scalar bundle and unspecified constants differ from HDblast. |
| Jeffrey B. Rauch and Frank J. Massey III, *Differentiability of solutions to hyperbolic initial-boundary value problems*, Transactions AMS 189 (1974), 303–318, [primary AMS PDF](https://www.ams.org/journals/tran/1974-189-00/S0002-9947-1974-0340832-0/S0002-9947-1974-0340832-0.pdf). | PDF pages 1–4 of 16: pp.303–306, assumptions (1.2)–(1.3), source-dependent initial recurrence (1.4), corner conditions (1.5), Theorem 3.1 and estimate (3.1). Later proof sections were not retrieved. | A primary regularity theorem explicitly retains source jets at the initial time. It assumes a smooth first-order system with nonsingular normal matrix and primal/adjoint L2 estimates. Those assumptions do not automatically follow from a standard first-order reduction of this mixed wave system. |
| Stefano Baiguera, Troels Harmark, Yang Lei and Nico Wintergerst, *Symmetry structure of the interactions in near-BPS corners of N=4 super-Yang-Mills*, [arXiv:2012.08532v2](https://arxiv.org/html/2012.08532v2), 7 April 2021; related DOI [10.1007/JHEP04(2021)029](https://doi.org/10.1007/JHEP04(2021)029). | Appendix A.1, especially (A.1)–(A.8), and the product expansion (A.14)–(A.15). | S3 scalar product overlaps factor into left/right SU(2) Clebsch–Gordan coefficients; vector contractions have their own coefficients. Only these harmonic identities are used, not any Yang–Mills dynamical conclusion. |

The published version related to 2106.12257 has the changed title *Stability and Lorentzian geometry for an inverse problem of a semilinear wave equation*, Analysis & PDE 18 (2025), 1065–1118, [DOI 10.2140/apde.2025.18.1065](https://doi.org/10.2140/apde.2025.18.1065). A publisher PDF retrieval covered only its first ten physical pages, which do not reach §2; therefore the precise mathematical passages above are pinned to the fully retrieved arXiv v1, not represented as a check of the published theorem. The partial publisher retrieval is recorded in the metadata.

## Source order: spatial H12 alone is incomplete

The retained source requirement at target order m=13 has the form

    S13(F) =
      Σ_(j=0)^12 ||∂t^j F||_(L1_t H^(12−j))
      + Σ_(j=0)^11 ||∂t^j F(0)||_(H^(11−j)).

This is the retained shell proof's interface, not a formula quoted from a different boundary theorem. It includes the top time-commuted L2 source, every lower mixed spatial norm, and initial source jets used to initialize the differentiated equation.

The external comparison is precise: [2106.12257, §2](https://arxiv.org/html/2106.12257v1#S2) defines E^s by controlling ∂t^j u in H^(s−j), charges compatible initial/boundary data and obtains one solution order above the source. Thus an E^12-type source is the relevant degree count for an order-thirteen second wave. Its proposition cannot be used to erase the Neumann rows, replace the retained norm, or assign a numerical shell constant.

For the complete force

    F2[U] = −C2(U,U) −2 E2[U] D U −2 D2(DU,DU)
                                         +2 Q1[U]:H U,

the potentially limiting term remains U D²U, including all Hessian/order corrections. Suppose a new proof establishes the actual mixed first-wave budget

    W14(U)=sup_t Σ_(j=0)^14 ||∂t^j U(t)||_(H^(14−j)).

Then the following **derivative allocation is a model-specific inference**, to be combined with the authenticated coefficient maps and the retained explicit H3-to-L∞ inequality.

For j+k≤12, a term of ∂t^j D_sp^k F2 has at most fourteen total field derivatives distributed across two factors. Put the factor with fewer total derivatives in L∞. It has at most seven derivatives, so appending the three spatial derivatives required for that embedding stays within W14. The other factor has at most fourteen derivatives and is controlled in L2. Coefficient derivatives must be charged separately through the actual mixed order used; they are not fields to which the same embedding argument is applied.

This degree count supports S13(F2) from W14 and finite coefficient bounds. It is not the evaluated product constant, and a bound on N14(t) alone must first be connected to all the mixed time derivatives in W14. The full noncommuted and time-commuted force must be the same retained equation, with its factor-two convention.

## Corner jets depend on the forcing

For an equation solved for u_tt, differentiating its time-jet recurrence introduces the corresponding derivatives of F(0). Smoothness of F by itself does not make these jets compatible with a homogeneous boundary. This is explicit in [Rauch–Massey, (1.4)–(1.5) and Theorem 3.1](https://www.ams.org/journals/tran/1974-189-00/S0002-9947-1974-0340832-0/S0002-9947-1974-0340832-0.pdf), and in the second-order Dirichlet example [2106.12257, (17)](https://arxiv.org/html/2106.12257v1#S2.E17).

For the retained model, the useful route is its actual initial flat collars: prove that every first-wave time jet used in the quadratic source is flat on those collars by the local background recurrence. The differentiated quadratic force then has flat initial collar jets as well. Combine these with the completed second initial values/velocities in its own forced recurrence.

This proves initial corner compatibility at the used finite orders. It does not impose a boundary condition on F2(t) for all later times, and it does not impose homogeneous rows on normal derivatives of V2.

## Quadratic angular products have a larger sector

The retained first wave has group type R=(1,1) under SU(2)_L×SU(2)_R, with positive natural SO(4) Casimir eight. The [S3 overlap formulas (A.1), (A.7), (A.8)](https://arxiv.org/html/2012.08532v2#A1.SS1) provide a primary harmonic-product comparison.

A conservative **representation-theoretic inference** is

    R⊗R = direct_sum_(a,b=0,1,2) (a,b),
    C_(a,b)=2[a(a+1)+b(b+1)] ≤24.

Covariant derivatives, tensor contractions and the complete background coefficient maps preserve natural equivariance when all derivative slots are included. After checking this for the actual force, its quadratic angular support therefore lies in Casimir eigenvalues at most twenty-four. In particular it need not remain in eigenvalue eight; Y² already has a nonzero constant component.

For invariant L2 products this gives the conditional ordered-word bound

    Σ_(|I|=r) ||R_I F2||² ≤24^r ||F2||².

It can simplify an independent forcing audit. It does not restrict the second wave's angular support unless the actual completed second Cauchy data are separately checked. In particular the earlier selected finite-angular correction must not be substituted for the actual canonical second derivative.

## Consequence for the continuation

The sources support the planned mixed-source hierarchy and forced corner recurrence, but none removes their numerical evaluation. The direct application obstruction is concrete: the modern energy-space comparison has only Dirichlet data, while the classical first-order theorem assumes a noncharacteristic normal matrix and given energy estimates. Neither is already the retained coupled mixed tensor problem.

No inspected source changes the remaining need for an actual-family metric-patch/amplitude bootstrap or the observer-clock, cut-area and positive-denominator estimates. No physical confirmation or novelty claim follows from this method review. Retrieval details and limits are in [LITERATURE_SOURCES.json](LITERATURE_SOURCES.json).

