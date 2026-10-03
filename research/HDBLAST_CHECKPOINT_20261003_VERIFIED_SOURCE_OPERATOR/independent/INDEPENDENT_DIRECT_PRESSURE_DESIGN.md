# Independent direct-pressure enclosure design

Status: **DRAFT_METHOD_AND_PROOF_ONLY; NOT_READY_FOR_PHYSICAL_EXECUTION**.

Prepared October 3, 2026 UTC (October 2, 2026 Pacific). This design follows a known diagnostic result. It introduces no new physical result, new blind study, fundamental mathematics claim, or interval certificate for the completed checkpoint. No retained array has been decoded and no registered physical source has been evaluated in preparing this document.

## Exact conditional target

Preserve both sources, both inherited momentum rules, K=64/128/256, the twelve-case universe, incoming unrestricted complex modes at a=-4.5, midpoint -4, endpoint b=-3.5, and all contact conventions. Regard each saved momentum, weight, incoming real/imaginary mode component, producer epsilon and producer pi as its exact represented rational. In particular, pi means 14488038916154245685/4611686018427387904; epsilon means 3777893186295716171/37778931862957161709568. Positive momenta/weights, exact cutoff memberships, and opaque input byte/shape pins must be authenticated before decoding.

The new real analytic source is B(z)=exp(1-1/(1-z²)), h=B or zB, z=eta+4, |z|<1; zero outside. Derivatives are the actual consistent derivatives, not retained native approximations. L=-1/eta and g=4L²h-2Lh'-h''. On the audited central interval |z|≤1/2, no flat endpoint is crossed. Retained native source jets, source phases and endpoints are historical computed values being diagnosed, not exact oracles for the new real source.

The principal independent target is the directly integrated **discrete-contact/matched** action integrand:

    I_d = sum_j omega_j integral_a^b F_j(t) dt,
    omega_j = weight_j k_j²/(2 pi²),
    F_j = L(R_j-3P_j)-3h'(R0_j+P0_j).

R and P are defined separately by varying the density and pressure operators and full fixed-reference W0/W2/W4 subtraction. Write U=u/epsilon, W=w/epsilon. The direct bare contributions are:

    r=[(2k²+3L²)Re(U)-L Re(W)-k Im(W)+Lh'+2L²h]/(2k),
    p=[(2k²/3-L²)Re(U)-L Re(W)-k Im(W)+Lh'-2L²h]/(2k),
    R=r-deltaS_R-4h R0,
    P=p-deltaS_P-4h P0,
    R0=R0_bare-S_R,  P0=P0_bare-S_P.

The proposed implementation derives deltaS_R and deltaS_P independently with pair-valued directional algebra applied to the preserved formal inventory. It must not import either old producer's generated contact evaluator, use density conservation to define pressure, normalize incoming modes, or define the integral by a measured endpoint density difference.

For efficiency, an exact algebraic expansion of the independently derived action formulas gives the modal component:

    F_mode = A(t) Re(u)+B(t) Re(w)+C(t) Im(w),
    A=3L³/(k epsilon), B=L²/(k epsilon), C=L/epsilon.

This simplification is an action-algebra identity, not a Ward-derived pressure definition. The explicit work -3h'(R0+P0) remains separately evaluated, enclosed and reported. The contact integrand remains L(C_R-3C_P)+work with both independent contact operators retained.

Preserve a second raw target I_A, using the registered analytic finite-band contacts, together with M_d=sum omega/(2k), M_A=K²/(8pi²), J=integral Lg and all signed matching corrections. For exact consistent sources:

    I_d = I_A + Delta(C_R^d-C_R^A) + (M_d-M_A)J.

This is a comparison identity, never the algorithm used to compute I_d. Each term needs an independent enclosure. Native endpoint source/operator discrepancies must be explicitly labeled; the formula does not license replacing original stored contacts by new exact-source contacts without reporting the discrepancy. The analytic pressure term -M_A*g/3 and all finite-cutoff source-work terms must remain in I_A.

## Distinct rigorous route: a posteriori polynomial ODE defects

Use a uniform provisional mesh of 64 panels, length H=1/64. All panel centers lie inside the original interval, and half-width r=1/128. Represent real/complex polynomial coefficients by explicit integer dyadic rationals. A separately implemented outward fixed-point interval layer provides bounds; no floating-point sin/cos, unverified complex exponential, mpmath precision agreement, or native phase evaluation enters the certificate.

On each panel with local coordinate x=t-t_left, choose a dyadic forcing polynomial p_g and a proved uniform bound |g-p_g|≤R_g. A source degree24 model on complex disks radius1/8 is feasible in analysis: the independent source-domain proof establishes |g|<64 on those disks, so the exact coefficient Taylor truncation on half-width1/128 is at most 1/(15*2^90). A dyadic coefficient approximation contributes an additional explicit polynomial rounding bound. Translating center coefficients to a left-coordinate polynomial is exact rational algebra or explicitly bounded outward arithmetic.

Construct a complex polynomial W_N, provisionally degree120, with a separately rounded dyadic recurrence approximating

    (n+1)w_(n+1)=i*(2k)*w_n-epsilon*g_n,
    U_(N+1)'=W_N exactly.

U has the preserved initial value and coefficient u_(n+1)=w_n/(n+1), represented as exact rationals or accompanied by its own exact residual. **Compute the actual residual polynomial from the delivered coefficient bytes**:

    R(x)=W_N'(x)-i*(2k)*W_N(x)+epsilon*p_g(x).

Do not infer a zero residual because a recurrence was intended. Coefficient quantization, reciprocal rounding, source coefficient uncertainty and any missing leading term must appear in R or in separately enclosed coefficient errors. For a chosen point polynomial, all these residual coefficients can be exact rational values. A polynomial with interval coefficients requires a specified point representative plus a separately bounded radius; an arbitrary independent choice from every interval cannot be treated as one solution.

Let r_w,r_u be rigorous incoming complex norm-error radii. Since k is real, the homogeneous flow exp(i*2k*x) is unitary. Set

    q_n >= abs(R_n)     (complex modulus or conservative |Re|+|Im|),
    Q(x)=sum_n q_n*x^n + |epsilon|*R_g.

Then for 0≤x≤H:

    |w-W_N| ≤ r_w + integral_0^x Q(s) ds,
    |u-U| ≤ r_u + x*r_w + integral_0^x (x-s)Q(s) ds.

The endpoint radii are obtained by substituting H. These formulas propagate an absolute envelope, not a conservation invariant. In particular a small signed invariant is never substituted for either radius. The proof follows by Duhamel applied to the residual equation; it works for arbitrarily small positive k without dividing by k in a phase recurrence. Exact action coefficients still contain 1/k, so weighted acceptance must bound that dependence separately.

The corresponding integrated error envelopes are exact elementary moments:

    E_W = H*r_w + sum_n q_n*H^(n+2)/[(n+1)(n+2)]
          + |epsilon|*R_g*H²/2,
    E_U = H*r_u + H²*r_w/2
          + sum_n q_n*H^(n+3)/[(n+1)(n+2)(n+3)]
          + |epsilon|*R_g*H³/6.

For action modal coefficients A,B,C, with a_max ≥ sup|A| and d_max ≥ sup(|B|+|C|), the modal integral error on that panel is at most a_max*E_U+d_max*E_W. This bound covers nested forcing moments and accumulated anchor error without testing an endpoint density cancellation. A tighter complex coefficient norm may be independently proved, but is optional. State/source and roundoff radii must remain separately reported.

## Geometry and direct integral

L varies on the panel. It cannot be frozen in the modal coefficients. For q=1,2,3, expand L(t)^q about a negative real center c. The Taylor coefficients are exact rational functions of c; the magnitude of coefficient n is |c|^(-q-n)*binomial(q+n-1,n). With rho=r/|c|<1 and degree p, the omitted tail is bounded by

    |c|^-q * binomial(q+p,p+1)*rho^(p+1) / (1-alpha),
    alpha=rho*(q+p+1)/(p+2) < 1.

This follows because the successive positive tail-term ratios decrease. An implementation must reject disks intersecting eta=0 and alpha≥1. Degree24 is a provisional geometry setting, not an established runtime/width result.

Expand A_P,B_P,C_P independently from these geometry polynomials. Integrate the actual product polynomial A_P*Re(U)+B_P*Re(W_N)+C_P*Im(W_N) **coefficient by coefficient** using exact monomial moments over the declared panel. No endpoint G or density primitive is used. A geometry error radius contributes

    R_A*integral|U| + (R_B+R_C)*integral|W_N|

in addition to the state/source error. Integrals of polynomial absolute norms can be safely bounded by the coefficient triangle using exact monomial moments. This is conservative and allows unresolved width outcomes.

The phase is represented by the ODE polynomial, rather than a standalone transcendental phase evaluator. The residual certifies the whole forced response. Synthetic high-phase tests must verify that the residual grows when polynomial degree is insufficient. High-precision final summation alone cannot repair a bad residual. The primary route's endpoint-primitive reduction is not imported or compared until both independent enclosures are complete.

## Source and contact certification

The contact implementation needs a separate interval time-Taylor algebra coupled to the metric-direction pair algebra. Source jets through order4 occur in the action contact; the polynomial integrand's Taylor degree requires additional derivative coefficients and a proved remainder. This must be implemented and reviewed before an executable registration. A provisional contact Taylor degree32 is deliberately separate from the source forcing degree24.

On a complex radius1/8 disk centered inside [-4.5,-3.5], eta never vanishes, 1-z² never vanishes, and Re(L²)>0. Thus k²+2L² remains in the open right half-plane for every real k>0. The principal square-root frequency is analytic and never zero. A simple metadata-independent lower bound is Re(L²)>1/64, hence |sqrt(k²+2L²)|>1/sqrt(32)>1/8. Tighter k-dependent lower bounds can reduce overestimation after authentication. The source polynomial recurrence provides explicit rational majorants for the needed h derivatives; none may be replaced by sampled maxima.

For each full contact/work function, certify a holomorphic magnitude bound M_C on the full complex disk by interval/polynomial majorants of the independently derived W0/W2/W4 inventory. The point Taylor coefficients and M_C give an actual contact remainder M_C*(r/R)^(p+1)/(1-r/R). A degree label without M_C is not a certificate. Compound pair algebra must respect the stated fixed physical mass, metric prefactors, operator variation and subtraction reference. Time-coefficient generation, branch selection, transcendental source values, inverse frequencies and finite-band asinh/sqrt functions all require enclosed arithmetic.

The signed external/source baseline work is directly integrated as its own Taylor model. No source-work cancellation is assumed. Finite-band analytic contacts use separately implemented formulas and branch/domain checks. Assembling a raw analytic contact after integrating a discrete target must not be represented as direct raw pressure integration.

A useful narrower prerequisite, if the complete direct ledger is computationally infeasible, is a separately registered uniform operator/source-remainder certificate. Its conclusions would concern those operators on the declared domains, and all twelve full ledger results would remain unresolved. It must not be substituted for a successful direct ledger after seeing physical results.

## Momentum weights, cancellation and output

Authenticate weights and nodes exactly. Compute every omega_j as an outward interval around the exact represented rational expression, or exactly as a rational. Preserve every inherited momentum node; use no interpolated common state and no continuum replacement. For nonnegative exact weights, total absolute error is bounded by sum_j omega_j*e_j; if a future rule has signed weights, use sum_j |omega_j|*e_j. Never replace that sum by a signed projection. All twelve cases, both sources and all declared cutoffs must be reported.

Keep enclosures for raw/matched integrals, source work, contacts, modal terms and all signed differences. Preserve native D_S and stored endpoint differences as exact historical computed targets. Evaluate interval attribution conservatively: |D_S|>2e-6; upper|DeltaR-[I_d]|≤0.1|D_S|; lower|[I_d]-S|≥0.9|D_S|. A crossing of any threshold is unresolved. The nominal total integral half-width2e-8 is a fixed acceptance condition, not an achieved bound.

The four proposed 5e-9 error allocations must cover source/mode propagation, dense/direct-integrand remainder, contact/work remainder, and arithmetic/output. Radius components cannot be double counted or silently spent twice. Preserve the existing2e-7 empirical and1e-12 bookkeeping gates as distinct checks; those checks alone do not establish interval inclusion.

Primary candidate arithmetic is integer dyadic outward rounding with provisional precision256bits. Complex norm upper bounds may use |Re|+|Im| to avoid unverified square roots. Rational coefficients and exact source inputs have unambiguous encodings. Internal interval endpoints are explicit integer numerators over powers of two. Publish endpoints as integer/power pairs or exact decimal expansions; if shortened decimal strings are also emitted, round lower endpoints down and upper endpoints up with integer arithmetic, then parse and verify containment. Ordinary format() of a midpoint is not a certificate.

## Independence and pre-freeze controls

The future route should have its own coefficient/residual generator, dyadic interval arithmetic, pair/Taylor contact implementation and validator. Sharing immutable input bytes and the mathematical action is required; sharing generated numerical mode/source arrays, pressure evaluators, primary endpoint primitives, quadrature nodes, or reduction intermediates would defeat the intended algorithm distinction. Any shared formal expression or proof must be explicitly identified and separately checked.

Required fabricated controls include constant/polynomial/oscillatory forcing, nonzero arbitrary complex incoming modes, very small positive k, maximum registered phase scale, deliberately insufficient polynomial degree, biased forcing, wrong phase sign, wrong work sign, omitted pressure/source/contact term, altered momentum weight, forged remainder, endpoint-equality substitution, inward decimal output, and exact metadata-contract/entrypoint fixtures under normal Python and -O. The tests must use explicit fabricated input, never registered source identifiers or retained arrays.

## Feasibility and registration readiness

64panels, source24, mode120, geometry24, contact32, 256bits are candidate engineering choices, not a settled registered production route. A naive Python Fraction implementation can have denominator growth and excessive runtime; a custom bounded dyadic coefficient implementation needs a measured complete fabricated universe. Contact geometry may dominate runtime. No full twelve-case resource receipt exists for this design, so **do not inherit900seconds/262144KiB or label it feasible**.

Before a new physical freeze, implement the actual full route, prove coefficient/remainder bounds, run fabricated full-shape resource benchmarks with both8192/16384-node rules, every cutoff, two sources, output and independent validation, settle a stated resource budget, and freeze/read back the exact files/schema/parameters. If absolute bound or resource feasibility fails, record the attempt. A prospective empirical fallback or narrower uniform-operator study needs its own explicit registration, rather than a change chosen after inspecting physical output.

The method controls a conditional anchored finite-momentum real-source problem. It does not certify incoming-state accuracy, continuum momentum convergence, the original whole-history trajectory, coupled Einstein evolution, reheating, or a higher-dimensional cause of the Big Bang.
