# C1 — Energy envelope for one produced particle cohort

**Status:** analytical inequality under explicit assumptions; exploratory numerical application to an archived fixed background. It is not a universal no-go theorem, a renormalized production calculation, or a newly coupled five-dimensional solution.

## Definitions and units

Take one real scalar degree of freedom with comoving occupation f(k)=F exp(-pi k²/q), using comoving scale factor a_star=1 at production. Then

n_star = F q^(3/2)/(8 pi³),  <k²> = 3q/(2pi).

The occupation amplitude F and the Gaussian shape are assumptions or archived ansatz inputs. No spin or species degeneracy is silently included. Hatted masses and momenta are measured in H0; densities in H0^4. The gravitationally normalized radiation is Rbar = kappa5² R/H0 = b Rhat, with b=(H0/M5)^3 for the convention used in the September analysis. It is not the radiation density R_model stored by the old solver: R_model=(b/rb) Rhat and Rbar=rb R_model.

For the fixed endpoint, sigmahat=rb sigma(phi). Its radiation contribution to H²/H0² is

rad = sigmahat Rbar/18 + Rbar²/36.

## Assumptions

Between production and the endpoint, a is monotone increasing. Particle mass remains nonnegative and at most M_max. Each particle either survives or decays once into freely redshifting radiation, with nonnegative probabilities whose sum is at most one. The supplied distribution bounds the specified cohort's comoving population. No additional production or incoming energy is attributed to that cohort.

Initial radiation is counted separately. Multiple production events require separately summed bounds. A contracting history, unbounded particle mass, or an unbounded extra population violates the assumptions.

## Inequality and proof

Let a_rel=a_end/a_star. For the specified cohort,

Rhat_cohort + Xhat <= n_star a_rel^(-3) sqrt[M_max² + 3q/(2pi a_rel²)],

where Xhat is surviving particle energy and Rhat_cohort is decay radiation energy at the endpoint.

For a particle with comoving momentum k, survival energy is sqrt(k²/a_rel²+m_end²). If it decays at a_d, its daughter radiation energy at the endpoint is

(a_d/a_end) sqrt(k²/a_d²+m_d²)
= sqrt(k²/a_end²+m_d² a_d²/a_end²).

Since a_d<=a_end and m_d<=M_max, each possible outcome is bounded by sqrt(k²/a_end²+M_max²). Multiply by the nonnegative outcome probabilities, whose sum is at most one, and integrate the number measure. Concavity of sqrt gives the stated Gaussian-moment bound.

No monotonic mass evolution is needed. A delayed decay may yield a greater budget than an early decay, but it still obeys this envelope. This does not guarantee that a given mass history is affordable: the scalar/background must pay its work through the appropriate energy ledger.

Add Rhat_initial a_rel^(-4) when initial radiation is present. The bound cannot be interpreted as a maximum over all future endpoints: Weyl, vacuum, mass, tension and expansion must be re-evaluated at a later endpoint.

## Fixed archived endpoint

Use the fine A1 Y=0 history at H0 tau=7.353355110666192:
phi_end=1.0361021913225144, ln(a_end)=2.0958245770395325,
Weyl/H0²=0.003666076408891769.
The tuned model has rb=7.838285538073243, delta=0.1,
c=0.5975949350280132, d=-3.106933495673783,
giving sigmahat=5.208194722502616.

The rows use the archived corrected population F=number_factor, the full-history characteristic/peak-mass screen b<=Lambda_full^(-3), and post-crossing mass ceiling Lambda_post. Values in this table are approximate lower ratios or upper budgets, not rounded exact rational bounds.

| phi_star | G | Local production screen | Rbar upper cap | Approximate lower r=Weyl/rad | Additional budget factor needed for r<=0.1 |
|---|---:|---|---:|---:|---:|
| 0.5 | 100 | Pass | 2.46307550e-4 | 51.4397455 | 508.300 |
| 0.5 | 1000 | Pass | 7.90636070e-5 | 160.2532642 | 1583.512 |
| 0.9 | 100 | Fail; formal Gaussian budget only | 9.54684805e-6 | 1327.1695109 | 13114.085 |
| 0.9 | 1000 | Fail; formal Gaussian budget only | 3.02897582e-6 | 4183.0289624 | 41333.500 |

The optimistic radiation target is obtained by inverting the full linear-plus-quadratic radiation contribution:

Rbar_required(r_target)
= 36 Weyl/r_target / [sqrt(sigmahat²+36 Weyl/r_target)+sigmahat].

The numerator is positive for this endpoint. Surviving matter and all decay radiation are already included in the optimistic envelope, so merely changing the decay time cannot raise this specified cohort beyond the bound.

The local production screen requires G>=43.184788 for phi_star=0.5 and G>=2189.805900 for phi_star=0.9. Thus the latter two rows do not validate particle production at G=100 or 1000. They are retained to expose the ansatz limitation. The narrower-window estimates give approximately 43.252847 and 2191.675334 respectively.

A Gaussian has arbitrarily large momentum tails. Screening its characteristic momentum and peak mass is not a complete ultraviolet/EFT certification. Integrating the whole Gaussian is optimistic relative to truncating this very distribution; it says nothing about an unspecified additional ultraviolet source.

## Conditional species scaling

For N copied, independent, identical scalar species with the same couplings, initial state and geometry, total source terms in the classical mode-plus-backreaction equations depend on b_eff=N b. This mapping assumes identical histories and does not apply to arbitrary interacting sectors.

If one additionally assumes a species cutoff Lambda_sp=C M5 N^(-1/3), the screen b N Lambda_dyn³<=C³ becomes b_eff Lambda_dyn³<=C³. This is a conditional mapping, not an established cutoff for brane-localized fields or a universal species no-gain theorem. A different cutoff exponent alpha scales an optimistic multiplicity budget as N^(1-3alpha); different masses, interactions, histories and repeated crossings are outside the copied-species mapping.

## Verification and next step

code/cohort_verify.js reproduces the endpoint rows, independently integrates the Gaussian moments, tests bounded positive decay histories including nonmonotone masses, and includes contraction and omitted-kinetic counterexamples. Results are saved in outputs/cohort_verify.json. The [coupling-corridor extension](C1_COUPLING_CORRIDOR.md) was explored after registration and has separate assumptions.

The physical next step is an action-derived source on an actual history, with renormalized stress/current and background work accounted for. September B1 stays INCONCLUSIVE; this bound diagnoses a conditional budget, not a revised simulation verdict.
