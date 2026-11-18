# Primary reference for the five-dimensional geometric mass budget

Verified by direct primary-body reading on 2026-11-18 UTC:
H. Maeda and M. Nozawa, *Generalized Misner-Sharp quasi-local mass in
Einstein-Gauss-Bonnet gravity*, [arXiv:0709.1199v4](https://arxiv.org/abs/0709.1199v4),
revised **2008-02-20 18:41:40 UTC**, Phys. Rev. D77, 064031 (2008),
[DOI:10.1103/PhysRevD.77.064031](https://doi.org/10.1103/PhysRevD.77.064031).
Read the [versioned PDF](https://arxiv.org/pdf/0709.1199v4), §§II–IV.A.
The arXiv record is correct; v4 is the latest version in the retrieved primary history.

## Definition and Einstein limit

Section II, Eqs.(2.6)–(2.7), assumes a two-dimensional Lorentzian base and
(n−2)-dimensional constant-sectional-curvature fibers:
g=g_ab dx^a dx^b+R²γ_ij dz^i dz^j, with k=−1,0,+1.
The paper explicitly assumes compact fibers and writes their unit-metric volume V_(n−2)^k.

Set the Gauss–Bonnet coupling α=0 and n=5 in Eq.(2.7):
m = [3V_3^k/(2κ_5²)] [R²{k−(DR)²}−ΛR⁴/6].

Therefore the requested purely geometric quantity is
μ_0=R²[k−(DR)²], and the explicitly Λ-subtracted normalized mass is
μ_Λ=μ_0−ΛR⁴/6=2κ_5²m/(3V_3^k).
For k=0, μ_0=−R²(DR)². No spherical k=1 term should be inserted.
The geometric normalization alone is not a physical-unit mass calibration.

## Work law and Kodama current

Section III, Eq.(3.1), defines K^μ=−ε^μν∇_νR.
Eqs.(3.9)–(3.13) construct the conserved current; on the field equations
J^μ=−T^μ_νK^ν and ∇_μJ^μ=0.
Eqs.(3.20)–(3.21) identify its quasi-local charge with m and define the areal volume.

Section IV.A, Eqs.(4.1)–(4.3):
P=−T^a_a/2,   ψ_a=T_a^bD_bR+P D_aR,
dm=Aψ_a dx^a+P dV,   A=V_3^k R³,   V=V_3^k R⁴/4.
The trace defining P is over the **two-dimensional base**, not the full5D stress.
V is areal volume; it is not a generic proper spatial-volume integral.
Both the flux term and the work term belong in a mass-budget comparison.

Equivalently, for the source convention G+Λg=κ_5²T and fixed orbit normalization,
D_a μ_Λ=(2κ_5²/3)R³[T_a^bD_bR−T^b_bD_aR].
This follows algebraically from the cited first law.
Use the registered Einstein-equation stress normalization before importing κ_5² factors.

## Flat-orbit and potential bookkeeping limits

A compact flat quotient, such as a specified torus, has finite V_3^0 and fits the
paper's compact-fiber volume convention.
For noncompact flatR³, total V_3^0 is infinite. The local identity still yields a
mass **per fixed fiducial coordinate orbit volume**, or one may explicitly choose
a compact quotient; that density interpretation is an application of the local
formula, not the paper's finite-total-charge theorem for a noncompact plane.
Fix γ and its coordinate normalization across all runs/comparisons.

If the scalar potential is retained in matter T, choose explicitΛ=0 for this formula:
μ_0 then includes that potential's contribution and its work law must use the fullT.
Alternatively split off one fixed constant vacuum contribution as explicitΛ and use
μ_Λ with the remainingT. These are equivalent bookkeeping choices if done consistently.
Do not subtract a variable U(φ)R⁴/6 as if U were a cosmological constant:
its derivative introduces an additional R⁴ dU/6 term and changes the budget.
For pure constant-curvature vacuum, μ_0=ΛR⁴/6 while μ_Λ=0; hence rawμ_0 need not
be positive in an AdS background.
General positivity/monotonicity claims in the paper have additional energy,
trapping and boundary hypotheses and are not inferred for the instrumented run.

This is a bulk warped-product mass reference. Brane/boundary exchange, moving endpoints,
and any distributional surface energy still require their own oriented balance terms.
It does not convert bulk mass loss into thermal brane radiation without a specified channel.
