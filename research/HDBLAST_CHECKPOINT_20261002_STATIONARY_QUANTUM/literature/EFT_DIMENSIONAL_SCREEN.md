# Independent dimensional screen for the stationary family

This note is an explicit dimensional analysis informed by the literature review, not a cutoff derivation from a UV completion. It distinguishes a formal stationary solution from a controlled low-energy interpretation. Use the actual geometry and scalar normalization when evaluating the diagnostics.

## Definitions and bulk screen

With the bulk Einstein action normalized by 1/(2 kappa5^2), define M5^3=1/kappa5^2. Then [kappa5^2]=length^3, [r]=[x]=length^(-2), and

gamma=N kappa5^2/L^3=N/(M5 L)^3,

M5 L=(N/gamma)^(1/3).

For an invariant-derived bulk curvature energy E_bulk and e_bulk=E_bulk L,

E_bulk/M5=e_bulk (gamma/N)^(1/3).

This ratio is independent of choosing meters, seconds, or electronvolts for L. A small value is a necessary Planck hierarchy for a conventional gravitational derivative expansion with cutoff no greater than M5; it is not sufficient evidence of control when Wilson coefficients or lower UV scales are unknown. An order-one or larger value removes that hierarchy. It is not a theorem locating a sharp failure threshold.

Define and report the curvature convention explicitly. Possibilities include the largest curvature eigenvalue's square root or a suitable fourth root of a positive curvature-square invariant. Numerical factors differ. In a general Lorentzian geometry a single scalar invariant can vanish despite nonzero curvature, so it should not be treated as a universal bound on all components. The current implementation uses the magnitudes of its tangential and normal curvature scales,

kT=U/6-w^2/12,  kN=U/6+w^2/4,

e_curv=max_sampled sqrt(max(|kT|,|kN|)),

in the declared dimensionless bulk variables. This convention is preferable to guessing the bulk scale from the much smaller shell H. It is a sampled diagnostic over the reported solved radial interval, not a certified supremum over an unexamined bulk extension. Include any independently large gradient/curvature scales when they are known. A tenfold hierarchy convention e_curv/M5hat <= 0.1 can be reported as a conservative screening threshold; it is not a derived universal cutoff, a radiative-stability proof, or a causal-stability result.

For illustration only, if the sampled geometry has e_bulk=0.16 throughout the relevant scan, the N=1 family gives:

| gamma | M5 L | E_bulk/M5 | kappa5^2 E_bulk^3/(24 pi^3) |
| ---: | ---: | ---: | ---: |
| 0.01 | 4.64159 | 0.03447 | 0.0000000550 |
| 1 | 1 | 0.16000 | 0.0000055043 |
| 100 | 0.215443 | 0.74265 | 0.00055043 |
| 10000 | 0.0464159 | 3.44710 | 0.055043 |
| 1000000 | 0.01 | 16.0000 | 5.5043 |

If L is specifically the radius of an included asymptotic AdS region, that region has |kT|=|kN|=1/L^2 and e_curv=1 in this convention. The illustrative 0.16 value is therefore not a bound on an unexamined asymptotic extension.

The final column is a conventional five-dimensional naive-dimensional-analysis loop proxy, not a coefficient extracted from any of the five primary PDFs. The factor 24 pi^3 is conventional and scheme/process dependent. The stronger requirement of a parametrically small curvature-to-UV-scale ratio is not supplied by that numerical denominator. In particular, a small loop proxy does not establish the absence of unsuppressed higher-dimension operators when E_bulk/M5 exceeds one.

The row labeled gamma=0 is defined here as the disabled-loop computational classical control. It is not a claim to realize kappa5=0 in a physical dimensionful action. The relation M5 L=gamma^(-1/3) is used only for the positive-gamma N=1 models. Planck ratios for the zero-loop control should be marked not applicable unless an independent retained kappa5 is specified.

For N=1, large gamma provides no large-N suppression of metric fluctuations. For variable N, the same gamma can arise either from many species or large gravitational coupling, with different loop hierarchies. Always record N and kappa5 independently.

## Conditional four-dimensional brane screen

If a normalizable graviton mode, geometry, and boundary conditions establish an effective four-dimensional regime with M4^2 approximately L/kappa5^2, and the relevant energy satisfies E L << 1, a dimensional matter-loop proxy is

epsilon_brane,4(E) approximately N E^2/(16 pi^2 M4^2) approximately gamma (E L)^2/(16 pi^2).

The physical M4 normalization may contain geometry factors or independent brane Einstein terms; it must be derived from this model before using the equality quantitatively. At E^2=x,

alpha_x=gamma x L^2/(16 pi^2),

alpha_r=gamma r L^2/(16 pi^2)=2 gamma (H0 L)^2/(16 pi^2),

alpha_x=alpha_r [1+b(eta-eta_ref)/2]^2.

Thus gamma*r/(16 pi^2) is dimensionless only when r denotes rhat=r L^2 or L=1 units have explicitly been adopted. The reference r is not itself an additional physical propagating species or energy. At x=r it coincides with the physical mass squared; away from that point physical mass, curvature, external frequencies, and thresholds matter. Use both sqrt(x) and H when screening the stationary solution, with threshold treatment appropriate to their ratio.

This proxy differs from the coefficient-level four-dimensional scalar result in Bento–Melo. Their massless species scaling requires m below the inferred scale and sufficiently large N. For N=1 it is more informative to state the gravitational hierarchy directly than to claim a large-species mechanism. Its value cannot certify five-dimensional bulk control, localization, or the source-field loop expansion.

## Distinct comparisons and species localization

For a typical brane vacuum contribution delta rho of order N H^4/(16 pi^2),

delta rho/(M4^2 H^2) approximately gamma (H L)^2/(16 pi^2),

delta rho/T_AdS approximately gamma (H L)^4/(16 pi^2), with T_AdS approximately 1/(kappa5^2 L).

These compare the loop source to different classical quantities. A small fractional change in the large background tension can give a significant fractional change in its much smaller curvature detuning. If instead x >> H^2 and a finite matched vacuum contribution of natural magnitude N x^2/(16 pi^2) remains, its ratio to M4^2 H^2 scales as gamma (x L^2)^2/[16 pi^2(H L)^2]. The actual finite action, logarithms, and matching conditions determine its value; these estimates do not overwrite those conditions.

A five-dimensional bulk gravitational loop scales as kappa5^2 E^3 times a loop factor. N four-dimensional brane fields do not automatically count as N bulk species. In an established four-dimensional graviton regime the usual species scaling is M4/sqrt(N), up to factors and thresholds. An N^(-1/3) scaling directly applies to N bulk species in five dimensions. A brane polarization coupled to a five-dimensional brane-to-brane graviton propagator can also scale as N kappa5^2 E^3, but establishing that regime requires the actual propagator and any induced brane Einstein term. Neither extrapolation follows from gamma alone.

## Reporting rule

For each root, report the stationary/junction residual separately from these prospective EFT diagnostics. List max(E_bulk/M5), H/M5, sqrt(x)/M5, the bulk curvature convention, and, only if localization is established, E/M4 and the brane loop proxy. Record any additional UV scale and the canonical source normalization when known. Large-gamma rows can remain valid mathematical closure stress tests while their controlled physical EFT interpretation is unestablished. No physical measurement or universal cutoff is inferred from this screen.
