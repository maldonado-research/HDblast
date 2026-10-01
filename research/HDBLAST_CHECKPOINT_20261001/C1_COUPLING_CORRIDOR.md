# Exploratory extension: a conditional coupling corridor

This extension was proposed and evaluated after REGISTRATION.md was committed. It is exploratory. It strengthens the fixed-history Gaussian budget, not the physical validity of the production ansatz.

At a fixed archived crossing set q=G v_star, Lambda_full=G L_full, Lambda_post=G L_post, and a_rel=a_end/a_star. If the occupation amplitude is bounded by a specified F_max, the characteristic cutoff screen b <= (G L_full)^(-3) gives

Rbar_cap <= [F_max v_star^(3/2)/(8 pi^3 a_rel^3 L_full^3)] sqrt[L_post^2/G + 3v_star/(2pi a_rel^2 G^2)].

The expression in the square root is decreasing for G>0 when v_star>0. Thus a local screen G>=G_min bounds the entire assumed coupling corridor by evaluation at G_min. A larger coupling alone cannot increase this optimistic single-cohort cap when its full-history peak mass forces b down by this rule.

| Crossing | Local G_min | Assumed F_max | Approximate Rbar upper cap | Approximate lower r=W/(sigmahat Rbar/18+Rbar²/36) |
|---|---:|---:|---:|---:|
| 0.5 | 43.184788 | 1 | 3.834048762e-4 | 33.0455654 |
| 0.5 | 43.184788 | 1.3 | 4.984263391e-4 | 25.4193850 |
| 0.9 | 2189.805900 | 1 | 2.050332759e-6 | 6179.6285117 |
| 0.9 | 2189.805900 | 1.3 | 2.665432586e-6 | 4753.5601129 |

For comparison, the target used in the earlier scan was r<=0.1. Under this fixed-history, bounded-Gaussian-population prescription, these optimistic budgets remain well short of that endpoint target throughout the screened corridor.

F_max=1 is an assumed bare Gaussian ceiling. F_max=1.3 is a sensitivity assumption, not a proven upper bound on an arbitrary physical spectrum. This corridor argument requires the archived v_star and mass-history coefficients, a single cohort, monotone expansion, fixed endpoint Weyl and tension, and the particular cutoff screen. It does not exclude repeated events, other distributions, other histories, different couplings, contractions, or unspecified ultraviolet production.

The phi_star=0.9 rows at G=100 and 1000 in C1 do not pass this local screen; the present corridor begins above G=2189.8059. Passing this screen is not proof of a controlled production calculation. Gaussian tails are not a hard momentum cutoff.

These formulae are elementary consequences of the stated envelope and scaling. External mathematical novelty is not claimed.
