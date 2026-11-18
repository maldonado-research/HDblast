# How much can completing the B1 decay accomplish?

18 November 2026. Numerical consequences of the archived B1 runs, with an exact algebraic conversion envelope. This concerns the earlier homogeneous, tuned `delta=0.1` B1 model; it does not import results from the later local packet model.

At every retained state in the specified windows, converting **all currently present massive-particle energy** into radiation would still leave the registered radiation-dominance diagnostics far from their targets. Across all eight registered cutoff cells and both grids, the largest fraction of the recorded `H^2` attributable to radiation after conversion is `0.01853141`, and the smallest Weyl/radiation ratio after this generous conversion is `56.7030`. The registered targets are `Omega_r>=0.9` and `|r|<=0.1`.

The nonlinear Friedmann formula distinguishes two deficits. The required radiative **contribution** is at least `567.03` times the available envelope; the required all-radiation **density** is at least `557.67` times the currently recorded total matter density. Merely attaining the `Omega_r>=0.9` criterion, without the Weyl-ratio criterion, requires at least `48.49` times that density. These are fixed-state numerical requirements, not predictions for a dynamically modified run.

## Exact envelope and denominator

Use the archived solver's variables:

* `R=kappa5^2 rho_r` is real decay-radiation density and `X=kappa5^2 rho_chi` is the full massive-species energy density, including kinetic energy.
* `sigma` is the current shell tension; `rho_b=H0^-1` is the run's background conversion scale.
* `h=H/H0`, `W=mathcal W/H0^2`, and the stored `rad` and `chi` are contributions to `h^2`, not matter densities.

At a fixed stored state, set `S=R+X`. The exact radiation contribution for an allocation `r` of this available density is

    g(r)=rho_b^2 [sigma*r/18+r^2/36].                         (1)

All selected states have `sigma>0` and `S>0`, so `g` is increasing on `[0,S]`. Allocating every unit of the existing matter density to radiation gives the maximum

    g_all=g(S)=rad+chi,
    rad=rho_b^2[sigma R/18+R^2/36],
    chi=rho_b^2[sigma X/18+(X^2+2RX)/36].                     (2)

The `2RX` cross term is retained. Treating the estimate as merely linear in density, or summing two separately evaluated radiation polynomials, would omit part of the nonlinear brane term.

The diagnostics use the **original recorded denominator** `h^2` and original `W`:

    Omega_all=g_all/h^2,    |r|_best=|W|/g_all.               (3)

No new expansion rate, fitted vacuum term, fitted Weyl term or radiation-only denominator is substituted. All material energy may be radiation under this allocation while its contribution to `h^2` remains small. `Omega_all` is not the material fraction `R/(R+X)`.

The full archived identity is

    h^2=vac+rad+chi+kin+fric+W,
    fric=-rho_b^2(2 sigma_prime J+J^2)/48.                    (4)

The replay checks (2) and (4) against the raw columns. Because the solver defines `W` through the same geometric identity, the closure check is a normalization/bookkeeping check, not an independent test of the PDE. The separately recorded near-shell Hamiltonian and momentum residuals provide the numerical constraint diagnostics below.

Equation (2) is a generous **frozen-state allocation envelope**, not a completed instantaneous-decay solution of the full boundary problem. Physically converting chi changes its pressure and scalar source `J`. The temporal Israel row and scalar-normal junction then change. One cannot generally hold every metric derivative, normal scalar derivative, `h`, `W`, and source fixed while imposing the new matter state. The result says that even before resolving those dynamical changes, merely reallocating the density currently present cannot meet the same-state targets. Future scalar work, new production, transport, changed Weyl geometry or later dilution are outside this envelope.

## Exact density requirement

To satisfy both registered instantaneous diagnostics at the same `h,W,sigma`, an all-radiation density must obey

    g(S_required)>=T,
    T=max(0.9 h^2, |W|/0.1).

Inverting the full quadratic polynomial gives the stable positive root

    S_required = 36T /
      {rho_b^2 [sqrt(sigma^2+36T/rho_b^2)+sigma]}.             (5)

The additional density requirement is `max(0,S_required-S)`. The replay substitutes (5) back into (1). It handles an empty reservoir explicitly: positive target with `S=0` has an infinite multiplier; zero target and zero reservoir needs no energy but its density ratio is undefined. No such zero case occurs in the selected states.

Throughout the retained samples the Weyl-ratio target is stricter: `10|W|/(0.9h^2)>11.25`. Thus the joint requirement equals the Weyl requirement here. The density multiplier is **not** obtained by multiplying the present density by `T/g(S)`: its quadratic correction is included in (5). Nor can this density multiplier be read as a permissible multiplier of `b`, a coupling, a production yield, or a late-time threshold; a new run would change the scalar and geometry through backreaction.

## Fixed selection and reliability

The source is `HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source`, including its registration and post-audit corrections. All eight registered `b_cut` cells are used:

    phi_star in {0.5,0.9}, G in {100,1000}, y in {0.1,1},
    dz_fine in {10^-3,5*10^-4}.

The registered rounded cutoff values are `6.49e-6`, `6.49e-9`, `1.39e-6`, `1.39e-9` for `(phi_star,G)=(.5,100),(.5,1000),(.9,100),(.9,1000)`, respectively. Calling these the registered cutoff cells does not establish a new exact cutoff calculation.

Before computing this conversion diagnostic, the new retrospective windows were fixed as `s=H0 tau` in `[3,4]`, `[4,5]`, `[5,6]`, and the illustrative state was fixed as the native stored record nearest `s=5`. These are not the original preregistered plateau windows. No window or cell was moved or removed after its conversion result was seen.

The export contains all 1,360 native stored states in `[3,6]`: 85 per run, 16 runs. Source indices are consecutive, with stored records strictly outside both outer endpoints. Every selected state has `H>0`. All three windows are strictly inside every run's **measured-constraint prefix**; its final measured anchor lies at `s=6.94317`–`6.96007`, whereas these runs stop later, at about `s=7.35`–`7.44`. No post-stop extrapolation or last-run-state classification is used.

All measured near-shell Hamiltonian and momentum residuals throughout the admitted prefixes are below the registered `0.05` ceiling. Their maxima are `4.236e-4` on the coarse grid and `1.164e-4` on the fine grid. Constraint values are stored at discrete diagnostic times, so this is not a validated continuous-time residual enclosure. The export includes all measured anchors in the admitted prefixes, including anchors strictly after each window endpoint.

The comparison of corresponding coarse/fine **window extrema** changes the reported quantities by at most `6.334e-5` relatively (`0.006334%`). This is two-grid sensitivity, not a rigorous discretization-error bound. The results below concern the sampled states; the calculation does not enclose every intervening continuous state.

| Fixed window in s | Stored states, all 16 runs | Smallest `|W|/g_all` | Largest `Omega_all` | Smallest `S_required/S` |
|---|---:|---:|---:|---:|
| [3,4] | 624 | 66.0679 | 0.0154329 | 640.359 |
| [4,5] | 432 | 60.8058 | 0.0169967 | 594.659 |
| [5,6] | 304 | 56.7030 | 0.0185315 | 557.672 |

The most favorable state is the fine-grid `(phi_star,G,y)=(.5,100,.1)` run at `s=5.9668838665`. It has contribution multiplier `567.0302968` and density multiplier `557.6723558`. These distinct numbers should not be interchanged.

## Illustrative states near s=5

This table uses the finer grid only; the machine report includes the corresponding coarser-grid values and all window extrema. Native record times are shown so no interpolation is mistaken for a stored state.

| phi_star | G | y | Actual s | `|W|/g_all` | `Omega_all` | `S_required/S` |
|---:|---:|---:|---:|---:|---:|---:|
| 0.5 | 100 | 0.1 | 5.019704 | 60.6013 | 0.0170655 | 592.842 |
| 0.5 | 100 | 1 | 5.019737 | 94.0434 | 0.0110629 | 919.937 |
| 0.5 | 1000 | 0.1 | 4.978320 | 208.675 | 0.00501149 | 2040.44 |
| 0.5 | 1000 | 1 | 4.978946 | 673.845 | 0.00155720 | 6588.59 |
| 0.9 | 100 | 0.1 | 4.979392 | 2924.15 | 0.000359260 | 28590.6 |
| 0.9 | 100 | 1 | 4.979389 | 3333.78 | 0.000315129 | 32595.8 |
| 0.9 | 1000 | 0.1 | 4.979433 | 5013.00 | 0.000209591 | 49014.1 |
| 0.9 | 1000 | 1 | 4.979439 | 8473.24 | 0.000124010 | 82846.1 |

For example, at the first row about 94.38% of the currently recorded matter density is still in chi. Completing that decay greatly improves its radiation readout, but even the generous envelope leaves a Weyl/radiation ratio above 60. This separates a slow-decay deficit from the much larger total-current-reservoir deficit.

## What is established

The current particle reservoir is insufficient to meet the registered instantaneous dominance diagnostics at any retained sampled state, even with perfect immediate conversion and the exact quadratic radiation contribution. This strengthens the earlier linear all-chi-to-radiation estimate with a full-density, nonlinear calculation and a direct density requirement.

It does **not** change the original registered aggregate verdict, **INCONCLUSIVE**, because the required `Delta ln a=0.5` plateau was not reached. It does not disprove later dominance, alternative production mechanisms, access to scalar/bulk energy, different initial states, or other parameter choices. A claim that new matter production succeeds must evolve the extra energy transfer and backreaction; the present numbers cannot replace that calculation.

## Reproduction

Run `python conversion_replay.py`. It uses only Python's standard library, authenticates `STATES.csv`, `CONSTRAINTS.csv`, and `INPUTS.json`, verifies the complete 8-by-2 selection and strict window containment, and emits `CONVERSION.json`. Normal Python and `python -O` give identical output; mathematical checks use explicit runtime conditions. Input metadata pins each original NPZ and summary and the solver/analyzer/registration. `extract_inputs.py` documents the read-only export from those archived NPZ files and needs NumPy only for that extraction. The replay does not rerun a PDE solver, fit a decay model, or modify any prior result.
