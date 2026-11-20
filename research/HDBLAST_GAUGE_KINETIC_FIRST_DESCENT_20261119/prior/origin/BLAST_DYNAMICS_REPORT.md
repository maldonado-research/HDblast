# HDBLAST Chat 9 - the "blast": nonlinear roll-off of the shell modulus in the 4D effective theory

16 September 2026. Everything here is computed inside the two-derivative Hamilton-Jacobi EFT (orchestrator result R3),
t -> 0. Nothing here is a 5D nonlinear result. Grades: **[EFT]** derived/numerically integrated within the EFT;
**[EST]** order-of-magnitude estimate; **[UNK]** unknown. Machine-readable numbers: `BLAST_RESULTS.json`
(merged), `BLAST_HOMOGENEOUS_RESULTS.json`, `BLAST_PARTICLE_RESULTS.json` (full spectra), `BLAST_REHEATING_RESULTS.json`,
`BLAST_TABLE_CHECKS.json`. Scripts: `blast_common.py` (tables), `blast_homogeneous.py`, `blast_particles.py`,
`blast_reheating.py`, `blast_summary.py` (numpy only, each < 2 min).

## 0. Set-up and checks

* Einstein frame, M_pl = 1: `g_E = f g_J`, `V_E = t v(Theta)`, `v = (1+c phi_b)/f^2`. With `tau = H_top t_E`,
  `H_top^2 = t v0/3`, `v0 = 0.2330302`, `u = v/v0`:
  `Theta'' + 3h Theta' + 3u_Theta = 0`, `h^2 = Theta'^2/6 + u - k/a^2`, `a''/a = u - Theta'^2/3`.
  **The classical dynamics is t-independent in Hubble units**; t enters only through the initial quantum displacement
  `dTheta = H_top/(2 pi)` and through 5D-unit diagnostics (which scale as sqrt t). For t = 1e-3: `H_top/M_pl = 8.813e-3`, `dTheta = 1.403e-3`.
* Tables rebuilt with a cancellation-free system (`I' = -D`, `D' = (2W/3)D - (2/3)W_phi^2 I`, `D = 1 - 2WI/3`); they reproduce the
  orchestrator's `eft_landscape.npz` (Theta, f, V) to <= 5e-8, `I(0) = I_plus` to 1e-13, hilltop `mu^2 = -7.719795` from the table.
* Measured linear growth rate (k=0, pure growing mode): s = 1.653 (+ side), 1.662 (throat side) vs 1.6574 predicted. Friedmann-constraint drift <= 9e-8 (+ side), 9e-6 (throat side). Energy ledger closes to 7e-8 / 9e-6.

## 1. Two structural findings about the field space (both [EFT], new relative to R5)

**(a) The field space does NOT end at phi_b = -1 (Theta = 0.5479).** y_b is a bad coordinate there, phi_b is a good one:
`Z(phi_b -> -1) = 9/23 = 0.391304` (numerically 0.3913044) is finite and equals the norm of the brane-localised scalar zero mode
`2 int_0^inf e^{-(10/9)s} e^{-4s} ds`; `dphi_b/dTheta = -2.1447` and `dlnV_E/dTheta = -3.185` are continuous across phi_b = -1.
The EFT continues smoothly onto the second throat-regular BPS branch `phi = -coth y` (phi_b < -1, bulk field rising from phi_b to -1
in the throat). On it `Z_E > 0`, `f > 0` at least down to phi_b = -6.7; `V_E = 0` at `phi_b = -1/c = -1.6734` (Theta = 0.8341), and
`V_E < 0` beyond (minimum `V_E = -0.0895 t` near phi_b ~ -2.9; W < 0, i.e. negative BPS tension, for phi_b < -2.104).

**(b) The phi_b -> +1 side ends at finite distance too, Theta_end = -3.46529**, but there it is a genuine end (shell at the AdS_+
boundary, y_b -> -inf, `DeltaTheta ~ e^{y_b/9}`): `f = 9 - (3/2) DeltaTheta^2`, `V_E = t (1+c)/81 * (1 + DeltaTheta^2/3)`, i.e. a
quadratic minimum AT the end point with **m^2 = 2 H^2 exactly** (numerically 1.99999998) - overdamped (s = 1, 2).

## 2. Homogeneous roll-off, registered linear detuning, t = 1e-3 (all [EFT])

Start: dTheta = +-H/(2 pi); k=0 on the growing mode, k=+1 at the time-symmetric bounce (a = 1/H, adot = 0, Theta' = 0).
(The H/2pi amplitude for a field with |m^2| = 7.7 H^2 is itself an [EST].)

| event (Einstein-frame e-folds N_E, time tau in 1/H_top) | + side k=0 | + side k=1 | throat k=0 | throat k=1 |
|---|---|---|---|---|
| V_E down 1% | N=2.55 (tau 2.55)| 1.74 (2.43) | 2.45 (2.46) | 1.65 (2.34) |
| V_E down 10% ("roll-off") | 3.37 (3.39) | 2.56 (3.26) | 3.07 (3.08) | 2.25 (2.95) |
| V_E down 50% | 4.29 (4.43) | 3.47 | 3.44 (3.48) | 2.62 (3.35) |
| peak kinetic energy K/V0 | 0.0858 at N=4.30, phi_b=0.903 | 0.0860 | 0.544 (phi_b=-2.38, beyond the throat point) | 0.546 |
| max eps_H | 0.444 (w <= -0.62): **acceleration never ends** | 0.445 | eps_H=1 at N=3.45, phi_b=-1.09 | N=2.64 |

* e-folds before roll-off (k=0): `N ~ 3.37 + ln(sqrt(1e-3/t))/1.657`: 2.68 (t=1e-2), 3.37 (1e-3), 5.46 (1e-6), 7.54 (1e-9), 9.63 (1e-12).
  This is Linde's fast-roll regime (hep-th/0110195): N = (1/s) ln(Theta_end/dTheta). The roll-off itself lasts ~1 Hubble time
  (FWHM of dphi_b/dtau_J = 1.46/H_topJ, equivalent sech^2 duration tau_p = 0.826/H_topJ; from the Hdot_J pulse tau_p = 0.767/H_topJ), i.e. ~1.4-1.6 linear growth times 1/s.
* **(i) toward phi_b -> +1: the modulus settles; no oscillation, no boundary hit, no "blast".** K never exceeds 15% of rho;
  100% of the potential drop (0.9154 V0) is removed by Hubble friction (ledger residual 7e-8). Late time: `DeltaTheta ~ e^{-h tau}`
  (measured dln DeltaTheta/dtau = -0.99997 h), `h_E -> 0.290927` (predicted sqrt(v_inf/v0) = 0.290927). Brane frame: H_J falls monotonically
  from H_topJ to 0.606400 H_topJ = sqrt((1+c)t/27), the exact Randall-Sundrum dS brane in AdS_+ with sigma = 2W(1)+t(1+c) at O(t).
  Bulk velocity: `sinh psi = dy_b/dtau_J -> -2.18922 sqrt(t)` vs the exact AdS kinematic identity `sinh psi = -l_+ H_J = -2.18924 sqrt(t)` (H = alpha sinh psi
  of the moving-shell theorem with beta = 0): the shell recedes toward the AdS_+ boundary at constant non-relativistic rapidity (|sinh psi| = 0.069 at t = 1e-3, the maximum over the whole run).
  The end state is again de Sitter, with V_E lower by a factor 11.8 - an inflating RS2 brane, not a hot universe.
* **(ii) toward the throat.** Requested diagnostic: the rapidity relative to the phi-leaves, `sinh psi = (dphi_b/dtau_J)/W_phi(phi_b)`, exceeds
  0.1 / 0.3 / 1 at 1+phi_b = 0.126 / 0.045 / 0.0130 (y_b = 1.35 / 1.89 / 2.52; t = 1e-3), in general `|sinh psi| = 1 at 1+phi_b = 0.44 sqrt t`.
  In the MOVING_SHELL theorem's frozen-bulk junction (S), `cosh psi * phi'(y_b) = -sigma_t'/2`, the same thing appears as cosh psi - 1 = tc/(2|W_phi|).
  **However this is a statement about the y_b chart, not about the EFT expansion parameters** [argued, not proved]: as W_phi -> 0 the bulk is AdS_-,
  whose boosts are isometries, the phi-leaves carry no energy, and the physical degree of freedom is the localised scalar zero mode with finite norm 9/23.
  The invariant small parameters at the crossing are `|dphi_b/dtau_J| = 0.0278` and `H_J l_- = 0.0182` (5D units, t = 1e-3; both ~ sqrt t). At the crossing
  (tau = 3.46, N_E = 3.42) K/V0 = 0.200, V_E/V0 = 0.532. Note that near the brane grad(phi) becomes timelike there (|phidot_b| > |sigma_t'|/2 = tc/2), so a frozen-bulk/leaf description is what fails.
  EFT continuation: V_E = 0 at tau = 3.67 (phi_b = -1.673), brane-frame Hubble rate H_J = 0 at tau = 3.73 (phi_b = -1.94) followed by **brane-frame collapse**
  (H_J = -31 H_topJ at phi_b = -6.7, because f grows steeply), eps_H up to 8, w >> 1 (negative V_E; cf. Felder-Frolov-Kofman-Linde hep-th/0202017).
  Genuine EFT failure for t = 1e-3: |H_J| > 0.1 (5D units) at phi_b = -5.1, |dphi_b/dtau_J| > 0.3 at phi_b = -6.0; for t = 1e-6 neither threshold is reached inside the table (phi_b >= -6.7).
  The Einstein-frame scale factor is still expanding at the end of the table (h_E = 0.069).
  **[UNK]** the 5D fate. The EFT says: through the throat point, onto the -coth branch, into negative tension and a collapsing induced metric heading for the
  singular end of the BPS flow. Whether the full 5D evolution (bulk radiation into the throat, the cone/horizon, higher-derivative terms whose coefficients near W_phi = 0 were not computed) does this is not known. The negative-V_E region lies outside the field range of the registered solution (phi in (-1,0]) and depends on the cubic W being unbounded below.
* **Inhomogeneity caveat [EST]:** with mu^2 = -7.7 every mode with k/a <~ |m| = 2.8 H grows, and the sign of dTheta is random per Hubble patch. A real shell
  fragments into (i)-type and (ii)-type regions separated by walls where phi_b ~ 0 (tachyonic/spinodal decomposition, cf. Felder et al. hep-ph/0012142). The homogeneous runs describe single patches.

## 3. Brane (Jordan) frame [EFT]

`a_J = a_E/sqrt f`, `dtau_J = dt_E/sqrt f`, `H_J = sqrt f (H_E - fdot/(2f))`. + side: f rises 2.07 -> 9, so N_J < N_E (3.29 vs 3.37 at roll-off); H_J/H_topJ: 1 -> 0.844 (peak of phidot) -> 0.6064, never negative; peak -Hdot_J = 0.251 H_topJ^2.
Throat side: N_J slightly > N_E until phi_b = -1 (f falls 2.07 -> 1.8), H_J = 0.797 H_topJ at the crossing, then turnaround as above. In 5D units H_topJ = 0.012686 (t = 1e-3), matching the registered h = 1.6093e-4.

## 4. Particle production along the + side trajectory (k=0), brane frame [EFT background + standard QFT in curved space]

Method: `v = a^{3/2} chi`, `Om^2 = k^2/a_J^2 + m_eff^2 - (9/4)H_J^2 - (3/2)Hdot_J`, Bogoliubov ODEs in the WKB basis with first-order adiabatic in/out numbers,
in-state = Bunch-Davies in the hilltop dS (9 extra e-folds prepended on the same growing-mode trajectory; modes start at k/a >= 100 m). Vectorised over 40 k's.
Code validated on the Chat 8 exact sech^2 result (0.18882258 vs 0.18882259; 7.49775e-3 vs 7.49775e-3; 0.8442564 vs 0.8442563; integer lambda: 3e-16) and on the
Bernard-Duncan tanh step (6.952117e-4 both). Independent check: the UV plateau reproduces the final-dS steady production `1/(e^{2 pi mu_f}-1)`: 1.191e-6 vs 1.193e-6 (m=1.6), 9.52e-9 vs 9.66e-9 (m=2).
Units H_topJ = 1; kappa = k/(a_* m), a_* at the peak of phidot_b. "Excess" = N_k minus the final-dS plateau; densities at a_*.

| case | N_IR (kappa=0.02) | N_max | n_excess/H^3 | rho_excess/H^4 | sech^2 envelope at kappa->0 (tau_p=0.767) |
|---|---|---|---|---|---|
| (a) minimal, m=1.6 | 0.155 (initial dS plateau 0.031) | 0.179 | 5.0e-4 | 9.3e-4 | 1.8e-3 |
| (a) m=2 | 2.9e-3 (2.5e-4) | 2.9e-3 | 7.5e-5 | 1.8e-4 | 2.6e-4 |
| (a) m=3 | 4.8e-6 (8.1e-8) | 7.3e-6 | 9.1e-7 | 3.3e-6 | 2.1e-6 |
| (b) m=2, g=10, phi_*=0.5 | 0.228 (KLS 0.254) | 0.228 | 1.42e-2 | 8.1e-2 | n/a |
| (b) m=2, g=30, phi_*=0.5 | 0.605 (KLS 0.634) | 0.605 | 0.127 | 1.95 | n/a |
| (b) m=2, g=5, phi_*=phi_inf=1 | 8.1e-3 | 8.1e-3 | 6.5e-4 | 1.8e-3 | 2.6e-4 |

* (a): production is exponentially small for m > 2H and spectrally cut off at k/a_* ~ m (N_k falls below 1e-3 N_IR by kappa ~ 2-3). The benchmark envelope `1/sinh^2(pi tau_p omega)` with the fitted
  tau_p is the right order of magnitude for m >= 2 (within a factor 2-10) but the conditions of the benchmark (|H| tau_p << 1, equal in/out masses) are violated here: H tau_p ~ 0.8, and the
  IR modes are dominated by the dS -> dS change of the effective frequency, not by a pulse. rho_excess <= 9.3e-4 H^4 sits below the Chat 8 ceiling 9.3775e-4/tau_p^4 = 2.0e-3 (tau_p=0.826) to 2.7e-3 H^4 (tau_p=0.767).
* (b) with phi_* inside the trajectory is a Kofman-Linde-Starobinsky crossing (hep-ph/9704452), `N_k = exp[-pi (k^2/a^2 + m^2 - 9H^2/4 - 3Hdot/2)/(g |phidot|)]`, `g|phidot_J| = 0.579 g H`; numerics agree to 5-10% in the IR
  and n = (g phidot)^{3/2}/(8 pi^3) e^{-...} = 0.0143 vs 0.01423. This channel is NOT bounded by the sech^2 ceiling (different profile: m_eff passes through a minimum and ends heavier).
* **Energy fraction.** rho_tot = 3 M_pl^2 H^2, so `rho_prod/rho_tot = (C/3)(H/M_pl)^2` with C = rho_excess/H^4 and `(H/M_pl)^2 = t v0/3` = 7.77e-5 at t = 1e-3:
  minimal coupling C <= 1e-3 -> <= 2.4e-8 (t=1e-3), 2.4e-11 (t=1e-6); g = 30 H: C ~ 2 -> 5e-5 (t=1e-3). Reaching O(1) needs g ~ 1e3 H_J at t = 1e-3 (g ~ 14 in 5D units) and more for smaller t [EST, no backreaction].
  On side (i) the produced quanta are in any case diluted by the continuing de Sitter expansion. Side (ii) has no asymptotic out-region inside the EFT; not computed.

## 5. Reheating ledger for modified detunings (NOT the registered tension)

* **(A) the suggested V = t(1-phi)(1+b phi), b = 1+c [EFT]:** hilltop mu^2 = -29.9 (s = 4.17, only 1.2 e-folds to roll-off at t=1e-3). Because `1-phi_b ~ e^{2y_b}` while `DeltaTheta ~ e^{y_b/9}`,
  `V_E ~ DeltaTheta^18` at the end point (fitted exponent 17.6 on the last 1e-7 of phi): **the modulus mass at phi_b = +1 is zero and it does not oscillate there.** Acceleration ends (eps_H = 1 at tau = 1.55, peak K/V0 = 0.256), the
  modulus goes into kination (w = 1) and reaches the field-space boundary in finite time (tau = 14.06) with diverging rapidity: |sinh psi| > 0.3 at tau = 13.2 (y_b = -33) for t = 1e-3, tau = 14.03 (y_b = -65) for t = 1e-6. There the EFT fails [UNK beyond].
  Reheating would have to be gravitational (Ford 1987; Peebles-Vilenkin astro-ph/9810509), rho_r ~ 1e-3..1e-2 H^4 [EST].
* **(B) interior Minkowski minimum, V = t(1-phi/phi_m)^2(1+b phi), b = c + 2/phi_m (hilltop kept at phi_b = 0) [EFT + EST]:** `m^2 = V''(phi_m)/(f^2 Z_E)` (table vs closed form agree to 1e-3), alpha = -(1/2) dln f/dTheta.

| phi_m | mu^2 hilltop | m/H_top | alpha | Gamma/H_top | T_rh/M_pl | BBN (T_rh > 4 MeV): H_top >, m >, t > |
|---|---|---|---|---|---|---|
| 0.5 | -207 | 7.68 | 0.300 | 40.8 (H/M_pl)^2 | 3.19 (H/M_pl)^{3/2} | 1.6e4 GeV, 1.2e5 GeV, 5e-28 |
| -0.5 | -141 | 16.5 | 0.126 | 71.9 (H/M_pl)^2 | 4.24 (H/M_pl)^{3/2} | 1.3e4 GeV, 2.1e5 GeV, 4e-28 |
| 0.9 | -77 | 1.08 | 0.312 | 0.12 (H/M_pl)^2 | 0.17 (H/M_pl)^{3/2} | overshoots the shallow minimum to phi_b -> 1; not a viable example |

  With Gamma = alpha^2 m^3/M_pl^2 (O(1) prefactor not computed) and T_rh = 0.5 sqrt(Gamma M_pl): t = 1e-3 gives T_rh ~ 6e15 GeV (phi_m = 0.5), t = 1e-10: 4e10 GeV-class, t = 1e-20: ~1e3 GeV, t = 1e-28: fails BBN.
  This is the usual cosmological moduli bound m >~ 1e5 GeV (Coughlan et al. 1983; Banks-Kaplan-Nelson hep-ph/9308292; de Carlos et al. hep-ph/9308325; T_rh >~ 4 MeV: Hannestad astro-ph/0403291, de Salas et al. 1511.00672, Hasegawa et al. 1908.10189).
  The modulus oscillates (eps_H up to 3, i.e. kinetic domination at each crossing) for (2/3) ln(H/Gamma) ~ 4-42 e-folds of matter-like expansion before decaying. If the same phase had to produce the CMB perturbations, H < 1.9e-5 M_pl (r < 0.036, BICEP/Keck 2110.00483) would need t < 5e-9 - but it cannot, see below.
* **What is NOT solved [EFT, firm]:** there is no slow-roll phase. Accelerated expansion before roll-off is 2.7-9.6 e-folds for t = 1e-2..1e-12 (registered detuning), 0.5-1.2 for the modified potentials; ~60 are needed.
  Horizon and flatness problems are not solved and no near-scale-invariant spectrum is produced (|mu^2| >= 7.7 gives a strongly blue/tachyonic spectrum). For the closed (k=+1) shell nucleated at a = 1/H: on side (i) of the registered detuning it
  does not recollapse only because it keeps inflating in the residual dS (which is not our universe); with a V = 0 minimum the curvature term wins almost immediately: numerically a_max/a_bounce = 1.25 (phi_m = 0.5), 2.06 (phi_m = -0.5), total lifetime ~ 4-6.5/H_top
  (dust estimate a_max = rho a^3/(3 H_top^2) checked against a direct integration: 1.1017 vs 1.1037). A long slow-roll stage (the 1%-tuned quadratic detuning of R5) or some other mechanism is indispensable.

## 6. Grading summary

* [EFT, numerically integrated, checks pass] sections 0-3, the tables, the ledger, case A/B dynamics, closed-shell recollapse.
* [EFT + standard QFT, validated code] section 4 spectra; [EST] energy fractions for large g, H/2pi initial amplitude, inhomogeneous fragmentation, Gamma and T_rh prefactors.
* [argued, not proved] regularity of the EFT expansion through phi_b = -1 beyond two derivatives (only the two-derivative coefficients Z, f, V were shown regular).
* [UNK] the 5D nonlinear fate on the throat side; the fate at the AdS_+ boundary in case A; all effects of the non-local (CFT/KK continuum) part of the effective action; O(t) corrections (R4 suggests ~2e-3 relative at t = 1e-3).
* Prior art: moduli-space/low-energy braneworld actions (Brax-van de Bruck-Davis-Rhodes, Kanno-Soda), fast-roll inflation (Linde hep-th/0110195), negative-potential cosmology (hep-th/0202017), KLS preheating (hep-ph/9704452), Bernard-Duncan/Birrell-Davies exact Bogoliubov problems, moduli problem literature. Nothing here is claimed as externally novel; the explicit numbers for the registered model and findings 1(a), 1(b) are new to this project. References quoted from memory, arXiv numbers not re-fetched in this session.
