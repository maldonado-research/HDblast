# HDBLAST Chat 9 - tensor sector, 4D gravity calibration, KK gap, shell-inflation test

Status tags: [P] hand-derived, [N] floating-point numerics (RK4, step-halving; not interval arithmetic), [C] conjecture/estimate.
All numbers: registered shell, t = 1e-3, c = c_star, kappa_5 = 1, unit length L_0. Scripts and raw outputs are in this folder.

## 1. TT graviton sector  (tensor_spectrum.py, tensor_continuum.py, tensor_amplitude_tscan.py)

[P] ds^2 = dy^2 + rho^2 (gamma_mn + h_mn), h TT on unit dS_4, (Box_gamma - 2) h = m^2 h (m in units of H = 1/rho_b):

    h'' + 4 (rho'/rho) h' + (m^2/rho^2) h = 0,     h'(y_b) = 0  (Z2 shell, tension terms cancel against the background junction).

Cone exponents h ~ y^a, a = -3/2 +- sqrt(9/4 - m^2); normalisable (int rho^2 h^2 dy < inf) only for the + branch, m^2 < 9/4.
Continuum threshold m^2 = 9H^2/4 (same threshold as the scalar sector R1). With dz = dy/rho, u = rho^{3/2} h:

    -u_zz + V1 u = m^2 u,   V1 = (9/4) rho'^2 + (3/2) rho rho'' = 9/4 + rho^2 [(15/4) F - phi'^2/2],  F = phi'^2/12 - U/6,
    -d_z^2 + V1 = Q^+ Q,  Q = d_z - (3/2) rho',  shell condition = (Q u)(z_b) = 0.

Hence m^2 = |Qu|^2/|u|^2 >= 0: no tensor tachyon, the zero mode h = const is the ground state. [P]
Massive bound states of Q^+Q map to Dirichlet bound states of the partner QQ^+ with
    V2 = 9/4 + rho^2 [ (9/16) phi'^2 - U/8 ]   ( -U/8 = W^2/12 - W_phi^2/16 ).
[N] On the registered background V2 - 9/4 >= 2.3e-5 > 0 everywhere (minimum at the first grid point y = 0.01, growing to 3622 at the shell)
=> no discrete massive tensor state below the continuum [P given the numerical positivity]. Direct shooting confirms it:
D(m^2) = rho_b h'/h at the shell is strictly monotone on m^2 in [-60, 2.2499] (670 samples), single zero at m^2 = 0, h nodeless, and
dD/dm^2|_0 = -0.013115 = -I_kept/rho_b (exact identity from (rho^4 h')' = -m^2 rho^2 h).
Note V1 is NOT a positive volcano here: the shell sits in the wall core and V1(shell) = -515 (thick-wall well) plus the delta well -3 rho'_b = -78.9.

Zero mode / Planck mass. [P] Reducing (1/2) int sqrt(-g) R on both Z2 copies with induced metric g_4 = rho_b^2 gamma gives
    M_4^2 = 2 int_0^{y_b} (rho/rho_b)^2 dy = 2 I_kept      (kappa_5 = 1).
[N]  I_kept = 1.03386639 (dv = 2e-4 vs 1e-4 agree to 1.5e-9),  M_4^2 = 2.0677328,  2 I_+ = 2.0715425,  8 pi G_N = 1/(2 I_kept) = 0.483621 (Theorem 4: 0.48362).
Theorem-4 sum rule re-evaluated independently: dlam/6 = 1.66668937e-4 = h I_kept (1.66380127e-4) + quadratic (2.8881e-7), residual 1e-13.
Remark: h = dlam/(6 I_+) x 1.001842 x (1 - 1.733e-3) - the zero-mode enhancement and the quadratic sum-rule term almost cancel, which is why eta_M462 t is only 1.2e-4.

KK continuum [N]: w(k) = |u_k(shell)|^2/|u_0(shell)|^2 per dk, m = H sqrt(9/4 + k^2). Monotone, no resonance (d ln w/d ln k in [0.02, 2.08]):
w -> 7.89e-5 k^2 at threshold (strong suppression of light KK modes), RS-like w ~ l_eff^2 H^2 k only in a narrow window (l_eff^2 peaks at 1.28 L_0^2 at k ~ 7;
the window H << m << k_- barely exists at H/k_- = 0.0228), w -> (2/pi) I_kept H = 8.35e-3 (flat 5D) for m L_0 >> 1. Step-halving differences < 1e-6 for k <= 150.

## 2. High-energy tensor-amplitude correction (LMW)

P_T = 2 H^2/(pi^2 M_4^2(H)), so relative to the low-energy limit F^2 = I_+/I_kept. [P for pure AdS: I_kept = (1/2k)[sqrt(1+x^2) - x^2 asinh(1/x)], x = H/k, which is exactly LMW.]
[N] registered shell: F^2 = 1.001842;  LMW with k_- = 5/9 (kept-side throat): 1.002076;  LMW with k_+: 1.0319 (wrong side, for reference only).
t-scan (t = 1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5): (F^2 - 1)_shell/(F^2 - 1)_LMW(k_-) = 0.894, 0.890, 0.888, 0.885, 0.884, 0.883 (t = 1e-5: float Newton solver fails).
[C, numerically supported] thick-wall asymptotics  F^2 - 1 = H^2/(2 k_-^3 I_+) [ ln(1/H L_0) + C ],  1/(2 k_-^3 I_+) = 2.8153 (pure AdS: 1/k^2 = 3.24, C = ln 2k - 1/2);
fitted C = -0.268, -0.291, -0.301, -0.305, -0.307, -0.301 (C ~ -0.30 to -0.31; the last point is limited by the accuracy of the float shell solve, so the asymptotic coefficient 2.8153 is supported, not proved). The tensor correction is a 0.18% effect: negligible.

## 3. Physical calibration (physical_calibration.py)

Units restored: S = M_5^3 [ int (R/2 - (d phi)^2/2 - U/L_0^2) - int sigma_t/L_0 ].  M_Pl^2 = 2 I_kept M_5^3 L_0,  H = 0.01268582/L_0 (t = 1e-3 fixed)
=>  L_0 = 0.0126858/H,  M_5^3 = M_Pl^2 H/(2 I_kept x 0.0126858).  Fixed ratios: H/k_- = 0.02283, H/k_+ = 0.1142, KK gap 1.5 H, tachyon |m| = 2.778 H, growth time 1/(sH) = 0.6034/H.

| case | H [GeV] | L_0 [m] | M_5 [GeV] | M_5 L_0 | k_- [GeV] | k_+ [GeV] | tension^(1/4) [GeV] | detuning^(1/4) [GeV] | growth time [s] |
|---|---|---|---|---|---|---|---|---|---|
| A | 1e13 | 2.50e-31 | 1.31e17 | 166 | 4.38e14 | 8.76e13 | 4.35e16 | 6.50e15 | 3.97e-38 |
| B | 1e10 | 2.50e-28 | 1.31e16 | 1.67e4 | 4.38e11 | 8.76e10 | 1.37e15 | 2.06e14 | 3.97e-35 |
| C | 1e5 | 2.50e-23 | 2.83e14 | 3.59e7 | 4.38e6 | 8.76e5 | 4.35e12 | 6.50e11 | 3.97e-30 |
| D: l_- = 50 um | 9.01e-14 | 2.78e-5 | 2.73e8 | 3.8e19 | 3.95e-12 | 7.9e-13 | 4.1e3 | 617 | 4.4e-12 |
| D': l_- = 30 um | 1.50e-13 | 1.67e-5 | 3.24e8 | 2.7e19 | 6.58e-12 | 1.3e-12 | 5.3e3 | 796 | 2.6e-12 |

Detuning energy density t(1 + c phi_b) M_5^3/L_0 = 3 H^2 M_Pl^2 (1 + 1.7e-3) in every row (check column in the script output).
Laboratory bound: gravitational-strength Yukawa excluded above 38.6 um (Lee et al., arXiv:2002.11761, abstract verified). With l_- = 1.8 L_0
this requires L_0 < ~ 17-28 um, i.e. for the registered t: H > ~ 1e-13 GeV and M_5 > ~ 3e8 GeV (the usual RS2 number). Cases A-C satisfy it by > 17 orders of magnitude.
Caveats: (i) which length controls the 1/r^3 correction for this thick asymmetric wall (l_-, l_eff = 1.68 L_0, ...) is an O(1) ambiguity not resolved here;
(ii) the lab bound concerns today's universe, whereas H above is the (unstable) shell's Hubble rate; today's H_0 would require t ~ 6 I_+ (H_0 L_0)^2, i.e. the usual cosmological-constant tuning;
(iii) M_5 L_0 >> 1 and k/M_5 << 1 in all rows, so classical 5D gravity is self-consistent; (iv) the low-energy condition H << k_- holds (0.0228) but H/k_+ = 0.114.

## 4. Shell-modulus inflation with the quadratically tuned detuning (shell_inflation.py, shell_inflation_B.py, inflection_scan.py)

EFT of R3, V_E/t = (1 + c phi + d phi^2/2)/f^2, canonical Theta (NB: in eft_landscape.npz Theta > 0 is the AdS_- throat side, opposite to the comment in eft_landscape.py).
Exact homogeneous equation in e-folds, Hubble-flow parameters, second-order n_s.

(A) Registered family c = c_star, d near d0 = 1.1134966 (orchestrator's proposal d = 1.105924, eta_V = -0.0175).  [N]
  * The hilltop is strongly asymmetric: gamma = V'''/V = +7.67 (canonical field). Toward the throat V has a dS minimum at Theta = 2|eta|/gamma = 0.0045 (trapped);
    the roll goes toward phi -> +1. With the cubic term, n_s - 1 = -2|eta| coth(|eta| N/2) <= -4/N  [P within slow roll], NOT 1 + 2 eta.
  * Numerics: d = 1.105924: n_s = 0.910 / 0.918 / 0.924 at N = 50/55/60, r = 7e-8, alpha_s = -1.2e-3. Best case d -> d0: n_s = 0.915 / 0.9235 / 0.930.
    => excluded by Planck n_s = 0.9649 +- 0.0042 at 8-10 sigma (N = 60-55). The R5 statement "n_s ~ 0.965 at d = 1.105924" is not correct once V''' is included.
  * Inflation never ends: max eps_H = 0.50 (max eps_V = 0.57) on the phi -> +1 side, where the shell becomes a detuned RS brane in AdS_+ (V_E -> t(1+c+d/2)/81 > 0).
    N above is counted back from the point of maximal eps_H (a nominal end).
(B) Inflection variant: c = c_star + dc, d = d0. V'/V = lam1 + gamma Theta^2/2, n_s - 1 = -4a cot(aN), a = sqrt(lam1 gamma/2).  [N]
  * n_s(55) = 0.9649 needs dc = -9.44e-5 (relative 1.6e-4; lam1 = 1.43e-4, N_tot ~ 133). Then n_s = 0.952/0.965/0.977 at N = 50/55/60 (strong N-dependence),
    r = 2.1e-7 (<< 0.036), alpha_s = -2.6e-3 (Planck: -0.0045 +- 0.0067), H = 1.2e11 GeV.
  * Planck 1 sigma corresponds to dc within +-10% (absolute window ~ 1e-5): a ~2e-5 tuning of the slope. d is NOT separately critical (d0 +- 3e-3 moves n_s by 1e-3);
    inflection_scan.py shows a whole curve of exact inflection points (c, d)(phi_i), always with gamma > 0, always rolling to phi -> +1, max eps_V ~ 0.5-0.6: no graceful exit anywhere in the quadratic family.
  * A_s = 2.1e-9 fixes only t/(M_5 L_0)^3 = 2.94e-14 (N = 55) (the overall action normalisation is the second free parameter), and H/M_Pl = 4.8e-8.
    With M_Pl^2 = f_* M_5^3 L_0 (f_* = 2.074; O(1) ambiguous because f keeps evolving to 9):
      t = 1e-3: M_5 L_0 = 3241, M_5 = 2.97e16 GeV, L_0 = 1.09e-13 GeV^-1 = 2.15e-29 m, k_- = 5.1e12 GeV;   t = 1e-2: M_5 = 2.0e16;  t = 1e-4: M_5 = 4.4e16;  t = 1e-6: M_5 = 9.4e16 GeV.
  * Honest verdict: a tuned, exit-less toy. O(t) corrections to the EFT (R4: delta mu^2 ~ 1.9 t, i.e. delta eta ~ 6e-4 at t = 1e-3) are comparable to the tuning width, so (c, d) would have to be re-tuned in the full 5D problem. No reheating mechanism. n_s numerics carry ~3e-4 noise from interpolated third derivatives.

## 5. Prior art
Zero mode + gap 3H/2 + F(H/k): Langlois-Maartens-Wands hep-th/0006007 (abstract verified); gap for dS branes: Garriga-Sasaki hep-th/9912118 (recalled, not re-fetched);
SUSY-QM factorisation of the graviton operator: DeWolfe-Freedman-Gubser-Karch hep-th/9909134 (flat walls; recalled); dS thick branes: Kobayashi-Koyama-Soda hep-th/0107025 (recalled).
Inflection-point n_s = 1 - 4a cot(aN): standard (e.g. Baumann et al. 0706.0360, recalled, not re-verified). Nothing in Sections 1-3 is claimed as new physics; the numbers are model-specific.
