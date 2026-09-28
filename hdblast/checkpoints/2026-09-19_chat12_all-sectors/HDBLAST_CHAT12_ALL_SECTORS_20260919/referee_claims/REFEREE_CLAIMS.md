# REFEREE_CLAIMS - Chat 12 (all sectors), prior art and safe wording
Date 2026-09-19. Scope: prior-art and claims check only. I read the three docstrings and the task summary; I did NOT re-run the SymPy scripts
(another referee's job). I re-derived identity T3 by hand (below).

## 1. Prior art
Legend: [V] = abstract page checked online today (arxiv.org/abs); [M] = from memory, details not re-verified.

| Ingredient | Prior art | Status |
|---|---|---|
| SUSY-QM factorisation of the graviton operator in an arbitrary 5D scalar-gravity background, "graviton is the susy ground state, no tachyons" | DeWolfe-Freedman-Gubser-Karch hep-th/9909134 (PRD 62 046008) - abstract says exactly this | [V] |
| Universal Schroedinger form for TT modes on thick branes | Csaki-Erlich-Hollowood-Shirman hep-th/0001033 | [V] (abstract; SUSY-QM details [M]) |
| Gap m = 3H/2 between zero mode and KK continuum on a dS brane | Garriga-Sasaki hep-th/9912118 (abstract states the gap); Langlois-Maartens-Wands hep-th/0006007 (abstract: massless mode + gap + continuum) | [V] |
| Same gap for thick dS branes with bulk scalar (V -> 9H^2/4 asymptotically) | e.g. Wang hep-th/0201051 (PRD 66 024024); later thick-dS-brane literature (0709.3552, 0911.0269, 1009.1684). Some of these models DO have a massive bound state below the gap, so its absence is model-dependent, not universal | [V] via search snippets only |
| Bulk-scalar fluctuations on dS brane | Kobayashi-Koyama-Soda hep-th/0009160 (bulk inflaton; abstract does not mention the gap) | [V] abstract; relevance to gap [M] |
| Karch-Randall hep-th/0011156 | concerns AdS4 branes (massive bound graviton), NOT dS. Cite only as a contrast, not for the 9H^2/4 gap | [V] |
| 4D-covariant (dS-slicing) decomposition; scalar/vector/tensor share the same KK spectrum on dS branes; single positive-tension dS brane has no radion; two-brane radion has negative m^2 ~ H^2 | Gen-Sasaki gr-qc/0011078 (abstract confirms all of this; the explicit value -4H^2 is not in the abstract: [M]) | [V] |
| Wall-displacement field with tachyonic mass, l = 1 modes = rigid translations of the dS wall (m^2 = -(N)H^2 for an N-dim dS worldsheet, i.e. -4H^2 for dS_4; translations are not instabilities) | Garriga-Vilenkin PRD 44 1007 (1991): abstract confirms covariant scalar-field description with tachyonic mass; the specific "-4H^2 / l=1 = translation" reading is [M] (also Garriga-Vilenkin PRD 45 3469, [M]) | [V] abstract / [M] details |
| No dynamical vector modes without brane matter sources (graviphoton is projected out by Z2 / is constrained; vector perturbations need anisotropic-stress source on the brane) | standard braneworld perturbation lore: Kodama-Ishibashi-Seto hep-th/0004160, Bridgman-Malik-Wands hep-th/0010133, Maartens Living Rev. gr-qc/0312059 | [M], not checked today |

Is the identity V2 - 9/4 = rho^2 (9 phi'^2/16 - U/8) or "U < 0 along the bulk => no massive graviton bound state below the gap" in the literature?
I did not find it stated (one web search + memory). That is NOT evidence of novelty: it is a two-line consequence of the background equations
(hand check: V2 - 9/4 = (3/4)(rho'^2 - 1) + rho^2 phi'^2/2 using rho rho'' = rho'^2 - 1 - rho^2 phi'^2/3, then the constraint
rho'^2 - 1 = rho^2 (phi'^2/12 - U/6) gives 9 phi'^2/16 - U/8. Identity CONFIRMED by hand.) The method (partner Hamiltonian bounds the excited
spectrum) is textbook SUSY-QM and is how DFGK-type arguments are routinely used. Treat as "elementary corollary, not located in the literature by us".

## 2. Technical remarks on the wording of the scripts (minor, but fix before publishing)
a. T4 says "m2 >= min V2 > 9/4". False as written: rho -> 0 at the cone (z -> -infinity), so inf V2 = 9/4 exactly; there is no positive gap of V2 above 9/4.
   Correct statement: V2 > 9/4 pointwise for y > 0, hence for any normalisable v != 0, <v,QQ^+ v> >= <v,V2 v> > (9/4)<v,v>, so every L2 tensor
   eigenvalue other than m2 = 0 satisfies m2 > 9/4 strictly. This excludes bound states AT or BELOW threshold; it says nothing (and need say nothing) about
   the continuum or embedded states above 9/4. Needs the (standard) facts that Qu in L2 with Dirichlet data at the shell and that boundary terms vanish at the cone end.
b. Special sector: the argument uses B != 0, not B > 0. The claim sentence "none, since B > 0 certified" should read "none, since B != 0 (certified B = +5.3e-4)".
   The registered shell (B < 0) has no special mode either - so this sector does not discriminate stable from unstable shells. Also B is small: a zero of B
   along a family of shells gives a marginal l=1 mode; say so.
c. "No mode" in vector/special sectors is an absence-of-degrees-of-freedom statement under stated assumptions (cone regularity, no brane matter, delta phi = 0 in the
   vector sector by symmetry, Z2). It is not a dynamical stability estimate.
d. The vector docstring works in gauge Fv = 0 with one explicit harmonic family; the generalisation to all transverse harmonics is by covariance (argued, not machine-checked).
e. "Continuum above 9/4" for tensors relies on V1 -> 9/4 at the cone end; fine, but it is units of the dS_4 slicing radius (H_4 = 1), not the physical shell Hubble rate - state units.
f. Whole result is conditional on: Chat 10 interval certificate (existence, phi in (-1,0), B), Chat 11 symbolic scalar equations, and the SymPy checks here. Tensor result additionally only needs U < 0 on the phi-range (exact Sturm count).

## 3. SAFE wording for a Zenodo note
"For the interval-certified shell solution S_8/5 of the 5D Einstein-scalar model (W = 1 - phi + phi^3/3, Z2 shell with tension sigma_t = 2W + t(1 + c phi + d phi^2/2),
d = 8/5), we examined linear perturbations decomposed in harmonics of the dS_4 slices, assuming regularity at the cone apex, Z2 symmetry and no matter on the shell.
(i) Scalar sector, generic harmonics: no normalisable mode below the continuum threshold 9/4 (in units of the slicing curvature), by an interval-arithmetic certificate (Chat 10), with the
linearised equations verified symbolically (Chat 11). (ii) Special (mu^2 = -4, l = 1) scalar harmonics: all bulk perturbations are pure gauge (the regular one is a rigid translation of the apex),
and the junction conditions admit a shell displacement only if B = 0; since B = +5.3187e-4 != 0 (certified), there is no physical mode. (iii) Vector sector: the field equations force
Bv = C/rho^2, singular at the apex unless C = 0; no regular vector mode. (iv) Tensor sector: the TT operator factorises as Q^+Q (no tachyon, in line with DeWolfe et al. hep-th/9909134);
the partner potential obeys V2 - 9/4 = rho^2(9 phi'^2/16 - U/8), and U < 0 on the certified field range (exact Sturm count), so the only normalisable tensor mode with m^2 <= 9/4 is the massless graviton.
Hence we find no exponentially growing or sub-threshold normalisable mode in any sector of this decomposition. Symbolic steps were checked with SymPy; existence and the scalar-sector count are computer-assisted
(interval arithmetic). The decomposition, the 9H^2/4 gap, the SUSY-QM structure and the interpretation of l = 1 modes as translations are standard (Garriga-Vilenkin 1991; Garriga-Sasaki hep-th/9912118;
Langlois-Maartens-Wands hep-th/0006007; Gen-Sasaki gr-qc/0011078; DeWolfe et al. hep-th/9909134; Csaki et al. hep-th/0001033). This is a mode-stability statement for one parameter point; it is not a proof of linear or nonlinear stability."

## 4. Statements that would be OVERCLAIMS (do not use)
1. "Linear stability proved" / "S_8/5 is stable". Mode stability (no bad normalisable separable modes) != linear stability (needs completeness of the mode expansion, a well-posed initial-value
   formulation in the bulk with the cone and shell, and decay/boundedness estimates for the continuum) and says nothing about nonlinear stability.
2. "Stable universe", "stable cosmology", "the HDBLAST universe is stable", or anything about the Big Bang, inflation, reheating, or observations. The object is a vacuum dS_4-sliced shell with no matter;
   Chat 9 already found the roll-off is not a hot Big Bang.
3. "Rigorous theorem / fully proved". It is computer-assisted (interval arithmetic + SymPy), conditional on code correctness, not independently peer reviewed; vector-sector generality is argued by covariance.
4. Any novelty claim: "first", "new mechanism", "new identity". The identity is an elementary corollary of the background equations; we merely did not locate it in print.
5. "Mass gap min V2 > 9/4" or "gap above 9/4 for massive gravitons": inf V2 = 9/4; the continuum starts at 9/4.
6. "Special modes absent because B > 0" as if the sign mattered, or "B > 0 stabilises the l = 1 sector". Only B != 0 is used; the unstable registered shell also has no such mode.
7. "All shells with d > 1.1135 are mode-stable" - only the single point S_8/5 is certified in all sectors (tensor/vector arguments do extend to any background with phi range where U < 0, but scalar certification does not).
8. "4D gravity is recovered / Newton's law on the shell" - a normalisable massless graviton exists; no coupling strength, no KK correction, no phenomenology was computed.
9. "No ghosts / quantum-mechanically stable / semiclassically stable" - nothing about kinetic-sign analysis of the full quadratic action, tunnelling, or the conical apex as a quantum object was done.
10. "Stable with matter on the shell" - brane matter sources vector and anisotropic-stress perturbations; excluded by assumption.
11. Citing Karch-Randall hep-th/0011156 for the dS gap (it is the AdS4 case).
