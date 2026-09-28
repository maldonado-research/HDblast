# Program map workstream (HDBLAST checkpoint, 27 September 2026)

## Question

What is the current state of the whole HDBLAST program, stated accurately enough to serve a new reader and the final synthesis? This covers:

- the hypothesis;
- the registered model;
- a dated chronology with Zenodo records;
- which claims are established, conditional, withdrawn or open;
- the registered predictions and falsification tests;
- where the key files and scripts are.

## Method

1. **Reading.** I read the overviews of the latest checkpoint (22 Sept 2026) and its reproduction, claim-review and summary files, and the frozen solver. I also read the overviews of Chats 9–14, the registered-shell controls package, the Chat 8 packages (folders 138–144), the v20/v21 Zenodo kits, the M462R1 root report, the M489G notes and the APS package.
2. **Targeted search.** I used grep and find over `new-files/D-Blast 3` to locate:
   - Zenodo record numbers;
   - the 29 Aug 2026 Zenodo provenance audit;
   - the PTA-knee contracts and gate memos;
   - GPD-site text naming record 22922928.
   
   I did not read the very large archives wholesale.
3. **Web check.** One web search was made for record 22922928. It returned no matching index entry; zenodo.org itself is blocked here.
4. **Verification.** `verify_program_map.py` recomputes every number in the map that can be recomputed cheaply. It checks date and phrase anchors in the source files, rebuilds the Zenodo chain, and checks that every quoted path exists. It includes deliberate controls that must fail. Nothing under `new-files/` was modified.

## Results

| Output | Content |
|---|---|
| `PROGRAM_MAP.md` | Hypothesis, model definition, chronology in 3 eras, claims ledger, predictions and tests, Zenodo table, file index, known inconsistencies |
| `program_map_checks.json` | 34 checks, 0 failed (27 Sept 2026); 34 chronology anchors; 21-version Zenodo chain parsed; 63 quoted paths all present |

Selected verified items:

| Check | Value | Label |
|---|---|---|
| c = 2/I₊ − 4/3 | 0.5975949350280132 (difference 0) | recomputed |
| Chat 9 closed-form μ²(t→0), plus the O(t) slope at t = 10⁻³ | −7.719796 and −7.717872, inside the certified (−7.71788, −7.71786) | recomputed |
| Growth rate from μ² = −7.7178716 | 1.6571936 | recomputed |
| Static +1 branch: residual slope of φ_b = 1 − 9cδ/64 | 2.0005 (the wrong coefficient 9c/32 gives 0.999) | recomputed from saved JSON |
| Static +1 branch: residual slope of the H² two-term expansion | 3.0004 (the wrong coefficient −c²/192 gives 1.99) | recomputed from saved JSON |
| H/H₀ = ρ₀/ρ_b | 0.6067217320; correction to H is −7.8685 ppm (the source quotes −7.868927 using its own benchmark value) | recomputed |
| Chat 14 turnaround (three grids) | H₀τ = 5.8944–5.8947, φ_b = −1.9417 to −1.9412 | read from the 22 Sept recompute JSON |
| Determinant intervals | v23 [2.1987, 2.4615]; v24 [2.1228, 2.5374]; Chat 8 forward [2.3290685, 2.3311463]; all positive and mutually consistent | containment check |
| K₂ knee filter x50 | 1.3391391 (ln x50 = 0.292027); a K₁-type filter gives 0.761 | recomputed |
| Fabricated record 22922929, fabricated path | not found / reported missing | controls pass |

Main findings of the map (details and sources are in `PROGRAM_MAP.md`):

- **Established within scope:**
  - the certified registered shell root (M462R1);
  - the certified tachyon of the registered shell;
  - the exact linearized equations;
  - the certified stable shell S₈⁄₅, mode-stable in all sectors;
  - the certified rank-two linear response (v23, v24, Chat 8 forward), which is conditional;
  - the exact constraint-transport identity and the initial-data obstruction;
  - the exact matter-junction algebra, as a proposal.
- **Negative:**
  - no radiation-dominated branch in the registered model's roll-off (Chats 9, 13, 14);
  - shell-modulus inflation fails (n_s ≤ 0.93);
  - the constant φ = 1 endpoint fails the scalar junction;
  - the v7plus Bessel-K PTA pilots FAIL/FAIL (Mar 2026);
  - the single-pulse quantum heating channel is excluded under the stated assumptions.
- **Withdrawn or demoted:**
  - the PHYS-M309 radiative TT certification was withdrawn and M311's 1.33765% KK energy fraction demoted (Zenodo v20, 21367209, 14 Jul 2026);
  - several draft statements in Chats 9, 13 and 14 were corrected.
- **Open:**
  - stability and attraction of the new static +1 branch;
  - the late-time endpoint;
  - constraint-controlled evolution;
  - matter production, thermalization and a radiation era;
  - perturbations and observables of the 5D model;
  - any derivation linking the 5D model to the PTA knee;
  - authoritative PTA data (`WAIT_FOR_RETURN`).

## Limitations

- M-series checkpoints (PHYS-M51 to M489G) are summarized from their titles and status lines. They were not re-audited.
- Publication evidence for v22, v23 and v24 comes from local transcripts and notes; for 22922928 it comes from site data. No live Zenodo API check was possible. The version label, concept membership and file list of 22922928 are unverified.
- The chronology dates for M313–M400 and M408–M453D are inferred from the bracketing dated files.
- Quoted results inherit the scope and premises of their sources. The verifier checks consistency and cheap identities; it does not re-run any certificate or PDE.
- Era-1 screening numbers are reported only as context. They carry the provenance caveats their own sources list.

## Reproduction

```sh
cd /home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927/program_map
python3 -B verify_program_map.py      # writes program_map_checks.json; about 1 s; needs numpy and scipy
```

The script was tested with Python 3 and NumPy 2.3.5 / SciPy 1.16.3. It reads only the source paths named inside it and PROGRAM_MAP.md. Expected output: `{"n_checks": 34, "n_failed": 0, "all_passed": true, ...}`.
