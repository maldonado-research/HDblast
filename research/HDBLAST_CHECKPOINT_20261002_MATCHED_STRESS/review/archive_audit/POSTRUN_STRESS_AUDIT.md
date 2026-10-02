# Independent post-run audit of saved stress modes and histories

Both normal and optimized Python audits passed. No source was resampled and
no response, forcing history or mode evolution was rerun. The audit imports
no producer. Its anchor is public freeze
`4a5dad8dda6d0a57a9cf88cc8b4a9b8feeaea45d`, with full registration
`4a1533fa7bf9df41cc635d76fa61a22b3aa119d8fbf215349fe09cf3ee54c893`.

The producer's four runs took 274.055 seconds and used a peak 102016 KiB.
Their actual longdouble/clongdouble archives retain the evolved u,w amplitudes
without Wronskian projection. This second audit verifies the frozen input and
actual artifact hashes before using them.

The independent operator reconstruction cancels the common plane-wave phase
analytically. It calculates

    A=Re(u)/k, A'=Re(w)/k,
    delta|Dv|^2=(k^2+L^2)A-Im(w)-L A',

then constructs BOTH minimal stress bilinears and their explicit mass
contacts. It never replaces Im(w) by 2kRe(u). This is different arithmetic
from the producer's full complex phase multiplication. It uses the hash-bound
complete subtraction arrays, which the separately frozen validator also
reconstructs independently; this second audit does not rederive those arrays.

All 24 mode snapshots pass. The audit reconstructs 576 quantity sums (eight
quantities, three cutoffs, six observations, four runs), checks 98,304 GL16
polynomial moments, and recomputes all 72 registered coarse/fine Ward
endpoints and 36 ledger refinements. The latter sums disjoint two-step
Simpson panels and then their prefixes, independently of the producer's global
odd/even grouping.

| Comparison | Largest discrepancy |
| --- | ---: |
| Phase-free bare operator versus complex reconstruction | 2.221e-16 |
| Reconstructed combined stress integrand | 1.373e-16 |
| Reconstructed momentum sums versus recorded values | 7.283e-17 |
| Independently grouped Simpson endpoints | 3.952e-19 |
| Fine registered Ward endpoint versus direct density | 1.986e-10 |
| Observation amplitude Wronskian / epsilon | 5.612e-17 |
| Canonical first-order energy variation integral / epsilon | 3.154e-17 |

The fine Ward endpoint discrepancy is below the frozen 2e-6 gate; the largest
coarse/fine ledger difference is 1.270e-7, below its 1e-6 gate. The audit also
reports all even history endpoints descriptively. Their largest fine
discrepancy is 2.462e-8, compared with 4.329e-7 at coarse resolution. Those
extra times are post-run diagnostics, not new registered acceptance points.

The source jets, stress subtraction variations and anomaly vanish exactly
after the pulse, while the direct stress can remain nonzero. For example,
the saved fine K=256 results are:

| Source | eta | a^4 delta_rho/epsilon | a^4 delta_p/epsilon |
| --- | ---: | ---: | ---: |
| positive_B | -2.5 | -0.00439299880 | -0.00337684550 |
| positive_B | -1.5 | -0.00506845118 | 0.000128626007 |
| signed_uB | -2.5 | -0.000712982419 | -0.00117442816 |
| signed_uB | -1.5 | -0.000388137621 | -0.0000993866541 |

These are finite-cutoff coherent linear responses. The separately derived
tail bounds are needed for continuum statements. The near-zero canonical
first-order energy variation does not make the physical minimal stress zero,
and these results do not establish positive radiation energy or heating.

Arithmetic consistency checks introduced by this post-run audit use a 2e-11
absolute allowance plus an extended-arithmetic rounding allowance. They do
not replace the frozen scientific gates. Finite integration accuracy remains
empirical; only the separately derived omitted UV band has an analytic
enclosure. Intermediate histories store directly integrated stresses rather
than all intermediate modes; independent mode reconstruction therefore covers
the six saved observations per run, while global Wronskian maxima remain
recorded producer evidence.

The portable CLI is:

```sh
python audit_saved_stress.py \
  --checkpoint /path/to/checkpoint \
  --modes /path/to/mode-output/results.json \
  --output /path/to/new-audit.json
```

Python 3.12 with NumPy 2.2.6 and longdouble fraction mantissa >=63 is required.
The actual Python patch is recorded. Outputs are never overwritten. The
normal and `-O` executions both passed without failures.

`IR_COHERENCE_AUDIT.md` is a separate analytic-only note: the positive pulse's
quadratic occupation-only variance has an infrared logarithm, even though its
canonical excitation energy is infrared finite. Coherence cancels the leading
variance singularity at fixed finite time. That note supplies no new numerical
second-order response calculation.
