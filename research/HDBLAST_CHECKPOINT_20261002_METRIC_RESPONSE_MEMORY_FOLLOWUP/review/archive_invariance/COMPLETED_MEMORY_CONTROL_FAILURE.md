# Completed memory followup: scientific failure retained

All four independent evolutions completed within the unchanged resource budgets, and the unchanged internal Ward gates then failed. The expected-failure wrapper successfully authenticated that outcome; the registered metric calibration remains **FAIL**, with cross-route calibration **NOT ESTABLISHED**. This read-only audit verified all 370 frozen replay-copy inputs against their registered SHA-256 and pinned Git-object blobs, verified all 32 captured fresh output hashes, and independently reconstructed the saved gate failure list. It did not query GitHub remote publication or recompute raw physical operators.

The checkpoint concerns homogeneous linear metric response at `H=1`, physical mass squared `r=2`, `xi=0`, fixed physical scalar and incoming BD state. It does not establish an actual shifted-root retarded matrix or coupled stability.

## Completed coverage and resources

- Frozen plan: 32 commands; attempted: 28; successful: 27, comprising 26 pure checks and the primary producer.
- Primary: all 12 source/time rows and 36 finite-cutoff points completed in 35.668402274 seconds.
- Independent: `positive_B` and `signed_uB`, each coarse/fine; four full evolutions, six observations and 18 finite-cutoff points per run; all 12 combined source/time rows saved.
- Independent producer elapsed time: 322.835271128 seconds, below the 900-second budget.
- Independent peak RSS: 103,400 KiB, below the 262,144-KiB budget.
- Failing command: `independent_modes`, ordinary exit 1; no timeout or resource-budget failure.
- Not reached after this failure: `validation`, `validation_optimized`, `summary`, `figures`. No cross-route calibration, summary/report or plot acceptance is asserted.

The prior high-precision round exceeded the same memory budget at 306,088 KiB after completing its two positive-source runs. The current lower recorded peak accompanies completion of all four runs. It is an observed storage-engineering improvement, supported by unchanged numerical-source ASTs and exact archived-value invariance; this audit does not infer a measured allocation profile or a physical change from those resource records.

## Unchanged internal gate failures

Saved-data gate recomputation produces **59 failures**: **29 Ward endpoint** failures against `2e-6`, and **30 Ward refinement** failures against `1e-6`. There are no observable refinement failures. All saved observables are finite; canonical and physical Wronskian maxima are respectively 2.4776127297494561e-17 and 3.4201904396460039e-17, below the unchanged `1e-10` gate.

The maximum Ward endpoint difference is 1.9552444892255006e-05, including the historical witness `positive_B`, fine, `eta=-1.5`, `K=256` (about 9.776222 times its gate). Maximum coarse/fine Ward-ledger difference is 0.021054270792825219 (about 21054.271 times its gate). This establishes failure of the registered numerical Ward diagnostic. It does not establish physical instability, heating, an inconsistency of general relativity or evidence for the blast hypothesis. A raw-operator diagnostic and analytic Ward review remain separate work.

| Observable | Maximum saved coarse/fine difference |
| --- | ---: |
| `Q0` | 1.0512424264419451e-15 |
| `current` | 5.9674487573602164e-16 |
| `p` | 2.3257507031360092e-11 |
| `p0` | 4.476430077310356e-13 |
| `q` | 5.9674487573602164e-16 |
| `q_prime` | 1.6675549829869851e-13 |
| `q_second` | 1.3959600142499085e-10 |
| `rho` | 3.1464414407267327e-14 |
| `rho0` | 4.0081219593313122e-14 |

Every one of the 36 combined source/time/cutoff gate points, every failed label, every verified input and every captured output pin is preserved in `COMPLETED_CONTROL_AUDIT.json`. The adjacent archival report documents exact equality of all 1,262 preserved positive-source member pairs and separates unused long-double padding differences from value equality.

## JSON consistency and audit history

The original post-run audit imposed `1e-18` absolute consistency when recomputing a saved Ward endpoint difference from separately serialized binary64 operands. That audit-only check failed: the maximum serialization discrepancy is about `1.054e-17`. Its initial source and failure receipt are retained as `audit_completed_memory_control.initial_too_strict.py.txt` and `READ_ONLY_AUDIT_INITIAL_FAILURE.json`. The final audit uses the frozen expected-failure wrapper's consistency rule (`1e-15` absolute, `1e-12` relative), records each residual and additionally requires that every scientific Ward gate classification agree before and after operand recomputation. No producer, physical value or scientific gate was changed. The retained calibration status remains **FAIL** throughout.

## Pins and reproduction

- Public source freeze: `19fde76912af6e2f30d6f55e066b27d88cc7e34b`.
- Registration SHA-256: `1c5bc9b21b34b1d036b7a7d12ccdb6766ac2ba7835ea04182dba868843866607`.
- Independent manifest SHA-256: `71830eb8aaf31bccd2cb628efeb9735490e41a7d7e5b5bffc527177f7c70e075`.
- Completed-control audit SHA-256: `5b4778b6bbabe52aba5635a498cac7d5786062060a6f6fcf154698f4c49985db`.
- Expected-failure completion receipt SHA-256: `59039c7f268a6be3a816a683fa410eebf5b2a5c2ce54ebcae117869ae4fdb1a3`.
- Underlying failed replay ledger SHA-256: `905a14532d96fa7a4833c9f6c876d809cad4e11385a3f69b861882fbc22bc3db`.
- Independent failed producer JSON SHA-256: `6ffd98ab0a18ba954cee7044d7d4f0581721eb8da6085c949be51f3785e16785`.

```bash
/workspace/hdblast-cloud-setup/venv-frw/bin/python audit_completed_memory_control.py \
  --control /workspace/hdblast-research-work/metric-memory-complete-control-001 \
  --repository /workspace/HDblast \
  --checkpoint-relative research/HDBLAST_CHECKPOINT_20261002_METRIC_RESPONSE_MEMORY_FOLLOWUP \
  --output /tmp/COMPLETED_CONTROL_AUDIT.json
```

This audit only streams saved file and Git-object bytes and evaluates arithmetic on already saved JSON values. No new physical evaluation occurs.
