# Read-only audit of the original registered ledger replay

**204/204 recorded-result checks pass.** The original replay completed all four commands: primary, independent, normal validator, and optimized validator. Both routes independently report `LEDGER_ERROR_DEMONSTRATED`, with the same eight attribution witnesses among the twelve fixed cases and no consistency failures. Normal and optimized validation reports agree exactly after excluding the optimization-mode field. The earlier metric status remains `FAIL`.

The audit binds public freeze `06984aa6b142499c592850e6b4afd49c28ced394` and registration SHA256 `e7a8fe5f6991bad9304440b037ce2eefafdfbaf9c011e177eb9b364b3239bba1`. It checks the registration against the immutable public Git blob and source hashes against the recorded registration and replay copy.

The standalone audit reads completed JSON diagnostic profiles/scalars, validation reports and execution receipts. It independently checks all twelve cases at both 80 and 100 digits in each route, full profile maxima, signed ledger identities, recurrence/direct gaps, triangle bounds, exact precision gaps, cross-route fields and classification criteria. Stream report bytes are authenticated against their recorded hashes; their physical values are not recalculated. No NPZ or physical input arrays are loaded, and neither physical routes nor validators are rerun.

| Measured evidence | Result |
| --- | ---: |
| Maximum cross-route common scalar/profile difference | approximately `7.65853e-73` |
| Maximum 100-digit R profile error | approximately `5.59349e-17` |
| Maximum 100-digit P profile error | approximately `5.00208e-17` |
| Maximum 100-digit absolute continuous defect | approximately `1.61797e-17` |
| Maximum 100-digit absolute projected flow defect | approximately `4.51246e-18` |
| Maximum absolute stored ledger defect | approximately `1.02863e-4` |
| Primary outer full-route execution | `92.754913` seconds, `56,832` KiB peak RSS |
| Independent outer full-route execution | `265.775931` seconds, `57,128` KiB peak RSS |

All consistency criteria use the registered `2e-7` science threshold. Arithmetic/precision/cross-route checks retain `1e-12`, and attribution retains `2e-6`. Each complete route stays within its original 900-second and 262,144-KiB outer budgets. The normal and optimized validators also complete within those limits.

Attribution witnesses:

| Source | Resolution | K |
| --- | --- | --- |
| positive_B | coarse | 64, 128, 256 |
| positive_B | fine | 64, 128 |
| signed_uB | coarse | 64, 128 |
| signed_uB | fine | 64 |

`FINAL_RECORDED_RESULT_AUDIT.json` contains all checks, per-case exact rational values and the report/source SHA256 inventory. Its SHA256 is `6a300482b5d85b7663f8f2eb3d3605caca039e6cb90c15ebb590eb0c39457784`. `SUMMARY.json` provides concise values and resources. The audit script SHA256 is `fcb97e1b5564e07726213e2aedbef1232b2bb7113588debbc0376d35a608fab0`.

Reproduce the read-only report audit:

```sh
/workspace/hdblast-cloud-setup/venv-frw/bin/python /workspace/hdblast-research-work/ledger-review-20261002/original-result-audit/audit_recorded_replay.py --output /tmp/ledger-original-recorded-result-audit.json
```

These results support the registered post-support stored-history ledger attribution. Report consistency does not independently establish raw source-data truth, repair the original metric calibration, establish full metric closure, or provide evidence of a higher-dimensional origin.
