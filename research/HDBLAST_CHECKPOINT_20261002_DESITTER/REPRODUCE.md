# Reproduce the de Sitter checkpoint

Use Python 3.12 and the pinned packages in this checkpoint's `requirements.txt`: NumPy 2.2.6, SciPy 1.15.3, SymPy 1.14.0 and mpmath 1.3.0. Choose a fresh output directory so original runs remain available for comparison.

The replay requires a full HDBLAST checkout containing this checkpoint and its inherited scientific inputs. A checkpoint-only ZIP is **not self-contained**: `--repo` must identify a checkout containing the pinned SMOOTH_FRW, JUNCTIONS and earlier classical background inputs. The source baseline is `0205cc651bfb614c32e39dfe833d93d957264229`; use the later checkpoint branch or release commit that also contains DESITTER. The replay checks pinned input bytes, rather than requiring the checkout HEAD to equal the older baseline.

From the root of that checkout:

```bash
python3.12 -m venv /tmp/hdblast-desitter-venv
/tmp/hdblast-desitter-venv/bin/python -m pip install -r research/HDBLAST_CHECKPOINT_20261002_DESITTER/requirements.txt
/tmp/hdblast-desitter-venv/bin/python research/HDBLAST_CHECKPOINT_20261002_DESITTER/code/replay_checkpoint.py \
  --repo /path/to/HDblast \
  --output /fresh/output
```

Replace the two absolute example paths with the checkout and an unused output location. The driver reruns the source benchmark, symbolic theory checks, independent proper-time implementation, exact-mode bridge, independent record audit, cross-method summary, and classical sensitivity calculation and verifier. Outputs and logs go to the supplied directory. The exact symbolic mode calculation and high-precision proper-time integrations can take several minutes.

## Expected scientific checks

| Replay component | Expected outcome |
|---|---|
| Primary source matrix | Four points, two resolutions; 57 passing gates; 3 negative controls detected |
| Symbolic theory | 37 assertions and 6 wrong-formula controls pass |
| Proper time | 8 evaluations and 12 passing refinement comparisons |
| Independent audit | 12 cross-method comparisons, 8 trace gates and 5 mutation controls pass |
| Exact-mode bridge | Separate energy, pressure and variance yield `11H⁴/(960pi²)`, `−11H⁴/(960pi²)`, `H²/(12pi²)` at `x=r=2H²` |
| Classical sensitivity | 19 runs, 11 numerical gates and 3 negative controls pass; positive `E1_ell≈9.63825e−14` |

Compare numerical values with the registered tolerances and inspect the reported status. Timestamps, execution paths, run-source hashes for portable wrappers, and platform-sensitive floating-point details need not match archived JSON byte for byte. Analytic truncation bounds and empirical numerical differences are different quantities; a passing replay does not create a global certified error enclosure.

## Provenance to retain

`PREREGISTRATION.json` pins the primary source files and inherited inputs. `REGISTRATION.md` specifies the fixed grid, precision sequence and gates. Do not run the producer's `--freeze` option when replaying an existing registration: that operation creates a new registration rather than reproducing the preserved one.

The primary local freeze occurred at `2026-10-02T06:17:06.731272Z`; public commit `85aea9955b99bba911a0e66e5869c184dd86a660` preceded the primary run beginning at `06:18:58.957725Z`. The review and sensitivity registrations were local records made before their own computations, which began before that public commit. A replay is a new execution; it does not recreate prospective timing.

The independent audit restores eight originally registered trace gates omitted from the original review aggregator's acceptance calculation. This is a retained post-run validation correction. The exact-mode bridge is also a post-registration diagnostic. `independent/INDEPENDENT_SCIENTIFIC_REVIEW.md` explains both.

The sensitivity directory retains the registered executed sources, the first failed symbolic verifier and its log, the corrected verifier, and `PORTABILITY.json`. Portable entry points provide explicit repository/output arguments; the numerical constants, solver sequence, gates and scientific functions retain the documented relationship to the registered sources. The failed verifier is historical evidence and is not part of the passing replay command.

The archive includes literature metadata and search provenance. Replaying the scientific calculations does not refetch papers or reproduce an external peer review. The prepared Zenodo metadata does not imply a deposit or DOI has been created.
