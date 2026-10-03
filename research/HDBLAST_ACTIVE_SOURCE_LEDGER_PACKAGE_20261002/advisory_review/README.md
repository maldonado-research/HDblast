# Independent advisory result review

These tools review completed serialized active-source diagnostic outputs after the public freeze. They execute no numerical producer, evaluate no physical source, and interpret no checkpoint arrays. The registered validator remains authoritative. This directory is external advisory tooling, not an amendment to the frozen numerical methods or acceptance gates.

`review_outputs.py` recomputes signed accounting identities, same-case witness margins, all 144 cross-control/context pairs, separated source and ledger controls, three discrete contact/baseline knots, and cancellation envelopes. Exact decimal strings are converted to rational numbers. Differences are reported without changing an authoritative scientific classification.

Supply the completed original results, both validator outputs, the outer execution receipt, and explicit freeze pins:

```bash
python review_outputs.py --review-after-public-freeze \
  --checkpoint-root /path/to/authenticated-checkpoint \
  --freeze-commit PUBLIC_COMMIT_40_HEX \
  --registration-sha256 REGISTRATION_64_HEX \
  --primary /path/to/primary.json \
  --independent /path/to/independent.json \
  --validation /path/to/validation-normal.json \
  --validation-optimized /path/to/validation-optimized.json \
  --execution /path/to/execution.json \
  --output-dir /path/to/new-external-review-directory
```

The output directory must be new and outside the checkpoint. It receives the original raw files, their hashes, `RESULT_REVIEW.json`, and `CASE_ATTRIBUTION_MARGINS.csv`. A consistent advisory audit exits zero; an accounting or replay difference exits one. Scientific `CONSISTENCY_FAILURE` can still have a consistent advisory audit: the audit should preserve and explain a valid negative result.

To compare a complete fresh packaged replay, add all three fresh JSON paths with `--fresh-primary`, `--fresh-independent`, and `--fresh-validation`. If the original independent package receipt is `None` and the fresh one is a manifest hash, also supply `--fresh-package-manifest` pointing to the actual package `MANIFEST.json`. The tool checks that manifest's bytes, every listed payload byte, and the registered source/input/freeze family. A valid-looking hash alone is insufficient. Root supplies separate evidence for the fresh outer execution and remote public verification; this tool does not establish either from result JSON alone.

The explicit science-only allowlist is accompanied by a stricter stable full-frame digest compatible with the inactive automation monitor's `HDBLAST_ACTIVE_SOURCE_SCIENTIFIC_PROJECTION_V1`. The latter excludes only measured route elapsed time, peak RSS, and the documented independent package-manifest receipt. Resource scopes, limits, frozen provenance, every scientific field and original list order remain retained. Raw original and fresh files preserve all omitted observations. Validator replay comparison additionally excludes its copied result hashes and Python optimization flag while retaining all classifications and failures.

Preparation tests use invented dictionaries and temporary invented package bytes only. `test_scientific_projection.py` checks mutations and digest compatibility; `test_accounting_review.py` checks signed accounting and actual manifest authentication. Their normal and `-O` receipts record zero physical evaluations. The latter is not a physically valid synthetic experiment and is not offered as scientific evidence.

See `REVIEW_CHECKLIST.md` for limits on interpretation. A successful active-interval ledger attribution concerns the registered finite-momentum diagnostic and empirical quadrature controls. It does not repair the historical `FAIL`, certify continuum accuracy or the initial quantum state, establish a complete cosmological model, or demonstrate a higher-dimensional origin of the Big Bang.
