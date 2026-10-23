# Bounded native initial-difference replay

The native producer uses the Python standard library. It evaluates96complete
closed cells of width1/128 on `[1,7/4]`, plus the analytic exact-zero region.
It does not execute a source campaign or B1. A typical96-cell execution takes
about2–3minutes here. Every output path must be outside the read-only input
capsules. Bytecode writing is disabled and the producer writes only its
explicit output file and progress text on stdout.

## Public dependency map

Create a JSON object with these six labels pointing at extracted public
capsule roots or the frozen branch report:

| Label | Complete capsule or input |
|---|---|
| P | `HDBLAST_COMMON_GAUGE_FORCING_20261016` |
| N | `HDBLAST_NUMERICAL_RESPONSE_INPUTS_20261016` |
| G | `HDBLAST_GEOMETRIC_BULK_SOURCES_20261016` |
| C | `HDBLAST_FORCED_RESPONSE_BOUNDS_20261023` |
| T | `HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023` |
| B | Frozen `BRANCH_COMPARISON_REPORT.json` of this new pilot |

P/N/G/C/T values must be the directories containing their MANIFEST.json.
B is the report filename. Paths may be absolute or resolved relative to the
process working directory. The producer's embedded fixed hashes authenticate
all P/N/G/C/T manifest closures before numerical imports and after evaluation.
The selected49-file T closure authenticates its copy of the old authenticator
before that module is loaded. No private planning or publication pin files
are required in public mode. The B report has its own embedded SHA256 pin;
the independent branch replay establishes its theorem dependency and proof.

From the directory containing the producer and checkers:

```
python -I -S -B direct_difference_pilot.py --case branch \
  --dependency-map deps.json --recenter --output fresh.json
python -I -S -B -O direct_difference_pilot.py --case branch \
  --dependency-map deps.json --recenter --output fresh-o.json
python -I -S -B check_difference_ledger.py --report fresh.json \
  --output ledger.json
python -I -S -B compare_difference_reports.py \
  RECENTERED_CANDIDATE_NORMAL.json fresh.json fresh-o.json \
  --output comparison.json
python -I -S -B check_difference_controls.py --output controls.json
python -I -S -B check_authentication_controls.py --dependency-map deps.json \
  --output authentication-controls.json --scratch-dir /tmp/fresh-owned-auth-copy
```

The difference controls write only their requested output. They check exact truncated
jet product/reciprocal identities and25scalar kernel derivative cases using
512-bit directed rational reference calculations. The authentication controls
mutate only owned copied T/N trees in a fresh nonexistent scratch directory,
reject six cache/extra/symlink cases before numerical imports, and remove
those copies. The ledger rejects eight
representative mutations. These controls supplement the analytic all-cell
proof; they are not a replacement for it.

The comparator excludes exactly `current_authentication.external_pin` and
`branch_input_certificate.path`. These location descriptions differ between
native external-pin authentication and public embedded-pin authentication.
All numerical fields, input hashes, candidate coefficients, error rows and
remaining provenance must match. Native normal/-O reports are byte-identical.

Omit `--recenter` to replay the fixed original candidate at the refined actual
beta interval. Its certified free-wave cap is about1.757028e30and FAILS1e28.
With the single registered background recentering, the cap is about8.518093e27
and PASSES1e28. The packet/cutoff radial component exceeds its illustrative
4e27design row; the exact recomputed free-wave expression is the accepted
criterion. Neither result evaluates B1, total U or an invariant observable.

Default mode, without `--dependency-map`, uses the recorded native research
paths and additionally checks their external publication/forward-plan pins.
It is for the retained research workspace. Public mode has no fallback to
those numerical dependency paths when a map has been supplied.
