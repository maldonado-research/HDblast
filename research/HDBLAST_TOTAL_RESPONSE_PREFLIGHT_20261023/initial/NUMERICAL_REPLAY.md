# Optional native initial-enclosure replay

Extract or copy the four published capsules at these fixed source versions. Keep every declared member and no extra regular files within their roots:

- P: [HDBLAST_COMMON_GAUGE_FORCING_20261016](https://github.com/maldonado-research/HDblast/tree/958bcff078e4ffb10884350fcd0bad9c5ab84386/research/HDBLAST_COMMON_GAUGE_FORCING_20261016)
- N: [HDBLAST_NUMERICAL_RESPONSE_INPUTS_20261016](https://github.com/maldonado-research/HDblast/tree/958bcff078e4ffb10884350fcd0bad9c5ab84386/research/HDBLAST_NUMERICAL_RESPONSE_INPUTS_20261016)
- G: [HDBLAST_GEOMETRIC_BULK_SOURCES_20261016](https://github.com/maldonado-research/HDblast/tree/958bcff078e4ffb10884350fcd0bad9c5ab84386/research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016)
- C: [HDBLAST_FORCED_RESPONSE_BOUNDS_20261023](https://github.com/maldonado-research/HDblast/tree/958bcff078e4ffb10884350fcd0bad9c5ab84386/research/HDBLAST_FORCED_RESPONSE_BOUNDS_20261023)

Copy DEPENDENCIES.example.json outside the frozen source and replace each value with its extracted capsule root. Run:

```
python -I -S -B initial/certify_initial_representation.py --dependency-map /path/to/deps.json --output /tmp/new-initial.json
python -I -S -B -O initial/certify_initial_representation.py --dependency-map /path/to/deps.json --output /tmp/new-initial-o.json
python -I -S -B initial/check_initial_gate_ledger.py --report /tmp/new-initial.json --output /tmp/new-ledger.json
```

The authenticator uses four fixed manifest hashes embedded outside the dependency capsules, then verifies their complete declared closures before numerical imports and after evaluation. Public-mode path/pin-source provenance differs from the retained local report; compare its exact representation, error_budget, point_preflight, partition and input_numeric_pins fields. No private planning files are needed. This bounded run covers the initial/free-wave candidate only. It does not evolve U, evaluate the observable or rerun earlier source campaigns.

`compare_initial_gate_reports.py --help` provides an exact numerical-field comparison between local and public-mode reports while allowing only the declared dependency-location/pin-source provenance differences.
