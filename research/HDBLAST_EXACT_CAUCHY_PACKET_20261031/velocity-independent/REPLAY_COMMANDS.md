# Independent velocity and conservative response replay

The inventory in MANIFEST.json is the essential portable package. Preserve its relative layout. INDEP denotes this package, REPO the authenticated repository containing the retained research sources, NATIVE the native velocity directory containing producer.py, FREEZE.json and INITIAL_NORMAL.json, ABS the independent absolute-refinement authority directory, CAUCHY the final native Cauchy response JSON, and OUT a fresh output directory. Each output filename and each control directory must be absent before a run.

Run each scientific checker once with python and once with python -O, using separate output names:

```sh
python "$INDEP/scripts/check_velocity_independent.py" --source "$REPO" --out "$OUT/INITIAL_NORMAL.json" --coverage
python -O "$INDEP/scripts/check_velocity_independent.py" --source "$REPO" --out "$OUT/INITIAL_OPTIMIZED.json" --coverage
python "$INDEP/scripts/check_velocity_regularity.py" --source-root "$REPO" --producer-root "$NATIVE" --pins "$INDEP/regularity-source-pins.json" --output "$OUT/REGULARITY_NORMAL.json"
python -O "$INDEP/scripts/check_velocity_regularity.py" --source-root "$REPO" --producer-root "$NATIVE" --pins "$INDEP/regularity-source-pins.json" --output "$OUT/REGULARITY_OPTIMIZED.json"
python "$INDEP/scripts/portable_velocity_audit.py" --source "$REPO" --own "$INDEP/INITIAL_NORMAL.json" --native "$NATIVE/INITIAL_NORMAL.json" --producer-root "$NATIVE" --output "$OUT/AUDIT_NORMAL.json"
python -O "$INDEP/scripts/portable_velocity_audit.py" --source "$REPO" --own "$INDEP/INITIAL_NORMAL.json" --native "$NATIVE/INITIAL_NORMAL.json" --producer-root "$NATIVE" --output "$OUT/AUDIT_OPTIMIZED.json"
python "$INDEP/scripts/portable_conservative_response.py" --source "$REPO" --absolute-dir "$ABS" --native-response "$CAUCHY" --output "$OUT/CONSERVATIVE_NORMAL.json"
python -O "$INDEP/scripts/portable_conservative_response.py" --source "$REPO" --absolute-dir "$ABS" --native-response "$CAUCHY" --output "$OUT/CONSERVATIVE_OPTIMIZED.json"
```

Each portable wrapper calls the independent rational checker in the same Python optimization mode and replaces provenance path keys with stable input roles. Its output excludes temporary filenames and private paths. The default independent directory and portable conservative freeze are relative to the wrapper; --independent-dir and --freeze can override them.

The conservative wrapper supplies these exact role assignments to the equation engine:

```text
own_initial=INDEP/INITIAL_NORMAL.json
own_regularity=INDEP/REGULARITY_FINAL_NORMAL.json
own_audit=INDEP/AUDIT_PORTABLE_FINAL_NORMAL.json
absolute=ABS/REFINEMENT.json
baseline=REPO/research/HDBLAST_EXACT_PACKET_INITIAL_20261031/reference/RESPONSE.json
baseline_initial=REPO/research/HDBLAST_EXACT_PACKET_INITIAL_20261031/reference/INITIAL.json
native_response=CAUCHY
```

The engine also accepts these directly as repeated --input role=path arguments with --freeze INDEP/CONSERVATIVE_RESPONSE_PORTABLE_FREEZE.json. Overrides only select a location; every frozen SHA remains required. The native response role is a consistency comparison with the smaller native bounds, while the conservative proof uses its own independently reconstructed larger velocity bounds.

Run each control driver once; each launches both normal and optimized rejection subprocesses internally:

```sh
python "$INDEP/scripts/run_velocity_controls.py" --source "$REPO" --producer-root "$NATIVE" --output-dir "$OUT/velocity-controls" --receipt "$OUT/VELOCITY_CONTROLS.json"
python "$INDEP/scripts/run_conservative_response_controls.py" --source "$REPO" --absolute-dir "$ABS" --native-response "$CAUCHY" --output-dir "$OUT/response-controls" --receipt "$OUT/RESPONSE_CONTROLS.json"
```

The declared checks are 48 active interval cells at 128 bits, plus the exact zero cell; four byte-identical normal/optimized report pairs; 14 initial/audit corruption cases with 28 rejection subprocesses; and 10 conservative-response corruption cases with 20 rejection subprocesses. Every rejection receipt records success_artifact_absent=true. The 192-bit fallback was not needed. The native initial proof has 49 leaves including its zero cell.
