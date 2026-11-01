# Independent weighted initial-monopole review

This component independently checks the source formulas, radial operators,
canonical conversions, compact homogeneous metric data, and strict negative
initial monopole moment in the accompanying
[weighted initial-data note](../WEIGHTED_MEAN_INITIAL_DATA.md).
The exact sign proof and conditional acceptance appear in [REVIEW.md](REVIEW.md).

The component contains 29 exact symbolic/rational checks, including two
controls that must reject deliberately incorrect formulas. It does not solve
the weighted boundary-value problem or evolve a shell response.

## Replay

Requirements are Python 3.9 or later and SymPy. The retained replay used
Python 3.12.14 and SymPy 1.14.0, also recorded in its JSON reports;
[requirements.txt](requirements.txt) pins the tested SymPy version. A source archive with
the six paths in [SOURCE_PINS.json](SOURCE_PINS.json) is sufficient; no Git
metadata, network access, producer checker, or private files are needed.

From this component directory, with `REPO_ROOT` naming that source archive:

```sh
python -I -B checks/independent_review.py --repo-root "$REPO_ROOT" --output "$FRESH_OUTPUT_DIR/NORMAL.json"
python -I -B -O checks/independent_review.py --repo-root "$REPO_ROOT" --output "$FRESH_OUTPUT_DIR/OPTIMIZED.json"
```

`FRESH_OUTPUT_DIR` must be outside both the source archive and this component.
Each output file must not already exist. Omitting `--output`, or using
`--output -`, writes JSON to standard output. The checker never writes into
its component or source tree.

The checker authenticates the six source files by exact byte count and
SHA-256 **before importing SymPy**, then authenticates them again after all
calculations. It also checks that its own code and manifest did not change
during the run. A failed authentication aborts without writing a result.

## Retained evidence

- [Normal replay](results/NORMAL.json) and [optimized replay](results/OPTIMIZED.json)
  were generated in fresh external files against a temporary source closure
  without Git metadata and copied here only after completion.
- [Packaging validation](PACKAGING_VALIDATION.json) records deliberate source
  authentication failures and rejected output locations. The import-order
  checks used isolated Python with site imports disabled, making SymPy
  unavailable; source failures were still reported before any symbolic import.
- [Source pins](SOURCE_PINS.json) list the complete repository-relative source
  closure and source commit provenance. The commit is descriptive; runtime
  verification depends on the pinned bytes.
- [Copy plan](COPY_PLAN.json) gives essential component files, final relative
  destinations, hashes, and reproduction commands for capsule integration.

The original review note hash is retained only as review history. Portable
packaging changes do not alter the 29 mathematical checks, their rejection
controls, or the exact rational sign constants. The input analytic inverse,
background enclosure, and evolution theorems remain conditional premises.
