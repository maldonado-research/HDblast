# Reproduce the continuous coordinate shell sources

The registered scientific programs use only the Python standard library.
The tested runtime is **Python 3.12.14**. The optional plotting recipe uses
the separate pinned `requirements-plot.txt`; plots are visualizations, not
the arithmetic authority for any inequality.

## Authenticated fresh replay

1. Obtain the source ZIP and its separate public
   [source-freeze receipt](https://github.com/maldonado-research/HDblast/blob/main/hdblast/continuation/2026.10.15-coordinate-shell-sources/SOURCE_FREEZE_RECEIPT.json).
   Compare the complete ZIP SHA-256 and byte count with that external receipt.
   Extract into a new directory, preserving every member.
2. Use the receipt's complete 64-character `manifest.sha256` as the expected
   source-manifest identity. From the extracted capsule directory, run:

   ```sh
   python -I -B reproduce.py --manifest-sha256 EXPECTED_64_CHARACTER_SHA256 --output-dir /tmp/hdblast-coordinate-fresh-replay
   ```

   Replace the uppercase placeholder with the actual external pin; do not
   pass the literal placeholder. The output directory must be fresh and
   outside the frozen source. Reading a source's own manifest hash without
   comparing an external pin establishes consistency, not its identity.
3. Inspect `REPRODUCTION_RECEIPT.json` in that output directory. The runner
   authenticates the complete source inventory before and after execution,
   copies every process into a separate fresh workspace, compares complete
   normal/optimized JSON goldens and every generated CSV byte, and requires
   each registered mathematical/input corruption to fail for its exact
   named reason. Process counts and the complete named inventories come
   from the frozen `REPLAY_ROSTER.json`.

The four programs registered in each mode are the primary generator, Astra
independent geometry/sign/norm checker, Sol independent complete-cell
reconstruction, and the independent full-tensor/source-scope checker.
The generator has no enumerated mathematical check groups; its complete
outputs are regenerated. Execution metadata may differ only at explicitly
declared runtime mode fields. A nonzero control exit by itself is not an
accepted mathematical rejection.

## Individual programs

For diagnosis, use fresh output names outside the capsule:

```sh
python -I -B compute_coordinate_sources.py --output-dir /tmp/hdblast-coordinate-generated
python -I -B astra/check_coordinate_sources.py --output /tmp/hdblast-coordinate-astra.json
python -I -B sol/check_coordinate_shell_sources.py --output /tmp/hdblast-coordinate-sol.json
python -I -B referee/check_full_source_formulation.py --output /tmp/hdblast-coordinate-referee.json
```

Repeat with Python's `-O` flag before the program name to verify optimized
execution. Explicit exceptions implement the checks; disabling Python assert
statements does not disable acceptance gates. Each independent checker has
CLI mathematical-error controls listed in the roster and its accompanying
note; the fresh runner also applies the registered copied-data corruptions.

The exact arithmetic authority is
`computed/COORDINATE_SHELL_SOURCE_ENCLOSURE.json` and its complete CSV.
Each row bounds the actual functions on its whole closed time interval.
The displayed decimals and enclosure midlines are not additional numerical
samples or exact trajectories. Intervals may be conservative and need not
prove the sign of a coefficient at every point.

## Optional figure

Install the plotting requirements in a separate environment if desired, and
use a writable cache directory without changing the scientific inputs:

```sh
MPLCONFIGDIR=/tmp/hdblast-coordinate-matplotlib XDG_CACHE_HOME=/tmp/hdblast-coordinate-cache python -I -B plot_coordinate_sources.py
```

This regenerates the PNG/PDF from the certified table and report. It writes
those figure files, so use a working copy if preserving a frozen capsule.
The plot shows interval bands and their midlines, source norm intervals and
a table of separately certified signs. It does not depict a solved U2.

The separate dated publication evidence records the actual fresh replay,
independent reviews, GitHub checks, merges and live readbacks. Those operational
observations do not certify the inherited analytic PDE hypotheses or replace
external peer review. No finite-amplitude nonlinear solution, effective
remainder, observational calibration or cosmological origin follows.
