"""Compare the frozen-wall tagged probe across Fourier and wall grids."""
import json
import sys
from pathlib import Path

base, refined = (json.loads(Path(arg).read_text()) for arg in sys.argv[1:])
windows = tuple(base["switches"])
bands = ("full", "screened")


def max_item(test):
    return max(test(grid, arm, row, w, band) for grid, data in base["grids"].items()
               for arm, row in data.items() for w in windows for band in bands)


def clock(key):
    return max_item(lambda g, a, row, w, b: abs(row["clocks"][key]))


def difference(q):
    return max_item(lambda g, a, row, w, b: abs(row[w][b][q]-refined["grids"][g][a][w][b][q]))


def contrast_delta(other, meshes, keys):
    diffs = []
    for g, og in meshes:
        for a, row in base["grids"][g].items():
            compare = other["grids"][og][a]
            for w in windows:
                for b in bands:
                    left, right = (r["contrasts_percent"][w][b] for r in (row,compare))
                    diffs.extend(abs(left[k]-right[k]) for k in keys)
    return max(diffs)


keys = list(base["flat"].keys() & {"number","energy","annihilation_number","stimulated_number",
                                   "annihilation_energy","stimulated_energy","net_parent_loss"})
result = dict(
    grids={k: [r["samples"],r["padding"]] for k,r in (("base",base),("refined",refined))},
    max_abs_fourier_grid_change={k:difference(k) for k in sorted(keys)},
    max_fourier_grid_contrast_change_pp=contrast_delta(refined, [(g,g) for g in base["grids"]], keys),
    max_wall_grid_inclusive_contrast_change_pp=contrast_delta(base,[("fine","coarse")],("number","energy")),
    max_abs_parseval_error=max_item(lambda g,a,r,w,b:abs(r[w]["time_identity"]["spectral_minus_time"])),
    max_abs_wronskian_error=clock("max_wronskian_error"),
    max_abs_conformal_clock_error=clock("final_eta_difference"),
    max_relative_mode_error_between_clocks=clock("relative_final_mode_error"),
    max_screened_daughter_momentum=max_item(lambda g,a,r,w,b:r[w]["max_screened_single_daughter_physical_momentum"]))
if (result["max_fourier_grid_contrast_change_pp"] > 1e-5 or result["max_abs_parseval_error"] > 2e-8
        or result["max_wall_grid_inclusive_contrast_change_pp"] > .01):
    raise ValueError("Resolution checks failed")
print(json.dumps(result, indent=2,allow_nan=False))
