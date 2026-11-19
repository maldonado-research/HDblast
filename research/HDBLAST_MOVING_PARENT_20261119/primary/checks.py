"""Diagnose the threshold, population and wall resolution in the moving-parent study."""
import argparse
import json
from pathlib import Path

ARMS = ("plus", "minus", "far")
BANDS = ("full", "screened")
WIN = ("wide","narrow")
LABELS = ("q0","q2","q4","q6","q8","uniform","tapered")
POPS = LABELS[-2:]
COEFFICIENTS = ("number","energy","ann_number","ann_energy","stim_number","stim_energy","net_parent_loss")
DENSITY = ("number_density","energy_density","isotropic_pressure")


def contrasts(grid, arm, window, label, band):
    src = grid[arm][window]["spectrum"][label][band]
    quiet = grid["quiet"][window]["spectrum"][label][band]
    return {k:100*(v/quiet[k]-1) for k,v in src.items()}


def compare(left,right,grid_pairs):
    maxima = {"max_relative_number":0.,"max_relative_energy":0.,"max_relative_density":0.,
              "max_contrast_change_percentage_points":0.,"max_inclusive_contrast_change_percentage_points":0.}
    for lg,rg in grid_pairs:
        g1,g2 = left["grids"][lg],right["grids"][rg]
        for arm in g1:
            for w in WIN:
                for label in LABELS:
                    for b in BANDS:
                        r1,r2 = [g[arm][w]["spectrum"][label][b] for g in (g1,g2)]
                        groups={"max_relative_number":("number","ann_number","stim_number","net_parent_loss"),
                                "max_relative_energy":("energy","ann_energy","stim_energy")}
                        if label in POPS:
                            groups["max_relative_density"]=DENSITY
                        for kind,keys in groups.items():
                            worst = max(abs(r1[k]/r2[k]-1) for k in keys)
                            maxima[kind]=max(maxima[kind],worst)
                        c1,c2 = [contrasts(g,arm,w,label,b) for g in (g1,g2)]
                        diff = max(abs(c1[k]-c2[k]) for k in c1)
                        maxima["max_contrast_change_percentage_points"] = max(maxima["max_contrast_change_percentage_points"],diff)
                        diff = max(abs(c1[k]-c2[k]) for k in ("number","energy"))
                        maxima["max_inclusive_contrast_change_percentage_points"] = max(maxima["max_inclusive_contrast_change_percentage_points"],diff)
    return maxima


def result(base,q_only,refined):
    same = (("fine","fine"),("coarse","coarse"))
    out = {
        "parent_quadrature":compare(base,q_only,same),
        "waveform_and_fourier_quadrature":compare(q_only,refined,same),
        "coarse_vs_fine_bulk":compare(refined,refined,[("coarse","fine")]),
    }
    if max(out[k]["max_inclusive_contrast_change_percentage_points"] for k in
           ("parent_quadrature","waveform_and_fourier_quadrature"))>2e-5:
        raise ValueError("Numerical refinement exceeds contrast tolerance")
    signs = {"plus":-1,"minus":1,"far":1}
    verified=0
    for grid in base["grids"].values():
        for arm, sign in signs.items():
            for w in WIN:
                for label in LABELS:
                    for b in BANDS:
                        diff = contrasts(grid,arm,w,label,b)
                        keys = ("number","energy")+(DENSITY[:2] if label in POPS else ())
                        if any(sign*diff[key] <= 0 for key in keys):
                            raise ValueError((arm,w,label,b))
                        verified += len(keys)
    out["pulse_sign_checks"]=verified
    rows=[r for grid in refined["grids"].values() for arm in grid.values() for r in arm.values()]
    for key in ("max_wronskian_error","max_parseval_error","max_clock_error","max_daughter_momentum_for_screen"):
        out[key]=max(r[key] for r in rows)
    out["rest_relative_discrepancy"] = refined["max_relative_error_vs_published_rest"]
    out["fine_population_contrasts_percent"] = {
        w: {label: {b: {arm: contrasts(refined["grids"]["fine"],arm,w,label,b) for arm in ARMS}
                    for b in BANDS} for label in POPS} for w in WIN}
    return out


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reports",nargs=3,type=Path,metavar="JSON")
    args=parser.parse_args()
    print(json.dumps(result(*(json.loads(p.read_bytes()) for p in args.reports)),indent=2,allow_nan=False))
